use crate::poseidon1;
use crate::xits;
use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_bluesky::from_const;
use starkom_ff::Field;
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, make_const, rvar, var,
};
use starkom_poseidon::{BlueSkyConfig3, BlueSkyConfig4};

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height `H`.
///
/// WARNING: `H` must be strictly less than 255. Do NOT use this chip if `H` spans the full BlueSky
/// range, as in that case the bit decomposition of the key would be UNSAFE! Use the
/// [`FullBinaryChip`] below instead.
///
/// The generic argument `L` is the number of lanes (parallel hash stages) used by the chip.
#[derive(Debug, Clone)]
pub struct BinaryChip<const H: usize, const L: usize> {
    decomposer: xits::BitDecomposerChip<H>,
    hasher_ir: poseidon1::PermutationChipIR<BlueSkyConfig3, 3>,
    hasher_er: [poseidon1::PermutationChipER<BlueSkyConfig3, 3>; H],
    path: [[Scalar; 2]; H],
}

impl<const H: usize, const L: usize> Default for BinaryChip<H, L> {
    fn default() -> Self {
        Self::new([[Scalar::ZERO; 2]; H])
    }
}

impl<const H: usize, const L: usize> BinaryChip<H, L> {
    const SELECTOR_HEIGHT: usize = 2;

    pub fn new(path: [[Scalar; 2]; H]) -> Self {
        assert!(L > 0, "need at least one lane");
        assert!(L <= H, "too many lanes");
        let hasher_ir = poseidon1::PermutationChipIR::default();
        let stage_width = hasher_ir.width() as isize;
        let stage_height = (Self::SELECTOR_HEIGHT + hasher_ir.height()) as isize;
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher_ir,
            hasher_er: std::array::from_fn(|i| {
                poseidon1::PermutationChipER::new(
                    ((i + 1) / L) as isize * -stage_height,
                    ((i + 1) % L) as isize * -stage_width,
                )
            }),
            path,
        }
    }

    fn stage_width(&self) -> usize {
        self.hasher_ir.width()
    }

    fn stage_height(&self) -> usize {
        Self::SELECTOR_HEIGHT + self.hasher_ir.height()
    }

    fn build_input_selector(
        &self,
        view: &mut impl CircuitView,
        hash: Option<Cell>,
        bit: Option<Cell>,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.connect(bit, view.cell(0, 2).into());
        view.add_gate(
            0,
            rvar(2, 0) * rvar(1, 0) + (make_const(1) - rvar(2, 0)) * rvar(0, 0) - rvar(0, 1),
        );
        view.add_gate(
            0,
            rvar(2, 0) * rvar(0, 0) + (make_const(1) - rvar(2, 0)) * rvar(1, 0) - rvar(1, 1),
        );
    }

    fn witness_input_selector(
        &self,
        view: &mut impl WitnessView,
        bits: &[CellOrUnconstrained],
        i: usize,
    ) {
        let bit = bits[i];
        let bit_value = match bit {
            CellOrUnconstrained::Cell(cell) => view.get(cell),
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        if bit_value != Scalar::ZERO {
            view.set(view.cell(0, 0), self.path[i][1]);
            view.set(view.cell(0, 1), self.path[i][0]);
        } else {
            view.set(view.cell(0, 0), self.path[i][0]);
            view.set(view.cell(0, 1), self.path[i][1]);
        }
        view.copy(bits[i], view.cell(0, 2));
        view.set(view.cell(1, 0), self.path[i][0]);
        view.set(view.cell(1, 1), self.path[i][1]);
    }
}

impl<const H: usize, const L: usize> PlonkChip<2, 1> for BinaryChip<H, L> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.stage_width() * L)
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.stage_height() * H.next_multiple_of(L) / L
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.build(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        {
            let bit = bits[0];
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, bit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..H {
            let bit = bits[i];
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, bit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 2],
    ) -> Result<[CellOrUnconstrained; 1]> {
        let [key, _] = inputs;
        let bits = self.decomposer.witness(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash;
        {
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &bits, 0)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..H {
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a ternary Sparse Merkle Tree of height `H`.
///
/// WARNING: `H` must be strictly less than 161. Do NOT use this chip if `H` spans the full BlueSky
/// range, as in that case the trit decomposition of the key would be UNSAFE! Use the
/// [`FullTernaryChip`] below instead.
#[derive(Debug, Clone)]
pub struct TernaryChip<const H: usize, const L: usize> {
    decomposer: xits::TritDecomposerChip<H>,
    hasher_ir: poseidon1::PermutationChipIR<BlueSkyConfig4, 4>,
    hasher_er: [poseidon1::PermutationChipER<BlueSkyConfig4, 4>; H],
    path: [[Scalar; 3]; H],
}

impl<const H: usize, const L: usize> Default for TernaryChip<H, L> {
    fn default() -> Self {
        Self::new([[Scalar::ZERO; 3]; H])
    }
}

impl<const H: usize, const L: usize> TernaryChip<H, L> {
    const SELECTOR_HEIGHT: usize = 2;

    pub fn new(path: [[Scalar; 3]; H]) -> Self {
        assert!(L > 0, "need at least one lane");
        assert!(L <= H, "too many lanes");
        let hasher_ir = poseidon1::PermutationChipIR::default();
        let stage_width = hasher_ir.width() as isize;
        let stage_height = (Self::SELECTOR_HEIGHT + hasher_ir.height()) as isize;
        Self {
            decomposer: xits::TritDecomposerChip::default(),
            hasher_ir,
            hasher_er: std::array::from_fn(|i| {
                poseidon1::PermutationChipER::new(
                    ((i + 1) / L) as isize * -stage_height,
                    ((i + 1) % L) as isize * -stage_width,
                )
            }),
            path,
        }
    }

    fn stage_width(&self) -> usize {
        self.hasher_ir.width()
    }

    fn stage_height(&self) -> usize {
        Self::SELECTOR_HEIGHT + self.hasher_ir.height()
    }

    fn build_input_selector(
        &self,
        view: &mut impl CircuitView,
        hash: Option<Cell>,
        trit: Option<Cell>,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.connect(trit, view.cell(0, 3).into());
        let l0 = ((var(3) ^ 2) - var(3) * 3 + 2) / 2;
        let l1 = var(3) * 2 - (var(3) ^ 2);
        let l2 = ((var(3) ^ 2) - var(3)) / 2;
        view.add_gate(
            0,
            l0.clone() * var(0) + l1.clone() * var(1) + l2.clone() * var(1) - rvar(0, 1),
        );
        view.add_gate(
            0,
            l0.clone() * var(1) + l1.clone() * var(0) + l2.clone() * var(2) - rvar(1, 1),
        );
        view.add_gate(0, l0 * var(2) + l1 * var(2) + l2 * var(0) - rvar(2, 1));
        view.add_gate(1, var(3));
    }

    fn witness_input_selector(
        &self,
        view: &mut impl WitnessView,
        trits: &[CellOrUnconstrained],
        i: usize,
    ) {
        let trit = trits[i];
        let trit_value = match trit {
            CellOrUnconstrained::Cell(cell) => view.get(cell),
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        const ZERO: Scalar = from_const(0);
        const ONE: Scalar = from_const(1);
        const TWO: Scalar = from_const(2);
        match trit_value {
            ZERO => {
                view.set(view.cell(0, 0), self.path[i][0]);
                view.set(view.cell(0, 1), self.path[i][1]);
                view.set(view.cell(0, 2), self.path[i][2]);
            }
            ONE => {
                view.set(view.cell(0, 0), self.path[i][1]);
                view.set(view.cell(0, 1), self.path[i][0]);
                view.set(view.cell(0, 2), self.path[i][2]);
            }
            TWO => {
                view.set(view.cell(0, 0), self.path[i][2]);
                view.set(view.cell(0, 1), self.path[i][0]);
                view.set(view.cell(0, 2), self.path[i][1]);
            }
            _ => panic!("invalid trit value {}", trit_value),
        };
        view.copy(trit, view.cell(0, 3).into());
        view.set(view.cell(1, 0), self.path[i][0]);
        view.set(view.cell(1, 1), self.path[i][1]);
        view.set(view.cell(1, 2), self.path[i][2]);
        view.set(view.cell(1, 3), Scalar::ZERO);
    }
}

impl<const H: usize, const L: usize> PlonkChip<2, 1> for TernaryChip<H, L> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.stage_width() * L)
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.stage_height() * H.next_multiple_of(L) / L
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let trits = self.decomposer.build(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        {
            let trit = trits[0];
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, trit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..H {
            let trit = trits[i];
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, trit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 2],
    ) -> Result<[CellOrUnconstrained; 1]> {
        let [key, _] = inputs;
        let trits = self.decomposer.witness(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash;
        {
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &trits, 0)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..H {
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &trits, i)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height 256.
///
/// The keys of such a tree span the full BlueSky range. Internally this chip uses a
/// [`xits::FullBitDecomposerChip`], making the 256-bit decomposition safe at the cost of some extra
/// constraints.
///
/// If you don't need 256- or 255-bit keys use [`BinaryChip`].
///
/// The generic argument `L` is the number of lanes (parallel hash stages) used by the chip.
#[derive(Debug, Clone)]
pub struct FullBinaryChip<const L: usize> {
    decomposer: xits::FullBitDecomposerChip,
    hasher_ir: poseidon1::PermutationChipIR<BlueSkyConfig3, 3>,
    hasher_er: [poseidon1::PermutationChipER<BlueSkyConfig3, 3>; 255],
    path: [[Scalar; 2]; 256],
}

impl<const L: usize> Default for FullBinaryChip<L> {
    fn default() -> Self {
        Self::new([[Scalar::ZERO; 2]; 256])
    }
}

impl<const L: usize> FullBinaryChip<L> {
    const SELECTOR_HEIGHT: usize = 2;

    pub fn new(path: [[Scalar; 2]; 256]) -> Self {
        assert!(L > 0, "need at least one lane");
        let hasher_ir = poseidon1::PermutationChipIR::default();
        let stage_width = hasher_ir.width() as isize;
        let stage_height = (Self::SELECTOR_HEIGHT + hasher_ir.height()) as isize;
        Self {
            decomposer: xits::FullBitDecomposerChip::default(),
            hasher_ir,
            hasher_er: std::array::from_fn(|i| {
                poseidon1::PermutationChipER::new(
                    ((i + 1) / L) as isize * -stage_height,
                    ((i + 1) % L) as isize * -stage_width,
                )
            }),
            path,
        }
    }

    fn stage_width(&self) -> usize {
        self.hasher_ir.width()
    }

    fn stage_height(&self) -> usize {
        Self::SELECTOR_HEIGHT + self.hasher_ir.height()
    }

    fn build_input_selector(
        &self,
        view: &mut impl CircuitView,
        hash: Option<Cell>,
        bit: Option<Cell>,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.connect(bit, view.cell(0, 2).into());
        view.add_gate(
            0,
            rvar(2, 0) * rvar(1, 0) + (make_const(1) - rvar(2, 0)) * rvar(0, 0) - rvar(0, 1),
        );
        view.add_gate(
            0,
            rvar(2, 0) * rvar(0, 0) + (make_const(1) - rvar(2, 0)) * rvar(1, 0) - rvar(1, 1),
        );
        view.add_gate(1, var(2));
    }

    fn witness_input_selector(
        &self,
        view: &mut impl WitnessView,
        bits: &[CellOrUnconstrained],
        i: usize,
    ) {
        let bit = bits[i];
        let bit_value = match bit {
            CellOrUnconstrained::Cell(cell) => view.get(cell),
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        if bit_value != Scalar::ZERO {
            view.set(view.cell(0, 0), self.path[i][1]);
            view.set(view.cell(0, 1), self.path[i][0]);
        } else {
            view.set(view.cell(0, 0), self.path[i][0]);
            view.set(view.cell(0, 1), self.path[i][1]);
        }
        view.copy(bits[i], view.cell(0, 2));
        view.set(view.cell(1, 0), self.path[i][0]);
        view.set(view.cell(1, 1), self.path[i][1]);
        view.set(view.cell(1, 2), Scalar::ZERO);
    }
}

impl<const L: usize> PlonkChip<2, 1> for FullBinaryChip<L> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.stage_width() * L)
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.stage_height() * 256usize.next_multiple_of(L) / L
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.build(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        {
            let bit = bits[0];
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, bit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..256 {
            let bit = bits[i];
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, bit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 2],
    ) -> Result<[CellOrUnconstrained; 1]> {
        let [key, _] = inputs;
        let bits = self.decomposer.witness(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash;
        {
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &bits, 0)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..256 {
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a ternary Sparse Merkle Tree of height 161.
///
/// The keys of such a tree span the full BlueSky range. Internally this chip uses a
/// [`xits::FullTritDecomposerChip`], making the 161-trit decomposition safe at the cost of some
/// extra constraints.
///
/// If you don't need 161-trit keys use [`TernaryChip`].
#[derive(Debug, Clone)]
pub struct FullTernaryChip<const L: usize> {
    decomposer: xits::FullTritDecomposerChip,
    hasher_ir: poseidon1::PermutationChipIR<BlueSkyConfig4, 4>,
    hasher_er: [poseidon1::PermutationChipER<BlueSkyConfig4, 4>; 160],
    path: [[Scalar; 3]; 161],
}

impl<const L: usize> Default for FullTernaryChip<L> {
    fn default() -> Self {
        Self::new([[Scalar::ZERO; 3]; 161])
    }
}

impl<const L: usize> FullTernaryChip<L> {
    const SELECTOR_HEIGHT: usize = 2;

    pub fn new(path: [[Scalar; 3]; 161]) -> Self {
        assert!(L > 0, "need at least one lane");
        let hasher_ir = poseidon1::PermutationChipIR::default();
        let stage_width = hasher_ir.width() as isize;
        let stage_height = (Self::SELECTOR_HEIGHT + hasher_ir.height()) as isize;
        Self {
            decomposer: xits::FullTritDecomposerChip::default(),
            hasher_ir,
            hasher_er: std::array::from_fn(|i| {
                poseidon1::PermutationChipER::new(
                    ((i + 1) / L) as isize * -stage_height,
                    ((i + 1) % L) as isize * -stage_width,
                )
            }),
            path,
        }
    }

    fn stage_width(&self) -> usize {
        self.hasher_ir.width()
    }

    fn stage_height(&self) -> usize {
        Self::SELECTOR_HEIGHT + self.hasher_ir.height()
    }

    fn build_input_selector(
        &self,
        view: &mut impl CircuitView,
        hash: Option<Cell>,
        trit: Option<Cell>,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.connect(trit, view.cell(0, 3).into());
        let l0 = ((var(3) ^ 2) - var(3) * 3 + 2) / 2;
        let l1 = var(3) * 2 - (var(3) ^ 2);
        let l2 = ((var(3) ^ 2) - var(3)) / 2;
        view.add_gate(
            0,
            l0.clone() * var(0) + l1.clone() * var(1) + l2.clone() * var(1) - rvar(0, 1),
        );
        view.add_gate(
            0,
            l0.clone() * var(1) + l1.clone() * var(0) + l2.clone() * var(2) - rvar(1, 1),
        );
        view.add_gate(0, l0 * var(2) + l1 * var(2) + l2 * var(0) - rvar(2, 1));
        view.add_gate(1, var(3));
    }

    fn witness_input_selector(
        &self,
        view: &mut impl WitnessView,
        trits: &[CellOrUnconstrained],
        i: usize,
    ) {
        let trit = trits[i];
        let trit_value = match trit {
            CellOrUnconstrained::Cell(cell) => view.get(cell),
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        const ZERO: Scalar = from_const(0);
        const ONE: Scalar = from_const(1);
        const TWO: Scalar = from_const(2);
        match trit_value {
            ZERO => {
                view.set(view.cell(0, 0), self.path[i][0]);
                view.set(view.cell(0, 1), self.path[i][1]);
                view.set(view.cell(0, 2), self.path[i][2]);
            }
            ONE => {
                view.set(view.cell(0, 0), self.path[i][1]);
                view.set(view.cell(0, 1), self.path[i][0]);
                view.set(view.cell(0, 2), self.path[i][2]);
            }
            TWO => {
                view.set(view.cell(0, 0), self.path[i][2]);
                view.set(view.cell(0, 1), self.path[i][0]);
                view.set(view.cell(0, 2), self.path[i][1]);
            }
            _ => panic!("invalid trit value {}", trit_value),
        };
        view.copy(trit, view.cell(0, 3).into());
        view.set(view.cell(1, 0), self.path[i][0]);
        view.set(view.cell(1, 1), self.path[i][1]);
        view.set(view.cell(1, 2), self.path[i][2]);
        view.set(view.cell(1, 3), Scalar::ZERO);
    }
}

impl<const L: usize> PlonkChip<2, 1> for FullTernaryChip<L> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.stage_width() * L)
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.stage_height() * 161usize.next_multiple_of(L) / L
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let trits = self.decomposer.build(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        {
            let trit = trits[0];
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, trit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..161 {
            let trit = trits[i];
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.build_input_selector(view, hash, trit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 2],
    ) -> Result<[CellOrUnconstrained; 1]> {
        let [key, _] = inputs;
        let trits = self.decomposer.witness(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash;
        {
            let mut view = view.sub(self.decomposer.height(), 0, stage_width);
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &trits, 0)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_ir, inputs)?;
        }
        for i in 1..161 {
            let mut view = view.sub(
                self.decomposer.height() + stage_height * (i / L),
                stage_width * (i % L),
                stage_width,
            );
            let inputs = std::array::from_fn(|i| view.cell(Self::SELECTOR_HEIGHT - 1, i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, stage_width, |view| {
                    self.witness_input_selector(view, &trits, i)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher_er[i - 1], inputs)?;
        }
        Ok([hash])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use primitive_types::{H256, U256};
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_ff::Field256;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon as poseidon1;
    use std::fmt::Debug;
    use std::sync::{Arc, LazyLock};

    const BLOWUP_LOG2: usize = 3;

    fn parse_hash(s: &'static str) -> H256 {
        s.parse().unwrap()
    }

    fn test_binary_smt<const H: usize, const L: usize>(
        key: u64,
        value: u64,
        path: [[Scalar; 2]; H],
        expected_root_hash: Scalar,
        circuit_commitment: H256,
    ) -> Result<()> {
        let key = Scalar::from(key);
        let value = Scalar::from(value);
        let chip = BinaryChip::<H, L>::new(path);
        assert_eq!(chip.width(), L * 6);
        assert_eq!(chip.height(), 1 + 196 * H.next_multiple_of(L) / L);
        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_rows([root_hash.unwrap().row()]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        let mut witness = circuit.make_witness();
        let inputs = [witness.cell(0, 0), witness.cell(0, 1)];
        witness.set(inputs[0], key);
        witness.set(inputs[1], value);
        let [root_hash] = witness.sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))?;
        let root_hash = match root_hash {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        };
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        let circuit = circuit.to_compressed::<Sha2Hash<Scalar>>(options);
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof)?;
        assert_eq!(openings[&root_hash], expected_root_hash);
        Ok(())
    }

    #[test]
    fn test_binary_smt_height_one_1() {
        let path = [[from_const(12), from_const(34)]];
        let root_hash =
            parse_scalar("0x45470d74563e5e49fe3bd2a161b36116e3c6a6a2f9c105bfe8c2599ff6116b06");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1, 1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1, 1>(1, 34, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_2() {
        let path = [[from_const(34), from_const(12)]];
        let root_hash =
            parse_scalar("0x6a6ca65c7ab651a6e7751e7a23df1d7ff66f745f1b09f4b39df2dfeb4e137422");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1, 1>(0, 34, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1, 1>(1, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78)]];
        let root_hash =
            parse_scalar("0x1ba4c686a3529d3bfc13890b2e1438b7adf780e2978cb2cabdd47653f402e8fe");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1, 1>(0, 56, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1, 1>(1, 78, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_one_lane_1() {
        let path = [
            [from_const(12), from_const(34)],
            [
                parse_scalar("0x45470d74563e5e49fe3bd2a161b36116e3c6a6a2f9c105bfe8c2599ff6116b06"),
                parse_scalar("0x1ba4c686a3529d3bfc13890b2e1438b7adf780e2978cb2cabdd47653f402e8fe"),
            ],
        ];
        let root_hash =
            parse_scalar("0x3f16169d0163139187336364cda1cac7f97b31dfbdabc4acba221d41792de5de");
        let c = parse_hash("0x4624a1fc0141a8d753764723b67f749af762216c30f212b4111c6efee396f361");
        assert!(test_binary_smt::<2, 1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<2, 1>(1, 34, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_one_lane_2() {
        let path = [
            [from_const(56), from_const(78)],
            [
                parse_scalar("0x45470d74563e5e49fe3bd2a161b36116e3c6a6a2f9c105bfe8c2599ff6116b06"),
                parse_scalar("0x1ba4c686a3529d3bfc13890b2e1438b7adf780e2978cb2cabdd47653f402e8fe"),
            ],
        ];
        let root_hash =
            parse_scalar("0x3f16169d0163139187336364cda1cac7f97b31dfbdabc4acba221d41792de5de");
        let c = parse_hash("0x4624a1fc0141a8d753764723b67f749af762216c30f212b4111c6efee396f361");
        assert!(test_binary_smt::<2, 1>(2, 56, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<2, 1>(3, 78, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_two_lanes_1() {
        let path = [
            [from_const(12), from_const(34)],
            [
                parse_scalar("0x45470d74563e5e49fe3bd2a161b36116e3c6a6a2f9c105bfe8c2599ff6116b06"),
                parse_scalar("0x1ba4c686a3529d3bfc13890b2e1438b7adf780e2978cb2cabdd47653f402e8fe"),
            ],
        ];
        let root_hash =
            parse_scalar("0x3f16169d0163139187336364cda1cac7f97b31dfbdabc4acba221d41792de5de");
        let c = parse_hash("0x98ebdec6c8e987a0e975b9dfcc881a41e6c31af1a96cc34c1b2df7942c24b331");
        assert!(test_binary_smt::<2, 2>(0, 12, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<2, 2>(1, 34, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_two_lanes_2() {
        let path = [
            [from_const(56), from_const(78)],
            [
                parse_scalar("0x45470d74563e5e49fe3bd2a161b36116e3c6a6a2f9c105bfe8c2599ff6116b06"),
                parse_scalar("0x1ba4c686a3529d3bfc13890b2e1438b7adf780e2978cb2cabdd47653f402e8fe"),
            ],
        ];
        let root_hash =
            parse_scalar("0x3f16169d0163139187336364cda1cac7f97b31dfbdabc4acba221d41792de5de");
        let c = parse_hash("0x98ebdec6c8e987a0e975b9dfcc881a41e6c31af1a96cc34c1b2df7942c24b331");
        assert!(test_binary_smt::<2, 2>(2, 56, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<2, 2>(3, 78, path, root_hash, c).is_ok());
    }

    fn test_ternary_smt<const H: usize, const L: usize>(
        key: u64,
        value: u64,
        path: [[Scalar; 3]; H],
        expected_root_hash: Scalar,
        circuit_commitment: H256,
    ) -> Result<()> {
        let key = Scalar::from(key);
        let value = Scalar::from(value);
        let chip = TernaryChip::<H, L>::new(path);
        assert_eq!(chip.width(), L * 8);
        assert_eq!(chip.height(), 1 + 196 * H.next_multiple_of(L) / L);
        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_rows([root_hash.unwrap().row()]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        let mut witness = circuit.make_witness();
        let inputs = [witness.cell(0, 0), witness.cell(0, 1)];
        witness.set(inputs[0], key);
        witness.set(inputs[1], value);
        let [root_hash] = witness.sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))?;
        let root_hash = match root_hash {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        };
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        let circuit = circuit.to_compressed::<Sha2Hash<Scalar>>(options);
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof)?;
        assert_eq!(openings[&root_hash], expected_root_hash);
        Ok(())
    }

    #[test]
    fn test_ternary_smt_height_one_1() {
        let path = [[from_const(12), from_const(34), from_const(56)]];
        let root_hash =
            parse_scalar("0x1125d1d7bcc64d065695f306f08db087abc90d214fd982461296e607de7d4d49");
        let c = parse_hash("0x512035d2b2db72ca18822ccd4cb4ce7d8455ffa95199f769260e89d31f4ca9f9");
        assert!(test_ternary_smt::<1, 1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1, 1>(1, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1, 1>(2, 56, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_2() {
        let path = [[from_const(34), from_const(56), from_const(12)]];
        let root_hash =
            parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3");
        let c = parse_hash("0x512035d2b2db72ca18822ccd4cb4ce7d8455ffa95199f769260e89d31f4ca9f9");
        assert!(test_ternary_smt::<1, 1>(0, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1, 1>(1, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1, 1>(2, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78), from_const(90)]];
        let root_hash =
            parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd");
        let c = parse_hash("0x512035d2b2db72ca18822ccd4cb4ce7d8455ffa95199f769260e89d31f4ca9f9");
        assert!(test_ternary_smt::<1, 1>(0, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1, 1>(1, 78, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1, 1>(2, 90, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_one_lane_1() {
        let path = [
            [from_const(12), from_const(34), from_const(56)],
            [
                parse_scalar("0x1125d1d7bcc64d065695f306f08db087abc90d214fd982461296e607de7d4d49"),
                parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3"),
                parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd"),
            ],
        ];
        let root_hash =
            parse_scalar("0x1bc60b83e94bbd9609c01954b66049bd4ba987570f4d5a68d89af45970a3930c");
        let c = parse_hash("0x7b12fdfa0fca6947c5e8c0a6273af575f003787182b65a04d47a972c168aa7bd");
        assert!(test_ternary_smt::<2, 1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 1>(1, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 1>(2, 56, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_one_lane_2() {
        let path = [
            [from_const(34), from_const(56), from_const(12)],
            [
                parse_scalar("0x1125d1d7bcc64d065695f306f08db087abc90d214fd982461296e607de7d4d49"),
                parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3"),
                parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd"),
            ],
        ];
        let root_hash =
            parse_scalar("0x1bc60b83e94bbd9609c01954b66049bd4ba987570f4d5a68d89af45970a3930c");
        let c = parse_hash("0x7b12fdfa0fca6947c5e8c0a6273af575f003787182b65a04d47a972c168aa7bd");
        assert!(test_ternary_smt::<2, 1>(3, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 1>(4, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 1>(5, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_one_lane_3() {
        let path = [
            [from_const(56), from_const(78), from_const(90)],
            [
                parse_scalar("0x1125d1d7bcc64d065695f306f08db087abc90d214fd982461296e607de7d4d49"),
                parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3"),
                parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd"),
            ],
        ];
        let root_hash =
            parse_scalar("0x1bc60b83e94bbd9609c01954b66049bd4ba987570f4d5a68d89af45970a3930c");
        let c = parse_hash("0x7b12fdfa0fca6947c5e8c0a6273af575f003787182b65a04d47a972c168aa7bd");
        assert!(test_ternary_smt::<2, 1>(6, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 1>(7, 78, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 1>(8, 90, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_two_lanes_1() {
        let path = [
            [from_const(12), from_const(34), from_const(56)],
            [
                parse_scalar("0x1125d1d7bcc64d065695f306f08db087abc90d214fd982461296e607de7d4d49"),
                parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3"),
                parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd"),
            ],
        ];
        let root_hash =
            parse_scalar("0x1bc60b83e94bbd9609c01954b66049bd4ba987570f4d5a68d89af45970a3930c");
        let c = parse_hash("0xabeda530f4fc805498b14b12f71239e35e6f4d0f1749aac803db2e76dfad7e4a");
        assert!(test_ternary_smt::<2, 2>(0, 12, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 2>(1, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 2>(2, 56, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_two_lanes_2() {
        let path = [
            [from_const(34), from_const(56), from_const(12)],
            [
                parse_scalar("0x1125d1d7bcc64d065695f306f08db087abc90d214fd982461296e607de7d4d49"),
                parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3"),
                parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd"),
            ],
        ];
        let root_hash =
            parse_scalar("0x1bc60b83e94bbd9609c01954b66049bd4ba987570f4d5a68d89af45970a3930c");
        let c = parse_hash("0xabeda530f4fc805498b14b12f71239e35e6f4d0f1749aac803db2e76dfad7e4a");
        assert!(test_ternary_smt::<2, 2>(3, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 2>(4, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 2>(5, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_two_lanes_3() {
        let path = [
            [from_const(56), from_const(78), from_const(90)],
            [
                parse_scalar("0x1125d1d7bcc64d065695f306f08db087abc90d214fd982461296e607de7d4d49"),
                parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3"),
                parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd"),
            ],
        ];
        let root_hash =
            parse_scalar("0x1bc60b83e94bbd9609c01954b66049bd4ba987570f4d5a68d89af45970a3930c");
        let c = parse_hash("0xabeda530f4fc805498b14b12f71239e35e6f4d0f1749aac803db2e76dfad7e4a");
        assert!(test_ternary_smt::<2, 2>(6, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 2>(7, 78, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2, 2>(8, 90, path, root_hash, c).is_ok());
    }

    trait Node: 'static + Debug + Send + Sync {
        fn hash(&self) -> Scalar;

        fn get_impl(&self, key: &U256) -> Scalar;

        fn get(&self, key: Scalar) -> Scalar {
            self.get_impl(&key.to_u256())
        }

        fn get_merkle_path_impl(&self, key: &U256) -> Vec<Vec<Scalar>>;

        fn get_merkle_path(&self, key: Scalar) -> Vec<Vec<Scalar>> {
            self.get_merkle_path_impl(&key.to_u256())
        }

        fn put_impl(self: Arc<Self>, key: &U256, value: Scalar) -> Arc<dyn Node>;

        fn put(self: Arc<Self>, key: Scalar, value: Scalar) -> Arc<dyn Node> {
            self.put_impl(&key.to_u256(), value)
        }
    }

    #[derive(Debug, Default, Copy, Clone)]
    struct Leaf(Scalar);

    impl Node for Leaf {
        fn hash(&self) -> Scalar {
            self.0
        }

        fn get_impl(&self, _key: &U256) -> Scalar {
            self.0
        }

        fn get_merkle_path_impl(&self, _key: &U256) -> Vec<Vec<Scalar>> {
            vec![]
        }

        fn put_impl(self: Arc<Self>, _key: &U256, value: Scalar) -> Arc<dyn Node> {
            Arc::new(Leaf(value))
        }
    }

    #[derive(Debug)]
    struct BinaryNode {
        level: usize,
        hash: Scalar,
        left: Arc<dyn Node>,
        right: Arc<dyn Node>,
    }

    impl BinaryNode {
        fn new(level: usize, left: Arc<dyn Node>, right: Arc<dyn Node>) -> Arc<dyn Node> {
            let hash = poseidon1::hash0::<poseidon1::BlueSkyConfig3, Scalar, 3, 2, 1>([
                left.hash(),
                right.hash(),
            ]);
            Arc::new(BinaryNode {
                level,
                hash,
                left,
                right,
            })
        }

        fn bit_at(&self, key: &U256) -> bool {
            (key >> self.level) & U256::one() != U256::zero()
        }
    }

    impl Node for BinaryNode {
        fn hash(&self) -> Scalar {
            self.hash
        }

        fn get_impl(&self, key: &U256) -> Scalar {
            if self.bit_at(key) {
                self.right.get_impl(key)
            } else {
                self.left.get_impl(key)
            }
        }

        fn get_merkle_path_impl(&self, key: &U256) -> Vec<Vec<Scalar>> {
            let mut path = if self.bit_at(key) {
                self.right.get_merkle_path_impl(key)
            } else {
                self.left.get_merkle_path_impl(key)
            };
            path.push(vec![self.left.hash(), self.right.hash()]);
            path
        }

        fn put_impl(self: Arc<Self>, key: &U256, value: Scalar) -> Arc<dyn Node> {
            if self.bit_at(key) {
                Self::new(
                    self.level,
                    self.left.clone(),
                    self.right.clone().put_impl(key, value),
                )
            } else {
                Self::new(
                    self.level,
                    self.left.clone().put_impl(key, value),
                    self.right.clone(),
                )
            }
        }
    }

    #[derive(Debug)]
    struct TernaryNode {
        level: usize,
        hash: Scalar,
        children: [Arc<dyn Node>; 3],
    }

    impl TernaryNode {
        fn new(level: usize, children: [Arc<dyn Node>; 3]) -> Arc<Self> {
            let hash = poseidon1::hash0::<poseidon1::BlueSkyConfig4, Scalar, 4, 3, 1>([
                children[0].hash(),
                children[1].hash(),
                children[2].hash(),
            ]);
            Arc::new(TernaryNode {
                level,
                hash,
                children,
            })
        }

        fn trit_at(&self, key: &U256) -> usize {
            let divisor = U256::from(3).pow(self.level.into());
            ((key / divisor) % 3).try_into().unwrap()
        }
    }

    impl Node for TernaryNode {
        fn hash(&self) -> Scalar {
            self.hash
        }

        fn get_impl(&self, key: &U256) -> Scalar {
            self.children[self.trit_at(key)].get_impl(key)
        }

        fn get_merkle_path_impl(&self, key: &U256) -> Vec<Vec<Scalar>> {
            let mut path = self.children[self.trit_at(key)].get_merkle_path_impl(key);
            path.push(self.children.iter().map(|child| child.hash()).collect());
            path
        }

        fn put_impl(self: Arc<Self>, key: &U256, value: Scalar) -> Arc<dyn Node> {
            let trit = self.trit_at(key);
            let child = self.children[trit].clone().put_impl(key, value);
            match trit {
                0 => Self::new(
                    self.level,
                    [child, self.children[1].clone(), self.children[2].clone()],
                ),
                1 => Self::new(
                    self.level,
                    [self.children[0].clone(), child, self.children[2].clone()],
                ),
                2 => Self::new(
                    self.level,
                    [self.children[0].clone(), self.children[1].clone(), child],
                ),
                _ => panic!(),
            }
        }
    }

    fn get_empty_binary_tree() -> Arc<dyn Node> {
        static TREE: LazyLock<Arc<dyn Node>> = LazyLock::new(|| {
            let mut node: Arc<dyn Node> = Arc::new(Leaf::default());
            for i in 0..256 {
                node = BinaryNode::new(i, node.clone(), node.clone());
            }
            node
        });
        TREE.clone()
    }

    fn get_empty_ternary_tree() -> Arc<dyn Node> {
        static TREE: LazyLock<Arc<dyn Node>> = LazyLock::new(|| {
            let mut node: Arc<dyn Node> = Arc::new(Leaf::default());
            for i in 0..161 {
                node = TernaryNode::new(i, [node.clone(), node.clone(), node.clone()]);
            }
            node
        });
        TREE.clone()
    }

    fn test_full_binary_smt_impl<I: IntoIterator<Item = (u64, u64)>>(
        entries: I,
        key: u64,
    ) -> Result<()> {
        const LANES: usize = 52;

        let tree = {
            let mut tree = get_empty_binary_tree();
            for (key, value) in entries {
                tree = tree.put(key.into(), value.into());
            }
            tree
        };
        let key = key.into();
        let value = tree.get(key);
        let path: [[Scalar; 2]; 256] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = tree.hash();

        let chip = FullBinaryChip::<LANES>::new(path);
        assert_eq!(chip.stage_width(), 6);
        assert_eq!(chip.stage_height(), 196);
        assert_eq!(chip.width(), std::cmp::max(257, 6 * LANES));
        assert_eq!(
            chip.height(),
            3 + chip.stage_height() * 256usize.next_multiple_of(LANES) / LANES
        );

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_rows([root_hash.unwrap().row()]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), chip.height() + 1);
        assert_eq!(circuit.num_columns(), chip.width());

        let mut witness = circuit.make_witness();
        let inputs = [witness.cell(0, 0), witness.cell(0, 1)];
        witness.set(inputs[0], key);
        witness.set(inputs[1], value);
        let [root_hash] = witness.sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))?;
        let root_hash = match root_hash {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        };

        circuit.check_witness(&witness).unwrap();

        let options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        let openings = circuit.verify(&proof, options)?;
        assert_eq!(openings[&root_hash], expected_root_hash);

        Ok(())
    }

    #[test]
    fn test_full_binary_smt_empty() {
        assert!(test_full_binary_smt_impl([], 0).is_ok());
        assert!(test_full_binary_smt_impl([], 1).is_ok());
        assert!(test_full_binary_smt_impl([], 2).is_ok());
        assert!(test_full_binary_smt_impl([], 3).is_ok());
        assert!(test_full_binary_smt_impl([], 4).is_ok());
        assert!(test_full_binary_smt_impl([], 5).is_ok());
    }

    #[test]
    fn test_full_binary_smt_one_entry() {
        let entries = [(12, 34)];
        assert!(test_full_binary_smt_impl(entries, 0).is_ok());
        assert!(test_full_binary_smt_impl(entries, 1).is_ok());
        assert!(test_full_binary_smt_impl(entries, 2).is_ok());
        assert!(test_full_binary_smt_impl(entries, 11).is_ok());
        assert!(test_full_binary_smt_impl(entries, 12).is_ok());
        assert!(test_full_binary_smt_impl(entries, 13).is_ok());
    }

    #[test]
    fn test_full_binary_smt_two_entries() {
        let entries = [(34, 56), (78, 12)];
        assert!(test_full_binary_smt_impl(entries, 0).is_ok());
        assert!(test_full_binary_smt_impl(entries, 1).is_ok());
        assert!(test_full_binary_smt_impl(entries, 2).is_ok());
        assert!(test_full_binary_smt_impl(entries, 33).is_ok());
        assert!(test_full_binary_smt_impl(entries, 34).is_ok());
        assert!(test_full_binary_smt_impl(entries, 35).is_ok());
        assert!(test_full_binary_smt_impl(entries, 77).is_ok());
        assert!(test_full_binary_smt_impl(entries, 78).is_ok());
        assert!(test_full_binary_smt_impl(entries, 79).is_ok());
    }

    fn test_full_ternary_smt_impl<I: IntoIterator<Item = (u64, u64)>>(
        entries: I,
        key: u64,
    ) -> Result<()> {
        const LANES: usize = 33;

        let tree = {
            let mut tree = get_empty_ternary_tree();
            for (key, value) in entries {
                tree = tree.put(key.into(), value.into());
            }
            tree
        };
        let key = key.into();
        let value = tree.get(key);
        let path: [[Scalar; 3]; 161] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = tree.hash();

        let chip = FullTernaryChip::<LANES>::new(path);
        assert_eq!(chip.stage_width(), 8);
        assert_eq!(chip.stage_height(), 196);
        assert_eq!(chip.width(), std::cmp::max(162, 8 * LANES));
        assert_eq!(
            chip.height(),
            4 + chip.stage_height() * 161usize.next_multiple_of(LANES) / LANES
        );

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_rows([root_hash.unwrap().row()]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), chip.height() + 1);
        assert_eq!(circuit.num_columns(), chip.width());

        let mut witness = circuit.make_witness();
        let inputs = [witness.cell(0, 0), witness.cell(0, 1)];
        witness.set(inputs[0], key);
        witness.set(inputs[1], value);
        let [root_hash] = witness.sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))?;
        let root_hash = match root_hash {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        };

        circuit.check_witness(&witness).unwrap();

        let options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        let openings = circuit.verify(&proof, options)?;
        assert_eq!(openings[&root_hash], expected_root_hash);

        Ok(())
    }

    #[test]
    fn test_full_ternary_smt_empty() {
        assert!(test_full_ternary_smt_impl([], 0).is_ok());
        assert!(test_full_ternary_smt_impl([], 1).is_ok());
        assert!(test_full_ternary_smt_impl([], 2).is_ok());
        assert!(test_full_ternary_smt_impl([], 3).is_ok());
        assert!(test_full_ternary_smt_impl([], 4).is_ok());
        assert!(test_full_ternary_smt_impl([], 5).is_ok());
    }

    #[test]
    fn test_full_ternary_smt_one_entry() {
        let entries = [(12, 34)];
        assert!(test_full_ternary_smt_impl(entries, 0).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 1).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 2).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 11).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 12).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 13).is_ok());
    }

    #[test]
    fn test_full_ternary_smt_two_entries() {
        let entries = [(34, 56), (78, 12)];
        assert!(test_full_ternary_smt_impl(entries, 0).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 1).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 2).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 33).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 34).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 35).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 77).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 78).is_ok());
        assert!(test_full_ternary_smt_impl(entries, 79).is_ok());
    }
}
