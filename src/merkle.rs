use crate::poseidon1;
use crate::xits;
use anyhow::Result;
use starkom_ff::{Field256, PrimeField256};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, make_const, rvar, var,
};
use starkom_poseidon::Config as PoseidonConfig;

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height `H`.
///
/// WARNING: `H` must be strictly less than 255. Do NOT use this chip if `H` spans the full BlueSky
/// range, as in that case the bit decomposition of the key would be UNSAFE! Use the
/// [`FullBinaryChip`] below instead.
///
/// The generic argument `L` is the number of lanes (parallel hash stages) used by the chip.
#[derive(Debug, Clone)]
pub struct BinaryChip<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 3>> {
    hasher: poseidon1::PermutationChip<F, C, 3>,
    path: [[F; 2]; H],
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 3>> Default for BinaryChip<F, H, C> {
    fn default() -> Self {
        Self::new([[F::ZERO; 2]; H])
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 3>> BinaryChip<F, H, C> {
    const SELECTOR_WIDTH: usize = 7;

    pub fn new(path: [[F; 2]; H]) -> Self {
        assert!(H < F::NUM_BITS);
        Self {
            hasher: poseidon1::PermutationChip::default(),
            path,
        }
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+----+
    /// | H1 | H2 | B  | S  | I1 | I2 | 0  |
    /// +----+----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = peer hash (unconstrained)
    /// B = key bit
    /// S = key bit sum
    /// I1 = left-hand-side input hash (either H1 or H2)
    /// I2 = right-hand-side input hash (either H1 or H2)
    /// 0 = a zero scalar used as input capacity
    ///
    /// The B column holds the decomposed bits of the key and the S column is used to reconstruct
    /// the original key by summing the decomposed bits weighted by the corresponding powers of two.
    /// Once reconstructed, the sum must be constrained to equal the original key.
    ///
    /// Note that the last three elements of the selector are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: Option<Cell>,
        i: usize,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        if i > 0 {
            view.add_gate(
                0,
                var(2) * (make_const(F::from(2u8)) ^ (i as isize)) + rvar(3, -1) - var(3),
            );
        } else {
            view.add_gate(0, var(2) - var(3));
        }
        view.add_gate(
            0,
            var(2) * var(1) + (make_const(F::ONE) - var(2)) * var(0) - var(4),
        );
        view.add_gate(
            0,
            var(2) * var(0) + (make_const(F::ONE) - var(2)) * var(1) - var(5),
        );
        view.add_gate(0, var(6));
    }

    /// See [`Self::build_input_selector`] for the layout.
    fn witness_input_selector(&self, view: &mut impl WitnessView<F>, bits: &[F], i: usize) {
        let bit = bits[i];
        if bit != F::ZERO {
            view.set(view.cell(0, 0), self.path[i][1]);
            view.set(view.cell(0, 1), self.path[i][0]);
        } else {
            view.set(view.cell(0, 0), self.path[i][0]);
            view.set(view.cell(0, 1), self.path[i][1]);
        }
        view.set(view.cell(0, 2), bit);
        if i > 0 {
            view.set(
                view.cell(0, 3),
                bit * F::from(2u8).pow_small(i) + view.get_at(view.cell(-1, 3)),
            );
        } else {
            view.copy(view.cell(0, 2).into(), view.cell(0, 3));
        }
        view.set(view.cell(0, 4), self.path[i][0]);
        view.set(view.cell(0, 5), self.path[i][1]);
        view.set(view.cell(0, 6), F::ZERO);
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 3>> PlonkChip<F, 2, 1>
    for BinaryChip<F, H, C>
{
    fn width(&self) -> usize {
        Self::SELECTOR_WIDTH + self.hasher.width()
    }

    fn height(&self) -> usize {
        H
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let width = self.width();
        let mut hash = value;
        for i in 0..H {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 4 + i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.build_input_selector(view, hash, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
        }
        view.connect(key, view.cell(H - 1, 3).into());
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 2],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        let [key, value] = inputs;
        let bits = xits::decompose_scalar_bits::<F, H>(view.get(key));
        let width = self.width();
        let mut hash = value;
        for i in 0..H {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 4 + i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
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
pub struct TernaryChip<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> {
    hasher: poseidon1::PermutationChip<F, C, 4>,
    path: [[F; 3]; H],
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> Default for TernaryChip<F, H, C> {
    fn default() -> Self {
        Self::new([[F::ZERO; 3]; H])
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> TernaryChip<F, H, C> {
    const SELECTOR_WIDTH: usize = 9;

    pub fn new(path: [[F; 3]; H]) -> Self {
        // TODO: assert that H is strictly less than the number of trits required to represent a
        // scalar.
        Self {
            hasher: poseidon1::PermutationChip::default(),
            path,
        }
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+----+----+----+
    /// | H1 | H2 | H3 | T  | S  | I1 | I2 | I3 | 0  |
    /// +----+----+----+----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = first peer hash (unconstrained)
    /// H3 = second peer hash (unconstrained)
    /// T = key trit
    /// S = key trit sum
    /// I1 = first input hash (one of H1, H2, or H3)
    /// I2 = second input hash (one of H1, H2, or H3)
    /// I3 = third input hash (one of H1, H2, or H3)
    /// 0 = a zero scalar used as input capacity
    ///
    /// The T column holds the decomposed trits of the key and the S column is used to reconstruct
    /// the original key by summing the decomposed trits weighted by the corresponding powers of
    /// three. Once reconstructed, the sum must be constrained to equal the original key.
    ///
    /// Note that the last four elements are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: Option<Cell>,
        i: usize,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        if i > 0 {
            view.add_gate(
                0,
                var(3) * (make_const(F::from(3u8)) ^ (i as isize)) + rvar(4, -1) - var(4),
            );
        } else {
            view.add_gate(0, var(3) - var(4));
        }
        let l0 = ((var(3) ^ 2) - var(3) * 3 + 2) / 2;
        let l1 = var(3) * 2 - (var(3) ^ 2);
        let l2 = ((var(3) ^ 2) - var(3)) / 2;
        view.add_gate(
            0,
            l0.clone() * var(0) + l1.clone() * var(1) + l2.clone() * var(1) - var(5),
        );
        view.add_gate(
            0,
            l0.clone() * var(1) + l1.clone() * var(0) + l2.clone() * var(2) - var(6),
        );
        view.add_gate(0, l0 * var(2) + l1 * var(2) + l2 * var(0) - var(7));
        view.add_gate(0, var(8));
    }

    /// See [`Self::build_input_selector`] for the layout.
    fn witness_input_selector(&self, view: &mut impl WitnessView<F>, trits: &[F], i: usize) {
        let trit = trits[i];
        if trit == F::from(0u8) {
            view.set(view.cell(0, 0), self.path[i][0]);
            view.set(view.cell(0, 1), self.path[i][1]);
            view.set(view.cell(0, 2), self.path[i][2]);
        } else if trit == F::from(1u8) {
            view.set(view.cell(0, 0), self.path[i][1]);
            view.set(view.cell(0, 1), self.path[i][0]);
            view.set(view.cell(0, 2), self.path[i][2]);
        } else if trit == F::from(2u8) {
            view.set(view.cell(0, 0), self.path[i][2]);
            view.set(view.cell(0, 1), self.path[i][0]);
            view.set(view.cell(0, 2), self.path[i][1]);
        } else {
            panic!("invalid trit value {}", trit);
        }
        view.set(view.cell(0, 3).into(), trit);
        if i > 0 {
            view.set(
                view.cell(0, 4),
                trit * F::from(3u8).pow_small(i) + view.get_at(view.cell(-1, 4)),
            );
        } else {
            view.copy(view.cell(0, 3).into(), view.cell(0, 4));
        }
        view.set(view.cell(0, 5), self.path[i][0]);
        view.set(view.cell(0, 6), self.path[i][1]);
        view.set(view.cell(0, 7), self.path[i][2]);
        view.set(view.cell(0, 8), F::ZERO);
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> PlonkChip<F, 2, 1>
    for TernaryChip<F, H, C>
{
    fn width(&self) -> usize {
        Self::SELECTOR_WIDTH + self.hasher.width()
    }

    fn height(&self) -> usize {
        H
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let width = self.width();
        let mut hash = value;
        for i in 0..H {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 5 + i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.build_input_selector(view, hash, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
        }
        view.connect(key, view.cell(H - 1, 4).into());
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 2],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        let [key, value] = inputs;
        let trits = xits::decompose_scalar_trits::<F, H>(view.get(key));
        let width = self.width();
        let mut hash = value;
        for i in 0..H {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 5 + i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.witness_input_selector(view, &trits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height 256.
///
/// The keys of such a tree span the full BlueSky range. Internally this chip uses a
/// [`xits::FullBitDecomposerChip256`], making the 256-bit decomposition safe at the cost of some
/// extra constraints.
///
/// If you don't need 256- or 255-bit keys use [`BinaryChip`].
///
/// The generic argument `L` is the number of lanes (parallel hash stages) used by the chip.
#[derive(Debug, Clone)]
pub struct FullBinaryChip<F: PrimeField256, C: PoseidonConfig<F, 3>> {
    decomposer: xits::FullBitDecomposerChip256<F>,
    hasher: poseidon1::PermutationChip<F, C, 3>,
    path: [[F; 2]; 256],
}

impl<F: PrimeField256, C: PoseidonConfig<F, 3>> Default for FullBinaryChip<F, C> {
    fn default() -> Self {
        Self::new([[F::ZERO; 2]; 256])
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 3>> FullBinaryChip<F, C> {
    const SELECTOR_WIDTH: usize = 6;

    pub fn new(path: [[F; 2]; 256]) -> Self {
        Self {
            decomposer: xits::FullBitDecomposerChip256::default(),
            hasher: poseidon1::PermutationChip::default(),
            path,
        }
    }

    fn stage_width(&self) -> usize {
        Self::SELECTOR_WIDTH + self.hasher.width()
    }

    fn stage_height(&self) -> usize {
        self.hasher.height()
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+
    /// | H1 | H2 | B  | I1 | I2 | 0  |
    /// +----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = peer hash
    /// B = key bit
    /// I1 = left-hand-side input hash (either H1 or H2)
    /// I2 = right-hand-side input hash (either H1 or H2)
    /// 0 = a zero scalar used as input capacity
    ///
    /// Note that the last three elements are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: Option<Cell>,
        bit: Option<Cell>,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.connect(bit, view.cell(0, 2).into());
        view.add_gate(
            0,
            var(2) * var(1) + (make_const(F::ONE) - var(2)) * var(0) - var(3),
        );
        view.add_gate(
            0,
            var(2) * var(0) + (make_const(F::ONE) - var(2)) * var(1) - var(4),
        );
        view.add_gate(0, var(5));
    }

    /// See [`Self::build_input_selector`] for the layout.
    fn witness_input_selector(
        &self,
        view: &mut impl WitnessView<F>,
        bits: &[CellOrUnconstrained<F>],
        i: usize,
    ) {
        let bit = bits[i];
        let bit_value = view.get(bit);
        if bit_value != F::ZERO {
            view.set(view.cell(0, 0), self.path[i][1]);
            view.set(view.cell(0, 1), self.path[i][0]);
        } else {
            view.set(view.cell(0, 0), self.path[i][0]);
            view.set(view.cell(0, 1), self.path[i][1]);
        }
        view.copy(bits[i], view.cell(0, 2));
        view.set(view.cell(0, 3), self.path[i][0]);
        view.set(view.cell(0, 4), self.path[i][1]);
        view.set(view.cell(0, 5), F::ZERO);
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 3>> PlonkChip<F, 2, 1> for FullBinaryChip<F, C> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.stage_width())
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.stage_height() * 256
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.build(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        for i in 0..256 {
            let bit = bits[i];
            let mut view = view.sub(
                self.decomposer.height() + stage_height * i,
                0,
                stage_width.into(),
                stage_height.into(),
            );
            let inputs = std::array::from_fn(|i| view.cell(0, 3 + i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), Some(1), |view| {
                    self.build_input_selector(view, hash, bit)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 2],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.witness(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        for i in 0..256 {
            let mut view = view.sub(
                self.decomposer.height() + stage_height * i,
                0,
                stage_width.into(),
                stage_height.into(),
            );
            let inputs = std::array::from_fn(|i| view.cell(0, 3 + i).into());
            [hash, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), Some(1), |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a ternary Sparse Merkle Tree of height 161.
///
/// The keys of such a tree span the full BlueSky range. Internally this chip uses a
/// [`xits::FullTritDecomposerChip256`], making the 161-trit decomposition safe at the cost of some
/// extra constraints.
///
/// If you don't need 161-trit keys use [`TernaryChip`].
#[derive(Debug, Clone)]
pub struct FullTernaryChip<F: PrimeField256, C: PoseidonConfig<F, 4>> {
    decomposer: xits::FullTritDecomposerChip256<F>,
    hasher: poseidon1::PermutationChip<F, C, 4>,
    path: [[F; 3]; 161],
}

impl<F: PrimeField256, C: PoseidonConfig<F, 4>> Default for FullTernaryChip<F, C> {
    fn default() -> Self {
        Self::new([[F::ZERO; 3]; 161])
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 4>> FullTernaryChip<F, C> {
    const SELECTOR_WIDTH: usize = 8;

    pub fn new(path: [[F; 3]; 161]) -> Self {
        Self {
            decomposer: xits::FullTritDecomposerChip256::default(),
            hasher: poseidon1::PermutationChip::default(),
            path,
        }
    }

    fn stage_width(&self) -> usize {
        Self::SELECTOR_WIDTH + self.hasher.width()
    }

    fn stage_height(&self) -> usize {
        self.hasher.height()
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+----+----+
    /// | H1 | H2 | H3 | T  | I1 | I2 | I3 | 0  |
    /// +----+----+----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = first peer hash
    /// H3 = second peer hash
    /// T = key trit
    /// I1 = first input hash (one of H1, H2, or H3)
    /// I2 = second input hash (one of H1, H2, or H3)
    /// I3 = third input hash (one of H1, H2, or H3)
    /// 0 = a zero scalar used as input capacity
    ///
    /// Note that the last four elements are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
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
            l0.clone() * var(0) + l1.clone() * var(1) + l2.clone() * var(1) - var(4),
        );
        view.add_gate(
            0,
            l0.clone() * var(1) + l1.clone() * var(0) + l2.clone() * var(2) - var(5),
        );
        view.add_gate(0, l0 * var(2) + l1 * var(2) + l2 * var(0) - var(6));
        view.add_gate(0, var(7));
    }

    /// See [`Self::build_input_selector`] for the layout.
    fn witness_input_selector(
        &self,
        view: &mut impl WitnessView<F>,
        trits: &[CellOrUnconstrained<F>],
        i: usize,
    ) {
        let trit = trits[i];
        let trit_value = view.get(trit);
        if trit_value == F::from(0u8) {
            view.set(view.cell(0, 0), self.path[i][0]);
            view.set(view.cell(0, 1), self.path[i][1]);
            view.set(view.cell(0, 2), self.path[i][2]);
        } else if trit_value == F::from(1u8) {
            view.set(view.cell(0, 0), self.path[i][1]);
            view.set(view.cell(0, 1), self.path[i][0]);
            view.set(view.cell(0, 2), self.path[i][2]);
        } else if trit_value == F::from(2u8) {
            view.set(view.cell(0, 0), self.path[i][2]);
            view.set(view.cell(0, 1), self.path[i][0]);
            view.set(view.cell(0, 2), self.path[i][1]);
        } else {
            panic!("invalid trit value {}", trit_value);
        }
        view.copy(trit, view.cell(0, 3).into());
        view.set(view.cell(0, 4), self.path[i][0]);
        view.set(view.cell(0, 5), self.path[i][1]);
        view.set(view.cell(0, 6), self.path[i][2]);
        view.set(view.cell(0, 7), F::ZERO);
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 4>> PlonkChip<F, 2, 1> for FullTernaryChip<F, C> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.stage_width())
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.stage_height() * 161
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let trits = self.decomposer.build(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        for i in 0..161 {
            let trit = trits[i];
            let mut view = view.sub(
                self.decomposer.height() + stage_height * i,
                0,
                stage_width.into(),
                stage_height.into(),
            );
            let inputs = std::array::from_fn(|i| view.cell(0, 4 + i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), Some(1), |view| {
                    self.build_input_selector(view, hash, trit)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 2],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        let [key, value] = inputs;
        let trits = self.decomposer.witness(view, [key])?;
        let stage_width = self.stage_width();
        let stage_height = self.stage_height();
        let mut hash = value;
        for i in 0..161 {
            let mut view = view.sub(
                self.decomposer.height() + stage_height * i,
                0,
                stage_width.into(),
                stage_height.into(),
            );
            let inputs = std::array::from_fn(|i| view.cell(0, 4 + i).into());
            [hash, _, _, _] = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), Some(1), |view| {
                    self.witness_input_selector(view, &trits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
        }
        Ok([hash])
    }
}

#[cfg(all(test, feature = "bluesky"))]
mod tests {
    use super::*;
    use primitive_types::{H256, U256};
    use starkom_bluesky::{Scalar, from_const, parse_scalar};
    use starkom_ff::Field;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon::{self as poseidon1, BlueSkyConfig3, BlueSkyConfig4};
    use std::collections::BTreeMap;
    use std::fmt::Debug;
    use std::sync::{Arc, LazyLock, Mutex};

    const BLOWUP_LOG2: usize = 3;

    fn parse_hash(s: &'static str) -> H256 {
        s.parse().unwrap()
    }

    fn test_binary_smt<const H: usize>(
        key: u64,
        value: u64,
        path: [[Scalar; 2]; H],
        expected_root_hash: Scalar,
        circuit_commitment: H256,
    ) -> Result<()> {
        let key = Scalar::from(key);
        let value = Scalar::from(value);
        let chip = BinaryChip::<Scalar, H, BlueSkyConfig3>::new(path);
        assert_eq!(chip.width(), 92);
        assert_eq!(chip.height(), H);
        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
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
        // assert_eq!(circuit.commitment(), circuit_commitment);  // TODO: re-enable
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
        assert!(test_binary_smt::<1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1>(1, 34, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_2() {
        let path = [[from_const(34), from_const(12)]];
        let root_hash =
            parse_scalar("0x6a6ca65c7ab651a6e7751e7a23df1d7ff66f745f1b09f4b39df2dfeb4e137422");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1>(0, 34, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1>(1, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78)]];
        let root_hash =
            parse_scalar("0x1ba4c686a3529d3bfc13890b2e1438b7adf780e2978cb2cabdd47653f402e8fe");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1>(0, 56, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1>(1, 78, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_1() {
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
        assert!(test_binary_smt::<2>(0, 12, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<2>(1, 34, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_2() {
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
        assert!(test_binary_smt::<2>(2, 56, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<2>(3, 78, path, root_hash, c).is_ok());
    }

    fn test_ternary_smt<const H: usize>(
        key: u64,
        value: u64,
        path: [[Scalar; 3]; H],
        expected_root_hash: Scalar,
        circuit_commitment: H256,
    ) -> Result<()> {
        let key = Scalar::from(key);
        let value = Scalar::from(value);
        let chip = TernaryChip::<Scalar, H, BlueSkyConfig4>::new(path);
        assert_eq!(chip.width(), 104);
        assert_eq!(chip.height(), H);
        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
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
        // assert_eq!(circuit.commitment(), circuit_commitment);  // TODO: re-enable
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
        assert!(test_ternary_smt::<1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(1, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(2, 56, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_2() {
        let path = [[from_const(34), from_const(56), from_const(12)]];
        let root_hash =
            parse_scalar("0x082b815a78ff9655cf614728ee7784b92be9d97086ccc0065b37cfa666efc2f3");
        let c = parse_hash("0x512035d2b2db72ca18822ccd4cb4ce7d8455ffa95199f769260e89d31f4ca9f9");
        assert!(test_ternary_smt::<1>(0, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(1, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(2, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78), from_const(90)]];
        let root_hash =
            parse_scalar("0x516b43041b6e111a7be5670972354589d8686593fbd2a994e14c53e55bb803cd");
        let c = parse_hash("0x512035d2b2db72ca18822ccd4cb4ce7d8455ffa95199f769260e89d31f4ca9f9");
        assert!(test_ternary_smt::<1>(0, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(1, 78, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(2, 90, path, root_hash, c).is_ok());
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
        assert!(test_ternary_smt::<2>(0, 12, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2>(1, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2>(2, 56, path, root_hash, c).is_ok());
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
        assert!(test_ternary_smt::<2>(3, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2>(4, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2>(5, 12, path, root_hash, c).is_ok());
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
        assert!(test_ternary_smt::<2>(6, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2>(7, 78, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<2>(8, 90, path, root_hash, c).is_ok());
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
            (key >> (self.level - 1)) & U256::one() != U256::zero()
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
        fn new(level: usize, children: [Arc<dyn Node>; 3]) -> Arc<dyn Node> {
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
            let divisor = U256::from(3).pow((self.level - 1).into());
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

    fn get_empty_binary_tree_locked(
        nodes_by_level: &mut BTreeMap<usize, Arc<dyn Node>>,
        level: usize,
    ) -> Arc<dyn Node> {
        match nodes_by_level.get_mut(&level) {
            Some(node) => node.clone(),
            None => {
                let node = if level > 0 {
                    let child = get_empty_binary_tree_locked(nodes_by_level, level - 1);
                    BinaryNode::new(level, child.clone(), child.clone())
                } else {
                    Arc::new(Leaf::default())
                };
                nodes_by_level.insert(level, node.clone());
                node
            }
        }
    }

    fn get_empty_binary_tree(level: usize) -> Arc<dyn Node> {
        static NODES_BY_LEVEL: LazyLock<Mutex<BTreeMap<usize, Arc<dyn Node>>>> =
            LazyLock::new(|| Mutex::new(BTreeMap::default()));
        let mut nodes_by_level = NODES_BY_LEVEL.lock().unwrap();
        get_empty_binary_tree_locked(&mut nodes_by_level, level)
    }

    fn get_empty_ternary_tree_locked(
        nodes_by_level: &mut BTreeMap<usize, Arc<dyn Node>>,
        level: usize,
    ) -> Arc<dyn Node> {
        match nodes_by_level.get_mut(&level) {
            Some(node) => node.clone(),
            None => {
                let node = if level > 0 {
                    let child = get_empty_ternary_tree_locked(nodes_by_level, level - 1);
                    TernaryNode::new(level, [child.clone(), child.clone(), child.clone()])
                } else {
                    Arc::new(Leaf::default())
                };
                nodes_by_level.insert(level, node.clone());
                node
            }
        }
    }

    fn get_empty_ternary_tree(level: usize) -> Arc<dyn Node> {
        static NODES_BY_LEVEL: LazyLock<Mutex<BTreeMap<usize, Arc<dyn Node>>>> =
            LazyLock::new(|| Mutex::new(BTreeMap::default()));
        let mut nodes_by_level = NODES_BY_LEVEL.lock().unwrap();
        get_empty_ternary_tree_locked(&mut nodes_by_level, level)
    }

    fn test_tall_binary_smt_impl<const H: usize>(
        entries: impl IntoIterator<Item = (u64, u64)>,
        key: u64,
    ) -> Result<()> {
        let tree = {
            let mut tree = get_empty_binary_tree(H);
            for (key, value) in entries {
                tree = tree.put(key.into(), value.into());
            }
            tree
        };
        let key = key.into();
        let value = tree.get(key);
        let path: [[Scalar; 2]; H] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = tree.hash();

        let chip = BinaryChip::<Scalar, H, BlueSkyConfig3>::new(path);
        assert_eq!(chip.width(), 92);
        assert_eq!(chip.height(), H);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
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
    fn test_tall_binary_smt_empty() {
        assert!(test_tall_binary_smt_impl::<20>([], 0).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 1).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 2).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 3).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 4).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 5).is_ok());
    }

    #[test]
    fn test_tall_binary_smt_one_entry() {
        let entries = [(12, 34)];
        assert!(test_tall_binary_smt_impl::<20>(entries, 0).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 1).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 2).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 11).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 12).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 13).is_ok());
    }

    #[test]
    fn test_tall_binary_smt_two_entries() {
        let entries = [(34, 56), (78, 12)];
        assert!(test_tall_binary_smt_impl::<20>(entries, 0).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 1).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 2).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 33).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 34).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 35).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 77).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 78).is_ok());
        assert!(test_tall_binary_smt_impl::<20>(entries, 79).is_ok());
    }

    #[test]
    #[ignore]
    fn test_taller_binary_smt_empty() {
        assert!(test_tall_binary_smt_impl::<160>([], 0).is_ok());
        assert!(test_tall_binary_smt_impl::<160>([], 1).is_ok());
        assert!(test_tall_binary_smt_impl::<160>([], 2).is_ok());
        assert!(test_tall_binary_smt_impl::<160>([], 3).is_ok());
        assert!(test_tall_binary_smt_impl::<160>([], 4).is_ok());
        assert!(test_tall_binary_smt_impl::<160>([], 5).is_ok());
    }

    #[test]
    #[ignore]
    fn test_taller_binary_smt_one_entry() {
        let entries = [(12, 34)];
        assert!(test_tall_binary_smt_impl::<160>(entries, 0).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 1).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 2).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 11).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 12).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 13).is_ok());
    }

    #[test]
    #[ignore]
    fn test_taller_binary_smt_two_entries() {
        let entries = [(34, 56), (78, 12)];
        assert!(test_tall_binary_smt_impl::<160>(entries, 0).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 1).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 2).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 33).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 34).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 35).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 77).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 78).is_ok());
        assert!(test_tall_binary_smt_impl::<160>(entries, 79).is_ok());
    }

    fn test_tall_ternary_smt_impl<const H: usize>(
        entries: impl IntoIterator<Item = (u64, u64)>,
        key: u64,
    ) -> Result<()> {
        let tree = {
            let mut tree = get_empty_ternary_tree(H);
            for (key, value) in entries {
                tree = tree.put(key.into(), value.into());
            }
            tree
        };
        let key = key.into();
        let value = tree.get(key);
        let path: [[Scalar; 3]; H] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = tree.hash();

        let chip = TernaryChip::<Scalar, H, BlueSkyConfig4>::new(path);
        assert_eq!(chip.width(), 104);
        assert_eq!(chip.height(), H);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
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
    fn test_tall_ternary_smt_empty() {
        assert!(test_tall_ternary_smt_impl::<13>([], 0).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>([], 1).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>([], 2).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>([], 3).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>([], 4).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>([], 5).is_ok());
    }

    #[test]
    fn test_tall_ternary_smt_one_entry() {
        let entries = [(12, 34)];
        assert!(test_tall_ternary_smt_impl::<13>(entries, 0).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 1).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 2).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 11).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 12).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 13).is_ok());
    }

    #[test]
    fn test_tall_ternary_smt_two_entries() {
        let entries = [(34, 56), (78, 12)];
        assert!(test_tall_ternary_smt_impl::<13>(entries, 0).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 1).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 2).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 33).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 34).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 35).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 77).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 78).is_ok());
        assert!(test_tall_ternary_smt_impl::<13>(entries, 79).is_ok());
    }

    #[test]
    #[ignore]
    fn test_taller_ternary_smt_empty() {
        assert!(test_tall_ternary_smt_impl::<101>([], 0).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>([], 1).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>([], 2).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>([], 3).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>([], 4).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>([], 5).is_ok());
    }

    #[test]
    #[ignore]
    fn test_taller_ternary_smt_one_entry() {
        let entries = [(12, 34)];
        assert!(test_tall_ternary_smt_impl::<101>(entries, 0).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 1).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 2).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 11).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 12).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 13).is_ok());
    }

    #[test]
    #[ignore]
    fn test_taller_ternary_smt_two_entries() {
        let entries = [(34, 56), (78, 12)];
        assert!(test_tall_ternary_smt_impl::<101>(entries, 0).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 1).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 2).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 33).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 34).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 35).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 77).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 78).is_ok());
        assert!(test_tall_ternary_smt_impl::<101>(entries, 79).is_ok());
    }

    fn test_full_binary_smt_impl<I: IntoIterator<Item = (u64, u64)>>(
        entries: I,
        key: u64,
    ) -> Result<()> {
        let tree = {
            let mut tree = get_empty_binary_tree(256);
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

        let chip = FullBinaryChip::<Scalar, BlueSkyConfig3>::new(path);
        assert_eq!(chip.stage_width(), 91);
        assert_eq!(chip.stage_height(), 1);
        assert_eq!(chip.width(), std::cmp::max(257, 91));
        assert_eq!(chip.height(), 259);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), chip.height() + 1);
        assert_eq!(circuit.degree_bound(), 512);
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
    #[ignore]
    fn test_full_binary_smt_empty() {
        assert!(test_full_binary_smt_impl([], 0).is_ok());
        assert!(test_full_binary_smt_impl([], 1).is_ok());
        assert!(test_full_binary_smt_impl([], 2).is_ok());
        assert!(test_full_binary_smt_impl([], 3).is_ok());
        assert!(test_full_binary_smt_impl([], 4).is_ok());
        assert!(test_full_binary_smt_impl([], 5).is_ok());
    }

    #[test]
    #[ignore]
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
    #[ignore]
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
        let tree = {
            let mut tree = get_empty_ternary_tree(161);
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

        let chip = FullTernaryChip::<Scalar, BlueSkyConfig4>::new(path);
        assert_eq!(chip.stage_width(), 103);
        assert_eq!(chip.stage_height(), 1);
        assert_eq!(chip.width(), std::cmp::max(163, 103));
        assert_eq!(chip.height(), 165);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), chip.height() + 1);
        assert_eq!(circuit.degree_bound(), 256);
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
    #[ignore]
    fn test_full_ternary_smt_empty() {
        assert!(test_full_ternary_smt_impl([], 0).is_ok());
        assert!(test_full_ternary_smt_impl([], 1).is_ok());
        assert!(test_full_ternary_smt_impl([], 2).is_ok());
        assert!(test_full_ternary_smt_impl([], 3).is_ok());
        assert!(test_full_ternary_smt_impl([], 4).is_ok());
        assert!(test_full_ternary_smt_impl([], 5).is_ok());
    }

    #[test]
    #[ignore]
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
    #[ignore]
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
