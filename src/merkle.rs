use crate::sponge;
use crate::xits;
use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, make_const, rvar,
};

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height `H`.
///
/// WARNING: `H` must be strictly less than 255. Do NOT use this chip if `H` spans the full BlueSky
/// range, as in that case the bit decomposition of the key would be UNSAFE! Use the
/// [`FullBinaryChip`] below instead.
#[derive(Debug, Clone)]
pub struct BinaryChip<const H: usize> {
    decomposer: xits::BitDecomposerChip<H>,
    hasher: sponge::Poseidon1ChipT3<2>,
    path: [[Scalar; 2]; H],
}

impl<const H: usize> BinaryChip<H> {
    pub fn new(path: [[Scalar; 2]; H]) -> Self {
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher: sponge::Poseidon1ChipT3::default(),
            path,
        }
    }
}

impl<const H: usize> Default for BinaryChip<H> {
    fn default() -> Self {
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher: sponge::Poseidon1ChipT3::default(),
            path: [[Scalar::ZERO; 2]; H],
        }
    }
}

impl<const H: usize> BinaryChip<H> {
    const STAGE_WIDTH: usize = 3;

    const SELECTOR_HEIGHT: usize = 2;

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

impl<const H: usize> PlonkChip<2, 1> for BinaryChip<H> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), Self::STAGE_WIDTH * H)
    }

    fn height(&self) -> usize {
        self.decomposer.height() + Self::SELECTOR_HEIGHT + self.hasher.height()
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.build(view, [key])?;
        let mut hash = value;
        for i in 0..H {
            let bit = bits[i];
            let mut view = view.sub(
                self.decomposer.height(),
                i * Self::STAGE_WIDTH,
                Self::STAGE_WIDTH,
            );
            let inputs = [view.cell(1, 0).into(), view.cell(1, 1).into()];
            [hash, _] = view
                .sub_fn(0, 0, Self::STAGE_WIDTH, |view| {
                    self.build_input_selector(view, hash, bit)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher, inputs)?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 2],
    ) -> Result<[CellOrUnconstrained; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.witness(view, [key])?;
        let mut hash = value;
        for i in 0..H {
            let mut view = view.sub(
                self.decomposer.height(),
                i * Self::STAGE_WIDTH,
                Self::STAGE_WIDTH,
            );
            let inputs = [view.cell(1, 0).into(), view.cell(1, 1).into()];
            [hash, _] = view
                .sub_fn(0, 0, Self::STAGE_WIDTH, |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher, inputs)?;
        }
        Ok([hash])
    }
}

// TODO

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};

    const BLOWUP_LOG2: usize = 3;

    fn test_binary_smt<const H: usize>(
        key: u64,
        value: u64,
        path: [[Scalar; 2]; H],
        expected_root_hash: Scalar,
    ) -> Result<()> {
        let key = Scalar::from(key);
        let value = Scalar::from(value);
        let chip = BinaryChip::<H>::new(path);
        assert_eq!(chip.width(), H * 3);
        assert_eq!(chip.height(), 200);
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
        let openings = circuit.verify(&proof, options)?;
        assert_eq!(openings[&root_hash], expected_root_hash);
        Ok(())
    }

    #[test]
    fn test_binary_smt_height_one_1() {
        let path = [[from_const(12), from_const(34)]];
        let root_hash =
            parse_scalar("0x45470d74563e5e49fe3bd2a161b36116e3c6a6a2f9c105bfe8c2599ff6116b06");
        assert!(test_binary_smt::<1>(0, 12, path, root_hash).is_ok());
        assert!(test_binary_smt::<1>(1, 34, path, root_hash).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_2() {
        let path = [[from_const(34), from_const(12)]];
        let root_hash =
            parse_scalar("0x6a6ca65c7ab651a6e7751e7a23df1d7ff66f745f1b09f4b39df2dfeb4e137422");
        assert!(test_binary_smt::<1>(0, 34, path, root_hash).is_ok());
        assert!(test_binary_smt::<1>(1, 12, path, root_hash).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78)]];
        let root_hash =
            parse_scalar("0x1ba4c686a3529d3bfc13890b2e1438b7adf780e2978cb2cabdd47653f402e8fe");
        assert!(test_binary_smt::<1>(0, 56, path, root_hash).is_ok());
        assert!(test_binary_smt::<1>(1, 78, path, root_hash).is_ok());
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
        assert!(test_binary_smt::<2>(0, 12, path, root_hash).is_ok());
        assert!(test_binary_smt::<2>(1, 34, path, root_hash).is_ok());
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
        assert!(test_binary_smt::<2>(2, 56, path, root_hash).is_ok());
        assert!(test_binary_smt::<2>(3, 78, path, root_hash).is_ok());
    }

    // TODO
}
