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
    const SELECTOR_HEIGHT: usize = 2;

    pub fn new(path: [[Scalar; 2]; H]) -> Self {
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher: sponge::Poseidon1ChipT3::default(),
            path,
        }
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

impl<const H: usize> PlonkChip<2, 1> for BinaryChip<H> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.hasher.width() * H)
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
                i * self.hasher.width(),
                self.hasher.width(),
            );
            let inputs = [view.cell(1, 0).into(), view.cell(1, 1).into()];
            [hash, _] = view
                .sub_fn(0, 0, self.hasher.width(), |view| {
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
                i * self.hasher.width(),
                self.hasher.width(),
            );
            let inputs = [view.cell(1, 0).into(), view.cell(1, 1).into()];
            [hash, _] = view
                .sub_fn(0, 0, self.hasher.width(), |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(Self::SELECTOR_HEIGHT, 0, &self.hasher, inputs)?;
        }
        Ok([hash])
    }
}

// TODO: TernaryChip

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height 256.
///
/// The keys of such a tree span the full BlueSky range. Internally this chip uses a
/// [`xits::FullBitDecomposerChip`], making the 256-bit decomposition safe at the cost of some extra
/// constraints.
///
/// If you don't need 256- or 255-bit keys use [`BinaryChip`].
#[derive(Debug, Clone)]
pub struct FullBinaryChip {
    decomposer: xits::FullBitDecomposerChip,
    hasher: sponge::Poseidon1ChipT3<2>,
    path: [[Scalar; 2]; 256],
}

impl Default for FullBinaryChip {
    fn default() -> Self {
        Self {
            decomposer: xits::FullBitDecomposerChip::default(),
            hasher: sponge::Poseidon1ChipT3::default(),
            path: [[Scalar::ZERO; 2]; 256],
        }
    }
}

impl FullBinaryChip {
    const SELECTOR_HEIGHT: usize = 2;

    pub fn new(path: [[Scalar; 2]; 256]) -> Self {
        Self {
            decomposer: xits::FullBitDecomposerChip::default(),
            hasher: sponge::Poseidon1ChipT3::default(),
            path,
        }
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

impl PlonkChip<2, 1> for FullBinaryChip {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.hasher.width() * 256)
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
        for i in 0..256 {
            let bit = bits[i];
            let mut view = view.sub(
                self.decomposer.height(),
                i * self.hasher.width(),
                self.hasher.width(),
            );
            let inputs = [view.cell(1, 0).into(), view.cell(1, 1).into()];
            [hash, _] = view
                .sub_fn(0, 0, self.hasher.width(), |view| {
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
        for i in 0..256 {
            let mut view = view.sub(
                self.decomposer.height(),
                i * self.hasher.width(),
                self.hasher.width(),
            );
            let inputs = [view.cell(1, 0).into(), view.cell(1, 1).into()];
            [hash, _] = view
                .sub_fn(0, 0, self.hasher.width(), |view| {
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
    use primitive_types::U256;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_ff::Field256;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon as poseidon1;
    use std::fmt::Debug;
    use std::sync::{Arc, LazyLock};

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

    fn test_full_binary_smt_impl<I: IntoIterator<Item = (u64, u64)>>(
        entries: I,
        key: u64,
    ) -> Result<()> {
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

        let chip = FullBinaryChip::new(path);
        assert_eq!(chip.width(), 256 * 3);
        assert_eq!(chip.height(), 202);
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

    // TODO
}
