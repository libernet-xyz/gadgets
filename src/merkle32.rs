//! Merkle lookups over Sparse Merkle Trees whose keys and values are elements of the 256-bit
//! extension of a 32-bit field.
//!
//! Every chip hashes the tree nodes with a Poseidon1 permutation whose state holds the children
//! followed by 16 capacity elements (~512 bits), warranting ~128 bits of preimage resistance in a
//! quantum setting. A child takes 8 elements, so the state width of the permutation is 32 in the
//! binary chip. The hash of a node is made of the first 8 elements of the output.

use crate::poseidon1;
use crate::xits;
use anyhow::{Result, anyhow};
use starkom_ff::{Field256, PrimeField32};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, make_const, rvar, var,
};
use starkom_poseidon::Config as PoseidonConfig;

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height `H`.
///
/// Key and values of the Merkle tree are defined over the 256-bit extension field and must be
/// decomposed into 32-bit base field elements off-circuit. Every node is hashed with the Poseidon1
/// permutation of state width 32: the input state holds the two children (8 elements each) followed
/// by 16 zeros.
#[derive(Debug, Clone)]
pub struct BinaryChip<F: PrimeField32, const H: usize, C: PoseidonConfig<F, 32>> {
    hasher: poseidon1::PermutationChip<F, C, 32>,
    path: [[F; 16]; H],
}

impl<F: PrimeField32, const H: usize, C: PoseidonConfig<F, 32>> Default for BinaryChip<F, H, C> {
    fn default() -> Self {
        Self::new([[F::ZERO; 16]; H])
    }
}

impl<F: PrimeField32, const H: usize, C: PoseidonConfig<F, 32>> BinaryChip<F, H, C> {
    const SELECTOR_WIDTH: usize = 51;

    pub fn new(path: [[F; 16]; H]) -> Self {
        assert!(H <= 256);
        Self {
            hasher: poseidon1::PermutationChip::default(),
            path,
        }
    }

    fn modulus_bit(i: usize) -> F {
        let modulus: u32 = F::MAX.to_u32() + 1;
        ((modulus >> i) & 1).try_into().unwrap()
    }

    fn modulus_has_bits_above(i: usize) -> bool {
        let modulus: u32 = F::MAX.to_u32() + 1;
        modulus.checked_shr(i as u32 + 1).unwrap_or(0) != 0
    }

    /// Selector layout:
    ///
    /// 0        8        16   17   18   19       27       35               51
    /// +--------+--------+----+----+----+--------+--------+----------------+
    /// |   H1   |   H2   | B  | C  | S  |   I1   |   I2   |       0        |
    /// +--------+--------+----+----+----+--------+--------+----------------+
    ///
    /// H1 = leaf-to-root path hash (8 digits)
    /// H2 = peer hash (8 digits, unconstrained)
    /// B = key bit (1 digit)
    /// C = partial result of comparison (1 digit)
    /// S = key bit sum (1 digit)
    /// I1 = left-hand-side input hash (8 digits, either H1 or H2)
    /// I2 = right-hand-side input hash (8 digits, either H1 or H2)
    /// 0 = the zero scalars used as input capacity (16 digits)
    ///
    /// Total [selector width](`Self::SELECTOR_WIDTH`): 51 field elements.
    ///
    /// The lookup key is an octet from the 256-bit extension field and it's decomposed into its
    /// 32-bit base field elements off-circuit, in little-endian order. The chip further decomposes
    /// each base field element (or "digit") into its 32 bits using the full decomposition pattern,
    /// similar to [`xits::FullBitDecomposerChip256`]. The B column contains all 256 decomposed bits
    /// in little-endian order, with each group of 32 bits representing a base field element (even
    /// if the base field is e.g. a 31-bit one). So for example: the first row of the chip contains
    /// the least significant bit of the least significant digit, the 31st row contains the most
    /// significant bit of the least significant digit, the 32nd row contains the least significant
    /// bit of the second-least significant digit, and so on up to the last (256th) row which
    /// contains the most significant bit of the most significant digit and of the whole key.
    ///
    /// The S column reconstructs the original base field elements of the key with Horner's method,
    /// from the most significant bit of each digit to the least significant one: every row holds
    /// twice the S value of the next row plus its own bit, except for the rows of the most
    /// significant bits of the digits (and the last row of the chip), which hold just their bit.
    /// So the 8 reconstructed digits are located at rows `32i`: 0, 32, 64, 96, 128, 160, 192, and
    /// 224, and each of them must be constrained to equal the corresponding digit of the original
    /// key. That way all rows but one in 32 (and the last one) have the same S constraint, so they
    /// share the same selector. If `H` is less than 256 the most significant digits of the key that
    /// don't have any rows in the chip are constrained to zero.
    ///
    /// Since the S sum only proves *congruency* rather than absolute equality, a malicious prover
    /// might inject key digits that only equal the provided lookup key after wraparound. To prevent
    /// that, this chip explicitly compares each key digit against the base field modulus bit by
    /// bit. The comparison is performed in the C column, which contains -1 if the comparison so far
    /// has determined that the digit is strictly less than the modulus, 0 if they're equal, and 1
    /// if the digit is strictly greater. `build_input_selector` constrains the C cell of the least
    /// significant bit of each digit to -1.
    ///
    /// NOTE: since the rows of the chip are ordered in little-endian (the first row calculates the
    /// leaf hash with the least significant bit of the key, the last row calculates the root hash
    /// with the most significant bit of the key), each C cell must refer to the *next* one
    /// (`rvar(17, 1)`) to determine whether or not the comparison has been already resolved by a
    /// higher bit. Likewise, each S cell refers to the next one (`rvar(18, 1)`).
    ///
    /// Note that the last 32 elements of the selector are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: [Option<Cell>; 8],
        i: usize,
    ) {
        for i in 0..8 {
            view.connect(hash[i], view.cell(0, i).into());
        }
        view.add_gate(0, var(16) * (var(16) - make_const(F::ONE)));

        let j = i & 31;
        if j != 31 && i < H - 1 {
            view.add_gate(
                0,
                rvar(17, 1)
                    + (make_const(F::ONE) - (rvar(17, 1) ^ 2))
                        * (var(16) - make_const(Self::modulus_bit(j)))
                    - var(17),
            );
            view.add_gate(0, rvar(18, 1) * 2 + var(16) - var(18));
        } else {
            if Self::modulus_has_bits_above(j) {
                view.add_gate(0, var(17) + make_const(F::ONE));
            } else {
                view.add_gate(0, var(16) - make_const(Self::modulus_bit(j)) - var(17));
            }
            view.add_gate(0, var(16) - var(18));
        }
        if j == 0 {
            view.add_gate(0, var(17) + make_const(F::ONE));
        }

        for j in 0..8 {
            view.add_gate(
                0,
                var(16) * var(8 + j) + (make_const(F::ONE) - var(16)) * var(j) - var(19 + j),
            );
            view.add_gate(
                0,
                var(16) * var(j) + (make_const(F::ONE) - var(16)) * var(8 + j) - var(27 + j),
            );
        }
        for j in 0..16 {
            view.add_gate(0, var(35 + j));
        }
    }

    /// See [`Self::build_input_selector`] for the layout.
    ///
    /// This function fills in all selector cells of the i-th row except the C and S columns, which
    /// cannot be determined without knowing the entries from the next row. The C and S columns are
    /// filled in separately by the caller.
    fn witness_input_selector(&self, view: &mut impl WitnessView<F>, bits: &[F], i: usize) {
        let bit = bits[i];
        if bit != F::ZERO {
            for j in 0..8 {
                view.set(view.cell(0, j), self.path[i][8 + j]);
                view.set(view.cell(0, 8 + j), self.path[i][j]);
            }
        } else {
            for j in 0..8 {
                view.set(view.cell(0, j), self.path[i][j]);
                view.set(view.cell(0, 8 + j), self.path[i][8 + j]);
            }
        }
        view.set(view.cell(0, 16), bit);
        for j in 0..8 {
            view.set(view.cell(0, 19 + j), self.path[i][j]);
            view.set(view.cell(0, 27 + j), self.path[i][8 + j]);
        }
        for j in 0..16 {
            view.set(view.cell(0, 35 + j), F::ZERO);
        }
    }
}

impl<F: PrimeField32, const H: usize, C: PoseidonConfig<F, 32>> PlonkChip<F, 16, 8>
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
        inputs: [Option<Cell>; 16],
    ) -> Result<[Option<Cell>; 8]> {
        let (key, value) = {
            let mut key = [None; 8];
            let mut value = [None; 8];
            key.copy_from_slice(&inputs[0..8]);
            value.copy_from_slice(&inputs[8..16]);
            (key, value)
        };
        let width = self.width();
        let mut hash = value;
        for i in 0..H {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 19 + i).into());
            let output = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.build_input_selector(view, hash, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
            hash.copy_from_slice(&output[0..8]);
        }
        for i in 0..8 {
            if 32 * i < H {
                view.connect(key[i], view.cell(32 * i, 18).into());
            } else {
                // The first capacity zero of the first row.
                view.connect(key[i], view.cell(0, 35).into());
            }
        }
        Ok(hash)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 16],
    ) -> Result<[CellOrUnconstrained<F>; 8]> {
        let key: [F; 8] = std::array::from_fn(|i| view.get(inputs[i]));
        let value: [CellOrUnconstrained<F>; 8] = std::array::from_fn(|i| inputs[8 + i]);
        let bits: [F; 256] = key
            .map(|digit| xits::decompose_scalar_bits::<F, 32>(digit))
            .concat()
            .try_into()
            .unwrap();
        for i in H..256 {
            if bits[i] != F::ZERO {
                return Err(anyhow!("key overflow"));
            }
        }
        let width = self.width();
        let mut hash = value;
        for i in 0..H {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 19 + i).into());
            let output = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?;
            hash.copy_from_slice(&output[0..8]);
        }
        let mut cmp = F::ZERO;
        let mut sum = F::ZERO;
        for i in (0..H).rev() {
            let j = i & 31;
            if j != 31 && i < H - 1 {
                cmp += (F::ONE - cmp.square()) * (bits[i] - Self::modulus_bit(j));
                sum = sum * F::from(2u8) + bits[i];
            } else {
                if Self::modulus_has_bits_above(j) {
                    cmp = -F::ONE;
                } else {
                    cmp = bits[i] - Self::modulus_bit(j);
                }
                sum = bits[i];
            }
            view.set(view.cell(i, 17), cmp);
            view.set(view.cell(i, 18), sum);
        }
        Ok(hash)
    }
}

// TODO: TernaryChip

#[cfg(all(test, feature = "koalabear"))]
mod tests {
    use super::*;
    use primitive_types::{H256, U256};
    use starkom_ff::Field;
    use starkom_koalabear::{KB, KB8};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{Circuit, CircuitBuilder, CircuitOptions, Witness};
    use starkom_poseidon::{self as poseidon1, KoalaBearConfig32};
    use std::collections::BTreeMap;
    use std::fmt::Debug;
    use std::str::FromStr;
    use std::sync::{Arc, LazyLock, Mutex};

    const SECURITY_BITS: usize = 128;
    const BLOWUP_LOG2: usize = 3;

    #[inline]
    const fn from_const(value: u32) -> KB {
        KB::from_const(value)
    }

    fn parse<V: FromStr<Err: Debug>>(s: &'static str) -> V {
        s.parse().unwrap()
    }

    fn from_base(digits: [KB; 8]) -> KB8 {
        let modulus: U256 = KB::MODULUS.parse().unwrap();
        let mut value = U256::zero();
        for i in 0..8 {
            value += digits[i].to_u256() * modulus.pow(i.into());
        }
        value.try_into().unwrap()
    }

    fn to_base(value: KB8) -> [KB; 8] {
        let modulus: U256 = KB::MODULUS.parse().unwrap();
        let mut value = value.to_u256();
        let mut digits = [KB::ZERO; 8];
        for i in 0..8 {
            digits[i] = KB::try_from(value % modulus).unwrap();
            value /= modulus;
        }
        digits
    }

    fn test_binary_smt<const H: usize>(
        key: u8,
        value: [KB; 8],
        path: [[KB; 16]; H],
        expected_root_hash: KB8,
        circuit_commitment: H256,
    ) -> Result<()> {
        let key = KB::from(key);
        let chip = BinaryChip::<KB, H, KoalaBearConfig32>::new(path);
        assert_eq!(chip.width(), 401);
        assert_eq!(chip.height(), H);
        let mut builder = CircuitBuilder::default();
        let inputs: [Option<Cell>; 16] = std::array::from_fn(|i| builder.cell(0, i).into());
        let root_hash = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells(root_hash.map(|digit| digit.unwrap()));
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: BLOWUP_LOG2,
                grinding_bits: 20,
                blind: true,
            })
            .unwrap();
        let mut witness = circuit.make_witness();
        let inputs: [Cell; 16] = std::array::from_fn(|i| witness.cell(0, i));
        witness.set(inputs[0], key);
        for i in 1..8 {
            witness.set(inputs[i], KB::ZERO);
        }
        for i in 0..8 {
            witness.set(inputs[8 + i], value[i]);
        }
        let root_hash = witness.sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))?;
        let root_hash = root_hash.map(|digit| match digit {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        });
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<KB8>>(witness, &circuit.default_proof_options())?;
        let circuit = circuit.to_compressed::<Sha2Hash<KB8>>();
        // assert_eq!(circuit.commitment(), circuit_commitment);  // TODO: re-enable
        let openings = circuit.verify(&proof)?;
        let expected_root_hash = to_base(expected_root_hash);
        for i in 0..8 {
            assert_eq!(openings[&root_hash[i]], expected_root_hash[i]);
        }
        Ok(())
    }

    #[test]
    fn test_binary_smt_height_one_1() {
        let value1 = [
            from_const(12),
            from_const(34),
            from_const(56),
            from_const(78),
            from_const(90),
            from_const(12),
            from_const(34),
            from_const(56),
        ];
        let value2 = [
            from_const(78),
            from_const(90),
            from_const(12),
            from_const(34),
            from_const(56),
            from_const(78),
            from_const(90),
            from_const(12),
        ];
        let path = [[value1, value2].concat().try_into().unwrap()];
        let root_hash = parse("0x00eb59cb6aea53a95e9ed4b4eeb46e5beb46e2e5ac005b31e81ec52d72fd41d9");
        let c = parse("0x003b5e8f62e8189a985dd38d93c152d4107e8ac4777a9768398ce204b46b68ae");
        assert!(test_binary_smt::<1>(0, value1, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1>(1, value2, path, root_hash, c).is_ok());
    }

    /// The height of [`BinarySmtFixture`], the smallest one spanning two digits of the key.
    const FIXTURE_HEIGHT: usize = 33;

    /// The circuit of a [`BinaryChip`] of height [`FIXTURE_HEIGHT`]. Building the circuit is slow,
    /// so the fixture is shared and only checks witnesses, without proving.
    struct BinarySmtFixture {
        circuit: Circuit<KB, KB8>,
    }

    impl BinarySmtFixture {
        fn new() -> Self {
            let chip = BinaryChip::<KB, FIXTURE_HEIGHT, KoalaBearConfig32>::default();
            let mut builder = CircuitBuilder::default();
            let inputs = std::array::from_fn(|i| builder.cell(0, i).into());
            let root_hash = builder.sub_chip(1, 0, &chip, inputs).unwrap();
            builder.declare_public_cells(root_hash.map(|digit| digit.unwrap()));
            let circuit = builder
                .build(CircuitOptions {
                    canonicalize_constraints: false,
                    security_bits: SECURITY_BITS,
                    blowup_log2: BLOWUP_LOG2,
                    grinding_bits: 20,
                    blind: true,
                })
                .unwrap();
            Self { circuit }
        }

        fn get() -> &'static Self {
            static FIXTURE: LazyLock<BinarySmtFixture> = LazyLock::new(BinarySmtFixture::new);
            &FIXTURE
        }

        /// Fills in a witness looking up `key` (whose digits are in the first 8 cells of the first
        /// row) and a zero value in a tree whose other nodes are all zero.
        ///
        /// The path follows the bits of the digits of the key like the chip does, which is not the
        /// order of the bits of the 256-bit key used by the [`Node`] trees.
        fn make_witness(&self, key: [KB; 8]) -> Result<Witness<KB>> {
            let bits: Vec<KB> = key
                .into_iter()
                .flat_map(xits::decompose_scalar_bits::<KB, 32>)
                .collect();
            let mut hash = [KB::ZERO; 8];
            let path: [[KB; 16]; FIXTURE_HEIGHT] = std::array::from_fn(|i| {
                let children = if bits[i] != KB::ZERO {
                    [[KB::ZERO; 8], hash]
                } else {
                    [hash, [KB::ZERO; 8]]
                };
                let children: [KB; 16] = children.concat().try_into().unwrap();
                let output =
                    poseidon1::hash::<KoalaBearConfig32, KB, 32, 16, 16>([KB::ZERO; 16], children);
                hash = std::array::from_fn(|j| output[j]);
                children
            });
            let chip = BinaryChip::<KB, FIXTURE_HEIGHT, KoalaBearConfig32>::new(path);
            let mut witness = self.circuit.make_witness();
            let inputs: [Cell; 16] = std::array::from_fn(|i| witness.cell(0, i));
            for i in 0..8 {
                witness.set(inputs[i], key[i]);
                witness.set(inputs[8 + i], KB::ZERO);
            }
            witness.sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))?;
            Ok(witness)
        }
    }

    /// A key spanning all 33 bits of [`BinarySmtFixture`]: the first digit is the largest one.
    const FIXTURE_KEY: [KB; 8] = [
        from_const(0x7f000000),
        from_const(1),
        from_const(0),
        from_const(0),
        from_const(0),
        from_const(0),
        from_const(0),
        from_const(0),
    ];

    #[test]
    fn test_binary_smt_witness() {
        let fixture = BinarySmtFixture::get();
        let witness = fixture.make_witness(FIXTURE_KEY).unwrap();
        assert!(fixture.circuit.check_witness(&witness).is_ok());
    }

    #[test]
    fn test_binary_smt_witness_key_overflow() {
        let fixture = BinarySmtFixture::get();
        assert!(
            fixture
                .make_witness([
                    from_const(0),
                    from_const(2),
                    from_const(0),
                    from_const(0),
                    from_const(0),
                    from_const(0),
                    from_const(0),
                    from_const(0),
                ])
                .is_err()
        );
    }

    fn test_binary_smt_witness_tampered_digit_impl(digit: usize, value: u32) {
        let fixture = BinarySmtFixture::get();
        let mut witness = fixture.make_witness(FIXTURE_KEY).unwrap();
        witness.set(witness.cell(0, digit), from_const(value));
        assert!(fixture.circuit.check_witness(&witness).is_err());
    }

    #[test]
    fn test_binary_smt_witness_tampered_digit() {
        test_binary_smt_witness_tampered_digit_impl(0, 0x7effffff);
        test_binary_smt_witness_tampered_digit_impl(1, 0);
        test_binary_smt_witness_tampered_digit_impl(2, 1);
        test_binary_smt_witness_tampered_digit_impl(7, 1);
    }

    // TODO

    trait Node: 'static + Debug + Send + Sync {
        fn hash(&self) -> KB8;

        fn get_impl(&self, key: &U256) -> KB8;

        fn get(&self, key: KB8) -> KB8 {
            self.get_impl(&key.to_u256())
        }

        fn get_merkle_path_impl(&self, key: &U256) -> Vec<Vec<KB8>>;

        fn get_merkle_path(&self, key: KB8) -> Vec<Vec<KB8>> {
            self.get_merkle_path_impl(&key.to_u256())
        }

        fn put_impl(self: Arc<Self>, key: &U256, value: KB8) -> Arc<dyn Node>;

        fn put(self: Arc<Self>, key: KB8, value: KB8) -> Arc<dyn Node> {
            self.put_impl(&key.to_u256(), value)
        }
    }

    #[derive(Debug, Default, Copy, Clone)]
    struct Leaf(KB8);

    impl Node for Leaf {
        fn hash(&self) -> KB8 {
            self.0
        }

        fn get_impl(&self, _key: &U256) -> KB8 {
            self.0
        }

        fn get_merkle_path_impl(&self, _key: &U256) -> Vec<Vec<KB8>> {
            vec![]
        }

        fn put_impl(self: Arc<Self>, _key: &U256, value: KB8) -> Arc<dyn Node> {
            Arc::new(Leaf(value))
        }
    }

    #[derive(Debug)]
    struct BinaryNode {
        level: usize,
        hash: KB8,
        left: Arc<dyn Node>,
        right: Arc<dyn Node>,
    }

    impl BinaryNode {
        fn new(level: usize, left: Arc<dyn Node>, right: Arc<dyn Node>) -> Arc<dyn Node> {
            let hash = poseidon1::hash::<KoalaBearConfig32, KB, 32, 16, 16>(
                std::array::from_fn(|_| KB::ZERO),
                [to_base(left.hash()), to_base(right.hash())].concat(),
            );
            let hash = from_base(std::array::from_fn(|i| hash[i]));
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
        fn hash(&self) -> KB8 {
            self.hash
        }

        fn get_impl(&self, key: &U256) -> KB8 {
            if self.bit_at(key) {
                self.right.get_impl(key)
            } else {
                self.left.get_impl(key)
            }
        }

        fn get_merkle_path_impl(&self, key: &U256) -> Vec<Vec<KB8>> {
            let mut path = if self.bit_at(key) {
                self.right.get_merkle_path_impl(key)
            } else {
                self.left.get_merkle_path_impl(key)
            };
            path.push(vec![self.left.hash(), self.right.hash()]);
            path
        }

        fn put_impl(self: Arc<Self>, key: &U256, value: KB8) -> Arc<dyn Node> {
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
        let path: [[KB; 16]; H] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| {
                entry
                    .into_iter()
                    .map(to_base)
                    .flatten()
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap()
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = to_base(tree.hash());

        let chip = BinaryChip::<KB, H, KoalaBearConfig32>::new(path);
        assert_eq!(chip.width(), 401);
        assert_eq!(chip.height(), H);

        let mut builder = CircuitBuilder::default();
        let inputs = std::array::from_fn(|i| builder.cell(0, i).into());
        let root_hash = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells(root_hash.map(|digit| digit.unwrap()));
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: BLOWUP_LOG2,
                grinding_bits: 20,
                blind: true,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), chip.height() + 1);
        assert_eq!(circuit.num_columns(), chip.width());

        let mut witness = circuit.make_witness();
        let inputs = std::array::from_fn(|i| witness.cell(0, i));
        let key = to_base(key);
        for i in 0..8 {
            witness.set(inputs[i], key[i]);
        }
        let value = to_base(value);
        for i in 0..8 {
            witness.set(inputs[8 + i], value[i]);
        }
        let root_hash = witness.sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))?;
        let root_hash = root_hash.map(|digit| match digit {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        });

        circuit.check_witness(&witness).unwrap();

        let proof = circuit.prove::<Sha2Hash<KB8>>(witness, &circuit.default_proof_options())?;
        let openings = circuit.verify(&proof)?;
        for i in 0..8 {
            assert_eq!(openings[&root_hash[i]], expected_root_hash[i]);
        }

        Ok(())
    }

    #[test]
    #[ignore]
    fn test_tall_binary_smt_empty() {
        assert!(test_tall_binary_smt_impl::<20>([], 0).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 1).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 2).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 3).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 4).is_ok());
        assert!(test_tall_binary_smt_impl::<20>([], 5).is_ok());
    }

    #[test]
    #[ignore]
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
    #[ignore]
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

    // TODO
}
