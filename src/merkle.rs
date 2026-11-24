//! Merkle lookups over Sparse Merkle Trees.
//!
//! Every chip hashes the tree nodes with a Poseidon1 permutation whose state holds the children
//! followed by two capacity elements, warranting ~128 bits of preimage resistance in a quantum
//! setting. The state width of the permutation is therefore 4 in the binary chips and 5 in the
//! ternary ones. The capacity elements of every permutation are initialized with the two domain
//! separator tags (DSTs) passed to the chip, and the hash of a node is the first element of the
//! output (same squeeze as [`hash0`](`starkom_poseidon::hash0`)).

use crate::poseidon1;
use crate::xits;
use anyhow::Result;
use primitive_types::U256;
use starkom_ff::{Field256, PrimeField256};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, make_const, rvar, var,
};
use starkom_poseidon::Config as PoseidonConfig;

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height `H`.
///
/// Every node is hashed with the Poseidon1 permutation of state width 4. The input state holds the
/// two children followed by the two DSTs of the chip.
///
/// WARNING: `H` must be strictly less than 256. Do NOT use this chip if `H` spans the full `F`
/// range, as in that case the bit decomposition of the key would be UNSAFE! Use the
/// [`FullBinaryChip`] below instead.
#[derive(Debug, Clone)]
pub struct BinaryChip<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> {
    hasher: poseidon1::PermutationChip<F, C, 4>,
    dst: [F; 2],
    path: [[F; 2]; H],
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> Default for BinaryChip<F, H, C> {
    fn default() -> Self {
        Self::new([F::ZERO; 2], [[F::ZERO; 2]; H])
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> BinaryChip<F, H, C> {
    const SELECTOR_WIDTH: usize = 8;

    /// Makes a chip that looks up a key through the Merkle `path`, initializing the capacity of
    /// every permutation with the domain separator tags `dst`.
    pub fn new(dst: [F; 2], path: [[F; 2]; H]) -> Self {
        assert!(H < F::NUM_BITS);
        Self {
            hasher: poseidon1::PermutationChip::default(),
            dst,
            path,
        }
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+----+----+
    /// | H1 | H2 | B  | S  | I1 | I2 | D1 | D2 |
    /// +----+----+----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = peer hash (unconstrained)
    /// B = key bit
    /// S = key bit sum
    /// I1 = left-hand-side input hash (either H1 or H2)
    /// I2 = right-hand-side input hash (either H1 or H2)
    /// D1, D2 = the two DSTs of the chip, used as input capacity
    ///
    /// The B column holds the decomposed bits of the key and the S column reconstructs the original
    /// key with Horner's method, from the most significant bit in the last row to the least
    /// significant one in the first row: every row but the last holds twice the S value of the next
    /// row plus its own bit, so the first row holds the whole key, which must be constrained to
    /// equal the original key. That way all rows but the last have the same constraints, so they
    /// share the same selector.
    ///
    /// Note that the last four elements of the selector are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: Option<Cell>,
        i: usize,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.add_gate(0, var(2) * (var(2) - make_const(F::ONE)));
        if i < H - 1 {
            view.add_gate(0, rvar(3, 1) * 2 + var(2) - var(3));
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
        view.add_gate(0, var(6) - make_const(self.dst[0]));
        view.add_gate(0, var(7) - make_const(self.dst[1]));
    }

    /// See [`Self::build_input_selector`] for the layout.
    ///
    /// This function fills in all selector cells of the i-th row except the S column, which depends
    /// on the next row. The S column is filled in separately by the caller.
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
        view.set(view.cell(0, 4), self.path[i][0]);
        view.set(view.cell(0, 5), self.path[i][1]);
        view.set(view.cell(0, 6), self.dst[0]);
        view.set(view.cell(0, 7), self.dst[1]);
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 4>> PlonkChip<F, 2, 1>
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
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.build_input_selector(view, hash, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        view.connect(key, view.cell(0, 3).into());
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
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        let mut sum = F::ZERO;
        for i in (0..H).rev() {
            sum = sum * F::from(2u8) + bits[i];
            view.set(view.cell(i, 3), sum);
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a ternary Sparse Merkle Tree of height `H`.
///
/// Every node is hashed with the Poseidon1 permutation of state width 5. The input state holds the
/// three children followed by the two DSTs of the chip.
///
/// WARNING: `H` must be strictly less than 162. Do NOT use this chip if `H` spans the full `F`
/// range, as in that case the trit decomposition of the key would be UNSAFE! Use the
/// [`FullTernaryChip`] below instead.
#[derive(Debug, Clone)]
pub struct TernaryChip<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 5>> {
    hasher: poseidon1::PermutationChip<F, C, 5>,
    dst: [F; 2],
    path: [[F; 3]; H],
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 5>> Default for TernaryChip<F, H, C> {
    fn default() -> Self {
        Self::new([F::ZERO; 2], [[F::ZERO; 3]; H])
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 5>> TernaryChip<F, H, C> {
    const SELECTOR_WIDTH: usize = 10;

    /// Makes a chip that looks up a key through the Merkle `path`, initializing the capacity of
    /// every permutation with the domain separator tags `dst`.
    pub fn new(dst: [F; 2], path: [[F; 3]; H]) -> Self {
        // TODO: assert that H is strictly less than the number of trits required to represent a
        // scalar.
        Self {
            hasher: poseidon1::PermutationChip::default(),
            dst,
            path,
        }
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+----+----+----+----+
    /// | H1 | H2 | H3 | T  | S  | I1 | I2 | I3 | D1 | D2 |
    /// +----+----+----+----+----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = first peer hash (unconstrained)
    /// H3 = second peer hash (unconstrained)
    /// T = key trit
    /// S = key trit sum
    /// I1 = first input hash (one of H1, H2, or H3)
    /// I2 = second input hash (one of H1, H2, or H3)
    /// I3 = third input hash (one of H1, H2, or H3)
    /// D1, D2 = the two DSTs of the chip, used as input capacity
    ///
    /// The T column holds the decomposed trits of the key and the S column reconstructs the
    /// original key with Horner's method, from the most significant trit in the last row to the
    /// least significant one in the first row: every row but the last holds three times the S value
    /// of the next row plus its own trit, so the first row holds the whole key, which must be
    /// constrained to equal the original key. That way all rows but the last have the same
    /// constraints, so they share the same selector.
    ///
    /// Note that the last five elements of the selector are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: Option<Cell>,
        i: usize,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.add_gate(
            0,
            var(3) * (var(3) - make_const(F::ONE)) * (var(3) - make_const(F::from(2u8))),
        );
        if i < H - 1 {
            view.add_gate(0, rvar(4, 1) * 3 + var(3) - var(4));
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
        view.add_gate(0, var(8) - make_const(self.dst[0]));
        view.add_gate(0, var(9) - make_const(self.dst[1]));
    }

    /// See [`Self::build_input_selector`] for the layout.
    ///
    /// This function fills in all selector cells of the i-th row except the S column, which
    /// depends on the next row. The S column is filled in separately by the caller.
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
        view.set(view.cell(0, 5), self.path[i][0]);
        view.set(view.cell(0, 6), self.path[i][1]);
        view.set(view.cell(0, 7), self.path[i][2]);
        view.set(view.cell(0, 8), self.dst[0]);
        view.set(view.cell(0, 9), self.dst[1]);
    }
}

impl<F: PrimeField256, const H: usize, C: PoseidonConfig<F, 5>> PlonkChip<F, 2, 1>
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
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.build_input_selector(view, hash, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        view.connect(key, view.cell(0, 4).into());
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
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.witness_input_selector(view, &trits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        let mut sum = F::ZERO;
        for i in (0..H).rev() {
            sum = sum * F::from(3u8) + trits[i];
            view.set(view.cell(i, 4), sum);
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
/// Every node is hashed like in [`BinaryChip`], with the Poseidon1 permutation of state width 4.
#[derive(Debug, Clone)]
pub struct FullBinaryChip<F: PrimeField256, C: PoseidonConfig<F, 4>> {
    hasher: poseidon1::PermutationChip<F, C, 4>,
    dst: [F; 2],
    path: [[F; 2]; 256],
}

impl<F: PrimeField256, C: PoseidonConfig<F, 4>> Default for FullBinaryChip<F, C> {
    fn default() -> Self {
        Self::new([F::ZERO; 2], [[F::ZERO; 2]; 256])
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 4>> FullBinaryChip<F, C> {
    const SELECTOR_WIDTH: usize = 9;

    /// Makes a chip that looks up a key through the Merkle `path`, initializing the capacity of
    /// every permutation with the domain separator tags `dst`.
    pub fn new(dst: [F; 2], path: [[F; 2]; 256]) -> Self {
        Self {
            hasher: poseidon1::PermutationChip::default(),
            dst,
            path,
        }
    }

    fn modulus_bit(i: usize) -> F {
        let modulus: U256 = F::MODULUS.parse().unwrap();
        ((modulus >> i) & U256::one()).as_u32().into()
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+----+----+----+
    /// | H1 | H2 | B  | C  | S  | I1 | I2 | D1 | D2 |
    /// +----+----+----+----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = peer hash (unconstrained)
    /// B = key bit
    /// C = partial result of comparison
    /// S = key bit sum
    /// I1 = left-hand-side input hash (either H1 or H2)
    /// I2 = right-hand-side input hash (either H1 or H2)
    /// D1, D2 = the two DSTs of the chip, used as input capacity
    ///
    /// The B column holds the decomposed bits of the key and the S column reconstructs the original
    /// key with Horner's method, like in [`BinaryChip`]: every row but the last holds twice the S
    /// value of the next row plus its own bit, so the first row holds the whole key, which must be
    /// constrained to equal the original key.
    ///
    /// Since this chip can handle large keys we need to explicitly compare the key against the
    /// field modulus bit by bit in order to prevent aliasing. Without the constraints enforced by
    /// this comparison, the sum in the S column might wrap around. The C column contains the
    /// partial result of the bit-by-bit comparison, with -1 indicating that the key is strictly
    /// less than the modulus, 0 that it's equal, and 1 that it's strictly greater. The first (ie.
    /// least significant) bit of C must be constrained to -1.
    ///
    /// NOTE: since the rows of the chip are ordered in little-endian (the first row calculates the
    /// leaf hash with the least significant bit of the key, the last row calculates the root hash
    /// with the most significant bit of the key), each C bit must refer to the *next* one
    /// (`rvar(3, +1)`) to determine whether or not the comparison has been already resolved by a
    /// higher bit.
    ///
    /// The last four elements of the selector are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: Option<Cell>,
        i: usize,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.add_gate(0, var(2) * (var(2) - make_const(F::ONE)));
        if i < 255 {
            view.add_gate(
                0,
                rvar(3, 1)
                    + (make_const(F::ONE) - (rvar(3, 1) ^ 2))
                        * (var(2) - make_const(Self::modulus_bit(i)))
                    - var(3),
            );
            view.add_gate(0, rvar(4, 1) * 2 + var(2) - var(4));
        } else {
            view.add_gate(0, var(2) - make_const(Self::modulus_bit(255)) - var(3));
            view.add_gate(0, var(2) - var(4));
        }
        view.add_gate(
            0,
            var(2) * var(1) + (make_const(F::ONE) - var(2)) * var(0) - var(5),
        );
        view.add_gate(
            0,
            var(2) * var(0) + (make_const(F::ONE) - var(2)) * var(1) - var(6),
        );
        view.add_gate(0, var(7) - make_const(self.dst[0]));
        view.add_gate(0, var(8) - make_const(self.dst[1]));
    }

    /// See [`Self::build_input_selector`] for the layout.
    ///
    /// This function fills in all selector cells of the i-th row except the C and S columns, which
    /// cannot be determined without knowing the entries from the next row. The C and S columns are
    /// filled in separately by the caller.
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
        view.set(view.cell(0, 5), self.path[i][0]);
        view.set(view.cell(0, 6), self.path[i][1]);
        view.set(view.cell(0, 7), self.dst[0]);
        view.set(view.cell(0, 8), self.dst[1]);
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 4>> PlonkChip<F, 2, 1> for FullBinaryChip<F, C> {
    fn width(&self) -> usize {
        Self::SELECTOR_WIDTH + self.hasher.width()
    }

    fn height(&self) -> usize {
        256
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let width = self.width();
        let mut hash = value;
        for i in 0..256 {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 5 + i).into());
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), Some(1), |view| {
                    self.build_input_selector(view, hash, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        view.add_gate(0, var(3) + make_const(F::ONE));
        view.connect(key, view.cell(0, 4).into());
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 2],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        let [key, value] = inputs;
        let bits = xits::decompose_scalar_bits::<F, 256>(view.get(key));
        let width = self.width();
        let mut hash = value;
        for i in 0..256 {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 5 + i).into());
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), Some(1), |view| {
                    self.witness_input_selector(view, &bits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        let mut cmp = bits[255] - Self::modulus_bit(255);
        view.set(view.cell(255, 3), cmp);
        for i in (0..255).rev() {
            cmp += (F::ONE - cmp.square()) * (bits[i] - Self::modulus_bit(i));
            view.set(view.cell(i, 3), cmp);
        }
        let mut sum = F::ZERO;
        for i in (0..256).rev() {
            sum = sum * F::from(2u8) + bits[i];
            view.set(view.cell(i, 4), sum);
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
///
/// Every node is hashed like in [`TernaryChip`], with the Poseidon1 permutation of state width 5.
#[derive(Debug, Clone)]
pub struct FullTernaryChip<F: PrimeField256, C: PoseidonConfig<F, 5>> {
    hasher: poseidon1::PermutationChip<F, C, 5>,
    dst: [F; 2],
    path: [[F; 3]; 162],
}

impl<F: PrimeField256, C: PoseidonConfig<F, 5>> Default for FullTernaryChip<F, C> {
    fn default() -> Self {
        Self::new([F::ZERO; 2], [[F::ZERO; 3]; 162])
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 5>> FullTernaryChip<F, C> {
    const SELECTOR_WIDTH: usize = 11;

    /// Makes a chip that looks up a key through the Merkle `path`, initializing the capacity of
    /// every permutation with the domain separator tags `dst`.
    pub fn new(dst: [F; 2], path: [[F; 3]; 162]) -> Self {
        Self {
            hasher: poseidon1::PermutationChip::default(),
            dst,
            path,
        }
    }

    fn modulus_trit(i: usize) -> F {
        let modulus: U256 = F::MODULUS.parse().unwrap();
        let three: U256 = 3.into();
        ((modulus / three.pow(i.into())) % three).as_u32().into()
    }

    /// Selector layout:
    ///
    /// +----+----+----+----+----+----+----+----+----+----+----+
    /// | H1 | H2 | H3 | T  | C  | S  | I1 | I2 | I3 | D1 | D2 |
    /// +----+----+----+----+----+----+----+----+----+----+----+
    ///
    /// H1 = leaf-to-root path hash
    /// H2 = first peer hash
    /// H3 = second peer hash
    /// T = key trit
    /// C = partial result of comparison
    /// S = key trit sum
    /// I1 = first input hash (one of H1, H2, or H3)
    /// I2 = second input hash (one of H1, H2, or H3)
    /// I3 = third input hash (one of H1, H2, or H3)
    /// D1, D2 = the two DSTs of the chip, used as input capacity
    ///
    /// The T column holds the decomposed trits of the key and the S column reconstructs the
    /// original key with Horner's method, like in [`TernaryChip`]: every row but the last holds
    /// three times the S value of the next row plus its own trit, so the first row holds the whole
    /// key, which must be constrained to equal the original key.
    ///
    /// Since this chip can handle large keys we need to explicitly compare the key against the
    /// field modulus trit by trit in order to prevent aliasing. Without the constraints enforced by
    /// this comparison, the sum in the S column might wrap around. The C column contains the
    /// partial result of the trit-by-trit comparison, with -1 indicating that the key is strictly
    /// less than the modulus, 0 that it's equal, and 1 that it's strictly greater. The first (ie.
    /// least significant) trit of C must be constrained to -1.
    ///
    /// NOTE: since the rows of the chip are ordered in little-endian (the first row calculates the
    /// leaf hash with the least significant trit of the key, the last row calculates the root hash
    /// with the most significant trit), each C bit must refer to the *next* one (`rvar(3, +1)`) to
    /// determine whether or not the comparison has been already resolved by a higher bit.
    ///
    /// The last five elements of the selector are the permutation input state vector.
    fn build_input_selector<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        hash: Option<Cell>,
        i: usize,
    ) {
        view.connect(hash, view.cell(0, 0).into());
        view.add_gate(
            0,
            var(3) * (var(3) - make_const(F::from(1u8))) * (var(3) - make_const(F::from(2u8))),
        );

        let trit = Self::modulus_trit(i);
        let cmp = ((var(3) - make_const(trit)) * make_const(F::from(7u8))
            - ((var(3) - make_const(trit)) ^ 3))
            / make_const(F::from(6u8));
        if i < 161 {
            view.add_gate(
                0,
                rvar(4, 1) + (make_const(F::ONE) - (rvar(4, 1) ^ 2)) * cmp - var(4),
            );
            view.add_gate(0, rvar(5, 1) * 3 + var(3) - var(5));
        } else {
            view.add_gate(0, cmp - var(4));
            view.add_gate(0, var(3) - var(5));
        }

        let l0 = ((var(3) ^ 2) - var(3) * 3 + 2) / 2;
        let l1 = var(3) * 2 - (var(3) ^ 2);
        let l2 = ((var(3) ^ 2) - var(3)) / 2;
        view.add_gate(
            0,
            l0.clone() * var(0) + l1.clone() * var(1) + l2.clone() * var(1) - var(6),
        );
        view.add_gate(
            0,
            l0.clone() * var(1) + l1.clone() * var(0) + l2.clone() * var(2) - var(7),
        );
        view.add_gate(0, l0 * var(2) + l1 * var(2) + l2 * var(0) - var(8));
        view.add_gate(0, var(9) - make_const(self.dst[0]));
        view.add_gate(0, var(10) - make_const(self.dst[1]));
    }

    /// See [`Self::build_input_selector`] for the layout.
    ///
    /// This function fills in all selector cells of the i-th row except the C and S columns, which
    /// cannot be determined without knowing the entries from the next row. The C and S columns are
    /// filled in separately by the caller.
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
        view.set(view.cell(0, 6), self.path[i][0]);
        view.set(view.cell(0, 7), self.path[i][1]);
        view.set(view.cell(0, 8), self.path[i][2]);
        view.set(view.cell(0, 9), self.dst[0]);
        view.set(view.cell(0, 10), self.dst[1]);
    }
}

impl<F: PrimeField256, C: PoseidonConfig<F, 5>> PlonkChip<F, 2, 1> for FullTernaryChip<F, C> {
    fn width(&self) -> usize {
        Self::SELECTOR_WIDTH + self.hasher.width()
    }

    fn height(&self) -> usize {
        162
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 2],
    ) -> Result<[Option<Cell>; 1]> {
        let [key, value] = inputs;
        let width = self.width();
        let mut hash = value;
        for i in 0..162 {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 6 + i).into());
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.build_input_selector(view, hash, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        view.add_gate(0, var(4) + make_const(F::ONE));
        view.connect(key, view.cell(0, 5).into());
        Ok([hash])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 2],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        let [key, value] = inputs;
        let trits = xits::decompose_scalar_trits::<F, 162>(view.get(key));
        let width = self.width();
        let mut hash = value;
        for i in 0..162 {
            let mut view = view.sub(i, 0, width.into(), Some(1));
            let inputs = std::array::from_fn(|i| view.cell(0, 6 + i).into());
            hash = view
                .sub_fn(0, 0, Self::SELECTOR_WIDTH.into(), None, |view| {
                    self.witness_input_selector(view, &trits, i)
                })
                .sub_chip(0, Self::SELECTOR_WIDTH, &self.hasher, inputs)?[0];
        }
        let diff = trits[161] - Self::modulus_trit(161);
        let mut cmp = (diff * F::from(7u8) - diff.cube()) / F::from(6u8);
        view.set(view.cell(161, 4), cmp);
        for i in (0..161).rev() {
            let diff = trits[i] - Self::modulus_trit(i);
            cmp += (F::ONE - cmp.square()) * ((diff * F::from(7u8) - diff.cube()) / F::from(6u8));
            view.set(view.cell(i, 4), cmp);
        }
        let mut sum = F::ZERO;
        for i in (0..162).rev() {
            sum = sum * F::from(3u8) + trits[i];
            view.set(view.cell(i, 5), sum);
        }
        Ok([hash])
    }
}

#[cfg(all(test, feature = "bluesky"))]
mod tests {
    use super::*;
    use crate::utils::hash_to_scalar;
    use primitive_types::{H256, U256};
    use starkom_bluesky::{Scalar, from_const, parse_scalar};
    use starkom_ff::Field;
    use starkom_pcs::hash::{Keccak256Hash, Sha2Hash};
    use starkom_plonk::{CircuitBuilder, CircuitOptions, ProofOptions};
    use starkom_poseidon::{
        self as poseidon1, BlueSkyConfigT4X5, BlueSkyConfigT4X7, BlueSkyConfigT5X5,
        BlueSkyConfigT5X7,
    };
    use std::any::TypeId;
    use std::collections::BTreeMap;
    use std::fmt::Debug;
    use std::marker::PhantomData;
    use std::sync::{Arc, LazyLock, Mutex};

    const SECURITY_BITS: usize = 128;
    const BLOWUP_LOG2: usize = 3;

    /// The domain separator tags of the chips under test and of the trees hashed like them.
    static DST: LazyLock<[Scalar; 2]> =
        LazyLock::new(|| [hash_to_scalar(b"test/dst/0"), hash_to_scalar(b"test/dst/1")]);

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
        let chip = BinaryChip::<Scalar, H, BlueSkyConfigT4X5>::new(*DST, path);
        assert_eq!(chip.width(), 103);
        assert_eq!(chip.height(), H);
        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: BLOWUP_LOG2,
                grinding_bits: 0,
                blind: true,
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
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, &circuit.default_proof_options())?;
        let circuit = circuit.to_compressed::<Sha2Hash<Scalar>>();
        // assert_eq!(circuit.commitment(), circuit_commitment);  // TODO: re-enable
        let openings = circuit.verify(&proof)?;
        assert_eq!(openings[&root_hash], expected_root_hash);
        Ok(())
    }

    #[test]
    fn test_binary_smt_height_one_1() {
        let path = [[from_const(12), from_const(34)]];
        let root_hash =
            parse_scalar("0x5c7e89842c3aaeea2f048e1d6de663af6ad11dce7d98d17b19cec47bf53c37b0");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1>(1, 34, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_2() {
        let path = [[from_const(34), from_const(12)]];
        let root_hash =
            parse_scalar("0x24f2dd48eeacd81ea136caab8718c4a9aa9447bef333bdb57b22cd3260923157");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1>(0, 34, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1>(1, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78)]];
        let root_hash =
            parse_scalar("0x76eb513de049a9f257d5e041ba17f3c7ce2e4a9341d39336df809ca1f1e23e1f");
        let c = parse_hash("0x9ea543dc5d7b98c872c7770f45442e1b682de56c0d8b739338ec877aab563285");
        assert!(test_binary_smt::<1>(0, 56, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<1>(1, 78, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_1() {
        let path = [
            [from_const(12), from_const(34)],
            [
                parse_scalar("0x5c7e89842c3aaeea2f048e1d6de663af6ad11dce7d98d17b19cec47bf53c37b0"),
                parse_scalar("0x76eb513de049a9f257d5e041ba17f3c7ce2e4a9341d39336df809ca1f1e23e1f"),
            ],
        ];
        let root_hash =
            parse_scalar("0x6f84a941a15c857ea77930bf2d3dddf58fbe61d29eb9d3da229d2158bdbd7445");
        let c = parse_hash("0x4624a1fc0141a8d753764723b67f749af762216c30f212b4111c6efee396f361");
        assert!(test_binary_smt::<2>(0, 12, path, root_hash, c).is_ok());
        assert!(test_binary_smt::<2>(1, 34, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_2() {
        let path = [
            [from_const(56), from_const(78)],
            [
                parse_scalar("0x5c7e89842c3aaeea2f048e1d6de663af6ad11dce7d98d17b19cec47bf53c37b0"),
                parse_scalar("0x76eb513de049a9f257d5e041ba17f3c7ce2e4a9341d39336df809ca1f1e23e1f"),
            ],
        ];
        let root_hash =
            parse_scalar("0x6f84a941a15c857ea77930bf2d3dddf58fbe61d29eb9d3da229d2158bdbd7445");
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
        let chip = TernaryChip::<Scalar, H, BlueSkyConfigT5X5>::new(*DST, path);
        assert_eq!(chip.width(), 115);
        assert_eq!(chip.height(), H);
        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: BLOWUP_LOG2,
                grinding_bits: 0,
                blind: true,
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
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, &circuit.default_proof_options())?;
        let circuit = circuit.to_compressed::<Sha2Hash<Scalar>>();
        // assert_eq!(circuit.commitment(), circuit_commitment);  // TODO: re-enable
        let openings = circuit.verify(&proof)?;
        assert_eq!(openings[&root_hash], expected_root_hash);
        Ok(())
    }

    #[test]
    fn test_ternary_smt_height_one_1() {
        let path = [[from_const(12), from_const(34), from_const(56)]];
        let root_hash =
            parse_scalar("0x475927bdf7f82c894db3074d32f1602b36df740779717e47b4cb035660f234ff");
        let c = parse_hash("0x512035d2b2db72ca18822ccd4cb4ce7d8455ffa95199f769260e89d31f4ca9f9");
        assert!(test_ternary_smt::<1>(0, 12, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(1, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(2, 56, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_2() {
        let path = [[from_const(34), from_const(56), from_const(12)]];
        let root_hash =
            parse_scalar("0x18ae6c0c541342ec2a98f155ad29bd056f7e10ca4aa384ff0fadd598d0e1a2a4");
        let c = parse_hash("0x512035d2b2db72ca18822ccd4cb4ce7d8455ffa95199f769260e89d31f4ca9f9");
        assert!(test_ternary_smt::<1>(0, 34, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(1, 56, path, root_hash, c).is_ok());
        assert!(test_ternary_smt::<1>(2, 12, path, root_hash, c).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78), from_const(90)]];
        let root_hash =
            parse_scalar("0x72297995c4fb54677035c53979a65dc027d50573af8489f3e9b1ac4b070a382a");
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
                parse_scalar("0x475927bdf7f82c894db3074d32f1602b36df740779717e47b4cb035660f234ff"),
                parse_scalar("0x18ae6c0c541342ec2a98f155ad29bd056f7e10ca4aa384ff0fadd598d0e1a2a4"),
                parse_scalar("0x72297995c4fb54677035c53979a65dc027d50573af8489f3e9b1ac4b070a382a"),
            ],
        ];
        let root_hash =
            parse_scalar("0x2d1ebbf130db703af5cc8013102a7da3df135585a8cc686d1134a1663909fdaf");
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
                parse_scalar("0x475927bdf7f82c894db3074d32f1602b36df740779717e47b4cb035660f234ff"),
                parse_scalar("0x18ae6c0c541342ec2a98f155ad29bd056f7e10ca4aa384ff0fadd598d0e1a2a4"),
                parse_scalar("0x72297995c4fb54677035c53979a65dc027d50573af8489f3e9b1ac4b070a382a"),
            ],
        ];
        let root_hash =
            parse_scalar("0x2d1ebbf130db703af5cc8013102a7da3df135585a8cc686d1134a1663909fdaf");
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
                parse_scalar("0x475927bdf7f82c894db3074d32f1602b36df740779717e47b4cb035660f234ff"),
                parse_scalar("0x18ae6c0c541342ec2a98f155ad29bd056f7e10ca4aa384ff0fadd598d0e1a2a4"),
                parse_scalar("0x72297995c4fb54677035c53979a65dc027d50573af8489f3e9b1ac4b070a382a"),
            ],
        ];
        let root_hash =
            parse_scalar("0x2d1ebbf130db703af5cc8013102a7da3df135585a8cc686d1134a1663909fdaf");
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

    /// A node of a binary tree hashed like the chips with the Poseidon configuration `C`, which
    /// must match the one of the chip under test.
    struct BinaryNode<C: PoseidonConfig<Scalar, 4> + 'static> {
        level: usize,
        hash: Scalar,
        left: Arc<dyn Node>,
        right: Arc<dyn Node>,
        _data: PhantomData<fn() -> C>,
    }

    impl<C: PoseidonConfig<Scalar, 4> + 'static> Debug for BinaryNode<C> {
        fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            formatter
                .debug_struct("BinaryNode")
                .field("level", &self.level)
                .field("hash", &self.hash)
                .field("left", &self.left)
                .field("right", &self.right)
                .finish()
        }
    }

    impl<C: PoseidonConfig<Scalar, 4> + 'static> BinaryNode<C> {
        fn new(level: usize, left: Arc<dyn Node>, right: Arc<dyn Node>) -> Arc<dyn Node> {
            let state = [left.hash(), right.hash(), DST[0], DST[1]];
            let hash = poseidon1::permutation::<C, Scalar, 4>(state)[0];
            Arc::new(BinaryNode::<C> {
                level,
                hash,
                left,
                right,
                _data: PhantomData,
            })
        }

        fn bit_at(&self, key: &U256) -> bool {
            (key >> (self.level - 1)) & U256::one() != U256::zero()
        }
    }

    impl<C: PoseidonConfig<Scalar, 4> + 'static> Node for BinaryNode<C> {
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

    /// A node of a ternary tree hashed like the chips with the Poseidon configuration `C`, which
    /// must match the one of the chip under test.
    struct TernaryNode<C: PoseidonConfig<Scalar, 5> + 'static> {
        level: usize,
        hash: Scalar,
        children: [Arc<dyn Node>; 3],
        _data: PhantomData<fn() -> C>,
    }

    impl<C: PoseidonConfig<Scalar, 5> + 'static> Debug for TernaryNode<C> {
        fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            formatter
                .debug_struct("TernaryNode")
                .field("level", &self.level)
                .field("hash", &self.hash)
                .field("children", &self.children)
                .finish()
        }
    }

    impl<C: PoseidonConfig<Scalar, 5> + 'static> TernaryNode<C> {
        fn new(level: usize, children: [Arc<dyn Node>; 3]) -> Arc<dyn Node> {
            let state = [
                children[0].hash(),
                children[1].hash(),
                children[2].hash(),
                DST[0],
                DST[1],
            ];
            let hash = poseidon1::permutation::<C, Scalar, 5>(state)[0];
            Arc::new(TernaryNode::<C> {
                level,
                hash,
                children,
                _data: PhantomData,
            })
        }

        fn trit_at(&self, key: &U256) -> usize {
            let divisor = U256::from(3).pow((self.level - 1).into());
            ((key / divisor) % 3).try_into().unwrap()
        }
    }

    impl<C: PoseidonConfig<Scalar, 5> + 'static> Node for TernaryNode<C> {
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

    fn get_empty_binary_tree_locked<C: PoseidonConfig<Scalar, 4> + 'static>(
        nodes_by_level: &mut BTreeMap<usize, Arc<dyn Node>>,
        level: usize,
    ) -> Arc<dyn Node> {
        match nodes_by_level.get_mut(&level) {
            Some(node) => node.clone(),
            None => {
                let node = if level > 0 {
                    let child = get_empty_binary_tree_locked::<C>(nodes_by_level, level - 1);
                    BinaryNode::<C>::new(level, child.clone(), child.clone())
                } else {
                    Arc::new(Leaf::default())
                };
                nodes_by_level.insert(level, node.clone());
                node
            }
        }
    }

    fn get_empty_binary_tree<C: PoseidonConfig<Scalar, 4> + 'static>(
        level: usize,
    ) -> Arc<dyn Node> {
        static NODES: LazyLock<Mutex<BTreeMap<TypeId, BTreeMap<usize, Arc<dyn Node>>>>> =
            LazyLock::new(|| Mutex::new(BTreeMap::default()));
        let mut nodes_by_config = NODES.lock().unwrap();
        let nodes_by_level = nodes_by_config.entry(TypeId::of::<C>()).or_default();
        get_empty_binary_tree_locked::<C>(nodes_by_level, level)
    }

    fn get_empty_ternary_tree_locked<C: PoseidonConfig<Scalar, 5> + 'static>(
        nodes_by_level: &mut BTreeMap<usize, Arc<dyn Node>>,
        level: usize,
    ) -> Arc<dyn Node> {
        match nodes_by_level.get_mut(&level) {
            Some(node) => node.clone(),
            None => {
                let node = if level > 0 {
                    let child = get_empty_ternary_tree_locked::<C>(nodes_by_level, level - 1);
                    TernaryNode::<C>::new(level, [child.clone(), child.clone(), child.clone()])
                } else {
                    Arc::new(Leaf::default())
                };
                nodes_by_level.insert(level, node.clone());
                node
            }
        }
    }

    fn get_empty_ternary_tree<C: PoseidonConfig<Scalar, 5> + 'static>(
        level: usize,
    ) -> Arc<dyn Node> {
        static NODES: LazyLock<Mutex<BTreeMap<TypeId, BTreeMap<usize, Arc<dyn Node>>>>> =
            LazyLock::new(|| Mutex::new(BTreeMap::default()));
        let mut nodes_by_config = NODES.lock().unwrap();
        let nodes_by_level = nodes_by_config.entry(TypeId::of::<C>()).or_default();
        get_empty_ternary_tree_locked::<C>(nodes_by_level, level)
    }

    fn test_tall_binary_smt_with_hasher<
        const W: usize,
        const H: usize,
        C: PoseidonConfig<Scalar, 4> + 'static,
    >(
        entries: impl IntoIterator<Item = (u64, u64)>,
        key: u64,
    ) -> Result<()> {
        let tree = {
            let mut tree = get_empty_binary_tree::<C>(H);
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

        let chip = BinaryChip::<Scalar, H, C>::new(*DST, path);
        assert_eq!(chip.width(), W);
        assert_eq!(chip.height(), H);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: BLOWUP_LOG2,
                grinding_bits: 0,
                blind: true,
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

        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, &circuit.default_proof_options())?;
        let openings = circuit.verify(&proof)?;
        assert_eq!(openings[&root_hash], expected_root_hash);

        Ok(())
    }

    fn test_tall_binary_smt_impl<const H: usize>(
        entries: impl IntoIterator<Item = (u64, u64)>,
        key: u64,
    ) -> Result<()> {
        test_tall_binary_smt_with_hasher::<103, H, BlueSkyConfigT4X5>(entries, key)
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

    #[test]
    fn test_tall_binary_smt_with_x7() {
        let entries = [(34, 56), (78, 12)];
        assert!(test_tall_binary_smt_with_hasher::<93, 20, BlueSkyConfigT4X7>(entries, 33).is_ok());
        assert!(test_tall_binary_smt_with_hasher::<93, 20, BlueSkyConfigT4X7>(entries, 34).is_ok());
        assert!(test_tall_binary_smt_with_hasher::<93, 20, BlueSkyConfigT4X7>(entries, 78).is_ok());
    }

    fn test_tall_ternary_smt_with_hasher<
        const W: usize,
        const H: usize,
        C: PoseidonConfig<Scalar, 5> + 'static,
    >(
        entries: impl IntoIterator<Item = (u64, u64)>,
        key: u64,
    ) -> Result<()> {
        let tree = {
            let mut tree = get_empty_ternary_tree::<C>(H);
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

        let chip = TernaryChip::<Scalar, H, C>::new(*DST, path);
        assert_eq!(chip.width(), W);
        assert_eq!(chip.height(), H);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: BLOWUP_LOG2,
                grinding_bits: 0,
                blind: true,
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

        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, &circuit.default_proof_options())?;
        let openings = circuit.verify(&proof)?;
        assert_eq!(openings[&root_hash], expected_root_hash);

        Ok(())
    }

    fn test_tall_ternary_smt_impl<const H: usize>(
        entries: impl IntoIterator<Item = (u64, u64)>,
        key: u64,
    ) -> Result<()> {
        test_tall_ternary_smt_with_hasher::<115, H, BlueSkyConfigT5X5>(entries, key)
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

    #[test]
    fn test_tall_ternary_smt_with_x7() {
        let entries = [(34, 56), (78, 12)];
        assert!(
            test_tall_ternary_smt_with_hasher::<105, 20, BlueSkyConfigT5X7>(entries, 33).is_ok()
        );
        assert!(
            test_tall_ternary_smt_with_hasher::<105, 20, BlueSkyConfigT5X7>(entries, 34).is_ok()
        );
        assert!(
            test_tall_ternary_smt_with_hasher::<105, 20, BlueSkyConfigT5X7>(entries, 78).is_ok()
        );
    }

    #[test]
    #[ignore]
    fn test_taller_ternary_smt_with_large_blowup() {
        let tree = {
            let mut tree = get_empty_ternary_tree::<BlueSkyConfigT5X7>(102);
            tree = tree.put(from_const(12), from_const(34));
            tree = tree.put(from_const(34), from_const(56));
            tree = tree.put(from_const(56), from_const(78));
            tree = tree.put(from_const(78), from_const(90));
            tree = tree.put(from_const(90), from_const(12));
            tree
        };
        let key = from_const(56);
        let value = from_const(78);
        let path: [[Scalar; 3]; 102] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = tree.hash();

        let chip = TernaryChip::<Scalar, 102, BlueSkyConfigT5X7>::new(*DST, path);
        assert_eq!(chip.width(), 105);
        assert_eq!(chip.height(), 102);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs).unwrap();
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: 8,
                grinding_bits: 16,
                blind: true,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), chip.height() + 1);
        assert_eq!(circuit.pcs_config().degree_bound, 128);
        assert_eq!(circuit.num_columns(), chip.width());

        let mut witness = circuit.make_witness();
        let inputs = [witness.cell(0, 0), witness.cell(0, 1)];
        witness.set(inputs[0], key);
        witness.set(inputs[1], value);
        let [root_hash] = witness
            .sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))
            .unwrap();
        let root_hash = match root_hash {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        };

        circuit.check_witness(&witness).unwrap();

        let proof = circuit
            .prove::<Keccak256Hash<Scalar>>(
                witness,
                &ProofOptions {
                    max_fri_rounds: 0,
                    ..circuit.default_proof_options()
                },
            )
            .unwrap();
        let openings = circuit.verify(&proof).unwrap();
        assert_eq!(openings[&root_hash], expected_root_hash);
    }

    fn test_full_binary_smt_impl<I: IntoIterator<Item = (u64, u64)>>(
        entries: I,
        key: u64,
    ) -> Result<()> {
        let tree = {
            let mut tree = get_empty_binary_tree::<BlueSkyConfigT4X5>(256);
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

        let chip = FullBinaryChip::<Scalar, BlueSkyConfigT4X5>::new(*DST, path);
        assert_eq!(chip.width(), 104);
        assert_eq!(chip.height(), 256);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
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
        assert_eq!(circuit.pcs_config().degree_bound, 512);
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

        let proof = circuit.prove::<Sha2Hash<Scalar>>(
            witness,
            &ProofOptions {
                max_fri_rounds: 0,
                ..circuit.default_proof_options()
            },
        )?;
        let openings = circuit.verify(&proof)?;
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
            let mut tree = get_empty_ternary_tree::<BlueSkyConfigT5X5>(162);
            for (key, value) in entries {
                tree = tree.put(key.into(), value.into());
            }
            tree
        };
        let key = key.into();
        let value = tree.get(key);
        let path: [[Scalar; 3]; 162] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = tree.hash();

        let chip = FullTernaryChip::<Scalar, BlueSkyConfigT5X5>::new(*DST, path);
        assert_eq!(chip.width(), 116);
        assert_eq!(chip.height(), 162);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs)?;
        builder.declare_public_cells([root_hash.unwrap()]);
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
        assert_eq!(circuit.pcs_config().degree_bound, 256);
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

        let proof = circuit.prove::<Sha2Hash<Scalar>>(
            witness,
            &ProofOptions {
                max_fri_rounds: 0,
                ..circuit.default_proof_options()
            },
        )?;
        let openings = circuit.verify(&proof)?;
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

    #[test]
    #[ignore]
    fn test_full_ternary_smt_with_large_blowup() {
        let tree = {
            let mut tree = get_empty_ternary_tree::<BlueSkyConfigT5X7>(162);
            tree = tree.put(from_const(12), from_const(34));
            tree = tree.put(from_const(34), from_const(56));
            tree = tree.put(from_const(56), from_const(78));
            tree = tree.put(from_const(78), from_const(90));
            tree = tree.put(from_const(90), from_const(12));
            tree
        };
        let key = from_const(56);
        let value = from_const(78);
        let path: [[Scalar; 3]; 162] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let expected_root_hash = tree.hash();

        let chip = FullTernaryChip::<Scalar, BlueSkyConfigT5X7>::new(*DST, path);
        assert_eq!(chip.width(), 106);
        assert_eq!(chip.height(), 162);

        let mut builder = CircuitBuilder::default();
        let inputs = [builder.cell(0, 0).into(), builder.cell(0, 1).into()];
        let [root_hash] = builder.sub_chip(1, 0, &chip, inputs).unwrap();
        builder.declare_public_cells([root_hash.unwrap()]);
        let circuit = builder
            .build(CircuitOptions {
                canonicalize_constraints: false,
                security_bits: SECURITY_BITS,
                blowup_log2: 8,
                grinding_bits: 16,
                blind: true,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), chip.height() + 1);
        assert_eq!(circuit.pcs_config().degree_bound, 256);
        assert_eq!(circuit.num_columns(), chip.width());

        let mut witness = circuit.make_witness();
        let inputs = [witness.cell(0, 0), witness.cell(0, 1)];
        witness.set(inputs[0], key);
        witness.set(inputs[1], value);
        let [root_hash] = witness
            .sub_chip(1, 0, &chip, inputs.map(CellOrUnconstrained::Cell))
            .unwrap();
        let root_hash = match root_hash {
            CellOrUnconstrained::Cell(cell) => cell,
            _ => panic!(),
        };

        circuit.check_witness(&witness).unwrap();

        let proof = circuit
            .prove::<Keccak256Hash<Scalar>>(
                witness,
                &ProofOptions {
                    max_fri_rounds: 0,
                    ..circuit.default_proof_options()
                },
            )
            .unwrap();
        let openings = circuit.verify(&proof).unwrap();
        assert_eq!(openings[&root_hash], expected_root_hash);
    }
}
