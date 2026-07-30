use anyhow::Result;
use primitive_types::U256;
use starkom_bluesky::{Scalar, from_const};
use starkom_ff::{Field, Field256, PrimeField};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, make_const,
    rvar, var,
};

/// Returns the smallest power of three that is >= n (returns 1 for n=0).
pub fn next_power_of_three(n: usize) -> usize {
    let mut pow = 1usize;
    while pow < n {
        pow *= 3;
    }
    pow
}

/// Checks if a number is a power of 3.
pub fn is_power_of_three(mut value: usize) -> bool {
    if value == 0 {
        return false;
    }
    while value > 1 {
        if value % 3 != 0 {
            return false;
        }
        value /= 3;
    }
    true
}

/// Computes the integer base-3 logarithm of `n`. For example, `ilog3(9) == 2`.
///
/// If `n` is not a power of 3 this function returns the logarithm rounded down to the nearest
/// integer, eg. `ilog3(8) == 1`.
pub fn ilog3(mut n: usize) -> usize {
    let mut c = 0;
    while n >= 3 {
        c += 1;
        n /= 3;
    }
    c
}

/// Returns the LSB of a scalar as a scalar.
pub fn and1(value: Scalar) -> Scalar {
    let lsb = value.to_le_bytes()[0];
    Scalar::from((lsb & 1) as u64)
}

/// Shifts the input [`Scalar`] to the right by `count` bits.
///
/// This is equivalent to the integer division by `2^count`.
pub fn shr(value: Scalar, count: usize) -> Scalar {
    (value.to_u256() >> U256::from(count)).try_into().unwrap()
}

/// Shifts the input [`Scalar`] to the right.
///
/// This is equivalent to the integer division by 2.
pub fn shr1(value: Scalar) -> Scalar {
    shr(value, 1)
}

/// Decomposes the input [`U256`] into `N` bits.
///
/// The `N` bits are represented as scalars and returned in little-endian order.
pub fn decompose_bits<const N: usize>(mut value: U256) -> [Scalar; N] {
    let mut bits = [Scalar::ZERO; N];
    for i in 0..N {
        bits[i] = if value & 1.into() != U256::zero() {
            Scalar::ONE
        } else {
            Scalar::ZERO
        };
        value >>= 1;
    }
    assert_eq!(value, U256::zero());
    bits
}

/// Decomposes the input [`Scalar`] into `N` bits.
///
/// The `N` bits are represented as scalars and returned in little-endian order.
pub fn decompose_scalar_bits<const N: usize>(value: Scalar) -> [Scalar; N] {
    decompose_bits::<N>(value.to_u256())
}

/// Decomposes the input signal into N bits.
///
/// The returned bits are in little-endian order.
///
/// WARNING: this chip is unsafe to use with 255 or 256 bits because it doesn't guard against
/// aliasing. Use the [`FullBitDecomposerChip`] for a full decomposition into 256 bits.
#[derive(Debug, Default, Clone)]
pub struct BitDecomposerChip<const N: usize> {}

impl<const N: usize> PlonkChip<1, N> for BitDecomposerChip<N> {
    fn width(&self) -> usize {
        N + 1
    }

    fn height(&self) -> usize {
        1
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; N]> {
        for i in 0..N {
            view.add_gate(0, var(i) * (make_const(1) - var(i)));
        }
        view.connect(inputs[0], view.cell(0, N).into());
        const TWO: Scalar = from_const(2);
        view.add_gate(
            0,
            var(N)
                - (0..N)
                    .map(|i| var(i) * TWO.pow_small(i))
                    .sum::<Constraint>(),
        );
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 1],
    ) -> Result<[CellOrUnconstrained; N]> {
        let value = match inputs[0] {
            CellOrUnconstrained::Cell(cell) => view.get(cell),
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        decompose_scalar_bits::<N>(value)
            .into_iter()
            .enumerate()
            .for_each(|(i, bit)| view.set(view.cell(0, i), bit));
        view.copy(inputs[0], view.cell(0, N));
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }
}

/// Compares the number represented by the input bits against a specified constant scalar.
///
/// The inputs bits must be provided in little-endian order.
///
/// The returned signal is:
///
///  * -1 if the input value is strictly less than the constant,
///  * 0 if the input value is equal to the constant,
///  * 1 if the input value is strictly greater than the constant.
#[derive(Debug, Default, Clone)]
pub struct ConstBitComparatorChip<const N: usize> {
    rhs: U256,
}

impl<const N: usize> ConstBitComparatorChip<N> {
    pub fn new(rhs: U256) -> Self {
        Self { rhs }
    }
}

impl<const N: usize> ConstBitComparatorChip<N> {
    fn get_rhs_bit(&self, i: usize) -> Scalar {
        ((self.rhs >> i) & 1.into()).try_into().unwrap()
    }
}

impl<const N: usize> PlonkChip<N, 1> for ConstBitComparatorChip<N> {
    fn width(&self) -> usize {
        N
    }

    fn height(&self) -> usize {
        2
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; N],
    ) -> Result<[Option<Cell>; 1]> {
        for i in 0..N {
            view.connect(inputs[i], view.cell(0, i).into());
        }
        view.add_gate(0, rvar(N - 1, 0) - self.get_rhs_bit(N - 1) - rvar(N - 1, 1));
        for i in (0..(N - 1)).rev() {
            let bit = self.get_rhs_bit(i);
            view.add_gate(
                0,
                (rvar(i + 1, 1) ^ 3) + (make_const(1) - (rvar(i + 1, 1) ^ 2)) * (rvar(i, 0) - bit)
                    - rvar(i, 1),
            );
        }
        Ok([view.cell(1, 0).into()])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; N],
    ) -> Result<[CellOrUnconstrained; 1]> {
        for i in 0..N {
            view.copy(inputs[i], view.cell(0, i));
        }
        view.set(
            view.cell(1, N - 1),
            view.get(view.cell(0, N - 1)) - self.get_rhs_bit(N - 1),
        );
        for i in (0..(N - 1)).rev() {
            let bit = self.get_rhs_bit(i);
            let cmp = view.get(view.cell(0, i)) - bit;
            let prev = view.get(view.cell(1, i + 1));
            view.set(
                view.cell(1, i),
                prev.cube() + (from_const(1) - prev.square()) * cmp,
            );
        }
        Ok([view.cell(1, 0).into()])
    }
}

/// Decomposes an input signal into 256 bits.
///
/// The returned bits are in little-endian order.
///
/// Note that the MSB will always be zero because BlueSky scalars don't cover the upper half of the
/// 256 bit range.
#[derive(Debug, Clone)]
pub struct FullBitDecomposerChip {
    decomposer: BitDecomposerChip<256>,
    comparator: ConstBitComparatorChip<256>,
}

impl Default for FullBitDecomposerChip {
    fn default() -> Self {
        Self {
            decomposer: BitDecomposerChip::default(),
            comparator: ConstBitComparatorChip::new(Scalar::MODULUS.parse().unwrap()),
        }
    }
}

impl PlonkChip<1, 256> for FullBitDecomposerChip {
    fn width(&self) -> usize {
        257
    }

    fn height(&self) -> usize {
        3
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 256]> {
        let bits = self.decomposer.build(view, inputs)?;
        let [cmp] = view.sub_chip(1, 0, &self.comparator, bits)?;
        view.add_gate(cmp.unwrap().row(), var(0) + 1);
        Ok(bits)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 1],
    ) -> Result<[CellOrUnconstrained; 256]> {
        let bits = self.decomposer.witness(view, inputs)?;
        view.sub_chip(1, 0, &self.comparator, bits)?;
        Ok(bits)
    }
}

/// Divides `value`, treated as an integer, by `3^exp`, rounding down.
pub fn div_pow3(value: Scalar, exp: usize) -> Scalar {
    let dividend = value.to_u256();
    let divisor = U256::from(3).pow(exp.into());
    (dividend / divisor).try_into().unwrap()
}

/// Divides `value`, treated as an integer, by 3, rounding down.
pub fn div3(value: Scalar) -> Scalar {
    let dividend = value.to_u256();
    (dividend / 3).try_into().unwrap()
}

/// Returns `value`, treated as an integer, modulo 3.
pub fn mod3(value: Scalar) -> Scalar {
    let value = value.to_u256();
    (value % 3).try_into().unwrap()
}

/// Decomposes `value` into its `N` base-3 digits (trits), least significant first.
///
/// Panics if `value` does not fit in `N` trits.
pub fn decompose_trits<const N: usize>(mut value: U256) -> [Scalar; N] {
    let mut trits = [Scalar::ZERO; N];
    for i in 0..N {
        trits[i] = Scalar::from((value % 3).as_u64());
        value /= 3;
    }
    assert_eq!(value, U256::zero());
    trits
}

/// Like [`decompose_trits`], but takes `value` as a [`Scalar`] rather than a [`U256`].
pub fn decompose_scalar_trits<const N: usize>(value: Scalar) -> [Scalar; N] {
    decompose_trits::<N>(value.to_u256())
}

/// Decomposes the input signal into N trits.
///
/// WARNING: this chip is unsafe to use with 160 or 161 trits because it doesn't guard against
/// aliasing. Use the [`FullTritDecomposerChip`] for a full decomposition into 161 trits.
#[derive(Debug, Default, Clone)]
pub struct TritDecomposerChip<const N: usize> {}

impl<const N: usize> PlonkChip<1, N> for TritDecomposerChip<N> {
    fn width(&self) -> usize {
        N + 1
    }

    fn height(&self) -> usize {
        1
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; N]> {
        for i in 0..N {
            view.add_gate(0, var(i) * (var(i) - 1) * (var(i) - 2));
        }
        view.connect(inputs[0], view.cell(0, N).into());
        const THREE: Scalar = from_const(3);
        view.add_gate(
            0,
            var(N)
                - (0..N)
                    .map(|i| var(i) * THREE.pow_small(i))
                    .sum::<Constraint>(),
        );
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; 1],
    ) -> Result<[CellOrUnconstrained; N]> {
        let value = match inputs[0] {
            CellOrUnconstrained::Cell(cell) => view.get(cell),
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        decompose_scalar_trits::<N>(value)
            .into_iter()
            .enumerate()
            .for_each(|(i, trit)| view.set(view.cell(0, i), trit));
        view.copy(inputs[0], view.cell(0, N));
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }
}

/// Compares the number represented by the input trits against a specified constant scalar.
///
/// The returned signal is:
///
///  * -1 if the input value is strictly less than the constant,
///  * 0 if the input value is equal to the constant,
///  * 1 if the input value is strictly greater than the constant.
#[derive(Debug, Default, Clone)]
pub struct ConstTritComparatorChip<const N: usize> {
    rhs: U256,
}

impl<const N: usize> ConstTritComparatorChip<N> {
    pub fn new(rhs: U256) -> Self {
        Self { rhs }
    }

    fn get_rhs_trit(&self, i: usize) -> Scalar {
        let three = U256::from(3);
        ((self.rhs / three.pow(i.into())) % three)
            .try_into()
            .unwrap()
    }
}

impl<const N: usize> PlonkChip<N, 1> for ConstTritComparatorChip<N> {
    fn width(&self) -> usize {
        N
    }

    fn height(&self) -> usize {
        2
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; N],
    ) -> Result<[Option<Cell>; 1]> {
        // TODO
        todo!()
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; N],
    ) -> Result<[CellOrUnconstrained; 1]> {
        // TODO
        todo!()
    }
}

// TODO

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::parse_scalar;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use std::cmp::Ordering;

    const BLOWUP_LOG2: usize = 1;

    #[inline]
    fn cell(row: usize, column: usize) -> Cell {
        Cell::new(row, column)
    }

    #[test]
    fn test_next_power_of_three() {
        assert_eq!(next_power_of_three(0), 1);
        assert_eq!(next_power_of_three(1), 1);
        assert_eq!(next_power_of_three(2), 3);
        assert_eq!(next_power_of_three(3), 3);
        assert_eq!(next_power_of_three(4), 9);
        assert_eq!(next_power_of_three(5), 9);
        assert_eq!(next_power_of_three(6), 9);
        assert_eq!(next_power_of_three(7), 9);
        assert_eq!(next_power_of_three(8), 9);
        assert_eq!(next_power_of_three(9), 9);
        assert_eq!(next_power_of_three(10), 27);
        assert_eq!(next_power_of_three(11), 27);
    }

    #[test]
    fn test_is_power_of_three() {
        assert!(!is_power_of_three(0));
        assert!(is_power_of_three(1));
        assert!(!is_power_of_three(2));
        assert!(is_power_of_three(3));
        assert!(!is_power_of_three(4));
        assert!(!is_power_of_three(5));
        assert!(!is_power_of_three(6));
        assert!(!is_power_of_three(7));
        assert!(!is_power_of_three(8));
        assert!(is_power_of_three(9));
        assert!(!is_power_of_three(10));
        assert!(!is_power_of_three(11));
    }

    #[test]
    fn test_ilog3() {
        assert_eq!(ilog3(0), 0);
        assert_eq!(ilog3(1), 0);
        assert_eq!(ilog3(2), 0);
        assert_eq!(ilog3(3), 1);
        assert_eq!(ilog3(4), 1);
        assert_eq!(ilog3(5), 1);
        assert_eq!(ilog3(6), 1);
        assert_eq!(ilog3(7), 1);
        assert_eq!(ilog3(8), 1);
        assert_eq!(ilog3(9), 2);
        assert_eq!(ilog3(10), 2);
        assert_eq!(ilog3(11), 2);
    }

    #[test]
    fn test_and1() {
        assert_eq!(and1(from_const(42)), from_const(0));
        assert_eq!(and1(from_const(43)), from_const(1));
        assert_eq!(and1(from_const(44)), from_const(0));
        assert_eq!(and1(from_const(45)), from_const(1));
    }

    #[test]
    fn test_and1_large() {
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            from_const(0)
        );
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f21"
            )),
            from_const(1)
        );
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f22"
            )),
            from_const(0)
        );
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f23"
            )),
            from_const(1)
        );
    }

    #[test]
    fn test_shr() {
        assert_eq!(
            shr(
                parse_scalar("0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"),
                4
            ),
            parse_scalar("0x00102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f2")
        );
    }

    #[test]
    fn test_shr1() {
        assert_eq!(
            shr1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            parse_scalar("0x008101820283038404850586068707880889098a0a8b0b8c0c8d0d8e0e8f0f90")
        );
    }

    #[test]
    fn test_decompose_bits_one() {
        assert_eq!(decompose_bits::<1>(0.into()), [from_const(0)]);
        assert_eq!(decompose_bits::<1>(1.into()), [from_const(1)]);
    }

    #[test]
    fn test_decompose_bits_two() {
        assert_eq!(
            decompose_bits::<2>(0.into()),
            [from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_bits::<2>(1.into()),
            [from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_bits::<2>(2.into()),
            [from_const(0), from_const(1)]
        );
        assert_eq!(
            decompose_bits::<2>(3.into()),
            [from_const(1), from_const(1)]
        );
    }

    #[test]
    fn test_decompose_bits_three() {
        assert_eq!(
            decompose_bits::<3>(0.into()),
            [from_const(0), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_bits::<3>(1.into()),
            [from_const(1), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_bits::<3>(2.into()),
            [from_const(0), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_bits::<3>(3.into()),
            [from_const(1), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_bits::<3>(4.into()),
            [from_const(0), from_const(0), from_const(1)]
        );
        assert_eq!(
            decompose_bits::<3>(5.into()),
            [from_const(1), from_const(0), from_const(1)]
        );
        assert_eq!(
            decompose_bits::<3>(6.into()),
            [from_const(0), from_const(1), from_const(1)]
        );
        assert_eq!(
            decompose_bits::<3>(7.into()),
            [from_const(1), from_const(1), from_const(1)]
        );
    }

    #[test]
    fn test_decompose_bits_large() {
        assert_eq!(
            decompose_bits::<64>(0xFFFFFFFFFFFFFFFFu64.into()),
            [from_const(1); 64]
        );
    }

    #[test]
    fn test_decompose_scalar_bits() {
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(0)),
            [from_const(0), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(1)),
            [from_const(1), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(2)),
            [from_const(0), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(3)),
            [from_const(1), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(4)),
            [from_const(0), from_const(0), from_const(1)]
        );
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(5)),
            [from_const(1), from_const(0), from_const(1)]
        );
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(6)),
            [from_const(0), from_const(1), from_const(1)]
        );
        assert_eq!(
            decompose_scalar_bits::<3>(from_const(7)),
            [from_const(1), from_const(1), from_const(1)]
        );
    }

    fn test_bit_decomposer_chip<const N: usize>(value: u64) {
        let chip = BitDecomposerChip::<N>::default();
        assert_eq!(chip.width(), N + 1);
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::default();
        chip.build(&mut builder, [None]).unwrap();
        builder.declare_public_rows([0]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.degree_bound(), 4);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let bits = chip
            .witness(&mut witness, [Scalar::from(value).into()])
            .unwrap()
            .map(|bit| match bit {
                CellOrUnconstrained::Cell(cell) => witness.get(cell),
                _ => panic!("the output bits must be constrained"),
            });
        assert_eq!(bits, decompose_bits::<N>(value.into())[0..N]);
        circuit.check_witness(&witness).unwrap();
        let proving_options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, proving_options.clone())
            .unwrap();
        let openings = circuit
            .to_compressed::<Sha2Hash<Scalar>>(proving_options)
            .verify(&proof)
            .unwrap();
        assert!((0..N).all(|i| openings[&cell(0, i)] == bits[i]));
    }

    #[test]
    fn test_bit_decomposer_chip_1() {
        test_bit_decomposer_chip::<1>(0);
        test_bit_decomposer_chip::<1>(1);
    }

    #[test]
    fn test_bit_decomposer_chip_2() {
        test_bit_decomposer_chip::<2>(0);
        test_bit_decomposer_chip::<2>(1);
        test_bit_decomposer_chip::<2>(2);
        test_bit_decomposer_chip::<2>(3);
    }

    #[test]
    fn test_bit_decomposer_chip_3() {
        test_bit_decomposer_chip::<3>(0);
        test_bit_decomposer_chip::<3>(1);
        test_bit_decomposer_chip::<3>(2);
        test_bit_decomposer_chip::<3>(3);
        test_bit_decomposer_chip::<3>(4);
        test_bit_decomposer_chip::<3>(5);
        test_bit_decomposer_chip::<3>(6);
        test_bit_decomposer_chip::<3>(7);
    }

    fn test_const_bit_comparator_chip<const N: usize>(lhs: u64, rhs: u64) {
        let mut builder = CircuitBuilder::default();
        let decomposer_chip = BitDecomposerChip::<N>::default();
        let bits = decomposer_chip.build(&mut builder, [None]).unwrap();
        let comparator_chip = ConstBitComparatorChip::<N>::new(rhs.into());
        assert_eq!(comparator_chip.width(), N);
        assert_eq!(comparator_chip.height(), 2);
        let [cmp] = builder
            .sub_chip(decomposer_chip.height(), 0, &comparator_chip, bits)
            .unwrap();
        builder.declare_public_rows([cmp.unwrap().row()]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(
            circuit.num_rows(),
            decomposer_chip.height() + comparator_chip.height()
        );
        assert_eq!(circuit.degree_bound(), 8);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let bits = decomposer_chip
            .witness(&mut witness, [Scalar::from(lhs).into()])
            .unwrap();
        assert!(
            witness
                .sub_chip(decomposer_chip.height(), 0, &comparator_chip, bits)
                .is_ok()
        );
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, options.clone())
            .unwrap();
        let openings = circuit
            .to_compressed::<Sha2Hash<Scalar>>(options)
            .verify(&proof)
            .unwrap();
        assert_eq!(
            openings[&cmp.unwrap()],
            match lhs.cmp(&rhs) {
                Ordering::Less => -from_const(1),
                Ordering::Equal => from_const(0),
                Ordering::Greater => from_const(1),
            }
        );
    }

    #[test]
    fn test_const_bit_comparator_chip_1() {
        test_const_bit_comparator_chip::<1>(0, 0);
        test_const_bit_comparator_chip::<1>(1, 0);
        test_const_bit_comparator_chip::<1>(0, 1);
        test_const_bit_comparator_chip::<1>(1, 1);
    }

    #[test]
    fn test_const_bit_comparator_chip_2() {
        test_const_bit_comparator_chip::<2>(0, 0);
        test_const_bit_comparator_chip::<2>(1, 0);
        test_const_bit_comparator_chip::<2>(2, 0);
        test_const_bit_comparator_chip::<2>(3, 0);
        test_const_bit_comparator_chip::<2>(0, 1);
        test_const_bit_comparator_chip::<2>(1, 1);
        test_const_bit_comparator_chip::<2>(2, 1);
        test_const_bit_comparator_chip::<2>(3, 1);
        test_const_bit_comparator_chip::<2>(0, 2);
        test_const_bit_comparator_chip::<2>(1, 2);
        test_const_bit_comparator_chip::<2>(2, 2);
        test_const_bit_comparator_chip::<2>(3, 2);
        test_const_bit_comparator_chip::<2>(0, 3);
        test_const_bit_comparator_chip::<2>(1, 3);
        test_const_bit_comparator_chip::<2>(2, 3);
        test_const_bit_comparator_chip::<2>(3, 3);
    }

    fn test_full_bit_decomposer_chip_impl(value: u64) {
        let chip = FullBitDecomposerChip::default();
        assert_eq!(chip.width(), 257);
        assert_eq!(chip.height(), 3);
        let mut builder = CircuitBuilder::default();
        chip.build(&mut builder, [None]).unwrap();
        builder.declare_public_rows([0]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), 3);
        assert_eq!(circuit.degree_bound(), 8);
        assert_eq!(circuit.num_columns(), 257);
        let mut witness = circuit.make_witness();
        let bits = chip
            .witness(&mut witness, [Scalar::from(value).into()])
            .unwrap()
            .map(|bit| match bit {
                CellOrUnconstrained::Cell(cell) => witness.get(cell),
                _ => panic!("the output bits must be constrained"),
            });
        assert_eq!(bits, decompose_bits::<256>(value.into())[0..256]);
        circuit.check_witness(&witness).unwrap();
        let proving_options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, proving_options.clone())
            .unwrap();
        let openings = circuit
            .to_compressed::<Sha2Hash<Scalar>>(proving_options)
            .verify(&proof)
            .unwrap();
        assert!((0..256).all(|i| openings[&cell(0, i)] == bits[i]));
    }

    #[test]
    fn test_full_bit_decomposer_chip_0() {
        test_full_bit_decomposer_chip_impl(0);
    }

    #[test]
    fn test_full_bit_decomposer_chip_1() {
        test_full_bit_decomposer_chip_impl(1);
    }

    #[test]
    fn test_full_bit_decomposer_chip_2() {
        test_full_bit_decomposer_chip_impl(2);
    }

    #[test]
    fn test_full_bit_decomposer_chip_3() {
        test_full_bit_decomposer_chip_impl(3);
    }

    #[test]
    fn test_full_bit_decomposer_chip_4() {
        test_full_bit_decomposer_chip_impl(4);
    }

    #[test]
    fn test_full_bit_decomposer_chip_5() {
        test_full_bit_decomposer_chip_impl(5);
    }

    #[test]
    fn test_full_bit_decomposer_chip_6() {
        test_full_bit_decomposer_chip_impl(6);
    }

    #[test]
    fn test_full_bit_decomposer_chip_7() {
        test_full_bit_decomposer_chip_impl(7);
    }

    #[test]
    fn test_div_pow3() {
        assert_eq!(
            div_pow3(
                parse_scalar("0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"),
                4
            ),
            parse_scalar("0x00032f71d3d0aac0e3aaca6871f05f0032c75591a1720ced55a4ab0058da7229")
        );
    }

    #[test]
    fn test_div3() {
        assert_eq!(
            div3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            parse_scalar("0x005601015702025803035904045a05055b06065c07075d08085e09095f0a0a60")
        );
    }

    #[test]
    fn test_mod3() {
        assert_eq!(mod3(from_const(42)), from_const(0));
        assert_eq!(mod3(from_const(43)), from_const(1));
        assert_eq!(mod3(from_const(44)), from_const(2));
        assert_eq!(mod3(from_const(45)), from_const(0));
        assert_eq!(mod3(from_const(46)), from_const(1));
        assert_eq!(mod3(from_const(47)), from_const(2));
    }

    #[test]
    fn test_mod3_large() {
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            from_const(0)
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f21"
            )),
            from_const(1)
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f22"
            )),
            from_const(2)
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f23"
            )),
            from_const(0)
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f24"
            )),
            from_const(1)
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f25"
            )),
            from_const(2)
        );
    }

    #[test]
    fn test_decompose_trits_one() {
        assert_eq!(decompose_trits::<1>(0.into()), [from_const(0)]);
        assert_eq!(decompose_trits::<1>(1.into()), [from_const(1)]);
        assert_eq!(decompose_trits::<1>(2.into()), [from_const(2)]);
    }

    #[test]
    fn test_decompose_trits_two() {
        assert_eq!(
            decompose_trits::<2>(0.into()),
            [from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<2>(1.into()),
            [from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<2>(2.into()),
            [from_const(2), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<2>(3.into()),
            [from_const(0), from_const(1)]
        );
        assert_eq!(
            decompose_trits::<2>(4.into()),
            [from_const(1), from_const(1)]
        );
        assert_eq!(
            decompose_trits::<2>(5.into()),
            [from_const(2), from_const(1)]
        );
        assert_eq!(
            decompose_trits::<2>(6.into()),
            [from_const(0), from_const(2)]
        );
        assert_eq!(
            decompose_trits::<2>(7.into()),
            [from_const(1), from_const(2)]
        );
        assert_eq!(
            decompose_trits::<2>(8.into()),
            [from_const(2), from_const(2)]
        );
    }

    #[test]
    fn test_decompose_trits_three() {
        assert_eq!(
            decompose_trits::<3>(0.into()),
            [from_const(0), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(1.into()),
            [from_const(1), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(2.into()),
            [from_const(2), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(3.into()),
            [from_const(0), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(4.into()),
            [from_const(1), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(5.into()),
            [from_const(2), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(6.into()),
            [from_const(0), from_const(2), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(7.into()),
            [from_const(1), from_const(2), from_const(0)]
        );
        assert_eq!(
            decompose_trits::<3>(8.into()),
            [from_const(2), from_const(2), from_const(0)]
        );
    }

    #[test]
    fn test_decompose_scalar_trits() {
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(0)),
            [from_const(0), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(1)),
            [from_const(1), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(2)),
            [from_const(2), from_const(0), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(3)),
            [from_const(0), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(4)),
            [from_const(1), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(5)),
            [from_const(2), from_const(1), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(6)),
            [from_const(0), from_const(2), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(7)),
            [from_const(1), from_const(2), from_const(0)]
        );
        assert_eq!(
            decompose_scalar_trits::<3>(from_const(8)),
            [from_const(2), from_const(2), from_const(0)]
        );
    }

    fn test_trit_decomposer_chip<const N: usize>(value: u64) {
        let chip = TritDecomposerChip::<N>::default();
        assert_eq!(chip.width(), N + 1);
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::default();
        chip.build(&mut builder, [None]).unwrap();
        builder.declare_public_rows([0]);
        let circuit = builder
            .build(CompilationOptions {
                canonicalize_constraints: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.degree_bound(), 4);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let trits = chip
            .witness(&mut witness, [Scalar::from(value).into()])
            .unwrap()
            .map(|trit| match trit {
                CellOrUnconstrained::Cell(cell) => witness.get(cell),
                _ => panic!("the output trits must be constrained"),
            });
        assert_eq!(trits, decompose_trits::<N>(value.into())[0..N]);
        circuit.check_witness(&witness).unwrap();
        let proving_options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, proving_options.clone())
            .unwrap();
        let openings = circuit
            .to_compressed::<Sha2Hash<Scalar>>(proving_options)
            .verify(&proof)
            .unwrap();
        assert!((0..N).all(|i| openings[&cell(0, i)] == trits[i]));
    }

    #[test]
    fn test_trit_decomposer_chip_1() {
        test_trit_decomposer_chip::<1>(0);
        test_trit_decomposer_chip::<1>(1);
        test_trit_decomposer_chip::<1>(2);
    }

    #[test]
    fn test_trit_decomposer_chip_2() {
        test_trit_decomposer_chip::<2>(0);
        test_trit_decomposer_chip::<2>(1);
        test_trit_decomposer_chip::<2>(2);
        test_trit_decomposer_chip::<2>(3);
        test_trit_decomposer_chip::<2>(4);
        test_trit_decomposer_chip::<2>(5);
        test_trit_decomposer_chip::<2>(6);
        test_trit_decomposer_chip::<2>(7);
        test_trit_decomposer_chip::<2>(8);
    }

    #[test]
    fn test_trit_decomposer_chip_3() {
        test_trit_decomposer_chip::<3>(0);
        test_trit_decomposer_chip::<3>(1);
        test_trit_decomposer_chip::<3>(2);
        test_trit_decomposer_chip::<3>(3);
        test_trit_decomposer_chip::<3>(4);
        test_trit_decomposer_chip::<3>(5);
        test_trit_decomposer_chip::<3>(6);
        test_trit_decomposer_chip::<3>(7);
        test_trit_decomposer_chip::<3>(8);
        test_trit_decomposer_chip::<3>(9);
        test_trit_decomposer_chip::<3>(10);
        test_trit_decomposer_chip::<3>(11);
        test_trit_decomposer_chip::<3>(12);
        test_trit_decomposer_chip::<3>(13);
        test_trit_decomposer_chip::<3>(14);
        test_trit_decomposer_chip::<3>(15);
        test_trit_decomposer_chip::<3>(16);
        test_trit_decomposer_chip::<3>(17);
        test_trit_decomposer_chip::<3>(18);
        test_trit_decomposer_chip::<3>(19);
        test_trit_decomposer_chip::<3>(20);
        test_trit_decomposer_chip::<3>(21);
        test_trit_decomposer_chip::<3>(22);
        test_trit_decomposer_chip::<3>(23);
        test_trit_decomposer_chip::<3>(24);
        test_trit_decomposer_chip::<3>(25);
        test_trit_decomposer_chip::<3>(26);
    }

    // TODO
}
