use anyhow::Result;
use primitive_types::U256;
use starkom_bluesky::Scalar;
use starkom_ff::{Field, Field256, PrimeField};
use starkom_plonk::{Chip as PlonkChip, CircuitBuilder, Wire, WireOrUnconstrained, Witness};

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

/// Calculates (1 - X^2).
///
/// Not an actual "logical NOT", it can only NOT the (-1,0,1) signals we use in the comparator
/// chips. `LogicalNotChip` is for internal use by those chips.
#[derive(Debug, Default, Clone)]
struct LogicalNotChip {}

impl PlonkChip<1, 1> for LogicalNotChip {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; 1]> {
        Ok([builder
            .add_unary_gate(
                Scalar::from_const(0),
                Scalar::from_const(0),
                -Scalar::from_const(1),
                -Scalar::from_const(1),
                Scalar::from_const(1),
                inputs[0],
            )
            .into()])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; 1]> {
        let input = inputs[0];
        let gate = witness.pop_gate();
        witness.copy(input.into(), Wire::LeftIn(gate));
        let input = witness.copy(input.into(), Wire::RightIn(gate));
        let out = Wire::Out(gate);
        witness.set(out, Scalar::from_const(1) - input.square());
        Ok([out.into()])
    }
}

pub fn and1(value: Scalar) -> Scalar {
    let lsb = value.to_le_bytes()[0];
    Scalar::from((lsb & 1) as u64)
}

pub fn shr(value: Scalar, count: usize) -> Scalar {
    (value.to_u256() >> U256::from(count)).try_into().unwrap()
}

pub fn shr1(value: Scalar) -> Scalar {
    shr(value, 1)
}

pub fn decompose_bits<const N: usize>(mut value: U256) -> [Scalar; N] {
    let mut bits = [Scalar::ZERO; N];
    for i in 0..N {
        bits[i] = Scalar::from((value & 1.into()).as_u64());
        value >>= 1;
    }
    assert_eq!(value, U256::zero());
    bits
}

pub fn decompose_scalar_bits<const N: usize>(value: Scalar) -> [Scalar; N] {
    decompose_bits::<N>(value.to_u256())
}

/// Decomposes the input signal into N bits.
///
/// WARNING: this chip is unsafe to use with 255 or 256 bits because it doesn't guard against
/// aliasing. Use the `FullBitDecomposerChip` for a full decomposition into 256 bits (note that the
/// MSB will always be zero because BLS12-381 scalars don't cover the upper half of the 256 bit
/// range).
#[derive(Debug, Default, Clone)]
pub struct BitDecomposerChip<const N: usize> {}

impl<const N: usize> PlonkChip<1, N> for BitDecomposerChip<N> {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; N]> {
        let mut sum = builder.add_const_gate(Scalar::ZERO);
        let mut power = Scalar::from_const(1);
        let bits = std::array::from_fn(|_| {
            sum =
                builder.add_linear_combination_gate(Scalar::from_const(1), sum.into(), power, None);
            power = power.double();
            let bit = Some(Wire::RightIn(sum.gate()));
            builder.add_bit_assertion_gate(bit);
            bit
        });
        if let Some(input) = inputs[0] {
            builder.connect(sum, input);
        }
        Ok(bits)
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; N]> {
        let mut input = match inputs[0] {
            WireOrUnconstrained::Wire(wire) => witness.get(wire),
            WireOrUnconstrained::Unconstrained(value) => value,
        };
        let mut sum = witness.assert_constant(Scalar::ZERO);
        let mut power = Scalar::from_const(1);
        let bits = std::array::from_fn(|_| {
            let bit = and1(input);
            input = shr1(input);
            sum = witness.combine(Scalar::from_const(1), sum.into(), power, bit.into());
            power = power.double();
            let bit = Wire::RightIn(sum.gate()).into();
            witness.assert_bit(bit);
            bit
        });
        Ok(bits)
    }
}

/// Compares the number represented by the input bits against a specified constant scalar.
///
/// The returned signal is:
///
///  * -1 if the input value is strictly less than the constant,
///  * 0 if the input value is equal to the constant,
///  * 1 if the input value is strictly greater than the constant.
#[derive(Debug, Clone)]
pub struct BitComparatorChip<const N: usize> {
    rhs: U256,
    logical_not: LogicalNotChip,
}

impl<const N: usize> BitComparatorChip<N> {
    pub fn new(rhs: U256) -> Self {
        Self {
            rhs,
            logical_not: LogicalNotChip::default(),
        }
    }
}

impl<const N: usize> BitComparatorChip<N> {
    fn get_rhs_bit(&self, i: usize) -> Scalar {
        ((self.rhs >> i) & 1.into()).try_into().unwrap()
    }
}

impl<const N: usize> PlonkChip<N, 1> for BitComparatorChip<N> {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; N],
    ) -> Result<[Option<Wire>; 1]> {
        assert!(N > 0);
        let mut cmp = builder.add_sub_const_gate(inputs[N - 1], self.get_rhs_bit(N - 1));
        for i in (0..(N - 1)).rev() {
            let cmp2 = builder.add_sub_const_gate(inputs[i], self.get_rhs_bit(i));
            let not = self.logical_not.build(builder, [cmp.into()])?[0];
            let rhs = builder.add_mul_gate(cmp2.into(), not);
            cmp = builder.add_sum_gate(cmp.into(), rhs.into());
        }
        Ok([Some(cmp)])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; N],
    ) -> Result<[WireOrUnconstrained; 1]> {
        assert!(N > 0);
        let mut cmp = witness.sub_const(inputs[N - 1], self.get_rhs_bit(N - 1));
        for i in (0..(N - 1)).rev() {
            let cmp2 = witness.sub_const(inputs[i], self.get_rhs_bit(i));
            let not = self.logical_not.witness(witness, [cmp.into()])?[0];
            let rhs = witness.mul(cmp2.into(), not);
            cmp = witness.add(cmp.into(), rhs.into());
        }
        Ok([cmp.into()])
    }
}

/// Decomposes an input signal into 256 bits.
#[derive(Debug, Clone)]
pub struct FullBitDecomposerChip {
    decomposer: BitDecomposerChip<256>,
    comparator: BitComparatorChip<256>,
}

impl Default for FullBitDecomposerChip {
    fn default() -> Self {
        Self {
            decomposer: BitDecomposerChip::default(),
            comparator: BitComparatorChip::new(Scalar::MODULUS.parse().unwrap()),
        }
    }
}

impl PlonkChip<1, 256> for FullBitDecomposerChip {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; 256]> {
        let bits = self.decomposer.build(builder, inputs)?;
        let cmp = self.comparator.build(builder, bits)?[0].unwrap();
        let c = builder.add_const_gate(-Scalar::from_const(1));
        builder.connect(cmp, c);
        Ok(bits)
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; 256]> {
        let bits = self.decomposer.witness(witness, inputs)?;
        self.comparator.witness(witness, bits)?;
        witness.assert_constant(-Scalar::from_const(1));
        Ok(bits)
    }
}

pub fn div_pow3(value: Scalar, exp: usize) -> Scalar {
    let dividend = value.to_u256();
    let divisor = U256::from(3).pow(exp.into());
    (dividend / divisor).try_into().unwrap()
}

pub fn div3(value: Scalar) -> Scalar {
    let dividend = value.to_u256();
    (dividend / 3).try_into().unwrap()
}

pub fn mod3(value: Scalar) -> Scalar {
    let value = value.to_u256();
    (value % 3).try_into().unwrap()
}

pub fn decompose_trits<const N: usize>(mut value: U256) -> [Scalar; N] {
    let mut trits = [Scalar::ZERO; N];
    for i in 0..N {
        trits[i] = Scalar::from((value % 3).as_u64());
        value /= 3;
    }
    assert_eq!(value, U256::zero());
    trits
}

pub fn decompose_scalar_trits<const N: usize>(value: Scalar) -> [Scalar; N] {
    decompose_trits::<N>(value.to_u256())
}

/// Decomposes the input signal into N trits.
///
/// WARNING: this chip is unsafe to use with 160 or 161 trits because it doesn't guard against
/// aliasing. Use the `FullTritDecomposerChip` for a full decomposition into 161 trits.
#[derive(Debug, Default, Clone)]
pub struct TritDecomposerChip<const N: usize> {}

impl<const N: usize> PlonkChip<1, N> for TritDecomposerChip<N> {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; N]> {
        let mut sum = builder.add_const_gate(Scalar::ZERO);
        let mut power = Scalar::from_const(1);
        let trits = std::array::from_fn(|_| {
            sum =
                builder.add_linear_combination_gate(Scalar::from_const(1), sum.into(), power, None);
            power = power.double() + power;
            let trit = Some(Wire::RightIn(sum.gate()));
            builder.add_trit_assertion_gate(trit);
            trit
        });
        if let Some(input) = inputs[0] {
            builder.connect(sum, input);
        }
        Ok(trits)
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; N]> {
        let mut input = match inputs[0] {
            WireOrUnconstrained::Wire(wire) => witness.get(wire),
            WireOrUnconstrained::Unconstrained(value) => value,
        };
        let mut sum = witness.assert_constant(Scalar::ZERO);
        let mut power = Scalar::from_const(1);
        let trits = std::array::from_fn(|_| {
            let trit = mod3(input);
            input = div3(input);
            sum = witness.combine(Scalar::from_const(1), sum.into(), power, trit.into());
            power = power.double() + power;
            let trit = Wire::RightIn(sum.gate()).into();
            witness.assert_trit(trit);
            trit
        });
        Ok(trits)
    }
}

/// Compares the number represented by the input trits against a specified constant scalar.
///
/// The returned signal is:
///
///  * -1 if the input value is strictly less than the constant,
///  * 0 if the input value is equal to the constant,
///  * 1 if the input value is strictly greater than the constant.
#[derive(Debug, Clone)]
pub struct TritComparatorChip<const N: usize> {
    rhs: U256,
    logical_not: LogicalNotChip,
}

impl<const N: usize> TritComparatorChip<N> {
    pub fn new(rhs: U256) -> Self {
        Self {
            rhs,
            logical_not: LogicalNotChip::default(),
        }
    }

    fn get_rhs_trit(&self, i: usize) -> Scalar {
        let three = U256::from(3);
        ((self.rhs / three.pow(i.into())) % three)
            .try_into()
            .unwrap()
    }

    fn build_compare_trits(builder: &mut CircuitBuilder, lhs: Option<Wire>, rhs: Scalar) -> Wire {
        let sub = builder.add_sub_const_gate(lhs, rhs);
        let square = builder.add_square_gate(sub.into());
        let cube = builder.add_mul_gate(sub.into(), square.into());
        builder.add_binary_gate(
            -Scalar::from_const(1),
            Scalar::from_const(7),
            -Scalar::from_const(6),
            Scalar::from_const(0),
            Scalar::from_const(0),
            cube.into(),
            sub.into(),
        )
    }

    fn witness_compare_trits(witness: &mut Witness, lhs: WireOrUnconstrained, rhs: Scalar) -> Wire {
        let sub = witness.sub_const(lhs, rhs);
        let square = witness.square(sub.into());
        let cube = witness.mul(sub.into(), square.into());
        let gate = witness.pop_gate();
        let lhs = witness.copy(cube.into(), Wire::LeftIn(gate).into());
        let rhs = witness.copy(sub.into(), Wire::RightIn(gate).into());
        let out = Wire::Out(gate);
        let div6 = Scalar::from_const(6).invert().into_option().unwrap();
        witness.set(out, (-lhs + rhs * Scalar::from_const(7)) * div6);
        out
    }
}

impl<const N: usize> PlonkChip<N, 1> for TritComparatorChip<N> {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; N],
    ) -> Result<[Option<Wire>; 1]> {
        assert!(N > 0);
        let mut cmp = Self::build_compare_trits(builder, inputs[N - 1], self.get_rhs_trit(N - 1));
        for i in (0..(N - 1)).rev() {
            let cmp2 = Self::build_compare_trits(builder, inputs[i], self.get_rhs_trit(i));
            let not = self.logical_not.build(builder, [cmp.into()])?[0];
            let rhs = builder.add_mul_gate(cmp2.into(), not);
            cmp = builder.add_sum_gate(cmp.into(), rhs.into());
        }
        Ok([Some(cmp)])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; N],
    ) -> Result<[WireOrUnconstrained; 1]> {
        assert!(N > 0);
        let mut cmp = Self::witness_compare_trits(witness, inputs[N - 1], self.get_rhs_trit(N - 1));
        for i in (0..(N - 1)).rev() {
            let cmp2 = Self::witness_compare_trits(witness, inputs[i], self.get_rhs_trit(i));
            let not = self.logical_not.witness(witness, [cmp.into()])?[0];
            let rhs = witness.mul(cmp2.into(), not);
            cmp = witness.add(cmp.into(), rhs.into());
        }
        Ok([cmp.into()])
    }
}

/// Decomposes an input signal into 161 trits.
#[derive(Debug, Clone)]
pub struct FullTritDecomposerChip {
    decomposer: TritDecomposerChip<161>,
    comparator: TritComparatorChip<161>,
}

impl Default for FullTritDecomposerChip {
    fn default() -> Self {
        Self {
            decomposer: TritDecomposerChip::default(),
            comparator: TritComparatorChip::new(Scalar::MODULUS.parse().unwrap()),
        }
    }
}

impl PlonkChip<1, 161> for FullTritDecomposerChip {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; 161]> {
        let trits = self.decomposer.build(builder, inputs)?;
        let cmp = self.comparator.build(builder, trits)?[0].unwrap();
        let c = builder.add_const_gate(-Scalar::from_const(1));
        builder.connect(cmp, c);
        Ok(trits)
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; 161]> {
        let trits = self.decomposer.witness(witness, inputs)?;
        self.comparator.witness(witness, trits)?;
        witness.assert_constant(-Scalar::from_const(1));
        Ok(trits)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::NUM_BLINDING_ROWS;
    use std::cmp::Ordering;

    const BLOWUP_LOG2: usize = 1;

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
        let mut builder = CircuitBuilder::default();
        let input = builder.add_const_gate(value.into());
        let chip = BitDecomposerChip::<N>::default();
        assert!(chip.build(&mut builder, [Some(input)]).is_ok());
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        let input = witness.assert_constant(value.into());
        let bits = chip
            .witness(&mut witness, [input.into()])
            .unwrap()
            .map(|bit| match bit {
                WireOrUnconstrained::Wire(wire) => witness.get(wire),
                _ => panic!("the output bits must be constrained"),
            });
        assert_eq!(bits, decompose_bits::<N>(value.into())[0..N]);
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, BLOWUP_LOG2)
            .unwrap();
        assert!(
            circuit
                .to_compressed::<Sha2Hash<Scalar>>(BLOWUP_LOG2)
                .verify(&proof)
                .is_ok()
        );
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

    fn test_bit_comparator_chip<const N: usize>(lhs: u64, rhs: u64) {
        let mut builder = CircuitBuilder::default();
        let input = builder.add_const_gate(lhs.into());
        let decomposer_chip = BitDecomposerChip::<N>::default();
        let bits = decomposer_chip.build(&mut builder, [input.into()]).unwrap();
        let comparator_chip = BitComparatorChip::<N>::new(rhs.into());
        let cmp = comparator_chip.build(&mut builder, bits).unwrap()[0].unwrap();
        builder.declare_public_gates([input.gate(), cmp.gate()]);
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        let input = witness.assert_constant(lhs.into());
        let bits = decomposer_chip
            .witness(&mut witness, [input.into()])
            .unwrap();
        assert!(comparator_chip.witness(&mut witness, bits).is_ok());
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, BLOWUP_LOG2)
            .unwrap();
        let openings = circuit
            .to_compressed::<Sha2Hash<Scalar>>(BLOWUP_LOG2)
            .verify(&proof)
            .unwrap();
        assert_eq!(openings[&input], lhs.into());
        assert_eq!(
            openings[&cmp],
            match lhs.cmp(&rhs) {
                Ordering::Less => -Scalar::from_const(1),
                Ordering::Equal => Scalar::from_const(0),
                Ordering::Greater => Scalar::from_const(1),
            }
        );
    }

    #[test]
    fn test_bit_comparator_chip_1() {
        test_bit_comparator_chip::<1>(0, 0);
        test_bit_comparator_chip::<1>(1, 0);
        test_bit_comparator_chip::<1>(0, 1);
        test_bit_comparator_chip::<1>(1, 1);
    }

    #[test]
    fn test_bit_comparator_chip_2() {
        test_bit_comparator_chip::<2>(0, 0);
        test_bit_comparator_chip::<2>(1, 0);
        test_bit_comparator_chip::<2>(2, 0);
        test_bit_comparator_chip::<2>(3, 0);
        test_bit_comparator_chip::<2>(0, 1);
        test_bit_comparator_chip::<2>(1, 1);
        test_bit_comparator_chip::<2>(2, 1);
        test_bit_comparator_chip::<2>(3, 1);
        test_bit_comparator_chip::<2>(0, 2);
        test_bit_comparator_chip::<2>(1, 2);
        test_bit_comparator_chip::<2>(2, 2);
        test_bit_comparator_chip::<2>(3, 2);
        test_bit_comparator_chip::<2>(0, 3);
        test_bit_comparator_chip::<2>(1, 3);
        test_bit_comparator_chip::<2>(2, 3);
        test_bit_comparator_chip::<2>(3, 3);
    }

    fn test_full_bit_decomposer_chip_impl(value: u64) {
        let mut builder = CircuitBuilder::default();
        let input = builder.add_const_gate(value.into());
        let chip = FullBitDecomposerChip::default();
        assert!(chip.build(&mut builder, [Some(input)]).is_ok());
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        let input = witness.assert_constant(value.into());
        let bits = chip
            .witness(&mut witness, [input.into()])
            .unwrap()
            .map(|bit| match bit {
                WireOrUnconstrained::Wire(wire) => witness.get(wire),
                _ => panic!("the output bits must be constrained"),
            });
        assert_eq!(bits, decompose_bits::<256>(value.into()));
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, BLOWUP_LOG2)
            .unwrap();
        assert!(
            circuit
                .to_compressed::<Sha2Hash<Scalar>>(BLOWUP_LOG2)
                .verify(&proof)
                .is_ok()
        );
    }

    #[test]
    fn test_full_bit_decomposer_chip() {
        test_full_bit_decomposer_chip_impl(0);
        test_full_bit_decomposer_chip_impl(1);
        test_full_bit_decomposer_chip_impl(2);
        test_full_bit_decomposer_chip_impl(3);
        test_full_bit_decomposer_chip_impl(4);
        test_full_bit_decomposer_chip_impl(5);
        test_full_bit_decomposer_chip_impl(6);
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
        let mut builder = CircuitBuilder::default();
        let input = builder.add_const_gate(value.into());
        let chip = TritDecomposerChip::<N>::default();
        assert!(chip.build(&mut builder, [Some(input)]).is_ok());
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        let input = witness.assert_constant(value.into());
        let trits = chip
            .witness(&mut witness, [input.into()])
            .unwrap()
            .map(|trit| match trit {
                WireOrUnconstrained::Wire(wire) => witness.get(wire),
                _ => panic!("the output trits must be constrained"),
            });
        assert_eq!(trits, decompose_trits::<N>(value.into())[0..N]);
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, BLOWUP_LOG2)
            .unwrap();
        assert!(
            circuit
                .to_compressed::<Sha2Hash<Scalar>>(BLOWUP_LOG2)
                .verify(&proof)
                .is_ok()
        );
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

    fn test_trit_comparator_chip<const N: usize>(lhs: u64, rhs: u64) {
        let mut builder = CircuitBuilder::default();
        let input = builder.add_const_gate(lhs.into());
        let decomposer_chip = TritDecomposerChip::<N>::default();
        let trits = decomposer_chip.build(&mut builder, [input.into()]).unwrap();
        let comparator_chip = TritComparatorChip::<N>::new(rhs.into());
        let cmp = comparator_chip.build(&mut builder, trits).unwrap()[0].unwrap();
        builder.declare_public_gates([input.gate(), cmp.gate()]);
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        let input = witness.assert_constant(lhs.into());
        let trits = decomposer_chip
            .witness(&mut witness, [input.into()])
            .unwrap();
        assert!(comparator_chip.witness(&mut witness, trits).is_ok());
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, BLOWUP_LOG2)
            .unwrap();
        let openings = circuit
            .to_compressed::<Sha2Hash<Scalar>>(BLOWUP_LOG2)
            .verify(&proof)
            .unwrap();
        assert_eq!(openings[&input], lhs.into());
        assert_eq!(
            openings[&cmp],
            match lhs.cmp(&rhs) {
                Ordering::Less => -Scalar::from_const(1),
                Ordering::Equal => Scalar::from_const(0),
                Ordering::Greater => Scalar::from_const(1),
            }
        );
    }

    #[test]
    fn test_trit_comparator_chip_1() {
        test_trit_comparator_chip::<1>(0, 0);
        test_trit_comparator_chip::<1>(0, 1);
        test_trit_comparator_chip::<1>(0, 2);
        test_trit_comparator_chip::<1>(1, 0);
        test_trit_comparator_chip::<1>(1, 1);
        test_trit_comparator_chip::<1>(1, 2);
        test_trit_comparator_chip::<1>(2, 0);
        test_trit_comparator_chip::<1>(2, 1);
        test_trit_comparator_chip::<1>(2, 2);
    }

    #[test]
    fn test_trit_comparator_chip_2() {
        for lhs in 0..9 {
            for rhs in 0..9 {
                test_trit_comparator_chip::<2>(lhs, rhs);
            }
        }
    }

    fn test_full_trit_decomposer_chip_impl(value: u64) {
        let mut builder = CircuitBuilder::default();
        let input = builder.add_const_gate(value.into());
        let chip = FullTritDecomposerChip::default();
        assert!(chip.build(&mut builder, [Some(input)]).is_ok());
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        let input = witness.assert_constant(value.into());
        let trits = chip
            .witness(&mut witness, [input.into()])
            .unwrap()
            .map(|trit| match trit {
                WireOrUnconstrained::Wire(wire) => witness.get(wire),
                _ => panic!("the output trits must be constrained"),
            });
        assert_eq!(trits, decompose_trits::<161>(value.into()));
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, BLOWUP_LOG2)
            .unwrap();
        assert!(
            circuit
                .to_compressed::<Sha2Hash<Scalar>>(BLOWUP_LOG2)
                .verify(&proof)
                .is_ok()
        );
    }

    #[test]
    fn test_full_trit_decomposer_chip_1() {
        test_full_trit_decomposer_chip_impl(0);
        test_full_trit_decomposer_chip_impl(1);
        test_full_trit_decomposer_chip_impl(2);
        test_full_trit_decomposer_chip_impl(3);
        test_full_trit_decomposer_chip_impl(4);
        test_full_trit_decomposer_chip_impl(5);
        test_full_trit_decomposer_chip_impl(6);
        test_full_trit_decomposer_chip_impl(7);
        test_full_trit_decomposer_chip_impl(8);
    }

    #[test]
    fn test_full_trit_decomposer_chip_2() {
        test_full_trit_decomposer_chip_impl(9);
        test_full_trit_decomposer_chip_impl(10);
        test_full_trit_decomposer_chip_impl(11);
        test_full_trit_decomposer_chip_impl(12);
        test_full_trit_decomposer_chip_impl(13);
        test_full_trit_decomposer_chip_impl(14);
        test_full_trit_decomposer_chip_impl(15);
        test_full_trit_decomposer_chip_impl(16);
        test_full_trit_decomposer_chip_impl(17);
    }

    #[test]
    fn test_full_trit_decomposer_chip_3() {
        test_full_trit_decomposer_chip_impl(18);
        test_full_trit_decomposer_chip_impl(19);
        test_full_trit_decomposer_chip_impl(20);
        test_full_trit_decomposer_chip_impl(21);
        test_full_trit_decomposer_chip_impl(22);
        test_full_trit_decomposer_chip_impl(23);
        test_full_trit_decomposer_chip_impl(24);
        test_full_trit_decomposer_chip_impl(25);
        test_full_trit_decomposer_chip_impl(26);
    }
}
