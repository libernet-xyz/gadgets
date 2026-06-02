use crate::utils;
use anyhow::Result;
use ff::{Field, PrimeField};
use primitive_types::U256;
use starkom_bluesky::Scalar;
use starkom_plonk::{Chip as PlonkChip, CircuitBuilder, Wire, WireOrUnconstrained, Witness};

/// Calculates (1 - X^2).
///
/// Not an actual "logical NOT", it can only NOT the (-1,0,1) signals we use in the comparator
/// chips. `LogicalNotChip` is for internal use by those chips.
#[derive(Debug, Default)]
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
    let lsb = value.to_little_endian()[0];
    Scalar::from((lsb & 1) as u64)
}

pub fn shr(value: Scalar, count: usize) -> Scalar {
    utils::u256_to_scalar(utils::scalar_to_u256(value) >> U256::from(count)).unwrap()
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
    decompose_bits::<N>(utils::scalar_to_u256(value))
}

/// Decomposes the input signal into N bits.
///
/// WARNING: this chip is unsafe to use with 255 or 256 bits because it doesn't guard against
/// aliasing. Use the `FullBitDecomposerChip` for a full decomposition into 256 bits (note that the
/// MSB will always be zero because BLS12-381 scalars don't cover the upper half of the 256 bit
/// range).
#[derive(Debug, Default)]
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
            sum = builder.add_linear_combination_gate(1.into(), sum.into(), power, None);
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
            sum = witness.combine(1.into(), sum.into(), power, bit.into());
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
#[derive(Debug)]
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
        utils::u256_to_scalar((self.rhs >> i) & 1.into()).unwrap()
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
#[derive(Debug)]
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

// TODO

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use super::*;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::NUM_BLINDING_ROWS;

    const BLOWUP_LOG2: usize = 1;

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

    // TODO
}
