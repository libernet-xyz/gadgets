use primitive_types::U256;
use starkom_bluesky::{Scalar, from_const};
use starkom_ff::{Field, Field256};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, make_const,
    var,
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

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; 1],
    ) -> anyhow::Result<[Option<Cell>; N]> {
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
    ) -> anyhow::Result<[CellOrUnconstrained; N]> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::parse_scalar;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};

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
        let chip = BitDecomposerChip::<N>::default();
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
        assert!(
            circuit
                .to_compressed::<Sha2Hash<Scalar>>(proving_options)
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

    // TODO
}
