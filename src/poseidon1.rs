use anyhow::Result;
use starkom_ff::{Field256, PrimeField};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, make_const,
    var,
};
use starkom_poseidon as poseidon1;
use std::fmt::{Debug, Formatter};
use std::marker::PhantomData;

/// Generic Poseidon1 permutation chip.
///
/// This chip achieves a full permutation in a single row using
/// degree-[alpha](`poseidon1::Config::alpha`) constraints. For example, this results in degree-5
/// constraints on BlueSky and degree-3 constraints on Schraderbrau.
pub struct PermutationChip<F: PrimeField, C: poseidon1::Config<F, T>, const T: usize> {
    _data: PhantomData<(F, C)>,
}

impl<F: PrimeField, C: poseidon1::Config<F, T>, const T: usize> Debug for PermutationChip<F, C, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PermutationChip").finish()
    }
}

impl<F: PrimeField, C: poseidon1::Config<F, T>, const T: usize> Default
    for PermutationChip<F, C, T>
{
    fn default() -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon1::Config<F, T>, const T: usize> Copy for PermutationChip<F, C, T> {}

impl<F: PrimeField, C: poseidon1::Config<F, T>, const T: usize> Clone for PermutationChip<F, C, T> {
    fn clone(&self) -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon1::Config<F, T>, const T: usize> PlonkChip<F, T, T>
    for PermutationChip<F, C, T>
{
    fn width(&self) -> usize {
        T * (2 + C::num_full_rounds_per_side() * 2) + C::num_partial_rounds() - 1
    }

    fn height(&self) -> usize {
        1
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; T],
    ) -> Result<[Option<Cell>; T]> {
        for i in 0..T {
            view.connect(inputs[i], view.cell(0, i).into());
        }

        let num_full_rounds_per_side = C::num_full_rounds_per_side();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        let c = C::get_round_constants();
        let a = C::alpha() as isize;
        let m = C::get_mds_matrix();
        let m0_inv = m[0].invert_unwrap();

        let mut offset = 0;
        for r in 0..num_full_rounds_per_side {
            for i in 0..T {
                view.add_gate(
                    0,
                    (0..T)
                        .map(|j| {
                            ((make_const(c[r * T + j]) + var(offset + j)) ^ a)
                                * make_const(m[i * T + j])
                        })
                        .sum::<Constraint<F>>()
                        - var(offset + T + i),
                );
            }
            offset += T;
        }

        let mut state: [Constraint<F>; T] = std::array::from_fn(|i| var(offset + i));
        offset += T;

        for r in num_full_rounds_per_side..(num_full_rounds_per_side + num_partial_rounds) {
            let linear: [Constraint<F>; T] = std::array::from_fn(|i| {
                (1..T)
                    .map(|j| {
                        (make_const(c[r * T + j]) + state[j].clone()) * make_const(m[i * T + j])
                    })
                    .sum()
            });
            view.add_gate(
                0,
                ((make_const(c[r * T]) + state[0].clone()) ^ a) * make_const(m[0])
                    + linear[0].clone()
                    - var(offset),
            );
            for i in 1..T {
                state[i] =
                    (var(offset) - linear[0].clone()) * make_const(m0_inv) * make_const(m[i * T])
                        + linear[i].clone();
            }
            state[0] = var(offset);
            offset += 1;
        }

        for i in 1..T {
            view.add_gate(0, state[i].clone() - var(offset + i - 1));
            state[i] = var(offset + i - 1);
        }
        offset += T - 1;

        for r in (num_full_rounds_per_side + num_partial_rounds)..num_total_rounds {
            for i in 0..T {
                view.add_gate(
                    0,
                    (0..T)
                        .map(|j| {
                            ((make_const(c[r * T + j]) + state[j].clone()) ^ a)
                                * make_const(m[i * T + j])
                        })
                        .sum::<Constraint<F>>()
                        - var(offset + i),
                );
            }
            state = std::array::from_fn(|i| var(offset + i));
            offset += T;
        }

        Ok(std::array::from_fn(|i| view.cell(0, offset - T + i).into()))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; T],
    ) -> Result<[CellOrUnconstrained<F>; T]> {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let mut state = inputs.map(|input| view.get(input));

        let num_full_rounds_per_side = C::num_full_rounds_per_side();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        let c = C::get_round_constants();
        let a = C::alpha();
        let m = C::get_mds_matrix();

        let mut offset = T;
        for r in 0..num_full_rounds_per_side {
            state = std::array::from_fn(|i| {
                (0..T)
                    .map(|j| (c[r * T + j] + state[j]).pow_small(a) * m[i * T + j])
                    .sum()
            });
            for i in 0..T {
                view.set(view.cell(0, offset + i), state[i]);
            }
            offset += T;
        }

        for r in num_full_rounds_per_side..(num_full_rounds_per_side + num_partial_rounds) {
            state = std::array::from_fn(|i| {
                std::iter::once((c[r * T] + state[0]).pow_small(a) * m[i * T])
                    .chain((1..T).map(|j| (c[r * T + j] + state[j]) * m[i * T + j]))
                    .sum()
            });
            view.set(view.cell(0, offset), state[0]);
            offset += 1;
        }

        for i in 1..T {
            view.set(view.cell(0, offset + i - 1), state[i]);
        }
        offset += T - 1;

        for r in (num_full_rounds_per_side + num_partial_rounds)..num_total_rounds {
            state = std::array::from_fn(|i| {
                (0..T)
                    .map(|j| (c[r * T + j] + state[j]).pow_small(a) * m[i * T + j])
                    .sum()
            });
            for i in 0..T {
                view.set(view.cell(0, offset + i), state[i]);
            }
            offset += T;
        }

        Ok(std::array::from_fn(|i| view.cell(0, offset - T + i).into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use primitive_types::H256;
    use starkom_ff::{PrimeField32, PrimeField64, PrimeField256};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, Options};
    use starkom_poseidon::Config;
    use std::fmt::Debug;
    use std::str::FromStr;

    #[cfg(feature = "bluesky")]
    use starkom_bluesky::Scalar as BS;

    #[cfg(feature = "goldilocks")]
    use starkom_goldilocks::{GL, GL4, from_const as gl};

    #[cfg(feature = "koalabear")]
    use starkom_koalabear::{KB, KB8, from_const as kb};

    #[cfg(feature = "schraderbrau")]
    use starkom_schraderbrau::Scalar as SB;

    const SECURITY_BITS: usize = 128;

    fn parse<T: FromStr<Err: Debug>>(s: &'static str) -> T {
        s.parse().unwrap()
    }

    fn test_permutation256_impl<F: PrimeField256, C: Config<F, T>, const T: usize>(
        inputs: [F; T],
        expected_output: [F; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChip::<F, C, T>::default();
        assert_eq!(
            chip.width(),
            T * (2 + C::num_full_rounds_per_side() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, F>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(Options {
            canonicalize_constraints: false,
            security_bits: SECURITY_BITS,
            blowup_log2,
            fri_fold_rate: 1,
            max_fri_rounds: usize::MAX,
            grinding_bits: 0,
            blind: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.num_columns(), chip.width());
        assert_eq!(circuit.get_max_gate_degree(), C::alpha());
        let degree_bound = circuit.pcs_config().degree_bound;
        assert!(degree_bound <= 256);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), degree_bound);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<F>>(witness)?;
        assert!(proof.pcs_config().satisfies(circuit.pcs_config()));
        let circuit = circuit.to_compressed::<Sha2Hash<F>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let public_inputs = circuit.verify(&proof)?;
        assert!(
            output
                .into_iter()
                .enumerate()
                .all(|(i, output)| match output {
                    CellOrUnconstrained::Cell(cell) => public_inputs[&cell],
                    CellOrUnconstrained::Unconstrained(value) => value,
                } == expected_output[i])
        );
        Ok(())
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_permutation_t3_bs() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into()];
        let outputs = [
            parse("0x7b68dcd80fa751ee8f2d76043bfd92c685601c79189393fc76e03c5214eed32b"),
            parse("0x0fbcb5720b463bf7e2ccabf373e77d2c10d27e6549f34cfa33eb2d06ea8b900a"),
            parse("0x26e03abfcc62da0101516b07aede8bc676a10c47299a57bedc6d9fe80484f3da"),
        ];
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT3X5, 3>(
                inputs,
                outputs,
                1,
                parse("0x9a85bd8655a910c9912d0b9faab4bd7a34f5d08c9df81f570c56a39fb89524b0")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT3X5, 3>(
                inputs,
                outputs,
                2,
                parse("0x49ab43bd49caca8544666800d372bcf77eab4d688565fbc6386cdd3fe7d0ffcd")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT3X5, 3>(
                inputs,
                outputs,
                3,
                parse("0x1d7224ff32d2bff94bb50b62febb55f48b9bd646010701dca887b61842b6c99f")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_permutation_t4_bs() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x12dde8a4c46760e349670d241e36ca7abacc991233039f8deaf6c58ce2230ef6"),
            parse("0x61e95d9456e9223b4d7926dabae10009da2b6fb9134ade8405f6ef1424e66aa1"),
            parse("0x2fcce25ab9efb3e26276f3b3aff1e02cdf82df48ce8d3eadbff900cfe015775b"),
            parse("0x2580707d57a8c1c0cad368e8d5705ffd96f269d66e1cd6f1433f93a3c66d9bf8"),
        ];
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT4X5, 4>(
                inputs,
                outputs,
                1,
                parse("0xd64c9568ca91f8745b55011856ee201ec2f58013525eb5b99453eb83e8018892")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT4X5, 4>(
                inputs,
                outputs,
                2,
                parse("0xebec1edb924ac7b6501aab1a4b38111e2b1bad6cfbb99cad20e4f924f325da21")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT4X5, 4>(
                inputs,
                outputs,
                3,
                parse("0xfc0b4307a220a99cb46dbe50aaa1d8f1db6fa0faab15ca3674055d5d7d532de7")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "schraderbrau")]
    #[test]
    fn test_permutation_t3_sb() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into()];
        let outputs = [
            parse("0x1531e124d8f663b8d26c56b9bf0b1b09ba15d4e20ea10a3d898cdcf1d2c41dba"),
            parse("0x1f5714cc13f8a33f4b32a07f1a85409de2d1b8d353aa269ca8f2bfa91071bd62"),
            parse("0x588c20c682f8b66c52049d910a4f0e7195dfd300a0d0cf3198b1383be1a4134a"),
        ];
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                1,
                parse("0x8a1d4149b1cc2149e06571e39bba4cbbf4e1049837afbacbf5cb7ddb61a8c47e")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                2,
                parse("0x096e2dc0e7dcc8604ac5b4ad4da6f34ea2c4fcd0ddd27fbc1325f2e6cedabfbc")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                3,
                parse("0xc0a5ad2226c570b1c95d3b4fa2f775754a89443514f20c7631e20d12d6b76c72")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "schraderbrau")]
    #[test]
    fn test_permutation_t4_sb() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x11de0b3702563747b0729abd2b93e720ec947ce067ca3f8c088b9829df1169d7"),
            parse("0x487ca87fa054ad960f6f571cf7aa6f5e075fd5ac3386e2e9fb22aeeb036fb1e1"),
            parse("0x28c775ddba9039531ee5deb0d4e2e6b1c9bb7d8be8da260a20d331bb041d0d38"),
            parse("0x10c93f76dca6f63507bc8ed1d2ea40647bb480ad43ac9922a3f10597b802f949"),
        ];
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                1,
                parse("0xf2298e36088f64b0161b3191a85282d514cab24c84301d87b990a9abe6b85cf8")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                2,
                parse("0xa1c2c4d6c53d5d83877973e64b357abd6ead82968190cece5a762b2731af2c31")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                3,
                parse("0x3e95734e48f65d4cccdd4a5a99f8de55b8ea93645dd70abcd7b94de3b6836f53")
            )
            .is_ok()
        );
    }

    fn test_permutation64_impl<
        F: PrimeField64,
        G: Field256<BaseField = F>,
        C: Config<F, T>,
        const T: usize,
    >(
        inputs: [F; T],
        expected_output: [F; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChip::<F, C, T>::default();
        assert_eq!(
            chip.width(),
            T * (2 + C::num_full_rounds_per_side() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, G>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(Options {
            canonicalize_constraints: false,
            security_bits: SECURITY_BITS,
            blowup_log2,
            fri_fold_rate: 1,
            max_fri_rounds: usize::MAX,
            grinding_bits: 0,
            blind: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.num_columns(), chip.width());
        assert_eq!(circuit.get_max_gate_degree(), C::alpha());
        let degree_bound = circuit.pcs_config().degree_bound;
        assert!(degree_bound <= 1024);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), degree_bound);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness)?;
        assert!(proof.pcs_config().satisfies(circuit.pcs_config()));
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let public_inputs = circuit.verify(&proof)?;
        assert!(
            output
                .into_iter()
                .enumerate()
                .all(|(i, output)| match output {
                    CellOrUnconstrained::Cell(cell) => public_inputs[&cell],
                    CellOrUnconstrained::Unconstrained(value) => value,
                } == expected_output[i])
        );
        Ok(())
    }

    #[cfg(feature = "goldilocks")]
    #[test]
    fn test_permutation_t12_gl() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            gl(0x056bda38ad308e78),
            gl(0x1f38944238b8ccd0),
            gl(0x80bef63a171f3156),
            gl(0x27bbc645b2a3198c),
            gl(0x9befae3f221509b3),
            gl(0xa1cfa54ae2c44c9e),
            gl(0xa1c876869f1c52f8),
            gl(0x7ffa21471eff65af),
            gl(0xdc565450ad52b99e),
            gl(0x4b8b1daf8e8ea3c6),
            gl(0xf866b42495e61984),
            gl(0x7af57b5f91f196fe),
        ];
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                1,
                parse("0x868fb44b8ae9ce0d6f436e53a318eb92742ead07fb10f863170cb405588b70b0")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0xc6d74c485785c83ab0a3f4e0ce9b9317a9a0926f5fa689ab56de4142e72fdf28")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0x6ddf85875d3315c9390c7a4c9fc17934d4fcdcd28983baf80a9f8ffb4a8c389b")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "goldilocks")]
    #[test]
    fn test_permutation_t16_gl() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            gl(0x6a84bf02be1f328d),
            gl(0xec14d274b936a21a),
            gl(0xc0539d7bd4eb66de),
            gl(0xb317ecf41fa8d55b),
            gl(0x80b0d36f66671f8a),
            gl(0x74a1592b9a16e832),
            gl(0x65e53afadfadc8c3),
            gl(0xa0007e5ee96ee4b2),
            gl(0x6dd5661a877003a8),
            gl(0xc36a09c2dc25cd6e),
            gl(0xcbda3d58f7cf85f4),
            gl(0x34cb1d63c35596cf),
            gl(0x4fcd09b24769e281),
            gl(0x6c514f906998c65d),
            gl(0xc447035d8d71952b),
            gl(0x591863454267826f),
        ];
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                1,
                parse("0x474950136c056c9432fec4220fe510a9ddeaa957194cf163695d6216cff304da")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x29581033b8ea2df4f437740366b199ef7e665456326a899e14456acdab080c15")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0x11e4887ff917a2aaefcb5487410bb32b36b7a0a70dda7604f52a62ee2243122f")
            )
            .is_ok()
        );
    }

    fn test_permutation32_impl<
        F: PrimeField32,
        G: Field256<BaseField = F>,
        C: Config<F, T>,
        const T: usize,
    >(
        inputs: [F; T],
        expected_output: [F; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChip::<F, C, T>::default();
        assert_eq!(
            chip.width(),
            T * (2 + C::num_full_rounds_per_side() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, G>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(Options {
            canonicalize_constraints: false,
            security_bits: SECURITY_BITS,
            blowup_log2,
            fri_fold_rate: 1,
            max_fri_rounds: usize::MAX,
            grinding_bits: 0,
            blind: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.num_columns(), chip.width());
        assert_eq!(circuit.get_max_gate_degree(), C::alpha());
        let degree_bound = circuit.pcs_config().degree_bound;
        assert!(degree_bound <= 2048);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), degree_bound);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness)?;
        assert!(proof.pcs_config().satisfies(circuit.pcs_config()));
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let public_inputs = circuit.verify(&proof)?;
        assert!(
            output
                .into_iter()
                .enumerate()
                .all(|(i, output)| match output {
                    CellOrUnconstrained::Cell(cell) => public_inputs[&cell],
                    CellOrUnconstrained::Unconstrained(value) => value,
                } == expected_output[i])
        );
        Ok(())
    }

    #[cfg(feature = "koalabear")]
    #[test]
    fn test_permutation_t24_kb() {
        let inputs = std::array::from_fn(|i| kb(i as u32));
        let outputs = [
            kb(0x20c420b9),
            kb(0x13ee7831),
            kb(0x1664a5fe),
            kb(0x61b8a9da),
            kb(0x62c0a36b),
            kb(0x6686aa16),
            kb(0x086b196f),
            kb(0x2ba9a06e),
            kb(0x48e393fa),
            kb(0x40e4aa30),
            kb(0x5afc6000),
            kb(0x5c9bd903),
            kb(0x3d866242),
            kb(0x44fb43ec),
            kb(0x4fd48650),
            kb(0x45c361f2),
            kb(0x5087cccd),
            kb(0x0c28329f),
            kb(0x58fc4353),
            kb(0x1fa8695a),
            kb(0x4a026a9d),
            kb(0x3e0c59ec),
            kb(0x3be00b5d),
            kb(0x687cd26e),
        ];
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig24, 24>(
                inputs,
                outputs,
                1,
                parse("0x040c69c30775c8ac0d552ae9df76126ba6b06cf2c133337d175238887da2a9d5")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig24, 24>(
                inputs,
                outputs,
                2,
                parse("0xc443e3c6da732ca4795b28c866701ba37261ed091ff680bc22743f59f975cd84")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig24, 24>(
                inputs,
                outputs,
                3,
                parse("0x381a39a0d92b37b572d4bd11dde9505b8657040d9a82dbf35f08c511d6d7586c")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "koalabear")]
    #[test]
    fn test_permutation_t32_kb() {
        let inputs = std::array::from_fn(|i| kb(i as u32));
        let outputs = [
            kb(0x2cd28abd),
            kb(0x4b7e0de6),
            kb(0x1c28010b),
            kb(0x341a4215),
            kb(0x2730f966),
            kb(0x6d5f6492),
            kb(0x2da41e53),
            kb(0x5b9553f8),
            kb(0x40f3bd78),
            kb(0x1847d418),
            kb(0x3b40f194),
            kb(0x61deec2d),
            kb(0x268de447),
            kb(0x4ff70519),
            kb(0x0276d981),
            kb(0x0ecd4498),
            kb(0x5852c99b),
            kb(0x0c0c5656),
            kb(0x1a5d72b5),
            kb(0x1556aa75),
            kb(0x05df8f05),
            kb(0x4e9c5eb7),
            kb(0x669108c9),
            kb(0x051bc5f5),
            kb(0x5d193fdf),
            kb(0x1fc5ddac),
            kb(0x478e26aa),
            kb(0x78fe0af3),
            kb(0x40b6dd16),
            kb(0x6cc4f791),
            kb(0x77e872fe),
            kb(0x6112c017),
        ];
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig32, 32>(
                inputs,
                outputs,
                1,
                parse("0x9eff29013a28e84baf95b1da936ddafae174fd1c5a918ae7e0aa8c78154b44a2")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig32, 32>(
                inputs,
                outputs,
                2,
                parse("0x23678bb6067f7132db3a1582bad083606272fd352d7fb0d5e01a3cac5e4ba46d")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig32, 32>(
                inputs,
                outputs,
                3,
                parse("0x18f81aac25ad1315f265866108ce641d549efaa4a335a2fd3f4c2d1631cdd02a")
            )
            .is_ok()
        );
    }
}
