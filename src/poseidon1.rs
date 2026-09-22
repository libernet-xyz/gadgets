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
/// degree-[`ALPHA`](`PrimeField::ALPHA`) constraints. For example, this results in degree-5
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
        T * (2 + C::num_full_rounds() * 2) + C::num_partial_rounds() - 1
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

        let num_full_rounds = C::num_full_rounds();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        let c = C::get_round_constants();
        let a = F::ALPHA as isize;
        let m = C::get_mds_matrix();
        let m0_inv = m[0].invert_unwrap();

        let mut offset = 0;
        for r in 0..num_full_rounds {
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

        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
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

        for r in (num_full_rounds + num_partial_rounds)..num_total_rounds {
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

        let num_full_rounds = C::num_full_rounds();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        let c = C::get_round_constants();
        let a = F::ALPHA;
        let m = C::get_mds_matrix();

        let mut offset = T;
        for r in 0..num_full_rounds {
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

        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
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

        for r in (num_full_rounds + num_partial_rounds)..num_total_rounds {
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
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
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
            T * (2 + C::num_full_rounds() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, F>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.degree_bound(), 4);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), 4);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<F>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 4);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 4 << blowup_log2);
        let circuit = circuit.to_compressed::<Sha2Hash<F>>(options);
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
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfig3, 3>(
                inputs,
                outputs,
                1,
                parse("0xf38f43b1031d10fc36e7b5506f486676298303ec7089a04b4f977fa7aae2faca")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfig3, 3>(
                inputs,
                outputs,
                2,
                parse("0x2c1468b45cd498a6fa8e297b51742d68a541e227f0fcedcc91201c7fd6ef78f5")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfig3, 3>(
                inputs,
                outputs,
                3,
                parse("0xd7cf2642d50b23d2a67a85798949e0d2a8aa7e4b6ad867f75cc75d324c18fa8d")
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
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfig4, 4>(
                inputs,
                outputs,
                1,
                parse("0xb5dafcffa8b4ec594f36a1ba0b101d80a82fbac5a4a96561583bc2bb98c96c91")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfig4, 4>(
                inputs,
                outputs,
                2,
                parse("0x39ccda5c4a841ee8fbf39403b4414d38cf8fa5aadc3aa554caeaf01019473716")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfig4, 4>(
                inputs,
                outputs,
                3,
                parse("0x4d60c26877c88738f171b4899c8c320c118da828eb7a65fead14fef87f6aed35")
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
                parse("0x6cbdfb134e3be4606502e7a77faf60f863e8951a2fc4ed2a4d1c1154e76cb072")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                2,
                parse("0xc822adf88e0162ff0321f5184887d77f5747d157d2ed74530e1cf87af913480d")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                3,
                parse("0x721332b6ca26baf6d8ad8c8ca675e6895e0c066b020c81fe6d7766cf9e587f83")
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
                parse("0xa5eba6d8a98b0c08160cd211c6dd0ca165ce7e69dd5b73ced08178bd4c3e167c")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                2,
                parse("0xa45a0ccc8a371133d2666365f34560bc5b0c6d1a42cfde3b3d403438a94718f2")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                3,
                parse("0xae2fd35f72c75ea4db8bd7a0e19d970a3b90eee6c8e7032d3f0af1a988b0522d")
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
            T * (2 + C::num_full_rounds() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, G>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.degree_bound(), 16);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), 16);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<G>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 16);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 16 << blowup_log2);
        let circuit = circuit.to_compressed::<Sha2Hash<G>>(options);
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
    fn test_permutation_t12_hw() {
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
                parse("0x96c37ae8089e392cf888bef6bce5ec9984429e3f6599f7ba94a35cdb9b274855")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0x821585f96b9af782405b3f1931c34c51eaa1bfed81e3b3ab95bdef9e831aa651")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0x943216ede1aca48240f963cffac6bb951f40154a18ef7cfe8e3af126329ada8b")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "goldilocks")]
    #[test]
    fn test_permutation_t16_hw() {
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
                parse("0xeb296b91f3415bc8eab2e18a3b4084c9f7f0a3d79f34af8528806e3288ce3c0b")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x580bbaf22cc0e8977704d58890e38607ca7cfdaf2842acaee9aa646a49bcce15")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0xafb2727920f5b9877476705b03137eb67c28d87190a7755fdf8d904b90375fd7")
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
            T * (2 + C::num_full_rounds() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, G>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.degree_bound(), 32);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), 32);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<G>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 32);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 32 << blowup_log2);
        let circuit = circuit.to_compressed::<Sha2Hash<G>>(options);
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
                parse("0xaa07e62436079cfefecfc42070bb63f8101726b2b87bf4d1a1ae860be0589d5f")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig24, 24>(
                inputs,
                outputs,
                2,
                parse("0xf9a0f4bc7cf733bf1933d16eab1ab75d957852bb5a9f619d9501faf0a1688ba2")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig24, 24>(
                inputs,
                outputs,
                3,
                parse("0x77d1af53fcb23c6ee510b23e03e8f3cdfbee7ffd8d78ef46137d8ed2ee292304")
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
                parse("0xead7a05ae16a0643a94214297c3565dc78a9513937433b83bc517ae2b3a4f7c4")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig32, 32>(
                inputs,
                outputs,
                2,
                parse("0x1c060a56f5fb0177be4c9ae8565905a981be15aa276a614cab0a5f2c21d6a6fb")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig32, 32>(
                inputs,
                outputs,
                3,
                parse("0x22b05bbdba463cd8f47a6e157fe8af9d7775ecf96771dccd8461515af77bdfac")
            )
            .is_ok()
        );
    }
}
