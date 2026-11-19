use anyhow::Result;
use starkom_ff::{Field256, PrimeField};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, make_const,
    var,
};
use starkom_poseidon2 as poseidon2;
use std::fmt::{Debug, Formatter};
use std::marker::PhantomData;

/// Generic Poseidon2 permutation chip.
///
/// This chip achieves a full permutation in a single row using
/// degree-[alpha](`poseidon2::Config::alpha`) constraints. For example, this results in degree-5
/// constraints on BlueSky and degree-3 constraints on Schraderbrau.
pub struct PermutationChip<F: PrimeField, C: poseidon2::Config<F, T>, const T: usize> {
    _data: PhantomData<(F, C)>,
}

impl<F: PrimeField, C: poseidon2::Config<F, T>, const T: usize> Debug for PermutationChip<F, C, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PermutationChip").finish()
    }
}

impl<F: PrimeField, C: poseidon2::Config<F, T>, const T: usize> Default
    for PermutationChip<F, C, T>
{
    fn default() -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon2::Config<F, T>, const T: usize> Copy for PermutationChip<F, C, T> {}

impl<F: PrimeField, C: poseidon2::Config<F, T>, const T: usize> Clone for PermutationChip<F, C, T> {
    fn clone(&self) -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon2::Config<F, T>, const T: usize> PlonkChip<F, T, T>
    for PermutationChip<F, C, T>
{
    fn width(&self) -> usize {
        T * (3 + C::num_full_rounds_per_side() * 2) + C::num_partial_rounds() - 1
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
        let external_matrix = C::get_external_matrix();
        let internal_matrix = C::get_internal_matrix();
        let internal_matrix_00_inverse = internal_matrix[0].invert_unwrap();

        for i in 0..T {
            view.add_gate(
                0,
                (0..T)
                    .map(|j| var(j) * make_const(external_matrix[i * T + j]))
                    .sum::<Constraint<F>>()
                    - var(T + i),
            );
        }
        let mut state: [Constraint<F>; T] = std::array::from_fn(|i| var(T + i));
        let mut offset = T * 2;

        for r in 0..num_full_rounds_per_side {
            let sbox_outputs: [Constraint<F>; T] =
                std::array::from_fn(|j| (make_const(c[r * T + j]) + state[j].clone()) ^ a);
            for i in 0..T {
                view.add_gate(
                    0,
                    (0..T)
                        .map(|j| sbox_outputs[j].clone() * make_const(external_matrix[i * T + j]))
                        .sum::<Constraint<F>>()
                        - var(offset + i),
                );
            }
            state = std::array::from_fn(|i| var(offset + i));
            offset += T;
        }

        for r in num_full_rounds_per_side..(num_full_rounds_per_side + num_partial_rounds) {
            let linear: [Constraint<F>; T] = std::array::from_fn(|i| {
                (1..T)
                    .map(|j| state[j].clone() * make_const(internal_matrix[i * T + j]))
                    .sum()
            });
            view.add_gate(
                0,
                ((make_const(c[r * T]) + state[0].clone()) ^ a) * make_const(internal_matrix[0])
                    + linear[0].clone()
                    - var(offset),
            );
            for i in 1..T {
                state[i] = (var(offset) - linear[0].clone())
                    * make_const(internal_matrix_00_inverse * internal_matrix[i * T])
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
            let sbox_outputs: [Constraint<F>; T] =
                std::array::from_fn(|j| (make_const(c[r * T + j]) + state[j].clone()) ^ a);
            for i in 0..T {
                view.add_gate(
                    0,
                    (0..T)
                        .map(|j| sbox_outputs[j].clone() * make_const(external_matrix[i * T + j]))
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
        let inputs = inputs.map(|input| view.get(input));

        let num_full_rounds_per_side = C::num_full_rounds_per_side();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        let c = C::get_round_constants();
        let a = C::alpha();
        let external_matrix = C::get_external_matrix();
        let internal_matrix = C::get_internal_matrix();

        let mut state: [F; T] =
            std::array::from_fn(|i| (0..T).map(|j| inputs[j] * external_matrix[i * T + j]).sum());
        for i in 0..T {
            view.set(view.cell(0, T + i), state[i]);
        }
        let mut offset = T * 2;

        for r in 0..num_full_rounds_per_side {
            let sbox_outputs: [F; T] =
                std::array::from_fn(|j| (c[r * T + j] + state[j]).pow_small(a));
            state = std::array::from_fn(|i| {
                (0..T)
                    .map(|j| sbox_outputs[j] * external_matrix[i * T + j])
                    .sum()
            });
            for i in 0..T {
                view.set(view.cell(0, offset + i), state[i]);
            }
            offset += T;
        }

        for r in num_full_rounds_per_side..(num_full_rounds_per_side + num_partial_rounds) {
            let sbox_output = (c[r * T] + state[0]).pow_small(a);
            state = std::array::from_fn(|i| {
                std::iter::once(sbox_output * internal_matrix[i * T])
                    .chain((1..T).map(|j| state[j] * internal_matrix[i * T + j]))
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
            let sbox_outputs: [F; T] =
                std::array::from_fn(|j| (c[r * T + j] + state[j]).pow_small(a));
            state = std::array::from_fn(|i| {
                (0..T)
                    .map(|j| sbox_outputs[j] * external_matrix[i * T + j])
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
    use starkom_ff::{PrimeField64, PrimeField256};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CircuitOptions};
    use starkom_poseidon2::Config;
    use std::fmt::Debug;
    use std::str::FromStr;

    #[cfg(feature = "bluesky")]
    use starkom_bluesky::Scalar as BS;

    #[cfg(feature = "goldilocks")]
    use starkom_goldilocks::{GL, GL4, from_const as gl};

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
            T * (3 + C::num_full_rounds_per_side() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, F>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(CircuitOptions {
            canonicalize_constraints: false,
            security_bits: SECURITY_BITS,
            blowup_log2,
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
        let proof = circuit.prove::<Sha2Hash<F>>(witness, &circuit.default_proof_options())?;
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
            parse("0x6f30582cde48a25b26015b7f718ba2fb359e93029caf04d8d0b3e66b1d46b941"),
            parse("0x5de8159372063ce76403529bb1a9725461b96467035d906400ff48d0937f9db6"),
            parse("0x3c88b37dc6d14d08960b6fe58344e09194d11a930ce9f60cc90294683fac4b9f"),
        ];
        assert!(
            test_permutation256_impl::<BS, poseidon2::BlueSkyConfig3, 3>(
                inputs,
                outputs,
                1,
                parse("0x1ce9cd7e81d1ebe39d506dbb3ce74f93c00ef134ed4ea8b64bee95801168ba81")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon2::BlueSkyConfig3, 3>(
                inputs,
                outputs,
                2,
                parse("0xe8906483f0a4b55088884fb7d219ca01a781a62c2653a133bbfc141bcf59772b")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon2::BlueSkyConfig3, 3>(
                inputs,
                outputs,
                3,
                parse("0x3c5bcdc92f2f42612b7e069e11bf33957629feafaef46e3b5d4fae5f237fe099")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_permutation_t4_bs() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x775049834d9decb40ec5a109116a27527fa9105a3521cee8a42777788fda1501"),
            parse("0x630ded08b39ceac4859c9ab6d14b548f48d01164ce1efada3a7a868f7d9248cb"),
            parse("0x14b47f414dececb9936dcbb89e2fdd8511c44acb30439d1d23e48119b1c03b4f"),
            parse("0x72de70292ce1ac7f30b859d04bbb6de5377288c1192a08863c34e11bc9269c4c"),
        ];
        assert!(
            test_permutation256_impl::<BS, poseidon2::BlueSkyConfig4, 4>(
                inputs,
                outputs,
                1,
                parse("0x75a7539c42fda544ea2302677627fc6b2f5b51feb778035409b248da57dd4640")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon2::BlueSkyConfig4, 4>(
                inputs,
                outputs,
                2,
                parse("0x9d8bee4dab499b1e3abffa2ee765eef3662295c32e0c40a6ace0ab63bc6e5a64")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon2::BlueSkyConfig4, 4>(
                inputs,
                outputs,
                3,
                parse("0xb8335edba106318875c1ceda0f6653f822c81ee3811ae75268caacdc9c0a1221")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "schraderbrau")]
    #[test]
    fn test_permutation_t3_sb() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into()];
        let outputs = [
            parse("0x4944db4e43f190b4efd2383ca68fc624f8c8f1f3cd6396aaccdc9f964011998b"),
            parse("0x48f659c353026326543ca20d6162afc5e2dc6871444e3271c500d9fd9a203ba3"),
            parse("0x264f179e410ea3505c03b2f726e8229a16bc796f0213ca5e325324f4407e1d43"),
        ];
        assert!(
            test_permutation256_impl::<SB, poseidon2::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                1,
                parse("0xf0092b7e94468f8c4a550a2f832519c344ca246b29b7c28d8bbc56a5d757a1dc")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon2::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                2,
                parse("0x93d8acdc272d92cccd79e9da307a15c271cb952686948ac3d318d38afdafb0f6")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon2::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                3,
                parse("0xa5bfe82585cd679270098804dd6760b4745a1fc4609181a766f6bb97d42bf728")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "schraderbrau")]
    #[test]
    fn test_permutation_t4_sb() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x0ed3f04e32c5a56b534177f0d5080e7808eb14b0e5a00087a8d9c1b31e15caa3"),
            parse("0x41cfd908539849048aa49853ffc331f6610ac53030220809863e41d3ccbd2db0"),
            parse("0x75b114742d553230c8d4aa77788c3c3eba96af41a79c6f43f1dbd500110facac"),
            parse("0x240692e1b40ee3a80e4ec79414f815d291866e817e59abee4d23d99cd547b8c7"),
        ];
        assert!(
            test_permutation256_impl::<SB, poseidon2::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                1,
                parse("0xd994e4f0e0dd0a0017760fa160dd2d22beded9727443810db12cb3721cb5c620")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon2::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                2,
                parse("0xef8608aa72dff40745fa5a5f38ba7b4ac7632a7dc0b1542ea089b5cc5ae73ed3")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon2::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                3,
                parse("0x58840a1be23fcfbcbdd88fd9ace76c041d58aabf73f64bc99f58f87df0a814df")
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
            T * (3 + C::num_full_rounds_per_side() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, G>::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(CircuitOptions {
            canonicalize_constraints: false,
            security_bits: SECURITY_BITS,
            blowup_log2,
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
        let proof = circuit.prove::<Sha2Hash<G>>(witness, &circuit.default_proof_options())?;
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
            gl(0x01eaef96bdf1c0c1),
            gl(0x1f0d2cc525b2540c),
            gl(0x6282c1dfe1e0358d),
            gl(0xe780d721f698e1e6),
            gl(0x280c0b6f753d833b),
            gl(0x1b942dd5023156ab),
            gl(0x43f0df3fcccb8398),
            gl(0xe8e8190585489025),
            gl(0x56bdbf72f77ada22),
            gl(0x7911c32bf9dcd705),
            gl(0xec467926508fbe67),
            gl(0x6a50450ddf85a6ed),
        ];
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon2::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                1,
                parse("0x5ee8dab7b3ea5466a20dc2655de4e732775564956dca02b73fcd0d1000ad8d29")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon2::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0x3a7423f485cc74ec20606e59ed49f4c84413fe74f64c95a358d13052e60933e9")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon2::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0xd68330f9ae77d844bec97aa3ef099a88e8b424cb015ce1d2fdd7f16c24954071")
            )
            .is_ok()
        );
    }

    #[cfg(feature = "goldilocks")]
    #[test]
    fn test_permutation_t16_gl() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            gl(0x85c54702470d9756),
            gl(0xaa53c7a7d52d9898),
            gl(0x285128096efb0dd7),
            gl(0xf3fde5edd3050ac8),
            gl(0xc7b65efd040df908),
            gl(0x4be3f6c467f57ae9),
            gl(0x274e9a67b41754fb),
            gl(0x0f7d39cd5de94dac),
            gl(0xd0224b9794d0b78c),
            gl(0x372f6139570042e1),
            gl(0xce6e8a93dc4ec26c),
            gl(0xace65e30a4daf7af),
            gl(0x016f2824cc1ba3db),
            gl(0x2e8f3af37c434dec),
            gl(0xc80831bb6e09da01),
            gl(0x3a7d670bf1a86ee8),
        ];
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon2::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                1,
                parse("0x6e301c0fe83b16b9235d34542060c9f536816cea1ff578061a0cdcea39a2912c")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon2::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0xa0a749eef26755980bd249d46d8c79e23600a54bff93f4d6362ddea1cee4711b")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon2::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0xe087520a3f11500eee7320e395063f3a7ed48ae84f72edae096c9420d16adb4d")
            )
            .is_ok()
        );
    }
}
