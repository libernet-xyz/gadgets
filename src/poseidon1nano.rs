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

        let mut offset = 0;
        for r in 0..num_full_rounds {
            let sboxed: [Constraint<F>; T] =
                std::array::from_fn(|j| (make_const(c[r * T + j]) + var(offset + j)) ^ a);
            for i in 0..T {
                view.add_gate(
                    0,
                    (0..T)
                        .map(|j| sboxed[j].clone() * make_const(m[i * T + j]))
                        .sum::<Constraint<F>>()
                        - var(offset + T + i),
                );
            }
            offset += T;
        }

        let mut state: [Constraint<F>; T] = std::array::from_fn(|i| var(offset + i));
        offset += T;
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            let sbox = (make_const(c[r * T]) + state[0].clone()) ^ a;
            let linear: [Constraint<F>; T] = std::array::from_fn(|i| {
                (1..T)
                    .map(|j| {
                        (make_const(c[r * T + j]) + state[j].clone()) * make_const(m[i * T + j])
                    })
                    .sum()
            });
            let next = var(offset);
            view.add_gate(
                0,
                sbox * make_const(m[0]) + linear[0].clone() - next.clone(),
            );
            // Keep the state linear by expressing the S-box output through the gate above.
            let sbox =
                (next.clone() - linear[0].clone()) * make_const(m[0].invert_vartime().unwrap());
            for i in 1..T {
                state[i] = sbox.clone() * make_const(m[i * T]) + linear[i].clone();
            }
            state[0] = next;
            offset += 1;
        }

        for i in 1..T {
            view.add_gate(0, state[i].clone() - var(offset + i - 1));
            state[i] = var(offset + i - 1);
        }
        offset += T - 1;

        for r in (num_full_rounds + num_partial_rounds)..num_total_rounds {
            let sboxed: [Constraint<F>; T] =
                std::array::from_fn(|j| (make_const(c[r * T + j]) + state[j].clone()) ^ a);
            for i in 0..T {
                view.add_gate(
                    0,
                    (0..T)
                        .map(|j| sboxed[j].clone() * make_const(m[i * T + j]))
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
            let sboxed: [F; T] = std::array::from_fn(|j| (c[r * T + j] + state[j]).pow_small(a));
            state = std::array::from_fn(|i| (0..T).map(|j| sboxed[j] * m[i * T + j]).sum());
            for i in 0..T {
                view.set(view.cell(0, offset + i), state[i]);
            }
            offset += T;
        }

        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            let sboxed: [F; T] = std::array::from_fn(|j| {
                let x = c[r * T + j] + state[j];
                if j == 0 { x.pow_small(a) } else { x }
            });
            state = std::array::from_fn(|i| (0..T).map(|j| sboxed[j] * m[i * T + j]).sum());
            view.set(view.cell(0, offset), state[0]);
            offset += 1;
        }

        for i in 1..T {
            view.set(view.cell(0, offset + i - 1), state[i]);
        }
        offset += T - 1;

        for r in (num_full_rounds + num_partial_rounds)..num_total_rounds {
            let sboxed: [F; T] = std::array::from_fn(|j| (c[r * T + j] + state[j]).pow_small(a));
            state = std::array::from_fn(|i| (0..T).map(|j| sboxed[j] * m[i * T + j]).sum());
            for i in 0..T {
                view.set(view.cell(0, offset + i), state[i]);
            }
            offset += T;
        }

        Ok(std::array::from_fn(|i| view.cell(0, offset - T + i).into()))
    }
}

#[cfg(all(test, feature = "schraderbrau"))]
mod tests {
    use super::*;
    use primitive_types::H256;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon::Config;
    use starkom_schraderbrau::Scalar;
    use std::fmt::Debug;
    use std::str::FromStr;

    fn parse<T: FromStr<Err: Debug>>(s: &'static str) -> T {
        s.parse().unwrap()
    }

    fn test_permutation_impl<C: Config<Scalar, T>, const T: usize>(
        inputs: [Scalar; T],
        expected_output: [Scalar; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChip::<Scalar, C, T>::default();
        assert_eq!(
            chip.width(),
            T * (2 + C::num_full_rounds() * 2) + C::num_partial_rounds() - 1
        );
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<Scalar, Scalar>::default();
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
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 4);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 4 << blowup_log2);
        let circuit = circuit.to_compressed::<Sha2Hash<Scalar>>(options);
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

    #[test]
    fn test_permutation_t3() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into()];
        let outputs = [
            parse("0x1531e124d8f663b8d26c56b9bf0b1b09ba15d4e20ea10a3d898cdcf1d2c41dba"),
            parse("0x1f5714cc13f8a33f4b32a07f1a85409de2d1b8d353aa269ca8f2bfa91071bd62"),
            parse("0x588c20c682f8b66c52049d910a4f0e7195dfd300a0d0cf3198b1383be1a4134a"),
        ];
        assert!(
            test_permutation_impl::<poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                1,
                parse("0x6cbdfb134e3be4606502e7a77faf60f863e8951a2fc4ed2a4d1c1154e76cb072")
            )
            .is_ok()
        );
        assert!(
            test_permutation_impl::<poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                2,
                parse("0xc822adf88e0162ff0321f5184887d77f5747d157d2ed74530e1cf87af913480d")
            )
            .is_ok()
        );
        assert!(
            test_permutation_impl::<poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                3,
                parse("0x721332b6ca26baf6d8ad8c8ca675e6895e0c066b020c81fe6d7766cf9e587f83")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x11de0b3702563747b0729abd2b93e720ec947ce067ca3f8c088b9829df1169d7"),
            parse("0x487ca87fa054ad960f6f571cf7aa6f5e075fd5ac3386e2e9fb22aeeb036fb1e1"),
            parse("0x28c775ddba9039531ee5deb0d4e2e6b1c9bb7d8be8da260a20d331bb041d0d38"),
            parse("0x10c93f76dca6f63507bc8ed1d2ea40647bb480ad43ac9922a3f10597b802f949"),
        ];
        assert!(
            test_permutation_impl::<poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                1,
                parse("0xa5eba6d8a98b0c08160cd211c6dd0ca165ce7e69dd5b73ced08178bd4c3e167c")
            )
            .is_ok()
        );
        assert!(
            test_permutation_impl::<poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                2,
                parse("0xa45a0ccc8a371133d2666365f34560bc5b0c6d1a42cfde3b3d403438a94718f2")
            )
            .is_ok()
        );
        assert!(
            test_permutation_impl::<poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                3,
                parse("0xae2fd35f72c75ea4db8bd7a0e19d970a3b90eee6c8e7032d3f0af1a988b0522d")
            )
            .is_ok()
        );
    }

    // TODO
}
