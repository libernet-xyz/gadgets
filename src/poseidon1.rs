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
            blowup_log2,
            blind: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.num_columns(), chip.width());
        assert_eq!(circuit.get_max_gate_degree(), C::alpha());
        let degree_bound = circuit.degree_bound();
        assert!(degree_bound <= 256);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), degree_bound);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<F>>(witness)?;
        assert_eq!(proof.degree_bound(), degree_bound);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), degree_bound << blowup_log2);
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
                parse("0xe064580fd786e7c56400181b1ee2c0229a32ec81f03289e95c914692ddaa4b44")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT3X5, 3>(
                inputs,
                outputs,
                2,
                parse("0xd9732b132ddb003e158479ec29c27393cc52788db512fa22987c5a2af5324501")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT3X5, 3>(
                inputs,
                outputs,
                3,
                parse("0x907529d1bb197d14b0d04495f4d292017b7f27450f260b321ebd2f1aa2cf3237")
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
                parse("0x9d283d00a407f7cc83745e9bb3c8d753b9f3d390c7e950b6cc8304429ca03e57")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT4X5, 4>(
                inputs,
                outputs,
                2,
                parse("0xf7b5e613e6fc926115d01a0549a29dadfd9d8cce92c883ea8f3396ee84bcf27f")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<BS, poseidon1::BlueSkyConfigT4X5, 4>(
                inputs,
                outputs,
                3,
                parse("0x4ba65c9f27631d2248e96286d44ab75ded518fa6a86e128e559e46b13be57941")
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
                parse("0x8337109609c02a68ee9a94503731841f85aae0cfb9e508eb4ceaff096709fc7d")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                2,
                parse("0xb205246f79389dcc2a9c9f0dbfcfed252cf7ce3bf39e2327ae505ac86acee8a0")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig3, 3>(
                inputs,
                outputs,
                3,
                parse("0x9813397aafe327c063ddf5d579771aa73bffd10083455086800855b84a6d42cc")
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
                parse("0xd94c3d97805b743004a60ca4c233ec737d99ab9fe96563b0a9379e87726391cb")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                2,
                parse("0xb288660aa911163029521a4bfde74c63724e78b3ed2800dc22e674c5300e6d1a")
            )
            .is_ok()
        );
        assert!(
            test_permutation256_impl::<SB, poseidon1::SchraderbrauConfig4, 4>(
                inputs,
                outputs,
                3,
                parse("0x7ec9036d87366b4bd4e0af9270564f09df5a3fe82a1c0f5fad291a871941bee8")
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
            blowup_log2,
            blind: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.num_columns(), chip.width());
        assert_eq!(circuit.get_max_gate_degree(), C::alpha());
        let degree_bound = circuit.degree_bound();
        assert!(degree_bound <= 1024);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), degree_bound);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness)?;
        assert_eq!(proof.degree_bound(), degree_bound);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), degree_bound << blowup_log2);
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
                parse("0x3af7c8adcd77c7142f11da66960357b1595a33a19fe5f0ca66d523013d3cd1ea")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0x4e151abfa9ea2b720f8fc19e4b2f870bae2abc2dbd4865bda7f753c5df1c4b54")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0x932192bd988f691129f4c1c19e78bc0089e2660ba7cab1b33df7473506f566b7")
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
                parse("0xbba319dc72544a48567389a2eff692d5dee291d8662221bba7b7fc0e7bbb7f53")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x745311090e041be5b2b52eae2181ca650b6dabab103e34a701780599df07a4a6")
            )
            .is_ok()
        );
        assert!(
            test_permutation64_impl::<GL, GL4, poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0x57857fffd3439733e5c49b9e0ba2435a21cd64920d9e493a2b1847938704e919")
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
            blowup_log2,
            blind: false,
        })?;
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.num_columns(), chip.width());
        assert_eq!(circuit.get_max_gate_degree(), C::alpha());
        let degree_bound = circuit.degree_bound();
        assert!(degree_bound <= 2048);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 1);
        assert_eq!(witness.degree_bound(), degree_bound);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness)?;
        assert_eq!(proof.degree_bound(), degree_bound);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), degree_bound << blowup_log2);
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
                parse("0x7a91806e225e9dcd82c86b808e2591c3e1c2b4867d8c817e47f677d17e92ce88")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig24, 24>(
                inputs,
                outputs,
                2,
                parse("0xce8daf5dcaacd272a22fb5a805aaa9b243aeb6aca8c45e4c945fdcb8281aaeeb")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig24, 24>(
                inputs,
                outputs,
                3,
                parse("0xf0ca7ad84c25e96e9e897980396da91f7a0d0269d9a0b40343529e0ef0112216")
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
                parse("0xea64ae6c3f834d4c182589adf0532c87dc51e9bfcb4a242b1037dd7cd3dc086a")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig32, 32>(
                inputs,
                outputs,
                2,
                parse("0x676f83f0699369e852f51e74c928b79d76a2c6f1410bb8549a456b8081c09163")
            )
            .is_ok()
        );
        assert!(
            test_permutation32_impl::<KB, KB8, poseidon1::KoalaBearConfig32, 32>(
                inputs,
                outputs,
                3,
                parse("0x358a1eaf7de64001f7086639fd7e21f28c8279bd1c8565b849411c0474da037d")
            )
            .is_ok()
        );
    }
}
