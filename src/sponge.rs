use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, rvar, var,
};

/// A hash with sponge construction.
#[derive(Debug, Default, Copy, Clone)]
pub struct Chip<const N: usize, P: PlonkChip<T, T>, const T: usize, const R: usize, const C: usize>
{
    permutation: P,
}

impl<const N: usize, P: PlonkChip<T, T>, const T: usize, const R: usize, const C: usize>
    Chip<N, P, T, R, C>
{
    fn build_absorb(
        &self,
        view: &mut impl CircuitView,
        state: [Option<Cell>; T],
        inputs: &mut impl Iterator<Item = Option<Cell>>,
    ) -> [Option<Cell>; T] {
        for i in 0..T {
            view.connect(state[i], view.cell(0, i).into());
        }
        for i in 0..R {
            match inputs.next() {
                Some(input) => {
                    view.connect(input, Some(view.cell(1, i)));
                }
                None => {
                    view.add_gate(1, var(i));
                }
            }
        }
        for i in 0..C {
            view.add_gate(1, var(R + i));
        }
        for i in 0..T {
            view.add_gate(1, rvar(i, -1) + rvar(i, 0) - rvar(i, 1));
        }
        std::array::from_fn(|i| Some(view.cell(2, i)))
    }

    fn witness_absorb(
        &self,
        view: &mut impl WitnessView,
        state: [CellOrUnconstrained; T],
        inputs: &mut impl Iterator<Item = CellOrUnconstrained>,
    ) -> [CellOrUnconstrained; T] {
        for i in 0..T {
            view.copy(state[i], view.cell(0, i));
        }
        for i in 0..R {
            match inputs.next() {
                Some(input) => {
                    view.copy(input, view.cell(1, i));
                }
                None => {
                    view.set(view.cell(1, i), Scalar::ZERO);
                }
            }
        }
        for i in 0..C {
            view.set(view.cell(1, R + i), Scalar::ZERO);
        }
        for i in 0..T {
            view.set(
                view.cell(2, i),
                view.get(view.cell(0, i)) + view.get(view.cell(1, i)),
            );
        }
        std::array::from_fn(|i| view.cell(2, i).into())
    }
}

impl<const N: usize, P: PlonkChip<T, T>, const T: usize, const R: usize, const C: usize>
    PlonkChip<N, R> for Chip<N, P, T, R, C>
{
    fn width(&self) -> usize {
        T * N.next_multiple_of(R) / R
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; N],
    ) -> Result<[Option<Cell>; R]> {
        let num_chunks = N.next_multiple_of(R) / R;
        assert!(num_chunks > 0);

        let mut state: [Option<Cell>; T] = std::array::from_fn(|i| {
            view.add_gate(0, var(i));
            Some(view.cell(0, i))
        });

        let mut input_it = inputs.iter().copied();

        state = self.build_absorb(view, state, &mut input_it);
        state = view.sub_chip(3, 0, &self.permutation, state)?;

        for c in 1..num_chunks {
            state = view
                .sub(0, c * T, T)
                .sub_fn(0, 0, T, |view| {
                    state = self.build_absorb(view, state, &mut input_it);
                })
                .sub_chip(0, 0, &self.permutation, state)?;
        }

        Ok(std::array::from_fn(|i| state[i]))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; N],
    ) -> Result<[CellOrUnconstrained; R]> {
        let num_chunks = N.next_multiple_of(R) / R;
        assert!(num_chunks > 0);

        let mut state: [CellOrUnconstrained; T] = std::array::from_fn(|i| {
            let cell = view.cell(0, i);
            view.set(cell, Scalar::ZERO);
            cell.into()
        });

        let mut input_it = inputs.iter().copied();

        state = self.witness_absorb(view, state, &mut input_it);
        state = view.sub_chip(3, 0, &self.permutation, state)?;

        for c in 1..num_chunks {
            state = view
                .sub(0, c * T, T)
                .sub_fn(0, 0, T, |view| {
                    state = self.witness_absorb(view, state, &mut input_it);
                })
                .sub_chip(0, 0, &self.permutation, state)?;
        }

        Ok(std::array::from_fn(|i| state[i]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poseidon1::PermutationChip as Poseidon1PermutationChip;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon as poseidon1;

    fn test_hash_v1<
        const N: usize,
        Cfg: poseidon1::Config<Scalar, T, R, C>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [Scalar; N],
        expected_output: [Scalar; R],
        blowup_log2: usize,
    ) -> Result<()> {
        let num_chunks = N.next_multiple_of(R) / R;
        let chip = Chip::<N, Poseidon1PermutationChip<Cfg, T, R, C>, T, R, C>::default();
        assert_eq!(chip.width(), T * num_chunks);
        let mut builder = CircuitBuilder::default();
        let output = chip.build(&mut builder, std::array::from_fn(|_| None))?;
        builder.declare_public_rows([output[0].unwrap().row()]);
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 197);
        assert_eq!(circuit.degree_bound(), 256);
        assert_eq!(circuit.num_columns(), T * num_chunks);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 197);
        assert_eq!(witness.degree_bound(), 256);
        assert_eq!(witness.num_columns(), T * num_chunks);
        let output = chip.witness(&mut witness, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 256);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 256 << blowup_log2);
        let public_inputs = circuit.verify(&proof, options)?;
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
    fn test_hash_v1_t3_1() {
        let inputs = [from_const(42)];
        let outputs = [
            parse_scalar("0x73952c443e4710be4a4c01e20046008b477f0d6fef5d87409cdebc4cdff3490c"),
            parse_scalar("0x05bf595cdacac4f9eba8679b69dcde4eeeca6db242005bf6b923fde28ea88a46"),
        ];
        assert!(test_hash_v1::<1, poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 1).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_2() {
        let inputs = [from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x28935bd3eba75f7b2d4f62babbd4e907b1ffcc28f73d1cae33654441a8a84023"),
            parse_scalar("0x339d0e485d8fdfb8c3391182d457fa3e73f043f566af1463ab05e57045122519"),
        ];
        assert!(test_hash_v1::<2, poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 1).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_3() {
        let inputs = [from_const(3), from_const(4), from_const(5)];
        let outputs = [
            parse_scalar("0x2bfc323795d99f44817eaa143a7db00103ff1eae1bd67ee3ab3f5a1006c7695d"),
            parse_scalar("0x5e7468521c84b23259b813d193017a2b3c7813ce82e94ce4cc74a8c527db0923"),
        ];
        assert!(test_hash_v1::<3, poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 1).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_4() {
        let inputs = [from_const(6), from_const(7), from_const(8), from_const(9)];
        let outputs = [
            parse_scalar("0x06ea9f66eddb8f036b0d6201dcf6a8c610b8aca9371e2bfc7fbd1deb1e5bb158"),
            parse_scalar("0x0bc4c477fdeee23bf2f139b12c2ea927d145f298e6204255cbad8461af9150c6"),
        ];
        assert!(test_hash_v1::<4, poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 1).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_5() {
        let inputs = [
            from_const(10),
            from_const(11),
            from_const(12),
            from_const(13),
            from_const(14),
        ];
        let outputs = [
            parse_scalar("0x05ae2c9b2bdbb5a64d4e838bd96b0b4c2366fc6d3cee4309793e01dfd2a589d1"),
            parse_scalar("0x67de663ef4d5db733c68cae13b6bb28aa97d0fc904dccdfa80f4c9fae36f51d0"),
        ];
        assert!(test_hash_v1::<5, poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 1).is_ok());
    }

    // TODO
}
