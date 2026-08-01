use crate::poseidon1::{
    PermutationChipER as Poseidon1PermutationChipER,
    PermutationChipHW as Poseidon1PermutationChipHW,
    PermutationChipIR as Poseidon1PermutationChipIR,
};
use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, rvar, var,
};
use starkom_poseidon as poseidon1;

/// A hash with sponge construction.
///
/// This chip implements only the sponge construction, while the PRP is outsourced to the `P` chip
/// specified in the generic arguments.
///
/// `N` is the number of input scalars, `T` is the state vector size, `R` is the ingestion rate, and
/// `C` is the number of capacity slots.
///
/// `T` must be equal to `R + C`. `N` can be any number and will be split in chunks of `R` scalars
/// each; if `N` is not a multiple of `R` the last chunk will be padded with zeros.
#[derive(Debug, Default, Copy, Clone)]
pub struct Chip<P: PlonkChip<T, T>, const T: usize, const R: usize, const C: usize, const N: usize>
{
    permutation: P,
}

impl<P: PlonkChip<T, T>, const T: usize, const R: usize, const C: usize, const N: usize>
    Chip<P, T, R, C, N>
{
    const ABSORB_HEIGHT: usize = 3;

    pub fn num_chunks() -> usize {
        N.next_multiple_of(R) / R
    }

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

impl<const T: usize, const R: usize, const C: usize, const N: usize> PlonkChip<N, R>
    for Chip<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig<T>, T>, T, R, C, N>
where
    poseidon1::BlueSkyConfig<T>: poseidon1::Config<Scalar, T>,
{
    fn width(&self) -> usize {
        self.permutation.width() * Self::num_chunks()
    }

    fn height(&self) -> usize {
        Self::ABSORB_HEIGHT + self.permutation.height()
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; N],
    ) -> Result<[Option<Cell>; R]> {
        let mut state: [Option<Cell>; T] = std::array::from_fn(|i| {
            view.add_gate(0, var(i));
            Some(view.cell(0, i))
        });

        let mut input_it = inputs.iter().copied();

        state = self.build_absorb(view, state, &mut input_it);
        state = view.sub_chip(Self::ABSORB_HEIGHT, 0, &self.permutation, state)?;

        let chunk_width = self.permutation.width();
        for c in 1..Self::num_chunks() {
            state = view
                .sub(0, c * chunk_width, chunk_width)
                .sub_fn(0, 0, T, |view| {
                    state = self.build_absorb(view, state, &mut input_it);
                })
                .sub_chip(Self::ABSORB_HEIGHT, 0, &self.permutation, state)?;
        }

        Ok(std::array::from_fn(|i| state[i]))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; N],
    ) -> Result<[CellOrUnconstrained; R]> {
        let mut state: [CellOrUnconstrained; T] = std::array::from_fn(|i| {
            let cell = view.cell(0, i);
            view.set(cell, Scalar::ZERO);
            cell.into()
        });

        let mut input_it = inputs.iter().copied();

        state = self.witness_absorb(view, state, &mut input_it);
        state = view.sub_chip(3, 0, &self.permutation, state)?;

        let chunk_width = self.permutation.width();
        for c in 1..Self::num_chunks() {
            state = view
                .sub(0, c * chunk_width, chunk_width)
                .sub_fn(0, 0, chunk_width, |view| {
                    state = self.witness_absorb(view, state, &mut input_it);
                })
                .sub_chip(Self::ABSORB_HEIGHT, 0, &self.permutation, state)?;
        }

        Ok(std::array::from_fn(|i| state[i]))
    }
}

/// T=3 Poseidon hash with [hard-wired round constants](`crate::poseidon1::RcModeHardWired`).
///
/// Compared to [`Poseidon1ChipT3IR`] and [`Poseidon1ChipT3ER`] this chip uses half the columns but
/// creates many more gate constraints (roughly 40 times more).
pub type Poseidon1ChipT3HW<const N: usize> =
    Chip<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>, 3, 2, 1, N>;

/// T=4 Poseidon hash with [hard-wired round constants](`crate::poseidon1::RcModeHardWired`).
///
/// Compared to [`Poseidon1ChipT4IR`] and [`Poseidon1ChipT4ER`] this chip uses half the columns but
/// creates many more gate constraints (roughly 40 times more).
pub type Poseidon1ChipT4HW<const N: usize> =
    Chip<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>, 4, 3, 1, N>;

/// T=3 Poseidon hash with
/// [internal ROM storage for round constants](`crate::poseidon1::RcModeInternalRom`).
///
/// The first chunk uses a [`PermutationChipIR` chip](`crate::poseidon1::PermutationChipIR`), while
/// subsequent chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`)
/// referring to the ROM of the first chip.
pub type Poseidon1ChipT3IR<const N: usize> =
    Chip<Poseidon1PermutationChipIR<poseidon1::BlueSkyConfig3, 3>, 3, 2, 1, N>;

/// T=4 Poseidon hash with
/// [internal ROM storage for round constants](`crate::poseidon1::RcModeInternalRom`).
///
/// The first chunk uses a [`PermutationChipIR` chip](`crate::poseidon1::PermutationChipIR`), while
/// subsequent chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`)
/// referring to the ROM of the first chip.
pub type Poseidon1ChipT4IR<const N: usize> =
    Chip<Poseidon1PermutationChipIR<poseidon1::BlueSkyConfig4, 4>, 4, 3, 1, N>;

/// T=3 Poseidon hash with
/// [external ROM storage for round constants](`crate::poseidon1::RcModeExternalRom`).
///
/// All chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`) referring to a
/// user-specified ROM area for the round constants, so a `Poseidon1ChipT3ER` can only be placed in
/// the circuit if a [`Poseidon1ChipT3IR`] is also present.
pub type Poseidon1ChipT3ER<'a, const N: usize> =
    Chip<Poseidon1PermutationChipER<'a, poseidon1::BlueSkyConfig3, 3>, 3, 2, 1, N>;

/// T=4 Poseidon hash with
/// [external ROM storage for round constants](`crate::poseidon1::RcModeExternalRom`).
///
/// All chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`) referring to a
/// user-specified ROM area for the round constants, so a `Poseidon1ChipT4ER` can only be placed in
/// the circuit if a [`Poseidon1ChipT4IR`] is also present.
pub type Poseidon1ChipT4ER<'a, const N: usize> =
    Chip<Poseidon1PermutationChipER<'a, poseidon1::BlueSkyConfig4, 4>, 4, 3, 1, N>;

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};

    const BLOWUP_LOG2: usize = 3;

    fn test_hash_v1<
        P: PlonkChip<T, T> + Default,
        const T: usize,
        const R: usize,
        const C: usize,
        const N: usize,
    >(
        inputs: [Scalar; N],
        expected_output: [Scalar; R],
    ) -> Result<()>
    where
        Chip<P, T, R, C, N>: PlonkChip<N, R>,
    {
        let num_chunks = N.next_multiple_of(R) / R;
        let chip = Chip::<P, T, R, C, N>::default();
        assert_eq!(chip.width(), T * num_chunks);
        assert_eq!(chip.height(), 197);
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
        let options = ProvingOptions {
            blowup_log2: BLOWUP_LOG2,
        };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 256);
        assert_eq!(proof.blowup_log2(), BLOWUP_LOG2);
        assert_eq!(proof.extended_domain_size(), 256 << BLOWUP_LOG2);
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
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>;
        assert!(test_hash_v1::<P, 3, 2, 1, 1>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_2() {
        let inputs = [from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x28935bd3eba75f7b2d4f62babbd4e907b1ffcc28f73d1cae33654441a8a84023"),
            parse_scalar("0x339d0e485d8fdfb8c3391182d457fa3e73f043f566af1463ab05e57045122519"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>;
        assert!(test_hash_v1::<P, 3, 2, 1, 2>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_3() {
        let inputs = [from_const(3), from_const(4), from_const(5)];
        let outputs = [
            parse_scalar("0x2bfc323795d99f44817eaa143a7db00103ff1eae1bd67ee3ab3f5a1006c7695d"),
            parse_scalar("0x5e7468521c84b23259b813d193017a2b3c7813ce82e94ce4cc74a8c527db0923"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>;
        assert!(test_hash_v1::<P, 3, 2, 1, 3>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_4() {
        let inputs = [from_const(6), from_const(7), from_const(8), from_const(9)];
        let outputs = [
            parse_scalar("0x06ea9f66eddb8f036b0d6201dcf6a8c610b8aca9371e2bfc7fbd1deb1e5bb158"),
            parse_scalar("0x0bc4c477fdeee23bf2f139b12c2ea927d145f298e6204255cbad8461af9150c6"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>;
        assert!(test_hash_v1::<P, 3, 2, 1, 4>(inputs, outputs).is_ok());
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
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>;
        assert!(test_hash_v1::<P, 3, 2, 1, 5>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_1() {
        let inputs = [from_const(42)];
        let outputs = [
            parse_scalar("0x2fdb574b84cca8f2c657ea588d8812bafbba305b7a9933728753de0fcf104c40"),
            parse_scalar("0x732f901b286e0f3575ab52e19494406c38f3db3e06169143f4c0369b3ba58ed9"),
            parse_scalar("0x24c8327a61a3bd811e04b11107609bd91b8916ab5cf53fe927edaa27a9e8d5da"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>;
        assert!(test_hash_v1::<P, 4, 3, 1, 1>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_2() {
        let inputs = [from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x33eaaa53f69ea75566e04bcb9318f965d5e74b68663bb4a09adfeeae27c752f4"),
            parse_scalar("0x292ad2994473be89dbfec5185888d85924bfa0f64b3be556609bbde3bad4360c"),
            parse_scalar("0x3f972105e69fcceafe6ce580dab417c50a34316d2de43d73a79f861ef55ca87a"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>;
        assert!(test_hash_v1::<P, 4, 3, 1, 2>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_3() {
        let inputs = [from_const(3), from_const(4), from_const(5)];
        let outputs = [
            parse_scalar("0x5220b264d93b85d22b4eb5a19c53ebfd08e1702e00dc76de14603165663006ea"),
            parse_scalar("0x664ca128f4f6f225f282a671b522c267389f30f01d858757a1f029941510d8ec"),
            parse_scalar("0x38ef442cd0ce47da5e7fdd912edfc2a95a36409b142fd0f94545267af135bcfa"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>;
        assert!(test_hash_v1::<P, 4, 3, 1, 3>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_4() {
        let inputs = [from_const(6), from_const(7), from_const(8), from_const(9)];
        let outputs = [
            parse_scalar("0x6cef40d837aeb6183356cf40d9818bc0ee109c557b17bffd80ab2905e4e2292f"),
            parse_scalar("0x688cf6e7f2aba6c399bc3253ce3827f7a003f8170fe679cbcc2b37e9ba65211e"),
            parse_scalar("0x1d14218c5f5ae32b4fc20b250b52ad8ec96a77627a6c103c8ecf3919290d6239"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>;
        assert!(test_hash_v1::<P, 4, 3, 1, 4>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_5() {
        let inputs = [
            from_const(10),
            from_const(11),
            from_const(12),
            from_const(13),
            from_const(14),
        ];
        let outputs = [
            parse_scalar("0x22c96d13097aa3b4782f9d2580dc2295378f87c85aaed5f47ee3f8b036faf8ee"),
            parse_scalar("0x0bb3130cba6d1aa9cd4ac577dd503905305ce7ccc08d04ec15d9a9700eb747a1"),
            parse_scalar("0x1cc1f59d0c8b31f60c5b10478b28db466bdcdefda0e8da296d96d5529177d621"),
        ];
        type P = Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>;
        assert!(test_hash_v1::<P, 4, 3, 1, 5>(inputs, outputs).is_ok());
    }
}
