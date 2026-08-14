use crate::poseidon1::{
    PermutationChipER as Poseidon1PermutationChipER,
    PermutationChipHW as Poseidon1PermutationChipHW,
    PermutationChipIR as Poseidon1PermutationChipIR,
};
use crate::poseidon2::{
    PermutationChipER as Poseidon2PermutationChipER,
    PermutationChipHW as Poseidon2PermutationChipHW,
    PermutationChipIR as Poseidon2PermutationChipIR,
};
use crate::sponge::internal::PrpMode;
use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, rvar, var,
};
use starkom_poseidon as poseidon1;
use starkom_poseidon2 as poseidon2;
use std::fmt::Debug;

mod internal {
    use super::*;

    /// Implements various layouts for the PRP (pseudo-random permutation) used by the sponge.
    ///
    /// The simplest mode is the [self-contained](`PrpModeSelfContained`) one, but it's also the
    /// most expensive in terms of number of gate constraints if you're using Poseidon because the
    /// self-contained Poseidon permutations (for both version 1 and 2) use hard-wired round
    /// constants, causing all gate constraints to be different from each other. See
    /// [`crate::poseidon1::RcModeHardWired`], [`crate::poseidon1::RcModeInternalRom`], and
    /// [`crate::poseidon1::RcModeExternalRom`] (or [`crate::poseidon2::RcModeHardWired`],
    /// [`crate::poseidon2::RcModeInternalRom`], and [`crate::poseidon2::RcModeExternalRom`]) for
    /// more information about the tradeoffs.
    pub trait PrpMode<const T: usize>: Debug + Copy + Clone {
        fn prp_width(&self) -> usize;
        fn prp_height(&self) -> usize;

        /// Builds the i-th permutation stage in the sponge.
        fn build_permute(
            &self,
            view: &mut impl CircuitView,
            index: usize,
            inputs: [Option<Cell>; T],
        ) -> Result<[Option<Cell>; T]>;

        /// Witnesses the i-th permutation stage in the sponge.
        fn witness_permute(
            &self,
            view: &mut impl WitnessView,
            index: usize,
            inputs: [CellOrUnconstrained; T],
        ) -> Result<[CellOrUnconstrained; T]>;
    }
}

/// Marker trait for self-contained PRP chips such as [`crate::poseidon1::PermutationChipHW`] and
/// [`crate::poseidon2::PermutationChipHW`].
pub trait SelfContainedPermutationChip<const T: usize>:
    PlonkChip<T, T> + Debug + Default + Copy + Clone
{
}

impl<Cfg: poseidon1::Config<Scalar, T>, const T: usize> SelfContainedPermutationChip<T>
    for Poseidon1PermutationChipHW<Cfg, T>
{
}

impl<Cfg: poseidon2::Config<Scalar, T>, const T: usize> SelfContainedPermutationChip<T>
    for Poseidon2PermutationChipHW<Cfg, T>
{
}

#[derive(Debug, Default, Copy, Clone)]
pub struct PrpModeSelfContained<P: PlonkChip<T, T> + Debug + Default + Copy + Clone, const T: usize>
{
    permutation: P,
}

impl<P: PlonkChip<T, T> + Debug + Default + Copy + Clone, const T: usize> internal::PrpMode<T>
    for PrpModeSelfContained<P, T>
{
    fn prp_width(&self) -> usize {
        self.permutation.width()
    }

    fn prp_height(&self) -> usize {
        self.permutation.height()
    }

    fn build_permute(
        &self,
        view: &mut impl CircuitView,
        _index: usize,
        inputs: [Option<Cell>; T],
    ) -> Result<[Option<Cell>; T]> {
        view.sub_chip(0, 0, &self.permutation, inputs)
    }

    fn witness_permute(
        &self,
        view: &mut impl WitnessView,
        _index: usize,
        inputs: [CellOrUnconstrained; T],
    ) -> Result<[CellOrUnconstrained; T]> {
        view.sub_chip(0, 0, &self.permutation, inputs)
    }
}

/// Marker trait for Poseidon permutation chips with internal ROMs
/// ([`crate::poseidon1::PermutationChipIR`] and [`crate::poseidon2::PermutationChipIR`]).
pub trait PoseidonPermutationChipIR<const T: usize>:
    PlonkChip<T, T> + Debug + Default + Copy + Clone
{
}

impl<C: poseidon1::Config<Scalar, T>, const T: usize> PoseidonPermutationChipIR<T>
    for Poseidon1PermutationChipIR<C, T>
{
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> PoseidonPermutationChipIR<T>
    for Poseidon2PermutationChipIR<C, T>
{
}

/// Marker trait for Poseidon permutation chips with external ROMs
/// ([`crate::poseidon1::PermutationChipER`] and [`crate::poseidon2::PermutationChipER`]).
pub trait PoseidonPermutationChipER<const T: usize>:
    PlonkChip<T, T> + Debug + Copy + Clone
{
    fn make_new(ir_chip_row_offset: isize, ir_chip_column_offset: isize) -> Self;
}

impl<C: poseidon1::Config<Scalar, T>, const T: usize> PoseidonPermutationChipER<T>
    for Poseidon1PermutationChipER<C, T>
{
    fn make_new(ir_chip_row_offset: isize, ir_chip_column_offset: isize) -> Self {
        Poseidon1PermutationChipER::new(ir_chip_row_offset, ir_chip_column_offset)
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> PoseidonPermutationChipER<T>
    for Poseidon2PermutationChipER<C, T>
{
    fn make_new(ir_chip_row_offset: isize, ir_chip_column_offset: isize) -> Self {
        Poseidon2PermutationChipER::new(ir_chip_row_offset, ir_chip_column_offset)
    }
}

#[derive(Debug, Copy, Clone)]
pub struct PrpModeInternalRom<
    IR: PoseidonPermutationChipIR<T>,
    ER: PoseidonPermutationChipER<T>,
    const T: usize,
    const N: usize,
> {
    ir_chip: IR,

    /// NOTE: we need `N` permutation chips in total of which 1 IR and N-1 ERs, but Rust doesn't
    /// allow NTT variable expressions so we can't instantiate `[ER; N-1]`. We have to instantiate
    /// `[ER; N]`, so there's an extra wasted ER chip instance. In order to keep the offset formulas
    /// simple we pick slot #0 as the extra wasted one, so permutation #1 is the ER chip at slot #1,
    /// permutation #2 is the ER chip at slot #2, and so on. Permutation #0 is of course the IR chip
    /// above.
    er_chip: [ER; N],
}

impl<
    IR: PoseidonPermutationChipIR<T>,
    ER: PoseidonPermutationChipER<T>,
    const T: usize,
    const N: usize,
> Default for PrpModeInternalRom<IR, ER, T, N>
{
    fn default() -> Self {
        let ir_chip = IR::default();
        let stage_width = ir_chip.width() as isize;
        Self {
            ir_chip,
            er_chip: std::array::from_fn(|i| ER::make_new(0, i as isize * -stage_width)),
        }
    }
}

impl<
    IR: PoseidonPermutationChipIR<T>,
    ER: PoseidonPermutationChipER<T>,
    const T: usize,
    const N: usize,
> internal::PrpMode<T> for PrpModeInternalRom<IR, ER, T, N>
{
    fn prp_width(&self) -> usize {
        self.ir_chip.width()
    }

    fn prp_height(&self) -> usize {
        self.ir_chip.height()
    }

    fn build_permute(
        &self,
        view: &mut impl CircuitView,
        index: usize,
        inputs: [Option<Cell>; T],
    ) -> Result<[Option<Cell>; T]> {
        if index > 0 {
            view.sub_chip(0, 0, &self.er_chip[index], inputs)
        } else {
            view.sub_chip(0, 0, &self.ir_chip, inputs)
        }
    }

    fn witness_permute(
        &self,
        view: &mut impl WitnessView,
        index: usize,
        inputs: [CellOrUnconstrained; T],
    ) -> Result<[CellOrUnconstrained; T]> {
        if index > 0 {
            view.sub_chip(0, 0, &self.er_chip[index], inputs)
        } else {
            view.sub_chip(0, 0, &self.ir_chip, inputs)
        }
    }
}

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
pub struct Chip<
    M: internal::PrpMode<T>,
    const T: usize,
    const R: usize,
    const C: usize,
    const N: usize,
> {
    mode: M,
}

impl<M: internal::PrpMode<T>, const T: usize, const R: usize, const C: usize, const N: usize>
    Chip<M, T, R, C, N>
{
    pub const ABSORB_HEIGHT: usize = 3;

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
                view.get_at(view.cell(0, i)) + view.get_at(view.cell(1, i)),
            );
        }
        std::array::from_fn(|i| view.cell(2, i).into())
    }
}

impl<
    P: SelfContainedPermutationChip<T>,
    const T: usize,
    const R: usize,
    const C: usize,
    const N: usize,
> PlonkChip<N, R> for Chip<PrpModeSelfContained<P, T>, T, R, C, N>
{
    fn width(&self) -> usize {
        self.mode.prp_width() * Self::num_chunks()
    }

    fn height(&self) -> usize {
        Self::ABSORB_HEIGHT + self.mode.prp_height()
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
        let chunk_width = self.mode.prp_width();
        let chunk_height = Self::ABSORB_HEIGHT + self.mode.prp_height();

        view.sub(0, 0, chunk_width.into(), chunk_height.into())
            .sub_fn(0, 0, None, None, |view| {
                state = self.build_absorb(view, state, &mut input_it)
            })
            .sub_fn_or(Self::ABSORB_HEIGHT, 0, None, None, |view| {
                state = self.mode.build_permute(view, 0, state)?;
                Ok(())
            })?;

        for c in 1..Self::num_chunks() {
            view.sub(0, c * chunk_width, chunk_width.into(), chunk_height.into())
                .sub_fn(0, 0, None, None, |view| {
                    state = self.build_absorb(view, state, &mut input_it);
                })
                .sub_fn_or(Self::ABSORB_HEIGHT, 0, None, None, |view| {
                    state = self.mode.build_permute(view, c, state)?;
                    Ok(())
                })?;
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
        let chunk_width = self.mode.prp_width();
        let chunk_height = Self::ABSORB_HEIGHT + self.mode.prp_height();

        view.sub(0, 0, chunk_width.into(), chunk_height.into())
            .sub_fn(0, 0, None, None, |view| {
                state = self.witness_absorb(view, state, &mut input_it)
            })
            .sub_fn_or(Self::ABSORB_HEIGHT, 0, None, None, |view| {
                state = self.mode.witness_permute(view, 0, state)?;
                Ok(())
            })?;

        for c in 1..Self::num_chunks() {
            view.sub(0, c * chunk_width, chunk_width.into(), chunk_height.into())
                .sub_fn(0, 0, None, None, |view| {
                    state = self.witness_absorb(view, state, &mut input_it);
                })
                .sub_fn_or(Self::ABSORB_HEIGHT, 0, None, None, |view| {
                    state = self.mode.witness_permute(view, c, state)?;
                    Ok(())
                })?;
        }

        Ok(std::array::from_fn(|i| state[i]))
    }
}

impl<
    IR: PoseidonPermutationChipIR<T>,
    ER: PoseidonPermutationChipER<T>,
    const T: usize,
    const R: usize,
    const C: usize,
    const N: usize,
> PlonkChip<N, R> for Chip<PrpModeInternalRom<IR, ER, T, N>, T, R, C, N>
{
    fn width(&self) -> usize {
        self.mode.prp_width() * Self::num_chunks()
    }

    fn height(&self) -> usize {
        Self::ABSORB_HEIGHT + self.mode.prp_height()
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; N],
    ) -> Result<[Option<Cell>; R]> {
        // TODO
        todo!()
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; N],
    ) -> Result<[CellOrUnconstrained; R]> {
        // TODO
        todo!()
    }
}

/// T=3 Poseidon hash with [hard-wired round constants](`crate::poseidon1::RcModeHardWired`).
///
/// Compared to [`Poseidon1ChipT3IR`] and [`Poseidon1ChipT3ER`] this chip uses half the columns but
/// creates many more gate constraints (roughly 40 times more).
pub type Poseidon1ChipT3HW<const N: usize> = Chip<
    PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>, 3>,
    3,
    2,
    1,
    N,
>;

/// T=4 Poseidon hash with [hard-wired round constants](`crate::poseidon1::RcModeHardWired`).
///
/// Compared to [`Poseidon1ChipT4IR`] and [`Poseidon1ChipT4ER`] this chip uses half the columns but
/// creates many more gate constraints (roughly 40 times more).
pub type Poseidon1ChipT4HW<const N: usize> = Chip<
    PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>, 4>,
    4,
    3,
    1,
    N,
>;

/// T=3 Poseidon2 hash with [hard-wired round constants](`crate::poseidon2::RcModeHardWired`).
pub type Poseidon2ChipT3HW<const N: usize> = Chip<
    PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig3, 3>, 3>,
    3,
    2,
    1,
    N,
>;

/// T=4 Poseidon2 hash with [hard-wired round constants](`crate::poseidon2::RcModeHardWired`).
pub type Poseidon2ChipT4HW<const N: usize> = Chip<
    PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig4, 4>, 4>,
    4,
    3,
    1,
    N,
>;

/// T=3 Poseidon hash with
/// [internal ROM storage for round constants](`crate::poseidon1::RcModeInternalRom`).
///
/// The first chunk uses a [`PermutationChipIR` chip](`crate::poseidon1::PermutationChipIR`), while
/// subsequent chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`)
/// referring to the ROM of the first chip.
pub type Poseidon1ChipT3IR<const N: usize> = Chip<
    PrpModeInternalRom<
        Poseidon1PermutationChipIR<poseidon1::BlueSkyConfig3, 3>,
        Poseidon1PermutationChipER<poseidon1::BlueSkyConfig3, 3>,
        3,
        N,
    >,
    3,
    2,
    1,
    N,
>;

/// T=4 Poseidon hash with
/// [internal ROM storage for round constants](`crate::poseidon1::RcModeInternalRom`).
///
/// The first chunk uses a [`PermutationChipIR` chip](`crate::poseidon1::PermutationChipIR`), while
/// subsequent chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`)
/// referring to the ROM of the first chip.
pub type Poseidon1ChipT4IR<const N: usize> = Chip<
    PrpModeInternalRom<
        Poseidon1PermutationChipIR<poseidon1::BlueSkyConfig4, 4>,
        Poseidon1PermutationChipER<poseidon1::BlueSkyConfig4, 4>,
        4,
        N,
    >,
    4,
    3,
    1,
    N,
>;

/// T=3 Poseidon2 hash with
/// [internal ROM storage for round constants](`crate::poseidon2::RcModeInternalRom`).
pub type Poseidon2ChipT3IR<const N: usize> = Chip<
    PrpModeInternalRom<
        Poseidon2PermutationChipIR<poseidon2::BlueSkyConfig3, 3>,
        Poseidon2PermutationChipER<poseidon2::BlueSkyConfig3, 3>,
        3,
        N,
    >,
    3,
    2,
    1,
    N,
>;

/// T=4 Poseidon2 hash with
/// [internal ROM storage for round constants](`crate::poseidon2::RcModeInternalRom`).
pub type Poseidon2ChipT4IR<const N: usize> = Chip<
    PrpModeInternalRom<
        Poseidon2PermutationChipIR<poseidon2::BlueSkyConfig4, 4>,
        Poseidon2PermutationChipER<poseidon2::BlueSkyConfig4, 4>,
        4,
        N,
    >,
    4,
    3,
    1,
    N,
>;

// /// T=3 Poseidon hash with
// /// [external ROM storage for round constants](`crate::poseidon1::RcModeExternalRom`).
// ///
// /// All chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`) referring to a
// /// user-specified ROM area for the round constants, so a `Poseidon1ChipT3ER` can only be placed in
// /// the circuit if a [`Poseidon1ChipT3IR`] is also present.
// pub type Poseidon1ChipT3ER<const N: usize> =
//     Chip<Poseidon1PermutationChipER<poseidon1::BlueSkyConfig3, 3>, 3, 2, 1, N>;

// /// T=4 Poseidon hash with
// /// [external ROM storage for round constants](`crate::poseidon1::RcModeExternalRom`).
// ///
// /// All chunks use [`PermutationChipER` chips](`crate::poseidon1::PermutationChipER`) referring to a
// /// user-specified ROM area for the round constants, so a `Poseidon1ChipT4ER` can only be placed in
// /// the circuit if a [`Poseidon1ChipT4IR`] is also present.
// pub type Poseidon1ChipT4ER<const N: usize> =
//     Chip<Poseidon1PermutationChipER<poseidon1::BlueSkyConfig4, 4>, 4, 3, 1, N>;

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};

    const BLOWUP_LOG2: usize = 3;

    fn test_hash<
        M: internal::PrpMode<T>,
        const T: usize,
        const R: usize,
        const C: usize,
        const N: usize,
    >(
        inputs: [Scalar; N],
        expected_output: [Scalar; R],
    ) -> Result<()>
    where
        Chip<M, T, R, C, N>: PlonkChip<N, R> + Default,
    {
        let num_chunks = N.next_multiple_of(R) / R;
        let chip = Chip::<M, T, R, C, N>::default();
        assert_eq!(chip.width(), T * num_chunks);
        assert_eq!(chip.height(), 197);
        let mut builder = CircuitBuilder::default();
        let output = builder.sub_chip(0, 0, &chip, std::array::from_fn(|_| None))?;
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
        let output = witness.sub_chip(0, 0, &chip, inputs.map(|input| input.into()))?;
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
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 1>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_2() {
        let inputs = [from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x28935bd3eba75f7b2d4f62babbd4e907b1ffcc28f73d1cae33654441a8a84023"),
            parse_scalar("0x339d0e485d8fdfb8c3391182d457fa3e73f043f566af1463ab05e57045122519"),
        ];
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 2>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_3() {
        let inputs = [from_const(3), from_const(4), from_const(5)];
        let outputs = [
            parse_scalar("0x2bfc323795d99f44817eaa143a7db00103ff1eae1bd67ee3ab3f5a1006c7695d"),
            parse_scalar("0x5e7468521c84b23259b813d193017a2b3c7813ce82e94ce4cc74a8c527db0923"),
        ];
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 3>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t3_4() {
        let inputs = [from_const(6), from_const(7), from_const(8), from_const(9)];
        let outputs = [
            parse_scalar("0x06ea9f66eddb8f036b0d6201dcf6a8c610b8aca9371e2bfc7fbd1deb1e5bb158"),
            parse_scalar("0x0bc4c477fdeee23bf2f139b12c2ea927d145f298e6204255cbad8461af9150c6"),
        ];
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 4>(inputs, outputs).is_ok());
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
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 5>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_1() {
        let inputs = [from_const(42)];
        let outputs = [
            parse_scalar("0x2fdb574b84cca8f2c657ea588d8812bafbba305b7a9933728753de0fcf104c40"),
            parse_scalar("0x732f901b286e0f3575ab52e19494406c38f3db3e06169143f4c0369b3ba58ed9"),
            parse_scalar("0x24c8327a61a3bd811e04b11107609bd91b8916ab5cf53fe927edaa27a9e8d5da"),
        ];
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 1>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_2() {
        let inputs = [from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x33eaaa53f69ea75566e04bcb9318f965d5e74b68663bb4a09adfeeae27c752f4"),
            parse_scalar("0x292ad2994473be89dbfec5185888d85924bfa0f64b3be556609bbde3bad4360c"),
            parse_scalar("0x3f972105e69fcceafe6ce580dab417c50a34316d2de43d73a79f861ef55ca87a"),
        ];
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 2>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_3() {
        let inputs = [from_const(3), from_const(4), from_const(5)];
        let outputs = [
            parse_scalar("0x5220b264d93b85d22b4eb5a19c53ebfd08e1702e00dc76de14603165663006ea"),
            parse_scalar("0x664ca128f4f6f225f282a671b522c267389f30f01d858757a1f029941510d8ec"),
            parse_scalar("0x38ef442cd0ce47da5e7fdd912edfc2a95a36409b142fd0f94545267af135bcfa"),
        ];
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 3>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v1_t4_4() {
        let inputs = [from_const(6), from_const(7), from_const(8), from_const(9)];
        let outputs = [
            parse_scalar("0x6cef40d837aeb6183356cf40d9818bc0ee109c557b17bffd80ab2905e4e2292f"),
            parse_scalar("0x688cf6e7f2aba6c399bc3253ce3827f7a003f8170fe679cbcc2b37e9ba65211e"),
            parse_scalar("0x1d14218c5f5ae32b4fc20b250b52ad8ec96a77627a6c103c8ecf3919290d6239"),
        ];
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 4>(inputs, outputs).is_ok());
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
        type M = PrpModeSelfContained<Poseidon1PermutationChipHW<poseidon1::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 5>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t3_1() {
        let inputs = [from_const(42)];
        let outputs = [
            parse_scalar("0x302e6d6d782c1367974698e051d9b55e18060b19393a4f0ac4b66f992bd5a5eb"),
            parse_scalar("0x26f778c0f82ffe3d4409ebb9d7e4611556ca89c6a3e1a77cf8b80528eb344777"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 1>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t3_2() {
        let inputs = [from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x2a24882111b586a835203bdeb7a97d8489e410eadf12a495624f49b729528873"),
            parse_scalar("0x69233d2461effb6b25dbec14086d466f3bf668ef2a38759fa5cb433bedf25778"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 2>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t3_3() {
        let inputs = [from_const(3), from_const(4), from_const(5)];
        let outputs = [
            parse_scalar("0x160be03feff499f1256ce2404ff9ee026fc378b6a91d434746bab98aafaecb63"),
            parse_scalar("0x14a259b91d964d8263af60fc1325c4874c68e8fd9caef509cc07622fc17718fe"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 3>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t3_4() {
        let inputs = [from_const(6), from_const(7), from_const(8), from_const(9)];
        let outputs = [
            parse_scalar("0x63d491b523ae737f62f117ef5affb8353996b67034ddaeb8586b574678ab440a"),
            parse_scalar("0x531109ee099551ccea55a61f6f7cab781bf0d3d0d0c4ba032476b65d1ebb9867"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 4>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t3_5() {
        let inputs = [
            from_const(10),
            from_const(11),
            from_const(12),
            from_const(13),
            from_const(14),
        ];
        let outputs = [
            parse_scalar("0x329255ad3db8a69a50a2a1f63fb4046d06d5bc6de30bf79bfe4138f4c93201df"),
            parse_scalar("0x6968c301186d76def97ee0d7bcc1f426b34df8f2e04a3afdeaa1acd8f9070d76"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig3, 3>, 3>;
        assert!(test_hash::<M, 3, 2, 1, 5>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t4_1() {
        let inputs = [from_const(42)];
        let outputs = [
            parse_scalar("0x109a9fd885b0047b036489dad6d0ca97749f6a9b21d9fc2c1cb7d25952e453a0"),
            parse_scalar("0x203e5346a31efe538f826a34e87c285ef6cfe0ce12a0316a25cbd4e2326abd29"),
            parse_scalar("0x4c53083849aedf3e11959d1dad010d2f1d2951adfdcce95f6a480666e63c5834"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 1>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t4_2() {
        let inputs = [from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x7c4e380d8a3935c0e8073420573f5b6aaf9ed2c727afc4da64f12401ab355faf"),
            parse_scalar("0x14b0dda71f3fb062cc99121629f080541891b8be0e65ba858906cf0b648042ac"),
            parse_scalar("0x5218f71044490008f3713824dfa6be57a708ad295ca8df9fb4176340d61fb681"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 2>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t4_3() {
        let inputs = [from_const(3), from_const(4), from_const(5)];
        let outputs = [
            parse_scalar("0x2582eca7bed4bca9d4326a9e2ca601e0b3779582bb5173318a4e19ab005e7495"),
            parse_scalar("0x0307db21c6063767f309fb09afcc4f250cbaed3f1d0870cd56d3276b18ced8d5"),
            parse_scalar("0x3024913685a18187c7ae50f1dcabe3ea2acb407fc3abd5d107bd79f3dbd2e90c"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 3>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t4_4() {
        let inputs = [from_const(6), from_const(7), from_const(8), from_const(9)];
        let outputs = [
            parse_scalar("0x6b13720a0ebd34f13327023c0232a3a3421f88d50b627bacfd114491ae48bfaa"),
            parse_scalar("0x235d36a318fbc8175bce6613d8b8812a1d8ab17c56a70565f8eb1253a248f5d0"),
            parse_scalar("0x379747ade21215413ccf1e1d91c7f367e1d5d8cd3e52b24da10e080dc8b25c43"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 4>(inputs, outputs).is_ok());
    }

    #[test]
    fn test_hash_v2_t4_5() {
        let inputs = [
            from_const(10),
            from_const(11),
            from_const(12),
            from_const(13),
            from_const(14),
        ];
        let outputs = [
            parse_scalar("0x4f07a42cf3cd73f35eeb9b42bff06b11e1c7ebe0fd8f65b7fab0dd5d551f1c6c"),
            parse_scalar("0x2507a1f641f6bab3bb2cc40cb6d14df149aeede849a53a4590d73b0d29af2d71"),
            parse_scalar("0x79f1da2d204aa97c4256d321055eac279959efeff119a32e106220ef69699d4f"),
        ];
        type M = PrpModeSelfContained<Poseidon2PermutationChipHW<poseidon2::BlueSkyConfig4, 4>, 4>;
        assert!(test_hash::<M, 4, 3, 1, 5>(inputs, outputs).is_ok());
    }
}
