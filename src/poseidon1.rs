use anyhow::Result;
use starkom_ff::{Field256, PrimeField};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, make_const,
    rvar, var,
};
use starkom_poseidon as poseidon;
use std::fmt::{Debug, Formatter};
use std::marker::PhantomData;
use std::ops::Mul;

mod internal {
    use super::*;

    /// Encodes the operations that differ between round constant modes
    /// ([hard-wired](`RcModeHardWired`), [internal ROM](`RcModeInternalRom`),
    /// [external ROM](`RcModeExternalRom`)).
    pub trait RcMode<F: PrimeField, const T: usize>: Debug + Copy + Clone {
        fn width(&self) -> usize;

        fn build_full_round<G: Field256 + From<F>>(
            &self,
            view: &mut impl CircuitView<F, G>,
            round: usize,
        ) where
            F: Mul<G, Output = G>,
            G: Mul<F, Output = G>;

        fn witness_full_round(&self, view: &mut impl WitnessView<F>, round: usize);

        fn build_partial_round<G: Field256 + From<F>>(
            &self,
            view: &mut impl CircuitView<F, G>,
            round: usize,
        ) where
            F: Mul<G, Output = G>,
            G: Mul<F, Output = G>;

        fn witness_partial_round(&self, view: &mut impl WitnessView<F>, round: usize);
    }
}

/// Hard-wired RC mode for the Poseidon1 permutation.
///
/// In this mode the chip takes exactly `T` columns, with `T` being the state vector size, but uses
/// a different gate for every ARC layer. You should use this mode only if the elevated number of
/// gates is not a concern.
///
/// This mode is best suited for circuits that run a small number of hashes, such as preimage
/// knowledge proofs and zkMAC signatures.
pub struct RcModeHardWired<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> {
    _data: PhantomData<(F, C)>,
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Debug for RcModeHardWired<F, C, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RcModeHardWired").finish()
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Default
    for RcModeHardWired<F, C, T>
{
    fn default() -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Copy for RcModeHardWired<F, C, T> {}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Clone for RcModeHardWired<F, C, T> {
    fn clone(&self) -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> internal::RcMode<F, T>
    for RcModeHardWired<F, C, T>
{
    fn width(&self) -> usize {
        T
    }

    fn build_full_round<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA as isize;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.add_gate(
                0,
                (0..T)
                    .map(|j| make_const(m[i * T + j]) * ((var(j) + make_const(c[j])) ^ a))
                    .sum::<Constraint<F>>()
                    - rvar(i, 1),
            );
        }
    }

    fn witness_full_round(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| {
                        m[i * T + j] * (view.get_at(view.cell(0, j)) + c[j]).pow_small_vartime(a)
                    })
                    .sum(),
            );
        }
    }

    fn build_partial_round<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA as isize;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.add_gate(
                0,
                std::iter::once(make_const(m[i * T + 0]) * ((var(0) + make_const(c[0])) ^ a))
                    .chain((1..T).map(|j| make_const(m[i * T + j]) * (var(j) + make_const(c[j]))))
                    .sum::<Constraint<F>>()
                    - rvar(i, 1),
            );
        }
    }

    fn witness_partial_round(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                std::iter::once(
                    m[i * T] * (view.get_at(view.cell(0, 0)) + c[0]).pow_small_vartime(a),
                )
                .chain((1..T).map(|j| m[i * T + j] * (view.get_at(view.cell(0, j)) + c[j])))
                .sum(),
            );
        }
    }
}

/// Internal ROM mode for the [`PermutationChip`].
///
/// In this mode the Poseidon round constants are stored in separate columns so that we can reuse
/// the same gate constraints for the round themselves. By contrast, [`RcModeHardWired`] causes the
/// gates of all rounds to be different from each other due to the hard-wired constants, and that
/// results in an elevated number of gates and increased proving cost.
///
/// NOTE: the overall number of gates instantiated by this mode is comparable to the hard-wired mode
/// because we still need to constrain the cells of the ROM area to have the correct values, but the
/// advantage of this mode is that you can subsequently instantiate an arbitrary number of cheap,
/// external ROM mode chips (see [`RcModeExternalRom`]). In other words: N hard-wired permutation
/// chips pay the gate cost of the round constants N times, whereas 1 IR chip + (N-1) ER chips pay
/// that cost only once (but still achieve N permutations).
pub struct RcModeInternalRom<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> {
    _data: PhantomData<(F, C)>,
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Debug
    for RcModeInternalRom<F, C, T>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RcModeInternalRom").finish()
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Default
    for RcModeInternalRom<F, C, T>
{
    fn default() -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Copy for RcModeInternalRom<F, C, T> {}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Clone
    for RcModeInternalRom<F, C, T>
{
    fn clone(&self) -> Self {
        Self { _data: PhantomData }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> internal::RcMode<F, T>
    for RcModeInternalRom<F, C, T>
{
    fn width(&self) -> usize {
        T * 2
    }

    fn build_full_round<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA as isize;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.add_gate(0, var(T + i) - make_const(c[i]));
            view.add_gate(
                0,
                (0..T)
                    .map(|j| make_const(m[i * T + j]) * ((var(j) + rvar(T + j, 0)) ^ a))
                    .sum::<Constraint<F>>()
                    - rvar(i, 1),
            );
        }
    }

    fn witness_full_round(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(view.cell(0, T + i), c[i]);
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| {
                        m[i * T + j] * (view.get_at(view.cell(0, j)) + c[j]).pow_small_vartime(a)
                    })
                    .sum(),
            );
        }
    }

    fn build_partial_round<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA as isize;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.add_gate(0, var(T + i) - make_const(c[i]));
            view.add_gate(
                0,
                std::iter::once(make_const(m[i * T]) * ((var(0) + var(T)) ^ a))
                    .chain((1..T).map(|j| make_const(m[i * T + j]) * (var(j) + var(T + j))))
                    .sum::<Constraint<F>>()
                    - rvar(i, 1),
            );
        }
    }

    fn witness_partial_round(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(view.cell(0, T + i), c[i]);
            view.set(
                view.cell(1, i),
                std::iter::once(
                    m[i * T] * (view.get_at(view.cell(0, 0)) + c[0]).pow_small_vartime(a),
                )
                .chain((1..T).map(|j| m[i * T + j] * (view.get_at(view.cell(0, j)) + c[j])))
                .sum(),
            );
        }
    }
}

/// External ROM mode for the [`PermutationChip`].
///
/// See [`RcModeInternalRom`] for more information about internal and external ROM modes.
pub struct RcModeExternalRom<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> {
    /// Row offset of the IR chip (ROM lender), relative to wherever this ER chip (ROM borrower)
    /// itself lands when it's built/witnessed. Added to this ER chip's own row to get the IR chip's
    /// row.
    ir_chip_row_offset: isize,

    /// Column offset of the IR chip (ROM lender), relative to wherever this ER chip (ROM borrower)
    /// itself lands when it's built/witnessed. Added to this ER chip's own column to get the IR
    /// chip's column.
    ir_chip_column_offset: isize,

    _data: PhantomData<(F, C)>,
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Debug
    for RcModeExternalRom<F, C, T>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RcModeExternalRom")
            .field("ir_chip_row_offset", &self.ir_chip_row_offset)
            .field("ir_chip_column_offset", &self.ir_chip_column_offset)
            .field("_data", &self._data)
            .finish()
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Copy for RcModeExternalRom<F, C, T> {}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> Clone
    for RcModeExternalRom<F, C, T>
{
    fn clone(&self) -> Self {
        Self {
            ir_chip_row_offset: self.ir_chip_row_offset,
            ir_chip_column_offset: self.ir_chip_column_offset,
            _data: PhantomData,
        }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> RcModeExternalRom<F, C, T> {
    fn new(ir_chip_row_offset: isize, ir_chip_column_offset: isize) -> Self {
        Self {
            ir_chip_row_offset,
            ir_chip_column_offset,
            _data: PhantomData,
        }
    }

    /// Returns the ROM cell from the remote IR chip that mirrors this ER chip's own
    /// `view.cell(1, T + i)`, i.e. the i-th constant of whatever ARC `view` is currently positioned
    /// at.
    ///
    /// Since the IR and ER chips are built/witnessed through the exact same sequence of `.sub()` /
    /// `.sub_fn()` calls, the IR chip's corresponding cell always sits at the same `(1, T + i)`
    /// local offset from `view`'s current position, shifted only by the constant offset between the
    /// two chips' own roots.
    fn remote_rom_cell<G: Field256 + From<F>>(
        &self,
        view: &impl CircuitView<F, G>,
        i: usize,
    ) -> Cell
    where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        view.cell(
            self.ir_chip_row_offset,
            self.ir_chip_column_offset + (T + i) as isize,
        )
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize> internal::RcMode<F, T>
    for RcModeExternalRom<F, C, T>
{
    fn width(&self) -> usize {
        T * 2
    }

    fn build_full_round<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        _round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let a = F::ALPHA as isize;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.connect(
                self.remote_rom_cell(view, i).into(),
                view.cell(0, T + i).into(),
            );
            view.add_gate(
                0,
                (0..T)
                    .map(|j| make_const(m[i * T + j]) * ((var(j) + rvar(T + j, 0)) ^ a))
                    .sum::<Constraint<F>>()
                    - rvar(i, 1),
            );
        }
    }

    fn witness_full_round(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(view.cell(0, T + i), c[i]);
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| {
                        m[i * T + j] * (view.get_at(view.cell(0, j)) + c[j]).pow_small_vartime(a)
                    })
                    .sum(),
            );
        }
    }

    fn build_partial_round<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        _round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let a = F::ALPHA as isize;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.connect(
                self.remote_rom_cell(view, i).into(),
                view.cell(0, T + i).into(),
            );
            view.add_gate(
                0,
                std::iter::once(make_const(m[i * T]) * ((var(0) + var(T)) ^ a))
                    .chain((1..T).map(|j| make_const(m[i * T + j]) * (var(j) + var(T + j))))
                    .sum::<Constraint<F>>()
                    - rvar(i, 1),
            );
        }
    }

    fn witness_partial_round(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = &C::get_round_constants()[(round * T)..((round + 1) * T)];
        let a = F::ALPHA;
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(view.cell(0, T + i), c[i]);
            view.set(
                view.cell(1, i),
                std::iter::once(
                    m[i * T] * (view.get_at(view.cell(0, 0)) + c[0]).pow_small_vartime(a),
                )
                .chain((1..T).map(|j| m[i * T + j] * (view.get_at(view.cell(0, j)) + c[j])))
                .sum(),
            );
        }
    }
}

/// Poseidon permutation chip.
///
/// You may want to use [`PermutationChipHW`], [`PermutationChipIR`], or [`PermutationChipER`]
/// rather than referring to this struct directly.
pub struct PermutationChip<
    F: PrimeField,
    C: poseidon::Config<F, T>,
    M: internal::RcMode<F, T>,
    const T: usize,
> {
    rc: M,
    _data: PhantomData<(F, C)>,
}

impl<F: PrimeField, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize> Debug
    for PermutationChip<F, C, M, T>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PermutationChip")
            .field("rc", &self.rc)
            .finish()
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, M: internal::RcMode<F, T> + Default, const T: usize>
    Default for PermutationChip<F, C, M, T>
{
    fn default() -> Self {
        Self {
            rc: M::default(),
            _data: PhantomData,
        }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize> Copy
    for PermutationChip<F, C, M, T>
{
}

impl<F: PrimeField, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize> Clone
    for PermutationChip<F, C, M, T>
{
    fn clone(&self) -> Self {
        Self {
            rc: self.rc.clone(),
            _data: PhantomData,
        }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, const T: usize>
    PermutationChip<F, C, RcModeExternalRom<F, C, T>, T>
{
    /// Constructs a `PermutationChipER` that borrows its round constant ROM from an IR chip located
    /// at the given row/column offsets relative to wherever this ER chip itself is later built or
    /// witnessed (i.e. the IR chip's coordinates are this ER chip's own coordinates plus the given
    /// offsets).
    ///
    /// See [`RcModeInternalRom`] for the rationale.
    pub fn new(ir_chip_row_offset: isize, ir_chip_column_offset: isize) -> Self {
        Self {
            rc: RcModeExternalRom::new(ir_chip_row_offset, ir_chip_column_offset),
            _data: PhantomData,
        }
    }
}

impl<F: PrimeField, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize>
    PlonkChip<F, T, T> for PermutationChip<F, C, M, T>
{
    fn width(&self) -> usize {
        self.rc.width()
    }

    fn height(&self) -> usize {
        C::num_total_rounds() + 1
    }

    fn build<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; T],
    ) -> Result<[Option<Cell>; T]>
    where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let num_full_rounds = C::num_full_rounds();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        assert_eq!(num_total_rounds, num_full_rounds * 2 + num_partial_rounds);
        let mut view = view.sub(0, 0, Some(T), None);
        for i in 0..T {
            view.connect(inputs[i], view.cell(0, i).into());
        }
        for r in 0..num_full_rounds {
            view.sub_fn(r, 0, None, Some(2), |view| {
                self.rc.build_full_round(view, r)
            });
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub_fn(r, 0, None, Some(2), |view| {
                self.rc.build_partial_round(view, r)
            });
        }
        for r in (num_full_rounds + num_partial_rounds)..num_total_rounds {
            view.sub_fn(r, 0, None, Some(2), |view| {
                self.rc.build_full_round(view, r)
            });
        }
        Ok(std::array::from_fn(|i| {
            Some(view.cell(num_total_rounds, i))
        }))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; T],
    ) -> Result<[CellOrUnconstrained<F>; T]> {
        let num_full_rounds = C::num_full_rounds();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        assert_eq!(num_total_rounds, num_full_rounds * 2 + num_partial_rounds);
        let mut view = view.sub(0, 0, Some(T), None);
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i).into());
        }
        for r in 0..num_full_rounds {
            view.sub_fn(r, 0, None, Some(2), |view| {
                self.rc.witness_full_round(view, r)
            });
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub_fn(r, 0, None, Some(2), |view| {
                self.rc.witness_partial_round(view, r)
            });
        }
        for r in (num_full_rounds + num_partial_rounds)..num_total_rounds {
            view.sub_fn(r, 0, None, Some(2), |view| {
                self.rc.witness_full_round(view, r)
            });
        }
        Ok(std::array::from_fn(|i| {
            view.cell(num_total_rounds, i).into()
        }))
    }
}

/// Poseidon permutation chip with [hard-wired round constants](`RcModeHardWired`).
pub type PermutationChipHW<F, C, const T: usize> =
    PermutationChip<F, C, RcModeHardWired<F, C, T>, T>;

/// Poseidon permutation chip with [internal ROM storage for round constants](`RcModeInternalRom`).
pub type PermutationChipIR<F, C, const T: usize> =
    PermutationChip<F, C, RcModeInternalRom<F, C, T>, T>;

/// Poseidon permutation chip with [external ROM storage for round constants](`RcModeExternalRom`).
pub type PermutationChipER<F, C, const T: usize> =
    PermutationChip<F, C, RcModeExternalRom<F, C, T>, T>;

#[cfg(all(test, feature = "bluesky", feature = "goldilocks"))]
mod tests {
    use super::*;
    use primitive_types::H256;
    use starkom_bluesky::Scalar as BS;
    use starkom_goldilocks::{GL, GL4};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon as poseidon1;
    use std::str::FromStr;

    fn parse<T: FromStr<Err: Debug>>(s: &'static str) -> T {
        s.parse().unwrap()
    }

    fn test_permutation_bluesky_impl<const T: usize>(
        chip: &impl PlonkChip<BS, T, T>,
        inputs: [BS; T],
        expected_output: [BS; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        assert_eq!(chip.height(), 65);
        let mut builder = CircuitBuilder::<BS, BS>::default();
        let output = builder.sub_chip(0, 0, chip, std::array::from_fn(|_| None))?;
        builder.declare_public_rows([output[0].unwrap().row()]);
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 65);
        assert_eq!(circuit.degree_bound(), 128);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 65);
        assert_eq!(witness.degree_bound(), 128);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<BS>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 128);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 128 << blowup_log2);
        let circuit = circuit.to_compressed::<Sha2Hash<BS>>(options);
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

    fn test_permutation_goldilocks_impl<const T: usize>(
        chip: &impl PlonkChip<GL, T, T>,
        inputs: [GL; T],
        expected_output: [GL; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        assert_eq!(chip.height(), 31);
        let mut builder = CircuitBuilder::<GL, GL4>::default();
        let output = builder.sub_chip(0, 0, chip, std::array::from_fn(|_| None))?;
        builder.declare_public_rows([output[0].unwrap().row()]);
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 31);
        assert_eq!(circuit.degree_bound(), 64);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 31);
        assert_eq!(witness.degree_bound(), 64);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<GL4>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 64);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 64 << blowup_log2);
        let circuit = circuit.to_compressed::<Sha2Hash<GL4>>(options);
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

    fn test_perm_bluesky_hw<
        Cfg: poseidon1::Config<BS, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [BS; T],
        expected_output: [BS; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipHW::<BS, Cfg, T>::default();
        assert_eq!(chip.width(), T);
        test_permutation_bluesky_impl::<T>(
            &chip,
            inputs,
            expected_output,
            blowup_log2,
            circuit_commitment,
        )
    }

    #[test]
    fn test_permutation_t3_hw() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into()];
        let outputs = [
            parse("0x7b68dcd80fa751ee8f2d76043bfd92c685601c79189393fc76e03c5214eed32b"),
            parse("0x0fbcb5720b463bf7e2ccabf373e77d2c10d27e6549f34cfa33eb2d06ea8b900a"),
            parse("0x26e03abfcc62da0101516b07aede8bc676a10c47299a57bedc6d9fe80484f3da"),
        ];
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse("0xea3b13b08ffab35ed39016acfacab9fe0b86978482aaa9b8032da18aadf3f5fd")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0x962586b041835c5dc7ec795301b13aebc9ed8ac3a03637890b079c3a7762a129")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0x74005e8550d520b3fb2eb9444d4e19fe604031bec42fc7c61fda7c564208f06c")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_hw() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x12dde8a4c46760e349670d241e36ca7abacc991233039f8deaf6c58ce2230ef6"),
            parse("0x61e95d9456e9223b4d7926dabae10009da2b6fb9134ade8405f6ef1424e66aa1"),
            parse("0x2fcce25ab9efb3e26276f3b3aff1e02cdf82df48ce8d3eadbff900cfe015775b"),
            parse("0x2580707d57a8c1c0cad368e8d5705ffd96f269d66e1cd6f1433f93a3c66d9bf8"),
        ];
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse("0x81c58afcd2fdd68feb646fd0a7a35d998491875f6e9c8a7b8cb22c91f63a5f79")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0x31ea4abe933234e0399a827743c008653a1708345485f9c3c1c5a7125147de76")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0x2ab529545598bee59f8dd5e4d02efb27812617a9321fdf37de928255ca9845dd")
            )
            .is_ok()
        );
    }

    fn test_perm_goldilocks_hw<Cfg: poseidon1::Config<GL, T>, const T: usize>(
        inputs: [GL; T],
        expected_output: [GL; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipHW::<GL, Cfg, T>::default();
        assert_eq!(chip.width(), T);
        test_permutation_goldilocks_impl::<T>(
            &chip,
            inputs,
            expected_output,
            blowup_log2,
            circuit_commitment,
        )
    }

    #[test]
    fn test_permutation_t12_hw() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            parse("0x056bda38ad308e78"),
            parse("0x1f38944238b8ccd0"),
            parse("0x80bef63a171f3156"),
            parse("0x27bbc645b2a3198c"),
            parse("0x9befae3f221509b3"),
            parse("0xa1cfa54ae2c44c9e"),
            parse("0xa1c876869f1c52f8"),
            parse("0x7ffa21471eff65af"),
            parse("0xdc565450ad52b99e"),
            parse("0x4b8b1daf8e8ea3c6"),
            parse("0xf866b42495e61984"),
            parse("0x7af57b5f91f196fe"),
        ];
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                1,
                parse("0xca22a19976890ccd04f8a4e4bddb12fef5a300fbba62856718a70ea978cb661d")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0xb6ba445d2ace00ab4359132cb6014da23a2caf3f09d1095a74005da51f364b48")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0x22cd7b09e43a265511c5337029510a511c6726f2ad632799deb0c56e4d47aefb")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t16_hw() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            parse("0x6a84bf02be1f328d"),
            parse("0xec14d274b936a21a"),
            parse("0xc0539d7bd4eb66de"),
            parse("0xb317ecf41fa8d55b"),
            parse("0x80b0d36f66671f8a"),
            parse("0x74a1592b9a16e832"),
            parse("0x65e53afadfadc8c3"),
            parse("0xa0007e5ee96ee4b2"),
            parse("0x6dd5661a877003a8"),
            parse("0xc36a09c2dc25cd6e"),
            parse("0xcbda3d58f7cf85f4"),
            parse("0x34cb1d63c35596cf"),
            parse("0x4fcd09b24769e281"),
            parse("0x6c514f906998c65d"),
            parse("0xc447035d8d71952b"),
            parse("0x591863454267826f"),
        ];
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                1,
                parse("0x0640901c4ac5acf1af213b200bd97c9969d4d29d09ba88a1b3186b69cd0bcf37")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x025f93dcd636a111f13b7ad1d00b95e767fd846916d87cd7491a5cbdf4d99c3d")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0x2023be32c15420fc5575734f9b46e6aec619e0cddb9243f03852f415f84f9004")
            )
            .is_ok()
        );
    }

    fn test_perm_bluesky_ir<
        Cfg: poseidon1::Config<BS, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [BS; T],
        expected_output: [BS; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipIR::<BS, Cfg, T>::default();
        assert_eq!(chip.width(), T * 2);
        test_permutation_bluesky_impl::<T>(
            &chip,
            inputs,
            expected_output,
            blowup_log2,
            circuit_commitment,
        )
    }

    #[test]
    fn test_permutation_t3_ir() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into()];
        let outputs = [
            parse("0x7b68dcd80fa751ee8f2d76043bfd92c685601c79189393fc76e03c5214eed32b"),
            parse("0x0fbcb5720b463bf7e2ccabf373e77d2c10d27e6549f34cfa33eb2d06ea8b900a"),
            parse("0x26e03abfcc62da0101516b07aede8bc676a10c47299a57bedc6d9fe80484f3da"),
        ];
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse("0xe9c445a6d8f5dd33274620bf39b91ebfbda1774040725d1e856fe7b3ba36db0f")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0xb1a33eb956c0eb259298bf3e850cb1d2f0f88a92d52d94fcd8c1b7b82010f6fb")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0xda6d4e1f7af1b92e8c416ce34929f3c384ea5f5ad129d31d55b747e09e25f559")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_ir() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x12dde8a4c46760e349670d241e36ca7abacc991233039f8deaf6c58ce2230ef6"),
            parse("0x61e95d9456e9223b4d7926dabae10009da2b6fb9134ade8405f6ef1424e66aa1"),
            parse("0x2fcce25ab9efb3e26276f3b3aff1e02cdf82df48ce8d3eadbff900cfe015775b"),
            parse("0x2580707d57a8c1c0cad368e8d5705ffd96f269d66e1cd6f1433f93a3c66d9bf8"),
        ];
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse("0xb05298f897fe2bda796bad79de1306b97c4d25b2a3fd2fac9c7e235a90ec4c3d")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0x8635b5e274847291a6d29f508b83d61ba96dedc174316c66c1b48745c02ff020")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0x6c549d26f46170ded362e705c62ee27e1cee45a117d67bb7a343396182e863ba")
            )
            .is_ok()
        );
    }

    fn test_perm_goldilocks_ir<Cfg: poseidon1::Config<GL, T>, const T: usize>(
        inputs: [GL; T],
        expected_output: [GL; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipIR::<GL, Cfg, T>::default();
        assert_eq!(chip.width(), T * 2);
        test_permutation_goldilocks_impl::<T>(
            &chip,
            inputs,
            expected_output,
            blowup_log2,
            circuit_commitment,
        )
    }

    #[test]
    fn test_permutation_t12_ir() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            parse("0x056bda38ad308e78"),
            parse("0x1f38944238b8ccd0"),
            parse("0x80bef63a171f3156"),
            parse("0x27bbc645b2a3198c"),
            parse("0x9befae3f221509b3"),
            parse("0xa1cfa54ae2c44c9e"),
            parse("0xa1c876869f1c52f8"),
            parse("0x7ffa21471eff65af"),
            parse("0xdc565450ad52b99e"),
            parse("0x4b8b1daf8e8ea3c6"),
            parse("0xf866b42495e61984"),
            parse("0x7af57b5f91f196fe"),
        ];
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                1,
                parse("0x1afa8b24d35602dfc572885c74efe4b078b98d89c54ec62b52fde8fae1f70caf")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0xaee97d8558f3198e79a0e5ee372b47a6a2876070293e916d4aa1129fa2204a54")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0xed8f71073d457d8e5a116bfb1b015a223bd630a5443ef5d686bc6c5e1206b101")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t16_ir() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            parse("0x6a84bf02be1f328d"),
            parse("0xec14d274b936a21a"),
            parse("0xc0539d7bd4eb66de"),
            parse("0xb317ecf41fa8d55b"),
            parse("0x80b0d36f66671f8a"),
            parse("0x74a1592b9a16e832"),
            parse("0x65e53afadfadc8c3"),
            parse("0xa0007e5ee96ee4b2"),
            parse("0x6dd5661a877003a8"),
            parse("0xc36a09c2dc25cd6e"),
            parse("0xcbda3d58f7cf85f4"),
            parse("0x34cb1d63c35596cf"),
            parse("0x4fcd09b24769e281"),
            parse("0x6c514f906998c65d"),
            parse("0xc447035d8d71952b"),
            parse("0x591863454267826f"),
        ];
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                1,
                parse("0x279cc1337eda951ccb13af87dc37f06b20194125ed7be83fa9e73ede1a9e7580")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x46f1d37cec11e7cce214db842eccff2d70dd9d82397da85cb5cd5dc6bbe874b8")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0x8af4f64c7bb731cf43c66453d2627852fe5e60152d66e40b1a49398553e0608f")
            )
            .is_ok()
        );
    }

    fn test_perm_bluesky_er<
        Cfg: poseidon1::Config<BS, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [BS; T],
        expected_output: [BS; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip_ir = PermutationChipIR::<BS, Cfg, T>::default();
        assert_eq!(chip_ir.width(), T * 2);
        assert_eq!(chip_ir.height(), 65);
        let ir_width = chip_ir.width();

        let chip_er = PermutationChipER::<BS, Cfg, T>::new(0, -(ir_width as isize));
        assert_eq!(chip_er.width(), T * 2);
        assert_eq!(chip_er.height(), 65);
        let er_width = chip_er.width();

        let mut builder = CircuitBuilder::<BS, BS>::default();
        let ir_output = builder.sub_chip(0, 0, &chip_ir, std::array::from_fn(|_| None))?;
        let er_output = builder.sub_chip(0, ir_width, &chip_er, std::array::from_fn(|_| None))?;

        for i in 0..T {
            builder.connect(ir_output[i], er_output[i]);
        }
        builder.declare_public_rows([ir_output[0].unwrap().row()]);

        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 65);
        assert_eq!(circuit.num_columns(), ir_width + er_width);

        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 65);
        assert_eq!(witness.num_columns(), ir_width + er_width);
        let ir_output = witness.sub_chip(0, 0, &chip_ir, inputs.map(|input| input.into()))?;
        let er_output =
            witness.sub_chip(0, ir_width, &chip_er, inputs.map(|input| input.into()))?;

        circuit.check_witness(&witness).unwrap();

        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<BS>>(witness, options.clone())?;

        let circuit = circuit.to_compressed::<Sha2Hash<BS>>(options);
        assert_eq!(circuit.commitment(), circuit_commitment);

        let public_inputs = circuit.verify(&proof)?;
        let get_value = |output: CellOrUnconstrained<BS>| match output {
            CellOrUnconstrained::Cell(cell) => public_inputs[&cell],
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        assert!(ir_output.into_iter().zip(er_output).enumerate().all(
            |(i, (ir_output, er_output))| get_value(ir_output) == expected_output[i]
                && get_value(er_output) == expected_output[i]
        ));
        Ok(())
    }

    #[test]
    fn test_permutation_t3_er() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into()];
        let outputs = [
            parse("0x7b68dcd80fa751ee8f2d76043bfd92c685601c79189393fc76e03c5214eed32b"),
            parse("0x0fbcb5720b463bf7e2ccabf373e77d2c10d27e6549f34cfa33eb2d06ea8b900a"),
            parse("0x26e03abfcc62da0101516b07aede8bc676a10c47299a57bedc6d9fe80484f3da"),
        ];
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse("0xa7d9e516c43ea92aa017c45fc367e84fbdc059344550a847777ad045699621ba")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0x6ae0d02555b056c67a58b64d0cfe01f13ed07902851b98795156c4fbf9d3595f")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0x63bac3890c17dff52faad356f47db5ee29f0a4f7e58d0ef78770cc384b396783")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_er() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x12dde8a4c46760e349670d241e36ca7abacc991233039f8deaf6c58ce2230ef6"),
            parse("0x61e95d9456e9223b4d7926dabae10009da2b6fb9134ade8405f6ef1424e66aa1"),
            parse("0x2fcce25ab9efb3e26276f3b3aff1e02cdf82df48ce8d3eadbff900cfe015775b"),
            parse("0x2580707d57a8c1c0cad368e8d5705ffd96f269d66e1cd6f1433f93a3c66d9bf8"),
        ];
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse("0x68629ac84ef302a47226cad2a29c563b8b6e0ce23495f7636938985e0a7b7ca7")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0x0286cfd6feb59ee7b00a099973c97dea3d96cbcdd02164bbb83194a9dc6a81fe")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0x069717fd1bf79551cada48ad780de13d486a6463bd474d39e355ce5747f6054a")
            )
            .is_ok()
        );
    }

    fn test_perm_goldilocks_er<Cfg: poseidon1::Config<GL, T>, const T: usize>(
        inputs: [GL; T],
        expected_output: [GL; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip_ir = PermutationChipIR::<GL, Cfg, T>::default();
        assert_eq!(chip_ir.width(), T * 2);
        assert_eq!(chip_ir.height(), 31);
        let ir_width = chip_ir.width();
        let ir_height = chip_ir.height();

        let chip_er = PermutationChipER::<GL, Cfg, T>::new(-(ir_height as isize), 0);
        assert_eq!(chip_er.width(), T * 2);
        assert_eq!(chip_er.height(), 31);

        let mut builder = CircuitBuilder::<GL, GL4>::default();
        let ir_output = builder.sub_chip(0, 0, &chip_ir, std::array::from_fn(|_| None))?;
        let er_output = builder.sub_chip(ir_height, 0, &chip_er, std::array::from_fn(|_| None))?;

        for i in 0..T {
            builder.connect(ir_output[i], er_output[i]);
        }
        builder.declare_public_rows([ir_output[0].unwrap().row(), er_output[0].unwrap().row()]);

        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 62);
        assert_eq!(circuit.num_columns(), ir_width);

        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 62);
        assert_eq!(witness.num_columns(), ir_width);
        let ir_output = witness.sub_chip(0, 0, &chip_ir, inputs.map(|input| input.into()))?;
        let er_output =
            witness.sub_chip(ir_height, 0, &chip_er, inputs.map(|input| input.into()))?;

        circuit.check_witness(&witness).unwrap();

        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<GL4>>(witness, options.clone())?;

        let circuit = circuit.to_compressed::<Sha2Hash<GL4>>(options);
        assert_eq!(circuit.commitment(), circuit_commitment);

        let public_inputs = circuit.verify(&proof)?;
        let get_value = |output: CellOrUnconstrained<GL>| match output {
            CellOrUnconstrained::Cell(cell) => public_inputs[&cell],
            CellOrUnconstrained::Unconstrained(value) => value,
        };
        assert!(ir_output.into_iter().zip(er_output).enumerate().all(
            |(i, (ir_output, er_output))| get_value(ir_output) == expected_output[i]
                && get_value(er_output) == expected_output[i]
        ));
        Ok(())
    }

    #[test]
    fn test_permutation_t12_er() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            parse("0x056bda38ad308e78"),
            parse("0x1f38944238b8ccd0"),
            parse("0x80bef63a171f3156"),
            parse("0x27bbc645b2a3198c"),
            parse("0x9befae3f221509b3"),
            parse("0xa1cfa54ae2c44c9e"),
            parse("0xa1c876869f1c52f8"),
            parse("0x7ffa21471eff65af"),
            parse("0xdc565450ad52b99e"),
            parse("0x4b8b1daf8e8ea3c6"),
            parse("0xf866b42495e61984"),
            parse("0x7af57b5f91f196fe"),
        ];
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                1,
                parse("0x7189de00db2ab4b59f2dd86e335a2967170cfc7fd2254bd1e106e38d03144c07")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0x10d09f7deea596a8eef52d5b4b6abf84f8cbdfde9ee09d994ea337745374d07d")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0x5d432700c3305b2d13a1f647fec239fdbd5c6f4b2e4c5a2c5b0ddb55eb2ea71c")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t16_er() {
        let inputs = std::array::from_fn(|i| i.try_into().unwrap());
        let outputs = [
            parse("0x6a84bf02be1f328d"),
            parse("0xec14d274b936a21a"),
            parse("0xc0539d7bd4eb66de"),
            parse("0xb317ecf41fa8d55b"),
            parse("0x80b0d36f66671f8a"),
            parse("0x74a1592b9a16e832"),
            parse("0x65e53afadfadc8c3"),
            parse("0xa0007e5ee96ee4b2"),
            parse("0x6dd5661a877003a8"),
            parse("0xc36a09c2dc25cd6e"),
            parse("0xcbda3d58f7cf85f4"),
            parse("0x34cb1d63c35596cf"),
            parse("0x4fcd09b24769e281"),
            parse("0x6c514f906998c65d"),
            parse("0xc447035d8d71952b"),
            parse("0x591863454267826f"),
        ];
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                1,
                parse("0x7e677f8d3fae6797b43acdd5597046655c0cf43c66b00a4bd228349b168a6351")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x7027d1db5706dd2949b15ac940f64f73803d27e0825b321046251b33c0363ab0")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0x9932e68664013aaed77de33ab431b7f7d0eed3e384cfe7e2c9efa9098e6220ce")
            )
            .is_ok()
        );
    }
}
