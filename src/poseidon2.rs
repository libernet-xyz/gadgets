use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, rvar, var,
};
use starkom_poseidon2 as poseidon2;
use std::fmt::{Debug, Formatter};
use std::marker::PhantomData;

mod internal {
    use super::*;

    /// Encodes the operations that differ between round constant modes
    /// ([hard-wired](`RcModeHardWired`), [internal ROM](`RcModeInternalRom`),
    /// [external ROM](`RcModeExternalRom`)).
    pub trait RcMode<C: poseidon2::Config<Scalar, T>, const T: usize>:
        Debug + Copy + Clone
    {
        fn width(&self) -> usize;

        fn build_first_fl_and_arc(&self, view: &mut impl CircuitView, inputs: [Option<Cell>; T]);

        fn witness_first_fl_and_arc(
            &self,
            view: &mut impl WitnessView,
            inputs: [CellOrUnconstrained; T],
        );

        fn build_linear_and_next_arc(
            &self,
            view: &mut impl CircuitView,
            round: usize,
            matrix: &[Scalar],
        );

        fn witness_linear_and_next_arc(
            &self,
            view: &mut impl WitnessView,
            round: usize,
            matrix: &[Scalar],
        );

        fn build_fl_and_next_arc(&self, view: &mut impl CircuitView, round: usize) {
            self.build_linear_and_next_arc(view, round, C::get_external_matrix());
        }

        fn witness_fl_and_next_arc(&self, view: &mut impl WitnessView, round: usize) {
            self.witness_linear_and_next_arc(view, round, C::get_external_matrix());
        }

        fn build_pl_and_next_arc(&self, view: &mut impl CircuitView, round: usize) {
            self.build_linear_and_next_arc(view, round, C::get_internal_matrix());
        }

        fn witness_pl_and_next_arc(&self, view: &mut impl WitnessView, round: usize) {
            self.witness_linear_and_next_arc(view, round, C::get_internal_matrix());
        }
    }
}

/// Hard-wired RC mode for the Poseidon2 permutation.
///
/// In this mode the chip takes exactly `T` columns, with `T` being the state vector size, but uses
/// a different gate for every ARC layer. You should use this mode only if the elevated number of
/// gates is not a concern.
///
/// This mode is best suited for circuits that run a small number of hashes, such as preimage
/// knowledge proofs and zkMAC signatures.
pub struct RcModeHardWired<C: poseidon2::Config<Scalar, T>, const T: usize> {
    _data: PhantomData<C>,
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Debug for RcModeHardWired<C, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RcModeHardWired").finish()
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Default for RcModeHardWired<C, T> {
    fn default() -> Self {
        Self {
            _data: Default::default(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Copy for RcModeHardWired<C, T> {}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Clone for RcModeHardWired<C, T> {
    fn clone(&self) -> Self {
        Self {
            _data: self._data.clone(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> internal::RcMode<C, T>
    for RcModeHardWired<C, T>
{
    fn width(&self) -> usize {
        T
    }

    fn build_first_fl_and_arc(&self, view: &mut impl CircuitView, inputs: [Option<Cell>; T]) {
        for i in 0..T {
            view.connect(inputs[i], Some(view.cell(0, i)));
        }
        let m = C::get_external_matrix();
        let c = C::get_round_constants();
        for i in 0..T {
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * m[i * T + j])
                    .sum::<Constraint>()
                    + c[i]
                    - rvar(i, 1),
            );
        }
    }

    fn witness_first_fl_and_arc(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; T],
    ) {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let m = C::get_external_matrix();
        let c = C::get_round_constants();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>()
                    + c[i],
            );
        }
    }

    fn build_linear_and_next_arc(&self, view: &mut impl CircuitView, round: usize, m: &[Scalar]) {
        let c = C::get_round_constants();
        for i in 0..T {
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * m[i * T + j])
                    .sum::<Constraint>()
                    + c[(round + 1) * T + i]
                    - rvar(i, 1),
            );
        }
    }

    fn witness_linear_and_next_arc(&self, view: &mut impl WitnessView, round: usize, m: &[Scalar]) {
        let c = C::get_round_constants();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>()
                    + c[(round + 1) * T + i],
            );
        }
    }
}

/// Internal ROM mode for the [`PermutationChip`].
///
/// In this mode the Poseidon2 round constants are stored in separate columns so that we can reuse
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
pub struct RcModeInternalRom<C: poseidon2::Config<Scalar, T>, const T: usize> {
    _data: PhantomData<C>,
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Debug for RcModeInternalRom<C, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RcModeInternalRom").finish()
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Default for RcModeInternalRom<C, T> {
    fn default() -> Self {
        Self {
            _data: Default::default(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Copy for RcModeInternalRom<C, T> {}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Clone for RcModeInternalRom<C, T> {
    fn clone(&self) -> Self {
        Self {
            _data: self._data.clone(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> internal::RcMode<C, T>
    for RcModeInternalRom<C, T>
{
    fn width(&self) -> usize {
        T * 2
    }

    fn build_first_fl_and_arc(&self, view: &mut impl CircuitView, inputs: [Option<Cell>; T]) {
        for i in 0..T {
            view.connect(inputs[i], Some(view.cell(0, i)));
        }
        let m = C::get_external_matrix();
        let c = C::get_round_constants();
        for i in 0..T {
            view.add_gate(1, var(T + i) - c[i]);
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * m[i * T + j])
                    .sum::<Constraint>()
                    + rvar(T + i, 1)
                    - rvar(i, 1),
            );
        }
    }

    fn witness_first_fl_and_arc(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; T],
    ) {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let m = C::get_external_matrix();
        let c = C::get_round_constants();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>()
                    + c[i],
            );
            view.set(view.cell(1, T + i), c[i]);
        }
    }

    fn build_linear_and_next_arc(&self, view: &mut impl CircuitView, round: usize, m: &[Scalar]) {
        let c = C::get_round_constants();
        for i in 0..T {
            view.add_gate(1, var(T + i) - c[(round + 1) * T + i]);
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * m[i * T + j])
                    .sum::<Constraint>()
                    + rvar(T + i, 1)
                    - rvar(i, 1),
            );
        }
    }

    fn witness_linear_and_next_arc(&self, view: &mut impl WitnessView, round: usize, m: &[Scalar]) {
        let c = C::get_round_constants();
        for i in 0..T {
            view.set(view.cell(1, T + i), c[(round + 1) * T + i]);
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>()
                    + c[(round + 1) * T + i],
            );
        }
    }
}

/// External ROM mode for the [`PermutationChip`].
///
/// See [`RcModeInternalRom`] for more information about internal and external ROM modes.
pub struct RcModeExternalRom<C: poseidon2::Config<Scalar, T>, const T: usize> {
    /// Row offset of the IR chip (ROM lender), relative to wherever this ER chip (ROM borrower)
    /// itself lands when it's built/witnessed. Added to this ER chip's own row to get the IR chip's
    /// row.
    ir_chip_row_offset: isize,

    /// Column offset of the IR chip (ROM lender), relative to wherever this ER chip (ROM borrower)
    /// itself lands when it's built/witnessed. Added to this ER chip's own column to get the IR
    /// chip's column.
    ir_chip_column_offset: isize,

    _data: PhantomData<C>,
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Debug for RcModeExternalRom<C, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RcModeExternalRom")
            .field("ir_chip_row_offset", &self.ir_chip_row_offset)
            .field("ir_chip_column_offset", &self.ir_chip_column_offset)
            .field("_data", &self._data)
            .finish()
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Copy for RcModeExternalRom<C, T> {}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> Clone for RcModeExternalRom<C, T> {
    fn clone(&self) -> Self {
        Self {
            ir_chip_row_offset: self.ir_chip_row_offset,
            ir_chip_column_offset: self.ir_chip_column_offset,
            _data: self._data.clone(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> RcModeExternalRom<C, T> {
    fn new(ir_chip_row_offset: isize, ir_chip_column_offset: isize) -> Self {
        Self {
            ir_chip_row_offset,
            ir_chip_column_offset,
            _data: PhantomData::default(),
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
    fn remote_rom_cell(&self, view: &impl CircuitView, i: usize) -> Cell {
        view.cell(
            self.ir_chip_row_offset + 1,
            self.ir_chip_column_offset + (T + i) as isize,
        )
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize> internal::RcMode<C, T>
    for RcModeExternalRom<C, T>
{
    fn width(&self) -> usize {
        T * 2
    }

    fn build_first_fl_and_arc(&self, view: &mut impl CircuitView, inputs: [Option<Cell>; T]) {
        for i in 0..T {
            view.connect(inputs[i], Some(view.cell(0, i)));
        }
        let m = C::get_external_matrix();
        for i in 0..T {
            view.connect(
                self.remote_rom_cell(view, i).into(),
                view.cell(1, T + i).into(),
            );
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * m[i * T + j])
                    .sum::<Constraint>()
                    + rvar(T + i, 1)
                    - rvar(i, 1),
            );
        }
    }

    fn witness_first_fl_and_arc(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; T],
    ) {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let m = C::get_external_matrix();
        let c = C::get_round_constants();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>()
                    + c[i],
            );
            view.set(view.cell(1, T + i), c[i]);
        }
    }

    fn build_linear_and_next_arc(&self, view: &mut impl CircuitView, _round: usize, m: &[Scalar]) {
        for i in 0..T {
            view.connect(
                self.remote_rom_cell(view, i).into(),
                view.cell(1, T + i).into(),
            );
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * m[i * T + j])
                    .sum::<Constraint>()
                    + rvar(T + i, 1)
                    - rvar(i, 1),
            );
        }
    }

    fn witness_linear_and_next_arc(&self, view: &mut impl WitnessView, round: usize, m: &[Scalar]) {
        let c = C::get_round_constants();
        for i in 0..T {
            view.set(view.cell(1, T + i), c[(round + 1) * T + i]);
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>()
                    + c[(round + 1) * T + i],
            );
        }
    }
}

/// Poseidon2 permutation chip.
///
/// You may want to use [`PermutationChipHW`], [`PermutationChipIR`], or [`PermutationChipER`]
/// rather than referring to this struct directly.
pub struct PermutationChip<
    C: poseidon2::Config<Scalar, T>,
    M: internal::RcMode<C, T>,
    const T: usize,
> {
    rc: M,
    _data: PhantomData<C>,
}

impl<C: poseidon2::Config<Scalar, T>, M: internal::RcMode<C, T>, const T: usize> Debug
    for PermutationChip<C, M, T>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PermutationChip")
            .field("rc", &self.rc)
            .finish()
    }
}

impl<C: poseidon2::Config<Scalar, T>, M: internal::RcMode<C, T> + Default, const T: usize> Default
    for PermutationChip<C, M, T>
{
    fn default() -> Self {
        Self {
            rc: M::default(),
            _data: PhantomData::default(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, M: internal::RcMode<C, T>, const T: usize> Copy
    for PermutationChip<C, M, T>
{
}

impl<C: poseidon2::Config<Scalar, T>, M: internal::RcMode<C, T>, const T: usize> Clone
    for PermutationChip<C, M, T>
{
    fn clone(&self) -> Self {
        Self {
            rc: self.rc.clone(),
            _data: self._data.clone(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, M: internal::RcMode<C, T>, const T: usize>
    PermutationChip<C, M, T>
{
    pub const FIRST_ARC_HEIGHT: usize = 2;
    pub const ROUND_HEIGHT: usize = 3;

    fn build_full_sbox(&self, view: &mut impl CircuitView) {
        for i in 0..T {
            view.add_gate(0, (rvar(i, -1) ^ 3) - rvar(i, 0));
            view.add_gate(0, (rvar(i, -1) ^ 2) * rvar(i, 0) - rvar(i, 1));
        }
    }

    fn witness_full_sbox(&self, view: &mut impl WitnessView) {
        for i in 0..T {
            let state = view.get_at(view.cell(-1, i));
            view.set(view.cell(0, i), state.cube());
            view.set(view.cell(1, i), state.square().square() * state);
        }
    }

    fn build_partial_sbox(&self, view: &mut impl CircuitView) {
        view.add_gate(0, (rvar(0, -1) ^ 3) - rvar(0, 0));
        view.add_gate(0, (rvar(0, -1) ^ 2) * rvar(0, 0) - rvar(0, 1));
        for i in 1..T {
            view.connect(Some(view.cell(-1, i)), Some(view.cell(1, i)));
        }
    }

    fn witness_partial_sbox(&self, view: &mut impl WitnessView) {
        let state = view.get_at(view.cell(-1, 0));
        view.set(view.cell(0, 0), state.cube());
        view.set(view.cell(1, 0), state.square().square() * state);
        for i in 1..T {
            view.copy(view.cell(-1, i).into(), view.cell(1, i));
        }
    }

    fn build_last_linear(&self, view: &mut impl CircuitView) {
        let m = C::get_external_matrix();
        for i in 0..T {
            view.add_gate(
                0,
                ((0..T)
                    .map(|j| rvar(j, 0) * m[i * T + j])
                    .sum::<Constraint>())
                    - rvar(i, 1),
            );
        }
    }

    fn witness_last_linear(&self, view: &mut impl WitnessView) {
        let m = C::get_external_matrix();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>(),
            );
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, const T: usize>
    PermutationChip<C, RcModeExternalRom<C, T>, T>
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
            _data: PhantomData::default(),
        }
    }
}

impl<C: poseidon2::Config<Scalar, T>, M: internal::RcMode<C, T>, const T: usize> PlonkChip<T, T>
    for PermutationChip<C, M, T>
{
    fn width(&self) -> usize {
        self.rc.width()
    }

    fn height(&self) -> usize {
        Self::FIRST_ARC_HEIGHT + Self::ROUND_HEIGHT * C::num_total_rounds()
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; T],
    ) -> Result<[Option<Cell>; T]> {
        let num_full_rounds = C::num_full_rounds();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        assert_eq!(num_total_rounds, num_full_rounds * 2 + num_partial_rounds);
        self.rc.build_first_fl_and_arc(view, inputs);
        let mut view = view.sub(
            Self::FIRST_ARC_HEIGHT,
            0,
            Some(self.width()),
            Some(Self::ROUND_HEIGHT * num_total_rounds),
        );
        for r in 0..num_full_rounds {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.build_full_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.build_fl_and_next_arc(view, r)
                });
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.build_partial_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.build_pl_and_next_arc(view, r)
                });
        }
        for r in (num_full_rounds + num_partial_rounds)..(num_total_rounds - 1) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.build_full_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.build_fl_and_next_arc(view, r)
                });
        }
        view.sub(
            (num_total_rounds - 1) * Self::ROUND_HEIGHT,
            0,
            None,
            Some(Self::ROUND_HEIGHT),
        )
        .sub_fn(0, 0, None, None, |view| self.build_full_sbox(view))
        .sub_fn(1, 0, None, None, |view| self.build_last_linear(view));
        Ok(std::array::from_fn(|i| {
            Some(view.cell(num_total_rounds * Self::ROUND_HEIGHT - 1, i))
        }))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; T],
    ) -> Result<[CellOrUnconstrained; T]> {
        let num_full_rounds = C::num_full_rounds();
        let num_partial_rounds = C::num_partial_rounds();
        let num_total_rounds = C::num_total_rounds();
        assert_eq!(num_total_rounds, num_full_rounds * 2 + num_partial_rounds);
        self.rc.witness_first_fl_and_arc(view, inputs);
        let mut view = view.sub(
            2,
            0,
            Some(self.width()),
            Some(Self::ROUND_HEIGHT * num_total_rounds),
        );
        for r in 0..num_full_rounds {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.witness_full_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.witness_fl_and_next_arc(view, r)
                });
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.witness_partial_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.witness_pl_and_next_arc(view, r)
                });
        }
        for r in (num_full_rounds + num_partial_rounds)..(num_total_rounds - 1) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.witness_full_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.witness_fl_and_next_arc(view, r)
                });
        }
        view.sub(
            (num_total_rounds - 1) * Self::ROUND_HEIGHT,
            0,
            None,
            Some(Self::ROUND_HEIGHT),
        )
        .sub_fn(0, 0, None, None, |view| self.witness_full_sbox(view))
        .sub_fn(1, 0, None, None, |view| self.witness_last_linear(view));
        Ok(std::array::from_fn(|i| {
            view.cell(num_total_rounds * Self::ROUND_HEIGHT - 1, i)
                .into()
        }))
    }
}

/// Poseidon2 permutation chip with [hard-wired round constants](`RcModeHardWired`).
pub type PermutationChipHW<C, const T: usize> = PermutationChip<C, RcModeHardWired<C, T>, T>;

/// Poseidon2 permutation chip with [internal ROM storage for round constants](`RcModeInternalRom`).
pub type PermutationChipIR<C, const T: usize> = PermutationChip<C, RcModeInternalRom<C, T>, T>;

/// Poseidon2 permutation chip with [external ROM storage for round constants](`RcModeExternalRom`).
pub type PermutationChipER<C, const T: usize> = PermutationChip<C, RcModeExternalRom<C, T>, T>;

#[cfg(test)]
mod tests {
    use super::*;
    use primitive_types::H256;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};

    fn parse_hash(s: &'static str) -> H256 {
        s.parse().unwrap()
    }

    fn test_permutation_impl<const T: usize>(
        chip: &impl PlonkChip<T, T>,
        inputs: [Scalar; T],
        expected_output: [Scalar; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        assert_eq!(chip.height(), 194);
        let mut builder = CircuitBuilder::default();
        let output = builder.sub_chip(0, 0, chip, std::array::from_fn(|_| None))?;
        builder.declare_public_rows([output[0].unwrap().row()]);
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 194);
        assert_eq!(circuit.degree_bound(), 256);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 194);
        assert_eq!(witness.degree_bound(), 256);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 256);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 256 << blowup_log2);
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

    fn test_perm_hw<
        Cfg: poseidon2::Config<Scalar, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [Scalar; T],
        expected_output: [Scalar; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipHW::<Cfg, T>::default();
        assert_eq!(chip.width(), T);
        test_permutation_impl::<T>(
            &chip,
            inputs,
            expected_output,
            blowup_log2,
            circuit_commitment,
        )
    }

    #[test]
    fn test_permutation_t3_hw() {
        let inputs = [from_const(0), from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x6f30582cde48a25b26015b7f718ba2fb359e93029caf04d8d0b3e66b1d46b941"),
            parse_scalar("0x5de8159372063ce76403529bb1a9725461b96467035d906400ff48d0937f9db6"),
            parse_scalar("0x3c88b37dc6d14d08960b6fe58344e09194d11a930ce9f60cc90294683fac4b9f"),
        ];
        assert!(
            test_perm_hw::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse_hash("0x469e99663f620884b2749346c83e523a3ccc726718618ff752257c8761b6be3d")
            )
            .is_ok()
        );
        assert!(
            test_perm_hw::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse_hash("0xef9bd00db1232aff06bf42592ecc6a494bdb0e9c0e159735112af67ac7699994")
            )
            .is_ok()
        );
        assert!(
            test_perm_hw::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse_hash("0xdf8e644ab3c7a89b596539436071b75a11eb5a3dd7f196cf137e7580eb50b03c")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_hw() {
        let inputs = [from_const(0), from_const(1), from_const(2), from_const(3)];
        let outputs = [
            parse_scalar("0x775049834d9decb40ec5a109116a27527fa9105a3521cee8a42777788fda1501"),
            parse_scalar("0x630ded08b39ceac4859c9ab6d14b548f48d01164ce1efada3a7a868f7d9248cb"),
            parse_scalar("0x14b47f414dececb9936dcbb89e2fdd8511c44acb30439d1d23e48119b1c03b4f"),
            parse_scalar("0x72de70292ce1ac7f30b859d04bbb6de5377288c1192a08863c34e11bc9269c4c"),
        ];
        assert!(
            test_perm_hw::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse_hash("0x27024531ca1f9f33e2c5985bfc0c5267b264ebc568e5aa498e93ebaf70f4ce1e")
            )
            .is_ok()
        );
        assert!(
            test_perm_hw::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse_hash("0xc5abc7a6f03c29acfbef015ad9fa00fb6c501cec83da197ef7f6f9521aec275b")
            )
            .is_ok()
        );
        assert!(
            test_perm_hw::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse_hash("0x2371fd1d96c3e9ff6d6a0cb0bf4b5e82f9c48a1a503fa3eaffebb0e56d342498")
            )
            .is_ok()
        );
    }

    fn test_perm_ir<
        Cfg: poseidon2::Config<Scalar, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [Scalar; T],
        expected_output: [Scalar; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipIR::<Cfg, T>::default();
        assert_eq!(chip.width(), T * 2);
        test_permutation_impl::<T>(
            &chip,
            inputs,
            expected_output,
            blowup_log2,
            circuit_commitment,
        )
    }

    #[test]
    fn test_permutation_t3_ir() {
        let inputs = [from_const(0), from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x6f30582cde48a25b26015b7f718ba2fb359e93029caf04d8d0b3e66b1d46b941"),
            parse_scalar("0x5de8159372063ce76403529bb1a9725461b96467035d906400ff48d0937f9db6"),
            parse_scalar("0x3c88b37dc6d14d08960b6fe58344e09194d11a930ce9f60cc90294683fac4b9f"),
        ];
        assert!(
            test_perm_ir::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse_hash("0x068d62a5c2772b353ab5d5fa68175fbde0583aa83e9bf4b9ac01fc37b661f565")
            )
            .is_ok()
        );
        assert!(
            test_perm_ir::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse_hash("0x80683e60d6c01ecb352e66d04509c10a5ab17b0b01ab3c2a7b238719c8c9933f")
            )
            .is_ok()
        );
        assert!(
            test_perm_ir::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse_hash("0x83a9a9d4fce7a821e763a7d8fa05724952586efe0e6f193013c75696abfcade2")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_ir() {
        let inputs = [from_const(0), from_const(1), from_const(2), from_const(3)];
        let outputs = [
            parse_scalar("0x775049834d9decb40ec5a109116a27527fa9105a3521cee8a42777788fda1501"),
            parse_scalar("0x630ded08b39ceac4859c9ab6d14b548f48d01164ce1efada3a7a868f7d9248cb"),
            parse_scalar("0x14b47f414dececb9936dcbb89e2fdd8511c44acb30439d1d23e48119b1c03b4f"),
            parse_scalar("0x72de70292ce1ac7f30b859d04bbb6de5377288c1192a08863c34e11bc9269c4c"),
        ];
        assert!(
            test_perm_ir::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse_hash("0xff432e4f6af66723bd76dd61676367a1c403e291c80478662ece1971d43e7551")
            )
            .is_ok()
        );
        assert!(
            test_perm_ir::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse_hash("0x6d69d24edddae3b2ca17a690e48edef60e9d645ac2ba1411169869ece1a977bc")
            )
            .is_ok()
        );
        assert!(
            test_perm_ir::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse_hash("0xf588bf93caeadebc8e4db49654fa36d2720b246c46c815d096f4851c30951918")
            )
            .is_ok()
        );
    }

    fn test_perm_er<
        Cfg: poseidon2::Config<Scalar, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [Scalar; T],
        expected_output: [Scalar; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip_ir = PermutationChipIR::<Cfg, T>::default();
        assert_eq!(chip_ir.width(), T * 2);
        assert_eq!(chip_ir.height(), 194);
        let ir_width = chip_ir.width();

        let chip_er = PermutationChipER::<Cfg, T>::new(0, -(ir_width as isize));
        assert_eq!(chip_er.width(), T * 2);
        assert_eq!(chip_er.height(), 194);
        let er_width = chip_er.width();

        let mut builder = CircuitBuilder::default();
        let ir_output = builder.sub_chip(0, 0, &chip_ir, std::array::from_fn(|_| None))?;
        let er_output = builder.sub_chip(0, ir_width, &chip_er, std::array::from_fn(|_| None))?;

        for i in 0..T {
            builder.connect(ir_output[i], er_output[i]);
        }
        builder.declare_public_rows([ir_output[0].unwrap().row()]);

        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 194);
        assert_eq!(circuit.num_columns(), ir_width + er_width);

        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 194);
        assert_eq!(witness.num_columns(), ir_width + er_width);
        let ir_output = witness.sub_chip(0, 0, &chip_ir, inputs.map(|input| input.into()))?;
        let er_output =
            witness.sub_chip(0, ir_width, &chip_er, inputs.map(|input| input.into()))?;

        circuit.check_witness(&witness).unwrap();

        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, options.clone())?;

        let circuit = circuit.to_compressed::<Sha2Hash<Scalar>>(options);
        assert_eq!(circuit.commitment(), circuit_commitment);

        let public_inputs = circuit.verify(&proof)?;
        let get_value = |output: CellOrUnconstrained| match output {
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
        let inputs = [from_const(0), from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x6f30582cde48a25b26015b7f718ba2fb359e93029caf04d8d0b3e66b1d46b941"),
            parse_scalar("0x5de8159372063ce76403529bb1a9725461b96467035d906400ff48d0937f9db6"),
            parse_scalar("0x3c88b37dc6d14d08960b6fe58344e09194d11a930ce9f60cc90294683fac4b9f"),
        ];
        assert!(
            test_perm_er::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse_hash("0xa85862db285e33a3f81876151ba4fa97a5df904030f64f90b6abf40299271549")
            )
            .is_ok()
        );
        assert!(
            test_perm_er::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse_hash("0xf68a8d79c16fc97cdbc79e3fede5f99b6e1dba598b1b3ecdf228b4ac5846da99")
            )
            .is_ok()
        );
        assert!(
            test_perm_er::<poseidon2::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse_hash("0xa87740b722f7a353a0a7a14681be872e97acae0409bf6a64fe7c8b92984ecc0c")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_er() {
        let inputs = [from_const(0), from_const(1), from_const(2), from_const(3)];
        let outputs = [
            parse_scalar("0x775049834d9decb40ec5a109116a27527fa9105a3521cee8a42777788fda1501"),
            parse_scalar("0x630ded08b39ceac4859c9ab6d14b548f48d01164ce1efada3a7a868f7d9248cb"),
            parse_scalar("0x14b47f414dececb9936dcbb89e2fdd8511c44acb30439d1d23e48119b1c03b4f"),
            parse_scalar("0x72de70292ce1ac7f30b859d04bbb6de5377288c1192a08863c34e11bc9269c4c"),
        ];
        assert!(
            test_perm_er::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse_hash("0x4bf1b953329fba8611cc69fca7b7843c64149edc8b52f0089111fda90cdd84b6")
            )
            .is_ok()
        );
        assert!(
            test_perm_er::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse_hash("0x412ce8a695f4fdbaa36c1f1cd7ba866a7bce77f1d1b6be0a752d538c8c8c904f")
            )
            .is_ok()
        );
        assert!(
            test_perm_er::<poseidon2::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse_hash("0x74e40b6c9f5527f604ccd1641c128217f4790f5634ec5c466185a88622958059")
            )
            .is_ok()
        );
    }
}
