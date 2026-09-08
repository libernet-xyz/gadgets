use crate::poseidon::Sbox;
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

        fn build_first_arc<G: Field256 + From<F>>(
            &self,
            view: &mut impl CircuitView<F, G>,
            inputs: [Option<Cell>; T],
        ) where
            F: Mul<G, Output = G>,
            G: Mul<F, Output = G>;

        fn witness_first_arc(
            &self,
            view: &mut impl WitnessView<F>,
            inputs: [CellOrUnconstrained<F>; T],
        );

        fn build_mds_and_next_arc<G: Field256 + From<F>>(
            &self,
            view: &mut impl CircuitView<F, G>,
            round: usize,
        ) where
            F: Mul<G, Output = G>,
            G: Mul<F, Output = G>;

        fn witness_mds_and_next_arc(&self, view: &mut impl WitnessView<F>, round: usize);
    }
}

/// Hard-wired RC mode for the Poseidon1 permutation.
///
/// In this mode the chip takes exactly `T` columns, with `T` being the state vector size, but uses
/// a different gate for every ARC layer. You should use this mode only if the elevated number of
/// gates is not a concern.
///
/// This mode is best suited for circuits that runs a small number of hashes, such as preimage
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

    fn build_first_arc<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; T],
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        for i in 0..T {
            view.connect(inputs[i], Some(view.cell(0, i)));
        }
        let c = C::get_round_constants();
        for i in 0..T {
            view.add_gate(0, rvar(i, 0) + make_const(c[i]) - rvar(i, 1));
        }
    }

    fn witness_first_arc(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; T],
    ) {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let c = C::get_round_constants();
        for i in 0..T {
            let state = view.get_at(view.cell(0, i));
            view.set(view.cell(1, i), state + c[i]);
        }
    }

    fn build_mds_and_next_arc<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let c = C::get_round_constants();
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * make_const(m[i * T + j]))
                    .sum::<Constraint<F>>()
                    + make_const(c[(round + 1) * T + i])
                    - rvar(i, 1),
            );
        }
    }

    fn witness_mds_and_next_arc(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = C::get_round_constants();
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<F>()
                    + c[(round + 1) * T + i],
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

    fn build_first_arc<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; T],
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        for i in 0..T {
            view.connect(inputs[i], Some(view.cell(0, i)));
        }
        let c = C::get_round_constants();
        for i in 0..T {
            view.add_gate(1, var(T + i) - make_const(c[i]));
            view.add_gate(0, rvar(i, 0) + rvar(T + i, 1) - rvar(i, 1));
        }
    }

    fn witness_first_arc(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; T],
    ) {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let c = C::get_round_constants();
        for i in 0..T {
            let state = view.get_at(view.cell(0, i));
            view.set(view.cell(1, i), state + c[i]);
            view.set(view.cell(1, T + i), c[i]);
        }
    }

    fn build_mds_and_next_arc<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let c = C::get_round_constants();
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.add_gate(1, var(T + i) - make_const(c[(round + 1) * T + i]));
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * make_const(m[i * T + j]))
                    .sum::<Constraint<F>>()
                    + rvar(T + i, 1)
                    - rvar(i, 1),
            );
        }
    }

    fn witness_mds_and_next_arc(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = C::get_round_constants();
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(view.cell(1, T + i), c[(round + 1) * T + i]);
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<F>()
                    + c[(round + 1) * T + i],
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
            self.ir_chip_row_offset + 1,
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

    fn build_first_arc<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; T],
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        for i in 0..T {
            view.connect(inputs[i], Some(view.cell(0, i)));
        }
        for i in 0..T {
            view.connect(
                self.remote_rom_cell(view, i).into(),
                view.cell(1, T + i).into(),
            );
            view.add_gate(0, rvar(i, 0) + rvar(T + i, 1) - rvar(i, 1));
        }
    }

    fn witness_first_arc(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; T],
    ) {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let c = C::get_round_constants();
        for i in 0..T {
            let state = view.get_at(view.cell(0, i));
            view.set(view.cell(1, i), state + c[i]);
            view.set(view.cell(1, T + i), c[i]);
        }
    }

    fn build_mds_and_next_arc<G: Field256 + From<F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        _round: usize,
    ) where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.connect(
                self.remote_rom_cell(view, i).into(),
                view.cell(1, T + i).into(),
            );
            view.add_gate(
                0,
                (0..T)
                    .map(|j| rvar(j, 0) * make_const(m[i * T + j]))
                    .sum::<Constraint<F>>()
                    + rvar(T + i, 1)
                    - rvar(i, 1),
            );
        }
    }

    fn witness_mds_and_next_arc(&self, view: &mut impl WitnessView<F>, round: usize) {
        let c = C::get_round_constants();
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(view.cell(1, T + i), c[(round + 1) * T + i]);
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<F>()
                    + c[(round + 1) * T + i],
            );
        }
    }
}

/// Poseidon permutation chip.
///
/// You may want to use [`PermutationChipHW`], [`PermutationChipIR`], or [`PermutationChipER`]
/// rather than referring to this struct directly.
pub struct PermutationChip<
    F: PrimeField + Sbox,
    C: poseidon::Config<F, T>,
    M: internal::RcMode<F, T>,
    const T: usize,
> {
    rc: M,
    _data: PhantomData<(F, C)>,
}

impl<F: PrimeField + Sbox, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize>
    Debug for PermutationChip<F, C, M, T>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PermutationChip")
            .field("rc", &self.rc)
            .finish()
    }
}

impl<
    F: PrimeField + Sbox,
    C: poseidon::Config<F, T>,
    M: internal::RcMode<F, T> + Default,
    const T: usize,
> Default for PermutationChip<F, C, M, T>
{
    fn default() -> Self {
        Self {
            rc: M::default(),
            _data: PhantomData,
        }
    }
}

impl<F: PrimeField + Sbox, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize>
    Copy for PermutationChip<F, C, M, T>
{
}

impl<F: PrimeField + Sbox, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize>
    Clone for PermutationChip<F, C, M, T>
{
    fn clone(&self) -> Self {
        Self {
            rc: self.rc.clone(),
            _data: PhantomData,
        }
    }
}

impl<F: PrimeField + Sbox, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize>
    PermutationChip<F, C, M, T>
{
    pub const FIRST_ARC_HEIGHT: usize = 2;
    pub const ROUND_HEIGHT: usize = 3;

    fn build_full_sbox<G: Field256 + From<F>>(&self, view: &mut impl CircuitView<F, G>)
    where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        for i in 0..T {
            view.sub_fn(0, i, Some(1), Some(2), |view| F::build_sbox::<G>(view));
        }
    }

    fn witness_full_sbox(&self, view: &mut impl WitnessView<F>) {
        for i in 0..T {
            view.sub_fn(0, i, Some(1), Some(2), |view| F::witness_sbox(view));
        }
    }

    fn build_partial_sbox<G: Field256 + From<F>>(&self, view: &mut impl CircuitView<F, G>)
    where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        view.sub_fn(0, 0, Some(1), Some(2), |view| F::build_sbox::<G>(view));
        for i in 1..T {
            view.connect(Some(view.cell(-1, i)), Some(view.cell(1, i)));
        }
    }

    fn witness_partial_sbox(&self, view: &mut impl WitnessView<F>) {
        view.sub_fn(0, 0, Some(1), Some(2), |view| F::witness_sbox(view));
        for i in 1..T {
            view.copy(view.cell(-1, i).into(), view.cell(1, i));
        }
    }

    fn build_last_mds<G: Field256 + From<F>>(&self, view: &mut impl CircuitView<F, G>)
    where
        F: Mul<G, Output = G>,
        G: Mul<F, Output = G>,
    {
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.add_gate(
                0,
                ((0..T)
                    .map(|j| rvar(j, 0) * make_const(m[i * T + j]))
                    .sum::<Constraint<F>>())
                    - rvar(i, 1),
            );
        }
    }

    fn witness_last_mds(&self, view: &mut impl WitnessView<F>) {
        let m = C::get_mds_matrix();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get_at(view.cell(0, j)) * m[i * T + j])
                    .sum::<F>(),
            );
        }
    }
}

impl<F: PrimeField + Sbox, C: poseidon::Config<F, T>, const T: usize>
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

impl<F: PrimeField + Sbox, C: poseidon::Config<F, T>, M: internal::RcMode<F, T>, const T: usize>
    PlonkChip<F, T, T> for PermutationChip<F, C, M, T>
{
    fn width(&self) -> usize {
        self.rc.width()
    }

    fn height(&self) -> usize {
        Self::FIRST_ARC_HEIGHT + Self::ROUND_HEIGHT * C::num_total_rounds()
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
        self.rc.build_first_arc(view, inputs);
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
                    self.rc.build_mds_and_next_arc(view, r)
                });
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.build_partial_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.build_mds_and_next_arc(view, r)
                });
        }
        for r in (num_full_rounds + num_partial_rounds)..(num_total_rounds - 1) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.build_full_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.build_mds_and_next_arc(view, r)
                });
        }
        view.sub(
            (num_total_rounds - 1) * Self::ROUND_HEIGHT,
            0,
            None,
            Some(Self::ROUND_HEIGHT),
        )
        .sub_fn(0, 0, None, None, |view| self.build_full_sbox(view))
        .sub_fn(1, 0, None, None, |view| self.build_last_mds(view));
        Ok(std::array::from_fn(|i| {
            Some(view.cell(num_total_rounds * Self::ROUND_HEIGHT - 1, i))
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
        self.rc.witness_first_arc(view, inputs);
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
                    self.rc.witness_mds_and_next_arc(view, r)
                });
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.witness_partial_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.witness_mds_and_next_arc(view, r)
                });
        }
        for r in (num_full_rounds + num_partial_rounds)..(num_total_rounds - 1) {
            view.sub(r * Self::ROUND_HEIGHT, 0, None, Some(Self::ROUND_HEIGHT))
                .sub_fn(0, 0, None, None, |view| self.witness_full_sbox(view))
                .sub_fn(1, 0, None, None, |view| {
                    self.rc.witness_mds_and_next_arc(view, r)
                });
        }
        view.sub(
            (num_total_rounds - 1) * Self::ROUND_HEIGHT,
            0,
            None,
            Some(Self::ROUND_HEIGHT),
        )
        .sub_fn(0, 0, None, None, |view| self.witness_full_sbox(view))
        .sub_fn(1, 0, None, None, |view| self.witness_last_mds(view));
        Ok(std::array::from_fn(|i| {
            view.cell(num_total_rounds * Self::ROUND_HEIGHT - 1, i)
                .into()
        }))
    }
}

/// Poseidon permutation chip with [hard-wired round constants](`RcModeHardWired`).
pub type PermutationChipHW<F, C, const T: usize> =
    PermutationChip<F, C, RcModeHardWired<F, C, T>, T>;

/// Poseidon permutation chip with [internal ROM storage for round constants](`RcModeInternalRom`).
pub type PermutationChipIR<F, C, const T: usize> =
    PermutationChip<F, C, RcModeInternalRom<F, C, T>, T>;

/// Poseidon permutation chip with [internal ROM storage for round constants](`RcModeExternalRom`).
pub type PermutationChipER<F, C, const T: usize> =
    PermutationChip<F, C, RcModeExternalRom<F, C, T>, T>;

#[cfg(all(test, feature = "bluesky", feature = "goldilocks"))]
mod tests {
    use super::*;
    use primitive_types::H256;
    use starkom_bluesky::Scalar as BS;
    use starkom_goldilocks::{GL, GL4, from_const};
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
        assert_eq!(chip.height(), 194);
        let mut builder = CircuitBuilder::<BS, BS>::default();
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
        let proof = circuit.prove::<Sha2Hash<BS>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 256);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 256 << blowup_log2);
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
        assert_eq!(chip.height(), 92);
        let mut builder = CircuitBuilder::<GL, GL4>::default();
        let output = builder.sub_chip(0, 0, chip, std::array::from_fn(|_| None))?;
        builder.declare_public_rows([output[0].unwrap().row()]);
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 92);
        assert_eq!(circuit.degree_bound(), 128);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 92);
        assert_eq!(witness.degree_bound(), 128);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<GL4>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 128);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 128 << blowup_log2);
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
                parse("0xcf0d11b627b09bbbddb45b72cbf61564c9966d4d0cd0392088810c72ef2d40f1")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0x801c8e968c58ef2e247315a2b0cc1cbfdeea10342fe4bec0ec827252b6f65f77")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0x581224f66337b1f48ae1938117fc8edb5345462614c68675661059a99df3f2ff")
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
                parse("0xf1732ed3f90394cc74ad7bb361ff1af683ebe10a18053c1f0b52df417f192edb")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0x0a0690dd5ee6ca4efa9505601fe3362d638735543e4196d49c0facf59a2b0f0f")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_hw::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0xff9c71a6794801b96710f41abaf4f2823ba29b153cc04dfddc11f712b167d342")
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
        let inputs = std::array::from_fn(|i| from_const(i as u64));
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
                parse("0xdbad34475f0d5f38cc11b9ab25936e3edbe3418dc2670745de140119f4987c82")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0x69adf9a8ffa5d637f4b3433d7d91f39102645754cc2cf4ea17e28894f0d468b8")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0xf2a96c0c58559129fa670ca2c6cbe6f7697cadded8c796e9c39167fb9edc6c0b")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t16_hw() {
        let inputs = std::array::from_fn(|i| from_const(i as u64));
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
                parse("0x5260ce7c2d728c1a41d19f53b5b2159326fadcb2db6c6ceca6b3addc3b37d70a")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x620c6a76c15f9bcf20d19bd8a3b8ada71898403ca9ba52c8d3b607c045d1a501")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_hw::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0x8236d782628179a9d407be3af27421e550a815252f000b17c8729f54eff665ce")
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
                parse("0xe4d769d9800824e77947b4bfc3245b40340160a133514a57b5843b0a205a4e29")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0x5b65a04b1aee48eced617acfc78e414b829dc9e6ff90e9964ec8b150387e67f7")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0xfc6b3292633d8898b39e28633c003db164ee8243628ccef32b3cd43137d94238")
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
                parse("0x24cd89be446211981f6968c0370a948c4382a9e8845db545462c6cbd2d2f9392")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0x0c80929bd4a4d878f9585c582adce9a7ffb5779afaee7902849e1f670cfe752b")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_ir::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0xb54ce41a7e8e1a408acb86c6b406c392c8731bb454c64d949c11838b707a28df")
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
        let inputs = std::array::from_fn(|i| from_const(i as u64));
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
                parse("0x0bcdce6774d6b947daceb0590ee91b0146bcfeadf7096b739604c0357e55f048")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0x96672484b088472576742f6fe48cef23e99cd239907177d5db546d619809c419")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0xa9b7d7a58feead08151d99649d559271796b23fb5466f86f7aa0a7376a12d937")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t16_ir() {
        let inputs = std::array::from_fn(|i| from_const(i as u64));
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
                parse("0xe65c79510dbbefbb2fde3583ae0f4aab4cd0fd4b8d12883297d0e097f58b39db")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0x483c454d97c79edc7b262d4b759db1a204cc779dd46ef2f0ad6b1d595d182b09")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_ir::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0xf2c27ea1dca488069739767022835c3f6cc6652c6f625d9ecf337771cdfd4a1c")
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
        assert_eq!(chip_ir.height(), 194);
        let ir_width = chip_ir.width();

        let chip_er = PermutationChipER::<BS, Cfg, T>::new(0, -(ir_width as isize));
        assert_eq!(chip_er.width(), T * 2);
        assert_eq!(chip_er.height(), 194);
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
                parse("0x97c79aee24f7fad057a23dd7b395072029836509b844becaa60512246429b1fc")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0xcefae940640b771aac1f0b5697a406a840f24a1df1b6e9f14dc239ee3932a02d")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0x6e065aa94c3348d7abb11c78e46b86f3baef0d24275fa78a0b93beeb677f3170")
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
                parse("0xd47f7c3a05ea2831ea1ad06a845a1cc32a825a5da1ffeab056a15cfd780a6b38")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0x0f346698864117f6ffceb1773a2c20d8581f3f12386af9c7360895c84ccc3abf")
            )
            .is_ok()
        );
        assert!(
            test_perm_bluesky_er::<poseidon1::BlueSkyConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0xd31a6155467a594029b6882e2c26c375b68116b26c821ad7010ab8324ce50f3c")
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
        assert_eq!(chip_ir.height(), 92);
        let ir_width = chip_ir.width();

        let chip_er = PermutationChipER::<GL, Cfg, T>::new(0, -(ir_width as isize));
        assert_eq!(chip_er.width(), T * 2);
        assert_eq!(chip_er.height(), 92);
        let er_width = chip_er.width();

        let mut builder = CircuitBuilder::<GL, GL4>::default();
        let ir_output = builder.sub_chip(0, 0, &chip_ir, std::array::from_fn(|_| None))?;
        let er_output = builder.sub_chip(0, ir_width, &chip_er, std::array::from_fn(|_| None))?;

        for i in 0..T {
            builder.connect(ir_output[i], er_output[i]);
        }
        builder.declare_public_rows([ir_output[0].unwrap().row()]);

        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 92);
        assert_eq!(circuit.num_columns(), ir_width + er_width);

        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 92);
        assert_eq!(witness.num_columns(), ir_width + er_width);
        let ir_output = witness.sub_chip(0, 0, &chip_ir, inputs.map(|input| input.into()))?;
        let er_output =
            witness.sub_chip(0, ir_width, &chip_er, inputs.map(|input| input.into()))?;

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
        let inputs = std::array::from_fn(|i| from_const(i as u64));
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
                parse("0x96d011ec1d620702853891a921ec6f14ec173ad3cf9afc4210e4ab097e57e6eb")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                2,
                parse("0xbfbfeefafa68d8b5377f1109a0892a2f4c7e3e0d4d5be7ebd0f9e3775176e31c")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig12, 12>(
                inputs,
                outputs,
                3,
                parse("0x32ad2105e5f4f221a1c4e06865b53950e66b3a069e00efcdd818f9eead1b3bb1")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t16_er() {
        let inputs = std::array::from_fn(|i| from_const(i as u64));
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
                parse("0x600e58e59bf0e664ac22446efbc4cadf712fd6cfc713cfaff4d9edfbd5ba2e6e")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                2,
                parse("0xd8f1b019ef5be35ddb0d44a24beb0fdac33a553ff118f95179499c41009a0dc0")
            )
            .is_ok()
        );
        assert!(
            test_perm_goldilocks_er::<poseidon1::GoldilocksConfig16, 16>(
                inputs,
                outputs,
                3,
                parse("0xebdb10fd8afcbf09701e3bb9586d93b670de4b4c92a22c6d280fd2d621d1f5e1")
            )
            .is_ok()
        );
    }
}
