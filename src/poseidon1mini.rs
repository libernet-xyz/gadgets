use starkom_ff::{Field256, PrimeField};
use starkom_plonk::{
    Chip as PlonkChip, CircuitView, Constraint, WitnessView, make_const, rvar, var,
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
                    .map(|j| make_const(m[i * T + j]) * ((var(j) + var(T + j)) ^ a))
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

/// Compact Poseidon permutation chip.
///
/// You may want to use [`PermutationChipHW`], [`PermutationChipIR`], or [`PermutationChipER`]
/// rather than referring to this struct directly.
///
/// The difference between this chip and [`crate::poseidon1::PermutationChip`] is that this one uses
/// a more compact layout but its constraints have degree [`F::ALPHA`](`PrimeField::ALPHA`), while
/// the latter uses three times as many rows but caps the constraint degree at 3.
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
        inputs: [Option<starkom_plonk::Cell>; T],
    ) -> anyhow::Result<[Option<starkom_plonk::Cell>; T]>
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
        inputs: [starkom_plonk::CellOrUnconstrained<F>; T],
    ) -> anyhow::Result<[starkom_plonk::CellOrUnconstrained<F>; T]> {
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

// /// Poseidon permutation chip with [external ROM storage for round constants](`RcModeExternalRom`).
// pub type PermutationChipER<F, C, const T: usize> =
//     PermutationChip<F, C, RcModeExternalRom<F, C, T>, T>;

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use primitive_types::H256;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CellOrUnconstrained, CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon as poseidon1;
    use starkom_schraderbrau::Scalar as SB;
    use std::str::FromStr;

    fn parse<T: FromStr<Err: Debug>>(s: &'static str) -> T {
        s.parse().unwrap()
    }

    fn test_permutation_schraderbrau_impl<const T: usize>(
        chip: &impl PlonkChip<SB, T, T>,
        inputs: [SB; T],
        expected_output: [SB; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let height = match T {
            3 => 92,
            4 => 93,
            _ => unimplemented!(),
        };
        assert_eq!(chip.height(), height);
        let mut builder = CircuitBuilder::<SB, SB>::default();
        let output = builder.sub_chip(0, 0, chip, std::array::from_fn(|_| None))?;
        builder.declare_public_cells(output.into_iter().flatten());
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), height);
        assert_eq!(circuit.degree_bound(), 128);
        assert_eq!(circuit.num_columns(), chip.width());
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), height);
        assert_eq!(witness.degree_bound(), 128);
        assert_eq!(witness.num_columns(), chip.width());
        let output = witness.sub_chip(0, 0, chip, inputs.map(|input| input.into()))?;
        circuit.check_witness(&witness).unwrap();
        let options = ProvingOptions { blowup_log2 };
        let proof = circuit.prove::<Sha2Hash<SB>>(witness, options.clone())?;
        assert_eq!(proof.degree_bound(), 128);
        assert_eq!(proof.blowup_log2(), blowup_log2);
        assert_eq!(proof.extended_domain_size(), 128 << blowup_log2);
        let circuit = circuit.to_compressed::<Sha2Hash<SB>>(options);
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

    fn test_perm_schraderbrau_hw<
        Cfg: poseidon1::Config<SB, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [SB; T],
        expected_output: [SB; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipHW::<SB, Cfg, T>::default();
        assert_eq!(chip.width(), T);
        test_permutation_schraderbrau_impl::<T>(
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
            parse("0x1531e124d8f663b8d26c56b9bf0b1b09ba15d4e20ea10a3d898cdcf1d2c41dba"),
            parse("0x1f5714cc13f8a33f4b32a07f1a85409de2d1b8d353aa269ca8f2bfa91071bd62"),
            parse("0x588c20c682f8b66c52049d910a4f0e7195dfd300a0d0cf3198b1383be1a4134a"),
        ];
        assert!(
            test_perm_schraderbrau_hw::<poseidon1::SchraderbrauConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse("0xae1b4ac18fed8a2b8145f17b870343017f45b5e1242d93ac6c6f7170673029c6")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_hw::<poseidon1::SchraderbrauConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0xd57cecde6c2e2500dd273b0b0fb091cae6bf9bf65028319423d87a97386920ef")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_hw::<poseidon1::SchraderbrauConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0xb2f743a95c4b832dd1aa2f004c68d07cafa46e171f2bd17f23c24b70180eb741")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_hw() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x11de0b3702563747b0729abd2b93e720ec947ce067ca3f8c088b9829df1169d7"),
            parse("0x487ca87fa054ad960f6f571cf7aa6f5e075fd5ac3386e2e9fb22aeeb036fb1e1"),
            parse("0x28c775ddba9039531ee5deb0d4e2e6b1c9bb7d8be8da260a20d331bb041d0d38"),
            parse("0x10c93f76dca6f63507bc8ed1d2ea40647bb480ad43ac9922a3f10597b802f949"),
        ];
        assert!(
            test_perm_schraderbrau_hw::<poseidon1::SchraderbrauConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse("0x6f5311bfa48fa99f0941ee57070002f063e2989690df99add464a568992b2a34")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_hw::<poseidon1::SchraderbrauConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0x6c124ce9258ef91862e434f6fa47c98b0c051f991a6c6e0a1759d17fa9b6b734")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_hw::<poseidon1::SchraderbrauConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0xd161f3e6ba8cf5e8d73aa3645c44b9082254d8c3985944df254a47ac841defc0")
            )
            .is_ok()
        );
    }

    fn test_perm_schraderbrau_ir<
        Cfg: poseidon1::Config<SB, T>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [SB; T],
        expected_output: [SB; T],
        blowup_log2: usize,
        circuit_commitment: H256,
    ) -> Result<()> {
        let chip = PermutationChipIR::<SB, Cfg, T>::default();
        assert_eq!(chip.width(), T * 2);
        test_permutation_schraderbrau_impl::<T>(
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
            parse("0x1531e124d8f663b8d26c56b9bf0b1b09ba15d4e20ea10a3d898cdcf1d2c41dba"),
            parse("0x1f5714cc13f8a33f4b32a07f1a85409de2d1b8d353aa269ca8f2bfa91071bd62"),
            parse("0x588c20c682f8b66c52049d910a4f0e7195dfd300a0d0cf3198b1383be1a4134a"),
        ];
        assert!(
            test_perm_schraderbrau_ir::<poseidon1::SchraderbrauConfig3, 3, 2, 1>(
                inputs,
                outputs,
                1,
                parse("0xd6f0fc99e8190a2b51e44b37111d37963afaff4ae49bb2becb5aee59714ddf0e")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_ir::<poseidon1::SchraderbrauConfig3, 3, 2, 1>(
                inputs,
                outputs,
                2,
                parse("0x8bbb1d9305dc9d61d43ffe2bed090e19048deb4aa5da594df3d50528994d84f8")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_ir::<poseidon1::SchraderbrauConfig3, 3, 2, 1>(
                inputs,
                outputs,
                3,
                parse("0x87fec259ec7c97d39145625dc926c563e09f46e49e801d803713324f076fde48")
            )
            .is_ok()
        );
    }

    #[test]
    fn test_permutation_t4_ir() {
        let inputs = [0u8.into(), 1u8.into(), 2u8.into(), 3u8.into()];
        let outputs = [
            parse("0x11de0b3702563747b0729abd2b93e720ec947ce067ca3f8c088b9829df1169d7"),
            parse("0x487ca87fa054ad960f6f571cf7aa6f5e075fd5ac3386e2e9fb22aeeb036fb1e1"),
            parse("0x28c775ddba9039531ee5deb0d4e2e6b1c9bb7d8be8da260a20d331bb041d0d38"),
            parse("0x10c93f76dca6f63507bc8ed1d2ea40647bb480ad43ac9922a3f10597b802f949"),
        ];
        assert!(
            test_perm_schraderbrau_ir::<poseidon1::SchraderbrauConfig4, 4, 3, 1>(
                inputs,
                outputs,
                1,
                parse("0xd8873b7baf50db8c6b44f3e39b85863fa5667b0ace57bcd7bb5357a2bb40f850")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_ir::<poseidon1::SchraderbrauConfig4, 4, 3, 1>(
                inputs,
                outputs,
                2,
                parse("0xd7c92583699b214e60ac8eec659d47c22c251f25180eb309c946a484c6a9f218")
            )
            .is_ok()
        );
        assert!(
            test_perm_schraderbrau_ir::<poseidon1::SchraderbrauConfig4, 4, 3, 1>(
                inputs,
                outputs,
                3,
                parse("0x8136c3a1a9553f7c207a2cf8efbe45b962e7cc0a7bcaf71643f266601d9b05ea")
            )
            .is_ok()
        );
    }

    // TODO
}
