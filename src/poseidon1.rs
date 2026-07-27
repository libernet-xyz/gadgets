use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, rvar,
};
use starkom_poseidon as poseidon;
use std::marker::PhantomData;

#[derive(Debug, Copy, Clone)]
pub struct PermutationChip<
    Cfg: poseidon::Config<Scalar, T, R, C>,
    const T: usize,
    const R: usize,
    const C: usize,
> {
    _data: PhantomData<Cfg>,
}

impl<Cfg: poseidon::Config<Scalar, T, R, C>, const T: usize, const R: usize, const C: usize> Default
    for PermutationChip<Cfg, T, R, C>
{
    fn default() -> Self {
        Self {
            _data: PhantomData::default(),
        }
    }
}

impl<Cfg: poseidon::Config<Scalar, T, R, C>, const T: usize, const R: usize, const C: usize>
    PermutationChip<Cfg, T, R, C>
{
    const ROUND_HEIGHT: usize = 3;

    fn build_first_arc(&self, view: &mut impl CircuitView, inputs: [Option<Cell>; T]) {
        for i in 0..T {
            view.connect(inputs[i], Some(view.cell(0, i)));
        }
        let c = Cfg::get_round_constants();
        for i in 0..T {
            view.add_gate(0, rvar(i, 0) + c[i] - rvar(i, 1));
        }
    }

    fn witness_first_arc(&self, view: &mut impl WitnessView, inputs: [CellOrUnconstrained; T]) {
        for i in 0..T {
            view.copy(inputs[i], view.cell(0, i));
        }
        let c = Cfg::get_round_constants();
        for i in 0..T {
            let state = view.get(view.cell(0, i));
            view.set(view.cell(1, i), state + c[i]);
        }
    }

    fn build_full_sbox(&self, view: &mut impl CircuitView) {
        for i in 0..T {
            view.add_gate(0, (rvar(i, -1) ^ 3) - rvar(i, 0));
            view.add_gate(0, (rvar(i, -1) ^ 2) * rvar(i, 0) - rvar(i, 1));
        }
    }

    fn witness_full_sbox(&self, view: &mut impl WitnessView) {
        for i in 0..T {
            let state = view.get(view.cell(-1, i));
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
        let state = view.get(view.cell(-1, 0));
        view.set(view.cell(0, 0), state.cube());
        view.set(view.cell(1, 0), state.square().square() * state);
        for i in 1..T {
            view.copy(view.cell(-1, i).into(), view.cell(1, i));
        }
    }

    fn build_mds_and_next_arc(&self, view: &mut impl CircuitView, round: usize) {
        let c = Cfg::get_round_constants();
        let m = Cfg::get_mds_matrix();
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

    fn witness_mds_and_next_arc(&self, view: &mut impl WitnessView, round: usize) {
        let c = Cfg::get_round_constants();
        let m = Cfg::get_mds_matrix();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>()
                    + c[(round + 1) * T + i],
            );
        }
    }

    fn build_last_mds(&self, view: &mut impl CircuitView) {
        let m = Cfg::get_mds_matrix();
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

    fn witness_last_mds(&self, view: &mut impl WitnessView) {
        let m = Cfg::get_mds_matrix();
        for i in 0..T {
            view.set(
                view.cell(1, i),
                (0..T)
                    .map(|j| view.get(view.cell(0, j)) * m[i * T + j])
                    .sum::<Scalar>(),
            );
        }
    }
}

impl<Cfg: poseidon::Config<Scalar, T, R, C>, const T: usize, const R: usize, const C: usize>
    PlonkChip<T, T> for PermutationChip<Cfg, T, R, C>
{
    fn width(&self) -> usize {
        T
    }

    fn build(
        &self,
        view: &mut impl CircuitView,
        inputs: [Option<Cell>; T],
    ) -> Result<[Option<Cell>; T]> {
        let num_full_rounds = Cfg::num_full_rounds();
        let num_partial_rounds = Cfg::num_partial_rounds();
        let num_total_rounds = Cfg::num_total_rounds();
        assert_eq!(num_total_rounds, num_full_rounds * 2 + num_partial_rounds);
        self.build_first_arc(view, inputs);
        let mut view = view.sub(2, 0, T);
        for r in 0..num_full_rounds {
            view.sub(r * Self::ROUND_HEIGHT, 0, T)
                .sub_fn(0, 0, T, |view| self.build_full_sbox(view))
                .sub_fn(1, 0, T, |view| self.build_mds_and_next_arc(view, r));
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub(r * Self::ROUND_HEIGHT, 0, T)
                .sub_fn(0, 0, T, |view| self.build_partial_sbox(view))
                .sub_fn(1, 0, T, |view| self.build_mds_and_next_arc(view, r));
        }
        for r in (num_full_rounds + num_partial_rounds)..(num_total_rounds - 1) {
            view.sub(r * Self::ROUND_HEIGHT, 0, T)
                .sub_fn(0, 0, T, |view| self.build_full_sbox(view))
                .sub_fn(1, 0, T, |view| self.build_mds_and_next_arc(view, r));
        }
        view.sub((num_total_rounds - 1) * Self::ROUND_HEIGHT, 0, T)
            .sub_fn(0, 0, T, |view| self.build_full_sbox(view))
            .sub_fn(1, 0, T, |view| self.build_last_mds(view));
        Ok(std::array::from_fn(|i| {
            Some(view.cell(num_total_rounds * Self::ROUND_HEIGHT - 1, i))
        }))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView,
        inputs: [CellOrUnconstrained; T],
    ) -> Result<[CellOrUnconstrained; T]> {
        let num_full_rounds = Cfg::num_full_rounds();
        let num_partial_rounds = Cfg::num_partial_rounds();
        let num_total_rounds = Cfg::num_total_rounds();
        assert_eq!(num_total_rounds, num_full_rounds * 2 + num_partial_rounds);
        self.witness_first_arc(view, inputs);
        let mut view = view.sub(2, 0, T);
        for r in 0..num_full_rounds {
            view.sub(r * Self::ROUND_HEIGHT, 0, T)
                .sub_fn(0, 0, T, |view| self.witness_full_sbox(view))
                .sub_fn(1, 0, T, |view| self.witness_mds_and_next_arc(view, r));
        }
        for r in num_full_rounds..(num_full_rounds + num_partial_rounds) {
            view.sub(r * Self::ROUND_HEIGHT, 0, T)
                .sub_fn(0, 0, T, |view| self.witness_partial_sbox(view))
                .sub_fn(1, 0, T, |view| self.witness_mds_and_next_arc(view, r));
        }
        for r in (num_full_rounds + num_partial_rounds)..(num_total_rounds - 1) {
            view.sub(r * Self::ROUND_HEIGHT, 0, T)
                .sub_fn(0, 0, T, |view| self.witness_full_sbox(view))
                .sub_fn(1, 0, T, |view| self.witness_mds_and_next_arc(view, r));
        }
        view.sub((num_total_rounds - 1) * Self::ROUND_HEIGHT, 0, T)
            .sub_fn(0, 0, T, |view| self.witness_full_sbox(view))
            .sub_fn(1, 0, T, |view| self.witness_last_mds(view));
        Ok(std::array::from_fn(|i| {
            view.cell(num_total_rounds * Self::ROUND_HEIGHT - 1, i)
                .into()
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, CompilationOptions, ProvingOptions};
    use starkom_poseidon as poseidon1;

    fn test_permutation<
        Cfg: poseidon1::Config<Scalar, T, R, C>,
        const T: usize,
        const R: usize,
        const C: usize,
    >(
        inputs: [Scalar; T],
        expected_output: [Scalar; T],
        blowup_log2: usize,
    ) -> Result<()> {
        let chip = PermutationChip::<Cfg, T, R, C>::default();
        assert_eq!(chip.width(), T);
        let mut builder = CircuitBuilder::default();
        let output = chip.build(&mut builder, std::array::from_fn(|_| None))?;
        builder.declare_public_rows([output[0].unwrap().row()]);
        let circuit = builder.build(CompilationOptions {
            canonicalize_constraints: false,
        })?;
        assert_eq!(circuit.num_rows(), 194);
        assert_eq!(circuit.degree_bound(), 256);
        assert_eq!(circuit.num_columns(), T);
        let mut witness = circuit.make_witness();
        assert_eq!(witness.num_rows(), 194);
        assert_eq!(witness.degree_bound(), 256);
        assert_eq!(witness.num_columns(), T);
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
    fn test_permutation_t3() {
        let inputs = [from_const(0), from_const(1), from_const(2)];
        let outputs = [
            parse_scalar("0x7b68dcd80fa751ee8f2d76043bfd92c685601c79189393fc76e03c5214eed32b"),
            parse_scalar("0x0fbcb5720b463bf7e2ccabf373e77d2c10d27e6549f34cfa33eb2d06ea8b900a"),
            parse_scalar("0x26e03abfcc62da0101516b07aede8bc676a10c47299a57bedc6d9fe80484f3da"),
        ];
        assert!(test_permutation::<poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 1).is_ok());
        assert!(test_permutation::<poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 2).is_ok());
        assert!(test_permutation::<poseidon1::BlueSkyConfig3, 3, 2, 1>(inputs, outputs, 3).is_ok());
    }

    #[test]
    fn test_permutation_t4() {
        let inputs = [from_const(0), from_const(1), from_const(2), from_const(3)];
        let outputs = [
            parse_scalar("0x12dde8a4c46760e349670d241e36ca7abacc991233039f8deaf6c58ce2230ef6"),
            parse_scalar("0x61e95d9456e9223b4d7926dabae10009da2b6fb9134ade8405f6ef1424e66aa1"),
            parse_scalar("0x2fcce25ab9efb3e26276f3b3aff1e02cdf82df48ce8d3eadbff900cfe015775b"),
            parse_scalar("0x2580707d57a8c1c0cad368e8d5705ffd96f269d66e1cd6f1433f93a3c66d9bf8"),
        ];
        assert!(test_permutation::<poseidon1::BlueSkyConfig4, 4, 3, 1>(inputs, outputs, 1).is_ok());
        assert!(test_permutation::<poseidon1::BlueSkyConfig4, 4, 3, 1>(inputs, outputs, 2).is_ok());
        assert!(test_permutation::<poseidon1::BlueSkyConfig4, 4, 3, 1>(inputs, outputs, 3).is_ok());
    }
}
