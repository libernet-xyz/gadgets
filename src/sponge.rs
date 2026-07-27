use anyhow::Result;
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
                    view.connect(input, Some(view.cell(0, i)));
                }
                None => {
                    view.add_gate(0, var(i));
                }
            }
        }
        for i in 0..C {
            view.add_gate(0, var(R + i));
        }
        for i in 0..T {
            view.add_gate(1, rvar(i, -1) + rvar(i, 0) - rvar(i, 1));
        }
        std::array::from_fn(|i| Some(view.cell(2, i)))
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

        // TODO
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // TODO
}
