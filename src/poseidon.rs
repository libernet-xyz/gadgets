use starkom_ff::{Field, Field256, PrimeField};
use starkom_plonk::{CircuitView, WitnessView, rvar};
use std::ops::Mul;

pub trait Sbox: PrimeField {
    fn build_sbox<G: Field256 + From<Self>>(view: &mut impl CircuitView<Self, G>)
    where
        Self: Mul<G, Output = G>,
        G: Mul<Self, Output = G>;

    fn witness_sbox(view: &mut impl WitnessView<Self>);
}

#[cfg(feature = "bluesky")]
impl Sbox for starkom_bluesky::Scalar {
    fn build_sbox<G: Field256 + From<Self>>(view: &mut impl CircuitView<Self, G>)
    where
        Self: Mul<G, Output = G>,
        G: Mul<Self, Output = G>,
    {
        view.add_gate(0, (rvar(0, -1) ^ 3) - rvar(0, 0));
        view.add_gate(0, (rvar(0, -1) ^ 2) * rvar(0, 0) - rvar(0, 1));
    }

    fn witness_sbox(view: &mut impl WitnessView<Self>) {
        let state = view.get_at(view.cell(-1, 0));
        view.set(view.cell(0, 0), state.cube());
        view.set(view.cell(1, 0), state.square().square() * state);
    }
}

#[cfg(feature = "goldilocks")]
impl Sbox for starkom_goldilocks::GL {
    fn build_sbox<G: Field256 + From<Self>>(view: &mut impl CircuitView<Self, G>)
    where
        Self: Mul<G, Output = G>,
        G: Mul<Self, Output = G>,
    {
        view.add_gate(0, (rvar(0, -1) ^ 3) - rvar(0, 0));
        view.add_gate(0, (rvar(0, 0) ^ 2) * rvar(0, -1) - rvar(0, 1));
    }

    fn witness_sbox(view: &mut impl WitnessView<Self>) {
        let state = view.get_at(view.cell(-1, 0));
        let cube = state.cube();
        view.set(view.cell(0, 0), cube);
        view.set(view.cell(1, 0), cube.square() * state);
    }
}
