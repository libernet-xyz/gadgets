use anyhow::{Result, anyhow};
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{Chip as PlonkChip, CircuitBuilder, Wire, WireOrUnconstrained, Witness};
use starkom_poseidon2::{Config, bluesky::BlueSkyConfig, bluesky::BlueSkyConfig4};

/// PLONK chip for our Poseidon hash instance over the BlueSky field.
///
/// `T` is the state vector size and `I` is the number of input scalars.
#[derive(Debug, Default, Clone)]
pub struct Chip<const T: usize, const I: usize> {}

impl<const T: usize, const I: usize> Chip<T, I>
where
    BlueSkyConfig<T>: Config<Scalar, T>,
{
    fn build_absorb_first(
        &self,
        builder: &mut CircuitBuilder,
        chunk: &[Option<Wire>],
    ) -> [Option<Wire>; T] {
        let mut state: [Option<Wire>; T] = [None; T];
        for i in 0..chunk.len() {
            state[i] = chunk[i];
        }
        for i in chunk.len()..T {
            state[i] = Some(builder.add_const_gate(Scalar::ZERO));
        }
        state
    }

    fn witness_absorb_first(
        &self,
        witness: &mut Witness,
        chunk: &[WireOrUnconstrained],
    ) -> [WireOrUnconstrained; T] {
        let mut state = [WireOrUnconstrained::Unconstrained(Scalar::ZERO); T];
        for i in 0..chunk.len() {
            state[i] = chunk[i];
        }
        for i in chunk.len()..T {
            state[i] = witness.assert_constant(Scalar::ZERO).into();
        }
        state
    }

    fn build_absorb(
        &self,
        builder: &mut CircuitBuilder,
        mut state: [Wire; T],
        chunk: &[Option<Wire>],
    ) -> [Wire; T] {
        for i in 0..chunk.len() {
            state[i] = builder.add_sum_gate(state[i].into(), chunk[i]);
        }
        state
    }

    fn witness_absorb(
        &self,
        witness: &mut Witness,
        mut state: [Wire; T],
        chunk: &[WireOrUnconstrained],
    ) -> [Wire; T] {
        for i in 0..chunk.len() {
            state[i] = witness.add(state[i].into(), chunk[i]);
        }
        state
    }

    fn build_sbox(&self, builder: &mut CircuitBuilder, wire: Wire) -> Wire {
        let out = builder.add_square_gate(wire.into());
        let out = builder.add_square_gate(out.into());
        builder.add_mul_gate(out.into(), wire.into())
    }

    fn witness_sbox(&self, witness: &mut Witness, wire: Wire) -> Wire {
        let out = witness.square(wire.into());
        let out = witness.square(out.into());
        witness.mul(out.into(), wire.into())
    }

    fn build_external_linear_t3(
        &self,
        builder: &mut CircuitBuilder,
        state: [Option<Wire>; T],
    ) -> [Wire; T] {
        let sum = builder.add_sum_gate(state[0], state[1]);
        let sum = builder.add_sum_gate(sum.into(), state[2]);
        std::array::from_fn(|i| builder.add_sum_gate(state[i], sum.into()))
    }

    fn witness_external_linear_t3(
        &self,
        witness: &mut Witness,
        state: [WireOrUnconstrained; T],
    ) -> [Wire; T] {
        let sum = witness.add(state[0], state[1]);
        let sum = witness.add(sum.into(), state[2]);
        std::array::from_fn(|i| witness.add(state[i], sum.into()))
    }

    fn build_external_linear_t4(
        &self,
        builder: &mut CircuitBuilder,
        state: [Option<Wire>; T],
    ) -> [Wire; T] {
        let m = BlueSkyConfig4::get_external_matrix();
        std::array::from_fn(|i| {
            let lhs = builder.add_linear_combination_gate(
                m[i * T + 0],
                state[0].into(),
                m[i * T + 1],
                state[1].into(),
            );
            let rhs = builder.add_linear_combination_gate(
                m[i * T + 2],
                state[2].into(),
                m[i * T + 3],
                state[3].into(),
            );
            builder.add_sum_gate(lhs.into(), rhs.into())
        })
    }

    fn witness_external_linear_t4(
        &self,
        witness: &mut Witness,
        state: [WireOrUnconstrained; T],
    ) -> [Wire; T] {
        let m = BlueSkyConfig4::get_external_matrix();
        std::array::from_fn(|i| {
            let lhs = witness.combine(m[i * T + 0], state[0].into(), m[i * T + 1], state[1].into());
            let rhs = witness.combine(m[i * T + 2], state[2].into(), m[i * T + 3], state[3].into());
            witness.add(lhs.into(), rhs.into())
        })
    }

    fn build_external_linear(
        &self,
        builder: &mut CircuitBuilder,
        state: [Option<Wire>; T],
    ) -> [Wire; T] {
        match T {
            3 => self.build_external_linear_t3(builder, state),
            4 => self.build_external_linear_t4(builder, state),
            _ => unimplemented!(),
        }
    }

    fn witness_external_linear(
        &self,
        witness: &mut Witness,
        state: [WireOrUnconstrained; T],
    ) -> [Wire; T] {
        match T {
            3 => self.witness_external_linear_t3(witness, state),
            4 => self.witness_external_linear_t4(witness, state),
            _ => unimplemented!(),
        }
    }

    fn build_internal_linear_t3(
        &self,
        builder: &mut CircuitBuilder,
        mut state: [Wire; T],
    ) -> [Wire; T] {
        let sum = builder.add_sum_gate(state[0].into(), state[1].into());
        let sum = builder.add_sum_gate(sum.into(), state[2].into());
        state[0] = builder.add_sum_gate(state[0].into(), sum.into());
        state[1] = builder.add_sum_gate(state[1].into(), sum.into());
        state[2] = builder.add_linear_combination_gate(
            Scalar::from_const(2),
            state[2].into(),
            Scalar::from_const(1),
            sum.into(),
        );
        state
    }

    fn witness_internal_linear_t3(&self, witness: &mut Witness, mut state: [Wire; T]) -> [Wire; T] {
        let sum = witness.add(state[0].into(), state[1].into());
        let sum = witness.add(sum.into(), state[2].into());
        state[0] = witness.add(state[0].into(), sum.into());
        state[1] = witness.add(state[1].into(), sum.into());
        state[2] = witness.combine(
            Scalar::from_const(2),
            state[2].into(),
            Scalar::from_const(1),
            sum.into(),
        );
        state
    }

    fn build_internal_linear_t4(
        &self,
        builder: &mut CircuitBuilder,
        state: [Wire; T],
    ) -> [Wire; T] {
        let lhs = builder.add_sum_gate(state[0].into(), state[1].into());
        let rhs = builder.add_sum_gate(state[2].into(), state[3].into());
        let sum = builder.add_sum_gate(lhs.into(), rhs.into());
        let m = BlueSkyConfig4::get_internal_matrix();
        std::array::from_fn(|i| {
            builder.add_linear_combination_gate(
                m[i * 5] - Scalar::from_const(1),
                state[i].into(),
                Scalar::from_const(1),
                sum.into(),
            )
        })
    }

    fn witness_internal_linear_t4(&self, witness: &mut Witness, state: [Wire; T]) -> [Wire; T] {
        let lhs = witness.add(state[0].into(), state[1].into());
        let rhs = witness.add(state[2].into(), state[3].into());
        let sum = witness.add(lhs.into(), rhs.into());
        let m = BlueSkyConfig4::get_internal_matrix();
        std::array::from_fn(|i| {
            witness.combine(
                m[i * 5] - Scalar::from_const(1),
                state[i].into(),
                Scalar::from_const(1),
                sum.into(),
            )
        })
    }

    fn build_internal_linear(&self, builder: &mut CircuitBuilder, state: [Wire; T]) -> [Wire; T] {
        match T {
            3 => self.build_internal_linear_t3(builder, state),
            4 => self.build_internal_linear_t4(builder, state),
            _ => unimplemented!(),
        }
    }

    fn witness_internal_linear(&self, witness: &mut Witness, state: [Wire; T]) -> [Wire; T] {
        match T {
            3 => self.witness_internal_linear_t3(witness, state),
            4 => self.witness_internal_linear_t4(witness, state),
            _ => unimplemented!(),
        }
    }

    fn build_full_round(
        &self,
        builder: &mut CircuitBuilder,
        state: [Wire; T],
        r: usize,
    ) -> [Wire; T] {
        let c = BlueSkyConfig::<T>::get_round_constants();
        let mut state: [Wire; T] =
            std::array::from_fn(|i| builder.add_sum_with_const_gate(Some(state[i]), c[r * T + i]));
        for i in 0..T {
            state[i] = self.build_sbox(builder, state[i]);
        }
        self.build_external_linear(builder, state.map(|state| state.into()))
    }

    fn witness_full_round(&self, witness: &mut Witness, state: [Wire; T], r: usize) -> [Wire; T] {
        let c = BlueSkyConfig::<T>::get_round_constants();
        let mut state: [Wire; T] = std::array::from_fn(|i| {
            witness.add_const(WireOrUnconstrained::Wire(state[i]), c[r * T + i].into())
        });
        for i in 0..T {
            state[i] = self.witness_sbox(witness, state[i]);
        }
        self.witness_external_linear(witness, state.map(|state| state.into()))
    }

    fn build_partial_round(
        &self,
        builder: &mut CircuitBuilder,
        mut state: [Wire; T],
        r: usize,
    ) -> [Wire; T] {
        let c = BlueSkyConfig::<T>::get_round_constants();
        state[0] = builder.add_sum_with_const_gate(Some(state[0]), c[r * T]);
        state[0] = self.build_sbox(builder, state[0]);
        self.build_internal_linear(builder, state)
    }

    fn witness_partial_round(
        &self,
        witness: &mut Witness,
        mut state: [Wire; T],
        r: usize,
    ) -> [Wire; T] {
        let c = BlueSkyConfig::<T>::get_round_constants();
        state[0] = witness.add_const(WireOrUnconstrained::Wire(state[0]), c[r * T].into());
        state[0] = self.witness_sbox(witness, state[0]);
        self.witness_internal_linear(witness, state)
    }

    fn build_permutation(
        &self,
        builder: &mut CircuitBuilder,
        state: [Option<Wire>; T],
    ) -> [Wire; T] {
        let num_full_rounds = BlueSkyConfig::<T>::num_full_rounds();
        let num_partial_rounds = BlueSkyConfig::<T>::num_partial_rounds();
        let mut state = self.build_external_linear(builder, state);
        for i in 0..num_full_rounds {
            state = self.build_full_round(builder, state, i);
        }
        for i in 0..num_partial_rounds {
            state = self.build_partial_round(builder, state, num_full_rounds + i);
        }
        for i in 0..num_full_rounds {
            state = self.build_full_round(builder, state, num_full_rounds + num_partial_rounds + i);
        }
        state
    }

    fn witness_permutation(
        &self,
        witness: &mut Witness,
        state: [WireOrUnconstrained; T],
    ) -> [Wire; T] {
        let num_full_rounds = BlueSkyConfig::<T>::num_full_rounds();
        let num_partial_rounds = BlueSkyConfig::<T>::num_partial_rounds();
        let mut state = self.witness_external_linear(witness, state);
        for i in 0..num_full_rounds {
            state = self.witness_full_round(witness, state, i);
        }
        for i in 0..num_partial_rounds {
            state = self.witness_partial_round(witness, state, num_full_rounds + i);
        }
        for i in 0..num_full_rounds {
            state =
                self.witness_full_round(witness, state, num_full_rounds + num_partial_rounds + i);
        }
        state
    }
}

impl<const T: usize, const I: usize> PlonkChip<I, T> for Chip<T, I>
where
    BlueSkyConfig<T>: Config<Scalar, T>,
{
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; I],
    ) -> Result<[Option<Wire>; T]> {
        let mut chunks = inputs.chunks(T - 1);
        let state = self.build_absorb_first(
            builder,
            match chunks.next() {
                Some(chunk) => chunk,
                None => return Err(anyhow!("at least one input scalar is required")),
            },
        );
        let mut state = self.build_permutation(builder, state);
        while let Some(chunk) = chunks.next() {
            state = self.build_absorb(builder, state, chunk);
            state = self.build_permutation(builder, state.map(|wire| Some(wire)));
        }
        Ok(state.map(|state| Some(state)))
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; I],
    ) -> Result<[WireOrUnconstrained; T]> {
        let mut chunks = inputs.chunks(T - 1);
        let state = self.witness_absorb_first(
            witness,
            match chunks.next() {
                Some(chunk) => chunk,
                None => return Err(anyhow!("at least one input scalar is required")),
            },
        );
        let mut state = self.witness_permutation(witness, state);
        while let Some(chunk) = chunks.next() {
            state = self.witness_absorb(witness, state, chunk);
            state = self
                .witness_permutation(witness, state.map(|wire| WireOrUnconstrained::Wire(wire)));
        }
        Ok(state.map(|state| WireOrUnconstrained::Wire(state)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_pcs::hash::{Hash, Poseidon2Hash, Sha2Hash};
    use starkom_plonk::NUM_BLINDING_ROWS;
    use starkom_poseidon2::hash;

    const fn from_const(value: u64) -> Scalar {
        Scalar::from_const(value)
    }

    fn test_hash_chip_impl<H: Hash<Scalar>, const T: usize, const I: usize>(
        blowup_log2: usize,
        inputs: [Scalar; I],
        expected_circuit_size: usize,
    ) where
        BlueSkyConfig<T>: Config<Scalar, T>,
    {
        let outputs = hash::<BlueSkyConfig<T>, Scalar, T>(&inputs);
        let mut builder = CircuitBuilder::default();
        let chip = Chip::<T, I>::default();
        let input_wires = inputs.map(|input| builder.add_const_gate(input));
        let output_wires = chip
            .build(&mut builder, input_wires.map(|wire| Some(wire)))
            .unwrap()
            .map(|wire| wire.unwrap());
        builder.declare_public_gates(
            input_wires
                .map(|wire| wire.gate())
                .into_iter()
                .chain(output_wires.map(|wire| wire.gate())),
        );
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        for i in 0..I {
            witness.assert_constant(inputs[i]);
        }
        assert_eq!(
            chip.witness(&mut witness, input_wires.map(|wire| wire.into()))
                .unwrap(),
            output_wires.map(|wire| WireOrUnconstrained::Wire(wire))
        );
        assert_eq!(output_wires.map(|wire| witness.get(wire)), outputs);
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        assert_eq!(circuit.size(), expected_circuit_size);
        let proof = circuit.prove::<H>(witness, blowup_log2).unwrap();
        let compressed_circuit = circuit.to_compressed::<H>(blowup_log2);
        let openings = compressed_circuit.verify(&proof).unwrap();
        assert_eq!(openings.len(), (inputs.len() + T) * 3);
        for (i, wire) in input_wires.iter().enumerate() {
            assert_eq!(openings[wire], inputs[i]);
        }
        assert_eq!(output_wires.map(|wire| openings[&wire]), outputs);
    }

    fn test_hash_chip<const T: usize, const I: usize>(
        inputs: [Scalar; I],
        expected_circuit_size: usize,
    ) where
        BlueSkyConfig<T>: Config<Scalar, T>,
    {
        test_hash_chip_impl::<Sha2Hash<Scalar>, T, I>(1, inputs, expected_circuit_size);
        test_hash_chip_impl::<Poseidon2Hash<Scalar>, T, I>(1, inputs, expected_circuit_size);
        test_hash_chip_impl::<Sha2Hash<Scalar>, T, I>(2, inputs, expected_circuit_size);
        test_hash_chip_impl::<Poseidon2Hash<Scalar>, T, I>(2, inputs, expected_circuit_size);
        test_hash_chip_impl::<Sha2Hash<Scalar>, T, I>(3, inputs, expected_circuit_size);
        test_hash_chip_impl::<Poseidon2Hash<Scalar>, T, I>(3, inputs, expected_circuit_size);
    }

    #[test]
    fn test_hash_chip_t3_1() {
        test_hash_chip::<3, 1>([from_const(42)], 651);
    }

    #[test]
    fn test_hash_chip_t3_2() {
        test_hash_chip::<3, 2>([from_const(1), from_const(2)], 651);
    }

    #[test]
    fn test_hash_chip_t3_3() {
        test_hash_chip::<3, 3>([from_const(3), from_const(4), from_const(5)], 1298);
    }

    #[test]
    fn test_hash_chip_t3_4() {
        test_hash_chip::<3, 4>(
            [from_const(6), from_const(7), from_const(8), from_const(9)],
            1300,
        );
    }

    #[test]
    fn test_hash_chip_t3_5() {
        test_hash_chip::<3, 5>(
            [
                from_const(10),
                from_const(11),
                from_const(12),
                from_const(13),
                from_const(14),
            ],
            1947,
        );
    }

    #[test]
    fn test_hash_chip_t4_1() {
        test_hash_chip::<4, 1>([from_const(42)], 859);
    }

    #[test]
    fn test_hash_chip_t4_2() {
        test_hash_chip::<4, 2>([from_const(1), from_const(2)], 859);
    }

    #[test]
    fn test_hash_chip_t4_3() {
        test_hash_chip::<4, 3>([from_const(3), from_const(4), from_const(5)], 859);
    }

    #[test]
    fn test_hash_chip_t4_4() {
        test_hash_chip::<4, 4>(
            [from_const(6), from_const(7), from_const(8), from_const(9)],
            1713,
        );
    }

    #[test]
    fn test_hash_chip_t4_5() {
        test_hash_chip::<4, 5>(
            [
                from_const(10),
                from_const(11),
                from_const(12),
                from_const(13),
                from_const(14),
            ],
            1715,
        );
    }

    const BLOWUP_LOG2: usize = 2;

    fn test_preimage_chip<const T: usize, const I: usize>(
        inputs: [Scalar; I],
        expected_circuit_size: usize,
    ) where
        BlueSkyConfig<T>: Config<Scalar, T>,
    {
        let outputs = hash::<BlueSkyConfig<T>, Scalar, T>(&inputs);
        let mut builder = CircuitBuilder::default();
        let chip = Chip::<T, I>::default();
        let output_wires = chip
            .build(&mut builder, std::array::from_fn(|_| None))
            .unwrap()
            .map(|wire| wire.unwrap());
        builder.declare_public_gates(output_wires.map(|wire| wire.gate()));
        let mut witness = Witness::new(builder.len() + NUM_BLINDING_ROWS);
        assert_eq!(
            chip.witness(
                &mut witness,
                inputs.map(|input| WireOrUnconstrained::Unconstrained(input))
            )
            .unwrap(),
            output_wires.map(|wire| WireOrUnconstrained::Wire(wire))
        );
        assert_eq!(output_wires.map(|wire| witness.get(wire)), outputs);
        assert!(builder.check_witness(&witness).is_ok());
        let circuit = builder.build();
        assert_eq!(circuit.size(), expected_circuit_size);
        let proof = circuit
            .prove::<Sha2Hash<Scalar>>(witness, BLOWUP_LOG2)
            .unwrap();
        let compressed_circuit = circuit.to_compressed::<Sha2Hash<Scalar>>(BLOWUP_LOG2);
        let openings = compressed_circuit.verify(&proof).unwrap();
        assert_eq!(openings.len(), T * 3);
        assert_eq!(output_wires.map(|wire| openings[&wire]), outputs);
    }

    #[test]
    fn test_preimage_chip_t3_1() {
        test_preimage_chip::<3, 1>([from_const(42)], 650);
    }

    #[test]
    fn test_preimage_chip_t3_2() {
        test_preimage_chip::<3, 2>([from_const(1), from_const(2)], 649);
    }

    #[test]
    fn test_preimage_chip_t3_3() {
        test_preimage_chip::<3, 3>([from_const(3), from_const(4), from_const(5)], 1295);
    }

    #[test]
    fn test_preimage_chip_t3_4() {
        test_preimage_chip::<3, 4>(
            [from_const(6), from_const(7), from_const(8), from_const(9)],
            1296,
        );
    }

    #[test]
    fn test_preimage_chip_t3_5() {
        test_preimage_chip::<3, 5>(
            [
                from_const(10),
                from_const(11),
                from_const(12),
                from_const(13),
                from_const(14),
            ],
            1942,
        );
    }

    #[test]
    fn test_preimage_chip_t4_1() {
        test_preimage_chip::<4, 1>([from_const(42)], 858);
    }

    #[test]
    fn test_preimage_chip_t4_2() {
        test_preimage_chip::<4, 2>([from_const(1), from_const(2)], 857);
    }

    #[test]
    fn test_preimage_chip_t4_3() {
        test_preimage_chip::<4, 3>([from_const(3), from_const(4), from_const(5)], 856);
    }

    #[test]
    fn test_preimage_chip_t4_4() {
        test_preimage_chip::<4, 4>(
            [from_const(6), from_const(7), from_const(8), from_const(9)],
            1709,
        );
    }

    #[test]
    fn test_preimage_chip_t4_5() {
        test_preimage_chip::<4, 5>(
            [
                from_const(10),
                from_const(11),
                from_const(12),
                from_const(13),
                from_const(14),
            ],
            1710,
        );
    }
}
