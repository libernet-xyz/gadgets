use anyhow::Result;
use starkom_ff::{Field256, PrimeField};
use starkom_plonk::{Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, WitnessView, var};
use std::marker::PhantomData;

/// Generic Poseidon sponge chip.
///
/// Hashes a fixed-length message with the permutation chip `P` (e.g.
/// [`poseidon1::PermutationChip`](`crate::poseidon1::PermutationChip`) or
/// [`poseidon2::PermutationChip`](`crate::poseidon2::PermutationChip`)) using exactly the same
/// sponge construction as [`starkom_poseidon::hash`] and [`starkom_poseidon2::hash`]: `T` is the
/// state size, `R` is the absorption rate, and `C` is the capacity. `T` must be equal to `R+C`.
///
/// The chip has `I` inputs: the first `C` are domain separator tags and the remaining `I-C` are the
/// message, which must not be empty. The DSTs are used to initialize the capacity elements of the
/// state, while the rate elements start at zero. The message is absorbed in blocks of `R` elements
/// that are added to the rate elements, each block followed by a permutation. The `R` outputs of
/// the chip are the rate elements of the final state.
///
/// WARNING: the chip doesn't constrain the DST inputs in any way. All domain separator tags must be
/// fixed and predetermined by the protocol, so it's up to the caller to constrain them accordingly
/// (e.g. by connecting them to cells constrained to constants). If they're left unconstrained a
/// malicious prover can choose the whole initial state, which voids all security properties of the
/// hash.
///
/// WARNING: like the `hash` functions, this chip implicitly pads with zeros up to the next rate
/// boundary, so for example if the rate is 4 the message [1, 2, 3] trivially collides with
/// [1, 2, 3, 0]. The message length is fixed for a given chip type, so this is only a concern if
/// the protocol allows comparing hashes of messages of different lengths. In that case you need to
/// include the length in the message, e.g. by prepending it.
///
/// The chip has one block for every permutation, each block being as tall as the permutation chip.
/// The first row of every block has the following layout:
///
/// ```text
/// 0        R        2R       3R
/// +--------+--------+--------+-----------------------+
/// |   S    |   P    |   M    |   permutation chip    |
/// +--------+--------+--------+-----------------------+
/// ```
///
/// S = rate elements of the state after absorbing M, i.e. the input of the permutation
/// P = rate elements of the state before absorbing M, i.e. the output of the previous permutation
/// M = message block being absorbed
///
/// Each group has R elements and S = P + M, which is the `state[i] += value` step of the `hash`
/// functions. P is constrained to zero in the first block and connected to the previous output in
/// all other blocks. M is connected to the message, except for the elements past the end of the
/// message in the last block, which are constrained to zero (that's the zero padding). The capacity
/// elements of the permutation input are directly connected to the DSTs in the first block and to
/// the capacity elements of the previous output in all other blocks.
///
/// The layout is the same in every block, so if the message fits in a single block the P and M
/// columns are just a waste and it's cheaper to use the permutation chip directly.
#[derive(Debug, Clone, Copy)]
pub struct SpongeChip<
    F: PrimeField,
    P: PlonkChip<F, T, T>,
    const T: usize,
    const R: usize,
    const C: usize,
    const I: usize,
> {
    permutation: P,
    _data: PhantomData<F>,
}

impl<
    F: PrimeField,
    P: PlonkChip<F, T, T>,
    const T: usize,
    const R: usize,
    const C: usize,
    const I: usize,
> SpongeChip<F, P, T, R, C, I>
{
    /// Number of message elements, i.e. inputs that are not DSTs.
    const MESSAGE_LENGTH: usize = I - C;

    /// Number of blocks, which is the same as the number of permutations.
    const NUM_BLOCKS: usize = Self::MESSAGE_LENGTH.div_ceil(R);

    /// Creates a sponge based on the provided permutation chip.
    pub fn new(permutation: P) -> Self {
        const { assert!(T == R + C) };
        const { assert!(R > 0) };
        const { assert!(I > C) };
        Self {
            permutation,
            _data: PhantomData,
        }
    }
}

impl<
    F: PrimeField,
    P: PlonkChip<F, T, T> + Default,
    const T: usize,
    const R: usize,
    const C: usize,
    const I: usize,
> Default for SpongeChip<F, P, T, R, C, I>
{
    fn default() -> Self {
        Self::new(P::default())
    }
}

impl<
    F: PrimeField,
    P: PlonkChip<F, T, T>,
    const T: usize,
    const R: usize,
    const C: usize,
    const I: usize,
> PlonkChip<F, I, R> for SpongeChip<F, P, T, R, C, I>
{
    fn width(&self) -> usize {
        3 * R + self.permutation.width()
    }

    fn height(&self) -> usize {
        Self::NUM_BLOCKS * self.permutation.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; I],
    ) -> Result<[Option<Cell>; R]> {
        let width = self.width();
        let block_height = self.permutation.height();
        let mut state: [Option<Cell>; T] =
            std::array::from_fn(|i| if i < R { None } else { inputs[i - R] });
        for block in 0..Self::NUM_BLOCKS {
            let mut view = view.sub(block * block_height, 0, Some(width), Some(block_height));
            for i in 0..R {
                if block != 0 {
                    view.connect(state[i], view.cell(0, R + i).into());
                } else {
                    view.add_gate(0, var(R + i));
                }
                match inputs.get(C + block * R + i) {
                    Some(&message) => view.connect(message, view.cell(0, 2 * R + i).into()),
                    None => view.add_gate(0, var(2 * R + i)),
                }
                view.add_gate(0, var(R + i) + var(2 * R + i) - var(i));
                state[i] = view.cell(0, i).into();
            }
            state = view.sub_chip(0, 3 * R, &self.permutation, state)?;
        }
        Ok(std::array::from_fn(|i| state[i]))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; I],
    ) -> Result<[CellOrUnconstrained<F>; R]> {
        let width = self.width();
        let block_height = self.permutation.height();
        let mut state: [CellOrUnconstrained<F>; T] =
            std::array::from_fn(|i| if i < R { F::ZERO.into() } else { inputs[i - R] });
        for block in 0..Self::NUM_BLOCKS {
            let mut view = view.sub(block * block_height, 0, Some(width), Some(block_height));
            for i in 0..R {
                let previous = if block != 0 {
                    view.copy(state[i], view.cell(0, R + i))
                } else {
                    view.set(view.cell(0, R + i), F::ZERO);
                    F::ZERO
                };
                let message = match inputs.get(C + block * R + i) {
                    Some(&message) => view.copy(message, view.cell(0, 2 * R + i)),
                    None => {
                        view.set(view.cell(0, 2 * R + i), F::ZERO);
                        F::ZERO
                    }
                };
                view.set(view.cell(0, i), previous + message);
                state[i] = view.cell(0, i).into();
            }
            state = view.sub_chip(0, 3 * R, &self.permutation, state)?;
        }
        Ok(std::array::from_fn(|i| state[i]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poseidon1;
    use crate::poseidon2;
    use primitive_types::H256;
    use starkom_ff::Field;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{Circuit, CircuitBuilder, CircuitOptions, Witness};
    use std::fmt::Debug;
    use std::str::FromStr;

    #[cfg(feature = "bluesky")]
    use starkom_bluesky::Scalar as BS;

    #[cfg(feature = "goldilocks")]
    use starkom_goldilocks::{GL, GL4, from_const as gl};

    #[cfg(feature = "koalabear")]
    use starkom_koalabear::{KB, KB8, from_const as kb};

    #[cfg(feature = "schraderbrau")]
    use starkom_schraderbrau::Scalar as SB;

    #[cfg(feature = "bluesky")]
    type Poseidon1BlueSkyT3 =
        poseidon1::PermutationChip<BS, starkom_poseidon::BlueSkyConfigT3X5, 3>;

    #[cfg(feature = "bluesky")]
    type Poseidon1BlueSkyT4 =
        poseidon1::PermutationChip<BS, starkom_poseidon::BlueSkyConfigT4X5, 4>;

    #[cfg(feature = "schraderbrau")]
    type Poseidon1SchraderbrauT3 =
        poseidon1::PermutationChip<SB, starkom_poseidon::SchraderbrauConfig3, 3>;

    #[cfg(feature = "goldilocks")]
    type Poseidon1GoldilocksT12 =
        poseidon1::PermutationChip<GL, starkom_poseidon::GoldilocksConfig12, 12>;

    #[cfg(feature = "koalabear")]
    type Poseidon1KoalaBearT24 =
        poseidon1::PermutationChip<KB, starkom_poseidon::KoalaBearConfig24, 24>;

    #[cfg(feature = "bluesky")]
    type Poseidon2BlueSkyT3 = poseidon2::PermutationChip<BS, starkom_poseidon2::BlueSkyConfig3, 3>;

    #[cfg(feature = "bluesky")]
    type Poseidon2BlueSkyT4 = poseidon2::PermutationChip<BS, starkom_poseidon2::BlueSkyConfig4, 4>;

    #[cfg(feature = "schraderbrau")]
    type Poseidon2SchraderbrauT4 =
        poseidon2::PermutationChip<SB, starkom_poseidon2::SchraderbrauConfig4, 4>;

    #[cfg(feature = "goldilocks")]
    type Poseidon2GoldilocksT12 =
        poseidon2::PermutationChip<GL, starkom_poseidon2::GoldilocksConfig12, 12>;

    const SECURITY_BITS: usize = 128;
    const BLOWUP_LOG2: usize = 3;

    fn parse<V: FromStr<Err: Debug>>(s: &'static str) -> V {
        s.parse().unwrap()
    }

    /// Returns the DSTs `100, 101, ...` followed by the message `1, 2, 3, ...`.
    ///
    /// The expected outputs used in the tests below have been calculated off-circuit by running
    /// [`starkom_poseidon::hash`] and [`starkom_poseidon2::hash`] on these inputs.
    fn make_inputs<F: PrimeField, const C: usize, const I: usize>() -> [F; I] {
        std::array::from_fn(|i| {
            if i < C {
                F::from(100 + i as u16)
            } else {
                F::from((i - C + 1) as u16)
            }
        })
    }

    /// A circuit containing a [`SpongeChip`] and a witness for it.
    ///
    /// The inputs of the sponge are in row 0 and the sponge itself starts at row 1. Both the input
    /// and output cells of the sponge are public. The witness is generated honestly from the
    /// inputs returned by [`make_inputs`].
    struct SpongeFixture<
        F: PrimeField,
        G: Field256<BaseField = F>,
        P: PlonkChip<F, T, T>,
        const T: usize,
        const R: usize,
        const C: usize,
        const I: usize,
    > {
        chip: SpongeChip<F, P, T, R, C, I>,
        inputs: [F; I],
        input_cells: [Cell; I],
        output_cells: [Cell; R],
        circuit: Circuit<F, G>,
        witness: Witness<F>,
        witness_outputs: [CellOrUnconstrained<F>; R],
    }

    impl<
        F: PrimeField,
        G: Field256<BaseField = F>,
        P: PlonkChip<F, T, T> + Default,
        const T: usize,
        const R: usize,
        const C: usize,
        const I: usize,
    > SpongeFixture<F, G, P, T, R, C, I>
    {
        fn new() -> Self {
            let chip = SpongeChip::<F, P, T, R, C, I>::default();
            let inputs = make_inputs::<F, C, I>();
            let mut builder = CircuitBuilder::<F, G>::default();
            let input_cells: [Cell; I] = std::array::from_fn(|i| builder.cell(0, i));
            let output_cells = builder
                .sub_chip(1, 0, &chip, input_cells.map(Some))
                .unwrap()
                .map(Option::unwrap);
            builder.declare_public_cells(input_cells.into_iter().chain(output_cells));
            let circuit = builder
                .build(CircuitOptions {
                    canonicalize_constraints: false,
                    security_bits: SECURITY_BITS,
                    blowup_log2: BLOWUP_LOG2,
                    grinding_bits: 0,
                    blind: false,
                })
                .unwrap();
            let mut witness = circuit.make_witness();
            for i in 0..I {
                witness.set(input_cells[i], inputs[i]);
            }
            let witness_outputs = witness
                .sub_chip(1, 0, &chip, input_cells.map(CellOrUnconstrained::Cell))
                .unwrap();
            Self {
                chip,
                inputs,
                input_cells,
                output_cells,
                circuit,
                witness,
                witness_outputs,
            }
        }

        /// Simulates a malicious prover that adds `delta` to the cells at the specified `columns`
        /// in the first row of the last block, then recomputes the last permutation from the S
        /// cells so that all wire constraints keep holding.
        ///
        /// NOTE: this assumes that the permutation chip copies its i-th input to its cell (0, i),
        /// which is the case for both the Poseidon1 and the Poseidon2 permutation chips.
        fn forge_last_block(&mut self, columns: &[usize], delta: F) {
            let row = 1 + self.chip.height() - self.chip.permutation.height();
            for &column in columns {
                let cell = self.witness.cell(row, column);
                self.witness.set(cell, self.witness.get_at(cell) + delta);
            }
            let inputs = std::array::from_fn(|i| {
                if i < R {
                    self.witness.cell(row, i).into()
                } else {
                    self.witness.cell(row, 3 * R + i).into()
                }
            });
            self.witness
                .sub_chip(row, 3 * R, &self.chip.permutation, inputs)
                .unwrap();
        }
    }

    fn test_sponge_impl<
        F: PrimeField,
        G: Field256<BaseField = F>,
        P: PlonkChip<F, T, T> + Default,
        const T: usize,
        const R: usize,
        const C: usize,
        const I: usize,
    >(
        expected_width: usize,
        expected_height: usize,
        expected_output: [F; R],
        circuit_commitment: H256,
    ) {
        let SpongeFixture {
            chip,
            inputs,
            input_cells,
            output_cells,
            circuit,
            witness,
            witness_outputs,
        } = SpongeFixture::<F, G, P, T, R, C, I>::new();
        assert_eq!(chip.width(), expected_width);
        assert_eq!(chip.height(), expected_height);
        assert_eq!(circuit.num_rows(), 1 + expected_height);
        assert_eq!(circuit.num_columns(), expected_width);
        assert_eq!(witness_outputs, output_cells.map(CellOrUnconstrained::Cell));
        circuit.check_witness(&witness).unwrap();
        let proof = circuit
            .prove::<Sha2Hash<G>>(witness, &circuit.default_proof_options())
            .unwrap();
        assert!(proof.pcs_config().satisfies(circuit.pcs_config()));
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let public_inputs = circuit.verify(&proof).unwrap();
        assert_eq!(input_cells.map(|cell| public_inputs[&cell]), inputs);
        assert_eq!(
            output_cells.map(|cell| public_inputs[&cell]),
            expected_output
        );
    }

    fn test_tampered_input_impl<
        F: PrimeField,
        G: Field256<BaseField = F>,
        P: PlonkChip<F, T, T> + Default,
        const T: usize,
        const R: usize,
        const C: usize,
        const I: usize,
    >(
        index: usize,
    ) -> Result<()> {
        let mut fixture = SpongeFixture::<F, G, P, T, R, C, I>::new();
        fixture
            .witness
            .set(fixture.input_cells[index], fixture.inputs[index] + F::ONE);
        fixture.circuit.check_witness(&fixture.witness)
    }

    fn test_forged_last_block_impl<
        F: PrimeField,
        G: Field256<BaseField = F>,
        P: PlonkChip<F, T, T> + Default,
        const T: usize,
        const R: usize,
        const C: usize,
        const I: usize,
    >(
        columns: &[usize],
        delta: F,
    ) -> Result<()> {
        let mut fixture = SpongeFixture::<F, G, P, T, R, C, I>::new();
        fixture.forge_last_block(columns, delta);
        fixture.circuit.check_witness(&fixture.witness)
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_sponge_poseidon1_t3_bs() {
        test_sponge_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 2>(
            91,
            1,
            [
                parse("0x4381d1ea2b2e3e264aec0c7a7a6f316f536e8db66b69123a7a7a64884bd636d0"),
                parse("0x5a58da4b384aa6b74c06a18fe40a837b7b8f756e06064180c3c966859f75d049"),
            ],
            parse("0xabc0b16022941f200d5f07c247b0f771d0cda528a9be0cea49f679c4bf03ab66"),
        );
        test_sponge_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 3>(
            91,
            1,
            [
                parse("0x2211e25575a1c1217c237255c41496f2f879011377d70dd985b51651ff0b0900"),
                parse("0x3d4bceb13ed19d3399429f5ef3cc255eaffb6b3109837d7df932ddbc9720a227"),
            ],
            parse("0x49b46c0ab33cd1074e391377683b4920aeacf197c3d213336e78a30c98a7a635"),
        );
        test_sponge_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 4>(
            91,
            2,
            [
                parse("0x02496e885e73e6afdd957d9752873172f59ae0bce4b202f4de13c1caeedf589a"),
                parse("0x1efb18e5114aeb46d0ee3709e07bbe7f1687a463f3eab9ee83342f828742b3f3"),
            ],
            parse("0xd22a0c61f1d32c21be57a26136351d8eb49a8c24ed511b79b53d67259c5c4e39"),
        );
        test_sponge_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(
            91,
            2,
            [
                parse("0x00dcb52229b91f80158e426fb94e4462ff0c0e6a06b9ab677e30ccb05502aa42"),
                parse("0x78a06199c8eb975c59ab0640a1c2905a2d79bb3eae1ceda2c6bda0a6afe03357"),
            ],
            parse("0x2204c286c58e06c5a931526899c3921bddf02d5393dd436d09e328906054f136"),
        );
        test_sponge_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 6>(
            91,
            3,
            [
                parse("0x563dfce2a5d490c45e6d77bf122d8d565db010e05d375a5bebea46553ee74bd8"),
                parse("0x28da1f990bcad848c583da829fe711c332dd43e61bd0d53ab3482637b45d3ae8"),
            ],
            parse("0x3460727108ab0bbd76c31512c4317f4316dbf9e6df37772e912876a0c5973d9f"),
        );
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_sponge_poseidon1_t4_bs() {
        test_sponge_impl::<BS, BS, Poseidon1BlueSkyT4, 4, 2, 2, 5>(
            101,
            2,
            [
                parse("0x61d60c7699b70ef978a9234014f39b65762e7200a738f56b51b0a2298d72ce8a"),
                parse("0x4e4030d8711c3e5cfe5f06585aca8de02b09c29f8954c2377f963276f1687c57"),
            ],
            parse("0x44b4ba3c824b9b6135db4dcf111190f7798d6ff860904028aca59f138fef60ca"),
        );
    }

    #[cfg(feature = "schraderbrau")]
    #[test]
    fn test_sponge_poseidon1_t3_sb() {
        test_sponge_impl::<SB, SB, Poseidon1SchraderbrauT3, 3, 2, 1, 4>(
            118,
            2,
            [
                parse("0x11efa96295661bad96d7b61520d43a8157dfc14916d8cb1545e73cc543885685"),
                parse("0x49e6a9d265a33ef1d2f77e8713a54706f4d59f825c33c34ac2df4d82339462bc"),
            ],
            parse("0xa862d0f7ad8e40fb645ce0433a62582634f8847ca5010f4bead1ba0e97bdaca4"),
        );
    }

    #[cfg(feature = "goldilocks")]
    #[test]
    fn test_sponge_poseidon1_t12_gl() {
        test_sponge_impl::<GL, GL4, Poseidon1GoldilocksT12, 12, 8, 4, 13>(
            165,
            2,
            [
                gl(0xd7fe703ef169a131),
                gl(0x0cbd70ff7fbd65c6),
                gl(0x0c13435ae650c3ad),
                gl(0x0cbb855b321838bb),
                gl(0x2d1e7422c239921d),
                gl(0x2332aac1f114504c),
                gl(0xf200f2eaee2cb60b),
                gl(0xc02e3026731d3837),
            ],
            parse("0x166af0b0103850ec058eafc3e946473097c4e274ba74ad4a8a00433c7ad7f9e4"),
        );
    }

    #[cfg(feature = "koalabear")]
    #[test]
    fn test_sponge_poseidon1_t24_kb() {
        test_sponge_impl::<KB, KB8, Poseidon1KoalaBearT24, 24, 16, 8, 25>(
            310,
            2,
            [
                kb(0x5e332d95),
                kb(0x091b617b),
                kb(0x64b2f9ee),
                kb(0x40ec4a68),
                kb(0x196bccb0),
                kb(0x3eb1a91a),
                kb(0x40c501d0),
                kb(0x4e1b8058),
                kb(0x6e56d387),
                kb(0x40e1ad12),
                kb(0x20745e9b),
                kb(0x467f10c7),
                kb(0x675d6591),
                kb(0x58755628),
                kb(0x7051ce17),
                kb(0x2b5a384d),
            ],
            parse("0xe9d1a20d3c6fafdbfa8d58d5c211931091f298098287bcd7c38f3619508104a3"),
        );
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_sponge_poseidon2_t3_bs() {
        test_sponge_impl::<BS, BS, Poseidon2BlueSkyT3, 3, 2, 1, 4>(
            94,
            2,
            [
                parse("0x12f29346e161ced9e035d89992527158f50addde6bfb49b9f9cb8910ca5af618"),
                parse("0x0d56c56ab54e77773ef8f7f19a99057d5782435bb32a34ca4934598569025701"),
            ],
            parse("0x9f2d2bd953763f78df14e0e119f960dfedc2bc0fd4a83c05166e11381eafca42"),
        );
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_sponge_poseidon2_t4_bs() {
        test_sponge_impl::<BS, BS, Poseidon2BlueSkyT4, 4, 3, 1, 5>(
            108,
            2,
            [
                parse("0x4a0f4b0ff78583837a67f960bcf4107507e2d0334260959a9f05ed046a19ab95"),
                parse("0x175f55b671cc19b43c465849f7a4afc0d06daf28383d7ca9495ab09ab09f3097"),
                parse("0x6fb73953a1eed14e87064341135048a1f5d78c4b89d971adedf1b8c23e9ed615"),
            ],
            parse("0xe2b887eb064998ebec1b6a1c010979182534642478f018eb63d595d8caaec4bf"),
        );
    }

    #[cfg(feature = "schraderbrau")]
    #[test]
    fn test_sponge_poseidon2_t4_sb() {
        test_sponge_impl::<SB, SB, Poseidon2SchraderbrauT4, 4, 3, 1, 4>(
            136,
            1,
            [
                parse("0x29e01c5113c0232788d5a2c0874c074c2aef0f3067a0246867709bd613fd467a"),
                parse("0x0e4c1a30cf610d310b4e5a76fb36106ce5e32db55004afdc79866acec608555a"),
                parse("0x70e8bd2734854bbd91473c7c6901482b9bb70855d1a6527b041d510aefc70783"),
            ],
            parse("0xd41585866a26c1cc4d9cd9eee978e9c2538cd727c98103bfea73b28c562ade9c"),
        );
    }

    #[cfg(feature = "goldilocks")]
    #[test]
    fn test_sponge_poseidon2_t12_gl() {
        test_sponge_impl::<GL, GL4, Poseidon2GoldilocksT12, 12, 8, 4, 5>(
            177,
            1,
            [
                gl(0xb815d2807191c099),
                gl(0x206b23ba36864cac),
                gl(0x2f7589808ab18f0a),
                gl(0xac18bd7a59939bf2),
                gl(0x146461900fe58aed),
                gl(0xc07463a4f3a8e90a),
                gl(0xf475eceec2c0cd6d),
                gl(0xd8aaf6a618e3f326),
            ],
            parse("0x3dcd5da359ac2ee5dd3788f5c24a7ef96e03011d73b5e89672964a86da8a24b6"),
        );
        test_sponge_impl::<GL, GL4, Poseidon2GoldilocksT12, 12, 8, 4, 13>(
            177,
            2,
            [
                gl(0x63412adbb21758bf),
                gl(0x90363f66fbf653b7),
                gl(0x94a891ffc19484d0),
                gl(0x6afcbe10527f96d6),
                gl(0x3ba03caf6fda74ff),
                gl(0x8e39404eca97d880),
                gl(0x04b0432cf916bb60),
                gl(0xcd9748d4dc482e4d),
            ],
            parse("0x77af805400a26e4445ffc0da0be16c6f790faf7da65c574e193a88b28c677479"),
        );
    }

    #[cfg(feature = "bluesky")]
    #[test]
    fn test_tampered_input() {
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(0).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(1).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(2).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(3).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(4).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT4, 4, 2, 2, 5>(0).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT4, 4, 2, 2, 5>(1).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT4, 4, 2, 2, 5>(2).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT4, 4, 2, 2, 5>(3).is_err());
        assert!(test_tampered_input_impl::<BS, BS, Poseidon1BlueSkyT4, 4, 2, 2, 5>(4).is_err());
    }

    /// Changes only S, so that S = P + M doesn't hold.
    #[cfg(feature = "bluesky")]
    #[test]
    fn test_forged_sum() {
        let zero = BS::ZERO;
        let one = BS::ONE;
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 2>(&[1], zero)
                .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 2>(&[1], one)
                .is_err()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 4>(&[0], zero)
                .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 4>(&[0], one)
                .is_err()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(&[1], zero)
                .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 5>(&[1], one)
                .is_err()
        );
    }

    /// Changes S and P consistently in a single-block sponge, so that only the zero constraint on
    /// the initial P can catch it. With R = 2 the S cells are in columns 0 and 1, and the P cells
    /// are in columns 2 and 3.
    #[cfg(feature = "bluesky")]
    #[test]
    fn test_forged_initial_state() {
        let zero = BS::ZERO;
        let one = BS::ONE;
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 3>(&[0, 2], zero)
                .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 3>(&[0, 2], one)
                .is_err()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 3>(&[1, 3], zero)
                .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 3>(&[1, 3], one)
                .is_err()
        );
    }

    /// Changes S and M consistently at a padding position of the last block, so that only the zero
    /// padding constraint can catch it. With R = 2 the S cells are in columns 0 and 1, and the M
    /// cells are in columns 4 and 5.
    #[cfg(feature = "bluesky")]
    #[test]
    fn test_forged_padding() {
        let zero = BS::ZERO;
        let one = BS::ONE;
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 2>(&[1, 5], zero)
                .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 2>(&[1, 5], one)
                .is_err()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 4>(&[1, 5], zero)
                .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<BS, BS, Poseidon1BlueSkyT3, 3, 2, 1, 4>(&[1, 5], one)
                .is_err()
        );
    }

    /// Same as `test_forged_padding` but with R = 8, so the last S cell is in column 7 and the last
    /// M cell is in column 23.
    #[cfg(feature = "goldilocks")]
    #[test]
    fn test_forged_padding_gl() {
        let zero = GL::ZERO;
        let one = GL::ONE;
        assert!(
            test_forged_last_block_impl::<GL, GL4, Poseidon2GoldilocksT12, 12, 8, 4, 5>(
                &[7, 23],
                zero
            )
            .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<GL, GL4, Poseidon2GoldilocksT12, 12, 8, 4, 5>(
                &[7, 23],
                one
            )
            .is_err()
        );
        assert!(
            test_forged_last_block_impl::<GL, GL4, Poseidon2GoldilocksT12, 12, 8, 4, 13>(
                &[7, 23],
                zero
            )
            .is_ok()
        );
        assert!(
            test_forged_last_block_impl::<GL, GL4, Poseidon2GoldilocksT12, 12, 8, 4, 13>(
                &[7, 23],
                one
            )
            .is_err()
        );
    }
}
