use crate::poseidon2;
use crate::xits;
use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::{Field, PrimeField};
use starkom_plonk::{Chip as PlonkChip, CircuitBuilder, Wire, WireOrUnconstrained, Witness};

/// Returns a boolean signal indicating whether an input trit is equal to 0.
///
/// For internal use by ternary SMTs.
///
/// The implementation uses the formula:
///
///   x^2 / 2 + 1 - 3x
///
/// derived from the Lagrange basis that activates on 0 and remains 0 on 1 and 2.
#[derive(Debug, Default, Clone)]
struct TernarySelector0 {}

impl PlonkChip<1, 1> for TernarySelector0 {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; 1]> {
        let [x] = inputs;
        let y = builder.add_unary_gate(
            -Scalar::from_const(3) * Scalar::TWO_INV,
            Scalar::from_const(0),
            -Scalar::from_const(1),
            Scalar::TWO_INV,
            Scalar::from_const(1),
            x,
        );
        Ok([y.into()])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; 1]> {
        let [x] = inputs;
        let gate = witness.pop_gate();
        let lhs = Wire::LeftIn(gate);
        let rhs = Wire::RightIn(gate);
        let out = Wire::Out(gate);
        witness.copy(x, lhs);
        witness.copy(x, rhs);
        let trit = witness.get(lhs).try_to_u8().unwrap() as usize;
        let value: u64 = [1, 0, 0][trit];
        witness.set(out, value.into());
        Ok([out.into()])
    }
}

/// Returns a boolean signal indicating whether an input trit is equal to 1.
///
/// For internal use by ternary SMTs.
///
/// The implementation uses the formula:
///
///   2x - x^2
///
/// derived from the Lagrange basis that activates on 1 and remains 0 on 0 and 2.
#[derive(Debug, Default, Clone)]
struct TernarySelector1 {}

impl PlonkChip<1, 1> for TernarySelector1 {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; 1]> {
        let [x] = inputs;
        let y = builder.add_unary_gate(
            Scalar::from_const(2),
            Scalar::from_const(0),
            -Scalar::from_const(1),
            -Scalar::from_const(1),
            Scalar::from_const(0),
            x,
        );
        Ok([y.into()])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; 1]> {
        let [x] = inputs;
        let gate = witness.pop_gate();
        let lhs = Wire::LeftIn(gate);
        let rhs = Wire::RightIn(gate);
        let out = Wire::Out(gate);
        witness.copy(x, lhs);
        witness.copy(x, rhs);
        let trit = witness.get(lhs).try_to_u8().unwrap() as usize;
        let value: u64 = [0, 1, 0][trit];
        witness.set(out, value.into());
        Ok([out.into()])
    }
}

/// Returns a boolean signal indicating whether an input trit is equal to 2.
///
/// For internal use by ternary SMTs.
///
/// The implementation uses the formula:
///
///   (x^2 - x) / 2
///
/// derived from the Lagrange basis that activates on 2 and remains 0 on 0 and 1.
#[derive(Debug, Default, Clone)]
struct TernarySelector2 {}

impl PlonkChip<1, 1> for TernarySelector2 {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; 1]> {
        let [x] = inputs;
        let y = builder.add_unary_gate(
            -Scalar::TWO_INV,
            Scalar::from_const(0),
            -Scalar::from_const(1),
            Scalar::TWO_INV,
            Scalar::from_const(0),
            x,
        );
        Ok([y.into()])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; 1]> {
        let [x] = inputs;
        let gate = witness.pop_gate();
        let lhs = Wire::LeftIn(gate);
        let rhs = Wire::RightIn(gate);
        let out = Wire::Out(gate);
        witness.copy(x, lhs);
        witness.copy(x, rhs);
        let trit = witness.get(lhs).try_to_u8().unwrap() as usize;
        let value: u64 = [0, 0, 1][trit];
        witness.set(out, value.into());
        Ok([out.into()])
    }
}

#[derive(Debug, Default, Clone)]
struct ThreeWayDemux {
    selector0: TernarySelector0,
    selector1: TernarySelector1,
    selector2: TernarySelector2,
}

impl PlonkChip<1, 3> for ThreeWayDemux {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 1],
    ) -> Result<[Option<Wire>; 3]> {
        let [out0] = self.selector0.build(builder, inputs)?;
        let [out1] = self.selector1.build(builder, inputs)?;
        let [out2] = self.selector2.build(builder, inputs)?;
        Ok([out0, out1, out2])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 1],
    ) -> Result<[WireOrUnconstrained; 3]> {
        let [out0] = self.selector0.witness(witness, inputs)?;
        let [out1] = self.selector1.witness(witness, inputs)?;
        let [out2] = self.selector2.witness(witness, inputs)?;
        Ok([out0, out1, out2])
    }
}

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height `H`.
///
/// WARNING: `H` must be strictly less than 255. Do NOT use this chip if `H` spans the full BlueSky
/// range, as in that case the bit decomposition of the key would be UNSAFE! Use the
/// [`FullBinaryChip`] below instead.
#[derive(Debug, Clone)]
pub struct BinaryChip<const H: usize> {
    decomposer: xits::BitDecomposerChip<H>,
    hasher: poseidon2::Chip<3, 2>,
    path: [[Scalar; 2]; H],
}

impl<const H: usize> BinaryChip<H> {
    pub fn new(path: [[Scalar; 2]; H]) -> Self {
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher: poseidon2::Chip::default(),
            path,
        }
    }
}

impl<const H: usize> Default for BinaryChip<H> {
    fn default() -> Self {
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher: poseidon2::Chip::default(),
            path: [[Scalar::ZERO; 2]; H],
        }
    }
}

impl<const H: usize> PlonkChip<2, 1> for BinaryChip<H> {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 2],
    ) -> Result<[Option<Wire>; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.build(builder, [key])?;
        let mut hash = value;
        for i in 0..H {
            let bit = bits[i];
            let not = builder.add_not_gate(bit).into();
            let lhs = {
                let lhs = builder.add_mul_gate(bit, None);
                let rhs = builder.add_mul_gate(not, hash.into());
                builder.add_sum_gate(lhs.into(), rhs.into())
            };
            let rhs = {
                let lhs = builder.add_mul_gate(bit, hash.into());
                let rhs = builder.add_mul_gate(not, None);
                builder.add_sum_gate(lhs.into(), rhs.into())
            };
            [hash, _, _] = self.hasher.build(builder, [lhs.into(), rhs.into()])?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 2],
    ) -> Result<[WireOrUnconstrained; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.witness(witness, [key])?;
        let mut hash = value;
        for i in 0..H {
            let bit = bits[i];
            let not = witness.not(bit).into();
            let lhs = {
                let lhs = witness.mul(bit, self.path[i][0].into());
                let rhs = witness.mul(not, hash);
                witness.add(lhs.into(), rhs.into())
            };
            let rhs = {
                let lhs = witness.mul(bit, hash);
                let rhs = witness.mul(not, self.path[i][1].into());
                witness.add(lhs.into(), rhs.into())
            };
            [hash, _, _] = self.hasher.witness(witness, [lhs.into(), rhs.into()])?;
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a ternary Sparse Merkle Tree of height `H`.
///
/// WARNING: `H` must be strictly less than 161. Do NOT use this chip if `H` spans the full BlueSky
/// range, as in that case the trit decomposition of the key would be UNSAFE! Use the
/// [`FullTernaryChip`] below instead.
#[derive(Debug, Clone)]
pub struct TernaryChip<const H: usize> {
    decomposer: xits::TritDecomposerChip<H>,
    demux: ThreeWayDemux,
    hasher: poseidon2::Chip<4, 3>,
    path: [[Scalar; 3]; H],
}

impl<const H: usize> TernaryChip<H> {
    pub fn new(path: [[Scalar; 3]; H]) -> Self {
        Self {
            decomposer: xits::TritDecomposerChip::default(),
            demux: ThreeWayDemux::default(),
            hasher: poseidon2::Chip::default(),
            path,
        }
    }
}

impl<const H: usize> Default for TernaryChip<H> {
    fn default() -> Self {
        Self {
            decomposer: xits::TritDecomposerChip::default(),
            demux: ThreeWayDemux::default(),
            hasher: poseidon2::Chip::default(),
            path: [[Scalar::ZERO; 3]; H],
        }
    }
}

impl<const H: usize> PlonkChip<2, 1> for TernaryChip<H> {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 2],
    ) -> Result<[Option<Wire>; 1]> {
        let [key, value] = inputs;
        let trits = self.decomposer.build(builder, [key])?;
        let mut hash = value;
        for i in 0..H {
            let trit = trits[i];
            let [selector0, selector1, selector2] = self.demux.build(builder, [trit])?;
            let inverted0 = builder.add_not_gate(selector0).into();
            let inverted1 = builder.add_not_gate(selector1).into();
            let inverted2 = builder.add_not_gate(selector2).into();
            let input0 = {
                let lhs = builder.add_mul_gate(selector0, hash.into());
                let rhs = builder.add_mul_gate(inverted0, None);
                builder.add_sum_gate(lhs.into(), rhs.into()).into()
            };
            let input1 = {
                let lhs = builder.add_mul_gate(selector1, hash.into());
                let rhs = builder.add_mul_gate(inverted1, None);
                builder.add_sum_gate(lhs.into(), rhs.into()).into()
            };
            let input2 = {
                let lhs = builder.add_mul_gate(selector2, hash.into());
                let rhs = builder.add_mul_gate(inverted2, None);
                builder.add_sum_gate(lhs.into(), rhs.into()).into()
            };
            [hash, _, _, _] = self.hasher.build(builder, [input0, input1, input2])?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 2],
    ) -> Result<[WireOrUnconstrained; 1]> {
        let [key, value] = inputs;
        let trits = self.decomposer.witness(witness, [key])?;
        let mut hash = value;
        for i in 0..H {
            let trit = trits[i];
            let [selector0, selector1, selector2] = self.demux.witness(witness, [trit])?;
            let inverted0 = witness.not(selector0).into();
            let inverted1 = witness.not(selector1).into();
            let inverted2 = witness.not(selector2).into();
            let input0 = {
                let lhs = witness.mul(selector0, hash);
                let rhs = witness.mul(inverted0, self.path[i][0].into());
                witness.add(lhs.into(), rhs.into()).into()
            };
            let input1 = {
                let lhs = witness.mul(selector1, hash);
                let rhs = witness.mul(inverted1, self.path[i][1].into());
                witness.add(lhs.into(), rhs.into()).into()
            };
            let input2 = {
                let lhs = witness.mul(selector2, hash);
                let rhs = witness.mul(inverted2, self.path[i][2].into());
                witness.add(lhs.into(), rhs.into()).into()
            };
            [hash, _, _, _] = self.hasher.witness(witness, [input0, input1, input2])?;
        }
        Ok([hash])
    }
}

/// Runs a Merkle lookup over a binary Sparse Merkle Tree of height 256.
///
/// The keys of such a tree span the full BlueSky range. Internally this chip uses a
/// [`xits::FullBitDecomposerChip`], making the 256-bit decomposition safe at the cost of some extra
/// constraints.
///
/// If you don't need 256- or 255-bit keys use [`BinaryChip`].
#[derive(Debug, Clone)]
pub struct FullBinaryChip {
    decomposer: xits::FullBitDecomposerChip,
    hasher: poseidon2::Chip<3, 2>,
    path: [[Scalar; 2]; 256],
}

impl FullBinaryChip {
    pub fn new(path: [[Scalar; 2]; 256]) -> Self {
        Self {
            decomposer: xits::FullBitDecomposerChip::default(),
            hasher: poseidon2::Chip::default(),
            path,
        }
    }
}

impl Default for FullBinaryChip {
    fn default() -> Self {
        Self {
            decomposer: xits::FullBitDecomposerChip::default(),
            hasher: poseidon2::Chip::default(),
            path: [[Scalar::ZERO; 2]; 256],
        }
    }
}

impl PlonkChip<2, 1> for FullBinaryChip {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; 2],
    ) -> Result<[Option<Wire>; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.build(builder, [key])?;
        let mut hash = value;
        for i in 0..256 {
            let bit = bits[i];
            let not = builder.add_not_gate(bit).into();
            let lhs = {
                let lhs = builder.add_mul_gate(bit, None);
                let rhs = builder.add_mul_gate(not, hash.into());
                builder.add_sum_gate(lhs.into(), rhs.into())
            };
            let rhs = {
                let lhs = builder.add_mul_gate(bit, hash.into());
                let rhs = builder.add_mul_gate(not, None);
                builder.add_sum_gate(lhs.into(), rhs.into())
            };
            [hash, _, _] = self.hasher.build(builder, [lhs.into(), rhs.into()])?;
        }
        Ok([hash])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; 2],
    ) -> Result<[WireOrUnconstrained; 1]> {
        let [key, value] = inputs;
        let bits = self.decomposer.witness(witness, [key])?;
        let mut hash = value;
        for i in 0..256 {
            let bit = bits[i];
            let not = witness.not(bit).into();
            let lhs = {
                let lhs = witness.mul(bit, self.path[i][0].into());
                let rhs = witness.mul(not, hash);
                witness.add(lhs.into(), rhs.into())
            };
            let rhs = {
                let lhs = witness.mul(bit, hash);
                let rhs = witness.mul(not, self.path[i][1].into());
                witness.add(lhs.into(), rhs.into())
            };
            [hash, _, _] = self.hasher.witness(witness, [lhs.into(), rhs.into()])?;
        }
        Ok([hash])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use primitive_types::U256;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_ff::Field256;
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{Circuit, CompressedCircuit};
    use std::fmt::Debug;
    use std::sync::{Arc, LazyLock};

    fn test_binary_smt<const H: usize>(
        key: u64,
        value: u64,
        path: [[Scalar; 2]; H],
        expected_root_hash: Scalar,
    ) -> Result<()> {
        let key = Scalar::from(key);
        let value = Scalar::from(value);
        let mut builder = CircuitBuilder::default();
        let inputs = builder.add_nop_gate(None, None, None);
        let key_wire = Wire::LeftIn(inputs);
        let value_wire = Wire::RightIn(inputs);
        let expected_root_hash_wire = Wire::Out(inputs);
        let chip = BinaryChip::<H>::new(path);
        let root_hash_wire =
            chip.build(&mut builder, [key_wire.into(), value_wire.into()])?[0].unwrap();
        builder.connect(root_hash_wire, expected_root_hash_wire);
        builder.declare_public_gates([inputs]);
        let circuit = builder.build();
        let mut witness = circuit.make_witness();
        witness.nop(key.into(), value.into(), expected_root_hash.into());
        chip.witness(&mut witness, [key.into(), value.into()])?;
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, 2)?;
        let public_inputs = circuit.verify(&proof)?;
        assert_eq!(public_inputs[&key_wire], key);
        assert_eq!(public_inputs[&value_wire], value);
        assert_eq!(public_inputs[&expected_root_hash_wire], expected_root_hash);
        Ok(())
    }

    #[test]
    fn test_binary_smt_height_one_1() {
        let path = [[from_const(12), from_const(34)]];
        let root_hash =
            parse_scalar("0x165e74be18ef4be6de5e232cd3480dcc38176807ac918b904576964612c5b6de");
        assert!(test_binary_smt::<1>(0, 12, path, root_hash).is_ok());
        assert!(test_binary_smt::<1>(1, 34, path, root_hash).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_2() {
        let path = [[from_const(34), from_const(12)]];
        let root_hash =
            parse_scalar("0x5edee3683c6686f34afdcc9bdabd06cc43c2d5a28d6509a6301eceb33e255e72");
        assert!(test_binary_smt::<1>(0, 34, path, root_hash).is_ok());
        assert!(test_binary_smt::<1>(1, 12, path, root_hash).is_ok());
    }

    #[test]
    fn test_binary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78)]];
        let root_hash =
            parse_scalar("0x38e7bb7b6ccae0c74031423877db058f4ab3a284964e2d91bf97497851eca5db");
        assert!(test_binary_smt::<1>(0, 56, path, root_hash).is_ok());
        assert!(test_binary_smt::<1>(1, 78, path, root_hash).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_1() {
        let path = [
            [from_const(12), from_const(34)],
            [
                parse_scalar("0x165e74be18ef4be6de5e232cd3480dcc38176807ac918b904576964612c5b6de"),
                parse_scalar("0x38e7bb7b6ccae0c74031423877db058f4ab3a284964e2d91bf97497851eca5db"),
            ],
        ];
        let root_hash =
            parse_scalar("0x2945d8bf64346bb2085bb6ec93520cd641a3c77f1501dba2cbcb8982eb8dbaa1");
        assert!(test_binary_smt::<2>(0, 12, path, root_hash).is_ok());
        assert!(test_binary_smt::<2>(1, 34, path, root_hash).is_ok());
    }

    #[test]
    fn test_binary_smt_height_two_2() {
        let path = [
            [from_const(56), from_const(78)],
            [
                parse_scalar("0x165e74be18ef4be6de5e232cd3480dcc38176807ac918b904576964612c5b6de"),
                parse_scalar("0x38e7bb7b6ccae0c74031423877db058f4ab3a284964e2d91bf97497851eca5db"),
            ],
        ];
        let root_hash =
            parse_scalar("0x2945d8bf64346bb2085bb6ec93520cd641a3c77f1501dba2cbcb8982eb8dbaa1");
        assert!(test_binary_smt::<2>(2, 56, path, root_hash).is_ok());
        assert!(test_binary_smt::<2>(3, 78, path, root_hash).is_ok());
    }

    fn test_ternary_selector_impl<C: Default + PlonkChip<1, 1>>(
        input: u8,
        output: u8,
    ) -> Result<()> {
        let mut builder = CircuitBuilder::default();
        let chip = C::default();
        let trit = builder.add_const_gate(input.into());
        let out = chip.build(&mut builder, [trit.into()])?[0].unwrap();
        let out = builder.add_nop_gate(None, None, out.into());
        builder.declare_public_gates([out]);
        let circuit = builder.build();
        let mut witness = circuit.make_witness();
        let trit = witness.assert_constant(input.into());
        let out = chip.witness(&mut witness, [trit.into()])?[0];
        let out = witness.nop(Scalar::ZERO.into(), Scalar::ZERO.into(), out);
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, 2)?;
        let public_inputs = circuit.verify(&proof)?;
        assert_eq!(public_inputs[&Wire::Out(out)], output.into());
        Ok(())
    }

    #[test]
    fn test_ternary_selector_0() {
        assert!(test_ternary_selector_impl::<TernarySelector0>(0, 1).is_ok());
        assert!(test_ternary_selector_impl::<TernarySelector0>(1, 0).is_ok());
        assert!(test_ternary_selector_impl::<TernarySelector0>(2, 0).is_ok());
    }

    #[test]
    fn test_ternary_selector_1() {
        assert!(test_ternary_selector_impl::<TernarySelector1>(0, 0).is_ok());
        assert!(test_ternary_selector_impl::<TernarySelector1>(1, 1).is_ok());
        assert!(test_ternary_selector_impl::<TernarySelector1>(2, 0).is_ok());
    }

    #[test]
    fn test_ternary_selector_2() {
        assert!(test_ternary_selector_impl::<TernarySelector2>(0, 0).is_ok());
        assert!(test_ternary_selector_impl::<TernarySelector2>(1, 0).is_ok());
        assert!(test_ternary_selector_impl::<TernarySelector2>(2, 1).is_ok());
    }

    fn test_three_way_demux_impl(input: u8, expected: [u8; 3]) -> Result<()> {
        let mut builder = CircuitBuilder::default();
        let chip = ThreeWayDemux::default();
        let trit = builder.add_const_gate(input.into());
        let [out0, out1, out2] = chip.build(&mut builder, [trit.into()])?;
        let out = builder.add_nop_gate(out0, out1, out2);
        builder.declare_public_gates([out]);
        let circuit = builder.build();
        let mut witness = circuit.make_witness();
        let trit = witness.assert_constant(input.into());
        let [out0, out1, out2] = chip.witness(&mut witness, [trit.into()])?;
        let out = witness.nop(out0, out1, out2);
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, 2)?;
        let public_inputs = circuit.verify(&proof)?;
        assert_eq!(public_inputs[&Wire::LeftIn(out)], expected[0].into());
        assert_eq!(public_inputs[&Wire::RightIn(out)], expected[1].into());
        assert_eq!(public_inputs[&Wire::Out(out)], expected[2].into());
        Ok(())
    }

    #[test]
    fn test_three_way_demux() {
        assert!(test_three_way_demux_impl(0, [1, 0, 0]).is_ok());
        assert!(test_three_way_demux_impl(1, [0, 1, 0]).is_ok());
        assert!(test_three_way_demux_impl(2, [0, 0, 1]).is_ok());
    }

    fn test_ternary_smt<const H: usize>(
        key: u64,
        value: u64,
        path: [[Scalar; 3]; H],
        expected_root_hash: Scalar,
    ) -> Result<()> {
        let key = Scalar::from(key);
        let value = Scalar::from(value);
        let mut builder = CircuitBuilder::default();
        let inputs = builder.add_nop_gate(None, None, None);
        let key_wire = Wire::LeftIn(inputs);
        let value_wire = Wire::RightIn(inputs);
        let expected_root_hash_wire = Wire::Out(inputs);
        let chip = TernaryChip::<H>::new(path);
        let root_hash_wire =
            chip.build(&mut builder, [key_wire.into(), value_wire.into()])?[0].unwrap();
        builder.connect(root_hash_wire, expected_root_hash_wire);
        builder.declare_public_gates([inputs]);
        let circuit = builder.build();
        let mut witness = circuit.make_witness();
        witness.nop(key.into(), value.into(), expected_root_hash.into());
        chip.witness(&mut witness, [key.into(), value.into()])?;
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, 2)?;
        let public_inputs = circuit.verify(&proof)?;
        assert_eq!(public_inputs[&key_wire], key);
        assert_eq!(public_inputs[&value_wire], value);
        assert_eq!(public_inputs[&expected_root_hash_wire], expected_root_hash);
        Ok(())
    }

    #[test]
    fn test_ternary_smt_height_one_1() {
        let path = [[from_const(12), from_const(34), from_const(56)]];
        let root_hash =
            parse_scalar("0x236092ebefc7e6565e0e75414d8fdce1ce2e19bb59002d36b794b9c3111bb9cd");
        assert!(test_ternary_smt::<1>(0, 12, path, root_hash).is_ok());
        assert!(test_ternary_smt::<1>(1, 34, path, root_hash).is_ok());
        assert!(test_ternary_smt::<1>(2, 56, path, root_hash).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_2() {
        let path = [[from_const(34), from_const(56), from_const(12)]];
        let root_hash =
            parse_scalar("0x33c425faba18725cb4ffa039bbd4dade2c5b47a61edbd416cf984541f6956581");
        assert!(test_ternary_smt::<1>(0, 34, path, root_hash).is_ok());
        assert!(test_ternary_smt::<1>(1, 56, path, root_hash).is_ok());
        assert!(test_ternary_smt::<1>(2, 12, path, root_hash).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_one_3() {
        let path = [[from_const(56), from_const(78), from_const(90)]];
        let root_hash =
            parse_scalar("0x2fa39a3a76d0cf8220bd6f9899b209110ad1cca7b0bdc2b340661fa7063f2ba0");
        assert!(test_ternary_smt::<1>(0, 56, path, root_hash).is_ok());
        assert!(test_ternary_smt::<1>(1, 78, path, root_hash).is_ok());
        assert!(test_ternary_smt::<1>(2, 90, path, root_hash).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_1() {
        let path = [
            [from_const(12), from_const(34), from_const(56)],
            [
                parse_scalar("0x236092ebefc7e6565e0e75414d8fdce1ce2e19bb59002d36b794b9c3111bb9cd"),
                parse_scalar("0x33c425faba18725cb4ffa039bbd4dade2c5b47a61edbd416cf984541f6956581"),
                parse_scalar("0x2fa39a3a76d0cf8220bd6f9899b209110ad1cca7b0bdc2b340661fa7063f2ba0"),
            ],
        ];
        let root_hash =
            parse_scalar("0x27228a5e7d694f88f1b5643dd325ddcc6497f4afaa807ddb64e197742fee8cb4");
        assert!(test_ternary_smt::<2>(0, 12, path, root_hash).is_ok());
        assert!(test_ternary_smt::<2>(1, 34, path, root_hash).is_ok());
        assert!(test_ternary_smt::<2>(2, 56, path, root_hash).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_2() {
        let path = [
            [from_const(34), from_const(56), from_const(12)],
            [
                parse_scalar("0x236092ebefc7e6565e0e75414d8fdce1ce2e19bb59002d36b794b9c3111bb9cd"),
                parse_scalar("0x33c425faba18725cb4ffa039bbd4dade2c5b47a61edbd416cf984541f6956581"),
                parse_scalar("0x2fa39a3a76d0cf8220bd6f9899b209110ad1cca7b0bdc2b340661fa7063f2ba0"),
            ],
        ];
        let root_hash =
            parse_scalar("0x27228a5e7d694f88f1b5643dd325ddcc6497f4afaa807ddb64e197742fee8cb4");
        assert!(test_ternary_smt::<2>(3, 34, path, root_hash).is_ok());
        assert!(test_ternary_smt::<2>(4, 56, path, root_hash).is_ok());
        assert!(test_ternary_smt::<2>(5, 12, path, root_hash).is_ok());
    }

    #[test]
    fn test_ternary_smt_height_two_3() {
        let path = [
            [from_const(56), from_const(78), from_const(90)],
            [
                parse_scalar("0x236092ebefc7e6565e0e75414d8fdce1ce2e19bb59002d36b794b9c3111bb9cd"),
                parse_scalar("0x33c425faba18725cb4ffa039bbd4dade2c5b47a61edbd416cf984541f6956581"),
                parse_scalar("0x2fa39a3a76d0cf8220bd6f9899b209110ad1cca7b0bdc2b340661fa7063f2ba0"),
            ],
        ];
        let root_hash =
            parse_scalar("0x27228a5e7d694f88f1b5643dd325ddcc6497f4afaa807ddb64e197742fee8cb4");
        assert!(test_ternary_smt::<2>(6, 56, path, root_hash).is_ok());
        assert!(test_ternary_smt::<2>(7, 78, path, root_hash).is_ok());
        assert!(test_ternary_smt::<2>(8, 90, path, root_hash).is_ok());
    }

    trait Node: 'static + Debug + Send + Sync {
        fn hash(&self) -> Scalar;
        fn get(&self, key: Scalar) -> Scalar;
        fn get_merkle_path(&self, key: Scalar) -> Vec<Vec<Scalar>>;
        fn put(self: Arc<Self>, key: Scalar, value: Scalar) -> Arc<dyn Node>;
    }

    #[derive(Debug, Default, Copy, Clone)]
    struct Leaf(Scalar);

    impl Node for Leaf {
        fn hash(&self) -> Scalar {
            self.0
        }

        fn get(&self, key: Scalar) -> Scalar {
            assert_eq!(key, Scalar::ZERO);
            self.0
        }

        fn get_merkle_path(&self, key: Scalar) -> Vec<Vec<Scalar>> {
            assert_eq!(key, Scalar::ZERO);
            vec![]
        }

        fn put(self: Arc<Self>, key: Scalar, value: Scalar) -> Arc<dyn Node> {
            assert_eq!(key, Scalar::ZERO);
            Arc::new(Leaf(value))
        }
    }

    #[derive(Debug)]
    struct BinaryNode {
        hash: Scalar,
        left: Arc<dyn Node>,
        right: Arc<dyn Node>,
    }

    impl BinaryNode {
        fn new(left: Arc<dyn Node>, right: Arc<dyn Node>) -> Arc<dyn Node> {
            use starkom_poseidon2::{bluesky::BlueSkyConfig3, hash0};
            let hash = hash0::<BlueSkyConfig3, Scalar, 3>(&[left.hash(), right.hash()]);
            Arc::new(BinaryNode { hash, left, right })
        }
    }

    impl Node for BinaryNode {
        fn hash(&self) -> Scalar {
            self.hash
        }

        fn get(&self, key: Scalar) -> Scalar {
            let key = key.to_u256();
            let next_key = Scalar::try_from(key >> 1).unwrap();
            if key & 1.into() != U256::zero() {
                self.right.get(next_key)
            } else {
                self.left.get(next_key)
            }
        }

        fn get_merkle_path(&self, key: Scalar) -> Vec<Vec<Scalar>> {
            let key = key.to_u256();
            let next_key = Scalar::try_from(key >> 1).unwrap();
            let mut path = if key & 1.into() != U256::zero() {
                self.right.get_merkle_path(next_key)
            } else {
                self.left.get_merkle_path(next_key)
            };
            path.push(vec![self.left.hash(), self.right.hash()]);
            path
        }

        fn put(self: Arc<Self>, key: Scalar, value: Scalar) -> Arc<dyn Node> {
            let key = key.to_u256();
            let next_key = Scalar::try_from(key >> 1).unwrap();
            if key & 1.into() != U256::zero() {
                Self::new(self.left.clone(), self.right.clone().put(next_key, value))
            } else {
                Self::new(self.left.clone().put(next_key, value), self.right.clone())
            }
        }
    }

    fn get_empty_tree() -> Arc<dyn Node> {
        static TREE: LazyLock<Arc<dyn Node>> = LazyLock::new(|| {
            let mut node: Arc<dyn Node> = Arc::new(Leaf::default());
            for _ in 0..256 {
                node = BinaryNode::new(node.clone(), node.clone());
            }
            node
        });
        TREE.clone()
    }

    fn get_full_binary_prover() -> &'static Circuit {
        static PROVER: LazyLock<Circuit> = LazyLock::new(|| {
            let chip = FullBinaryChip::default();
            let mut builder = CircuitBuilder::default();
            let inputs = builder.add_nop_gate(None, None, None);
            let key_wire = Wire::LeftIn(inputs);
            let value_wire = Wire::RightIn(inputs);
            let expected_root_hash_wire = Wire::Out(inputs);
            let [root_hash_wire] = chip
                .build(&mut builder, [key_wire.into(), value_wire.into()])
                .unwrap();
            builder.connect(root_hash_wire.unwrap(), expected_root_hash_wire);
            builder.declare_public_gates([inputs]);
            builder.build()
        });
        &*PROVER
    }

    fn get_full_binary_verifier() -> &'static CompressedCircuit {
        static VERIFIER: LazyLock<CompressedCircuit> =
            LazyLock::new(|| get_full_binary_prover().compress::<Sha2Hash<Scalar>>(2));
        &*VERIFIER
    }

    fn test_full_binary_smt_impl<I: IntoIterator<Item = (u64, u64)>>(
        entries: I,
        key: u64,
    ) -> Result<()> {
        let tree = {
            let mut tree = get_empty_tree();
            for (key, value) in entries {
                tree = tree.put(key.into(), value.into());
            }
            tree
        };
        let key = key.into();
        let value = tree.get(key);
        let path: [[Scalar; 2]; 256] = tree
            .get_merkle_path(key.into())
            .into_iter()
            .map(|entry| entry.try_into().unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let root_hash = tree.hash();
        let chip = FullBinaryChip::new(path);
        let prover = get_full_binary_prover();
        let mut witness = prover.make_witness();
        let inputs = witness.nop(key.into(), value.into(), root_hash.into());
        let key_wire = Wire::LeftIn(inputs);
        let value_wire = Wire::RightIn(inputs);
        let expected_root_hash_wire = Wire::Out(inputs);
        chip.witness(&mut witness, [key.into(), value.into()])?;
        let proof = prover.prove::<Sha2Hash<Scalar>>(witness, 2)?;
        let public_inputs = get_full_binary_verifier().verify(&proof)?;
        assert_eq!(public_inputs[&key_wire], key);
        assert_eq!(public_inputs[&value_wire], value);
        assert_eq!(public_inputs[&expected_root_hash_wire], root_hash);
        Ok(())
    }

    #[test]
    fn test_full_binary_smt_empty() {
        assert!(test_full_binary_smt_impl([], 0).is_ok());
        assert!(test_full_binary_smt_impl([], 1).is_ok());
        assert!(test_full_binary_smt_impl([], 2).is_ok());
        assert!(test_full_binary_smt_impl([], 3).is_ok());
        assert!(test_full_binary_smt_impl([], 4).is_ok());
        assert!(test_full_binary_smt_impl([], 5).is_ok());
    }

    // TODO
}
