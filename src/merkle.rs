use crate::poseidon2;
use crate::xits;
use anyhow::Result;
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_plonk::{Chip as PlonkChip, CircuitBuilder, Wire, WireOrUnconstrained, Witness};

/// Runs a Merkle lookup over a Sparse Merkle Tree of height `H`.
///
/// `W` is the arity of the tree and must be either 2 or 3.
///
/// WARNING: do NOT use this chip if `H` spans the full BlueSky range, as in that case the bit or
/// trit decomposition of the key would be UNSAFE! Use the [`FullChip`] below instead.
#[derive(Debug)]
pub struct Chip<const W: usize, const H: usize> {
    decomposer: xits::BitDecomposerChip<H>,
    hasher: poseidon2::Chip<3, 2>,
    path: [[Scalar; W]; H],
}

impl<const W: usize, const H: usize> Chip<W, H> {
    /// Constructs a Merkle chip from the given Merkle path.
    pub fn new(path: [[Scalar; W]; H]) -> Self {
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher: poseidon2::Chip::default(),
            path,
        }
    }
}

impl<const W: usize, const H: usize> Default for Chip<W, H> {
    fn default() -> Self {
        Self {
            decomposer: xits::BitDecomposerChip::default(),
            hasher: poseidon2::Chip::default(),
            path: [[Scalar::ZERO; W]; H],
        }
    }
}

impl<const H: usize> PlonkChip<2, 1> for Chip<2, H> {
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
            hash = self.hasher.build(builder, [lhs.into(), rhs.into()])?[0];
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
            hash = self.hasher.witness(witness, [lhs.into(), rhs.into()])?[0];
        }
        Ok([hash])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::{from_const, parse_scalar};
    use starkom_pcs::hash::Sha2Hash;

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
        let chip = Chip::<2, H>::new(path);
        let root_hash_wire =
            chip.build(&mut builder, [key_wire.into(), value_wire.into()])?[0].unwrap();
        builder.connect(root_hash_wire, expected_root_hash_wire);
        builder.declare_public_gates([inputs]);
        let circuit = builder.build();
        let mut witness = circuit.make_witness();
        witness.nop(key.into(), value.into(), expected_root_hash.into());
        chip.witness(&mut witness, [key.into(), value.into()])?;
        let proof = circuit.prove::<Sha2Hash<Scalar>>(witness, 1)?;
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

    // TODO
}
