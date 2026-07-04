use crate::poseidon2;
use crate::utils;
use anyhow::{Result, anyhow};
use starkom_bluesky::Scalar;
use starkom_ff::Field;
use starkom_pcs::hash::{Hash, Poseidon2Hash, Sha2Hash};
use starkom_plonk::{
    Chip as PlonkChip, Circuit, CircuitBuilder, CompressedCircuit, Proof as PlonkProof, Wire,
    WireOrUnconstrained, Witness,
};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::sync::{Arc, LazyLock, Mutex, OnceLock};

/// Domain separator tag for the identity commitment.
static IDENTITY_DST: LazyLock<Scalar> =
    LazyLock::new(|| utils::hash_to_scalar(b"starkom/hmac/identity"));

/// Domain separator tag for the HMAC.
static AUTHENTICATION_DST: LazyLock<Scalar> =
    LazyLock::new(|| utils::hash_to_scalar(b"starkom/hmac/authentication"));

/// Poseidon2 T=3 instance over the BlueSky field.
///
/// We use this instance for hashing both the identity commitment and the HMAC.
fn hash_t3(inputs: &[Scalar]) -> Scalar {
    starkom_poseidon2::hash::<starkom_poseidon2::bluesky::BlueSkyConfig3, Scalar, 3>(inputs)[0]
}

/// Returns the identity commitment of a secret key.
///
/// The commitment is calculated by hashing of the secret key with a domain separator tag.
pub fn identity_commitment(secret_key: Scalar) -> Scalar {
    hash_t3(&[*IDENTITY_DST, secret_key])
}

/// Proves the HMAC of a message.
///
/// NOTE: the generic parameter `N` is the length of the message and the generic parameter `M` MUST
/// be equal to `N+3`. `M` is used as the size of the internal HMAC hash and is equal to the size of
/// the signed message plus 3 scalars for: the DST, the secret key, and the message length.
#[derive(Debug, Default, Clone)]
pub struct Chip<const N: usize, const M: usize> {
    secret_key: Scalar,
    identity_hasher: poseidon2::Chip<3, 2>,
    authentication_hasher: poseidon2::Chip<3, M>,
}

impl<const N: usize, const M: usize> Chip<N, M> {
    pub fn new(secret_key: Scalar) -> Self {
        Self {
            secret_key,
            identity_hasher: poseidon2::Chip::default(),
            authentication_hasher: poseidon2::Chip::default(),
        }
    }
}

impl<const N: usize, const M: usize> PlonkChip<N, 2> for Chip<N, M> {
    fn build(
        &self,
        builder: &mut CircuitBuilder,
        inputs: [Option<Wire>; N],
    ) -> Result<[Option<Wire>; 2]> {
        let identity_dst = builder.add_const_gate(*IDENTITY_DST);
        let authentication_dst = builder.add_const_gate(*AUTHENTICATION_DST);
        let secret_key = Wire::Out(builder.add_nop_gate(None, None, None));
        let length = builder.add_const_gate(Scalar::from_const(N as u64));
        let identity_commitment = self
            .identity_hasher
            .build(builder, [identity_dst.into(), secret_key.into()])?[0];
        let signature = self.authentication_hasher.build(
            builder,
            [authentication_dst.into(), secret_key.into(), length.into()]
                .into_iter()
                .chain(inputs.into_iter())
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
        )?[0];
        Ok([identity_commitment, signature])
    }

    fn witness(
        &self,
        witness: &mut Witness,
        inputs: [WireOrUnconstrained; N],
    ) -> Result<[WireOrUnconstrained; 2]> {
        let identity_dst = witness.assert_constant(*IDENTITY_DST);
        let authentication_dst = witness.assert_constant(*AUTHENTICATION_DST);
        let secret_key = Wire::Out(witness.nop(
            Scalar::ZERO.into(),
            Scalar::ZERO.into(),
            self.secret_key.into(),
        ));
        let length = witness.assert_constant(Scalar::from_const(N as u64));
        let identity_commitment = self
            .identity_hasher
            .witness(witness, [identity_dst.into(), secret_key.into()])?[0];
        let signature = self.authentication_hasher.witness(
            witness,
            [authentication_dst.into(), secret_key.into(), length.into()]
                .into_iter()
                .chain(inputs.into_iter())
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
        )?[0];
        Ok([identity_commitment, signature])
    }
}

/// A PLONK circuit to sign messages with a secret key using zkMAC.
#[derive(Debug, Clone)]
pub struct SignerCircuit<const N: usize, const M: usize> {
    chip: Chip<N, M>,
    inner: Circuit,
    identity_commitment_wire: Wire,
    message_wires: [Wire; N],
    signature_wire: Wire,
}

impl<const N: usize, const M: usize> SignerCircuit<N, M> {
    fn make(secret_key: Option<Scalar>) -> Self {
        let mut builder = CircuitBuilder::default();
        let message: [Option<Wire>; N] =
            std::array::from_fn(|_| Some(Wire::Out(builder.add_nop_gate(None, None, None))));
        let chip = match secret_key {
            Some(secret_key) => Chip::<N, M>::new(secret_key),
            None => Chip::<N, M>::default(),
        };
        let [identity_commitment, signature] = chip.build(&mut builder, message).unwrap();
        builder.declare_public_gates(
            message
                .into_iter()
                .chain([identity_commitment, signature])
                .map(|wire| wire.unwrap().gate()),
        );
        Self {
            chip,
            inner: builder.build(),
            identity_commitment_wire: identity_commitment.unwrap(),
            message_wires: message.map(|wire| wire.unwrap()),
            signature_wire: signature.unwrap(),
        }
    }

    pub fn new(secret_key: Scalar) -> Self {
        Self::make(Some(secret_key))
    }

    pub fn to_verifier<H: Hash<Scalar>>(self, blowup_log2: usize) -> VerifierCircuit<H, N> {
        VerifierCircuit {
            inner: self.inner.to_compressed::<H>(blowup_log2),
            identity_commitment_wire: self.identity_commitment_wire,
            message_wires: self.message_wires,
            signature_wire: self.signature_wire,
            _data: Default::default(),
        }
    }
}

impl<const N: usize, const M: usize> Default for SignerCircuit<N, M> {
    fn default() -> Self {
        Self::make(None)
    }
}

/// Dyn trait for a generic signer circuit.
///
/// This trait abstracts the generic parameters `N` and `M` away so that signers for all sizes can
/// be cached in the same data structure.
pub trait AbstractSignerCircuit<H: Hash<Scalar>>: Debug + Send + Sync {
    /// Signs the provided message and returns the signature scalar and the zkSTARK proof.
    fn witness(&self, message: &[Scalar], blowup_log2: usize) -> (Scalar, PlonkProof<H>);
}

impl<H: Hash<Scalar>, const N: usize, const M: usize> AbstractSignerCircuit<H>
    for SignerCircuit<N, M>
{
    fn witness(&self, message: &[Scalar], blowup_log2: usize) -> (Scalar, PlonkProof<H>) {
        let mut witness = self.inner.make_witness();
        let message = std::array::from_fn(|i| {
            Wire::Out(witness.nop(Scalar::ZERO.into(), Scalar::ZERO.into(), message[i].into()))
                .into()
        });
        self.chip.witness(&mut witness, message).unwrap();
        let signature = witness.get(self.signature_wire);
        let proof = self.inner.prove::<H>(witness, blowup_log2).unwrap();
        (signature, proof)
    }
}

#[derive(Debug, Clone)]
pub struct VerifierCircuit<H: Hash<Scalar>, const N: usize> {
    inner: CompressedCircuit,
    identity_commitment_wire: Wire,
    message_wires: [Wire; N],
    signature_wire: Wire,
    _data: PhantomData<H>,
}

/// Dyn trait for a generic verifier circuit.
///
/// This trait abstracts the generic parameter `N` away so that verifiers for all sizes can be
/// cached in the same data structure.
pub trait AbstractVerifierCircuit<H: Hash<Scalar>>: Debug + Send + Sync {
    fn verify(
        &self,
        identity_commitment: Scalar,
        message: &[Scalar],
        signature: Scalar,
        proof: &PlonkProof<H>,
    ) -> Result<()>;
}

impl<H: Hash<Scalar> + Debug + Send + Sync, const N: usize> AbstractVerifierCircuit<H>
    for VerifierCircuit<H, N>
{
    fn verify(
        &self,
        identity_commitment: Scalar,
        message: &[Scalar],
        signature: Scalar,
        proof: &PlonkProof<H>,
    ) -> Result<()> {
        let openings = self.inner.verify(proof)?;
        if openings[&self.identity_commitment_wire] != identity_commitment {
            return Err(anyhow!(
                "identity commitment mismatch (got {}, want {})",
                openings[&self.identity_commitment_wire],
                identity_commitment
            ));
        }
        if self
            .message_wires
            .iter()
            .enumerate()
            .any(|(i, value)| openings[value] != message[i])
        {
            return Err(anyhow!("signed message mismatch"));
        }
        if openings[&self.signature_wire] != signature {
            return Err(anyhow!(
                "signature mismatch (got {}, want {})",
                openings[&self.signature_wire],
                signature
            ));
        }
        Ok(())
    }
}

/// Caches signer circuits.
#[derive(Debug)]
struct SignerCache<H: Hash<Scalar>> {
    /// The first component of the key is the secret key used to sign; the second one is the size of
    /// the signed input, `N`.
    circuits: Mutex<BTreeMap<(Scalar, usize), Arc<OnceLock<Arc<dyn AbstractSignerCircuit<H>>>>>>,
    _data: PhantomData<H>,
}

impl<H: Hash<Scalar>> SignerCache<H> {
    fn make<const N: usize, const M: usize>(
        secret_key: Scalar,
    ) -> Arc<dyn AbstractSignerCircuit<H>> {
        Arc::new(SignerCircuit::<N, M>::new(secret_key))
    }

    fn make_monomorphic(secret_key: Scalar, n: usize) -> Arc<dyn AbstractSignerCircuit<H>> {
        match n {
            1 => Self::make::<1, 4>(secret_key),
            2 => Self::make::<2, 5>(secret_key),
            3 => Self::make::<3, 6>(secret_key),
            4 => Self::make::<4, 7>(secret_key),
            5 => Self::make::<5, 8>(secret_key),
            6 => Self::make::<6, 9>(secret_key),
            7 => Self::make::<7, 10>(secret_key),
            8 => Self::make::<8, 11>(secret_key),
            9 => Self::make::<9, 12>(secret_key),
            10 => Self::make::<10, 13>(secret_key),
            11 => Self::make::<11, 14>(secret_key),
            12 => Self::make::<12, 15>(secret_key),
            13 => Self::make::<13, 16>(secret_key),
            14 => Self::make::<14, 17>(secret_key),
            15 => Self::make::<15, 18>(secret_key),
            16 => Self::make::<16, 19>(secret_key),
            17 => Self::make::<17, 20>(secret_key),
            18 => Self::make::<18, 21>(secret_key),
            19 => Self::make::<19, 22>(secret_key),
            20 => Self::make::<20, 23>(secret_key),
            _ => unimplemented!(),
        }
    }

    fn get_or_make(
        &self,
        secret_key: Scalar,
        n: usize,
    ) -> Arc<OnceLock<Arc<dyn AbstractSignerCircuit<H>>>> {
        let key = (secret_key, n);
        let mut circuits = self.circuits.lock().unwrap();
        if !circuits.contains_key(&key) {
            circuits.insert(key, Arc::default());
        }
        circuits.get(&key).unwrap().clone()
    }
}

impl SignerCache<Sha2Hash<Scalar>> {
    fn get(secret_key: Scalar, n: usize) -> Arc<dyn AbstractSignerCircuit<Sha2Hash<Scalar>>> {
        static CACHE: LazyLock<SignerCache<Sha2Hash<Scalar>>> = LazyLock::new(|| SignerCache {
            circuits: Mutex::default(),
            _data: Default::default(),
        });
        let once_lock = CACHE.get_or_make(secret_key, n);
        once_lock
            .get_or_init(|| Self::make_monomorphic(secret_key, n))
            .clone()
    }
}

impl SignerCache<Poseidon2Hash<Scalar>> {
    fn get(secret_key: Scalar, n: usize) -> Arc<dyn AbstractSignerCircuit<Poseidon2Hash<Scalar>>> {
        static CACHE: LazyLock<SignerCache<Poseidon2Hash<Scalar>>> =
            LazyLock::new(|| SignerCache {
                circuits: Mutex::default(),
                _data: Default::default(),
            });
        let once_lock = CACHE.get_or_make(secret_key, n);
        once_lock
            .get_or_init(|| Self::make_monomorphic(secret_key, n))
            .clone()
    }
}

/// Caches verifier circuits.
#[derive(Debug)]
struct VerifierCache<H: Hash<Scalar> + 'static> {
    /// The first component of the key is the size of the signed input, `N`; the second one is the
    /// log2 of the blowup factor for the circuit.
    circuits: Mutex<BTreeMap<(usize, usize), Arc<OnceLock<Arc<dyn AbstractVerifierCircuit<H>>>>>>,
    _data: PhantomData<H>,
}

impl<H: Hash<Scalar> + Debug + Send + Sync> VerifierCache<H> {
    fn make<const N: usize, const M: usize>(
        blowup_log2: usize,
    ) -> Arc<dyn AbstractVerifierCircuit<H>> {
        Arc::new(SignerCircuit::<N, M>::default().to_verifier(blowup_log2))
    }

    fn make_monomorphic(n: usize, blowup_log2: usize) -> Arc<dyn AbstractVerifierCircuit<H>> {
        match n {
            1 => Self::make::<1, 4>(blowup_log2),
            2 => Self::make::<2, 5>(blowup_log2),
            3 => Self::make::<3, 6>(blowup_log2),
            4 => Self::make::<4, 7>(blowup_log2),
            5 => Self::make::<5, 8>(blowup_log2),
            6 => Self::make::<6, 9>(blowup_log2),
            7 => Self::make::<7, 10>(blowup_log2),
            8 => Self::make::<8, 11>(blowup_log2),
            9 => Self::make::<9, 12>(blowup_log2),
            10 => Self::make::<10, 13>(blowup_log2),
            11 => Self::make::<11, 14>(blowup_log2),
            12 => Self::make::<12, 15>(blowup_log2),
            13 => Self::make::<13, 16>(blowup_log2),
            14 => Self::make::<14, 17>(blowup_log2),
            15 => Self::make::<15, 18>(blowup_log2),
            16 => Self::make::<16, 19>(blowup_log2),
            17 => Self::make::<17, 20>(blowup_log2),
            18 => Self::make::<18, 21>(blowup_log2),
            19 => Self::make::<19, 22>(blowup_log2),
            20 => Self::make::<20, 23>(blowup_log2),
            _ => unimplemented!(),
        }
    }

    fn get_or_make(
        &self,
        n: usize,
        blowup_log2: usize,
    ) -> Arc<OnceLock<Arc<dyn AbstractVerifierCircuit<H>>>> {
        let key = (n, blowup_log2);
        let mut circuits = self.circuits.lock().unwrap();
        if !circuits.contains_key(&key) {
            circuits.insert(key, Arc::default());
        }
        circuits.get(&key).unwrap().clone()
    }
}

impl VerifierCache<Sha2Hash<Scalar>> {
    fn get(n: usize, blowup_log2: usize) -> Arc<dyn AbstractVerifierCircuit<Sha2Hash<Scalar>>> {
        static CACHE: LazyLock<VerifierCache<Sha2Hash<Scalar>>> = LazyLock::new(|| VerifierCache {
            circuits: Mutex::default(),
            _data: Default::default(),
        });
        let once_lock = CACHE.get_or_make(n, blowup_log2);
        once_lock
            .get_or_init(|| Self::make_monomorphic(n, blowup_log2))
            .clone()
    }
}

impl VerifierCache<Poseidon2Hash<Scalar>> {
    fn get(
        n: usize,
        blowup_log2: usize,
    ) -> Arc<dyn AbstractVerifierCircuit<Poseidon2Hash<Scalar>>> {
        static CACHE: LazyLock<VerifierCache<Poseidon2Hash<Scalar>>> =
            LazyLock::new(|| VerifierCache {
                circuits: Mutex::default(),
                _data: Default::default(),
            });
        let once_lock = CACHE.get_or_make(n, blowup_log2);
        once_lock
            .get_or_init(|| Self::make_monomorphic(n, blowup_log2))
            .clone()
    }
}

/// Utility trait to fetch a circuit from the cache for all known hash backends.
///
/// For internal use only. It's declared as public because it needs to surface in the signature of
/// `sign`.
pub trait CachedSigner: Hash<Scalar> {
    fn get_signer_circuit(secret_key: Scalar, n: usize) -> Arc<dyn AbstractSignerCircuit<Self>>;
}

impl CachedSigner for Sha2Hash<Scalar> {
    fn get_signer_circuit(secret_key: Scalar, n: usize) -> Arc<dyn AbstractSignerCircuit<Self>> {
        SignerCache::<Sha2Hash<Scalar>>::get(secret_key, n)
    }
}

impl CachedSigner for Poseidon2Hash<Scalar> {
    fn get_signer_circuit(secret_key: Scalar, n: usize) -> Arc<dyn AbstractSignerCircuit<Self>> {
        SignerCache::<Poseidon2Hash<Scalar>>::get(secret_key, n)
    }
}

/// Signs a message with the given secret key using the zkMAC scheme, which consists of making a
/// zkSTARK proof of an HMAC (hash-based message authentication code).
///
/// The secret key remains a private input of the zkSTARK circuit. Other parties can verify the
/// zkSTARK proof, and therefore the signature, without ever seeing the secret key.
///
/// The returned pair contains the signature scalar and the zkSTARK proof.
pub fn sign<H: Hash<Scalar> + CachedSigner, const N: usize>(
    secret_key: Scalar,
    message: &[Scalar; N],
    blowup_log2: usize,
) -> (Scalar, PlonkProof<H>) {
    let circuit = H::get_signer_circuit(secret_key, N);
    circuit.witness(message, blowup_log2)
}

/// Utility trait to fetch a circuit from the cache for all known hash backends.
///
/// For internal use only. It's declared as public because it needs to surface in the signature of
/// `verify`.
pub trait CachedVerifier: Hash<Scalar> {
    fn get_verifier_circuit(n: usize, blowup_log2: usize)
    -> Arc<dyn AbstractVerifierCircuit<Self>>;
}

impl CachedVerifier for Sha2Hash<Scalar> {
    fn get_verifier_circuit(
        n: usize,
        blowup_log2: usize,
    ) -> Arc<dyn AbstractVerifierCircuit<Self>> {
        VerifierCache::<Sha2Hash<Scalar>>::get(n, blowup_log2)
    }
}

impl CachedVerifier for Poseidon2Hash<Scalar> {
    fn get_verifier_circuit(
        n: usize,
        blowup_log2: usize,
    ) -> Arc<dyn AbstractVerifierCircuit<Self>> {
        VerifierCache::<Poseidon2Hash<Scalar>>::get(n, blowup_log2)
    }
}

/// Verifies a message signed with our zkMAC scheme.
///
/// The `identity_commitment` is the same value returned by the `identity_commitment()` function
/// above for the secret key used to sign the message, and is one of the public inputs of the
/// circuit. This function verifies the zkSTARK proof and checks that the identity commitment and
/// message match the public inputs.
pub fn verify<H: Hash<Scalar> + CachedVerifier, const N: usize>(
    identity_commitment: Scalar,
    message: &[Scalar; N],
    signature: Scalar,
    proof: &PlonkProof<H>,
) -> Result<()> {
    let circuit = H::get_verifier_circuit(N, proof.blowup_log2());
    circuit.verify(identity_commitment, message, signature, proof)
}

#[cfg(test)]
mod tests {
    use super::*;
    use starkom_bluesky::{from_const, parse_scalar};

    const BLOWUP_LOG2: usize = 1;

    #[test]
    fn test_identity_commitment() {
        assert_eq!(
            identity_commitment(parse_scalar(
                "0x57833b8d7c2e4b4fa73b9b701496c153628a211d802f02e3f948437627123680"
            )),
            parse_scalar("0x300a3af1081f62e3bcb20da3798a2b50ffce22d2ac24b9731aeca3bf8ac88837")
        );
        assert_eq!(
            identity_commitment(parse_scalar(
                "0x75b9be52955cd0933248b6324139bb988e9442fcf657c391b227de67f7d84561"
            )),
            parse_scalar("0x65d90105b752f620a9fcee49e878df81ef3b7e9a53c4752b9b16bba2d3b04f93")
        );
    }

    #[test]
    fn test_signature_one_scalar_sha2() {
        let secret_key =
            parse_scalar("0x57833b8d7c2e4b4fa73b9b701496c153628a211d802f02e3f948437627123680");
        let (signature, proof) =
            sign::<Sha2Hash<Scalar>, 1>(secret_key, &[from_const(42)], BLOWUP_LOG2);
        assert_eq!(
            signature,
            parse_scalar("0x232869fb537eb95b880e872fe7683145c887b26fedd0427791d6480048d57b67")
        );
        assert!(
            verify::<Sha2Hash<Scalar>, 1>(
                identity_commitment(secret_key),
                &[from_const(42)],
                signature,
                &proof
            )
            .is_ok()
        );
    }

    #[test]
    fn test_signature_one_scalar_poseidon2() {
        let secret_key =
            parse_scalar("0x57833b8d7c2e4b4fa73b9b701496c153628a211d802f02e3f948437627123680");
        let (signature, proof) =
            sign::<Poseidon2Hash<Scalar>, 1>(secret_key, &[from_const(42)], BLOWUP_LOG2);
        assert_eq!(
            signature,
            parse_scalar("0x232869fb537eb95b880e872fe7683145c887b26fedd0427791d6480048d57b67")
        );
        assert!(
            verify::<Poseidon2Hash<Scalar>, 1>(
                identity_commitment(secret_key),
                &[from_const(42)],
                signature,
                &proof
            )
            .is_ok()
        );
    }

    #[test]
    fn test_signature_one_scalar_different_key_sha2() {
        let secret_key =
            parse_scalar("0x75b9be52955cd0933248b6324139bb988e9442fcf657c391b227de67f7d84561");
        let (signature, proof) =
            sign::<Sha2Hash<Scalar>, 1>(secret_key, &[from_const(42)], BLOWUP_LOG2);
        assert_eq!(
            signature,
            parse_scalar("0x20d84c0bb6bb0ee23f73141017a95539a4a6e68fdf5b79ff8d993e1adaa148fb")
        );
        assert!(
            verify::<Sha2Hash<Scalar>, 1>(
                identity_commitment(secret_key),
                &[from_const(42)],
                signature,
                &proof
            )
            .is_ok()
        );
    }

    #[test]
    fn test_signature_one_scalar_different_key_poseidon2() {
        let secret_key =
            parse_scalar("0x75b9be52955cd0933248b6324139bb988e9442fcf657c391b227de67f7d84561");
        let (signature, proof) =
            sign::<Poseidon2Hash<Scalar>, 1>(secret_key, &[from_const(42)], BLOWUP_LOG2);
        assert_eq!(
            signature,
            parse_scalar("0x20d84c0bb6bb0ee23f73141017a95539a4a6e68fdf5b79ff8d993e1adaa148fb")
        );
        assert!(
            verify::<Poseidon2Hash<Scalar>, 1>(
                identity_commitment(secret_key),
                &[from_const(42)],
                signature,
                &proof
            )
            .is_ok()
        );
    }

    #[test]
    fn test_signature_another_scalar_sha2() {
        let secret_key =
            parse_scalar("0x57833b8d7c2e4b4fa73b9b701496c153628a211d802f02e3f948437627123680");
        let (signature, proof) =
            sign::<Sha2Hash<Scalar>, 1>(secret_key, &[from_const(123)], BLOWUP_LOG2);
        assert_eq!(
            signature,
            parse_scalar("0x7395bb057b232be412d3db8aa35e54f4fa44b140d7cfe05d01853fe81b7401ef")
        );
        assert!(
            verify::<Sha2Hash<Scalar>, 1>(
                identity_commitment(secret_key),
                &[from_const(123)],
                signature,
                &proof
            )
            .is_ok()
        );
    }

    #[test]
    fn test_signature_another_scalar_poseidon2() {
        let secret_key =
            parse_scalar("0x57833b8d7c2e4b4fa73b9b701496c153628a211d802f02e3f948437627123680");
        let (signature, proof) =
            sign::<Poseidon2Hash<Scalar>, 1>(secret_key, &[from_const(123)], BLOWUP_LOG2);
        assert_eq!(
            signature,
            parse_scalar("0x7395bb057b232be412d3db8aa35e54f4fa44b140d7cfe05d01853fe81b7401ef")
        );
        assert!(
            verify::<Poseidon2Hash<Scalar>, 1>(
                identity_commitment(secret_key),
                &[from_const(123)],
                signature,
                &proof
            )
            .is_ok()
        );
    }

    #[test]
    fn test_signature_two_scalars_sha2() {
        let secret_key =
            parse_scalar("0x57833b8d7c2e4b4fa73b9b701496c153628a211d802f02e3f948437627123680");
        let (signature, proof) = sign::<Sha2Hash<Scalar>, 2>(
            secret_key,
            &[from_const(123), from_const(456)],
            BLOWUP_LOG2,
        );
        assert_eq!(
            signature,
            parse_scalar("0x48983c3c58b69a1354453ccc59a6f828c397627513937db41d2f32c4401723b6")
        );
        assert!(
            verify::<Sha2Hash<Scalar>, 2>(
                identity_commitment(secret_key),
                &[from_const(123), from_const(456)],
                signature,
                &proof
            )
            .is_ok()
        );
    }

    #[test]
    fn test_signature_two_scalars_poseidon2() {
        let secret_key =
            parse_scalar("0x57833b8d7c2e4b4fa73b9b701496c153628a211d802f02e3f948437627123680");
        let (signature, proof) = sign::<Poseidon2Hash<Scalar>, 2>(
            secret_key,
            &[from_const(123), from_const(456)],
            BLOWUP_LOG2,
        );
        assert_eq!(
            signature,
            parse_scalar("0x48983c3c58b69a1354453ccc59a6f828c397627513937db41d2f32c4401723b6")
        );
        assert!(
            verify::<Poseidon2Hash<Scalar>, 2>(
                identity_commitment(secret_key),
                &[from_const(123), from_const(456)],
                signature,
                &proof
            )
            .is_ok()
        );
    }
}
