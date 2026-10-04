use anyhow::Result;
use primitive_types::U256;
use starkom_ff::{Field, Field32, Field64, Field128, Field256};
use starkom_plonk::{
    Cell, CellOrUnconstrained, Chip as PlonkChip, CircuitView, Constraint, WitnessView, make_const,
    rvar, var,
};
use std::marker::PhantomData;

/// Returns the smallest power of three that is >= n (returns 1 for n=0).
pub const fn next_power_of_three(n: usize) -> usize {
    let mut pow = 1usize;
    while pow < n {
        pow *= 3;
    }
    pow
}

/// Checks if a number is a power of 3.
pub const fn is_power_of_three(mut value: usize) -> bool {
    if value == 0 {
        return false;
    }
    while value > 1 {
        if value % 3 != 0 {
            return false;
        }
        value /= 3;
    }
    true
}

/// Computes the integer base-3 logarithm of `n`. For example, `ilog3(9) == 2`.
///
/// If `n` is not a power of 3 this function returns the logarithm rounded down to the nearest
/// integer, eg. `ilog3(8) == 1`.
pub const fn ilog3(mut n: usize) -> usize {
    let mut c = 0;
    while n >= 3 {
        c += 1;
        n /= 3;
    }
    c
}

/// Returns the LSB of a scalar as a scalar.
pub fn and1<F: Field>(value: F) -> F {
    value.is_odd().unwrap_u8().into()
}

/// Shifts the input [`Scalar`] to the right by `count` bits.
///
/// This is equivalent to the integer division by `2^count`.
pub fn shr<F: Field>(value: F, count: usize) -> F {
    (value.to_u256() >> count).try_into().unwrap()
}

/// Shifts the input [`Scalar`] to the right.
///
/// This is equivalent to the integer division by 2.
#[inline]
pub fn shr1<F: Field>(value: F) -> F {
    shr(value, 1)
}

/// Decomposes the input [`U256`] into `N` bits.
///
/// The `N` bits are represented as scalars and returned in little-endian order.
pub fn decompose_bits<F: Field, const N: usize>(mut value: U256) -> [F; N] {
    let mut bits = [F::ZERO; N];
    for i in 0..N {
        bits[i] = if value & 1.into() != U256::zero() {
            F::ONE
        } else {
            F::ZERO
        };
        value >>= 1;
    }
    assert_eq!(value, U256::zero());
    bits
}

/// Decomposes the input [`Scalar`] into `N` bits.
///
/// The `N` bits are represented as scalars and returned in little-endian order.
pub fn decompose_scalar_bits<F: Field, const N: usize>(value: F) -> [F; N] {
    decompose_bits::<F, N>(value.to_u256())
}

/// Decomposes the input signal into N bits.
///
/// The returned bits are in little-endian order.
///
/// WARNING: this chip is unsafe to use with a number of bits close to the field's own bit length
/// because it doesn't guard against aliasing. Use [`FullBitDecomposerChip32`],
/// [`FullBitDecomposerChip64`], [`FullBitDecomposerChip128`], or [`FullBitDecomposerChip256`] for a
/// full decomposition instead.
#[derive(Debug, Default, Clone)]
pub struct BitDecomposerChip<F: Field, const N: usize> {
    _data: PhantomData<F>,
}

impl<F: Field, const N: usize> PlonkChip<F, 1, N> for BitDecomposerChip<F, N> {
    fn width(&self) -> usize {
        N + 1
    }

    fn height(&self) -> usize {
        1
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; N]> {
        for i in 0..N {
            view.add_gate(0, var(i) * (make_const(1u8.into()) - var(i)));
        }
        view.connect(inputs[0], view.cell(0, N).into());
        let mut pow = F::ONE;
        view.add_gate(
            0,
            var(N)
                - (0..N)
                    .map(|i| {
                        let value = var(i) * make_const(pow);
                        pow = pow.double();
                        value
                    })
                    .sum::<Constraint<F>>(),
        );
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; N]> {
        let value = view.get(inputs[0]);
        decompose_scalar_bits::<F, N>(value)
            .into_iter()
            .enumerate()
            .for_each(|(i, bit)| view.set(view.cell(0, i), bit));
        view.copy(inputs[0], view.cell(0, N));
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }
}

/// Compares the number represented by the input bits against a specified constant scalar.
///
/// The inputs bits must be provided in little-endian order.
///
/// The returned signal is:
///
///  * -1 if the input value is strictly less than the constant,
///  * 0 if the input value is equal to the constant,
///  * 1 if the input value is strictly greater than the constant.
#[derive(Debug, Default, Clone)]
pub struct ConstBitComparatorChip<F: Field, const N: usize> {
    rhs: U256,
    _data: PhantomData<F>,
}

impl<F: Field, const N: usize> ConstBitComparatorChip<F, N> {
    pub fn new(rhs: U256) -> Self {
        Self {
            rhs,
            _data: PhantomData,
        }
    }
}

impl<F: Field, const N: usize> ConstBitComparatorChip<F, N> {
    fn get_rhs_bit(&self, i: usize) -> F {
        ((self.rhs >> i) & 1.into()).try_into().unwrap()
    }
}

impl<F: Field, const N: usize> PlonkChip<F, N, 1> for ConstBitComparatorChip<F, N> {
    fn width(&self) -> usize {
        N
    }

    fn height(&self) -> usize {
        2
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; N],
    ) -> Result<[Option<Cell>; 1]> {
        for i in 0..N {
            view.connect(inputs[i], view.cell(0, i).into());
        }
        view.add_gate(
            0,
            rvar(N - 1, 0) - make_const(self.get_rhs_bit(N - 1)) - rvar(N - 1, 1),
        );
        for i in (0..(N - 1)).rev() {
            let bit = self.get_rhs_bit(i);
            view.add_gate(
                0,
                (rvar(i + 1, 1) ^ 3)
                    + (make_const(1u8.into()) - (rvar(i + 1, 1) ^ 2))
                        * (rvar(i, 0) - make_const(bit))
                    - rvar(i, 1),
            );
        }
        Ok([view.cell(1, 0).into()])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; N],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        for i in 0..N {
            view.copy(inputs[i], view.cell(0, i));
        }
        view.set(
            view.cell(1, N - 1),
            view.get_at(view.cell(0, N - 1)) - self.get_rhs_bit(N - 1),
        );
        for i in (0..(N - 1)).rev() {
            let bit = self.get_rhs_bit(i);
            let cmp = view.get_at(view.cell(0, i)) - bit;
            let prev = view.get_at(view.cell(1, i + 1));
            view.set(
                view.cell(1, i),
                prev.cube() + (F::ONE - prev.square()) * cmp,
            );
        }
        Ok([view.cell(1, 0).into()])
    }
}

/// Decomposes an input signal into `N` bits, covering the full range of a field whose modulus fits
/// in `N` bits.
///
/// The returned bits are in little-endian order.
///
/// This is the shared implementation behind [`FullBitDecomposerChip32`], [`FullBitDecomposerChip64`],
/// [`FullBitDecomposerChip128`], and [`FullBitDecomposerChip256`]; use one of those instead of this type
/// directly.
#[derive(Debug, Clone)]
struct FullBitDecomposer<F: Field, const N: usize> {
    decomposer: BitDecomposerChip<F, N>,
    comparator: ConstBitComparatorChip<F, N>,
}

impl<F: Field, const N: usize> Default for FullBitDecomposer<F, N> {
    fn default() -> Self {
        Self {
            decomposer: BitDecomposerChip::default(),
            comparator: ConstBitComparatorChip::new(F::MODULUS.parse().unwrap()),
        }
    }
}

impl<F: Field, const N: usize> FullBitDecomposer<F, N> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.comparator.width())
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.comparator.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; N]> {
        let bits = view.sub_chip(0, 0, &self.decomposer, inputs)?;
        let mut view = view.sub(self.decomposer.height(), 0, None, None);
        view.sub_chip(0, 0, &self.comparator, bits)?;
        view.add_gate(self.comparator.height() - 1, var(0) + 1);
        Ok(bits)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; N]> {
        let bits = view.sub_chip(0, 0, &self.decomposer, inputs)?;
        view.sub_chip(1, 0, &self.comparator, bits)?;
        Ok(bits)
    }
}

/// Decomposes an input signal into 32 bits, covering the full range of any [`Field32`].
///
/// The returned bits are in little-endian order.
///
/// Note that the MSB will always be zero for fields whose modulus doesn't cover the upper half of
/// the 32 bit range, such as KoalaBear.
#[derive(Debug, Default, Clone)]
pub struct FullBitDecomposerChip32<F: Field32>(FullBitDecomposer<F, 32>);

impl<F: Field32> PlonkChip<F, 1, 32> for FullBitDecomposerChip32<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 32]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 32]> {
        self.0.witness(view, inputs)
    }
}

/// Decomposes an input signal into 64 bits, covering the full range of any [`Field64`].
///
/// The returned bits are in little-endian order.
#[derive(Debug, Default, Clone)]
pub struct FullBitDecomposerChip64<F: Field64>(FullBitDecomposer<F, 64>);

impl<F: Field64> PlonkChip<F, 1, 64> for FullBitDecomposerChip64<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 64]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 64]> {
        self.0.witness(view, inputs)
    }
}

/// Decomposes an input signal into 128 bits, covering the full range of any [`Field128`].
///
/// The returned bits are in little-endian order.
#[derive(Debug, Default, Clone)]
pub struct FullBitDecomposerChip128<F: Field128>(FullBitDecomposer<F, 128>);

impl<F: Field128> PlonkChip<F, 1, 128> for FullBitDecomposerChip128<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 128]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 128]> {
        self.0.witness(view, inputs)
    }
}

/// Decomposes an input signal into 256 bits, covering the full range of any [`Field256`].
///
/// The returned bits are in little-endian order.
///
/// Note that the MSB will always be zero for fields whose modulus doesn't cover the upper half of
/// the 256 bit range, such as BlueSky and Schraderbrau.
#[derive(Debug, Default, Clone)]
pub struct FullBitDecomposerChip256<F: Field256>(FullBitDecomposer<F, 256>);

impl<F: Field256> PlonkChip<F, 1, 256> for FullBitDecomposerChip256<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 256]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 256]> {
        self.0.witness(view, inputs)
    }
}

/// Divides `value`, treated as an integer, by `3^exp`, rounding down.
pub fn div_pow3<F: Field>(value: F, exp: usize) -> F {
    let dividend = value.to_u256();
    let divisor = U256::from(3).pow(exp.into());
    (dividend / divisor).try_into().unwrap()
}

/// Divides `value`, treated as an integer, by 3, rounding down.
pub fn div3<F: Field>(value: F) -> F {
    let dividend = value.to_u256();
    (dividend / 3).try_into().unwrap()
}

/// Returns `value`, treated as an integer, modulo 3.
pub fn mod3<F: Field>(value: F) -> F {
    let value = value.to_u256();
    (value % 3).try_into().unwrap()
}

/// Decomposes `value` into its `N` base-3 digits (trits), least significant first.
///
/// Panics if `value` does not fit in `N` trits.
pub fn decompose_trits<F: Field, const N: usize>(mut value: U256) -> [F; N] {
    let mut trits = [F::ZERO; N];
    for i in 0..N {
        trits[i] = F::from((value % 3).as_u32() as u8);
        value /= 3;
    }
    assert_eq!(value, U256::zero());
    trits
}

/// Like [`decompose_trits`], but takes `value` as a [`Scalar`] rather than a [`U256`].
pub fn decompose_scalar_trits<F: Field, const N: usize>(value: F) -> [F; N] {
    decompose_trits::<F, N>(value.to_u256())
}

/// Decomposes the input signal into N trits.
///
/// WARNING: this chip is unsafe to use with a number of trits close to the field's own bit length
/// because it doesn't guard against aliasing. Use [`FullTritDecomposerChip32`],
/// [`FullTritDecomposerChip64`], [`FullTritDecomposerChip128`], or [`FullTritDecomposerChip256`]
/// for a full decomposition instead, picking the one that matches the field's size.
#[derive(Debug, Default, Clone)]
pub struct TritDecomposerChip<F: Field, const N: usize> {
    _data: PhantomData<F>,
}

impl<F: Field, const N: usize> PlonkChip<F, 1, N> for TritDecomposerChip<F, N> {
    fn width(&self) -> usize {
        N + 1
    }

    fn height(&self) -> usize {
        1
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; N]> {
        for i in 0..N {
            view.add_gate(0, var(i) * (var(i) - 1) * (var(i) - 2));
        }
        view.connect(inputs[0], view.cell(0, N).into());
        let three: F = 3u8.into();
        let mut pow = F::ONE;
        view.add_gate(
            0,
            var(N)
                - (0..N)
                    .map(|i| {
                        let value = var(i) * make_const(pow);
                        pow *= three;
                        value
                    })
                    .sum::<Constraint<F>>(),
        );
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; N]> {
        let value = view.get(inputs[0]);
        decompose_scalar_trits::<F, N>(value)
            .into_iter()
            .enumerate()
            .for_each(|(i, trit)| view.set(view.cell(0, i), trit));
        view.copy(inputs[0], view.cell(0, N));
        Ok(std::array::from_fn(|i| view.cell(0, i).into()))
    }
}

/// Compares the number represented by the input trits against a specified constant scalar.
///
/// The returned signal is:
///
///  * -1 if the input value is strictly less than the constant,
///  * 0 if the input value is equal to the constant,
///  * 1 if the input value is strictly greater than the constant.
#[derive(Debug, Default, Clone)]
pub struct ConstTritComparatorChip<F: Field, const N: usize> {
    rhs: U256,
    _data: PhantomData<F>,
}

impl<F: Field, const N: usize> ConstTritComparatorChip<F, N> {
    pub fn new(rhs: U256) -> Self {
        Self {
            rhs,
            _data: PhantomData,
        }
    }

    fn get_rhs_trit(&self, i: usize) -> F {
        let three = U256::from(3);
        ((self.rhs / three.pow(i.into())) % three)
            .try_into()
            .unwrap()
    }
}

impl<F: Field, const N: usize> PlonkChip<F, N, 1> for ConstTritComparatorChip<F, N> {
    fn width(&self) -> usize {
        N
    }

    fn height(&self) -> usize {
        3
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; N],
    ) -> Result<[Option<Cell>; 1]> {
        for i in 0..N {
            view.connect(inputs[i], view.cell(0, i).into());
            let trit = self.get_rhs_trit(i);
            view.add_gate(
                0,
                ((rvar(i, 0) - make_const(trit)) * 7 - ((rvar(i, 0) - make_const(trit)) ^ 3)) / 6
                    - rvar(i, 1),
            );
        }
        view.connect(view.cell(1, N - 1).into(), view.cell(2, N - 1).into());
        for i in (0..(N - 1)).rev() {
            view.add_gate(
                1,
                rvar(i + 1, 1) + (make_const(1u8.into()) - (rvar(i + 1, 1) ^ 2)) * rvar(i, 0)
                    - rvar(i, 1),
            );
        }
        Ok([view.cell(2, 0).into()])
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; N],
    ) -> Result<[CellOrUnconstrained<F>; 1]> {
        for i in 0..N {
            view.copy(inputs[i], view.cell(0, i).into());
            let diff = view.get_at(view.cell(0, i)) - self.get_rhs_trit(i);
            view.set(
                view.cell(1, i),
                (diff * F::from(7u8) - diff.cube()) / F::from(6u8),
            );
        }
        view.copy(view.cell(1, N - 1).into(), view.cell(2, N - 1));
        for i in (0..(N - 1)).rev() {
            let cmp = view.get_at(view.cell(1, i));
            let prev = view.get_at(view.cell(2, i + 1));
            view.set(view.cell(2, i), prev + (F::ONE - prev.square()) * cmp);
        }
        Ok([view.cell(2, 0).into()])
    }
}

/// Decomposes an input signal into `N` trits, covering the full range of a field whose modulus fits
/// in `N` trits.
///
/// The returned trits are in little-endian order.
///
/// This is the shared implementation behind [`FullTritDecomposerChip32`],
/// [`FullTritDecomposerChip64`], [`FullTritDecomposerChip128`], and [`FullTritDecomposerChip256`];
/// use one of those instead of this type directly.
#[derive(Debug, Clone)]
struct FullTritDecomposer<F: Field, const N: usize> {
    decomposer: TritDecomposerChip<F, N>,
    comparator: ConstTritComparatorChip<F, N>,
}

impl<F: Field, const N: usize> Default for FullTritDecomposer<F, N> {
    fn default() -> Self {
        Self {
            decomposer: TritDecomposerChip::default(),
            comparator: ConstTritComparatorChip::new(F::MODULUS.parse().unwrap()),
        }
    }
}

impl<F: Field, const N: usize> FullTritDecomposer<F, N> {
    fn width(&self) -> usize {
        std::cmp::max(self.decomposer.width(), self.comparator.width())
    }

    fn height(&self) -> usize {
        self.decomposer.height() + self.comparator.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; N]> {
        let trits = view.sub_chip(0, 0, &self.decomposer, inputs)?;
        let mut view = view.sub(self.decomposer.height(), 0, None, None);
        view.sub_chip(0, 0, &self.comparator, trits)?;
        view.add_gate(self.comparator.height() - 1, var(0) + 1);
        Ok(trits)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; N]> {
        let trits = view.sub_chip(0, 0, &self.decomposer, inputs)?;
        view.sub_chip(self.decomposer.height(), 0, &self.comparator, trits)?;
        Ok(trits)
    }
}

/// Decomposes an input signal into 21 trits, covering the full range of any [`Field32`].
///
/// The returned trits are in little-endian order.
#[derive(Debug, Default, Clone)]
pub struct FullTritDecomposerChip32<F: Field32>(FullTritDecomposer<F, 21>);

impl<F: Field32> PlonkChip<F, 1, 21> for FullTritDecomposerChip32<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 21]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 21]> {
        self.0.witness(view, inputs)
    }
}

/// Decomposes an input signal into 41 trits, covering the full range of any [`Field64`].
///
/// The returned trits are in little-endian order.
#[derive(Debug, Default, Clone)]
pub struct FullTritDecomposerChip64<F: Field64>(FullTritDecomposer<F, 41>);

impl<F: Field64> PlonkChip<F, 1, 41> for FullTritDecomposerChip64<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 41]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 41]> {
        self.0.witness(view, inputs)
    }
}

/// Decomposes an input signal into 81 trits, covering the full range of any [`Field128`].
///
/// The returned trits are in little-endian order.
#[derive(Debug, Default, Clone)]
pub struct FullTritDecomposerChip128<F: Field128>(FullTritDecomposer<F, 81>);

impl<F: Field128> PlonkChip<F, 1, 81> for FullTritDecomposerChip128<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 81]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 81]> {
        self.0.witness(view, inputs)
    }
}

/// Decomposes an input signal into 162 trits, covering the full range of any [`Field256`].
///
/// The returned trits are in little-endian order.
#[derive(Debug, Default, Clone)]
pub struct FullTritDecomposerChip256<F: Field256>(FullTritDecomposer<F, 162>);

impl<F: Field256> PlonkChip<F, 1, 162> for FullTritDecomposerChip256<F> {
    fn width(&self) -> usize {
        self.0.width()
    }

    fn height(&self) -> usize {
        self.0.height()
    }

    fn build<G: Field256<BaseField = F>>(
        &self,
        view: &mut impl CircuitView<F, G>,
        inputs: [Option<Cell>; 1],
    ) -> Result<[Option<Cell>; 162]> {
        self.0.build(view, inputs)
    }

    fn witness(
        &self,
        view: &mut impl WitnessView<F>,
        inputs: [CellOrUnconstrained<F>; 1],
    ) -> Result<[CellOrUnconstrained<F>; 162]> {
        self.0.witness(view, inputs)
    }
}

#[cfg(all(
    test,
    feature = "bluesky",
    feature = "goldilocks",
    feature = "koalabear"
))]
mod tests {
    use super::*;
    use primitive_types::H256;
    use starkom_bluesky::{Scalar as BS, parse_scalar};
    use starkom_goldilocks::{GL, GL4};
    use starkom_koalabear::{KB, KB8};
    use starkom_pcs::hash::Sha2Hash;
    use starkom_plonk::{CircuitBuilder, Options};
    use std::cmp::Ordering;

    const BLOWUP_LOG2: usize = 1;

    #[inline]
    fn cell(row: usize, column: usize) -> Cell {
        Cell::new(row, column)
    }

    fn parse_hash(s: &'static str) -> H256 {
        s.parse().unwrap()
    }

    #[test]
    fn test_next_power_of_three() {
        assert_eq!(next_power_of_three(0), 1);
        assert_eq!(next_power_of_three(1), 1);
        assert_eq!(next_power_of_three(2), 3);
        assert_eq!(next_power_of_three(3), 3);
        assert_eq!(next_power_of_three(4), 9);
        assert_eq!(next_power_of_three(5), 9);
        assert_eq!(next_power_of_three(6), 9);
        assert_eq!(next_power_of_three(7), 9);
        assert_eq!(next_power_of_three(8), 9);
        assert_eq!(next_power_of_three(9), 9);
        assert_eq!(next_power_of_three(10), 27);
        assert_eq!(next_power_of_three(11), 27);
    }

    #[test]
    fn test_is_power_of_three() {
        assert!(!is_power_of_three(0));
        assert!(is_power_of_three(1));
        assert!(!is_power_of_three(2));
        assert!(is_power_of_three(3));
        assert!(!is_power_of_three(4));
        assert!(!is_power_of_three(5));
        assert!(!is_power_of_three(6));
        assert!(!is_power_of_three(7));
        assert!(!is_power_of_three(8));
        assert!(is_power_of_three(9));
        assert!(!is_power_of_three(10));
        assert!(!is_power_of_three(11));
    }

    #[test]
    fn test_ilog3() {
        assert_eq!(ilog3(0), 0);
        assert_eq!(ilog3(1), 0);
        assert_eq!(ilog3(2), 0);
        assert_eq!(ilog3(3), 1);
        assert_eq!(ilog3(4), 1);
        assert_eq!(ilog3(5), 1);
        assert_eq!(ilog3(6), 1);
        assert_eq!(ilog3(7), 1);
        assert_eq!(ilog3(8), 1);
        assert_eq!(ilog3(9), 2);
        assert_eq!(ilog3(10), 2);
        assert_eq!(ilog3(11), 2);
    }

    #[test]
    fn test_and1() {
        assert_eq!(and1(BS::from(42u8)), 0u8.into());
        assert_eq!(and1(BS::from(43u8)), 1u8.into());
        assert_eq!(and1(BS::from(44u8)), 0u8.into());
        assert_eq!(and1(BS::from(45u8)), 1u8.into());
        assert_eq!(and1(GL::from(42u8)), 0u8.into());
        assert_eq!(and1(GL::from(43u8)), 1u8.into());
        assert_eq!(and1(GL::from(44u8)), 0u8.into());
        assert_eq!(and1(GL::from(45u8)), 1u8.into());
    }

    #[test]
    fn test_and1_large() {
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            0u8.into()
        );
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f21"
            )),
            1u8.into()
        );
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f22"
            )),
            0u8.into()
        );
        assert_eq!(
            and1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f23"
            )),
            1u8.into()
        );
    }

    #[test]
    fn test_shr() {
        assert_eq!(
            shr(
                parse_scalar("0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"),
                4
            ),
            parse_scalar("0x00102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f2")
        );
    }

    #[test]
    fn test_shr1() {
        assert_eq!(
            shr1(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            parse_scalar("0x008101820283038404850586068707880889098a0a8b0b8c0c8d0d8e0e8f0f90")
        );
    }

    #[test]
    fn test_decompose_bits_one() {
        assert_eq!(decompose_bits::<BS, 1>(0.into()), [0u8.into()]);
        assert_eq!(decompose_bits::<BS, 1>(1.into()), [1u8.into()]);
    }

    #[test]
    fn test_decompose_bits_two() {
        assert_eq!(decompose_bits::<BS, 2>(0.into()), [0u8.into(), 0u8.into()]);
        assert_eq!(decompose_bits::<BS, 2>(1.into()), [1u8.into(), 0u8.into()]);
        assert_eq!(decompose_bits::<BS, 2>(2.into()), [0u8.into(), 1u8.into()]);
        assert_eq!(decompose_bits::<BS, 2>(3.into()), [1u8.into(), 1u8.into()]);
    }

    #[test]
    fn test_decompose_bits_three() {
        assert_eq!(
            decompose_bits::<BS, 3>(0.into()),
            [0u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_bits::<BS, 3>(1.into()),
            [1u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_bits::<BS, 3>(2.into()),
            [0u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_bits::<BS, 3>(3.into()),
            [1u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_bits::<BS, 3>(4.into()),
            [0u8.into(), 0u8.into(), 1u8.into()]
        );
        assert_eq!(
            decompose_bits::<BS, 3>(5.into()),
            [1u8.into(), 0u8.into(), 1u8.into()]
        );
        assert_eq!(
            decompose_bits::<BS, 3>(6.into()),
            [0u8.into(), 1u8.into(), 1u8.into()]
        );
        assert_eq!(
            decompose_bits::<BS, 3>(7.into()),
            [1u8.into(), 1u8.into(), 1u8.into()]
        );
    }

    #[test]
    fn test_decompose_bits_large() {
        assert_eq!(
            decompose_bits::<BS, 64>(0xFFFFFFFFFFFFFFFFu64.into()),
            [1u8.into(); 64]
        );
    }

    #[test]
    fn test_decompose_scalar_bits() {
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(0u8.into()),
            [0u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(1u8.into()),
            [1u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(2u8.into()),
            [0u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(3u8.into()),
            [1u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(4u8.into()),
            [0u8.into(), 0u8.into(), 1u8.into()]
        );
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(5u8.into()),
            [1u8.into(), 0u8.into(), 1u8.into()]
        );
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(6u8.into()),
            [0u8.into(), 1u8.into(), 1u8.into()]
        );
        assert_eq!(
            decompose_scalar_bits::<BS, 3>(7u8.into()),
            [1u8.into(), 1u8.into(), 1u8.into()]
        );
    }

    fn test_bit_decomposer_chip<F: Field, G: Field256<BaseField = F>, const N: usize>(
        value: u8,
        expected_degree_bound: usize,
        circuit_commitment: H256,
    ) {
        let chip = BitDecomposerChip::<F, N>::default();
        assert_eq!(chip.width(), N + 1);
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, G>::default();
        assert!(builder.sub_chip(0, 0, &chip, [None]).is_ok());
        builder.declare_public_cells((0..N).map(|i| cell(0, i)));
        let circuit = builder
            .build(Options {
                canonicalize_constraints: false,
                blowup_log2: BLOWUP_LOG2,
                blind: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.degree_bound(), expected_degree_bound);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let bits = witness
            .sub_chip(0, 0, &chip, [F::from(value).into()])
            .unwrap()
            .map(|bit| match bit {
                CellOrUnconstrained::Cell(cell) => witness.get_at(cell),
                _ => panic!("the output bits must be constrained"),
            });
        assert_eq!(bits, decompose_bits::<F, N>(value.into())[0..N]);
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness).unwrap();
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof).unwrap();
        assert!((0..N).all(|i| openings[&cell(0, i)] == bits[i]));
    }

    #[test]
    fn test_bit_decomposer_chip_bluesky_1() {
        let c = parse_hash("0xdc14c67633d5b395003630393a38431efede5b70f66a9ba7f2545751a95787b7");
        test_bit_decomposer_chip::<BS, BS, 1>(0, 2, c);
        test_bit_decomposer_chip::<BS, BS, 1>(1, 2, c);
    }

    #[test]
    fn test_bit_decomposer_chip_goldilocks_1() {
        let c = parse_hash("0xba7f20698c555604a69fe1d3359b4606d37dfa6eb4a452b7b198ae7eadfae6df");
        test_bit_decomposer_chip::<GL, GL4, 1>(0, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 1>(1, 2, c);
    }

    #[test]
    fn test_bit_decomposer_chip_bluesky_2() {
        let c = parse_hash("0x06a21a9eb71dffbc38b809ac6267c6519e3ce0beeac41ab8522b9b6a41d67a04");
        test_bit_decomposer_chip::<BS, BS, 2>(0, 2, c);
        test_bit_decomposer_chip::<BS, BS, 2>(1, 2, c);
        test_bit_decomposer_chip::<BS, BS, 2>(2, 2, c);
        test_bit_decomposer_chip::<BS, BS, 2>(3, 2, c);
    }

    #[test]
    fn test_bit_decomposer_chip_goldilocks_2() {
        let c = parse_hash("0x07c2e9a289c5a99fc0a6f6fced028b06f6faca261e7d5fcf5054fde074d4e1bf");
        test_bit_decomposer_chip::<GL, GL4, 2>(0, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 2>(1, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 2>(2, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 2>(3, 2, c);
    }

    #[test]
    fn test_bit_decomposer_chip_bluesky_3() {
        let c = parse_hash("0x3a7ca4209d00270b379b76458c0dc3f5179e166b438a946802da520e287c6d11");
        test_bit_decomposer_chip::<BS, BS, 3>(0, 2, c);
        test_bit_decomposer_chip::<BS, BS, 3>(1, 2, c);
        test_bit_decomposer_chip::<BS, BS, 3>(2, 2, c);
        test_bit_decomposer_chip::<BS, BS, 3>(3, 2, c);
        test_bit_decomposer_chip::<BS, BS, 3>(4, 2, c);
        test_bit_decomposer_chip::<BS, BS, 3>(5, 2, c);
        test_bit_decomposer_chip::<BS, BS, 3>(6, 2, c);
        test_bit_decomposer_chip::<BS, BS, 3>(7, 2, c);
    }

    #[test]
    fn test_bit_decomposer_chip_goldilocks_3() {
        let c = parse_hash("0xae3660f5d31368727052219f7bce7844dd72c94027f08043abe466be958a3d2d");
        test_bit_decomposer_chip::<GL, GL4, 3>(0, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 3>(1, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 3>(2, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 3>(3, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 3>(4, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 3>(5, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 3>(6, 2, c);
        test_bit_decomposer_chip::<GL, GL4, 3>(7, 2, c);
    }

    fn test_const_bit_comparator_chip<F: Field, G: Field256<BaseField = F>, const N: usize>(
        lhs: u8,
        rhs: u8,
        expected_degree_bound: usize,
        circuit_commitment: H256,
    ) {
        let mut builder = CircuitBuilder::<F, G>::default();
        let decomposer_chip = BitDecomposerChip::<F, N>::default();
        let bits = builder.sub_chip(0, 0, &decomposer_chip, [None]).unwrap();
        let comparator_chip = ConstBitComparatorChip::<F, N>::new(rhs.into());
        assert_eq!(comparator_chip.width(), N);
        assert_eq!(comparator_chip.height(), 2);
        let [cmp] = builder
            .sub_chip(decomposer_chip.height(), 0, &comparator_chip, bits)
            .unwrap();
        builder.declare_public_cells([cmp.unwrap()]);
        let circuit = builder
            .build(Options {
                canonicalize_constraints: false,
                blowup_log2: BLOWUP_LOG2,
                blind: false,
            })
            .unwrap();
        assert_eq!(
            circuit.num_rows(),
            decomposer_chip.height() + comparator_chip.height()
        );
        assert_eq!(circuit.degree_bound(), expected_degree_bound);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let bits = witness
            .sub_chip(0, 0, &decomposer_chip, [F::from(lhs).into()])
            .unwrap();
        assert!(
            witness
                .sub_chip(decomposer_chip.height(), 0, &comparator_chip, bits)
                .is_ok()
        );
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness).unwrap();
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof).unwrap();
        assert_eq!(
            openings[&cmp.unwrap()],
            match lhs.cmp(&rhs) {
                Ordering::Less => -F::from(1u8),
                Ordering::Equal => F::from(0u8),
                Ordering::Greater => F::from(1u8),
            }
        );
    }

    #[test]
    fn test_const_bit_comparator_chip_bluesky_1() {
        let c = parse_hash("0x3284d10bd08d027345d4deb12d2bb761e6607e8201c813967b55137b4f31ea4e");
        test_const_bit_comparator_chip::<BS, BS, 1>(0, 0, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 1>(1, 0, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 1>(0, 1, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 1>(1, 1, 4, c);
    }

    #[test]
    fn test_const_bit_comparator_chip_goldilocks_1() {
        let c = parse_hash("0x2f225916217af558e21e57a4478c83153d3f16c6a20ce09d766542a5e172f467");
        test_const_bit_comparator_chip::<GL, GL4, 1>(0, 0, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 1>(1, 0, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 1>(0, 1, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 1>(1, 1, 4, c);
    }

    #[test]
    fn test_const_bit_comparator_chip_bluesky_2() {
        let c = parse_hash("0x95e03b1867d1810cea10c732bde9688357ed443b4074ff465634fd38ce8c275f");
        test_const_bit_comparator_chip::<BS, BS, 2>(0, 0, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(1, 0, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(2, 0, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(3, 0, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(0, 1, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(1, 1, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(2, 1, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(3, 1, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(0, 2, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(1, 2, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(2, 2, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(3, 2, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(0, 3, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(1, 3, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(2, 3, 4, c);
        test_const_bit_comparator_chip::<BS, BS, 2>(3, 3, 4, c);
    }

    #[test]
    fn test_const_bit_comparator_chip_goldilocks_2() {
        let c = parse_hash("0x31aac105af8e62d8564c533afb0113d24f389436eb7b1d2fa99eb6a5f7a7ff51");
        test_const_bit_comparator_chip::<GL, GL4, 2>(0, 0, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(1, 0, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(2, 0, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(3, 0, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(0, 1, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(1, 1, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(2, 1, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(3, 1, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(0, 2, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(1, 2, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(2, 2, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(3, 2, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(0, 3, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(1, 3, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(2, 3, 4, c);
        test_const_bit_comparator_chip::<GL, GL4, 2>(3, 3, 4, c);
    }

    fn test_full_bit_decomposer_chip_impl<
        F: Field,
        G: Field256<BaseField = F>,
        C: PlonkChip<F, 1, N> + Default,
        const N: usize,
    >(
        value: u8,
        expected_degree_bound: usize,
        circuit_commitment: H256,
    ) {
        let chip = C::default();
        assert_eq!(chip.width(), N + 1);
        assert_eq!(chip.height(), 3);
        let mut builder = CircuitBuilder::<F, G>::default();
        assert!(builder.sub_chip(0, 0, &chip, [None]).is_ok());
        builder.declare_public_cells((0..N).map(|i| cell(0, i)));
        let circuit = builder
            .build(Options {
                canonicalize_constraints: false,
                blowup_log2: BLOWUP_LOG2,
                blind: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), 3);
        assert_eq!(circuit.degree_bound(), expected_degree_bound);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let bits = witness
            .sub_chip(0, 0, &chip, [F::from(value).into()])
            .unwrap()
            .map(|bit| match bit {
                CellOrUnconstrained::Cell(cell) => witness.get_at(cell),
                _ => panic!("the output bits must be constrained"),
            });
        assert_eq!(bits, decompose_bits::<F, N>(value.into()));
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness).unwrap();
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof).unwrap();
        assert!((0..N).all(|i| openings[&cell(0, i)] == bits[i]));
    }

    #[test]
    fn test_full_bit_decomposer_chip_bluesky() {
        type Chip = FullBitDecomposerChip256<BS>;
        let c = parse_hash("0x12d8ab57282ec992b7e59c232c3ab2efeea3f1169d77a02eac3692f19e99eeaa");
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(0, 4, c);
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(1, 4, c);
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(2, 4, c);
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(3, 4, c);
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(4, 4, c);
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(5, 4, c);
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(6, 4, c);
        test_full_bit_decomposer_chip_impl::<BS, BS, Chip, 256>(7, 4, c);
    }

    #[test]
    fn test_full_bit_decomposer_chip_goldilocks() {
        type Chip = FullBitDecomposerChip64<GL>;
        let c = parse_hash("0x41a664b29157ba96b68dc23569008ebd300bb8c5b3f24245f72f30df991b47e9");
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(0, 4, c);
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(1, 4, c);
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(2, 4, c);
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(3, 4, c);
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(4, 4, c);
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(5, 4, c);
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(6, 4, c);
        test_full_bit_decomposer_chip_impl::<GL, GL4, Chip, 64>(7, 4, c);
    }

    #[test]
    fn test_full_bit_decomposer_chip_koalabear() {
        type Chip = FullBitDecomposerChip32<KB>;
        let c = parse_hash("0xb3cd374e95df821c23027f238801e77c7cfd7b0d5094a93197d61f5e7891ce5c");
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(0, 4, c);
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(1, 4, c);
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(2, 4, c);
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(3, 4, c);
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(4, 4, c);
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(5, 4, c);
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(6, 4, c);
        test_full_bit_decomposer_chip_impl::<KB, KB8, Chip, 32>(7, 4, c);
    }

    #[test]
    fn test_div_pow3() {
        assert_eq!(
            div_pow3(
                parse_scalar("0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"),
                4
            ),
            parse_scalar("0x00032f71d3d0aac0e3aaca6871f05f0032c75591a1720ced55a4ab0058da7229")
        );
    }

    #[test]
    fn test_div3() {
        assert_eq!(
            div3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            parse_scalar("0x005601015702025803035904045a05055b06065c07075d08085e09095f0a0a60")
        );
    }

    #[test]
    fn test_mod3() {
        assert_eq!(mod3(BS::from(42u8)), 0u8.into());
        assert_eq!(mod3(BS::from(43u8)), 1u8.into());
        assert_eq!(mod3(BS::from(44u8)), 2u8.into());
        assert_eq!(mod3(BS::from(45u8)), 0u8.into());
        assert_eq!(mod3(BS::from(46u8)), 1u8.into());
        assert_eq!(mod3(BS::from(47u8)), 2u8.into());
        assert_eq!(mod3(GL::from(42u8)), 0u8.into());
        assert_eq!(mod3(GL::from(43u8)), 1u8.into());
        assert_eq!(mod3(GL::from(44u8)), 2u8.into());
        assert_eq!(mod3(GL::from(45u8)), 0u8.into());
        assert_eq!(mod3(GL::from(46u8)), 1u8.into());
        assert_eq!(mod3(GL::from(47u8)), 2u8.into());
    }

    #[test]
    fn test_mod3_large() {
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f20"
            )),
            0u8.into()
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f21"
            )),
            1u8.into()
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f22"
            )),
            2u8.into()
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f23"
            )),
            0u8.into()
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f24"
            )),
            1u8.into()
        );
        assert_eq!(
            mod3(parse_scalar(
                "0x0102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f25"
            )),
            2u8.into()
        );
    }

    #[test]
    fn test_decompose_trits_one() {
        assert_eq!(decompose_trits::<BS, 1>(0.into()), [0u8.into()]);
        assert_eq!(decompose_trits::<BS, 1>(1.into()), [1u8.into()]);
        assert_eq!(decompose_trits::<BS, 1>(2.into()), [2u8.into()]);
    }

    #[test]
    fn test_decompose_trits_two() {
        assert_eq!(decompose_trits::<BS, 2>(0.into()), [0u8.into(), 0u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(1.into()), [1u8.into(), 0u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(2.into()), [2u8.into(), 0u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(3.into()), [0u8.into(), 1u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(4.into()), [1u8.into(), 1u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(5.into()), [2u8.into(), 1u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(6.into()), [0u8.into(), 2u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(7.into()), [1u8.into(), 2u8.into()]);
        assert_eq!(decompose_trits::<BS, 2>(8.into()), [2u8.into(), 2u8.into()]);
    }

    #[test]
    fn test_decompose_trits_three() {
        assert_eq!(
            decompose_trits::<BS, 3>(0.into()),
            [0u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(1.into()),
            [1u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(2.into()),
            [2u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(3.into()),
            [0u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(4.into()),
            [1u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(5.into()),
            [2u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(6.into()),
            [0u8.into(), 2u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(7.into()),
            [1u8.into(), 2u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_trits::<BS, 3>(8.into()),
            [2u8.into(), 2u8.into(), 0u8.into()]
        );
    }

    #[test]
    fn test_decompose_scalar_trits() {
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(0u8.into()),
            [0u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(1u8.into()),
            [1u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(2u8.into()),
            [2u8.into(), 0u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(3u8.into()),
            [0u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(4u8.into()),
            [1u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(5u8.into()),
            [2u8.into(), 1u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(6u8.into()),
            [0u8.into(), 2u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(7u8.into()),
            [1u8.into(), 2u8.into(), 0u8.into()]
        );
        assert_eq!(
            decompose_scalar_trits::<BS, 3>(8u8.into()),
            [2u8.into(), 2u8.into(), 0u8.into()]
        );
    }

    fn test_trit_decomposer_chip<F: Field, G: Field256<BaseField = F>, const N: usize>(
        value: u8,
        expected_degree_bound: usize,
        circuit_commitment: H256,
    ) {
        let chip = TritDecomposerChip::<F, N>::default();
        assert_eq!(chip.width(), N + 1);
        assert_eq!(chip.height(), 1);
        let mut builder = CircuitBuilder::<F, G>::default();
        assert!(builder.sub_chip(0, 0, &chip, [None]).is_ok());
        builder.declare_public_cells((0..N).map(|i| cell(0, i)));
        let circuit = builder
            .build(Options {
                canonicalize_constraints: false,
                blowup_log2: BLOWUP_LOG2,
                blind: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), 1);
        assert_eq!(circuit.degree_bound(), expected_degree_bound);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let trits = witness
            .sub_chip(0, 0, &chip, [F::from(value).into()])
            .unwrap()
            .map(|trit| match trit {
                CellOrUnconstrained::Cell(cell) => witness.get_at(cell),
                _ => panic!("the output trits must be constrained"),
            });
        assert_eq!(trits, decompose_trits::<F, N>(value.into())[0..N]);
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness).unwrap();
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof).unwrap();
        assert!((0..N).all(|i| openings[&cell(0, i)] == trits[i]));
    }

    #[test]
    fn test_trit_decomposer_chip_bluesky_1() {
        let c = parse_hash("0xdc14c67633d5b395003630393a38431efede5b70f66a9ba7f2545751a95787b7");
        test_trit_decomposer_chip::<BS, BS, 1>(0, 2, c);
        test_trit_decomposer_chip::<BS, BS, 1>(1, 2, c);
        test_trit_decomposer_chip::<BS, BS, 1>(2, 2, c);
    }

    #[test]
    fn test_trit_decomposer_chip_goldilocks_1() {
        let c = parse_hash("0xba7f20698c555604a69fe1d3359b4606d37dfa6eb4a452b7b198ae7eadfae6df");
        test_trit_decomposer_chip::<GL, GL4, 1>(0, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 1>(1, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 1>(2, 2, c);
    }

    #[test]
    fn test_trit_decomposer_chip_bluesky_2() {
        let c = parse_hash("0x06a21a9eb71dffbc38b809ac6267c6519e3ce0beeac41ab8522b9b6a41d67a04");
        test_trit_decomposer_chip::<BS, BS, 2>(0, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(1, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(2, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(3, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(4, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(5, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(6, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(7, 2, c);
        test_trit_decomposer_chip::<BS, BS, 2>(8, 2, c);
    }

    #[test]
    fn test_trit_decomposer_chip_goldilocks_2() {
        let c = parse_hash("0x07c2e9a289c5a99fc0a6f6fced028b06f6faca261e7d5fcf5054fde074d4e1bf");
        test_trit_decomposer_chip::<GL, GL4, 2>(0, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(1, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(2, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(3, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(4, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(5, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(6, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(7, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 2>(8, 2, c);
    }

    #[test]
    fn test_trit_decomposer_chip_bluesky_3() {
        let c = parse_hash("0x3a7ca4209d00270b379b76458c0dc3f5179e166b438a946802da520e287c6d11");
        test_trit_decomposer_chip::<BS, BS, 3>(0, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(1, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(2, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(3, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(4, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(5, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(6, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(7, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(8, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(9, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(10, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(11, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(12, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(13, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(14, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(15, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(16, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(17, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(18, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(19, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(20, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(21, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(22, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(23, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(24, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(25, 2, c);
        test_trit_decomposer_chip::<BS, BS, 3>(26, 2, c);
    }

    #[test]
    fn test_trit_decomposer_chip_goldilocks_3() {
        let c = parse_hash("0xae3660f5d31368727052219f7bce7844dd72c94027f08043abe466be958a3d2d");
        test_trit_decomposer_chip::<GL, GL4, 3>(0, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(1, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(2, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(3, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(4, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(5, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(6, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(7, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(8, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(9, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(10, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(11, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(12, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(13, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(14, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(15, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(16, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(17, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(18, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(19, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(20, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(21, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(22, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(23, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(24, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(25, 2, c);
        test_trit_decomposer_chip::<GL, GL4, 3>(26, 2, c);
    }

    fn test_const_trit_comparator_chip<F: Field, G: Field256<BaseField = F>, const N: usize>(
        lhs: u8,
        rhs: u8,
        expected_degree_bound: usize,
        circuit_commitment: H256,
    ) {
        let mut builder = CircuitBuilder::<F, G>::default();
        let decomposer_chip = TritDecomposerChip::<F, N>::default();
        let trits = builder.sub_chip(0, 0, &decomposer_chip, [None]).unwrap();
        let comparator_chip = ConstTritComparatorChip::<F, N>::new(rhs.into());
        assert_eq!(comparator_chip.width(), N);
        assert_eq!(comparator_chip.height(), 3);
        let [cmp] = builder
            .sub_chip(decomposer_chip.height(), 0, &comparator_chip, trits)
            .unwrap();
        builder.declare_public_cells([cmp.unwrap()]);
        let circuit = builder
            .build(Options {
                canonicalize_constraints: false,
                blowup_log2: BLOWUP_LOG2,
                blind: false,
            })
            .unwrap();
        assert_eq!(
            circuit.num_rows(),
            decomposer_chip.height() + comparator_chip.height()
        );
        assert_eq!(circuit.degree_bound(), expected_degree_bound);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let trits = witness
            .sub_chip(0, 0, &decomposer_chip, [F::from(lhs).into()])
            .unwrap();
        assert!(
            witness
                .sub_chip(decomposer_chip.height(), 0, &comparator_chip, trits)
                .is_ok()
        );
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness).unwrap();
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof).unwrap();
        assert_eq!(
            openings[&cmp.unwrap()],
            match lhs.cmp(&rhs) {
                Ordering::Less => -F::from(1u8),
                Ordering::Equal => F::from(0u8),
                Ordering::Greater => F::from(1u8),
            }
        );
    }

    #[test]
    fn test_const_trit_comparator_chip_bluesky_1() {
        let c = parse_hash("0x4b359bece2aeaa0d815ed51d66798bb4c8eb99922585e8470f4a4cd4e650fa18");
        test_const_trit_comparator_chip::<BS, BS, 1>(0, 0, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(1, 0, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(2, 0, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(0, 1, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(1, 1, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(2, 1, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(0, 2, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(1, 2, 4, c);
        test_const_trit_comparator_chip::<BS, BS, 1>(2, 2, 4, c);
    }

    #[test]
    fn test_const_trit_comparator_chip_goldilocks_1() {
        let c = parse_hash("0x95feb533cc15d44e29e6e5c019881adec19d72736bf5a53a17e0742393878c18");
        test_const_trit_comparator_chip::<GL, GL4, 1>(0, 0, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(1, 0, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(2, 0, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(0, 1, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(1, 1, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(2, 1, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(0, 2, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(1, 2, 4, c);
        test_const_trit_comparator_chip::<GL, GL4, 1>(2, 2, 4, c);
    }

    #[test]
    fn test_const_trit_comparator_chip_bluesky_2() {
        let c = parse_hash("0xe71bafef60cfd6df00da4a434dbb72d443fa8f56cf0ff26a52e6fa99a3383017");
        for i in 0..9 {
            for j in 0..9 {
                test_const_trit_comparator_chip::<BS, BS, 2>(i, j, 4, c);
            }
        }
    }

    #[test]
    fn test_const_trit_comparator_chip_goldilocks_2() {
        let c = parse_hash("0x4e6ac69f2c6649d6217cceb7d203b54cf34bd08c7202e21b84a1ea6c947eb81a");
        for i in 0..9 {
            for j in 0..9 {
                test_const_trit_comparator_chip::<GL, GL4, 2>(i, j, 4, c);
            }
        }
    }

    fn test_full_trit_decomposer_chip_impl<
        F: Field,
        G: Field256<BaseField = F>,
        C: PlonkChip<F, 1, N> + Default,
        const N: usize,
    >(
        value: u8,
        expected_degree_bound: usize,
        circuit_commitment: H256,
    ) {
        let chip = C::default();
        assert_eq!(chip.width(), N + 1);
        assert_eq!(chip.height(), 4);
        let mut builder = CircuitBuilder::<F, G>::default();
        assert!(builder.sub_chip(0, 0, &chip, [None]).is_ok());
        builder.declare_public_cells((0..N).map(|i| cell(0, i)));
        let circuit = builder
            .build(Options {
                canonicalize_constraints: false,
                blowup_log2: BLOWUP_LOG2,
                blind: false,
            })
            .unwrap();
        assert_eq!(circuit.num_rows(), 4);
        assert_eq!(circuit.degree_bound(), expected_degree_bound);
        assert_eq!(circuit.num_columns(), N + 1);
        let mut witness = circuit.make_witness();
        let trits = witness
            .sub_chip(0, 0, &chip, [F::from(value).into()])
            .unwrap()
            .map(|trit| match trit {
                CellOrUnconstrained::Cell(cell) => witness.get_at(cell),
                _ => panic!("the output trits must be constrained"),
            });
        assert_eq!(trits, decompose_trits::<F, N>(value.into()));
        circuit.check_witness(&witness).unwrap();
        let proof = circuit.prove::<Sha2Hash<G>>(witness).unwrap();
        let circuit = circuit.to_compressed::<Sha2Hash<G>>();
        assert_eq!(circuit.commitment(), circuit_commitment);
        let openings = circuit.verify(&proof).unwrap();
        assert!((0..N).all(|i| openings[&cell(0, i)] == trits[i]));
    }

    #[test]
    fn test_full_trit_decomposer_chip_bluesky() {
        type Chip = FullTritDecomposerChip256<BS>;
        let c = parse_hash("0x74254fa4297fe0e848462886f58e08f32f0b20bc1d73965b9d67c014499f5ea9");
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(0, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(1, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(2, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(3, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(4, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(5, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(6, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(7, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(8, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(9, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(10, 4, c);
        test_full_trit_decomposer_chip_impl::<BS, BS, Chip, 162>(11, 4, c);
    }

    #[test]
    fn test_full_trit_decomposer_chip_goldilocks() {
        type Chip = FullTritDecomposerChip64<GL>;
        let c = parse_hash("0xda770da42421fb74bc9d87d907cdded42a7aa795f0eaa1a31d9cd54bbd04d084");
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(0, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(1, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(2, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(3, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(4, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(5, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(6, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(7, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(8, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(9, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(10, 4, c);
        test_full_trit_decomposer_chip_impl::<GL, GL4, Chip, 41>(11, 4, c);
    }

    #[test]
    fn test_full_trit_decomposer_chip_koalabear() {
        type Chip = FullTritDecomposerChip32<KB>;
        let c = parse_hash("0xa2ad1f596c4fc73cebf75ed02aaf9d3950dea1b704c63917bd17d94ebd8562b0");
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(0, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(1, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(2, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(3, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(4, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(5, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(6, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(7, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(8, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(9, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(10, 4, c);
        test_full_trit_decomposer_chip_impl::<KB, KB8, Chip, 21>(11, 4, c);
    }
}
