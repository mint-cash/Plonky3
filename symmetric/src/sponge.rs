use alloc::string::String;
use core::marker::PhantomData;

// use itertools::Itertools;

use alloc::vec::Vec;
use core::iter::IntoIterator;


use p3_field::{reduce_32, Field, PrimeField, PrimeField32};

use crate::hasher::CryptographicHasher;
use crate::permutation::CryptographicPermutation;

/// A padding-free, overwrite-mode sponge function.
///
/// `WIDTH` is the sponge's rate plus the sponge's capacity.
#[derive(Clone, Debug)]
pub struct PaddingFreeSponge<P, const WIDTH: usize, const RATE: usize, const OUT: usize> {
    permutation: P,
}

impl<P, const WIDTH: usize, const RATE: usize, const OUT: usize>
    PaddingFreeSponge<P, WIDTH, RATE, OUT>
{
    pub const fn new(permutation: P) -> Self {
        Self { permutation }
    }
}

impl<T, P, const WIDTH: usize, const RATE: usize, const OUT: usize> CryptographicHasher<T, [T; OUT]>
    for PaddingFreeSponge<P, WIDTH, RATE, OUT>
where
    T: Default + Copy,
    P: CryptographicPermutation<[T; WIDTH]>,
{
    fn hash_iter<I>(&self, input: I) -> [T; OUT]
    where
        I: IntoIterator<Item = T>,
    {
        let mut state = [T::default(); WIDTH];
        let mut input_iter = input.into_iter();

        loop {
            let mut absorbed = 0;
            for s in state.iter_mut().take(RATE) {
                if let Some(i) = input_iter.next() {
                    *s = i;
                    absorbed += 1;
                } else {
                    break;
                }
            }
            if absorbed == 0 {
                break;
            }
            state = self.permutation.permute(state);
        }

        state[..OUT].try_into().unwrap()
    }
}

/// A padding-free, overwrite-mode sponge function that operates natively over PF but accepts elements
/// of F: PrimeField32.
///
/// `WIDTH` is the sponge's rate plus the sponge's capacity.
#[derive(Clone, Debug)]
pub struct MultiField32PaddingFreeSponge<
    F,
    PF,
    P,
    const WIDTH: usize,
    const RATE: usize,
    const OUT: usize,
> {
    permutation: P,
    num_f_elms: usize,
    _phantom: PhantomData<(F, PF)>,
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize, const OUT: usize>
    MultiField32PaddingFreeSponge<F, PF, P, WIDTH, RATE, OUT>
where
    F: PrimeField32,
    PF: Field,
{
    pub fn new(permutation: P) -> Result<Self, String> {
        if F::order() >= PF::order() {
            return Err(String::from("F::order() must be less than PF::order()"));
        }

        let num_f_elms = PF::bits() / F::bits();
        Ok(Self {
            permutation,
            num_f_elms,
            _phantom: PhantomData,
        })
    }
}

impl<F, PF, P, const WIDTH: usize, const RATE: usize, const OUT: usize>
    CryptographicHasher<F, [PF; OUT]> for MultiField32PaddingFreeSponge<F, PF, P, WIDTH, RATE, OUT>
where
    F: PrimeField32,
    PF: PrimeField + Default + Copy,
    P: CryptographicPermutation<[PF; WIDTH]>,
{
    fn hash_iter<I>(&self, input: I) -> [PF; OUT]
    where
        I: IntoIterator<Item = F>,
    {
        let mut state = [PF::default(); WIDTH];
        let mut input_iter = input.into_iter();
        let mut temp_chunk = Vec::with_capacity(self.num_f_elms);

        'outer: loop {
            for chunk_id in 0..RATE {
                temp_chunk.clear();
                for _ in 0..self.num_f_elms {
                    if let Some(item) = input_iter.next() {
                        temp_chunk.push(item);
                    } else {
                        if chunk_id == 0 && temp_chunk.is_empty() {
                            break 'outer;
                        }
                        break;
                    }
                }
                if !temp_chunk.is_empty() {
                    state[chunk_id] = reduce_32(&temp_chunk);
                }
            }
            state = self.permutation.permute(state);
        }

        state[..OUT].try_into().unwrap()
    }
}