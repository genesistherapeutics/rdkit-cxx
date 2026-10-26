// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Python's `random.Random`, the Mersenne Twister RDKit's pure-Python
//! modules draw from, so that their ports draw the same numbers from the
//! same seed (CONVENTIONS.md R2).
//!
//! Python's modules use the `random` module's global generator or a
//! `random.Random` they make; their ports take a [`Random`] instead.

use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hasher};

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

/// A Mersenne Twister (MT19937) seeded as Python seeds one, which draws the
/// numbers Python's `random.Random` draws. Mirrors Python's
/// `random.Random`.
#[derive(Clone)]
pub struct Random {
    mt: Box<[u32; N]>,
    index: usize,
}

impl std::fmt::Debug for Random {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Random").finish_non_exhaustive()
    }
}

impl Random {
    /// A generator seeded with `seed`, as `random.Random(seed)` seeds one
    /// for a nonnegative `int` (and for a negative one, with its absolute
    /// value).
    pub fn new(seed: u64) -> Random {
        // Python splits the seed into 32-bit words, least significant
        // first, and seeds with at least one.
        let mut key = vec![seed as u32];
        if seed >> 32 != 0 {
            key.push((seed >> 32) as u32);
        }
        Random::from_key(&key)
    }

    /// A generator seeded from the operating system's randomness, as
    /// `random.Random()` is.
    pub fn from_entropy() -> Random {
        let mut key = Vec::with_capacity(8);
        for i in 0..4u64 {
            let mut hasher = RandomState::new().build_hasher();
            hasher.write_u64(i);
            let x = hasher.finish();
            key.push(x as u32);
            key.push((x >> 32) as u32);
        }
        Random::from_key(&key)
    }

    /// `init_by_array` of the reference implementation, which Python's
    /// `random_seed` calls with the seed's 32-bit words.
    fn from_key(key: &[u32]) -> Random {
        let mut r = Random::from_u32(19_650_218);
        let mt = &mut r.mt;
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..N.max(key.len()) {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..N - 1 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_566_083_941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
        }
        mt[0] = 0x8000_0000;
        r
    }

    /// `init_genrand` of the reference implementation.
    fn from_u32(s: u32) -> Random {
        let mut mt = Box::new([0u32; N]);
        mt[0] = s;
        for i in 1..N {
            mt[i] = 1_812_433_253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Random { mt, index: N }
    }

    /// The next 32 random bits, `genrand_uint32` of the reference
    /// implementation.
    fn next_u32(&mut self) -> u32 {
        if self.index >= N {
            let mt = &mut self.mt;
            for kk in 0..N {
                let y = (mt[kk] & UPPER_MASK) | (mt[(kk + 1) % N] & LOWER_MASK);
                let mag = if y & 1 == 0 { 0 } else { MATRIX_A };
                mt[kk] = mt[(kk + M) % N] ^ (y >> 1) ^ mag;
            }
            self.index = 0;
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// A float in [0, 1) with 53 random bits. Mirrors `Random.random`.
    pub fn random(&mut self) -> f64 {
        let a = self.next_u32() >> 5;
        let b = self.next_u32() >> 6;
        (f64::from(a) * 67_108_864.0 + f64::from(b)) * (1.0 / 9_007_199_254_740_992.0)
    }

    /// An integer of `k` random bits, as 32-bit words, least significant
    /// first: one word for `k` up to 32 (none for 0). Mirrors
    /// `Random.getrandbits(k)`.
    pub fn getrandbits(&mut self, k: u32) -> Vec<u32> {
        if k == 0 {
            return Vec::new();
        }
        let words = (k - 1) / 32 + 1;
        let mut left = k;
        (0..words)
            .map(|_| {
                let mut r = self.next_u32();
                if left < 32 {
                    r >>= 32 - left;
                }
                left = left.saturating_sub(32);
                r
            })
            .collect()
    }

    /// A random integer in `[0, n)`, for `n` below 2^64. Mirrors Python's
    /// `Random._randbelow`, which `randrange`, `choice`, `shuffle`, and
    /// `sample` use.
    pub fn randbelow(&mut self, n: u64) -> u64 {
        let k = u64::BITS - n.leading_zeros();
        loop {
            let words = self.getrandbits(k);
            let r = words
                .iter()
                .rev()
                .fold(0u64, |acc, &w| (acc << 32) | u64::from(w));
            if r < n {
                return r;
            }
        }
    }

    /// `k` distinct items of `population`, in the order Python's
    /// `Random.sample(population, k)` picks them. `k` larger than the
    /// population panics, where Python raises `ValueError`.
    pub fn sample<T: Clone>(&mut self, population: &[T], k: usize) -> Vec<T> {
        let n = population.len();
        assert!(k <= n, "Sample larger than population");
        let mut setsize = 21usize;
        if k > 5 {
            // Python computes ceil(log(3k, 4)) in floating point.
            let exp = ((k as f64 * 3.0).ln() / 4f64.ln()).ceil() as u32;
            setsize += 4usize.pow(exp);
        }
        let mut result = Vec::with_capacity(k);
        if n <= setsize {
            let mut pool: Vec<T> = population.to_vec();
            for i in 0..k {
                let j = self.randbelow((n - i) as u64) as usize;
                result.push(pool[j].clone());
                pool[j] = pool[n - i - 1].clone();
            }
        } else {
            let mut selected = std::collections::HashSet::new();
            for _ in 0..k {
                let mut j = self.randbelow(n as u64) as usize;
                while selected.contains(&j) {
                    j = self.randbelow(n as u64) as usize;
                }
                selected.insert(j);
                result.push(population[j].clone());
            }
        }
        result
    }

    /// Shuffles `x` in place as `rdkit.RDRandom.shuffle(x,
    /// random=random.random)` does: Python 3.10's `random.shuffle` with a
    /// float function, which RDKit's modules call and Python 3.11 removed.
    pub fn shuffle_with_random<T>(&mut self, x: &mut [T]) {
        for i in (1..x.len()).rev() {
            let j = (self.random() * (i + 1) as f64).floor() as usize;
            x.swap(i, j);
        }
    }
}

/// Python's `hash` of a tuple of small nonnegative integers and of tuples of
/// them, as 3.8 and later compute it (xxHash's mixing) on 64-bit platforms:
/// what `EnumerateStereoisomers` seeds its default generator with.
pub(crate) fn hash_tuple(items: &[u64]) -> i64 {
    const P1: u64 = 11_400_714_785_074_694_791;
    const P2: u64 = 14_029_467_366_897_019_727;
    const P5: u64 = 2_870_177_450_012_600_261;
    let mut acc = P5;
    for &lane in items {
        acc = acc.wrapping_add(lane.wrapping_mul(P2));
        acc = acc.rotate_left(31);
        acc = acc.wrapping_mul(P1);
    }
    acc = acc.wrapping_add((items.len() as u64) ^ (P5 ^ 3_527_539));
    if acc == u64::MAX {
        return 1_546_275_796;
    }
    acc as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    // Values from Python 3.12's random module.
    #[test]
    fn matches_python() {
        let mut r = Random::new(23);
        assert_eq!(r.random(), 0.9248652516259452);
        assert_eq!(r.random(), 0.9486057779931771);
        let mut r = Random::new(0xf00d);
        assert_eq!(r.getrandbits(20), [485_582]);
        assert_eq!(r.getrandbits(40), [3_965_503_819, 32]);
        let mut r = Random::new(1 << 40);
        assert_eq!(r.random(), 0.1036394562812375);
        let mut r = Random::new(0);
        assert_eq!(r.random(), 0.8444218515250481);
        let mut r = Random::new(42);
        assert_eq!(
            r.sample(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 10),
            [1, 0, 4, 9, 6, 5, 8, 2, 3, 7]
        );
        let mut r = Random::new(42);
        let mut x: Vec<u32> = (0..10).collect();
        r.shuffle_with_random(&mut x);
        assert_eq!(x, [9, 7, 8, 5, 3, 4, 1, 2, 0, 6]);
        assert_eq!(
            hash_tuple(&[hash_tuple(&[1, 6]) as u64, hash_tuple(&[2, 8]) as u64]),
            7_099_509_676_313_867_053
        );
        // Thirty of a hundred, from Python's pool; five of a thousand, by sets.
        let mut r = Random::new(7);
        let population: Vec<u32> = (0..100).collect();
        assert_eq!(
            r.sample(&population, 30),
            [
                41, 19, 50, 83, 6, 9, 68, 12, 46, 74, 7, 64, 27, 4, 11, 55, 53, 8, 30, 85, 70, 54,
                89, 72, 15, 28, 77, 97, 95, 90
            ]
        );
        let population: Vec<u32> = (0..1000).collect();
        assert_eq!(
            Random::new(7).sample(&population, 5),
            [331, 970, 154, 404, 666]
        );
    }
}
