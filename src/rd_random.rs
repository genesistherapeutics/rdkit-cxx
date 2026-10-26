// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.RDRandom`, a Python module: the random numbers RDKit's
//! Python code draws, which come from Python's `random` module.
//!
//! [`Random`] is Python's generator (the Mersenne Twister of
//! `_randommodule.c`, seeded as `random.seed` seeds it from an integer), so
//! the same seed gives the same numbers as Python. Python's module-level
//! functions share one generator; [`seed`], [`random`], and [`shuffle`] do
//! too, behind a lock, seeded from the operating system's randomness until
//! [`seed`] is called. Code that needs a sequence no other thread can draw
//! from holds its own [`Random`].
//!
//! ```no_run
//! use rdkit_cxx::rd_random::Random;
//!
//! let mut rng = Random::new(23);
//! let x = rng.random();
//! ```

use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hasher};
use std::sync::{Mutex, MutexGuard};

use crate::{Error, Result};

const N: usize = 624;
const M: usize = 397;

/// Python's random number generator, the Mersenne Twister MT19937, with
/// `random.Random`'s seeding and floats.
#[derive(Clone)]
pub struct Random {
    state: [u32; N],
    index: usize,
}

impl std::fmt::Debug for Random {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Random").finish_non_exhaustive()
    }
}

impl Random {
    /// A generator seeded as `random.Random(seed)` seeds one. Python seeds
    /// with the absolute value of an integer, so `-5` seeds as `5`.
    pub fn new(seed: i64) -> Random {
        let mut rng = Random {
            state: [0; N],
            index: N,
        };
        rng.seed(seed);
        rng
    }

    /// Reseeds the generator, as `random.seed(seed)` does for an integer.
    pub fn seed(&mut self, seed: i64) {
        let n = seed.unsigned_abs();
        // Python's key is the integer's 32-bit words, least significant
        // first, and one zero word for zero.
        let key: Vec<u32> = if n >> 32 != 0 {
            vec![n as u32, (n >> 32) as u32]
        } else {
            vec![n as u32]
        };
        self.init_by_array(&key);
    }

    fn init_genrand(&mut self, s: u32) {
        self.state[0] = s;
        for i in 1..N {
            let prev = self.state[i - 1];
            self.state[i] = 1812433253u32
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(i as u32);
        }
        self.index = N;
    }

    fn init_by_array(&mut self, key: &[u32]) {
        self.init_genrand(19650218);
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..N.max(key.len()) {
            let prev = self.state[i - 1];
            self.state[i] = (self.state[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1664525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                self.state[0] = self.state[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..N - 1 {
            let prev = self.state[i - 1];
            self.state[i] = (self.state[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1566083941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                self.state[0] = self.state[N - 1];
                i = 1;
            }
        }
        self.state[0] = 0x80000000;
        self.index = N;
    }

    fn genrand_u32(&mut self) -> u32 {
        const MAG01: [u32; 2] = [0, 0x9908b0df];
        const UPPER: u32 = 0x80000000;
        const LOWER: u32 = 0x7fffffff;
        if self.index >= N {
            for kk in 0..N {
                let y = (self.state[kk] & UPPER) | (self.state[(kk + 1) % N] & LOWER);
                self.state[kk] = self.state[(kk + M) % N] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
            }
            self.index = 0;
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^ (y >> 18)
    }

    /// The next float in `[0, 1)`, with 53 random bits, as `random.random()`
    /// returns it.
    pub fn random(&mut self) -> f64 {
        let a = self.genrand_u32() >> 5;
        let b = self.genrand_u32() >> 6;
        (a as f64 * 67108864.0 + b as f64) * (1.0 / 9007199254740992.0)
    }

    /// `k` random bits, as Python's `getrandbits(k)` draws them: 32-bit
    /// words, least significant first, the last cut to the bits left.
    fn getrandbits(&mut self, k: u32) -> u64 {
        let mut r = 0u64;
        let mut left = k;
        let mut shift = 0;
        while left > 0 {
            let mut word = self.genrand_u32();
            if left < 32 {
                word >>= 32 - left;
            }
            r |= u64::from(word) << shift;
            shift += 32;
            left = left.saturating_sub(32);
        }
        r
    }

    /// A random integer below `n`, as Python's `random._randbelow(n)` draws
    /// it (by rejection, from as many bits as `n` has). `n` is positive.
    fn randbelow(&mut self, n: u64) -> u64 {
        let k = 64 - n.leading_zeros();
        loop {
            let r = self.getrandbits(k);
            if r < n {
                return r;
            }
        }
    }

    /// Shuffles `x` in place, as [`shuffle`] does, drawing from this
    /// generator.
    pub fn shuffle<T>(&mut self, x: &mut [T]) {
        for i in (1..x.len()).rev() {
            let j = self.randbelow(i as u64 + 1);
            x.swap(i, j as usize);
        }
    }
}

fn entropy_seed() -> i64 {
    let mut hasher = RandomState::new().build_hasher();
    hasher.write_u64(0x5eed);
    hasher.finish() as i64
}

static GENERATOR: Mutex<Option<Random>> = Mutex::new(None);

/// The shared generator, seeded from the operating system's randomness on
/// first use, as Python's is on import.
fn generator() -> MutexGuard<'static, Option<Random>> {
    let mut guard = GENERATOR.lock().unwrap_or_else(|e| e.into_inner());
    if guard.is_none() {
        *guard = Some(Random::new(entropy_seed()));
    }
    guard
}

fn with_generator<T>(f: impl FnOnce(&mut Random) -> T) -> T {
    f(generator().as_mut().expect("the generator is seeded"))
}

/// Reseeds the shared generator. Mirrors `RDRandom.seed` (Python's
/// `random.seed`) for an integer seed.
pub fn seed(a: i64) {
    with_generator(|rng| rng.seed(a));
}

/// The shared generator's next float in `[0, 1)`. Mirrors `RDRandom.random`
/// (Python's `random.random`).
pub fn random() -> f64 {
    with_generator(Random::random)
}

/// Shuffles `x` in place with the shared generator, swapping each item,
/// from the last, with one at a random index at or before it. Mirrors
/// `RDRandom.shuffle(x)`.
pub fn shuffle<T>(x: &mut [T]) {
    with_generator(|rng| rng.shuffle(x));
}

/// Shuffles `x` in place with the random floats `random` returns, as
/// `RDRandom.shuffle(x, random)` does: item `i`, from the last, swaps with
/// item `floor(random() * (i + 1))`. RDKit's Python code passes
/// `RDRandom.random`, as in `shuffle_with_random(x, rd_random::random)`,
/// which draws differently from [`shuffle`].
///
/// `Error::Value` when `random` returns a float outside `[0, 1)` that
/// makes an index outside `x`; the items may then be partly shuffled. (Python
/// raises `IndexError` there, or wraps a negative index.)
pub fn shuffle_with_random<T>(x: &mut [T], mut random: impl FnMut() -> f64) -> Result<()> {
    for i in (1..x.len()).rev() {
        let j = (random() * (i + 1) as f64).floor();
        if !(0.0..x.len() as f64).contains(&j) {
            return Err(Error::Value {
                message: format!("random() gave the index {j}, outside the sequence"),
            });
        }
        x.swap(i, j as usize);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::Random;

    #[test]
    fn wide_draws_match_python() {
        // random.Random(9)'s _randbelow(2**40 + 3), _randbelow(2**33), and
        // getrandbits(64), in Python 3.12.
        let mut rng = Random::new(9);
        assert_eq!(rng.randbelow((1 << 40) + 3), 585718914800);
        assert_eq!(rng.randbelow(1 << 33), 595022250);
        assert_eq!(rng.getrandbits(64), 12481870810729348949);
    }
}
