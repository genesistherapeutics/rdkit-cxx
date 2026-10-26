// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Checks the port of rdkit.RDRandom against Python's random module in
// Python 3.12, which RDKit 2026.09.1's RDRandom draws from.

use rdkit_cxx::Error;
use rdkit_cxx::rd_random::{self, Random};

#[test]
fn generator_matches_python() {
    let mut rng = Random::new(0);
    assert_eq!(rng.random(), 0.8444218515250481);
    assert_eq!(rng.random(), 0.7579544029403025);
    assert_eq!(Random::new(-3).random(), 0.23796462709189137);
    assert_eq!(Random::new(-3).random(), Random::new(3).random());
    assert_eq!(Random::new((1 << 40) + 7).random(), 0.6137037779936511);
    rng.seed(0);
    assert_eq!(rng.random(), 0.8444218515250481);
}

#[test]
fn shuffles_match_python() {
    // `random.shuffle(x)`: RDRandom.shuffle(x) itself raises AttributeError
    // in RDKit 2026.09.1, reading a `random._randbelow` that Python 3.12
    // lacks.
    let mut x: Vec<u32> = (0..12).collect();
    Random::new(5).shuffle(&mut x);
    assert_eq!(x, [10, 2, 11, 7, 1, 3, 6, 0, 8, 5, 4, 9]);
    // `RDRandom.shuffle(x, random=RDRandom.random)`.
    let mut rng = Random::new(5);
    let mut x: Vec<u32> = (0..12).collect();
    rd_random::shuffle_with_random(&mut x, || rng.random()).unwrap();
    assert_eq!(x, [9, 4, 1, 3, 2, 0, 6, 5, 10, 11, 8, 7]);
    let mut x = [1, 2, 3];
    assert!(matches!(
        rd_random::shuffle_with_random(&mut x, || 1.5),
        Err(Error::Value { .. })
    ));
}

#[test]
fn shared_generator() {
    // The only test in this binary to draw from the shared generator.
    rd_random::seed(5);
    let mut x: Vec<u32> = (0..12).collect();
    rd_random::shuffle(&mut x);
    assert_eq!(x, [10, 2, 11, 7, 1, 3, 6, 0, 8, 5, 4, 9]);
    rd_random::seed(0);
    assert_eq!(rd_random::random(), 0.8444218515250481);
}

#[test]
fn draws_for_the_chem_ports_match_python() {
    // getrandbits, sample, and RDRandom's float shuffle, as BRICS,
    // EnumerateStereoisomers, and Randomize draw them.
    let mut r = Random::new(23);
    assert_eq!(r.random(), 0.9248652516259452);
    assert_eq!(r.random(), 0.9486057779931771);
    let mut r = Random::new(0xf00d);
    assert_eq!(r.getrandbits_words(20), [485_582]);
    assert_eq!(r.getrandbits_words(40), [3_965_503_819, 32]);
    let mut r = Random::new(0xf00d);
    assert_eq!(r.getrandbits(20), 485_582);
    assert_eq!(r.getrandbits(40), (32 << 32) | 3_965_503_819);
    assert!(Random::new(1).getrandbits_words(0).is_empty());
    assert_eq!(Random::new(1 << 40).random(), 0.1036394562812375);
    let mut r = Random::new(42);
    assert_eq!(
        r.sample(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 10),
        [1, 0, 4, 9, 6, 5, 8, 2, 3, 7]
    );
    let mut r = Random::new(42);
    let mut x: Vec<u32> = (0..10).collect();
    rd_random::shuffle_with_random(&mut x, || r.random()).unwrap();
    assert_eq!(x, [9, 7, 8, 5, 3, 4, 1, 2, 0, 6]);
    // Thirty of a hundred, from Python's pool; five of a thousand, by sets.
    let mut r = Random::new(7);
    let population: Vec<u32> = (0..100).collect();
    assert_eq!(
        r.sample(&population, 30),
        [
            41, 19, 50, 83, 6, 9, 68, 12, 46, 74, 7, 64, 27, 4, 11, 55, 53, 8, 30, 85, 70, 54, 89,
            72, 15, 28, 77, 97, 95, 90
        ]
    );
    let population: Vec<u32> = (0..1000).collect();
    assert_eq!(
        Random::new(7).sample(&population, 5),
        [331, 970, 154, 404, 666]
    );
}
