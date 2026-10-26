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
