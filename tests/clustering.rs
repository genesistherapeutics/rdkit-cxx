// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::ml::cluster::clustering::*;

// Expected values come from Python's rdkit.ML.Cluster.Clustering.

type Merges = (Vec<i64>, Vec<i64>, Vec<f64>);

fn assert_merges(actual: &Merges, ia: &[i64], ib: &[i64], crit: &[f64]) {
    assert_eq!(actual.0, ia, "{actual:?}");
    assert_eq!(actual.1, ib, "{actual:?}");
    assert_eq!(actual.2.len(), crit.len());
    for (a, e) in actual.2.iter().zip(crit) {
        assert!((a - e).abs() <= 1e-12 * e.abs().max(1.0), "{actual:?}");
    }
}

fn points() -> Vec<Vec<f64>> {
    vec![
        vec![0.0, 0.0],
        vec![1.0, 0.5],
        vec![4.0, 4.0],
        vec![5.0, 4.5],
        vec![10.0, 0.0],
        vec![0.5, 9.0],
    ]
}

const DISTS: [f64; 10] = [1.0, 4.0, 3.5, 2.0, 6.0, 5.5, 7.0, 7.5, 2.5, 3.0];

#[test]
fn points_by_criterion() {
    let expected: [(&[i64], &[i64], &[f64]); 7] = [
        (
            &[1, 3, 3, 1, 1, 0],
            &[2, 4, 6, 3, 5, 0],
            &[
                0.625,
                0.625,
                25.708333333333332,
                45.94166666666667,
                62.80833333333333,
                0.0,
            ],
        ),
        (
            &[1, 3, 1, 1, 1, 0],
            &[2, 4, 3, 6, 5, 0],
            &[1.25, 1.25, 21.25, 37.25, 45.25, 0.0],
        ),
        (
            &[1, 3, 3, 1, 1, 0],
            &[2, 4, 6, 3, 5, 0],
            &[1.25, 1.25, 40.5, 81.25, 171.25, 0.0],
        ),
        (
            &[1, 3, 1, 1, 1, 0],
            &[2, 4, 3, 6, 5, 0],
            &[1.25, 1.25, 32.625, 57.875, 89.95, 0.0],
        ),
        (
            &[1, 3, 1, 1, 1, 0],
            &[2, 4, 3, 6, 5, 0],
            &[1.25, 1.25, 32.625, 57.875, 120.4375, 0.0],
        ),
        (
            &[1, 3, 1, 1, 1, 0],
            &[2, 4, 3, 6, 5, 0],
            &[1.25, 1.25, 32.0, 49.5625, 103.890625, 0.0],
        ),
        (
            &[1, 3, 1, 1, 1, 0],
            &[2, 4, 3, 6, 5, 0],
            &[1.25, 1.25, 32.0, 49.5625, 75.37, 0.0],
        ),
    ];
    for (option, (ia, ib, crit)) in (1..=7).zip(expected) {
        assert_merges(
            &murtagh_cluster(&points(), 6, 2, option).unwrap(),
            ia,
            ib,
            crit,
        );
    }
}

#[test]
fn distances_by_criterion() {
    let expected: [(&[i64], &[i64], &[f64]); 7] = [
        (
            &[1, 3, 3, 1, 0],
            &[2, 5, 4, 3, 0],
            &[0.5, 1.25, 2.4166666666666665, 4.233333333333333, 0.0],
        ),
        (
            &[1, 1, 3, 1, 0],
            &[2, 4, 5, 3, 0],
            &[1.0, 2.0, 2.5, 3.0, 0.0],
        ),
        (
            &[1, 3, 3, 1, 0],
            &[2, 5, 4, 3, 0],
            &[1.0, 2.5, 5.5, 7.5, 0.0],
        ),
        (
            &[1, 3, 1, 1, 0],
            &[2, 5, 4, 3, 0],
            &[1.0, 2.5, 4.0, 5.083333333333333, 0.0],
        ),
        (
            &[1, 3, 1, 1, 0],
            &[2, 5, 4, 3, 0],
            &[1.0, 2.5, 4.0, 4.875, 0.0],
        ),
        (
            &[1, 3, 3, 1, 0],
            &[2, 5, 4, 3, 0],
            &[1.0, 2.5, 3.625, 3.28125, 0.0],
        ),
        (
            &[1, 3, 3, 1, 0],
            &[2, 5, 4, 3, 0],
            &[1.0, 2.5, 3.625, 3.527777777777778, 0.0],
        ),
    ];
    for (option, (ia, ib, crit)) in (1..=7).zip(expected) {
        assert_merges(
            &murtagh_dist_cluster(&DISTS, 5, option).unwrap(),
            ia,
            ib,
            crit,
        );
    }
    // As Python, values past n_pts' are ignored.
    assert_merges(
        &murtagh_dist_cluster(&[1.0, 2.0, 3.0, 99.0], 3, 2).unwrap(),
        &[1, 1, 0],
        &[2, 3, 0],
        &[1.0, 2.0, 0.0],
    );
}

#[test]
fn python_errors() {
    // Python's own checks raise RuntimeError.
    assert!(matches!(
        murtagh_cluster(&points(), 5, 2, 1),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        murtagh_cluster(&points(), 6, 3, 1),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        murtagh_dist_cluster(&[], 3, 1),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        murtagh_cluster(&[vec![0.0], vec![1.0, 2.0]], 2, 1, 1),
        Err(Error::Value { .. })
    ));
}

#[test]
fn guards() {
    // Fewer than 2 items: out-of-bounds writes in hc_ (heap corruption in
    // Python).
    for n in [-1, 0, 1] {
        assert!(matches!(
            murtagh_dist_cluster(&[1.0], n, 1),
            Err(Error::Value { .. })
        ));
    }
    assert!(matches!(
        murtagh_cluster(&[vec![0.0, 1.0]], 1, 2, 1),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        murtagh_cluster(&[], 0, 2, 1),
        Err(Error::Value { .. })
    ));
    // A matrix shorter than n_pts' (read past its end in Python).
    assert!(matches!(
        murtagh_dist_cluster(&[1.0, 2.0], 3, 1),
        Err(Error::Value { .. })
    ));
    // No criterion outside 1 to 7.
    for option in [0, 8, -3] {
        assert!(matches!(
            murtagh_dist_cluster(&DISTS, 5, option),
            Err(Error::Value { .. })
        ));
    }
    // Distances that never compare below the sentinel (out-of-bounds merges
    // in Python), and NaN.
    for bad in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        let mut dists = DISTS;
        dists[3] = bad;
        assert!(matches!(
            murtagh_dist_cluster(&dists, 5, 2),
            Err(Error::Value { .. })
        ));
    }
    // Squared distances that overflow.
    assert!(matches!(
        murtagh_cluster(&[vec![0.0], vec![1e200]], 2, 1, 1),
        Err(Error::Value { .. })
    ));
}

#[test]
fn large_distances_are_scaled() {
    // Distances at and past hc_'s 1e20 sentinel merge as the same matrix
    // scaled down by a power of two does, with crit scaled back exactly.
    let scale = 2f64.powi(80);
    for option in 1..=7 {
        let small = murtagh_dist_cluster(&DISTS, 5, option).unwrap();
        let large: Vec<f64> = DISTS.iter().map(|d| d * scale).collect();
        let large = murtagh_dist_cluster(&large, 5, option).unwrap();
        assert_eq!(large.0, small.0);
        assert_eq!(large.1, small.1);
        let scaled: Vec<f64> = small.2.iter().map(|c| c * scale).collect();
        assert_eq!(large.2, scaled);
    }
    // The largest finite distances.
    let huge = vec![f64::MAX, f64::MAX / 2.0, f64::MAX / 4.0];
    let merges = murtagh_dist_cluster(&huge, 3, 3).unwrap();
    assert_eq!(merges.0, [2, 1, 0]);
    assert_eq!(merges.1, [3, 2, 0]);
    assert_eq!(merges.2, [f64::MAX / 4.0, f64::MAX, 0.0]);
}

#[test]
fn threads() {
    // hc_ keeps its state in statics, so calls are serialized.
    let expected = murtagh_cluster(&points(), 6, 2, 1).unwrap();
    std::thread::scope(|scope| {
        for t in 0..8 {
            let expected = &expected;
            scope.spawn(move || {
                for _ in 0..200 {
                    if t % 2 == 0 {
                        assert_eq!(&murtagh_cluster(&points(), 6, 2, 1).unwrap(), expected);
                    } else {
                        let merges = murtagh_dist_cluster(&DISTS, 5, 3).unwrap();
                        assert_eq!(merges.0, [1, 3, 3, 1, 0]);
                    }
                }
            });
        }
    });
}
