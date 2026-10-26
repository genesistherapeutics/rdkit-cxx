// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_dist_geom::{
    GetMoleculeBoundsMatrixParams, get_molecule_bounds_matrix_with_params,
};
use rdkit_cxx::distance_geometry::*;

// Expected values below come from Python's DistGeom. Floats are compared
// with a tolerance, since they differ in the last bits between platforms,
// and embedded coordinates, which differ more, by the bounds they meet.

fn assert_matrix_close(actual: &[Vec<f64>], expected: &[Vec<f64>]) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert_eq!(a.len(), e.len());
        for (x, y) in a.iter().zip(e) {
            assert!((x - y).abs() < 1e-9, "{actual:?} != {expected:?}");
        }
    }
}

/// Asserts that the points meet the bounds, to a tolerance.
fn assert_within_bounds(points: &[rdkit_cxx::geometry::Point3D], bounds: &[Vec<f64>]) {
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            let d = points[i].distance(&points[j]);
            assert!(
                bounds[j][i] - 0.1 < d && d < bounds[i][j] + 0.1,
                "{i}-{j}: {d}"
            );
        }
    }
}

fn butane_bounds() -> Vec<Vec<f64>> {
    get_molecule_bounds_matrix_with_params(
        &mol_from_smiles("CCCC").unwrap(),
        &GetMoleculeBoundsMatrixParams {
            do_triangle_smoothing: false,
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn defaults() {
    assert_eq!(
        EmbedBoundsMatrixParams::default(),
        EmbedBoundsMatrixParams {
            max_iters: 10,
            randomize_on_failure: false,
            num_zero_fail: 2,
            weights: Vec::new(),
            random_seed: -1,
        }
    );
}

#[test]
fn triangle_smoothing() {
    let mut bounds = butane_bounds();
    assert!((bounds[2][0] - 2.432790629029564).abs() < 1e-9);
    assert!((bounds[0][3] - 3.810719484563623).abs() < 1e-9);
    let raw = bounds.clone();
    assert!(do_triangle_smoothing(&mut bounds).unwrap());
    assert_matrix_close(&bounds, &raw);

    let mut bad = vec![
        vec![0.0, 1.0, 1.0],
        vec![0.5, 0.0, 1.0],
        vec![3.0, 0.5, 0.0],
    ];
    assert!(!do_triangle_smoothing(&mut bad).unwrap());
    assert_matrix_close(
        &bad,
        &[
            vec![0.0, 1.0, 1.0],
            vec![0.5, 0.0, 1.0],
            vec![3.0, 2.0, 0.0],
        ],
    );
    let mut close = vec![
        vec![0.0, 1.0, 1.0],
        vec![0.5, 0.0, 1.0],
        vec![2.02, 0.5, 0.0],
    ];
    assert!(!do_triangle_smoothing_with_tol(&mut close, 0.05).unwrap());
    assert_matrix_close(
        &close,
        &[
            vec![0.0, 1.0, 1.0],
            vec![0.5, 0.0, 1.02],
            vec![2.02, 1.02, 0.0],
        ],
    );
    let value = |result: Result<bool, Error>| matches!(result, Err(Error::Value { .. }));
    assert!(value(do_triangle_smoothing(&mut [])));
    assert!(value(do_triangle_smoothing(&mut [vec![0.0, 1.0]])));
}

#[test]
fn embedding() {
    let seeded = EmbedBoundsMatrixParams {
        random_seed: 42,
        ..Default::default()
    };
    let others = [
        seeded.clone(),
        EmbedBoundsMatrixParams {
            weights: vec![(0, 3, 10.0)],
            ..seeded.clone()
        },
        EmbedBoundsMatrixParams {
            max_iters: 3,
            randomize_on_failure: true,
            num_zero_fail: 1,
            ..seeded.clone()
        },
    ];
    let bounds = butane_bounds();
    for params in others {
        let coords = embed_bounds_matrix_with_params(&bounds, &params).unwrap();
        assert_eq!(coords.len(), 4);
        assert_within_bounds(&coords, &bounds);
    }
    // A seed repeats.
    let first = embed_bounds_matrix_with_params(&bounds, &seeded).unwrap();
    let again = embed_bounds_matrix_with_params(&bounds, &seeded).unwrap();
    for (p, q) in first.iter().zip(&again) {
        assert!(p.distance(q) < 1e-12);
    }
    let unseeded = embed_bounds_matrix(&bounds).unwrap();
    assert_eq!(unseeded.len(), 4);
    assert_within_bounds(&unseeded, &bounds);
    // RDKit adds iter * 999 to the seed in int, which overflows here.
    let overflowing = EmbedBoundsMatrixParams {
        random_seed: i32::MAX - 500,
        max_iters: 3,
        ..Default::default()
    };
    let _ = embed_bounds_matrix_with_params(&butane_bounds(), &overflowing);
}

#[test]
fn embedding_inputs_rdkit_trusts() {
    // A lower bound above its upper bound fails RDKit's invariant.
    assert!(matches!(
        embed_bounds_matrix(&[vec![0.0, 1.0], vec![2.0, 0.0]]),
        Err(Error::Runtime { .. })
    ));
    let value = |result: rdkit_cxx::Result<_>| matches!(result, Err(Error::Value { .. }));
    assert!(value(embed_bounds_matrix(&[])));
    assert!(value(embed_bounds_matrix(&[vec![0.0, 1.0]])));
    assert!(value(embed_bounds_matrix(&[
        vec![0.0, f64::INFINITY],
        vec![1.0, 0.0]
    ])));
    // RDKit casts the summed squares of the distances it picks to int.
    let n = 12;
    let long: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i < j { 1e6 } else { 0.0 }).collect())
        .collect();
    assert!(value(embed_bounds_matrix(&long)));
}
