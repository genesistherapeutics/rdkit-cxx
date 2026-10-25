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

// Expected values below come from Python's DistGeom.

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
    assert_eq!(bounds[2][0], 2.4257022902910914);
    assert_eq!(bounds[0][3], 3.8753832408954354);
    let raw = bounds.clone();
    assert!(do_triangle_smoothing(&mut bounds).unwrap());
    assert_eq!(bounds, raw);

    let mut bad = vec![
        vec![0.0, 1.0, 1.0],
        vec![0.5, 0.0, 1.0],
        vec![3.0, 0.5, 0.0],
    ];
    assert!(!do_triangle_smoothing(&mut bad).unwrap());
    assert_eq!(
        bad,
        vec![
            vec![0.0, 1.0, 1.0],
            vec![0.5, 0.0, 1.0],
            vec![3.0, 2.0, 0.0]
        ]
    );
    let mut close = vec![
        vec![0.0, 1.0, 1.0],
        vec![0.5, 0.0, 1.0],
        vec![2.02, 0.5, 0.0],
    ];
    assert!(!do_triangle_smoothing_with_tol(&mut close, 0.05).unwrap());
    assert_eq!(
        close,
        vec![
            vec![0.0, 1.0, 1.0],
            vec![0.5, 0.0, 1.02],
            vec![2.02, 1.02, 0.0]
        ]
    );
    let value = |result: Result<bool, Error>| matches!(result, Err(Error::Value { .. }));
    assert!(value(do_triangle_smoothing(&mut [])));
    assert!(value(do_triangle_smoothing(&mut [vec![0.0, 1.0]])));
}

#[test]
fn embedding() {
    let expected = [
        (
            -1.5209907982334718,
            -0.5716127504791275,
            -0.15517330685885306,
        ),
        (-0.6739436193085492, 0.5639549369723209, 0.3497242904897715),
        (0.6660322905746084, 0.5808228503019247, -0.3443470543822283),
        (1.528902126967413, -0.5731650367951179, 0.149796070751309),
    ];
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
    for params in others {
        let coords = embed_bounds_matrix_with_params(&butane_bounds(), &params).unwrap();
        assert_eq!(coords.len(), 4);
        for (point, e) in coords.iter().zip(expected) {
            assert!((point.x - e.0).abs() < 1e-9);
            assert!((point.y - e.1).abs() < 1e-9);
            assert!((point.z - e.2).abs() < 1e-9);
        }
    }
    assert_eq!(embed_bounds_matrix(&butane_bounds()).unwrap().len(), 4);
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
