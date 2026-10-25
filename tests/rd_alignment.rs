// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::geometry::Point3D;
use rdkit_cxx::numerics::rd_alignment::*;

// Expected values below come from Python's rdAlignment.

fn points(coords: &[(f64, f64, f64)]) -> Vec<Point3D> {
    coords
        .iter()
        .map(|&(x, y, z)| Point3D::new(x, y, z))
        .collect()
}

fn assert_rows(actual: &[Vec<f64>], expected: [[f64; 4]; 4]) {
    for (row, e) in actual.iter().zip(expected) {
        for (a, b) in row.iter().zip(e) {
            assert!((a - b).abs() < 1e-12, "{actual:?}");
        }
    }
}

#[test]
fn alignment_transforms() {
    assert_eq!(
        GetAlignmentTransformParams::default(),
        GetAlignmentTransformParams {
            weights: None,
            reflect: false,
            max_iterations: 50,
        }
    );
    let reference = points(&[
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.5),
    ]);
    let probe = points(&[
        (1.0, 1.0, 1.0),
        (1.0, 2.0, 1.0),
        (0.0, 1.0, 1.0),
        (1.2, 1.0, 2.4),
    ]);
    let (ssd, trans) = get_alignment_transform(&reference, &probe).unwrap();
    assert!((ssd - 0.00991970967898581).abs() < 1e-12);
    assert_rows(
        &trans,
        [
            [
                -0.01367855473628183,
                0.9996035741669639,
                0.02460877187016158,
                -1.0217834659443974,
            ],
            [
                -0.9943985691024161,
                -0.016178727140344136,
                0.10444967474906626,
                0.9247352032961236,
            ],
            [
                0.10480640680508982,
                -0.023042206941814535,
                0.9942256654260424,
                -1.022247015091961,
            ],
            [0.0, 0.0, 0.0, 1.0],
        ],
    );
    let params = GetAlignmentTransformParams {
        weights: Some(vec![1.0, 2.0, 3.0, 4.0]),
        reflect: true,
        max_iterations: 10,
    };
    let (ssd, trans) = get_alignment_transform_with_params(&reference, &probe, &params).unwrap();
    assert!((ssd - 1.3731858933562364).abs() < 1e-12);
    assert_rows(
        &trans,
        [
            [
                0.7793447160713263,
                0.1441923796633452,
                -0.609778952718709,
                0.3703354321095372,
            ],
            [
                -0.23315997340929243,
                -0.8365449276609681,
                -0.4958114669957561,
                2.2591845809657896,
            ],
            [
                0.5815997251817455,
                -0.5285840913721557,
                0.6183369777207172,
                -0.18395256123949338,
            ],
            [0.0, 0.0, 0.0, 1.0],
        ],
    );

    let value = |result: rdkit_cxx::Result<_>| matches!(result, Err(Error::Value { .. }));
    assert!(value(get_alignment_transform(&[], &[])));
    assert!(value(get_alignment_transform(&reference, &probe[..3])));
    let one = GetAlignmentTransformParams {
        weights: Some(vec![1.0]),
        ..Default::default()
    };
    assert!(value(get_alignment_transform_with_params(
        &reference, &probe, &one
    )));
    let negative = GetAlignmentTransformParams {
        weights: Some(vec![1.0, -1.0, 1.0, 1.0]),
        ..Default::default()
    };
    assert!(matches!(
        get_alignment_transform_with_params(&reference, &probe, &negative),
        Err(Error::Runtime { .. })
    ));
}
