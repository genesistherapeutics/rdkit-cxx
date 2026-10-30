// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{assert_close, ethanol, ethanol_two_confs};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_shape_helpers::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::data_structs::DiscreteValueType;
use rdkit_cxx::geometry::{Point3D, UniformGrid3D, UniformGrid3DParams};

// Expected values below come from Python's rdShapeHelpers on the same
// molecules. The comparisons encode both shapes in the first conformer's
// principal axes, whose signs can differ between platforms, so their
// values are compared to 0.02.

fn assert_point(actual: Point3D, expected: (f64, f64, f64)) {
    assert_close(actual.x, expected.0);
    assert_close(actual.y, expected.1);
    assert_close(actual.z, expected.2);
}

fn translation() -> Vec<Vec<f64>> {
    let mut trans = vec![vec![0.0; 4]; 4];
    for (i, row) in trans.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    trans[0][3] = 1.0;
    trans[1][3] = -2.0;
    trans
}

fn grid() -> UniformGrid3D {
    UniformGrid3D::new_with_params(
        10.0,
        10.0,
        10.0,
        &UniformGrid3DParams {
            off_set: Some(Point3D::new(-5.0, -5.0, -5.0)),
            ..Default::default()
        },
    )
    .unwrap()
}

/// Asserts a grid's total value, to 1% for points on a sphere's edge that
/// another platform's rounding may place differently.
fn assert_total(grid: &UniformGrid3D, expected: u32) {
    let total = grid.get_occupancy_vect().get_total_val();
    assert!(
        total.abs_diff(expected) * 100 <= expected,
        "{total} != {expected}"
    );
}

fn value<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

#[test]
fn params_defaults() {
    assert_eq!(
        EncodeShapeParams::default(),
        EncodeShapeParams {
            conf_id: -1,
            trans: None,
            vdw_scale: 0.8,
            step_size: 0.25,
            max_layers: -1,
            ignore_hs: true,
        }
    );
    assert_eq!(
        ShapeCompareParams::default(),
        ShapeCompareParams {
            conf_id1: -1,
            conf_id2: -1,
            grid_spacing: 0.5,
            bits_per_point: DiscreteValueType::Twobitvalue,
            vdw_scale: 0.8,
            step_size: 0.25,
            max_layers: -1,
            ignore_hs: true,
        }
    );
    assert!(ShapeProtrudeDistParams::default().allow_reordering);
    assert_eq!(
        ComputeConfBoxParams::default(),
        ComputeConfBoxParams {
            trans: None,
            padding: 2.0,
        }
    );
}

#[test]
fn boxes() {
    let mol = ethanol_two_confs();
    let conf = mol.get_conformer().unwrap();
    let (lower, upper) = compute_conf_box(&conf).unwrap();
    assert_point(lower, (-3.5833, -3.4177, -3.23));
    assert_point(upper, (3.8339, 3.032, 3.0296));
    let shifted = ComputeConfBoxParams {
        trans: Some(translation()),
        padding: 0.0,
    };
    let (lower, upper) = compute_conf_box_with_params(&conf, &shifted).unwrap();
    assert_point(lower, (-0.5833, -3.4177, -1.23));
    assert_point(upper, (2.8339, -0.968, 1.0296));

    let (dims, offset) = compute_conf_dims_and_offset(&conf).unwrap();
    assert_point(dims, (7.4172, 6.4497, 6.2596));
    assert_point(offset, (-3.5833, -3.4177, -3.23));
    let padded = ComputeConfBoxParams {
        trans: Some(translation()),
        padding: 1.0,
    };
    let (dims, offset) = compute_conf_dims_and_offset_with_params(&conf, &padded).unwrap();
    assert_point(dims, (5.4172, 4.4497, 4.2596));
    assert_point(offset, (-1.5833, -4.4177, -2.23));

    let box1 = compute_conf_box(&conf).unwrap();
    let box2 = compute_conf_box_with_params(
        &mol.get_conformer_with_id(1).unwrap(),
        &ComputeConfBoxParams {
            trans: Some(translation()),
            ..Default::default()
        },
    )
    .unwrap();
    let (lower, upper) = compute_union_box(box1, box2);
    assert_point(lower, (-3.5833, -5.4966, -3.23));
    assert_point(upper, (4.8331, 3.032, 3.0296));

    // A conformer outside any molecule has positions enough.
    let free = Conformer::from(conf);
    assert_eq!(compute_conf_box(&free).unwrap(), box1);
    let bad = ComputeConfBoxParams {
        trans: Some(vec![vec![1.0; 3]; 3]),
        ..Default::default()
    };
    assert!(value(compute_conf_box_with_params(&free, &bad)));
}

#[test]
fn encoding() {
    let mol = ethanol_two_confs();
    let mut g = grid();
    encode_shape(&mol, &mut g).unwrap();
    assert_total(&g, 917);
    let origin = g.get_grid_point_index(&Point3D::new(0.0, 0.0, 0.0));
    assert_eq!(g.get_val(origin as u32).unwrap(), 3);
    let mut g = grid();
    let params = EncodeShapeParams {
        conf_id: 1,
        trans: Some(translation()),
        vdw_scale: 1.0,
        step_size: 0.5,
        max_layers: 1,
        ignore_hs: false,
    };
    encode_shape_with_params(&mol, &mut g, &params).unwrap();
    assert_total(&g, 1290);

    assert!(value(encode_shape_with_params(
        &mol,
        &mut grid(),
        &EncodeShapeParams {
            conf_id: 2,
            ..Default::default()
        }
    )));
}

#[test]
fn encoding_inputs_rdkit_trusts() {
    let mol = ethanol();
    // RDKit casts each center in grid units to int; centers that far away
    // are outside the grid, which it skips.
    let mut far = translation();
    far[0][3] = 1e300;
    let mut g = grid();
    let params = EncodeShapeParams {
        trans: Some(far),
        ..Default::default()
    };
    encode_shape_with_params(&mol, &mut g, &params).unwrap();
    assert_eq!(g.get_occupancy_vect().get_total_val(), 0);
    // And each sphere's extent; a negative step casts a negative layer to
    // unsigned int. Python gets 0 and 57 from those casts here. The grid is
    // unchanged.
    for params in [
        EncodeShapeParams {
            vdw_scale: 1e12,
            ..Default::default()
        },
        EncodeShapeParams {
            step_size: -0.25,
            ..Default::default()
        },
    ] {
        let mut g = grid();
        assert!(value(encode_shape_with_params(&mol, &mut g, &params)));
        assert_eq!(g.get_occupancy_vect().get_total_val(), 0);
    }
}

#[test]
fn comparisons() {
    let mol = ethanol_two_confs();
    let other = ShapeCompareParams {
        conf_id1: 0,
        conf_id2: 1,
        ..Default::default()
    };
    let near = |actual: f64, expected: f64| {
        assert!((actual - expected).abs() < 0.02, "{actual} != {expected}");
    };
    assert_eq!(shape_tanimoto_dist(&mol, &mol).unwrap(), 0.0);
    near(
        shape_tanimoto_dist_with_params(&mol, &mol, &other).unwrap(),
        0.1510204081632653,
    );
    let changed = ShapeCompareParams {
        grid_spacing: 0.3,
        bits_per_point: DiscreteValueType::Onebitvalue,
        vdw_scale: 1.0,
        step_size: 0.5,
        max_layers: 1,
        ignore_hs: false,
        ..other.clone()
    };
    near(
        shape_tanimoto_dist_with_params(&mol, &mol, &changed).unwrap(),
        0.23928258967629046,
    );
    assert_eq!(shape_tversky_index(&mol, &mol, 1.0, 1.0).unwrap(), 1.0);
    near(
        shape_tversky_index_with_params(&mol, &mol, 0.7, 0.3, &other).unwrap(),
        0.9154929577464789,
    );
    assert_eq!(shape_protrude_dist(&mol, &mol).unwrap(), 0.0);
    let protrude = ShapeProtrudeDistParams {
        compare: other.clone(),
        ..Default::default()
    };
    let expected = 0.07452725250278086;
    near(
        shape_protrude_dist_with_params(&mol, &mol, &protrude).unwrap(),
        expected,
    );
    // RDKit's protrude distance always uses a spacing of 0.5.
    let spaced = ShapeProtrudeDistParams {
        compare: ShapeCompareParams {
            grid_spacing: 0.3,
            ..other.clone()
        },
        ..Default::default()
    };
    assert_eq!(
        shape_protrude_dist_with_params(&mol, &mol, &spaced).unwrap(),
        shape_protrude_dist_with_params(&mol, &mol, &protrude).unwrap()
    );
    let unordered = ShapeProtrudeDistParams {
        compare: ShapeCompareParams {
            conf_id1: 1,
            conf_id2: 0,
            ..Default::default()
        },
        allow_reordering: false,
    };
    near(
        shape_protrude_dist_with_params(&mol, &mol, &unordered).unwrap(),
        0.0784313725490196,
    );
}

#[test]
fn comparison_inputs_rdkit_trusts() {
    let mol = ethanol();
    // initGrid casts the grid's size in points to unsigned int; Python
    // raises an IndexError for a negative index from the overflow.
    let fine = ShapeCompareParams {
        grid_spacing: 1e-6,
        ..Default::default()
    };
    assert!(value(shape_tanimoto_dist_with_params(&mol, &mol, &fine)));
    // setSphereOccupancy casts each sphere's extent to int.
    for params in [
        ShapeCompareParams {
            vdw_scale: 1e12,
            ..Default::default()
        },
        ShapeCompareParams {
            step_size: -1.0,
            ..Default::default()
        },
    ] {
        assert!(value(shape_tversky_index_with_params(
            &mol, &mol, 1.0, 1.0, &params
        )));
    }
    let protrude = ShapeProtrudeDistParams {
        compare: ShapeCompareParams {
            vdw_scale: 1e12,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(value(shape_protrude_dist_with_params(
        &mol, &mol, &protrude
    )));
    // No conformer.
    let flat = mol_from_smiles("CCO").unwrap();
    assert!(value(shape_tanimoto_dist(&flat, &mol)));
}

/// The comparisons orient both shapes along the first's principal axes,
/// which Eigen's eigensolver finds, so the coordinates must be finite.
#[test]
fn comparisons_need_finite_coordinates() {
    let good = ethanol();
    let mut bad = ethanol();
    bad.get_conformer_mut()
        .unwrap()
        .set_atom_position(
            0,
            Point3D {
                x: f64::NAN,
                y: 0.0,
                z: 0.0,
            },
        )
        .unwrap();
    assert!(value(shape_tanimoto_dist(&bad, &good)));
    assert!(value(shape_tanimoto_dist(&good, &bad)));
    assert_close(shape_tanimoto_dist(&good, &good).unwrap(), 0.0);
}
