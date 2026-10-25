// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{ETHANOL_CONF2, assert_close, assert_close_to, ethanol, ethanol_two_confs};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_gaussian_shape::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Expected values come from Python's rdGaussianShape on the same molecules.
// Alignments optimize, so their scores are compared to 0.02 and their poses
// by invariants; shapes and scores of shapes left in place are compared
// closely.

/// Phenol with hydrogens, as RDKit's embedder placed it with seed 42.
const PHENOL: &str = "
     RDKit          3D

 13 13  0  0  0  0  0  0  0  0999 V2000
    2.4811   -0.3094    0.1041 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0852   -0.1409    0.0818 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2270   -1.2295    0.0635 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1497   -1.0540    0.0416 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6498    0.2281    0.0383 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8016    1.3255    0.0564 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5787    1.1380    0.0783 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0277   -0.4385   -0.7257 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6702   -2.2198    0.0670 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7622   -1.9434    0.0281 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7220    0.3350    0.0210 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2416    2.3272    0.0530 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2570    1.9817    0.0927 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  2  1  0
  1  8  1  0
  3  9  1  0
  4 10  1  0
  5 11  1  0
  6 12  1  0
  7 13  1  0
M  END
";

fn phenol() -> Mol {
    let params = MolFromMolBlockParams {
        remove_hs: false,
        ..Default::default()
    };
    mol_from_mol_block_with_params(PHENOL, &params).unwrap()
}

/// Ethanol in its second conformer, as its only one.
fn ethanol2() -> Mol {
    let mut mol = ethanol();
    let positions: Vec<Point3D> = ETHANOL_CONF2
        .iter()
        .map(|&(x, y, z)| Point3D::new(x, y, z))
        .collect();
    mol.get_conformer_mut()
        .unwrap()
        .set_positions(&positions)
        .unwrap();
    mol
}

fn shape_with(mol: &Mol, shape_opt: ShapeInputOptions) -> rdkit_cxx::Result<ShapeInput> {
    ShapeInput::new_with_params(
        mol,
        &ShapeInputParams {
            shape_opt,
            ..Default::default()
        },
    )
}

fn assert_value_error<T>(result: rdkit_cxx::Result<T>) {
    assert!(
        matches!(result, Err(Error::Value { .. })),
        "expected Error::Value"
    );
}

fn assert_index_error<T>(result: rdkit_cxx::Result<T>) {
    assert!(
        matches!(result, Err(Error::Index { .. })),
        "expected Error::Index"
    );
}

fn assert_triple(actual: (f64, f64, f64), expected: (f64, f64, f64), tol: f64) {
    assert_close_to(actual.0, expected.0, tol);
    assert_close_to(actual.1, expected.1, tol);
    assert_close_to(actual.2, expected.2, tol);
}

fn assert_five(actual: (f64, f64, f64, f64, f64), expected: (f64, f64, f64, f64, f64)) {
    assert_triple(
        (actual.0, actual.1, actual.2),
        (expected.0, expected.1, expected.2),
        1e-6,
    );
    assert_close_to(actual.3, expected.3, 1e-6);
    assert_close_to(actual.4, expected.4, 1e-6);
}

/// Checks that a transformation's 16 values, by rows, are a rotation and a
/// translation.
fn assert_rigid(matrix: &[f64]) {
    assert_eq!(matrix.len(), 16);
    for i in 0..3 {
        for j in 0..3 {
            let dot: f64 = (0..3).map(|k| matrix[4 * i + k] * matrix[4 * j + k]).sum();
            assert_close_to(dot, if i == j { 1.0 } else { 0.0 }, 1e-6);
        }
        assert!(matrix[4 * i + 3].is_finite());
    }
    assert_eq!(&matrix[12..], [0.0, 0.0, 0.0, 1.0]);
}

const ETH_PH: (f64, f64, f64) = (0.3526061582312596, 0.6202436414764871, 0.08496867498603206);
const ETH_ETH2_SCORES: (f64, f64, f64, f64, f64) = (
    0.8999097395930902,
    0.9778488248013629,
    0.8219706543848173,
    98.44286330395616,
    9.592410956014085,
);

#[test]
fn option_defaults() {
    let opts = ShapeInputOptions::default();
    assert!(opts.use_colors);
    assert!(opts.all_carbon_radii);
    assert!(opts.atom_subset.is_empty());
    assert!(opts.custom_features.is_empty());
    assert!(opts.atom_radii.is_empty());
    assert_eq!(opts.shape_prune_threshold, -1.0);
    assert!(opts.sort_shapes);
    assert!(opts.include_dummies);
    let overlay = ShapeOverlayOptions::default();
    assert_eq!(overlay.start_mode, StartMode::ALaPubchem);
    assert_eq!(overlay.optim_mode, OptimMode::ShapePlusColorScore);
    assert_eq!((overlay.sim_alpha, overlay.sim_beta), (1.0, 1.0));
    assert_eq!((overlay.opt_param, overlay.n_steps), (0.5, 100));
    assert!(overlay.normalize);
    assert!(overlay.use_dist_cutoff);
    assert_eq!(overlay.dist_cutoff, 4.5);
    assert_eq!(overlay.shape_convergence_criterion, 0.001);
    let params = ShapeInputParams::default();
    assert_eq!(params.conf_id, -1);
    assert_eq!(params.shape_opt, opts);
    assert_eq!(params.overlay_opts, overlay);
    assert_eq!(
        ShapeToMolParams::default(),
        ShapeToMolParams {
            include_colors: false,
            with_bonds: true
        }
    );
    assert_eq!(BestSimilarityParams::default().threshold, -1.0);
    let align = AlignMolParams::default();
    assert_eq!((align.ref_conf_id, align.fit_conf_id), (-1, -1));
    assert_eq!(align.ref_opts, opts);
    assert_eq!(align.overlay_opts, overlay);
    assert_eq!(AlignMolShapeInputParams::default().fit_conf_id, -1);
    let score = ScoreMolParams::default();
    assert_eq!((score.ref_conf_id, score.fit_conf_id), (-1, -1));
    assert_eq!(score.fit_opts, opts);
    assert_eq!(ScoreMolShapeInputParams::default().fit_conf_id, -1);
    assert_eq!(
        ScoreMoleculeAllConformersParams::default().overlay_opts,
        overlay
    );
    assert_eq!(StartMode::Rotate45Fragment as i32, 6);
    assert_eq!(OptimMode::ShapePlusColor as i32, 2);
}

#[test]
fn makes_shapes() {
    let shape = ShapeInput::new(&ethanol()).unwrap();
    assert_eq!(
        (shape.num_atoms(), shape.num_features(), shape.num_shapes()),
        (3, 2, 1)
    );
    assert_close(shape.shape_volume(), 100.10615354135277);
    assert_close(shape.color_volume(), 10.63121364091602);
    assert_eq!(shape.get_smiles().unwrap(), "[H]OC([H])([H])C([H])([H])[H]");
    assert_eq!(shape.get_active_shape(), 0);

    let shape = ShapeInput::new(&phenol()).unwrap();
    assert_eq!(
        (shape.num_atoms(), shape.num_features(), shape.num_shapes()),
        (7, 3, 1)
    );
    assert_close(shape.shape_volume(), 307.64213726160807);
    assert_close(shape.color_volume(), 15.94682046137403);
    assert_eq!(
        shape.get_smiles().unwrap(),
        "[H]Oc1c([H])c([H])c([H])c([H])c1[H]"
    );

    let mut shape = ShapeInput::new(&ethanol_two_confs()).unwrap();
    assert_eq!(shape.num_shapes(), 2);
    assert_close(shape.shape_volume(), 100.10615354135277);
    shape.set_active_shape(1).unwrap();
    assert_eq!(shape.get_active_shape(), 1);
    assert_close(shape.shape_volume(), 99.00959580299264);
    assert_close(shape.color_volume(), 10.63121364091602);
    assert_index_error(shape.set_active_shape(2));
    let copy = shape.clone();
    assert_eq!(copy.get_active_shape(), 1);
    assert_close(copy.shape_volume(), 99.00959580299264);

    let one = ShapeInput::new_with_params(
        &ethanol_two_confs(),
        &ShapeInputParams {
            conf_id: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(one.num_shapes(), 1);
    assert_close(one.shape_volume(), 99.00959580299264);
}

#[test]
fn shape_options() {
    let ph = phenol();
    let shape = shape_with(
        &ph,
        ShapeInputOptions {
            all_carbon_radii: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_close(shape.shape_volume(), 294.7414891377527);
    let shape = shape_with(
        &ph,
        ShapeInputOptions {
            use_colors: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(shape.num_features(), 0);
    assert_eq!(shape.color_volume(), 0.0);
    let shape = shape_with(
        &ph,
        ShapeInputOptions {
            atom_subset: vec![0, 1, 2, 3],
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!((shape.num_atoms(), shape.num_features()), (4, 1));
    assert_close(shape.shape_volume(), 144.28723045556862);
    assert_eq!(shape.get_smiles().unwrap(), "C:C:CO");

    // Unlike Python, the radii are not left on the molecule.
    let eth = ethanol();
    let shape = shape_with(
        &eth,
        ShapeInputOptions {
            atom_radii: vec![(0, 2.0)],
            all_carbon_radii: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_close(shape.shape_volume(), 98.11902443903753);
    assert!(!eth.get_atom_with_idx(0).unwrap().has_prop("BespokeRadius"));

    let custom = vec![vec![
        (1, Point3D::new(0.0, 0.0, 0.0), 1.0, vec![]),
        (2, Point3D::new(1.0, 0.0, 0.0), 1.0, vec![2]),
    ]];
    let shape = shape_with(
        &eth,
        ShapeInputOptions {
            custom_features: custom,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(shape.num_features(), 2);
    assert_close(shape.color_volume(), 8.377580409597229);

    let pruned = shape_with(
        &ethanol_two_confs(),
        ShapeInputOptions {
            shape_prune_threshold: 0.5,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(pruned.num_shapes(), 1);
}

#[test]
fn shapes_to_molecules() {
    let shape = ShapeInput::new(&phenol()).unwrap();
    let mol = shape.shape_to_mol().unwrap();
    assert_eq!((mol.get_num_atoms(), mol.get_num_bonds()), (7, 7));
    assert_eq!(mol_to_smiles(&mol).unwrap(), "Oc1ccccc1");
    let first = mol.get_conformer().unwrap().get_atom_position(0).unwrap();
    assert_close(first.x, 2.4811);
    let colors = ShapeToMolParams {
        include_colors: true,
        with_bonds: false,
    };
    let mol = shape.shape_to_mol_with_params(&colors).unwrap();
    assert_eq!((mol.get_num_atoms(), mol.get_num_bonds()), (10, 0));

    // Normalizing centers the shape's atoms on the origin.
    let mut shape = shape;
    shape.normalize_coords();
    let positions = shape
        .shape_to_mol()
        .unwrap()
        .get_conformer()
        .unwrap()
        .get_positions();
    for axis in 0..3 {
        let center: f64 = positions.iter().map(|p| [p.x, p.y, p.z][axis]).sum::<f64>() / 7.0;
        assert_close_to(center, 0.0, 1e-9);
    }
    assert_close(shape.shape_volume(), 307.64213726160807);
}

#[test]
fn aligns_molecules() {
    let eth = ethanol();
    let mut fit = ethanol2();
    assert_triple(
        align_mol(&eth, &mut fit).unwrap(),
        (0.9964825027986293, 0.9991791152669989, 0.9937858903302595),
        0.02,
    );
    // The fit's heavy atoms land near the reference's.
    let (aligned, target) = (
        fit.get_conformer().unwrap().get_positions(),
        eth.get_conformer().unwrap().get_positions(),
    );
    for i in 0..3 {
        let (a, b) = (aligned[i], target[i]);
        let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
        assert!(d < 0.3, "atom {i} is {d} away");
    }

    let mut fit = phenol();
    assert_triple(align_mol(&eth, &mut fit).unwrap(), ETH_PH, 0.02);
    let cases = [
        (
            StartMode::ALaPubchem,
            OptimMode::ShapeOnly,
            (0.3493199960520435, 0.6231351738172328, 0.07550481828685415),
        ),
        (
            StartMode::ALaPubchem,
            OptimMode::ShapePlusColor,
            (0.4480012787207643, 0.5266480239857415, 0.36935453345578706),
        ),
        (
            StartMode::Rotate0,
            OptimMode::ShapePlusColorScore,
            (0.3493199960520435, 0.6231351738172328, 0.07550481828685415),
        ),
        (
            StartMode::Rotate180,
            OptimMode::ShapePlusColorScore,
            (0.3493199960520435, 0.6231351738172328, 0.07550481828685415),
        ),
        (
            StartMode::Rotate45Fragment,
            OptimMode::ShapePlusColorScore,
            (0.3878778451542281, 0.6003738630900795, 0.17538182721837667),
        ),
    ];
    for (start_mode, optim_mode, expected) in cases {
        let params = AlignMolParams {
            overlay_opts: ShapeOverlayOptions {
                start_mode,
                optim_mode,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut fit = phenol();
        assert_triple(
            align_mol_with_params(&eth, &mut fit, &params).unwrap(),
            expected,
            0.02,
        );
    }
}

#[test]
fn aligns_to_shapes() {
    let reference = ShapeInput::new(&ethanol()).unwrap();
    let mut fit = phenol();
    let before = fit.get_conformer().unwrap().get_positions();
    assert_triple(
        align_mol_shape_input(&reference, &mut fit).unwrap(),
        ETH_PH,
        0.02,
    );
    // The fit moved rigidly.
    let after = fit.get_conformer().unwrap().get_positions();
    assert_ne!(after, before);
    let dist = |p: &[Point3D], i: usize, j: usize| {
        ((p[i].x - p[j].x).powi(2) + (p[i].y - p[j].y).powi(2) + (p[i].z - p[j].z).powi(2)).sqrt()
    };
    for (i, j) in [(0, 1), (0, 4), (2, 12)] {
        assert_close_to(dist(&after, i, j), dist(&before, i, j), 1e-6);
    }
    let mut fit = phenol();
    let params = AlignMolShapeInputParams {
        fit_conf_id: 0,
        ..Default::default()
    };
    assert_triple(
        align_mol_shape_input_with_params(&reference, &mut fit, &params).unwrap(),
        ETH_PH,
        0.02,
    );

    let mut fit_shape = ShapeInput::new(&phenol()).unwrap();
    let (combo, shape, color, matrix) = align_shapes(&reference, &mut fit_shape).unwrap();
    assert_triple((combo, shape, color), ETH_PH, 0.02);
    assert_rigid(&matrix);
    let colorless = ShapeOverlayOptions {
        optim_mode: OptimMode::ShapeOnly,
        ..Default::default()
    };
    let mut fit_shape = ShapeInput::new(&phenol()).unwrap();
    let (combo, ..) =
        align_shapes_with_overlay_opts(&reference, &mut fit_shape, &colorless).unwrap();
    assert_close_to(combo, 0.3493199960520435, 0.02);
}

#[test]
fn scores_in_place() {
    let (eth, eth2, ph) = (ethanol(), ethanol2(), phenol());
    assert_five(score_mol(&eth, &eth2).unwrap(), ETH_ETH2_SCORES);
    assert_five(
        score_mol(&eth, &ph).unwrap(),
        (
            0.29853368565780236,
            0.5644258975624138,
            0.032641473753190896,
            147.1106399955384,
            0.8401233386532209,
        ),
    );
    let params = ScoreMolParams {
        ref_conf_id: 0,
        fit_conf_id: 0,
        ..Default::default()
    };
    assert_five(
        score_mol_with_params(&eth, &eth2, &params).unwrap(),
        ETH_ETH2_SCORES,
    );
    let reference = ShapeInput::new(&eth).unwrap();
    assert_five(
        score_mol_shape_input(&reference, &eth2).unwrap(),
        ETH_ETH2_SCORES,
    );
    assert_five(
        score_mol_shape_input_with_params(&reference, &eth2, &ScoreMolShapeInputParams::default())
            .unwrap(),
        ETH_ETH2_SCORES,
    );
    let fit = ShapeInput::new(&eth2).unwrap();
    assert_five(score_shape(&reference, &fit), ETH_ETH2_SCORES);
    let (combo, shape, color, _, color_overlap) = score_shape_with_overlay_opts(
        &reference,
        &fit,
        &ShapeOverlayOptions {
            optim_mode: OptimMode::ShapeOnly,
            ..Default::default()
        },
    );
    assert_close_to(shape, 0.9778488248013629, 1e-6);
    assert_close_to(combo, shape, 1e-6);
    assert_eq!((color, color_overlap), (0.0, 0.0));
    assert_value_error(score_mol_with_params(
        &eth,
        &eth2,
        &ScoreMolParams {
            ref_conf_id: 3,
            ..Default::default()
        },
    ));
}

#[test]
fn compares_conformers() {
    let two = ethanol_two_confs();
    let ph = phenol();
    let before = ph.get_conformer().unwrap().get_positions();
    let (scores, best_ref, best_fit, matrix) = score_molecule_all_conformers(&two, &ph).unwrap();
    assert_eq!((scores.len(), scores[0].len()), (2, 1));
    assert_close_to(scores[0][0], 0.3526061582312596, 0.02);
    assert_close_to(scores[1][0], 0.35665361038718346, 0.02);
    assert_eq!(best_fit, 0);
    assert!(best_ref == 0 || best_ref == 1);
    assert_rigid(&matrix);
    assert_eq!(ph.get_conformer().unwrap().get_positions(), before);
    let (scores, best_ref, best_fit, _) =
        score_molecule_all_conformers_with_params(&ph, &two, &Default::default()).unwrap();
    assert_eq!((scores.len(), scores[0].len()), (1, 2));
    assert_close_to(scores[0][0], 0.35264885636108234, 0.02);
    assert_close_to(scores[0][1], 0.3567082565391058, 0.02);
    assert_eq!(best_ref, 0);
    assert!(best_fit == 0 || best_fit == 1);

    let mut shapes = ShapeInput::new(&two).unwrap();
    let fit = ShapeInput::new(&ph).unwrap();
    let ((combo, shape, color), this_shape, fit_shape, matrix) =
        shapes.best_similarity(&fit).unwrap();
    assert_triple(
        (combo, shape, color),
        (0.3566536103871839, 0.6194420298721873, 0.09386519090218054),
        0.02,
    );
    assert!(this_shape < 2);
    assert_eq!(fit_shape, 0);
    assert_rigid(&matrix);
    assert_eq!(shapes.get_active_shape(), 0);
    let params = BestSimilarityParams {
        threshold: 0.99,
        ..Default::default()
    };
    let (scores, this_shape, fit_shape, _) =
        shapes.best_similarity_with_params(&fit, &params).unwrap();
    assert_eq!(scores, (-1.0, -1.0, -1.0));
    assert_eq!((this_shape, fit_shape), (u32::MAX, u32::MAX));
    assert_close(shapes.max_possible_similarity(&fit), 0.49603236154907154);
    let colorless = ShapeOverlayOptions {
        opt_param: 0.0,
        ..Default::default()
    };
    assert_close(
        shapes.max_possible_similarity_with_overlay_opts(&fit, &colorless),
        0.32539805643147646,
    );
}

// ShapeInput's constructor trusts its options' atom indices: atomSubset's
// write into a bitset unchecked (past its end in Python, silently).
#[test]
fn atom_indices_are_checked() {
    let eth = ethanol();
    assert_index_error(shape_with(
        &eth,
        ShapeInputOptions {
            atom_subset: vec![0, 100],
            ..Default::default()
        },
    ));
    assert_index_error(shape_with(
        &eth,
        ShapeInputOptions {
            atom_radii: vec![(9, 1.5)],
            ..Default::default()
        },
    ));
    assert_value_error(ShapeInput::new_with_params(
        &eth,
        &ShapeInputParams {
            conf_id: 1,
            ..Default::default()
        },
    ));
    assert_value_error(ShapeInput::new(&mol_from_smiles("CCO").unwrap()));
}

// Custom features are counted from the first conformer, and the others'
// read past their ends (a segfault in Python).
#[test]
fn custom_features_are_checked() {
    let two = ethanol_two_confs();
    let feature = |x: f64, radius: f64| (1, Point3D::new(x, 0.0, 0.0), radius, vec![]);
    let uneven = ShapeInputOptions {
        custom_features: vec![vec![feature(0.0, 1.0), feature(1.0, 1.0)], vec![]],
        ..Default::default()
    };
    assert_value_error(shape_with(&two, uneven.clone()));
    // One conformer's features are all that one conformer needs.
    let one = ShapeInput::new_with_params(
        &two,
        &ShapeInputParams {
            conf_id: 0,
            shape_opt: uneven,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(one.num_features(), 2);
    // Too few lists for the conformers is RDKit's own error.
    let short = ShapeInputOptions {
        custom_features: vec![vec![feature(0.0, 1.0)]],
        ..Default::default()
    };
    assert!(matches!(
        shape_with(&two, short),
        Err(Error::Runtime { .. })
    ));
    for bad in [
        feature(f64::NAN, 1.0),
        feature(0.0, 0.0),
        feature(0.0, f64::INFINITY),
    ] {
        assert_value_error(shape_with(
            &ethanol(),
            ShapeInputOptions {
                custom_features: vec![vec![bad]],
                ..Default::default()
            },
        ));
    }
}

// Volumes that are not numbers reach a std::sort, and shapes without volume
// score 0 / 0 in every pose (a NaN transformation in Python).
#[test]
fn shapes_need_volume() {
    for v in [f64::NAN, f64::INFINITY, 1e9] {
        let mut odd = ethanol();
        odd.get_conformer_mut()
            .unwrap()
            .set_atom_position(1, Point3D::new(v, 0.0, 0.0))
            .unwrap();
        assert_value_error(ShapeInput::new(&odd));
    }
    for radius in [0.0, 1e-4, 1e4, f64::NAN] {
        assert_value_error(shape_with(
            &ethanol(),
            ShapeInputOptions {
                atom_radii: vec![(0, radius)],
                all_carbon_radii: false,
                ..Default::default()
            },
        ));
    }
    let hydrogens = ShapeInputOptions {
        atom_subset: vec![3, 4],
        ..Default::default()
    };
    assert_value_error(shape_with(&ethanol(), hydrogens.clone()));
    let mut fit = ethanol2();
    let params = AlignMolParams {
        fit_opts: hydrogens,
        ..Default::default()
    };
    let before = fit.get_conformer().unwrap().get_positions();
    assert_value_error(align_mol_with_params(&ethanol(), &mut fit, &params));
    assert_eq!(fit.get_conformer().unwrap().get_positions(), before);
}

// The overlay sorts its starting poses by score and returns its best score
// uninitialized unless one is at least 0, and counts nSteps - 10 further
// steps unsigned.
#[test]
fn overlay_options_are_checked() {
    let reference = ShapeInput::new(&ethanol()).unwrap();
    let bad = [
        ShapeOverlayOptions {
            sim_alpha: 0.0,
            sim_beta: 0.0,
            ..Default::default()
        },
        ShapeOverlayOptions {
            sim_alpha: 1.5,
            ..Default::default()
        },
        ShapeOverlayOptions {
            sim_beta: f64::NAN,
            ..Default::default()
        },
        ShapeOverlayOptions {
            opt_param: -0.1,
            ..Default::default()
        },
        ShapeOverlayOptions {
            opt_param: f64::NAN,
            ..Default::default()
        },
        ShapeOverlayOptions {
            n_steps: 5,
            ..Default::default()
        },
    ];
    for overlay_opts in bad {
        let mut fit = ethanol2();
        let params = AlignMolParams {
            overlay_opts: overlay_opts.clone(),
            ..Default::default()
        };
        assert_value_error(align_mol_with_params(&ethanol(), &mut fit, &params));
        let mut fit_shape = ShapeInput::new(&ethanol2()).unwrap();
        assert_value_error(align_shapes_with_overlay_opts(
            &reference,
            &mut fit_shape,
            &overlay_opts,
        ));
        let mut shapes = reference.clone();
        let params = BestSimilarityParams {
            threshold: -1.0,
            overlay_opts: overlay_opts.clone(),
        };
        assert_value_error(shapes.best_similarity_with_params(&fit_shape, &params));
        let params = ScoreMoleculeAllConformersParams {
            overlay_opts,
            ..Default::default()
        };
        assert_value_error(score_molecule_all_conformers_with_params(
            &ethanol(),
            &ethanol2(),
            &params,
        ));
    }
    // A Tversky score weighted to one side is allowed.
    let mut fit = ethanol2();
    let params = AlignMolParams {
        overlay_opts: ShapeOverlayOptions {
            sim_alpha: 0.0,
            ..Default::default()
        },
        ..Default::default()
    };
    let (combo, ..) = align_mol_with_params(&ethanol(), &mut fit, &params).unwrap();
    assert_close_to(combo, 1.0010040051117843, 0.02);
}

// Feature patterns are matched under the control_c lock and the module's
// own; shapes made in threads agree.
#[test]
fn makes_shapes_in_threads() {
    let results: Vec<(f64, f64)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..4)
            .map(|_| {
                scope.spawn(|| {
                    let shape = ShapeInput::new(&phenol()).unwrap();
                    (shape.shape_volume(), shape.color_volume())
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for r in &results {
        assert_eq!(*r, results[0]);
    }
}
