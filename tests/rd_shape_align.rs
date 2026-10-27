// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{ETHANOL_CONF2, assert_close, assert_close_to, ethanol};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_depictor::compute_2d_coords;
use rdkit_cxx::chem::rd_shape_align::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Expected values come from RDKit 2025.03.6's Python rdShapeAlign on the
// same molecules. Alignments optimize, so their scores are compared to 0.02
// and their poses by invariants; shapes are compared closely.

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

fn assert_value_error<T>(result: rdkit_cxx::Result<T>) {
    assert!(
        matches!(result, Err(Error::Value { .. })),
        "expected Error::Value"
    );
}

fn assert_pair(actual: (f64, f64), expected: (f64, f64), tol: f64) {
    assert_close_to(actual.0, expected.0, tol);
    assert_close_to(actual.1, expected.1, tol);
}

#[test]
fn option_defaults() {
    let opts = ShapeInputOptions::default();
    assert!(opts.use_colors);
    assert!(!opts.include_dummies);
    assert_close(opts.dummy_radius, 2.16);
    assert!(opts.atom_subset.is_empty());
    assert!(opts.not_color_atoms.is_empty());
    assert!(opts.atom_radii.is_empty());
    let prepare = PrepareConformerParams::default();
    assert_eq!(prepare.conf_id, -1);
    assert_eq!(prepare.opts, opts);
    let align = AlignMolParams::default();
    assert_eq!(
        (
            align.ref_conf_id,
            align.probe_conf_id,
            align.use_colors,
            align.opt_param
        ),
        (-1, -1, true, 1.0)
    );
    assert_eq!((align.max_preiters, align.max_postiters), (10, 30));
    let align = AlignMolShapeInputParams::default();
    assert_eq!(
        (
            align.probe_conf_id,
            align.use_colors,
            align.opt_param,
            align.apply_ref_shift
        ),
        (-1, true, 1.0, false)
    );
    assert_eq!((align.max_preiters, align.max_postiters), (10, 30));
    let align = AlignShapesParams::default();
    assert_eq!(
        (align.opt_param, align.max_preiters, align.max_postiters),
        (1.0, 10, 30)
    );
}

#[test]
fn prepares_shapes() {
    let shape = prepare_conformer(&ethanol()).unwrap();
    let alphas = [
        0.8111517942578272,
        0.8111517942578272,
        1.0146419171594185,
        1.9999715393635706,
        1.9999715393635706,
    ];
    for (alpha, expected) in shape.alpha_vector().iter().zip(alphas) {
        assert_close(*alpha, expected);
    }
    assert_eq!(shape.atom_type_vector(), [0, 0, 0, 1, 2]);
    assert_eq!(shape.volume_atom_index_vector(), [0, 1, 2]);
    assert_close_to(shape.sov(), 85.41806213981589, 1e-6);
    assert_close_to(shape.sof(), 10.14849439791215, 1e-6);
    let shift = shape.shift();
    for (v, expected) in shift
        .iter()
        .zip([-0.2546, -0.016066666666666674, -0.20553333333333335])
    {
        assert_close_to(*v, expected, 1e-6);
    }
    // The shape is centered on its atoms.
    let coord = shape.coord();
    assert_eq!(coord.len(), 15);
    assert_close_to(coord[0], -1.2186000347137451, 1e-6);
    for axis in 0..3 {
        let center: f64 = (0..3).map(|atom| coord[3 * atom + axis]).sum::<f64>() / 3.0;
        assert_close_to(center, 0.0, 1e-6);
    }
    let copy = shape.clone();
    assert_eq!(copy.coord(), coord);

    let shape = prepare_conformer(&phenol()).unwrap();
    assert_eq!(shape.atom_type_vector(), [0, 0, 0, 0, 0, 0, 0, 1, 2, 3]);
    assert_eq!(shape.volume_atom_index_vector(), [0, 1, 2, 3, 4, 5, 6]);
    assert_close_to(shape.sov(), 286.72416514419785, 1e-6);
    assert_close_to(shape.sof(), 15.222741596868225, 1e-6);
    assert_close_to(shape.coord()[0], 2.370971441268921, 1e-6);

    let no_colors = PrepareConformerParams {
        opts: ShapeInputOptions {
            use_colors: false,
            ..Default::default()
        },
        ..Default::default()
    };
    let shape = prepare_conformer_with_params(&phenol(), &no_colors).unwrap();
    assert_eq!(shape.alpha_vector().len(), 7);
    assert_eq!(shape.sof(), 0.0);
}

#[test]
fn shape_options() {
    let mol = ethanol();
    let prepare = |opts: ShapeInputOptions| {
        prepare_conformer_with_params(&mol, &PrepareConformerParams { conf_id: -1, opts }).unwrap()
    };
    let shape = prepare(ShapeInputOptions {
        atom_radii: vec![(0, 2.0)],
        not_color_atoms: vec![2],
        ..Default::default()
    });
    assert_eq!(shape.atom_type_vector(), [0, 0, 0]);
    assert_close(shape.alpha_vector()[0], 0.5860571713512801);
    assert_close_to(shape.sov(), 107.23191337296689, 1e-6);
    let shape = prepare(ShapeInputOptions {
        atom_subset: vec![0, 1],
        ..Default::default()
    });
    assert_eq!(shape.atom_type_vector(), [0, 0]);
    assert_close_to(shape.sov(), 54.617890545023215, 1e-6);
    assert_close_to(shape.shift()[0], 0.22995, 1e-6);

    // PubChem's own features, from an SD property.
    let mut featured = ethanol();
    featured.set_prop("PUBCHEM_PHARMACOPHORE_FEATURES", "1\n1 3 donor\n");
    let shape = prepare_conformer(&featured).unwrap();
    assert_eq!(shape.atom_type_vector(), [0, 0, 0, 4]);
    assert_close_to(shape.sof(), 5.074247198956075, 1e-6);
}

#[test]
fn aligns_molecules() {
    let reference = ethanol();
    let mut probe = ethanol2();
    let scores = align_mol(&reference, &mut probe).unwrap();
    assert_pair(scores, (0.9962729073701257, 0.9758271088349679), 0.02);
    assert_close(
        probe.get_double_prop("shape_align_shape_tanimoto").unwrap(),
        scores.0,
    );
    assert_close(
        probe.get_double_prop("shape_align_color_tanimoto").unwrap(),
        scores.1,
    );
    // The probe's heavy atoms land near the reference's.
    let (aligned, target) = (
        probe.get_conformer().unwrap().get_positions(),
        reference.get_conformer().unwrap().get_positions(),
    );
    for i in 0..3 {
        let (a, b) = (aligned[i], target[i]);
        let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
        assert!(d < 0.3, "atom {i} is {d} away");
    }

    let mut probe = ethanol2();
    let no_colors = AlignMolParams {
        use_colors: false,
        ..Default::default()
    };
    assert_pair(
        align_mol_with_params(&reference, &mut probe, &no_colors).unwrap(),
        (0.9962729073701257, 0.0),
        0.02,
    );
    let mut probe = phenol();
    assert_pair(
        align_mol(&reference, &mut probe).unwrap(),
        (0.6013704215976929, 3.907557049983984e-7),
        0.02,
    );
    let mut probe = phenol();
    let half = AlignMolParams {
        opt_param: 0.5,
        ..Default::default()
    };
    assert_pair(
        align_mol_with_params(&reference, &mut probe, &half).unwrap(),
        (0.599743900555401, 0.08672677808897558),
        0.02,
    );
}

#[test]
fn aligns_to_shapes() {
    let reference = prepare_conformer(&ethanol()).unwrap();
    let mut probe = phenol();
    assert_pair(
        align_mol_shape_input(&reference, &mut probe).unwrap(),
        (0.6013704215976929, 3.907557049983984e-7),
        0.02,
    );
    // With the reference's shift, the probe moves to the reference
    // molecule's frame: its center to ethanol's heavy-atom center.
    let mut probe = ethanol2();
    let params = AlignMolShapeInputParams {
        apply_ref_shift: true,
        ..Default::default()
    };
    align_mol_shape_input_with_params(&reference, &mut probe, &params).unwrap();
    let center = |mol: &Mol| {
        let p = mol.get_conformer().unwrap().get_positions();
        (0..3).map(|i| p[i].x).sum::<f64>() / 3.0
    };
    assert_close_to(center(&probe), center(&ethanol()), 0.05);

    let mut shape = prepare_conformer(&phenol()).unwrap();
    let before = shape.coord();
    let (st, ct, matrix) = align_shapes(&reference, &mut shape).unwrap();
    assert_close_to(st, 0.6013704215976929, 0.02);
    assert!(ct.abs() < 0.02);
    assert_eq!(matrix.len(), 12);
    assert_ne!(shape.coord(), before);
    // The matrix's rotation is orthonormal.
    for i in 0..3 {
        let norm: f64 = (0..3).map(|k| matrix[3 * i + k].powi(2)).sum();
        assert_close_to(norm, 1.0, 1e-4);
    }
    // Moving the conformer the probe shape came from as the shape moved.
    let mut mol = phenol();
    let probe = prepare_conformer(&phenol()).unwrap();
    let mut conf = mol.get_conformer_mut().unwrap();
    transform_conformer(&[0.0; 3], &matrix, &probe, &mut conf).unwrap();
    let moved = conf.get_positions();
    let coord = shape.coord();
    // The shape's points are the heavy atoms in order.
    for (i, p) in moved.iter().take(7).enumerate() {
        assert_close_to(p.x, coord[3 * i], 1e-3);
        assert_close_to(p.y, coord[3 * i + 1], 1e-3);
        assert_close_to(p.z, coord[3 * i + 2], 1e-3);
    }
    // The probe shape is unchanged.
    assert_eq!(probe.coord(), prepare_conformer(&phenol()).unwrap().coord());
    // As in Python, the sizes are checked.
    assert_value_error(transform_conformer(&[0.0; 2], &matrix, &probe, &mut conf));
    assert_value_error(transform_conformer(
        &[0.0; 3],
        &matrix[1..],
        &probe,
        &mut conf,
    ));
    let mut huge = matrix.clone();
    huge[0] = 1e39;
    assert_value_error(transform_conformer(&[0.0; 3], &huge, &probe, &mut conf));
}

#[test]
fn conformers_are_checked() {
    let missing = PrepareConformerParams {
        conf_id: 2,
        ..Default::default()
    };
    assert_value_error(prepare_conformer_with_params(&ethanol(), &missing));
    let mut flat = mol_from_smiles("CCO").unwrap();
    compute_2d_coords(&mut flat).unwrap();
    assert_value_error(prepare_conformer(&flat));
    // PubChem stores the coordinates as float.
    for v in [f64::NAN, f64::INFINITY, 1e9] {
        let mut odd = ethanol();
        odd.get_conformer_mut()
            .unwrap()
            .set_atom_position(1, Point3D::new(v, 0.0, 0.0))
            .unwrap();
        assert_value_error(prepare_conformer(&odd));
    }
}

// RDKit allocates the counts in PUBCHEM_PHARMACOPHORE_FEATURES before it
// checks them: bad_alloc, or the OOM killer, in Python.
#[test]
fn feature_counts_are_checked() {
    for features in [
        "4000000000\n",
        "1\n4000000000 1 donor\n",
        "2\n1 3 donor\n",
        "1\n2 3 donor\n",
    ] {
        let mut mol = ethanol();
        mol.set_prop("PUBCHEM_PHARMACOPHORE_FEATURES", features);
        assert_value_error(prepare_conformer(&mol));
    }
    let mut mol = ethanol();
    mol.set_prop("PUBCHEM_PHARMACOPHORE_FEATURES", "1\n1 3 metal\n");
    assert_value_error(prepare_conformer(&mol));
}

// PubChem's optimizer reads past an empty shape's coordinates (a segfault
// in Python), and returns uninitialized memory as the pose when no pose
// overlaps.
#[test]
fn alignments_need_overlap() {
    let reference = ethanol();
    let hydrogens = ShapeInputOptions {
        atom_subset: vec![3, 4],
        ..Default::default()
    };
    let params = PrepareConformerParams {
        conf_id: -1,
        opts: hydrogens,
    };
    let empty = prepare_conformer_with_params(&reference, &params).unwrap();
    assert!(empty.alpha_vector().is_empty());
    let full = prepare_conformer(&reference).unwrap();
    let mut probe = empty.clone();
    assert_value_error(align_shapes(&full, &mut probe));
    let mut probe = full.clone();
    assert_value_error(align_shapes(&empty, &mut probe));
    let mut probe = ethanol2();
    assert_value_error(align_mol_shape_input(&empty, &mut probe));

    // Shapes far apart overlap in no starting pose. PrepareConformer
    // centers the shapes, so the far one is a shape whose atoms are spread
    // far apart: dihydrogen's two atoms, 100 apart, as dummies.
    let mut far = mol_from_smiles("*.*").unwrap();
    let mut conf = Conformer::new(2);
    conf.set_atom_position(1, Point3D::new(100.0, 0.0, 0.0))
        .unwrap();
    conf.set_3d(true);
    far.add_conformer(&conf).unwrap();
    let dummies = PrepareConformerParams {
        conf_id: -1,
        opts: ShapeInputOptions {
            include_dummies: true,
            ..Default::default()
        },
    };
    let mut shape = prepare_conformer_with_params(&far, &dummies).unwrap();
    let coord = shape.coord();
    assert_value_error(align_shapes(&full, &mut shape));
    assert_eq!(shape.coord(), coord);

    // With colors in common, the color weight must lie in (0, 1].
    for opt_param in [0.0, -0.5, 1.5, f64::NAN] {
        let mut probe = ethanol2();
        let params = AlignMolParams {
            opt_param,
            ..Default::default()
        };
        assert_value_error(align_mol_with_params(&reference, &mut probe, &params));
    }
    // Without colors it only breaks ties, so any value works.
    let mut probe = ethanol2();
    let params = AlignMolParams {
        opt_param: 0.0,
        use_colors: false,
        ..Default::default()
    };
    assert_pair(
        align_mol_with_params(&reference, &mut probe, &params).unwrap(),
        (0.9962729073701257, 0.0),
        0.02,
    );

    // Radii so small that the overlaps' exponential approximation breaks.
    let tiny = PrepareConformerParams {
        conf_id: -1,
        opts: ShapeInputOptions {
            atom_radii: vec![(0, 0.01)],
            ..Default::default()
        },
    };
    let tiny = prepare_conformer_with_params(&reference, &tiny).unwrap();
    let mut probe = tiny.clone();
    assert_value_error(align_shapes(&tiny, &mut probe));
}

// The global state PubChem's code writes is locked; alignments in threads
// take turns and agree.
#[test]
fn aligns_in_threads() {
    let results: Vec<(f64, f64)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..4)
            .map(|_| {
                scope.spawn(|| {
                    let reference = ethanol();
                    let mut probe = phenol();
                    align_mol(&reference, &mut probe).unwrap()
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for r in &results {
        assert_eq!(*r, results[0]);
    }
}

// An FMCS final match check runs inside a substructure match, holding
// control_c_mutex, so it cannot wait for this module's lock, which comes
// before it (CONVENTIONS.md U4): calls that take it fail instead.
#[test]
fn shapes_in_final_match_checks_fail() {
    use rdkit_cxx::chem::rd_fmcs::{MCSFinalMatchCheck, MCSParameters, find_mcs_with_params};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    let calls = Arc::new(AtomicUsize::new(0));
    let reference = prepare_conformer(&ethanol()).unwrap();
    let check = {
        let calls = Arc::clone(&calls);
        MCSFinalMatchCheck::new(move |mol1, _, _, _, _| {
            assert!(matches!(
                prepare_conformer(mol1),
                Err(Error::Runtime { .. })
            ));
            let mut probe = reference.clone();
            assert!(matches!(
                align_shapes(&reference, &mut probe),
                Err(Error::Runtime { .. })
            ));
            calls.fetch_add(1, Ordering::Relaxed);
            Ok(true)
        })
    };
    let params = MCSParameters {
        final_match_checker: Some(check),
        ..MCSParameters::default()
    };
    find_mcs_with_params(&[ethanol(), ethanol2()], &params).unwrap();
    assert!(calls.load(Ordering::Relaxed) > 0);
    prepare_conformer(&ethanol()).unwrap();
}
