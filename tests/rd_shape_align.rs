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

// Expected values come from RDKit 2024.09.6's Python rdShapeAlign on the
// same molecules, after an alignment has set PubChem's distance cutoff.
// Alignments optimize, so their scores are compared to 0.02; shapes, which
// are only centered, are compared closely.

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

fn no_colors() -> PrepareConformerParams {
    PrepareConformerParams {
        use_colors: false,
        ..Default::default()
    }
}

#[test]
fn prepares_shapes() {
    let shape = prepare_conformer(&ethanol()).unwrap();
    // The heavy atoms, then a donor and an acceptor at the oxygen.
    assert_eq!(shape.atom_type_vector(), [0, 0, 0, 1, 2]);
    assert_eq!(shape.volume_atom_index_vector(), [0, 1, 2]);
    assert_eq!(shape.alpha_vector().len(), 5);
    assert_close(shape.alpha_vector()[0], 0.8111517942578272);
    assert_close_to(shape.sov(), 85.41806213981589, 1e-6);
    assert_close_to(shape.sof(), 10.14849439791215, 1e-6);
    assert_close_to(shape.coord()[0], -1.2186000347137451, 1e-6);
    let shift = shape.shift();
    assert_close(shift[0], -0.2546);
    assert_close(shift[1], -0.016066666666666674);
    assert_close(shift[2], -0.20553333333333335);
    let clone = shape.clone();
    assert_eq!(clone.coord(), shape.coord());

    let plain = prepare_conformer_with_params(&ethanol(), &no_colors()).unwrap();
    assert_eq!(plain.atom_type_vector(), [0, 0, 0]);
    assert_close_to(plain.sov(), 85.41806213981589, 1e-6);
    assert_eq!(plain.sof(), 0.0);

    let shape = prepare_conformer(&phenol()).unwrap();
    assert_eq!(shape.atom_type_vector(), [0, 0, 0, 0, 0, 0, 0, 1, 2, 3]);
    assert_close(shape.alpha_vector()[0], 1.0146419171594185);
    assert_close_to(shape.sov(), 286.72416514419785, 1e-6);
    assert_close_to(shape.sof(), 15.222741596868225, 1e-6);
    assert_close_to(shape.coord()[0], 2.370971441268921, 1e-6);
}

#[test]
fn aligns_molecules() {
    let reference = ethanol();
    let mut probe = phenol();
    let scores = align_mol(&reference, &mut probe).unwrap();
    assert_pair(scores, (0.6013704215976929, 3.907557049983984e-07), 0.02);
    assert_close(
        probe.get_double_prop("shape_align_shape_tanimoto").unwrap(),
        scores.0,
    );
    assert_close(
        probe.get_double_prop("shape_align_color_tanimoto").unwrap(),
        scores.1,
    );
    let mut probe = phenol();
    let plain = AlignMolParams {
        use_colors: false,
        ..Default::default()
    };
    assert_pair(
        align_mol_with_params(&reference, &mut probe, &plain).unwrap(),
        (0.6013704215976929, 0.0),
        0.02,
    );
    let mut probe = ethanol2();
    assert_pair(
        align_mol(&reference, &mut probe).unwrap(),
        (0.9962729073701257, 0.9758271088349679),
        0.02,
    );
    let mut probe = phenol();
    assert_pair(align_mol(&phenol(), &mut probe).unwrap(), (1.0, 1.0), 0.02);
}

#[test]
fn aligns_to_shapes() {
    let reference = prepare_conformer(&ethanol()).unwrap();
    let mut probe = phenol();
    assert_pair(
        align_mol_shape_input(&reference, &mut probe).unwrap(),
        (0.6013704215976929, 3.907557049983984e-07),
        0.02,
    );
    let mut probe = ethanol2();
    let params = AlignMolShapeInputParams {
        opt_param: 0.5,
        ..Default::default()
    };
    assert_pair(
        align_mol_shape_input_with_params(&reference, &mut probe, &params).unwrap(),
        (0.9962729073701257, 0.9758271088349679),
        0.02,
    );
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

// RDKit 2024.09 writes each heavy atom's coordinates at its atom index in
// an array of the heavy atoms and features: past the end for a heavy atom
// after the hydrogens.
#[test]
fn heavy_atoms_come_first() {
    let mol = phenol();
    let n = mol.get_num_atoms();
    // The last hydrogen first.
    let order: Vec<u32> = std::iter::once(n - 1).chain(0..n - 1).collect();
    let reordered = renumber_atoms(&mol, &order).unwrap();
    assert_value_error(prepare_conformer(&reordered));
    let reference = ethanol();
    let mut probe = renumber_atoms(&mol, &order).unwrap();
    assert_value_error(align_mol(&reference, &mut probe));
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
fn alignments_need_atoms() {
    let mut hydrogen = mol_from_smiles("[H][H]").unwrap();
    let mut conf = Conformer::new(2);
    conf.set_atom_position(1, Point3D::new(0.74, 0.0, 0.0))
        .unwrap();
    conf.set_3d(true);
    hydrogen.add_conformer(&conf).unwrap();
    let empty = prepare_conformer_with_params(&hydrogen, &no_colors()).unwrap();
    assert!(empty.alpha_vector().is_empty());
    let mut probe = ethanol2();
    assert_value_error(align_mol_shape_input(&empty, &mut probe));
    let mut probe = hydrogen.clone();
    assert_value_error(align_mol(&ethanol(), &mut probe));
    // With colors, only the features' shape would differ from the empty one.
    let reference = prepare_conformer(&ethanol()).unwrap();
    let mut probe = ethanol2();
    let bad = AlignMolShapeInputParams {
        opt_param: 2.0,
        ..Default::default()
    };
    assert_value_error(align_mol_shape_input_with_params(
        &reference, &mut probe, &bad,
    ));
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
            let mut probe = mol1.clone();
            assert!(matches!(
                align_mol_shape_input(&reference, &mut probe),
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
