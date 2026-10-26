// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::{assert_close, assert_close_to, ethanol};
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_free_sasa::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Expected values come from Python's rdFreeSASA on the same block: the
// heavy atoms of the first two residues of 1d3z's first model, from RDKit's
// External/FreeSASA/test_data.
const MET_GLN: &str = "\
ATOM      1  N   MET A   1      52.923 -90.016   8.509  1.00  9.67           N
ATOM      2  CA  MET A   1      51.653 -89.304   8.833  1.00 10.38           C
ATOM      3  C   MET A   1      50.851 -89.086   7.556  1.00  9.62           C
ATOM      4  O   MET A   1      51.414 -89.033   6.462  1.00  9.62           O
ATOM      5  CB  MET A   1      51.976 -87.958   9.485  1.00 13.77           C
ATOM      6  CG  MET A   1      52.864 -87.131   8.557  1.00 16.29           C
ATOM      7  SD  MET A   1      53.355 -85.606   9.403  1.00 17.17           S
ATOM      8  CE  MET A   1      53.359 -84.524   7.954  1.00 16.11           C
ATOM     20  N   GLN A   2      49.535 -88.953   7.698  1.00  9.27           N
ATOM     21  CA  GLN A   2      48.664 -88.733   6.544  1.00  9.07           C
ATOM     22  C   GLN A   2      48.393 -87.246   6.365  1.00  8.72           C
ATOM     23  O   GLN A   2      48.199 -86.506   7.332  1.00  8.22           O
ATOM     24  CB  GLN A   2      47.322 -89.445   6.739  1.00 14.46           C
ATOM     25  CG  GLN A   2      47.506 -90.960   6.658  1.00 17.01           C
ATOM     26  CD  GLN A   2      47.856 -91.369   5.234  1.00 20.10           C
ATOM     27  OE1 GLN A   2      47.605 -90.618   4.294  1.00 21.89           O
ATOM     28  NE2 GLN A   2      48.403 -92.531   5.017  1.00 19.49           N
END
";

const PROTOR: [f64; 17] = [
    1.64, 1.88, 1.61, 1.42, 1.88, 1.88, 1.77, 1.88, 1.64, 1.88, 1.61, 1.42, 1.88, 1.88, 1.61, 1.42,
    1.64,
];

fn met_gln() -> Mol {
    mol_from_pdb_block(MET_GLN).unwrap()
}

fn class_names(mol: &Mol) -> Vec<String> {
    (0..mol.get_num_atoms())
        .map(|i| {
            mol.get_atom_with_idx(i)
                .unwrap()
                .get_prop("SASAClassName")
                .unwrap()
        })
        .collect()
}

fn with_probe(probe_radius: f64) -> CalcSASAParams<'static> {
    CalcSASAParams {
        opts: SASAOpts::new_with_pr(
            SASAAlgorithm::LeeRichards,
            SASAClassifier::Protor,
            probe_radius,
        ),
        ..Default::default()
    }
}

fn assert_value_error<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) {
    assert!(matches!(result, Err(Error::Value { .. })), "{result:?}");
}

#[test]
fn sasa_opts_defaults() {
    let opts = SASAOpts::default();
    assert_eq!(opts.algorithm, SASAAlgorithm::LeeRichards);
    assert_eq!(opts.classifier, SASAClassifier::Protor);
    assert_close(opts.probe_radius, 1.4);
    let opts = SASAOpts::new(SASAAlgorithm::ShrakeRupley, SASAClassifier::Oons);
    assert_eq!(opts.algorithm, SASAAlgorithm::ShrakeRupley);
    assert_eq!(opts.classifier, SASAClassifier::Oons);
    assert_close(opts.probe_radius, 1.4);
    let params = CalcSASAParams::default();
    assert_eq!(params.conf_idx, -1);
    assert!(params.query.is_none());
    assert_eq!(params.opts, SASAOpts::default());
}

#[test]
fn classifies_residue_atoms() {
    let mut mol = met_gln();
    let radii = classify_atoms(&mut mol);
    assert_eq!(radii.len(), 17);
    for (radius, expected) in radii.iter().zip(PROTOR) {
        assert_close(*radius, expected);
    }
    let classes: Vec<i32> = (0..17)
        .map(|i| {
            mol.get_atom_with_idx(i)
                .unwrap()
                .get_int_prop("SASAClass")
                .unwrap()
        })
        .collect();
    assert_eq!(classes, [1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1]);
    assert_eq!(
        class_names(&mol)[..4],
        ["Polar", "Apolar", "Apolar", "Polar"]
    );

    let naccess = SASAOpts::new(SASAAlgorithm::LeeRichards, SASAClassifier::Naccess);
    let radii = classify_atoms_with_options(&mut mol, &naccess);
    assert_close(radii[2], 1.76);
    assert_close(radii[6], 1.85);
    assert_eq!(class_names(&mol)[6], "Apolar");
    let oons = SASAOpts::new(SASAAlgorithm::LeeRichards, SASAClassifier::Oons);
    let radii = classify_atoms_with_options(&mut mol, &oons);
    assert_close(radii[1], 2.0);
    assert_close(radii[2], 1.55);
    assert_eq!(class_names(&mol)[2], "Polar");
}

#[test]
fn atoms_without_residues_get_radius_zero() {
    let mut mol = ethanol();
    assert_eq!(classify_atoms(&mut mol), vec![0.0; 9]);
    assert!(class_names(&mol).iter().all(|name| name == "Unclassified"));
    assert_close_to(
        calc_sasa(&mut mol, &[1.7; 9]).unwrap(),
        249.67315676029267,
        1e-6,
    );
}

#[test]
fn computes_areas() {
    let mut mol = met_gln();
    let radii = classify_atoms(&mut mol);
    let lee_richards = calc_sasa(&mut mol, &radii).unwrap();
    assert_close_to(lee_richards, 467.7077888464553, 1e-6);
    assert_close_to(mol.get_double_prop("SASA").unwrap(), lee_richards, 1e-12);
    let atom0 = mol
        .get_atom_with_idx(0)
        .unwrap()
        .get_double_prop("SASA")
        .unwrap();
    assert_close_to(atom0, 45.462124314419874, 1e-6);

    let shrake = CalcSASAParams {
        opts: SASAOpts::new(SASAAlgorithm::ShrakeRupley, SASAClassifier::Protor),
        ..Default::default()
    };
    assert_close_to(
        calc_sasa_with_params(&mut mol, &radii, &shrake).unwrap(),
        472.9801948817594,
        1e-6,
    );
    assert_close_to(
        calc_sasa_with_params(&mut mol, &radii, &with_probe(2.0)).unwrap(),
        564.07589407925,
        1e-6,
    );
    let shrake2 = CalcSASAParams {
        opts: SASAOpts::new_with_pr(SASAAlgorithm::ShrakeRupley, SASAClassifier::Protor, 2.0),
        ..Default::default()
    };
    assert_close_to(
        calc_sasa_with_params(&mut mol, &radii, &shrake2).unwrap(),
        556.9531192557344,
        1e-6,
    );
}

#[test]
fn queries_limit_the_area() {
    let mut mol = met_gln();
    let radii = classify_atoms(&mut mol);
    let apolar_query = make_free_sasa_a_polar_atom_query();
    let polar_query = make_free_sasa_polar_atom_query();
    assert!(apolar_query.has_query());
    let apolar = CalcSASAParams {
        query: Some(&apolar_query),
        ..Default::default()
    };
    let polar = CalcSASAParams {
        query: Some(&polar_query),
        ..Default::default()
    };
    let apolar = calc_sasa_with_params(&mut mol, &radii, &apolar).unwrap();
    let polar = calc_sasa_with_params(&mut mol, &radii, &polar).unwrap();
    assert_close_to(apolar, 228.42560450066512, 1e-6);
    assert_close_to(polar, 239.28218434579014, 1e-6);

    // As in Python, the query must be a query atom.
    let carbon = Atom::from_num(6);
    let plain = CalcSASAParams {
        query: Some(&carbon),
        ..Default::default()
    };
    assert_value_error(calc_sasa_with_params(&mut mol, &radii, &plain));
}

#[test]
fn conformers_are_checked() {
    let mut mol = ethanol();
    let params = CalcSASAParams {
        conf_idx: 3,
        ..Default::default()
    };
    assert_value_error(calc_sasa_with_params(&mut mol, &[1.7; 9], &params));
    mol.remove_all_conformers();
    assert!(matches!(
        calc_sasa(&mut mol, &[1.7; 9]),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        calc_sasa(&mut Mol::default(), &[]),
        Err(Error::Runtime { .. })
    ));
}

// Python's CalcSASA reads a radius for every atom from the list it is
// given: it crashes for an empty list, and reads past a short one.
#[test]
fn radii_must_match_the_atoms() {
    let mut mol = ethanol();
    assert_value_error(calc_sasa(&mut mol, &[]));
    assert_value_error(calc_sasa(&mut mol, &[1.7; 3]));
    assert_value_error(calc_sasa(&mut mol, &[1.7; 10]));
}

// Each of these crashes Python's CalcSASA, in freesasa's cell lists.
#[test]
fn coordinates_and_radii_are_checked() {
    let mut mol = met_gln();
    let radii = classify_atoms(&mut mol);
    // Every radius plus the probe at most 0, or too small for coordinates
    // this large.
    assert_value_error(calc_sasa_with_params(
        &mut mol,
        &[0.0; 17],
        &with_probe(0.0),
    ));
    assert_value_error(calc_sasa_with_params(
        &mut mol,
        &[-1.0; 17],
        &with_probe(0.5),
    ));
    assert_value_error(calc_sasa_with_params(
        &mut mol,
        &[1e-300; 17],
        &with_probe(0.0),
    ));
    for probe in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_value_error(calc_sasa_with_params(&mut mol, &radii, &with_probe(probe)));
    }
    let mut huge = radii.clone();
    huge[0] = f64::INFINITY;
    assert_value_error(calc_sasa(&mut mol, &huge));
    for (x, y) in [
        (f64::NAN, 0.0),
        (f64::INFINITY, 0.0),
        (1e12, 0.0),
        (1e5, 1e5),
    ] {
        let mut moved = mol.clone();
        moved
            .get_conformer_mut()
            .unwrap()
            .set_atom_position(3, Point3D::new(x, y, y))
            .unwrap();
        assert_value_error(calc_sasa(&mut moved, &radii));
    }
    // A molecule spread over many cells, but fewer than i32 counts, works.
    let mut spread = mol.clone();
    spread
        .get_conformer_mut()
        .unwrap()
        .set_atom_position(3, Point3D::new(500.0, 0.0, 0.0))
        .unwrap();
    assert!(calc_sasa(&mut spread, &radii).unwrap() > 467.0);
}

// freesasa compares an uninitialized buffer with its classifier's names for
// a blank residue or atom name; the wrapper classifies those as unknown
// names, as Python gives for the name "?".
#[test]
fn blank_names_are_unknown() {
    let mut mol = met_gln();
    for (idx, atom_name, residue_name) in [(1, "", "MET"), (9, "CA", "   "), (12, " \t", "\0GLN")] {
        let mut info = mol
            .get_atom_with_idx(idx)
            .unwrap()
            .get_pdb_residue_info()
            .unwrap()
            .unwrap();
        info.set_name(atom_name);
        info.set_residue_name(residue_name);
        mol.get_atom_with_idx_mut(idx)
            .unwrap()
            .set_pdb_residue_info(&info);
    }
    assert!(classify_atoms(&mut mol).is_empty());
    let names = class_names(&mol);
    assert_eq!(names[1], "Unclassified");
    assert_eq!(names[9], "Unclassified");
    assert_eq!(names[12], "Unclassified");
    assert_eq!(names[0], "Polar");
    // The names are restored.
    let info = mol
        .get_atom_with_idx(12)
        .unwrap()
        .get_pdb_residue_info()
        .unwrap()
        .unwrap();
    assert_eq!(info.get_name().unwrap(), " \t");
    assert_eq!(info.get_residue_name().unwrap(), "\0GLN");
    let info = mol
        .get_atom_with_idx(9)
        .unwrap()
        .get_pdb_residue_info()
        .unwrap()
        .unwrap();
    assert_eq!(info.get_name().unwrap(), "CA");
    assert_eq!(info.get_residue_name().unwrap(), "   ");
}
