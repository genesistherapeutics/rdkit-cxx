// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Expected values are from Python's rdCoordGen. Coordinates depend on the
// platform's floating point, so the tests compare invariants with
// tolerances.

use std::collections::BTreeMap;
use std::path::PathBuf;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_coord_gen::*;
use rdkit_cxx::chem::rd_depictor::{compute_2d_coords, is_coord_gen_support_available};
use rdkit_cxx::chem::{Conformer, Mol, mol_from_smiles};
use rdkit_cxx::geometry::{Point2D, Point3D};

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn position(mol: &Mol, idx: u32) -> Point3D {
    mol.get_conformer().unwrap().get_atom_position(idx).unwrap()
}

fn distance(a: Point3D, b: Point3D) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

fn assert_bond_lengths(mol: &Mol, expected: f64, tol: f64) {
    for b in mol.get_bonds() {
        let length = distance(
            position(mol, b.get_begin_atom_idx()),
            position(mol, b.get_end_atom_idx()),
        );
        assert!(
            (length - expected).abs() < tol,
            "bond length {length}, expected {expected}"
        );
    }
}

fn is_value_error(result: Result<(), Error>, needle: &str) -> bool {
    matches!(result, Err(Error::Value { message }) if message.contains(needle))
}

// A benzene depiction from RDKit's depictor.
fn benzene() -> Mol {
    let mut mol = smiles("c1ccccc1");
    compute_2d_coords(&mut mol).unwrap();
    mol
}

#[test]
fn defaults_match_python() {
    let p = CoordGenParams::default();
    assert!(p.coord_map.is_empty() && p.template_mol.is_none());
    assert_eq!(p.coordgen_scaling, 50.0);
    assert_eq!(p.template_file_dir, PathBuf::new());
    assert_eq!(p.minimizer_precision, 0.01f32);
    assert!(p.dbg_use_constrained && !p.dbg_use_fixed);
    assert!(!p.treat_nonterminal_bonds_to_metal_as_zobs);
    assert_eq!(p.sketcher_coarse_precision(), 0.01f32);
    assert_eq!(p.sketcher_standard_precision(), 1.0);
    assert_eq!(p.sketcher_best_precision(), 3.0);
    assert_eq!(p.sketcher_quick_precision(), 0.2f32);
}

#[test]
fn add_coords_gives_one_2d_conformer() {
    assert!(is_coord_gen_support_available());
    let mut mol = smiles("c1ccccc1CCO");
    add_coords(&mut mol).unwrap();
    add_coords(&mut mol).unwrap();
    assert_eq!(mol.get_num_conformers(), 1);
    assert!(!mol.get_conformer().unwrap().is_3d());
    assert_bond_lengths(&mol, 1.0, 1e-2);

    // coordgenScaling is CoordGen's bond length in RDKit's units.
    let raw = CoordGenParams {
        coordgen_scaling: 1.0,
        ..Default::default()
    };
    add_coords_with_params(&mut mol, &raw).unwrap();
    assert_bond_lengths(&mol, 50.0, 0.5);

    let mut empty = Mol::default();
    add_coords(&mut empty).unwrap();
    assert_eq!(empty.get_num_conformers(), 1);
}

#[test]
fn templates_fix_atoms() {
    let template = benzene();
    let mut params = CoordGenParams {
        dbg_use_fixed: true,
        ..Default::default()
    };
    let map: BTreeMap<u32, Point2D> = (0..6)
        .map(|i| {
            let p = position(&template, i);
            (i, Point2D::new(p.x, p.y))
        })
        .collect();
    params.set_coord_map(&map);
    let mut mol = smiles("c1ccccc1CCO");
    add_coords_with_params(&mut mol, &params).unwrap();
    for i in 0..6 {
        assert!(distance(position(&mol, i), position(&template, i)) < 1e-3);
    }

    let mut params = CoordGenParams {
        dbg_use_fixed: true,
        ..Default::default()
    };
    params.set_template_mol(&template);
    let mut mol = smiles("c1ccccc1CCO");
    add_coords_with_params(&mut mol, &params).unwrap();
    for i in 0..6 {
        assert!(distance(position(&mol, i), position(&template, i)) < 1e-3);
    }

    // CoordGen looks up only the molecule's atoms in the map.
    let past_end = CoordGenParams {
        coord_map: BTreeMap::from([(100, Point2D::new(0.0, 0.0))]),
        ..Default::default()
    };
    let mut mol = smiles("CCO");
    add_coords_with_params(&mut mol, &past_end).unwrap();
}

#[test]
fn add_coords_refuses_what_coord_gen_would_mishandle() {
    let mut mol = smiles("c1ccccc1CCO");
    // CoordGen casts six times the precision to int; a negative one makes it
    // read the best of no solutions (Python raises MemoryError from the
    // garbage it reads).
    for precision in [f32::NAN, f32::INFINITY, -1.0, 4e8] {
        let params = CoordGenParams {
            minimizer_precision: precision,
            ..Default::default()
        };
        assert!(is_value_error(
            add_coords_with_params(&mut mol, &params),
            "minimizerPrecision"
        ));
    }
    for precision in [0.0, 3.0, 3e8] {
        let params = CoordGenParams {
            minimizer_precision: precision,
            ..Default::default()
        };
        add_coords_with_params(&mut mol, &params).unwrap();
    }

    for p in [Point2D::new(f64::NAN, 0.0), Point2D::new(0.0, 2e5)] {
        let params = CoordGenParams {
            coord_map: BTreeMap::from([(0, p)]),
            ..Default::default()
        };
        assert!(is_value_error(
            add_coords_with_params(&mut mol, &params),
            "template coordinates"
        ));
    }
    let mut template = benzene();
    template
        .get_conformer_mut()
        .unwrap()
        .set_atom_position(2, Point3D::new(f64::INFINITY, 0.0, 0.0))
        .unwrap();
    let mut params = CoordGenParams::default();
    params.set_template_mol(&template);
    assert!(is_value_error(
        add_coords_with_params(&mut mol, &params),
        "template coordinates"
    ));
    // A template with two conformers is ignored, coordinates and all.
    template.add_conformer(&Conformer::new(6)).unwrap();
    params.set_template_mol(&template);
    add_coords_with_params(&mut mol, &params).unwrap();
}

#[test]
fn morgan_scores_that_overflow_are_refused() {
    // canonicalOrdering sums Morgan scores in int: a 40-carbon chain
    // overflows them.
    let mut chain = smiles(&"C".repeat(45));
    assert!(is_value_error(add_coords(&mut chain), "Morgan"));
    let mut short = smiles(&"C".repeat(30));
    add_coords(&mut short).unwrap();
    assert_bond_lengths(&short, 1.0, 1e-2);

    // Zero-order bonds do not count. With treatNonterminalBondsToMetalAsZOBs
    // the iron splits this chain in two.
    let text = format!("{0}[Fe]{0}", "C".repeat(25));
    let mut split = smiles(&text);
    assert!(is_value_error(add_coords(&mut split), "Morgan"));
    let zobs = CoordGenParams {
        treat_nonterminal_bonds_to_metal_as_zobs: true,
        ..Default::default()
    };
    add_coords_with_params(&mut split, &zobs).unwrap();
    assert_eq!(split.get_num_conformers(), 1);
}

#[test]
fn template_directories() {
    set_default_template_file_dir("/nonexistent/templates").unwrap();
    let mut mol = smiles("CCO");
    add_coords(&mut mol).unwrap();
    let elsewhere = CoordGenParams {
        template_file_dir: PathBuf::from("/nonexistent/other"),
        ..Default::default()
    };
    add_coords_with_params(&mut mol, &elsewhere).unwrap();
    assert!(matches!(
        set_default_template_file_dir("a\0b"),
        Err(Error::Value { .. })
    ));
    set_default_template_file_dir("").unwrap();
}

#[test]
fn coord_gen_from_several_threads() {
    // CoordGen's templates and template directory are globals the wrappers
    // lock (U4).
    std::thread::scope(|scope| {
        for i in 0..4 {
            scope.spawn(move || {
                let params = CoordGenParams {
                    template_file_dir: PathBuf::from(format!("/nonexistent/{i}")),
                    ..Default::default()
                };
                for text in ["c1ccccc1CCO", "C1CCC2CCCCC2C1", "CC(C)(C)C(=O)N"] {
                    for _ in 0..5 {
                        let mut mol = smiles(text);
                        add_coords_with_params(&mut mol, &params).unwrap();
                        assert_bond_lengths(&mol, 1.0, 0.05);
                        compute_2d_coords(&mut mol).unwrap();
                    }
                }
            });
        }
    });
}
