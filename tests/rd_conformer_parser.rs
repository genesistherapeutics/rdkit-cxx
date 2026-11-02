// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_conformer_parser::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Expected values come from Python's rdConformerParser.

fn trajectory(name: &str) -> std::path::PathBuf {
    let path = common::temp_path(name);
    std::fs::write(&path, "title\n0 0 0 1.5 0 0 2 1 0\n0 0 1 1.5 0 1 2 1 1\n").unwrap();
    path
}

#[test]
fn amber_trajectories() {
    let path = trajectory("traj.crd");
    let mut mol = mol_from_smiles("CCO").unwrap();
    assert_eq!(
        add_conformers_from_amber_trajectory(&mut mol, &path).unwrap(),
        [0, 1]
    );
    assert_eq!(mol.get_num_conformers(), 2);
    assert_eq!(
        mol.get_conformer_with_id(1)
            .unwrap()
            .get_atom_position(2)
            .unwrap(),
        Point3D::new(2.0, 1.0, 1.0)
    );
    let one = AddConformersFromAmberTrajectoryParams {
        num_confs: 1,
        clear_confs: false,
    };
    assert_eq!(
        add_conformers_from_amber_trajectory_with_params(&mut mol, &path, &one).unwrap(),
        [2]
    );
    assert_eq!(mol.get_num_conformers(), 3);
    std::fs::remove_file(path).unwrap();
}

#[test]
fn bad_input() {
    let path = trajectory("bad.crd");
    let mut mol = mol_from_smiles("CCO").unwrap();
    add_conformers_from_amber_trajectory(&mut mol, &path).unwrap();
    assert!(matches!(
        add_conformers_from_amber_trajectory(&mut mol, common::temp_path("missing.crd")),
        Err(Error::Os { .. })
    ));
    // RDKit asserts that the file holds the frames asked for.
    let five = AddConformersFromAmberTrajectoryParams {
        num_confs: 5,
        ..Default::default()
    };
    assert!(matches!(
        add_conformers_from_amber_trajectory_with_params(&mut mol, &path, &five),
        Err(Error::Value { .. })
    ));
    // RDKit divides by three times the atoms.
    let mut empty = Mol::default();
    assert!(matches!(
        add_conformers_from_amber_trajectory(&mut empty, &path),
        Err(Error::Value { .. })
    ));
    let mut four = mol_from_smiles("CCCO").unwrap();
    assert!(matches!(
        add_conformers_from_amber_trajectory(&mut four, &path),
        Err(Error::Value { .. })
    ));
    std::fs::write(&path, "title\n0 0 zero\n").unwrap();
    assert!(matches!(
        add_conformers_from_amber_trajectory(&mut mol, &path),
        Err(Error::Value { .. })
    ));
    // A failed read leaves the conformers alone; Python clears them first.
    assert_eq!(mol.get_num_conformers(), 2);
    std::fs::remove_file(path).unwrap();
}
