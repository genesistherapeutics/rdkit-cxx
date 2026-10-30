// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports the doctests of Python's Chem.FindMolChiralCenters, with values from
// RDKit 2025.03.6. The file holds one test, so the one call that turns legacy
// stereo perception off runs with no other thread reading the environment.

use rdkit_cxx::Error;
use rdkit_cxx::chem::{
    FindMolChiralCentersParams, find_mol_chiral_centers, find_mol_chiral_centers_with_params,
    get_use_legacy_stereo_perception, mol_from_smiles, set_use_legacy_stereo_perception,
};

fn centers(smiles: &str, params: &FindMolChiralCentersParams) -> Vec<(u32, String)> {
    let mut mol = mol_from_smiles(smiles).unwrap();
    find_mol_chiral_centers_with_params(&mut mol, params).unwrap()
}

fn expected(codes: &[(u32, &str)]) -> Vec<(u32, String)> {
    codes
        .iter()
        .map(|&(idx, code)| (idx, code.to_owned()))
        .collect()
}

#[test]
fn find_mol_chiral_centers_matches_python() {
    // RDKit 2025.03 perceives stereo with its legacy code by default.
    assert!(get_use_legacy_stereo_perception());

    let mut mol = mol_from_smiles("[C@H](Cl)(F)Br").unwrap();
    assert_eq!(
        find_mol_chiral_centers(&mut mol).unwrap(),
        expected(&[(0, "R")])
    );
    let mut mol = mol_from_smiles("[C@@H](Cl)(F)Br").unwrap();
    assert_eq!(
        find_mol_chiral_centers(&mut mol).unwrap(),
        expected(&[(0, "S")])
    );
    let mut mol = mol_from_smiles("CCC").unwrap();
    assert_eq!(find_mol_chiral_centers(&mut mol).unwrap(), expected(&[]));

    // By default unassigned stereo centers are not reported, but they can be:
    let defaults = FindMolChiralCentersParams::default();
    let unassigned = FindMolChiralCentersParams {
        include_unassigned: true,
        ..Default::default()
    };
    assert_eq!(
        centers("C[C@H](F)C(F)(Cl)Br", &defaults),
        expected(&[(1, "S")])
    );
    assert_eq!(
        centers("C[C@H](F)C(F)(Cl)Br", &unassigned),
        expected(&[(1, "S"), (3, "?")])
    );

    // The legacy implementation mishandles dependent stereochemistry:
    let legacy = FindMolChiralCentersParams {
        include_unassigned: true,
        use_legacy_implementation: Some(true),
        ..Default::default()
    };
    assert_eq!(
        centers("C1CC(C)C(C)C(C)C1", &legacy),
        expected(&[(2, "?"), (6, "?")])
    );
    assert_eq!(
        centers("C1C[C@H](C)C(C)[C@H](C)C1", &legacy),
        expected(&[(2, "S"), (4, "?"), (6, "R")])
    );
    assert_eq!(
        centers("C[C@H]1CCC[C@@H](C)[C@@H]1Cl", &legacy),
        expected(&[(1, "S"), (5, "R"), (7, "R")])
    );

    // The new implementation gets it right, and para-stereochemistry too:
    let new = FindMolChiralCentersParams {
        use_legacy_implementation: Some(false),
        ..Default::default()
    };
    let new_unassigned = FindMolChiralCentersParams {
        include_unassigned: true,
        ..new.clone()
    };
    assert_eq!(
        centers("C1CC(C)C(C)C(C)C1", &new_unassigned),
        expected(&[(2, "?"), (4, "?"), (6, "?")])
    );
    assert_eq!(
        centers("C1C[C@H](C)[C@H](C)[C@H](C)C1", &new),
        expected(&[(2, "S"), (4, "s"), (6, "R")])
    );
    assert_eq!(
        centers("C[C@H]1CCC[C@@H](C)[C@@H]1Cl", &new_unassigned),
        expected(&[(1, "S"), (5, "R"), (7, "r")])
    );

    // Without CIP labels, the new implementation gives the descriptors:
    let no_cip = FindMolChiralCentersParams {
        include_cip: false,
        ..new.clone()
    };
    assert_eq!(
        centers("C1C[C@H](C)[C@H](C)[C@H](C)C1", &no_cip),
        expected(&[(2, "Tet_CCW"), (4, "Tet_CCW"), (6, "Tet_CCW")])
    );
    let no_cip_unassigned = FindMolChiralCentersParams {
        include_unassigned: true,
        ..no_cip
    };
    let mut mol = mol_from_smiles("CC(F)Cl").unwrap();
    assert_eq!(
        find_mol_chiral_centers_with_params(&mut mol, &no_cip_unassigned).unwrap(),
        expected(&[(1, "?")])
    );
    assert!(
        mol.get_atom_with_idx(1)
            .unwrap()
            .has_prop("_ChiralityPossible")
    );

    // Safe Rust cannot switch the process-wide setting for the call, as
    // Python does, so the legacy implementation needs it on.
    // SAFETY: this test is the only one in its process, so no other thread
    // reads the environment.
    unsafe { set_use_legacy_stereo_perception(false) };
    let mut mol = mol_from_smiles("C1CC(C)C(C)C(C)C1").unwrap();
    let result = find_mol_chiral_centers_with_params(&mut mol, &legacy);
    // The default follows the setting, now the new implementation.
    let mut mol = mol_from_smiles("C1CC(C)C(C)C(C)C1").unwrap();
    let default_unassigned = find_mol_chiral_centers_with_params(&mut mol, &unassigned);
    // SAFETY: as above.
    unsafe { set_use_legacy_stereo_perception(true) };
    assert!(matches!(result, Err(Error::Value { .. })));
    assert_eq!(
        default_unassigned.unwrap(),
        expected(&[(2, "?"), (4, "?"), (6, "?")])
    );
}
