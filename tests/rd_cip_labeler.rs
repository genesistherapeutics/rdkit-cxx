// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdCIPLabeler.

fn atom_labels(mol: &Mol) -> Vec<(u32, String)> {
    mol.get_atoms()
        .filter(|atom| atom.has_prop("_CIPCode"))
        .map(|atom| (atom.get_idx(), atom.get_prop("_CIPCode").unwrap()))
        .collect()
}

fn bond_labels(mol: &Mol) -> Vec<(u32, String)> {
    mol.get_bonds()
        .filter(|bond| bond.has_prop("_CIPCode"))
        .map(|bond| (bond.get_idx(), bond.get_prop("_CIPCode").unwrap()))
        .collect()
}

fn without_labels(smiles: &str) -> Mol {
    let mut mol = mol_from_smiles(smiles).unwrap();
    for idx in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(idx)
            .unwrap()
            .clear_prop("_CIPCode");
    }
    mol
}

#[test]
fn cip_labels() {
    let mut mol = without_labels("C[C@H](N)C(=O)O");
    assign_cip_labels(&mut mol).unwrap();
    assert_eq!(atom_labels(&mol), [(1, "S".to_string())]);

    let mut mol = without_labels("C/C=C/[C@H](F)Cl");
    assign_cip_labels_with_params(
        &mut mol,
        &AssignCIPLabelsParams {
            atoms_to_label: Some(vec![3]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(atom_labels(&mol), [(3, "R".to_string())]);
    assert!(bond_labels(&mol).is_empty());

    let mut mol = without_labels("C/C=C/[C@H](F)Cl");
    assign_cip_labels_with_params(
        &mut mol,
        &AssignCIPLabelsParams {
            bonds_to_label: Some(vec![1]),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(atom_labels(&mol).is_empty());
    assert_eq!(bond_labels(&mol), [(1, "E".to_string())]);
}

#[test]
fn iteration_limit() {
    // RDKit resolves easy centers in a first pass the limit does not bound,
    // so it takes RDKit's own example of a hard one.
    let mut hard = mol_from_mol_file_with_params(
        concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/cip_max_iterations.mol"
        ),
        &MolFromMolFileParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let limited = AssignCIPLabelsParams {
        max_recursive_iterations: 1,
        ..Default::default()
    };
    assert!(matches!(
        assign_cip_labels_with_params(&mut hard, &limited),
        Err(Error::Runtime { .. })
    ));
    let mut mol = without_labels("C[C@H](N)C(=O)O");
    assign_cip_labels_with_params(
        &mut mol,
        &AssignCIPLabelsParams {
            max_recursive_iterations: 100,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(atom_labels(&mol), [(1, "S".to_string())]);
}

#[test]
fn indices_that_rdkit_trusts() {
    // Python sets each index in a bitset unchecked.
    let mut mol = mol_from_smiles("C[C@H](N)C(=O)O").unwrap();
    for params in [
        AssignCIPLabelsParams {
            atoms_to_label: Some(vec![6]),
            ..Default::default()
        },
        AssignCIPLabelsParams {
            bonds_to_label: Some(vec![6]),
            ..Default::default()
        },
    ] {
        assert!(matches!(
            assign_cip_labels_with_params(&mut mol, &params),
            Err(Error::Index { .. })
        ));
    }
}

// RDKit 2025.09 reads past an empty vector for these bonds, and Python
// crashes; RDKit 2026.09 leaves them unlabeled, as here.
#[test]
fn bonds_without_carriers() {
    // The bare end's stereo atom is the bond's other atom. The other double
    // bond is still labeled, and the stale label cleared.
    let mut mol = mol_from_smiles("[C]=CC/C=C/C").unwrap();
    let mut bond = mol.get_bond_with_idx_mut(0).unwrap();
    bond.set_stereo_atoms(1, 2).unwrap();
    bond.set_stereo(BondStereo::Stereocis);
    bond.set_prop("_CIPCode", "Z");
    assign_cip_labels(&mut mol).unwrap();
    assert_eq!(bond_labels(&mol), [(3, "E".to_string())]);

    // An atropisomer bond whose ends have no other bonds.
    let mut mol = mol_from_smiles("CC").unwrap();
    let mut bond = mol.get_bond_with_idx_mut(0).unwrap();
    bond.set_stereo(BondStereo::Stereoatropcw);
    assign_cip_labels(&mut mol).unwrap();
    assert!(bond_labels(&mol).is_empty());
}

/// A chain of `n` stereocenters.
fn stereocenters(n: usize) -> Mol {
    mol_from_smiles(&format!("C{}C", "[C@H](F)".repeat(n))).unwrap()
}

// RDKit 2025.09 pairs the descriptors of a branch's stereo atoms and bonds
// in a 64-bit mask, by a shift that is undefined past 64.
#[test]
fn at_most_64_stereo_atoms_and_bonds() {
    let mut mol = stereocenters(64);
    assign_cip_labels(&mut mol).unwrap();
    assert_eq!(atom_labels(&mol).len(), 64);
    let mut mol = stereocenters(65);
    match assign_cip_labels(&mut mol) {
        Err(Error::Value { message }) => assert!(message.contains("past 64"), "{message}"),
        other => panic!("{other:?}"),
    }
    // Only the ones asked for count.
    let params = AssignCIPLabelsParams {
        atoms_to_label: Some((1..65).collect()),
        ..Default::default()
    };
    assign_cip_labels_with_params(&mut mol, &params).unwrap();
}
