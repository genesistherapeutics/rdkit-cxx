// Copyright (c) 2026 Genesis Molecular AI, Inc.
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

// RDKit 2025.03 reads past an empty vector for an atropisomer bond without
// carriers, and Python crashes; RDKit 2026.09 leaves it unlabeled, as here.
#[test]
fn bonds_without_carriers() {
    // A double bond whose bare end's stereo atom is the bond's other atom
    // is labeled from findStereoAtoms, as Python does.
    let mut mol = mol_from_smiles("[C]=CC/C=C/C").unwrap();
    let mut bond = mol.get_bond_with_idx_mut(0).unwrap();
    bond.set_stereo_atoms(1, 2).unwrap();
    bond.set_stereo(BondStereo::Stereocis);
    bond.set_prop("_CIPCode", "Z");
    assign_cip_labels(&mut mol).unwrap();
    assert_eq!(
        bond_labels(&mol),
        [(0, "Z".to_string()), (3, "E".to_string())]
    );

    // An atropisomer bond whose ends have no other bonds.
    let mut mol = mol_from_smiles("CC").unwrap();
    let mut bond = mol.get_bond_with_idx_mut(0).unwrap();
    bond.set_stereo(BondStereo::Stereoatropcw);
    assign_cip_labels(&mut mol).unwrap();
    assert!(bond_labels(&mol).is_empty());
}

// A chain of `centers` tetrahedral stereo atoms.
fn stereo_chain(centers: usize) -> String {
    format!("C{}C", "[C@H](Cl)".repeat(centers))
}

#[test]
fn labeling_refuses_more_than_64_descriptors() {
    // RDKit 2025.03's PairList shifts a 64-bit pairing by the number of
    // descriptors it holds, which is undefined past 64.
    let mut mol = without_labels(&stereo_chain(65));
    assert!(matches!(
        assign_cip_labels(&mut mol),
        Err(Error::Value { .. })
    ));
    // Asking for fewer of them labels those.
    assign_cip_labels_with_params(
        &mut mol,
        &AssignCIPLabelsParams {
            atoms_to_label: Some((1..=64).collect()),
            bonds_to_label: Some(vec![]),
            ..Default::default()
        },
    )
    .unwrap();
    let mut mol = without_labels(&stereo_chain(64));
    assign_cip_labels(&mut mol).unwrap();
    assert_eq!(atom_labels(&mol).len(), 64);
}
