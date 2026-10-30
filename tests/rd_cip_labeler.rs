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

// RDKit 2024.09's Sp2Bond reads no neighbor of its own, so the double bond
// whose stereo atom is its other end is labeled as Python labels it. Its
// AtropisomerBond has no carriers for a bond with a bare end and reads them
// once the other end has two other bonds (Python segfaults for "CC(C)O");
// every such bond is left out and its label cleared, as RDKit 2026.09 does.
#[test]
fn bonds_without_carriers() {
    let mut mol = mol_from_smiles("[C]=CC/C=C/C").unwrap();
    let mut bond = mol.get_bond_with_idx_mut(0).unwrap();
    bond.set_stereo_atoms(1, 2).unwrap();
    bond.set_stereo(BondStereo::Stereocis);
    bond.set_prop("_CIPCode", "X");
    assign_cip_labels(&mut mol).unwrap();
    assert_eq!(
        bond_labels(&mol),
        [(0, "Z".to_string()), (3, "E".to_string())]
    );

    for smiles in ["CC", "CCO", "CC(C)O", "CC(C)C", "OC(C)C"] {
        let mut mol = mol_from_smiles(smiles).unwrap();
        let idx = if smiles == "OC(C)C" { 2 } else { 0 };
        let mut bond = mol.get_bond_with_idx_mut(idx).unwrap();
        bond.set_stereo(BondStereo::Stereoatropcw);
        bond.set_prop("_CIPCode", "X");
        assign_cip_labels(&mut mol).unwrap();
        assert!(bond_labels(&mol).is_empty(), "{smiles}");
    }
}

/// RDKit 2024.09's PairList shifts a 64-bit pairing by 63 less the number of
/// descriptors it pairs, undefined past 64.
#[test]
fn at_most_64_configurations() {
    let chain = |n: usize| format!("C{}C", "[C@H](F)".repeat(n));
    let mut mol = mol_from_smiles(&chain(65)).unwrap();
    assert!(matches!(
        assign_cip_labels(&mut mol),
        Err(Error::Value { .. })
    ));
    let mut mol = mol_from_smiles(&chain(64)).unwrap();
    assign_cip_labels(&mut mol).unwrap();
    // Only those asked for count.
    let mut mol = mol_from_smiles(&chain(65)).unwrap();
    let params = AssignCIPLabelsParams {
        atoms_to_label: Some((1..=64).collect()),
        ..Default::default()
    };
    assign_cip_labels_with_params(&mut mol, &params).unwrap();
}
