// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! What the wrapper refuses where RDKit 2024.03's algorithms would reach
//! undefined behavior (CONVENTIONS.md U1), and that ordinary uses still
//! match Python's results (RDKit 2024.03.6).

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_cip_labeler::{AssignCIPLabelsParams, assign_cip_labels_with_params};
use rdkit_cxx::chem::rd_fmcs::{MCSParameters, find_mcs, find_mcs_with_params};
use rdkit_cxx::chem::rd_mmpa::fragment_mol_as_smiles;
use rdkit_cxx::chem::rd_mol_hash::{HashFunction, mol_hash};
use rdkit_cxx::chem::rd_rascal_mces::{RascalOptions, find_mces, find_mces_with_opts};
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

/// Without the legacy labels the SMILES parser assigns.
fn without_labels(text: &str) -> Mol {
    let mut mol = smiles(text);
    for idx in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(idx)
            .unwrap()
            .clear_prop("_CIPCode");
    }
    mol
}

fn is_value_error<T: std::fmt::Debug>(result: &Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

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

// MaximumCommonSubgraph::find dereferences the null query molecule StoreAll
// leaves when Verbose is set.
#[test]
fn fmcs_verbose_with_store_all() {
    let mols = [smiles("c1ccccc1C"), smiles("c1ccccc1CC")];
    let both = MCSParameters {
        verbose: true,
        store_all: true,
        ..MCSParameters::default()
    };
    assert!(is_value_error(&find_mcs_with_params(&mols, &both)));
    // Each alone is fine; StoreAll keeps the one MCS, as in Python.
    let store_all = MCSParameters {
        store_all: true,
        ..MCSParameters::default()
    };
    let result = find_mcs_with_params(&mols, &store_all).unwrap();
    assert_eq!(result.num_atoms(), 7);
    let keys: Vec<&String> = result.degenerate_smarts_query_mol_dict().keys().collect();
    assert_eq!(keys, ["[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]"]);
    let verbose = MCSParameters {
        verbose: true,
        ..MCSParameters::default()
    };
    assert_eq!(
        find_mcs_with_params(&mols, &verbose).unwrap().num_atoms(),
        7
    );
}

// Seed::grow shifts 1uLL by a seed's new bonds, undefined for exactly 64,
// which on the supported CPUs leaves the seed ungrown (README, Known
// limitations); queries of 64 or more bonds are searched, as in Python.
#[test]
fn fmcs_query_bonds() {
    let long = smiles(&"C".repeat(65));
    let result = find_mcs(&[long.clone(), long.clone()]).unwrap();
    assert_eq!((result.num_atoms(), result.num_bonds()), (65, 64));
    let chain = smiles(&"C".repeat(64));
    let result = find_mcs(&[chain.clone(), chain]).unwrap();
    assert_eq!((result.num_atoms(), result.num_bonds()), (64, 63));
    let mols = [smiles("CCC"), smiles("CCCC"), smiles(&"C".repeat(70))];
    let result = find_mcs(&mols).unwrap();
    assert_eq!(result.smarts_string(), "[#6]-[#6]-[#6]");
    let half = MCSParameters {
        threshold: 0.5,
        ..MCSParameters::default()
    };
    let result = find_mcs_with_params(&mols, &half).unwrap();
    assert_eq!((result.num_atoms(), result.num_bonds()), (4, 3));
    let mols = [smiles("CCC"), long.clone(), long];
    assert_eq!(find_mcs(&mols).unwrap().smarts_string(), "[#6]-[#6]-[#6]");
    let result = find_mcs_with_params(&mols, &half).unwrap();
    assert_eq!((result.num_atoms(), result.num_bonds()), (65, 64));
}

// checkTimeout reads front() of a vector it has just cleared.
#[test]
fn rascal_timeouts() {
    assert_eq!(RascalOptions::default().timeout, -1);
    let m1 = smiles("c1ccccc1CCO");
    let m2 = smiles("c1ccccc1CCN");
    for timeout in [60, 1, 0, -2, i32::MIN] {
        let opts = RascalOptions {
            timeout,
            ..RascalOptions::default()
        };
        assert!(
            is_value_error(&find_mces_with_opts(&m1, &m2, &opts)),
            "{timeout}"
        );
    }
    // As Python with timeout -1.
    let results = find_mces(&m1, &m2).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].smarts_string(), "c1:c:c:c:c:c:1-CC");
    assert!((results[0].similarity() - 0.7901234567901234).abs() < 1e-12);
}

// PairList::addAndPair shifts 0x1 by 31 - size, negative past 32
// descriptors.
#[test]
fn cip_labeler_descriptor_count() {
    let chain = |n: usize| format!("F{}Cl", "[C@H](O)".repeat(n));
    let mut mol = without_labels(&chain(33));
    let result = assign_cip_labels(&mut mol);
    assert!(is_value_error(&result), "{result:?}");
    // Labeling at most 32 of them is fine, as Python labels them.
    let some = AssignCIPLabelsParams {
        atoms_to_label: Some((1..64).step_by(2).collect()),
        ..AssignCIPLabelsParams::default()
    };
    assign_cip_labels_with_params(&mut mol, &some).unwrap();
    let labels: String = atom_labels(&mol)
        .into_iter()
        .map(|(_, label)| label)
        .collect();
    assert_eq!(labels, "SRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSS");
    // 32, as Python labels them.
    let mut mol = without_labels(&chain(32));
    assign_cip_labels(&mut mol).unwrap();
    let labels: String = atom_labels(&mol)
        .into_iter()
        .map(|(_, label)| label)
        .collect();
    assert_eq!(labels, "SRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSR");
    // Atoms and bonds both count; untagged ones do not.
    let mut mol = without_labels("C[C@H](O)/C=C/[C@@H](N)C");
    assign_cip_labels(&mut mol).unwrap();
    assert_eq!(
        atom_labels(&mol),
        [(1, "S".to_string()), (5, "S".to_string())]
    );
    assert_eq!(bond_labels(&mol), [(3, "E".to_string())]);
}

// RDKit 2024.03's Sp2Bond labels a double bond whose stereo atom is the
// other end of the bond (2026.03 read past an empty vector there).
#[test]
fn cip_labeler_stereo_atom_is_other_end() {
    for (stereo, expected) in [(BondStereo::Stereocis, "Z"), (BondStereo::Stereotrans, "E")] {
        let mut mol = smiles("FC=CF");
        let mut bond = mol.get_bond_with_idx_mut(1).unwrap();
        bond.set_stereo_atoms(2, 3).unwrap();
        bond.set_stereo(stereo);
        assign_cip_labels(&mut mol).unwrap();
        assert_eq!(bond_labels(&mol), [(1, expected.to_string())]);
        let bond = mol.get_bond_with_idx(1).unwrap();
        assert_eq!(bond.get_stereo_atoms(), [2, 3]);
        assert_eq!(bond.get_stereo(), stereo);
    }
}

// ExtendedMurckoScaffold allocates a byte of stack per atom with alloca.
#[test]
fn extended_murcko_atom_limit() {
    // Disconnected methanes keep RDKit's ring search short.
    let big = smiles(&vec!["C"; 4001].join("."));
    assert!(is_value_error(&mol_hash(
        &big,
        HashFunction::ExtendedMurcko
    )));
    let limit = smiles(&vec!["C"; 4000].join("."));
    assert!(mol_hash(&limit, HashFunction::ExtendedMurcko).is_ok());
    // As Python.
    assert_eq!(
        mol_hash(&smiles("c1ccccc1CCOC1CC1CCC"), HashFunction::ExtendedMurcko).unwrap(),
        "*C1CC1OCCc1ccccc1"
    );
}

// computeMorganCodeHash left-shifts a negative formal charge, which GCC and
// Clang define; anions are fragmented as in Python.
#[test]
fn mmpa_anion() {
    let mut pairs = fragment_mol_as_smiles(&smiles("CC(C)CC(=O)[O-]")).unwrap();
    pairs.sort();
    let expected = [
        ("", "CC(C)C[*:1].O=C([O-])[*:1]"),
        ("", "CC(C)[*:1].O=C([O-])C[*:1]"),
        ("", "CC(CC(=O)[O-])[*:1].C[*:1]"),
        ("C(C([*:1])[*:2])[*:3]", "C[*:1].C[*:2].O=C([O-])[*:3]"),
        ("C([*:1])([*:2])[*:3]", "C[*:1].C[*:2].O=C([O-])C[*:3]"),
        ("C([*:1])[*:2]", "CC(C)[*:2].O=C([O-])[*:1]"),
        ("CC(C[*:2])[*:1]", "C[*:1].O=C([O-])[*:2]"),
        ("CC([*:1])[*:2]", "C[*:1].O=C([O-])C[*:2]"),
        ("O=C([O-])CC([*:1])[*:2]", "C[*:1].C[*:2]"),
    ];
    let expected: Vec<(String, String)> = expected
        .iter()
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .collect();
    assert_eq!(pairs, expected);
}

// RDKit 2024.03's getImplicitValence returns -1 before the valences are
// computed, as Python does.
#[test]
fn implicit_valence_before_property_cache() {
    let mut mol = mol_from_smiles_with_params(
        "CC",
        &SmilesParserParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_implicit_valence()
            .unwrap(),
        -1
    );
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .update_property_cache()
        .unwrap();
    assert_eq!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_implicit_valence()
            .unwrap(),
        3
    );
    // Outside a molecule it is still RDKit's precondition.
    assert!(matches!(
        Atom::from_num(6).get_implicit_valence(),
        Err(Error::Runtime { .. })
    ));
}
