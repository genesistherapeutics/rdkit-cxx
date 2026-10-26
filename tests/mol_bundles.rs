// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn smarts(text: &str) -> Mol {
    mol_from_smarts(text).unwrap()
}

fn bundle(mols: &[Mol]) -> MolBundle {
    let mut bundle = MolBundle::default();
    for mol in mols {
        bundle.add_mol(mol);
    }
    bundle
}

#[test]
fn add_and_get() {
    let mut bundle = MolBundle::default();
    assert_eq!(bundle.add_mol(&smiles("CCO")), 1);
    assert_eq!(bundle.add_mol(&smiles("c1ccccc1O")), 2);
    assert_eq!(bundle.size(), 2);
    assert_eq!(
        mol_to_smiles(bundle.get_mol(1).unwrap()).unwrap(),
        "Oc1ccccc1"
    );
    assert_eq!(
        bundle.get_mol(5).unwrap_err(),
        Error::Index { index: 5, len: 2 }
    );
}

#[test]
fn substructure_queries() {
    let bundle = bundle(&[smiles("CCO"), smiles("c1ccccc1O")]);
    let query = smarts("cO");
    assert!(bundle.has_substruct_match(&query).unwrap());
    assert_eq!(bundle.get_substruct_match(&query).unwrap(), [5, 6]);
    assert_eq!(bundle.get_substruct_matches(&query).unwrap(), [vec![5, 6]]);
    let none = smarts("N");
    assert!(!bundle.has_substruct_match(&none).unwrap());
    assert!(bundle.get_substruct_match(&none).unwrap().is_empty());

    // The first query molecule that matches gives the match.
    let queries = self::bundle(&[smarts("N"), smarts("CO")]);
    assert!(bundle.has_substruct_match_mol_bundle(&queries).unwrap());
    assert_eq!(
        bundle.get_substruct_match_mol_bundle(&queries).unwrap(),
        [1, 2]
    );
    assert_eq!(
        bundle.get_substruct_matches_mol_bundle(&queries).unwrap(),
        [vec![1, 2]]
    );

    let mol = smiles("OCCN");
    assert!(mol.has_substruct_match_mol_bundle(&queries).unwrap());
    assert_eq!(mol.get_substruct_match_mol_bundle(&queries).unwrap(), [3]);
    assert_eq!(
        mol.get_substruct_matches_mol_bundle(&queries).unwrap(),
        [vec![3]]
    );

    let params = SubstructMatchParameters {
        max_matches: 1,
        ..Default::default()
    };
    let matches = bundle
        .get_substruct_matches_with_params(&smarts("[#6]"), &params)
        .unwrap();
    assert_eq!(matches.len(), 1);
    assert!(
        mol.has_substruct_match_mol_bundle_with_params(&queries, &params)
            .unwrap()
    );
}

#[test]
fn pickles() {
    assert!(mol_bundle_can_serialize());
    let bundle = bundle(&[smiles("CCO"), smiles("c1ccccc1O")]);
    // SAFETY: a pickle RDKit wrote.
    let read = unsafe { MolBundle::from_pkl_string(&bundle.to_binary().unwrap()) }.unwrap();
    assert_eq!(read.size(), 2);
    assert_eq!(
        mol_to_smiles(read.get_mol(1).unwrap()).unwrap(),
        "Oc1ccccc1"
    );
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { MolBundle::from_pkl_string(b"garbage") },
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn fixed_mol_size_bundles() {
    let mut bundle = FixedMolSizeMolBundle::default();
    assert_eq!(bundle.add_mol(&smiles("CCO")).unwrap(), 1);
    assert_eq!(bundle.add_mol(&smiles("CCN")).unwrap(), 2);
    let Err(Error::Value { message }) = bundle.add_mol(&smiles("CCCO")) else {
        panic!("added a molecule of another size");
    };
    assert!(message.contains("same number of atoms"), "{message}");
    let Err(Error::Value { message }) = bundle.add_mol(&smiles("C1CC1")) else {
        panic!("added a molecule of another size");
    };
    assert!(message.contains("same number of bonds"), "{message}");
    assert_eq!(bundle.size(), 2);

    // Everything but add_mol comes from MolBundle.
    assert_eq!(bundle.get_substruct_match(&smarts("N")).unwrap(), [2]);
    assert!(
        smiles("NCC")
            .has_substruct_match_mol_bundle(&bundle)
            .unwrap()
    );
}

#[test]
fn editable_mol() {
    let mut em = EditableMol::new(&smiles("CCC"));
    assert_eq!(em.add_atom(&Atom::from_num(8)), 3);
    assert_eq!(em.add_bond_with_order(0, 3, BondType::Single).unwrap(), 3);
    let mut mol = em.get_mol();
    sanitize_mol(&mut mol).unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "CCCO");

    em.remove_bond(0, 3).unwrap();
    em.remove_atom(3).unwrap();
    assert_eq!(mol_to_smiles(&em.get_mol()).unwrap(), "CCC");
    assert_eq!(
        em.remove_atom(9).unwrap_err(),
        Error::Index { index: 9, len: 3 }
    );

    em.begin_batch_edit().unwrap();
    em.remove_atom(0).unwrap();
    em.rollback_batch_edit();
    assert_eq!(em.get_mol().get_num_atoms(), 3);
    em.begin_batch_edit().unwrap();
    em.remove_atom(0).unwrap();
    em.commit_batch_edit().unwrap();
    assert_eq!(em.get_mol().get_num_atoms(), 2);

    em.replace_atom(0, &Atom::from_num(7)).unwrap();
    assert_eq!(mol_to_smiles(&em.get_mol()).unwrap(), "CN");
    let double = smiles("C=C");
    let bond = double.get_bond_with_idx(0).unwrap();
    em.replace_bond_with_preserve_props(0, &bond, true).unwrap();
    assert_eq!(mol_to_smiles(&em.get_mol()).unwrap(), "C=N");
    assert!(matches!(em.add_bond(0, 1), Err(Error::Runtime { .. })));
}
