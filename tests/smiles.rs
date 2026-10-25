// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::{Mol, mol_from_smiles, mol_to_smiles};

#[test]
fn round_trip_canonicalizes() {
    let mol = mol_from_smiles("OCC").unwrap();
    assert_eq!((mol.get_num_atoms(), mol.get_num_bonds()), (3, 2));
    assert_eq!(mol_to_smiles(&mol).unwrap(), "CCO");
}

#[test]
fn counts_with_implicit_hydrogens() {
    let mol = mol_from_smiles("OCC").unwrap();
    assert_eq!(mol.get_num_atoms_with_only_explicit(false).unwrap(), 9);
    assert_eq!(mol.get_num_bonds_with_only_heavy(false).unwrap(), 8);
    assert_eq!(mol.get_num_heavy_atoms(), 3);
}

#[test]
fn empty_smiles_is_an_empty_molecule() {
    let mol = mol_from_smiles("").unwrap();
    assert_eq!(mol.get_num_atoms(), 0);
    assert_eq!(mol_to_smiles(&mol).unwrap(), "");
}

#[test]
fn unparsable_smiles_is_a_parse_error() {
    assert!(matches!(mol_from_smiles("C1CC"),
                     Err(Error::Parse { message }) if message.contains("C1CC")));
}

#[test]
fn sanitization_errors_keep_their_atoms() {
    assert!(matches!(mol_from_smiles("C(C)(C)(C)(C)C"),
                     Err(Error::AtomValence { atom_idx: 0, message }) if message.contains("valence")));
    assert!(matches!(mol_from_smiles("c1cccc1"),
                     Err(Error::Kekulize { atom_indices, .. }) if atom_indices == [0, 1, 2, 3, 4]));
    assert!(matches!(
        mol_from_smiles("cC"),
        Err(Error::AtomKekulize { atom_idx: 0, .. })
    ));
}

#[test]
fn clone_outlives_the_original() {
    let mol = mol_from_smiles("OCC").unwrap();
    let copy = mol.clone();
    drop(mol);
    assert_eq!(mol_to_smiles(&copy).unwrap(), "CCO");
}

#[test]
fn default_is_empty() {
    assert_eq!(Mol::default().get_num_atoms(), 0);
}

#[test]
fn debug_shows_counts() {
    assert_eq!(
        format!("{:?}", mol_from_smiles("OCC").unwrap()),
        "Mol { num_atoms: 3, num_bonds: 2 }"
    );
}

#[test]
fn mol_is_send() {
    fn assert_send<T: Send>() {}
    assert_send::<Mol>();
}

#[test]
fn random_smiles_from_several_threads_parse_back() {
    // Random SMILES draw from RDKit's global generator, which the wrapper
    // locks (U4).
    use rdkit_cxx::chem::{SmilesWriteParams, mol_to_smiles_with_params};
    let params = SmilesWriteParams {
        do_random: true,
        canonical: false,
        ..SmilesWriteParams::default()
    };
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                let mol = mol_from_smiles("c1ccccc1C(=O)NCCO").unwrap();
                for _ in 0..50 {
                    let smiles = mol_to_smiles_with_params(&mol, &params).unwrap();
                    let back = mol_from_smiles(&smiles).unwrap();
                    assert_eq!(mol_to_smiles(&back).unwrap(), "O=C(NCCO)c1ccccc1");
                }
            });
        }
    });
}
