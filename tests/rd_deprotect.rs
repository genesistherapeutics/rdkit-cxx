// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_deprotect::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdDeprotect.

const BOC: &str = "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]";

fn protected() -> Mol {
    mol_from_smiles("N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C").unwrap()
}

#[test]
fn built_in_deprotections() {
    let all = get_deprotections();
    assert_eq!(all.len(), 25);
    let first = &all[0];
    assert_eq!(first.deprotection_class(), "alcohol");
    assert_eq!(first.abbreviation(), "TBDMS");
    assert_eq!(first.full_name(), "tert-butyldimethylsilyl");
    assert_eq!(
        first.reaction_smarts(),
        "CC(C)([Si](C)(C)[O;H0:1])C>>[O;H1:1]"
    );
    assert_eq!(first.example(), "CC(C)(C)[Si](C)(C)OC>>OC");
    assert!(first.is_valid());
    assert_eq!(first.clone().abbreviation(), "TBDMS");

    let out = deprotect(&protected()).unwrap();
    assert_eq!(mol_to_smiles(&out).unwrap(), "NCc1ccccc1N");
    assert_eq!(
        out.get_prop_auto_convert("DEPROTECTIONS").unwrap(),
        PropValue::StringVec(vec!["Boc".to_string(), "Boc".to_string()])
    );
    assert_eq!(out.get_int_prop("DEPROTECTION_COUNT").unwrap(), 2);

    let mut mol = protected();
    assert!(deprotect_in_place(&mut mol).unwrap());
    assert_eq!(mol_to_smiles(&mol).unwrap(), "NCc1ccccc1N");
    assert!(!deprotect_in_place(&mut mol).unwrap());
}

#[test]
fn custom_deprotections() {
    let boc = DeprotectData::new("amine", BOC, "Boc", "tert-butyloxycarbonyl").unwrap();
    assert!(boc.is_valid());
    assert_eq!(boc.example(), "");
    let out = deprotect_with_deprotections(&protected(), std::slice::from_ref(&boc)).unwrap();
    assert_eq!(mol_to_smiles(&out).unwrap(), "NCc1ccccc1N");
    let mut mol = protected();
    assert!(deprotect_in_place_with_deprotections(&mut mol, &[boc]).unwrap());

    let two_products = DeprotectData::new("x", "C>>C.C", "a", "b").unwrap();
    assert!(!two_products.is_valid());
    // Invalid deprotections are skipped.
    let out = deprotect_with_deprotections(&protected(), &[two_products]).unwrap();
    assert_eq!(out.get_num_atoms(), protected().get_num_atoms());
    assert!(matches!(
        DeprotectData::new("x", "C>>>C", "a", "b"),
        Err(Error::Parse { .. })
    ));
    // Running in place needs one reactant.
    let two_reactants = DeprotectData::new("x", "[N:1].[O:2]>>[N:1][O:2]", "a", "b").unwrap();
    assert!(two_reactants.is_valid());
    assert!(matches!(
        deprotect_in_place_with_deprotections(&mut protected(), &[two_reactants]),
        Err(Error::Value { .. })
    ));
}

#[test]
fn deprotecting_from_threads() {
    // The built-in list's reactions are shared by every thread.
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for _ in 0..5 {
                    let out = deprotect(&protected()).unwrap();
                    assert_eq!(mol_to_smiles(&out).unwrap(), "NCc1ccccc1N");
                }
            });
        }
    });
}
