// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_sln_parse::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdSLNParse.

/// conda-forge's macOS build of RDKit 2024.09.6 ships no SLN parser, and
/// every SLN function is then `Error::Value`; its Linux builds have one.
fn has_sln_parser() -> bool {
    let has = match mol_from_sln("C") {
        Ok(_) => true,
        Err(Error::Value { message }) if message.contains("no SLN parser") => false,
        other => panic!("{other:?}"),
    };
    assert!(has || !cfg!(target_os = "linux"));
    has
}

#[test]
fn molecules_and_queries() {
    if !has_sln_parser() {
        assert!(matches!(mol_from_query_sln("C"), Err(Error::Value { .. })));
        return;
    }
    let mol = mol_from_sln("CH3CH2OH").unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "CCO");
    let raw = mol_from_sln_with_params(
        "CH3CH2OH",
        &MolFromSLNParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(raw.get_num_atoms(), 9);
    assert!(matches!(
        mol_from_sln("CH3CH2(((OH"),
        Err(Error::Parse { .. })
    ));

    let query = mol_from_query_sln("C[1]:C:C:C:C:C:@1").unwrap();
    assert_eq!(query.get_num_atoms(), 6);
    assert!(
        mol_from_smiles("c1ccccc1O")
            .unwrap()
            .has_substruct_match(&query)
            .unwrap()
    );
    let unmerged = mol_from_query_sln_with_params(
        "CH3",
        &MolFromQuerySLNParams {
            merge_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(unmerged.get_num_atoms(), 1);
}

#[test]
fn parsing_from_threads() {
    if !has_sln_parser() {
        return;
    }
    // The parser's debug flag is a global the wrapper guards (U4).
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for _ in 0..20 {
                    assert_eq!(mol_from_sln("CH3CH2OH").unwrap().get_num_atoms(), 3);
                }
            });
        }
    });
}

#[test]
fn params_defaults() {
    assert_eq!(
        MolFromSLNParams::default(),
        MolFromSLNParams {
            sanitize: true,
            debug_parser: false,
        }
    );
    assert_eq!(
        MolFromQuerySLNParams::default(),
        MolFromQuerySLNParams {
            merge_hs: true,
            debug_parser: false,
        }
    );
}
