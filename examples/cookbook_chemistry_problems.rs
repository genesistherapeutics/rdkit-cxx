// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Detect Chemistry Problems (RDKitCB_14): find the problems sanitization
//! would report, without sanitizing.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#detect-chemistry-problems),
//! by Greg Landrum. Python's problems are exception objects; here they are
//! `Error` variants that keep the same atom indices and messages.

use rdkit_cxx::chem::{SmilesParserParams, detect_chemistry_problems, mol_from_smiles_with_params};
use rdkit_cxx::{Error, Result};

fn main() -> Result<()> {
    let unsanitized = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };

    let m = mol_from_smiles_with_params("CN(C)(C)C", &unsanitized)?;
    let problems = detect_chemistry_problems(&m)?;
    println!("{}", problems.len());
    println!("{:?}", problems[0]);
    assert_eq!(problems.len(), 1);
    let Error::AtomValence { atom_idx, message } = &problems[0] else {
        panic!("expected a valence problem, not {:?}", problems[0]);
    };
    assert_eq!(*atom_idx, 1);
    assert_eq!(
        message,
        "Explicit valence for atom # 1 N, 4, is greater than permitted"
    );

    let m2 = mol_from_smiles_with_params("c1cncc1", &unsanitized)?;
    let problems = detect_chemistry_problems(&m2)?;
    println!("{}", problems.len());
    println!("{:?}", problems[0]);
    assert_eq!(problems.len(), 1);
    let Error::Kekulize {
        atom_indices,
        message,
    } = &problems[0]
    else {
        panic!("expected a kekulization problem, not {:?}", problems[0]);
    };
    assert_eq!(*atom_indices, [0, 1, 2, 3, 4]);
    assert_eq!(message, "Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 4");
    Ok(())
}
