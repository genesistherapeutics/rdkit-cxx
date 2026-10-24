// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Helpers shared by the integration tests.

#![allow(dead_code)] // each test file uses some of these

use std::path::PathBuf;

use rdkit_cxx::chem::{mol_from_smiles, mol_to_mol_block};

/// A path in the temporary directory that no other test process uses.
pub fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rdkit-cxx-{}-{name}", std::process::id()))
}

/// Writes a mol file for methanol whose title line is Latin-1, not UTF-8.
pub fn latin1_mol_file() -> PathBuf {
    let block = mol_to_mol_block(&mol_from_smiles("CO").unwrap()).unwrap();
    let rest = block.split_once('\n').unwrap().1;
    let path = temp_path("latin1.mol");
    std::fs::write(&path, [b"\xe9t\xe9\n".as_slice(), rest.as_bytes()].concat()).unwrap();
    path
}
