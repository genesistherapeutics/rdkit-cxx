// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Macrocycles with SMARTS queries (RDKitCB_13): match atoms in rings of 12
//! or more atoms.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#macrocycles-with-smarts-queries),
//! by Ivan Tubert-Brohman, David Cosgrove, and Vincent Scalfani.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles};

fn main() -> Result<()> {
    let erythromycin = mol_from_smiles(
        "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
    )?;

    // A SMARTS pattern for ring size >= 12, an RDKit SMARTS extension
    let macro_ = mol_from_smarts("[r{12-}]")?;
    let matches = erythromycin.get_substruct_matches(&macro_)?;
    println!("{matches:?}");
    assert_eq!(
        matches,
        [
            [2],
            [3],
            [4],
            [5],
            [6],
            [8],
            [9],
            [10],
            [11],
            [12],
            [13],
            [14],
            [15],
            [17]
        ]
    );
    Ok(())
}
