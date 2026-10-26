// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Functional Group with SMARTS queries (RDKitCB_10): match primary and
//! secondary alcohols with SMARTS.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#functional-group-with-smarts-queries),
//! by Paolo Tosco.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles};

fn main() -> Result<()> {
    let sucrose =
        "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O";
    let sucrose_mol = mol_from_smiles(sucrose)?;

    let primary_alcohol = mol_from_smarts("[CH2][OH1]")?;
    let matches = sucrose_mol.get_substruct_matches(&primary_alcohol)?;
    println!("{matches:?}");
    assert_eq!(matches, [[0, 22], [13, 14], [17, 18]]);

    let secondary_alcohol = mol_from_smarts("[CH1][OH1]")?;
    let matches = sucrose_mol.get_substruct_matches(&secondary_alcohol)?;
    println!("{matches:?}");
    assert_eq!(matches, [[2, 21], [3, 20], [4, 19], [9, 16], [10, 15]]);
    Ok(())
}
