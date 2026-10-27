// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Explicit Valence and Number of Hydrogens (RDKitCB_11): an atom's explicit
//! valence, its explicit, implicit, and total hydrogens.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#explicit-valence-and-number-of-hydrogens),
//! by Michael Palmer and Greg Landrum.
//!
//! Explicit usually refers to atoms in the graph and implicit to atoms that
//! are not (hydrogens), so each atom of pyrrole has an explicit valence of 3
//! and one hydrogen. But a hydrogen written inside square brackets in SMILES
//! is "explicit" as far as `get_num_explicit_hs` is concerned.

use rdkit_cxx::Result;
use rdkit_cxx::chem::mol_from_smiles;

/// An atom's symbol, explicit valence, implicit, explicit, and total
/// hydrogens.
type Row = (String, i32, u32, u32, u32);

/// A row for each atom.
fn hydrogens(smiles: &str) -> Result<Vec<Row>> {
    let mol = mol_from_smiles(smiles)?;
    let mut rows = Vec::new();
    for atom in mol.get_atoms() {
        let row = (
            atom.get_symbol()?,
            atom.get_explicit_valence()?,
            atom.get_num_implicit_hs()?,
            atom.get_num_explicit_hs(),
            atom.get_total_num_hs()?,
        );
        println!("{} {} {} {} {}", row.0, row.1, row.2, row.3, row.4);
        rows.push(row);
    }
    Ok(rows)
}

fn main() -> Result<()> {
    let row = |symbol: &str, valence, implicit, explicit, total| {
        (symbol.to_owned(), valence, implicit, explicit, total)
    };

    let pyrrole = hydrogens("C1=CNC=C1")?;
    assert_eq!(
        pyrrole,
        [
            row("C", 3, 1, 0, 1),
            row("C", 3, 1, 0, 1),
            row("N", 3, 0, 1, 1),
            row("C", 3, 1, 0, 1),
            row("C", 3, 1, 0, 1),
        ]
    );
    println!();

    let mol1 = hydrogens("C1=CNCC1")?;
    assert_eq!(
        mol1,
        [
            row("C", 3, 1, 0, 1),
            row("C", 3, 1, 0, 1),
            row("N", 2, 1, 0, 1),
            row("C", 2, 2, 0, 2),
            row("C", 2, 2, 0, 2),
        ]
    );
    println!();

    let mol2 = hydrogens("C1=C[NH]CC1")?;
    assert_eq!(
        mol2,
        [
            row("C", 3, 1, 0, 1),
            row("C", 3, 1, 0, 1),
            row("N", 3, 0, 1, 1),
            row("C", 2, 2, 0, 2),
            row("C", 2, 2, 0, 2),
        ]
    );
    Ok(())
}
