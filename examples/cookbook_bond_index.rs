// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Include a Bond Index (RDKitCB_40): draw a molecule with its bond indices.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#include-a-bond-index),
//! by Jeremy Monat.

mod common;

use common::{MOL_SIZE, count_class, show, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::mol_from_smiles;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let console = MolDrawOptions::default();

    // Test in a kinase inhibitor
    let mol =
        mol_from_smiles("C1CC2=C3C(=CC=C2)C(=CN3C1)[C@H]4[C@@H](C(=O)NC4=O)C5=CNC6=CC=CC=C65")?;
    // Default
    let svg = show(&mol, MOL_SIZE, &console)?;
    write_image("bond_index_default.svg", &svg)?;
    assert_eq!(count_class(&svg, "note"), 0);

    // Add bond indices
    let console = MolDrawOptions {
        add_bond_indices: true,
        ..console
    };
    let svg = show(&mol, (350, 300), &console)?;
    write_image("bond_index.svg", &svg)?;
    // A note for each character of each index.
    let digits: usize = (0..mol.get_num_bonds()).map(|i| i.to_string().len()).sum();
    assert_eq!(count_class(&svg, "note"), digits);
    Ok(())
}
