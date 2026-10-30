// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Highlight a Substructure in a Molecule (RDKitCB_2): draw a molecule with
//! its substructure matches highlighted, as Jupyter shows it.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#highlight-a-substructure-in-a-molecule),
//! by Greg Landrum.
//!
//! In RDKit 2025.03, Jupyter shows a molecule whose substructure matches
//! were just found with the matched atoms highlighted (all matches' atoms
//! in one list, Python's `__sssAtoms`), and with them the bonds between
//! those atoms. The recipe then sets the atoms to highlight by hand.

mod common;

use common::{MOL_SIZE, show_with_highlights, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles};

/// The highlight color, `MolDrawOptions::highlight_colour`'s default.
const HIGHLIGHT: &str = "#FF7F7F";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let m = mol_from_smiles("c1cc(C(=O)O)c(OC(=O)C)cc1")?;
    let substructure = mol_from_smarts("C(=O)O")?;
    let matches = m.get_substruct_matches(&substructure)?;
    println!("{matches:?}");
    assert_eq!(matches, [[3, 4, 5], [8, 9, 7]]);

    let options = MolDrawOptions::default();
    let atoms: Vec<u32> = matches.concat();
    let svg = show_with_highlights(&m, MOL_SIZE, &options, &atoms)?;
    write_image("highlight_substructure.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), 20);

    // You can also set the atoms to highlight by hand:
    let svg = show_with_highlights(&m, MOL_SIZE, &options, &[0, 1, 2, 6, 11, 12])?;
    write_image("highlight_substructure_by_hand.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), 24);
    Ok(())
}
