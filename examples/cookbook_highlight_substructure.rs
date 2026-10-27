// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Highlight a Substructure in a Molecule (RDKitCB_2): draw a molecule with
//! its substructure matches highlighted, as Jupyter shows it.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#highlight-a-substructure-in-a-molecule),
//! by Greg Landrum.
//!
//! Jupyter shows a molecule whose substructure matches were just found with
//! the atoms of every match highlighted, and the bonds between them (RDKit
//! 2024.09's console keeps the matches' atoms in `mol.__sssAtoms` and draws
//! them with `Draw._moltoSVG`). The recipe then sets the atoms to highlight
//! by hand.

mod common;

use common::{MOL_SIZE, show_with_highlights, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles};

/// The highlight color, `MolDrawOptions::highlight_colour`'s default.
const HIGHLIGHT: &str = "#FF7F7F";

/// How often RDKit 2024.09.6's Python draws the highlight color for the
/// matches, and for the atoms set by hand: `Draw._moltoSVG(m, molSize,
/// highlightAtoms, "", True, drawOptions=MolDrawOptions())`.
const SSS_HIGHLIGHTS: usize = 20;
const BY_HAND_HIGHLIGHTS: usize = 24;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let m = mol_from_smiles("c1cc(C(=O)O)c(OC(=O)C)cc1")?;
    let substructure = mol_from_smarts("C(=O)O")?;
    let matches = m.get_substruct_matches(&substructure)?;
    println!("{matches:?}");
    assert_eq!(matches, [[3, 4, 5], [8, 9, 7]]);

    // The console's drawing options.
    let console = MolDrawOptions::default();
    let sss_atoms: Vec<u32> = matches.concat();
    let svg = show_with_highlights(&m, MOL_SIZE, &console, &sss_atoms)?;
    write_image("highlight_substructure.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), SSS_HIGHLIGHTS);

    // You can also set the atoms to highlight by hand:
    let sss_atoms = [0, 1, 2, 6, 11, 12];
    let svg = show_with_highlights(&m, MOL_SIZE, &console, &sss_atoms)?;
    write_image("highlight_substructure_by_hand.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), BY_HAND_HIGHLIGHTS);
    Ok(())
}
