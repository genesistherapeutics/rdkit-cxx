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
//! Jupyter shows a molecule whose substructure matches were just found with
//! the matches highlighted (`Draw.DrawMolWithMatches`, with the query, so
//! the matched bonds are highlighted too). The recipe then sets the atoms
//! to highlight by hand, but Jupyter keeps the query, which those atoms do
//! not fit, and in RDKit 2026.09 drawing them fails (Python raises
//! `AttributeError`; here the drawing is `Error::Value`). Drawn without the
//! query, as the recipe intends, the atoms are highlighted.

mod common;

use common::{MOL_SIZE, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::draw::{DrawMolWithMatchesParams, draw_mol_with_matches_as_svg_with_params};
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles};

/// The highlight color, `MolDrawOptions::highlight_colour`'s default.
const HIGHLIGHT: &str = "#FF7F7F";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let m = mol_from_smiles("c1cc(C(=O)O)c(OC(=O)C)cc1")?;
    let substructure = mol_from_smarts("C(=O)O")?;
    let matches = m.get_substruct_matches(&substructure)?;
    println!("{matches:?}");
    assert_eq!(matches, [[3, 4, 5], [8, 9, 7]]);

    let shown = DrawMolWithMatchesParams {
        mol_size: MOL_SIZE,
        qry: Some(substructure),
        options: MolDrawOptions::default(),
        ..Default::default()
    };
    let svg = draw_mol_with_matches_as_svg_with_params(&m, &matches, &shown)?;
    write_image("highlight_substructure.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), 20);

    // You can also set the atoms to highlight by hand:
    let atoms = vec![vec![0, 1, 2, 6, 11, 12]];
    assert!(draw_mol_with_matches_as_svg_with_params(&m, &atoms, &shown).is_err());
    let by_hand = DrawMolWithMatchesParams { qry: None, ..shown };
    let svg = draw_mol_with_matches_as_svg_with_params(&m, &atoms, &by_hand)?;
    write_image("highlight_substructure_by_hand.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), 12);
    Ok(())
}
