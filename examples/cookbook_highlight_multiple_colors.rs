// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Highlight Molecule with Multiple Colors (RDKitCB_39): highlight a
//! molecule's aromatic atoms and bonds in one color and the rest in
//! another.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#highlight-molecule-with-multiple-colors),
//! by Vincent Scalfani.

mod common;

use std::collections::BTreeMap;

use common::{png_size, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::{DrawColour, MolDraw2DCairo};
use rdkit_cxx::chem::mol_from_smiles;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mol = mol_from_smiles("CC1=CC(=CC=C1)NC(=O)CCC2=CC=CC=C2")?;
    let colors = [
        DrawColour::rgba(0.0, 0.0, 1.0, 0.1),
        DrawColour::rgba(1.0, 0.0, 0.0, 0.2),
    ];
    let mut athighlights: BTreeMap<u32, Vec<DrawColour>> = BTreeMap::new();
    let mut arads = BTreeMap::new();
    for a in mol.get_atoms() {
        let aid = a.get_idx();
        let color = if a.get_is_aromatic() {
            colors[0]
        } else {
            colors[1]
        };
        athighlights.entry(aid).or_default().push(color);
        arads.insert(aid, 0.3);
    }
    let mut bndhighlights: BTreeMap<u32, Vec<DrawColour>> = BTreeMap::new();
    for bond in mol.get_bonds() {
        let (aid1, aid2) = (bond.get_begin_atom_idx(), bond.get_end_atom_idx());
        let bid = mol
            .get_bond_between_atoms(aid1, aid2)?
            .expect("the bond joins its atoms")
            .get_idx();
        let color = if bond.get_is_aromatic() {
            colors[0]
        } else {
            colors[1]
        };
        bndhighlights.entry(bid).or_default().push(color);
    }
    // The two rings' twelve atoms and bonds are aromatic.
    let aromatic = |highlights: &BTreeMap<u32, Vec<DrawColour>>| {
        highlights.values().filter(|c| c[0] == colors[0]).count()
    };
    assert_eq!(aromatic(&athighlights), 12);
    assert_eq!(aromatic(&bndhighlights), 12);

    let mut d2d = MolDraw2DCairo::new(350, 400)?;
    d2d.draw_molecule_with_highlights(
        &mol,
        "",
        &athighlights,
        &bndhighlights,
        &arads,
        &BTreeMap::new(),
    )?;
    d2d.finish_drawing();
    let png = d2d.get_drawing_text()?;
    write_image("highlight_multiple_colors.png", &png)?;
    assert_eq!(png_size(&png), (350, 400));
    Ok(())
}
