// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Highlight Entire Molecule (RDKitCB_38): highlight all of a molecule's
//! atoms and bonds.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#highlight-entire-molecule),
//! by Vincent Scalfani.

mod common;

use common::{png_size, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::{
    DrawColour, MolDraw2DCairo, PrepareAndDrawMoleculeParams, prepare_and_draw_molecule_with_params,
};
use rdkit_cxx::chem::mol_from_smiles;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mol = mol_from_smiles("CC(C)CN1C(=O)COC2=C1C=CC(=C2)NC(=O)/C=C/C3=CC=CC=C3")?;
    let rgba_color = DrawColour::rgba(0.0, 0.0, 1.0, 0.1); // transparent blue
    let atoms: Vec<u32> = mol.get_atoms().map(|a| a.get_idx()).collect();
    let mut bonds = Vec::new();
    for bond in mol.get_bonds() {
        let aid1 = atoms[bond.get_begin_atom_idx() as usize];
        let aid2 = atoms[bond.get_end_atom_idx() as usize];
        let bond = mol
            .get_bond_between_atoms(aid1, aid2)?
            .expect("the bond joins its atoms");
        bonds.push(bond.get_idx());
    }
    assert_eq!(bonds.len() as u32, mol.get_num_bonds());

    let mut drawer = MolDraw2DCairo::new(350, 300)?;
    let mut options = drawer.draw_options();
    options.fill_highlights = true;
    options.set_highlight_colour(rgba_color)?;
    options.highlight_bond_width_multiplier = 20;
    options.clear_background = false;
    drawer.set_draw_options(&options)?;
    prepare_and_draw_molecule_with_params(
        &mut drawer,
        &mol,
        &PrepareAndDrawMoleculeParams {
            highlight_atoms: atoms,
            highlight_bonds: bonds,
            ..Default::default()
        },
    )?;
    let png = drawer.get_drawing_text()?;
    write_image("highlight_molecule.png", &png)?;
    assert_eq!(png_size(&png), (350, 300));
    Ok(())
}
