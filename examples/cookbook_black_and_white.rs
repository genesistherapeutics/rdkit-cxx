// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Black and White Molecules (RDKitCB_1): draw molecules and reactions in
//! black and white.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#black-and-white-molecules),
//! by Greg Landrum and Vincent Scalfani.

mod common;

use std::collections::BTreeSet;

use common::{colours, png_size, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::{MolDraw2DCairo, MolDraw2DCairoParams, MolDrawOptions};
use rdkit_cxx::chem::draw::{
    MolsToGridImageParams, mols_to_grid_image_use_svg, mols_to_grid_image_use_svg_with_params,
};
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ms = [
        "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12",
        "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na]",
    ]
    .iter()
    .map(|smiles| mol_from_smiles(smiles))
    .collect::<rdkit_cxx::Result<Vec<_>>>()?;

    // Oxygens are red, nitrogens blue, and sulfurs yellow.
    let svg = mols_to_grid_image_use_svg(&ms)?;
    write_image("black_and_white_colour.svg", &svg)?;
    let black_and_white = BTreeSet::from(["#000000".to_owned(), "#FFFFFF".to_owned()]);
    assert!(colours(&svg).is_superset(&black_and_white));
    assert!(colours(&svg).contains("#FF0000"));

    // With the black and white palette, only black on white.
    let mut console = MolDrawOptions::default();
    console.use_bw_atom_palette();
    let grid = MolsToGridImageParams {
        draw_options: Some(console),
        ..Default::default()
    };
    let svg = mols_to_grid_image_use_svg_with_params(&ms, &grid)?;
    write_image("black_and_white.svg", &svg)?;
    assert_eq!(colours(&svg), black_and_white);

    // Alternatively, use the rdMolDraw2D package
    let mut drawer = MolDraw2DCairo::new_with_params(
        500,
        180,
        &MolDraw2DCairoParams {
            panel_width: 200,
            panel_height: 180,
            ..Default::default()
        },
    )?;
    let mut options = drawer.draw_options();
    options.use_bw_atom_palette();
    drawer.set_draw_options(&options)?;
    drawer.draw_molecules(&ms)?;
    drawer.finish_drawing();
    let png = drawer.get_drawing_text()?;
    write_image("black_and_white_molecules.png", &png)?;
    assert_eq!(png_size(&png), (500, 180));

    // works for reactions too (the reaction is from
    // https://github.com/rdkit/UGM_2020/blob/master/Notebooks/Landrum_WhatsNew.ipynb):
    let rxn = reaction_from_smarts(concat!(
        "[cH:1]:1:[cH:2]:[cH:3]:[cH:4]:[cH:5](-[C:6]#[N:7]):[c:8]:1-[Cl].",
        "[cH:10]:1:[cH:11]:[cH:12](-[Cl:16]):[cH:13]:[cH:14]:[cH:15]:1-B(-O)-O>>",
        "[cH:1]:1:[cH:2]:[cH:3]:[cH:4]:[cH:5](-[C:6]#[N:7]):[c:8]:1-[cH:15]:1[cH:10]:[cH:11]:[cH:12](-[Cl:16]):[cH:13]:[cH:14]:1",
    ))?;
    let mut drawer = MolDraw2DCairo::new(700, 300)?;
    let mut options = drawer.draw_options();
    options.use_bw_atom_palette();
    drawer.set_draw_options(&options)?;
    drawer.draw_reaction(&rxn)?;
    drawer.finish_drawing();
    let png = drawer.get_drawing_text()?;
    write_image("black_and_white_reaction.png", &png)?;
    assert_eq!(png_size(&png), (700, 300));
    Ok(())
}
