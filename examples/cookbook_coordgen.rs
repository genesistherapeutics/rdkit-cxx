// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Using CoordGen Library (RDKitCB_37): draw a macrocycle with coordinates
//! from CoordGen, which lays out macrocycles better than RDKit's default.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#using-coordgen-library),
//! by Vincent Scalfani.

mod common;

use common::{show, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::{MolDraw2DSVG, MolDrawOptions};
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_coord_gen::add_coords;
use rdkit_cxx::chem::rd_depictor::set_prefer_coord_gen;

const CYCLOSPORIN: &str = "C/C=C/CC(C)C(O)C1C(=O)NC(CC)C(=O)N(C)CC(=O)N(C)C(CC(C)C)C(=O)NC(C(C)C)C(=O)N(C)C(CC(C)C)C(=O)NC(C)C(=O)NC(C)C(=O)N(C)C(CC(C)C)C(=O)N(C)C(CC(C)C)C(=O)N(C)C(C(C)C)C(=O)N1C";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let console = MolDrawOptions::default();
    let mol_size = (350, 300);

    // default drawing
    let mut mol = mol_from_smiles(CYCLOSPORIN)?;
    let default_svg = show(&mol, mol_size, &console)?;
    write_image("coordgen_default.svg", &default_svg)?;

    // with CoordGen
    add_coords(&mut mol)?;
    assert_eq!(mol.get_num_conformers(), 1);
    let svg = show(&mol, mol_size, &console)?;
    write_image("coordgen.svg", &svg)?;
    assert_ne!(svg, default_svg);

    // CoordGen can also lay out what MolDraw2D draws:
    set_prefer_coord_gen(true);
    let mol = mol_from_smiles(CYCLOSPORIN)?;
    let mut drawer = MolDraw2DSVG::new(300, 300)?;
    let mut options = drawer.draw_options();
    options.add_stereo_annotation = false;
    drawer.set_draw_options(&options)?;
    drawer.draw_molecule(&mol)?;
    drawer.finish_drawing();
    let svg = drawer.get_drawing_text()?;
    write_image("coordgen_moldraw2d.svg", &svg)?;
    assert!(svg.starts_with("<?xml") && svg.trim_end().ends_with("</svg>"));
    Ok(())
}
