// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Highlight Molecule Differences (RDKitCB_36): highlight the atoms outside
//! two molecules' maximum common substructure.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#highlight-molecule-differences),
//! by Takayuki Serizawa.

mod common;

use common::write_image;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::draw::{MolsToGridImageParams, mols_to_grid_image_use_svg_with_params};
use rdkit_cxx::chem::rd_depictor::set_prefer_coord_gen;
use rdkit_cxx::chem::rd_fmcs::find_mcs;
use rdkit_cxx::chem::{Mol, mol_from_smarts, mol_from_smiles};

/// The atoms of each molecule outside their maximum common substructure.
fn view_difference(mol1: &Mol, mol2: &Mol) -> rdkit_cxx::Result<(Vec<u32>, Vec<u32>)> {
    let mcs = find_mcs(&[mol1.clone(), mol2.clone()])?;
    let mcs_mol = mol_from_smarts(mcs.smarts_string())?;
    let outside = |mol: &Mol| -> rdkit_cxx::Result<Vec<u32>> {
        let matched = mol.get_substruct_match(&mcs_mol)?;
        Ok(mol
            .get_atoms()
            .map(|atom| atom.get_idx())
            .filter(|idx| !matched.contains(idx))
            .collect())
    };
    Ok((outside(mol1)?, outside(mol2)?))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    set_prefer_coord_gen(true);
    let console = MolDrawOptions {
        min_font_size: 20,
        ..Default::default()
    };

    let mol1 = mol_from_smiles("FC1=CC=C2C(=C1)C=NN2")?;
    let mol2 = mol_from_smiles("CCC1=C2NN=CC2=CC(Cl)=C1")?;
    let grid = MolsToGridImageParams {
        draw_options: Some(console),
        ..Default::default()
    };
    let mols = [mol1.clone(), mol2.clone()];
    let svg = mols_to_grid_image_use_svg_with_params(&mols, &grid)?;
    write_image("highlight_differences_plain.svg", &svg)?;

    let (target_atm1, target_atm2) = view_difference(&mol1, &mol2)?;
    println!("{target_atm1:?} {target_atm2:?}");
    assert_eq!(target_atm1, [0]);
    assert_eq!(target_atm2, [0, 1, 10]);
    let highlighted = MolsToGridImageParams {
        highlight_atom_lists: Some(vec![target_atm1, target_atm2]),
        ..grid
    };
    let svg = mols_to_grid_image_use_svg_with_params(&mols, &highlighted)?;
    write_image("highlight_differences.svg", &svg)?;
    assert_eq!(svg.matches("#FF7F7F").count(), 10);
    Ok(())
}
