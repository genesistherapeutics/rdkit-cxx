// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Include Stereo Annotations (RDKitCB_32): draw molecules with their CIP
//! labels.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#include-stereo-annotations),
//! by Greg Landrum.

mod common;

use common::{count_class, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::draw::{MolsToGridImageParams, mols_to_grid_image_use_svg_with_params};
use rdkit_cxx::chem::{Mol, assign_cip_labels, mol_from_smiles};

/// Each atom's CIP label, by index.
fn cip_codes(mol: &Mol) -> rdkit_cxx::Result<Vec<(u32, String)>> {
    mol.get_atoms()
        .filter(|atom| atom.has_prop("_CIPCode"))
        .map(|atom| Ok((atom.get_idx(), atom.get_prop("_CIPCode")?)))
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let console = MolDrawOptions {
        add_atom_indices: false,
        add_stereo_annotation: true,
        ..Default::default()
    };
    let grid = MolsToGridImageParams {
        sub_img_size: (250, 250),
        draw_options: Some(console),
        ..Default::default()
    };
    let label = |codes: &[(u32, &str)]| -> Vec<(u32, String)> {
        codes
            .iter()
            .map(|&(idx, code)| (idx, code.to_owned()))
            .collect()
    };

    // The default labels come from the legacy FindMolChiralCenters() code.
    let mut m1 = mol_from_smiles("C1CC1[C@H](F)C1CCC1")?;
    let mut m2 = mol_from_smiles("F[C@H]1CC[C@H](O)CC1")?;
    let mols = [m1.clone(), m2.clone()];
    let svg = mols_to_grid_image_use_svg_with_params(&mols, &grid)?;
    write_image("stereo_annotations_legacy.svg", &svg)?;
    assert_eq!(cip_codes(&m1)?, label(&[(3, "R")]));
    assert_eq!(cip_codes(&m2)?, label(&[]));
    // One note element for each character of "(R)": RDKit 2024.09 draws
    // CIP labels as atom notes.
    assert_eq!(count_class(&svg, "note"), 3);

    // The new stereochemistry code (2020.09) gives more accurate labels.
    assign_cip_labels(&mut m1)?;
    assign_cip_labels(&mut m2)?;
    let svg = mols_to_grid_image_use_svg_with_params(&[m1.clone(), m2.clone()], &grid)?;
    write_image("stereo_annotations.svg", &svg)?;
    println!("{:?} {:?}", cip_codes(&m1)?, cip_codes(&m2)?);
    assert_eq!(cip_codes(&m1)?, label(&[(3, "S")]));
    assert_eq!(cip_codes(&m2)?, label(&[(1, "r"), (4, "r")]));
    // "(S)", "(r)", and "(r)".
    assert_eq!(count_class(&svg, "note"), 9);
    Ok(())
}
