// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Include a Calculation (RDKitCB_23): draw a molecule with a calculated
//! value (here, its Gasteiger charges) by each atom.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#include-a-calculation),
//! by Greg Landrum.

mod common;

use common::{count_class, show, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_partial_charges::compute_gasteiger_charges;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let console = MolDrawOptions::default();
    let mol_size = (250, 250);

    let mut m = mol_from_smiles("c1ncncc1C(=O)[O-]")?;
    compute_gasteiger_charges(&mut m)?;
    let svg = show(&m, mol_size, &console)?;
    write_image("calculation_plain.svg", &svg)?;
    assert_eq!(count_class(&svg, "note"), 0);

    // Each atom's charge as its note.
    let mut m2 = m.clone();
    let mut lbls = Vec::new();
    for idx in 0..m2.get_num_atoms() {
        let mut at = m2.get_atom_with_idx_mut(idx)?;
        let lbl = format!("{:.2}", at.get_double_prop("_GasteigerCharge")?);
        at.set_prop("atomNote", &lbl);
        lbls.push(lbl);
    }
    println!("{lbls:?}");
    assert_eq!(
        lbls,
        [
            "0.04", "-0.24", "0.11", "-0.24", "0.04", "0.03", "0.07", "-0.54", "-0.54"
        ]
    );
    let svg = show(&m2, mol_size, &console)?;
    write_image("calculation.svg", &svg)?;
    // A note for each character of each label.
    assert_eq!(
        count_class(&svg, "note"),
        lbls.iter().map(String::len).sum::<usize>()
    );
    Ok(())
}
