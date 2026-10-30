// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Without Implicit Hydrogens (RDKitCB_17): draw atoms without their
//! implicit hydrogens.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#without-implicit-hydrogens),
//! by Greg Landrum.

mod common;

use common::{MOL_SIZE, count_class, show, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::mol_from_smiles;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let console = MolDrawOptions::default();
    // Each atom's label characters.
    let label_lengths = |svg: &str| -> Vec<usize> {
        (0..5)
            .map(|idx| count_class(svg, &format!("atom-{idx}")))
            .collect()
    };

    let mut m = mol_from_smiles("[Pt](Cl)(Cl)(N)N")?;
    // The nitrogens are labeled H3N and NH3.
    let svg = show(&m, MOL_SIZE, &console)?;
    write_image("without_implicit_hs_default.svg", &svg)?;
    assert_eq!(label_lengths(&svg), [2, 2, 2, 3, 3]);

    // An atom label replaces the symbol and its hydrogens.
    for idx in 0..m.get_num_atoms() {
        let mut atom = m.get_atom_with_idx_mut(idx)?;
        let symbol = atom.get_symbol()?;
        atom.set_prop("atomLabel", &symbol);
    }
    let svg = show(&m, MOL_SIZE, &console)?;
    write_image("without_implicit_hs.svg", &svg)?;
    assert_eq!(label_lengths(&svg), [2, 2, 2, 1, 1]);
    Ok(())
}
