// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Include an Atom Index (RDKitCB_0): draw a molecule with its atom indices.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#include-an-atom-index),
//! by Takayuki Serizawa.

mod common;

use std::collections::BTreeSet;

use common::{MOL_SIZE, count_class, labelled_atoms, show, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::{Mol, mol_from_smiles};

/// Labels each atom with its index, as an atom map number.
fn mol_with_atom_index(mol: &mut Mol) -> rdkit_cxx::Result<()> {
    for idx in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(idx)?.set_atom_map_num(idx as i32);
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let console = MolDrawOptions::default();

    // Test in a kinase inhibitor
    let smiles = "C1CC2=C3C(=CC=C2)C(=CN3C1)[C@H]4[C@@H](C(=O)NC4=O)C5=CNC6=CC=CC=C65";
    let mut mol = mol_from_smiles(smiles)?;
    // Default: only the heteroatoms have labels.
    let svg = show(&mol, MOL_SIZE, &console)?;
    write_image("atom_index_default.svg", &svg)?;
    assert_eq!(labelled_atoms(&svg), BTreeSet::from([10, 15, 16, 18, 21]));

    // With atom index: each atom is labelled with its map number, but atom
    // 0, whose map number is 0 (no map number), is not.
    mol_with_atom_index(&mut mol)?;
    let svg = show(&mol, MOL_SIZE, &console)?;
    write_image("atom_index_map_numbers.svg", &svg)?;
    assert_eq!(labelled_atoms(&svg), (1..mol.get_num_atoms()).collect());

    // A simpler way to add atom indices is to set the drawing options. The
    // indices are then drawn near the atoms rather than at them.
    let console = MolDrawOptions {
        add_atom_indices: true,
        ..console
    };
    let mol = mol_from_smiles(smiles)?;
    let svg = show(&mol, (300, 300), &console)?;
    write_image("atom_index_options.svg", &svg)?;
    // A note for each character of each index.
    let digits: usize = (0..mol.get_num_atoms()).map(|i| i.to_string().len()).sum();
    assert_eq!(count_class(&svg, "note"), digits);
    Ok(())
}
