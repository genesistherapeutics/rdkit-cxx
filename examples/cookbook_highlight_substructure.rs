// Copyright (c) 2026 Steven Kearnes
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
//! the matched atoms highlighted, and the bonds between them
//! (`IPythonConsole` keeps the matched atoms in `__sssAtoms` and passes no
//! bond highlights, so the drawer highlights the bonds that join two
//! highlighted atoms). The recipe then sets the atoms to highlight by hand.
//! `mol_to_svg_with_params` takes the bonds explicitly, so these pass the
//! bonds between the atoms too.

mod common;

use common::{MOL_SIZE, write_image};
use rdkit_cxx::chem::draw::{MolToSVGParams, mol_to_svg_with_params};
use rdkit_cxx::chem::{Mol, mol_from_smarts, mol_from_smiles};

/// The highlight color, `MolDrawOptions::highlight_colour`'s default.
const HIGHLIGHT: &str = "#FF7F7F";

/// Draws `mol` as Jupyter does with `atoms` highlighted.
fn highlighted(mol: &Mol, atoms: Vec<u32>) -> rdkit_cxx::Result<String> {
    let bonds = mol
        .get_bonds()
        .filter(|b| {
            atoms.contains(&b.get_begin_atom_idx()) && atoms.contains(&b.get_end_atom_idx())
        })
        .map(|b| b.get_idx())
        .collect();
    let params = MolToSVGParams {
        size: MOL_SIZE,
        highlight_atoms: atoms,
        highlight_bonds: bonds,
        ..Default::default()
    };
    mol_to_svg_with_params(mol, &params)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let m = mol_from_smiles("c1cc(C(=O)O)c(OC(=O)C)cc1")?;
    let substructure = mol_from_smarts("C(=O)O")?;
    let matches = m.get_substruct_matches(&substructure)?;
    println!("{matches:?}");
    assert_eq!(matches, [[3, 4, 5], [8, 9, 7]]);

    let svg = highlighted(&m, matches.concat())?;
    write_image("highlight_substructure.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), 20);

    // You can also set the atoms to highlight by hand:
    let svg = highlighted(&m, vec![0, 1, 2, 6, 11, 12])?;
    write_image("highlight_substructure_by_hand.svg", &svg)?;
    assert_eq!(svg.matches(HIGHLIGHT).count(), 24);
    Ok(())
}
