// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Writing Molecules" chapter.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    Mol, SmilesWriteParams, kekulize, mol_from_smiles, mol_to_smiles, mol_to_smiles_with_params,
};

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[
    ("RDKitCB_4", "Kekule SMILES", kekule_smiles),
    ("RDKitCB_5", "Isomeric SMILES without isotopes", smiles_without_isotopes),
];

/// Kekulizes a molecule and writes Kekulé SMILES. By Paolo Tosco.
fn kekule_smiles() -> Result<()> {
    let mut mol = mol_from_smiles("CN1C(NC2=NC=CC=C2)=CC=C1")?;
    println!("{}", mol_to_smiles(&mol)?);
    assert_eq!(mol_to_smiles(&mol)?, "Cn1cccc1Nc1ccccn1");

    kekulize(&mut mol)?;
    let params = SmilesWriteParams {
        do_kekule: true,
        ..Default::default()
    };
    let smiles = mol_to_smiles_with_params(&mol, &params)?;
    println!("{smiles}");
    assert_eq!(smiles, "CN1C=CC=C1NC1=CC=CC=N1");
    Ok(())
}

/// Writes isomeric SMILES with the isotopes cleared, then restores them.
fn mol_without_isotopes_to_smiles(mol: &mut Mol) -> Result<String> {
    let isotopes: Vec<u32> = mol.get_atoms().map(|atom| atom.get_isotope()).collect();
    for (idx, &isotope) in (0..).zip(&isotopes) {
        if isotope != 0 {
            mol.get_atom_with_idx_mut(idx)?.set_isotope(0);
        }
    }
    let smiles = mol_to_smiles(mol)?;
    for (idx, &isotope) in (0..).zip(&isotopes) {
        if isotope != 0 {
            mol.get_atom_with_idx_mut(idx)?.set_isotope(isotope);
        }
    }
    Ok(smiles)
}

/// Writes isomeric SMILES that keep stereochemistry but drop isotopes. By
/// Andrew Dalke. Isotopic hydrogens, and stereocenters that exist only
/// because of isotopes, need more care (see the Cookbook).
fn smiles_without_isotopes() -> Result<()> {
    let mut mol = mol_from_smiles("[19F][13C@H]([16OH])[35Cl]")?;
    let smiles = mol_without_isotopes_to_smiles(&mut mol)?;
    println!("{smiles}");
    assert_eq!(smiles, "O[C@@H](F)Cl");
    assert_eq!(mol_to_smiles(&mol)?, "[16OH][13C@@H]([19F])[35Cl]");
    Ok(())
}
