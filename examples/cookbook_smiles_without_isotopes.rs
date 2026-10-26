// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Isomeric SMILES without isotopes (RDKitCB_5): write isomeric SMILES with
//! stereochemistry but no isotopes.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#isomeric-smiles-without-isotopes),
//! by Andrew Dalke. The method has two limitations, isotopic hydrogens and
//! stereochemistry that needs recalculating; the source discussion explains
//! them.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{Mol, mol_from_smiles, mol_to_smiles};

fn mol_without_isotopes_to_smiles(mol: &mut Mol) -> Result<String> {
    let isotopes: Vec<u32> = mol.get_atoms().map(|atom| atom.get_isotope()).collect();
    for (idx, &isotope) in isotopes.iter().enumerate() {
        if isotope != 0 {
            mol.get_atom_with_idx_mut(idx as u32)?.set_isotope(0);
        }
    }
    let smiles = mol_to_smiles(mol)?;
    // Restore the original isotopes.
    for (idx, &isotope) in isotopes.iter().enumerate() {
        if isotope != 0 {
            mol.get_atom_with_idx_mut(idx as u32)?.set_isotope(isotope);
        }
    }
    Ok(smiles)
}

fn main() -> Result<()> {
    let mut mol = mol_from_smiles("[19F][13C@H]([16OH])[35Cl]")?;
    let smiles = mol_without_isotopes_to_smiles(&mut mol)?;
    println!("{smiles}");
    assert_eq!(smiles, "O[C@@H](F)Cl");
    assert_eq!(mol_to_smiles(&mol)?, "[16OH][13C@@H]([19F])[35Cl]");
    Ok(())
}
