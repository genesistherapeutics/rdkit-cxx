// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Create Fragments (RDKitCB_7): fragment molecules on bonds.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#create-fragments),
//! by Paolo Tosco. The recipe draws each result; this prints their SMILES.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    Mol, RWMol, fragment_on_bonds, fragment_on_some_bonds, get_mol_frags_as_mols, mol_from_smiles,
    mol_to_smiles,
};

fn smiles(mols: &[Mol]) -> Result<Vec<String>> {
    mols.iter().map(mol_to_smiles).collect()
}

fn main() -> Result<()> {
    // Explicit bonds in the SMILES make the bond indices easy to follow:
    let mol = mol_from_smiles("O-C-C-C-C-N")?;

    // fragment_on_bonds fragments all the bonds at once and returns a single
    // molecule with all the cuts applied. By default empty valences are
    // filled with dummy atoms, each labeled with the index of the atom it
    // replaces as its isotope (atom 0's label, 0, writes no isotope):
    let mol1_f = fragment_on_bonds(&mol, &[0, 2, 4])?;
    println!("{}", mol_to_smiles(&mol1_f)?);
    assert_eq!(mol_to_smiles(&mol1_f)?, "*CC[3*].[1*]O.[2*]CC[5*].[4*]N");

    // This molecule can be split into individual fragments with
    // get_mol_frags_as_mols:
    let frags = smiles(&get_mol_frags_as_mols(&mol1_f)?)?;
    println!("{frags:?}");
    assert_eq!(frags, ["[1*]O", "*CC[3*]", "[2*]CC[5*]", "[4*]N"]);

    // fragment_on_some_bonds fragments one bond at a time (Python's
    // numToBreak defaults to 1) and returns one molecule per cut:
    let mol2_f = smiles(&fragment_on_some_bonds(&mol, &[0, 2, 4])?)?;
    println!("{mol2_f:?}");
    assert_eq!(
        mol2_f,
        ["*CCCCN.[1*]O", "[2*]CCN.[3*]CCO", "[4*]N.[5*]CCCCO"]
    );

    // Finally, you can cut bonds by hand with RWMol::remove_bond. Python's
    // `with Chem.RWMol(mol) as rwmol:` makes the edits a batch, so the bond
    // indices hold until the batch is committed:
    let mut rwmol = RWMol::from(mol.clone());
    rwmol.begin_batch_edit()?;
    for b_idx in [0, 2, 4] {
        let b = rwmol.get_bond_with_idx(b_idx)?;
        let (begin, end) = (b.get_begin_atom_idx(), b.get_end_atom_idx());
        rwmol.remove_bond(begin, end)?;
    }
    rwmol.commit_batch_edit()?;
    // get_mol_frags_as_mols then returns sanitized fragments whose empty
    // valences are filled with implicit hydrogens:
    let frags = smiles(&get_mol_frags_as_mols(&rwmol)?)?;
    println!("{frags:?}");
    assert_eq!(frags, ["O", "CC", "CC", "N"]);
    Ok(())
}
