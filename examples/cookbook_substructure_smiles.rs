// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Returning Substructure Matches as SMILES (RDKitCB_18): write the SMILES of
//! a matched substructure.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#returning-substructure-matches-as-smiles),
//! by Andrew Dalke.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    Mol, MolFragmentToSmilesParams, mol_fragment_to_smiles, mol_fragment_to_smiles_with_params,
    mol_from_smarts, mol_from_smiles,
};

/// The molecule's bonds that the query's bonds match.
fn get_match_bond_indices(query: &Mol, mol: &Mol, match_atom_indices: &[u32]) -> Result<Vec<u32>> {
    let mut bond_indices = Vec::new();
    for query_bond in query.get_bonds() {
        let atom_index1 = match_atom_indices[query_bond.get_begin_atom_idx() as usize];
        let atom_index2 = match_atom_indices[query_bond.get_end_atom_idx() as usize];
        let bond = mol
            .get_bond_between_atoms(atom_index1, atom_index2)?
            .expect("the match has a bond for each query bond");
        bond_indices.push(bond.get_idx());
    }
    Ok(bond_indices)
}

fn main() -> Result<()> {
    let pat = mol_from_smarts("[NX1]#[CX2]")?; // matches nitrile
    let mol = mol_from_smiles("CCCC#N")?; // Butyronitrile
    let atom_indices = mol.get_substruct_match(&pat)?;
    println!("{atom_indices:?}");
    assert_eq!(atom_indices, [4, 3]);

    let nitrile = mol_fragment_to_smiles(&mol, &atom_indices)?; // returns the nitrile
    println!("{nitrile}");
    assert_eq!(nitrile, "C#N");

    // If only the atom indices are given, mol_fragment_to_smiles includes
    // every bond between those atoms, even bonds the SMARTS does not match:
    let pat = mol_from_smarts("*~*~*~*")?; // match 4 linear atoms
    let mol = mol_from_smiles("C1CCC1")?; // ring of size 4
    let atom_indices = mol.get_substruct_match(&pat)?;
    println!("{atom_indices:?}");
    assert_eq!(atom_indices, [0, 1, 2, 3]);

    let ring = mol_fragment_to_smiles(&mol, &atom_indices)?; // returns the ring
    println!("{ring}");
    assert_eq!(ring, "C1CCC1");

    // If this matters, pass the bond indices too, which the query graph's
    // bonds give in the molecule graph:
    let bond_indices = get_match_bond_indices(&pat, &mol, &atom_indices)?;
    println!("{bond_indices:?}");
    assert_eq!(bond_indices, [0, 1, 2]);

    let params = MolFragmentToSmilesParams {
        bonds_to_use: Some(bond_indices),
        ..Default::default()
    };
    let chain = mol_fragment_to_smiles_with_params(&mol, &atom_indices, &params)?;
    println!("{chain}");
    assert_eq!(chain, "CCCC");
    Ok(())
}
