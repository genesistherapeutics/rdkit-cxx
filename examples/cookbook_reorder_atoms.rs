// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Reorder Atoms (RDKitCB_28): renumber atoms in a canonical order that does
//! not depend on the input.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#reorder-atoms),
//! by Jeffrey Van Santen and Paolo Tosco. The recipe draws the molecules with
//! their atom indices; this prints them as atom map numbers.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    Mol, SmilesWriteParams, canonical_rank_atoms, mol_from_smiles, mol_to_smiles,
    mol_to_smiles_with_params, renumber_atoms,
};

/// The atoms in canonical order: the atom with rank 0 first.
fn canonical_order(mol: &Mol) -> Result<Vec<u32>> {
    let ranks = canonical_rank_atoms(mol)?;
    let mut order: Vec<u32> = (0..ranks.len() as u32).collect();
    order.sort_by_key(|&i| ranks[i as usize]);
    Ok(order)
}

/// Sets each atom's map number to its index.
fn add_atom_indices(mol: &mut Mol) -> Result<()> {
    for i in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(i)?.set_atom_map_num(i as i32);
    }
    Ok(())
}

/// The SMILES with atom map numbers, written in input order.
fn indexed_smiles(mol: &Mol) -> Result<String> {
    let params = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    mol_to_smiles_with_params(mol, &params)
}

fn main() -> Result<()> {
    let mut m = mol_from_smiles("c1([C@H](C)CC)cccc2ccccc12")?;
    let mut m1 = mol_from_smiles("c12ccccc1c(ccc2)[C@H](C)CC")?;
    println!("{}", mol_to_smiles(&m)? == mol_to_smiles(&m1)?);
    assert_eq!(mol_to_smiles(&m)?, mol_to_smiles(&m1)?);

    // The atoms are numbered differently:
    let m_neworder = canonical_order(&m)?;
    let m1_neworder = canonical_order(&m1)?;
    println!("{}", m_neworder == m1_neworder);
    assert_ne!(m_neworder, m1_neworder);

    add_atom_indices(&mut m)?;
    add_atom_indices(&mut m1)?;
    println!("{}", indexed_smiles(&m)?);
    println!("{}", indexed_smiles(&m1)?);

    // Renumbered in the same canonical order, they match atom for atom:
    let mut m_renum = renumber_atoms(&m, &m_neworder)?;
    let mut m1_renum = renumber_atoms(&m1, &m1_neworder)?;
    add_atom_indices(&mut m_renum)?;
    add_atom_indices(&mut m1_renum)?;
    println!("{}", indexed_smiles(&m_renum)?);
    println!("{}", indexed_smiles(&m1_renum)?);
    assert_eq!(indexed_smiles(&m_renum)?, indexed_smiles(&m1_renum)?);
    Ok(())
}
