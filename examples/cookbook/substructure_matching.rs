// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Substructure Matching" chapter.

use std::collections::BTreeSet;

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    Mol, MolFragmentToSmilesParams, get_mol_frags, mol_fragment_to_smiles,
    mol_fragment_to_smiles_with_params, mol_from_smarts, mol_from_smiles,
};

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[
    ("RDKitCB_10", "Functional Group with SMARTS queries", functional_group),
    ("RDKitCB_13", "Macrocycles with SMARTS queries", macrocycles),
    ("RDKitCB_18", "Returning Substructure Matches as SMILES", matches_as_smiles),
    ("RDKitCB_20", "Within the Same Fragment", within_the_same_fragment),
];

/// Matches primary and secondary alcohols. By Paolo Tosco.
fn functional_group() -> Result<()> {
    let sucrose = mol_from_smiles(
        "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O",
    )?;
    let primary_alcohol = mol_from_smarts("[CH2][OH1]")?;
    let matches = sucrose.get_substruct_matches(&primary_alcohol)?;
    println!("{matches:?}");
    assert_eq!(matches, [[0, 22], [13, 14], [17, 18]]);

    let secondary_alcohol = mol_from_smarts("[CH1][OH1]")?;
    let matches = sucrose.get_substruct_matches(&secondary_alcohol)?;
    println!("{matches:?}");
    assert_eq!(matches, [[2, 21], [3, 20], [4, 19], [9, 16], [10, 15]]);
    Ok(())
}

/// Matches the atoms of rings larger than 12 atoms. By Ivan Tubert-Brohman,
/// David Cosgrove, and Vincent Scalfani.
fn macrocycles() -> Result<()> {
    let erythromycin = mol_from_smiles(
        "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)\
         O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)\
         (C)O)C)C)O)(C)O",
    )?;
    // `r{12-}` is RDKit's SMARTS extension for a ring size range.
    let macro_ring = mol_from_smarts("[r{12-}]")?;
    let matches = erythromycin.get_substruct_matches(&macro_ring)?;
    println!("{matches:?}");
    let atoms: Vec<u32> = matches.concat();
    assert_eq!(atoms, [2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 17]);
    Ok(())
}

/// The bonds of `mol` that a match of `query` covers.
fn get_match_bond_indices(query: &Mol, mol: &Mol, match_atom_indices: &[u32]) -> Result<Vec<u32>> {
    let mut bond_indices = Vec::new();
    for query_bond in query.get_bonds() {
        let atom_index1 = match_atom_indices[query_bond.get_begin_atom_idx() as usize];
        let atom_index2 = match_atom_indices[query_bond.get_end_atom_idx() as usize];
        let bond = mol
            .get_bond_between_atoms(atom_index1, atom_index2)?
            .expect("a match maps query bonds to bonds");
        bond_indices.push(bond.get_idx());
    }
    Ok(bond_indices)
}

/// Writes the matched part of a molecule as SMILES. By Andrew Dalke.
fn matches_as_smiles() -> Result<()> {
    let pat = mol_from_smarts("[NX1]#[CX2]")?; // a nitrile
    let mol = mol_from_smiles("CCCC#N")?; // butyronitrile
    let atom_indices = mol.get_substruct_match(&pat)?;
    println!("{atom_indices:?}");
    assert_eq!(atom_indices, [4, 3]);
    let smiles = mol_fragment_to_smiles(&mol, &atom_indices)?;
    println!("{smiles}");
    assert_eq!(smiles, "C#N");

    // Given only atoms, MolFragmentToSmiles keeps every bond between them,
    // including bonds the pattern did not match.
    let pat = mol_from_smarts("*~*~*~*")?; // four atoms in a row
    let mol = mol_from_smiles("C1CCC1")?; // a four-membered ring
    let atom_indices = mol.get_substruct_match(&pat)?;
    println!("{atom_indices:?}");
    assert_eq!(atom_indices, [0, 1, 2, 3]);
    let smiles = mol_fragment_to_smiles(&mol, &atom_indices)?;
    println!("{smiles}");
    assert_eq!(smiles, "C1CCC1");

    // The query's bonds name the molecule's bonds that the match covers.
    let bond_indices = get_match_bond_indices(&pat, &mol, &atom_indices)?;
    println!("{bond_indices:?}");
    assert_eq!(bond_indices, [0, 1, 2]);
    let smiles = mol_fragment_to_smiles_with_params(
        &mol,
        &atom_indices,
        &MolFragmentToSmilesParams {
            bonds_to_use: Some(bond_indices),
            ..Default::default()
        },
    )?;
    println!("{smiles}");
    assert_eq!(smiles, "CCCC");
    Ok(())
}

/// The first match of `p` that lies within one fragment of `m`.
fn fragsearch(m: &Mol, p: &Mol) -> Result<Option<BTreeSet<u32>>> {
    let matches: Vec<BTreeSet<u32>> = m
        .get_substruct_matches(p)?
        .into_iter()
        .map(BTreeSet::from_iter)
        .collect();
    for frag in get_mol_frags(m) {
        let frag = BTreeSet::from_iter(frag);
        if let Some(found) = matches.iter().find(|found| found.is_subset(&frag)) {
            return Ok(Some(found.clone()));
        }
    }
    Ok(None)
}

/// Matches a disconnected pattern only within one fragment. By Greg Landrum.
fn within_the_same_fragment() -> Result<()> {
    let p = mol_from_smarts("O.N")?;
    let m1 = mol_from_smiles("OCCCN.CCC")?;
    let m2 = mol_from_smiles("OCCC.CCCN")?;
    // The dot in a SMARTS pattern lets its parts match anywhere.
    let has_match = [m1.has_substruct_match(&p)?, m2.has_substruct_match(&p)?];
    println!("{has_match:?}");
    assert_eq!(has_match, [true, true]);

    let found = [fragsearch(&m1, &p)?, fragsearch(&m2, &p)?];
    println!("{found:?}");
    assert_eq!(found, [Some(BTreeSet::from([0, 4])), None]);
    Ok(())
}
