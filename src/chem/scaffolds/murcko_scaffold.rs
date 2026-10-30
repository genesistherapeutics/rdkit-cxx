// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Scaffolds.MurckoScaffold`, a Python module: Murcko
//! scaffolds, the ring systems of a molecule and the linkers between them.
//!
//! Python builds its patterns and reactions once, on import; each function
//! here that returns one makes a new one, since neither a `Mol` nor a
//! `ChemicalReaction` can be shared between threads (CONVENTIONS.md S1, S4).
//! The scaffold itself comes from RDKit's `MurckoDecompose`, as in Python;
//! the patterns are those of Python's private fallback, written in Python, which
//! Python never calls.

use crate::Result;
use crate::chem::rd_chem_reactions::{ChemicalReaction, reaction_from_smarts};
use crate::chem::{
    BondType, ChiralType, Mol, SmilesWriteParams, get_symm_sssr, mol_from_smarts, mol_from_smiles,
    mol_to_smiles_with_params, murcko_decompose, remove_hs,
};
use crate::error::Error;

fn reaction(smarts: &str) -> ChemicalReaction {
    reaction_from_smarts(smarts).expect("the MurckoScaffold reaction SMARTS parse")
}

/// The reactions that replace a side chain's last atom with a hydrogen.
/// Mirrors `MurckoScaffold.murckoTransforms`.
pub fn murcko_transforms() -> Vec<ChemicalReaction> {
    [
        "[*:1]-[!#1;D1]>>[*:1][H]",
        "[*:1]-[!#1;D2]#[AD1]>>[*:1][H]",
        "[*:1]-[!#1;D2]=[AD1]>>[*:1][H]",
        "[*:1]-[!#1;D3](=[AD1])=[AD1]>>[*:1][H]",
    ]
    .into_iter()
    .map(reaction)
    .collect()
}

const MURCKO_PATTS: [&str; 3] = [
    "[!#1;D3;$([D3]-[!#1])](=[AD1])=[AD1]",
    "[!#1;D2;$([D2]-[!#1])]=,#[AD1]",
    "[!#1;D1;$([D1]-[!#1;!n])]",
];

/// The patterns of the atoms that are not part of a scaffold: an atom's
/// doubly bonded terminal pair, a doubly or triply bonded terminal atom,
/// and a terminal atom not on an aromatic nitrogen. Mirrors
/// `MurckoScaffold.murckoPatts`.
pub fn murcko_patts() -> Vec<Mol> {
    MURCKO_PATTS
        .into_iter()
        .map(|sma| mol_from_smarts(sma).expect("the MurckoScaffold SMARTS parse"))
        .collect()
}

/// One pattern for any of [`murcko_patts`]. Mirrors
/// `MurckoScaffold.murckoQ`.
pub fn murcko_q() -> Mol {
    let alternatives: Vec<String> = MURCKO_PATTS.iter().map(|x| format!("$({x})")).collect();
    mol_from_smarts(&format!("[{}]", alternatives.join(",")))
        .expect("the MurckoScaffold SMARTS parse")
}

/// The reaction that gives an aromatic nitrogen a hydrogen in place of its
/// terminal neighbor. Mirrors `MurckoScaffold.aromaticNTransform`.
pub fn aromatic_n_transform() -> ChemicalReaction {
    reaction("[n:1]-[D1]>>[nH:1]")
}

/// The scaffold made generic: every atom but hydrogen a carbon, every bond
/// single, with no aromaticity, isotopes, charges, chirality, or explicit
/// hydrogens; hydrogens are then removed. `Error::AtomValence` (or another
/// sanitization variant) where removing hydrogens fails, as in Python.
/// Mirrors `MurckoScaffold.MakeScaffoldGeneric`.
pub fn make_scaffold_generic(mol: &Mol) -> Result<Mol> {
    let mut res = mol.clone();
    for idx in 0..res.get_num_atoms() {
        let mut atom = res.get_atom_with_idx_mut(idx)?;
        if atom.get_atomic_num() != 1 {
            atom.set_atomic_num(6);
        }
        atom.set_is_aromatic(false);
        atom.set_isotope(0);
        atom.set_formal_charge(0);
        atom.set_chiral_tag(ChiralType::ChiUnspecified);
        atom.set_no_implicit(false);
        atom.set_num_explicit_hs(0);
    }
    for idx in 0..res.get_num_bonds() {
        let mut bond = res.get_bond_with_idx_mut(idx)?;
        bond.set_bond_type(BondType::Single);
        bond.set_is_aromatic(false);
    }
    remove_hs(&res)
}

/// The molecule's Murcko scaffold, with its property cache and rings
/// computed afresh. Mirrors `MurckoScaffold.GetScaffoldForMol`.
pub fn get_scaffold_for_mol(mol: &Mol) -> Result<Mol> {
    let mut res = murcko_decompose(mol)?;
    res.clear_computed_props();
    res.update_property_cache()?;
    get_symm_sssr(&mut res)?;
    Ok(res)
}

/// The SMILES of the scaffold of the molecule a SMILES describes, without
/// stereochemistry. Empty SMILES is `Error::Value`, as in Python, and
/// SMILES RDKit cannot parse is `Error::Parse`, where Python raises
/// `ValueError`. Mirrors `MurckoScaffold.MurckoScaffoldSmiles(smiles)`.
pub fn murcko_scaffold_smiles(smiles: &str) -> Result<String> {
    murcko_scaffold_smiles_with_include_chirality(smiles, false)
}

/// As [`murcko_scaffold_smiles`], with stereochemistry when
/// `include_chirality` is set. Mirrors
/// `MurckoScaffold.MurckoScaffoldSmiles(smiles, includeChirality=...)`.
pub fn murcko_scaffold_smiles_with_include_chirality(
    smiles: &str,
    include_chirality: bool,
) -> Result<String> {
    if smiles.is_empty() {
        return Err(Error::Value {
            message: "No molecule provided".to_owned(),
        });
    }
    murcko_scaffold_smiles_mol_with_include_chirality(&mol_from_smiles(smiles)?, include_chirality)
}

/// The SMILES of the molecule's scaffold, without stereochemistry. Mirrors
/// `MurckoScaffold.MurckoScaffoldSmiles(mol=mol)`.
pub fn murcko_scaffold_smiles_mol(mol: &Mol) -> Result<String> {
    murcko_scaffold_smiles_mol_with_include_chirality(mol, false)
}

/// As [`murcko_scaffold_smiles_mol`], with stereochemistry when
/// `include_chirality` is set. Mirrors
/// `MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=...)`.
pub fn murcko_scaffold_smiles_mol_with_include_chirality(
    mol: &Mol,
    include_chirality: bool,
) -> Result<String> {
    let scaffold = get_scaffold_for_mol(mol)?;
    mol_to_smiles_with_params(
        &scaffold,
        &SmilesWriteParams {
            do_isomeric_smiles: include_chirality,
            ..Default::default()
        },
    )
}

/// Mirrors `MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smiles)`, which is
/// [`murcko_scaffold_smiles`].
pub fn murcko_scaffold_smiles_from_smiles(smiles: &str) -> Result<String> {
    murcko_scaffold_smiles(smiles)
}

/// Mirrors `MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smiles,
/// includeChirality)`.
pub fn murcko_scaffold_smiles_from_smiles_with_include_chirality(
    smiles: &str,
    include_chirality: bool,
) -> Result<String> {
    murcko_scaffold_smiles_with_include_chirality(smiles, include_chirality)
}
