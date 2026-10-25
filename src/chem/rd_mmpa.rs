// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMMPA`: the fragmentation behind matched molecular
//! pair analysis.
//!
//! Each fragmentation is a core and its side chains, with dummy atoms
//! numbered by cut (atom map numbers). A single cut has no core. Python's
//! `resultsAsMols=False` gives SMILES instead (A4): the `_as_smiles`
//! forms, with an empty core for none, as in Python.

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdmolfiles::mol_to_smiles;
use crate::error::expect_ok;
use crate::ffi::mmpa::ffi;
use crate::ffi::molops::ffi::{MolList, mol_list_size, mol_list_take};

// Python's default pattern: acyclic single bonds from carbon.
const DEFAULT_PATTERN: &str = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]";

fn into_pairs(mut list: UniquePtr<MolList>) -> Vec<(Option<Mol>, Mol)> {
    let size = expect_ok(mol_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    (0..size / 2)
        .map(|i| {
            let core = expect_ok(mol_list_take(list.pin_mut(), 2 * i));
            let side_chains = expect_ok(mol_list_take(list.pin_mut(), 2 * i + 1));
            (
                (!core.is_null()).then(|| Mol::from_ffi(core)),
                Mol::from_ffi(side_chains),
            )
        })
        .collect()
}

fn as_smiles(pairs: Vec<(Option<Mol>, Mol)>) -> Result<Vec<(String, String)>> {
    pairs
        .iter()
        .map(|(core, side_chains)| {
            let core = match core {
                Some(core) => mol_to_smiles(core)?,
                None => String::new(),
            };
            Ok((core, mol_to_smiles(side_chains)?))
        })
        .collect()
}

/// Parameters for `fragment_mol_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentMolParams {
    /// The most bonds cut at once.
    pub max_cuts: u32,
    /// Nothing is cut when the pattern matches more bonds than this.
    pub max_cut_bonds: u32,
    /// SMARTS for the bonds to cut; its first two atoms are the bond.
    pub pattern: String,
}

impl Default for FragmentMolParams {
    fn default() -> FragmentMolParams {
        FragmentMolParams {
            max_cuts: 3,
            max_cut_bonds: 20,
            pattern: DEFAULT_PATTERN.to_string(),
        }
    }
}

/// Every way of cutting one to `max_cuts` of the bonds the pattern
/// matches; empty when it matches none or more than `max_cut_bonds`.
/// Mirrors `rdMMPA.FragmentMol(mol, maxCuts, maxCutBonds, pattern)`.
/// `Error::Parse` for a pattern that does not parse, and `Error::Value`
/// for one with fewer than two atoms (RDKit reads each match's first two).
pub fn fragment_mol(mol: &Mol) -> Result<Vec<(Option<Mol>, Mol)>> {
    fragment_mol_with_params(mol, &FragmentMolParams::default())
}

/// Mirrors `rdMMPA.FragmentMol(mol, maxCuts, maxCutBonds, pattern)`.
pub fn fragment_mol_with_params(
    mol: &Mol,
    params: &FragmentMolParams,
) -> Result<Vec<(Option<Mol>, Mol)>> {
    Ok(into_pairs(ffi::fragment_mol(
        mol.as_ffi(),
        1,
        params.max_cuts,
        params.max_cut_bonds,
        &params.pattern,
    )?))
}

/// Mirrors `rdMMPA.FragmentMol(mol, resultsAsMols=False)`.
pub fn fragment_mol_as_smiles(mol: &Mol) -> Result<Vec<(String, String)>> {
    as_smiles(fragment_mol(mol)?)
}

/// Mirrors `rdMMPA.FragmentMol(mol, maxCuts, maxCutBonds, pattern,
/// resultsAsMols=False)`.
pub fn fragment_mol_as_smiles_with_params(
    mol: &Mol,
    params: &FragmentMolParams,
) -> Result<Vec<(String, String)>> {
    as_smiles(fragment_mol_with_params(mol, params)?)
}

/// Every way of cutting `min_cuts` to `max_cuts` of the bonds the default
/// pattern matches. Mirrors `rdMMPA.FragmentMol(mol, minCuts, maxCuts,
/// maxCutBonds)`. `Error::Value` for `min_cuts` of 0 or above
/// `max_cuts`.
pub fn fragment_mol_min_cuts(
    mol: &Mol,
    min_cuts: u32,
    max_cuts: u32,
    max_cut_bonds: u32,
) -> Result<Vec<(Option<Mol>, Mol)>> {
    fragment_mol_min_cuts_with_pattern(mol, min_cuts, max_cuts, max_cut_bonds, DEFAULT_PATTERN)
}

/// Mirrors `rdMMPA.FragmentMol(mol, minCuts, maxCuts, maxCutBonds,
/// pattern)`.
pub fn fragment_mol_min_cuts_with_pattern(
    mol: &Mol,
    min_cuts: u32,
    max_cuts: u32,
    max_cut_bonds: u32,
    pattern: &str,
) -> Result<Vec<(Option<Mol>, Mol)>> {
    Ok(into_pairs(ffi::fragment_mol(
        mol.as_ffi(),
        min_cuts,
        max_cuts,
        max_cut_bonds,
        pattern,
    )?))
}

/// Mirrors `rdMMPA.FragmentMol(mol, minCuts, maxCuts, maxCutBonds,
/// resultsAsMols=False)`.
pub fn fragment_mol_min_cuts_as_smiles(
    mol: &Mol,
    min_cuts: u32,
    max_cuts: u32,
    max_cut_bonds: u32,
) -> Result<Vec<(String, String)>> {
    as_smiles(fragment_mol_min_cuts(
        mol,
        min_cuts,
        max_cuts,
        max_cut_bonds,
    )?)
}

/// Mirrors `rdMMPA.FragmentMol(mol, minCuts, maxCuts, maxCutBonds,
/// pattern, resultsAsMols=False)`.
pub fn fragment_mol_min_cuts_as_smiles_with_pattern(
    mol: &Mol,
    min_cuts: u32,
    max_cuts: u32,
    max_cut_bonds: u32,
    pattern: &str,
) -> Result<Vec<(String, String)>> {
    as_smiles(fragment_mol_min_cuts_with_pattern(
        mol,
        min_cuts,
        max_cuts,
        max_cut_bonds,
        pattern,
    )?)
}

/// Parameters for `fragment_mol_bonds_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentMolBondsParams {
    /// The fewest bonds cut at once.
    pub min_cuts: u32,
    /// The most bonds cut at once.
    pub max_cuts: u32,
}

impl Default for FragmentMolBondsParams {
    fn default() -> FragmentMolBondsParams {
        FragmentMolBondsParams {
            min_cuts: 1,
            max_cuts: 3,
        }
    }
}

/// Every way of cutting one to three of the given bonds. Mirrors
/// `rdMMPA.FragmentMol(mol, bondsToCut)`. `Error::Value` for no bonds or
/// a bond given twice, and `Error::Index` for a bond past the end.
pub fn fragment_mol_bonds(mol: &Mol, bonds_to_cut: &[u32]) -> Result<Vec<(Option<Mol>, Mol)>> {
    fragment_mol_bonds_with_params(mol, bonds_to_cut, &FragmentMolBondsParams::default())
}

/// Mirrors `rdMMPA.FragmentMol(mol, bondsToCut, minCuts, maxCuts)`.
pub fn fragment_mol_bonds_with_params(
    mol: &Mol,
    bonds_to_cut: &[u32],
    params: &FragmentMolBondsParams,
) -> Result<Vec<(Option<Mol>, Mol)>> {
    Ok(into_pairs(ffi::fragment_mol_bonds(
        mol.as_ffi(),
        bonds_to_cut,
        params.min_cuts,
        params.max_cuts,
    )?))
}

/// Mirrors `rdMMPA.FragmentMol(mol, bondsToCut, resultsAsMols=False)`.
pub fn fragment_mol_bonds_as_smiles(
    mol: &Mol,
    bonds_to_cut: &[u32],
) -> Result<Vec<(String, String)>> {
    as_smiles(fragment_mol_bonds(mol, bonds_to_cut)?)
}

/// Mirrors `rdMMPA.FragmentMol(mol, bondsToCut, minCuts, maxCuts,
/// resultsAsMols=False)`.
pub fn fragment_mol_bonds_as_smiles_with_params(
    mol: &Mol,
    bonds_to_cut: &[u32],
    params: &FragmentMolBondsParams,
) -> Result<Vec<(String, String)>> {
    as_smiles(fragment_mol_bonds_with_params(mol, bonds_to_cut, params)?)
}
