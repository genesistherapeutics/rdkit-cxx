// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdReducedGraphs`: extended reduced graphs and ErG
//! fingerprints.
//!
//! Python's `atomTypes` arguments raise unless they are `None`, so they
//! are not bound.

use crate::Result;
use crate::chem::Mol;
use crate::ffi::reduced_graphs::ffi;

/// The molecule's extended reduced graph: pharmacophore feature types in
/// each atom's `_ErGAtomTypes` property, and a dummy atom for each ring
/// smaller than eight. Mirrors
/// `rdReducedGraphs.GenerateMolExtendedReducedGraph`. `Error::Runtime`
/// for a molecule without ring information.
pub fn generate_mol_extended_reduced_graph(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::generate_mol_extended_reduced_graph(
        mol.as_ffi(),
    )?))
}

/// Parameters for the ErG fingerprint functions' `_with_params` forms,
/// with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetErGFingerprintParams {
    /// The amount added to the neighboring distance bins.
    pub fuzz_increment: f64,
    /// The shortest path, in bonds, counted.
    pub min_path: i32,
    /// The longest path counted; `Error::Value` unless
    /// `0 <= min_path < max_path` (RDKit counts them as unsigned).
    pub max_path: i32,
}

impl Default for GetErGFingerprintParams {
    fn default() -> GetErGFingerprintParams {
        GetErGFingerprintParams {
            fuzz_increment: 0.3,
            min_path: 1,
            max_path: 15,
        }
    }
}

/// The ErG fingerprint of a molecule that is already a reduced graph (from
/// `generate_mol_extended_reduced_graph`). Mirrors
/// `rdReducedGraphs.GenerateErGFingerprintForReducedGraph`. `Error::Key`
/// for an atom without `_ErGAtomTypes`; `Error::Value` for a type that is
/// not an ErG type (RDKit computes the bin from it unchecked), or above
/// 46340 atoms (U3).
pub fn generate_erg_fingerprint_for_reduced_graph(mol: &Mol) -> Result<Vec<f64>> {
    generate_erg_fingerprint_for_reduced_graph_with_params(mol, &GetErGFingerprintParams::default())
}

/// Mirrors `rdReducedGraphs.GenerateErGFingerprintForReducedGraph(mol,
/// atomTypes, fuzzIncrement, minPath, maxPath)`.
pub fn generate_erg_fingerprint_for_reduced_graph_with_params(
    mol: &Mol,
    params: &GetErGFingerprintParams,
) -> Result<Vec<f64>> {
    Ok(ffi::generate_erg_fingerprint_for_reduced_graph(
        mol.as_ffi(),
        params.fuzz_increment,
        params.min_path,
        params.max_path,
    )?)
}

/// The molecule's ErG fingerprint. Mirrors
/// `rdReducedGraphs.GetErGFingerprint`.
pub fn get_erg_fingerprint(mol: &Mol) -> Result<Vec<f64>> {
    get_erg_fingerprint_with_params(mol, &GetErGFingerprintParams::default())
}

/// Mirrors `rdReducedGraphs.GetErGFingerprint(mol, atomTypes,
/// fuzzIncrement, minPath, maxPath)`.
pub fn get_erg_fingerprint_with_params(
    mol: &Mol,
    params: &GetErGFingerprintParams,
) -> Result<Vec<f64>> {
    Ok(ffi::get_erg_fingerprint(
        mol.as_ffi(),
        params.fuzz_increment,
        params.min_path,
        params.max_path,
    )?)
}
