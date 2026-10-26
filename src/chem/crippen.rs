// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Crippen`, a Python module: Wildman and Crippen's
//! atom-based logP and molar refractivity (S. A. Wildman and G. M. Crippen,
//! JCICS 39, 868-873 (1999)).
//!
//! Python's functions call `rdMolDescriptors.CalcCrippenDescriptors`, as
//! these do. Its private pattern reader and Python reference
//! implementations, and its `force` arguments, which only steer RDKit's
//! descriptor caches (CONVENTIONS.md U8), are not ported.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rd_mol_descriptors::calc_crippen_descriptors_with_include_hs;

/// Wildman and Crippen's logP, with hydrogens. Mirrors `Crippen.MolLogP`.
pub fn mol_log_p(mol: &Mol) -> Result<f64> {
    mol_log_p_with_include_hs(mol, true)
}

/// Without implicit hydrogens unless `include_hs`. Mirrors
/// `Crippen.MolLogP(mol, includeHs)`.
pub fn mol_log_p_with_include_hs(mol: &Mol, include_hs: bool) -> Result<f64> {
    Ok(calc_crippen_descriptors_with_include_hs(mol, include_hs)?.0)
}

/// Wildman and Crippen's molar refractivity, with hydrogens. Mirrors
/// `Crippen.MolMR`.
pub fn mol_mr(mol: &Mol) -> Result<f64> {
    mol_mr_with_include_hs(mol, true)
}

/// Without implicit hydrogens unless `include_hs`. Mirrors
/// `Crippen.MolMR(mol, includeHs)`.
pub fn mol_mr_with_include_hs(mol: &Mol, include_hs: bool) -> Result<f64> {
    Ok(calc_crippen_descriptors_with_include_hs(mol, include_hs)?.1)
}
