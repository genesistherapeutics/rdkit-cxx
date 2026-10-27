// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdCIPLabeler`: CIP labels from a full
//! implementation of the CIP rules.

use crate::Result;
use crate::chem::Mol;
use crate::ffi::cip_labeler::ffi;

/// Parameters for `assign_cip_labels_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AssignCIPLabelsParams {
    /// The atoms to label. With neither list, every atom and bond is
    /// labeled. Python treats an empty list as `None`.
    pub atoms_to_label: Option<Vec<u32>>,
    /// The bonds to label. Python treats an empty list as `None`.
    pub bonds_to_label: Option<Vec<u32>>,
    /// The labeler's iteration limit, 0 for none. Exceeding it is
    /// `Error::Runtime` (Python's `RuntimeError`).
    pub max_recursive_iterations: u32,
}

/// Labels the molecule's stereocenters and stereo bonds, in their
/// `_CIPCode` properties. Mirrors `rdCIPLabeler.AssignCIPLabels`. An
/// atropisomer bond with an end of no other bonds whose other end has two
/// other bonds is left out, with any label it had: RDKit 2024.09 (and
/// Python) reads past its empty carriers for one.
pub fn assign_cip_labels(mol: &mut Mol) -> Result<()> {
    assign_cip_labels_with_params(mol, &AssignCIPLabelsParams::default())
}

/// Mirrors `rdCIPLabeler.AssignCIPLabels(mol, atomsToLabel, bondsToLabel,
/// maxRecursiveIterations)`. `Error::Index` for an atom or bond past the
/// end.
pub fn assign_cip_labels_with_params(mol: &mut Mol, params: &AssignCIPLabelsParams) -> Result<()> {
    Ok(ffi::assign_cip_labels(
        mol.pin_mut(),
        params.atoms_to_label.is_some(),
        params.atoms_to_label.as_deref().unwrap_or_default(),
        params.bonds_to_label.is_some(),
        params.bonds_to_label.as_deref().unwrap_or_default(),
        params.max_recursive_iterations,
    )?)
}
