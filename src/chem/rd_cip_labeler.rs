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
/// atropisomer bond with an end of no other bonds is left unlabeled: RDKit
/// 2024.03 (and Python) reads past an empty vector for one. RDKit 2024.03
/// pairs stereo descriptors in a 32-bit mask, by a shift that is undefined
/// past 32 (`PairList`), so labeling more than 32 stereo atoms and bonds
/// (those with a tetrahedral tag or a double-bond or atropisomer stereo,
/// among the ones asked for) is `Error::Value`. Within that limit, a fused
/// ring system that repeats one stereo atom in a branch of RDKit's digraph
/// many times could still reach the shift; this is not checked, since a
/// bound that covered it would refuse common polycycles such as steroids.
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
