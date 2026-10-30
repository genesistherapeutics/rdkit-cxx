// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdPartialCharges`: Gasteiger partial charges.

use crate::Result;
use crate::chem::Mol;
use crate::ffi::partial_charges::ffi;

/// Parameters for `compute_gasteiger_charges_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ComputeGasteigerChargesParams {
    /// The number of iterations; `Error::Value` when negative (RDKit counts
    /// them as unsigned).
    pub n_iter: i32,
    /// Returns `Error::Value` for an atom with no parameters, rather than
    /// giving it a charge of NaN.
    pub throw_on_param_failure: bool,
}

impl Default for ComputeGasteigerChargesParams {
    fn default() -> ComputeGasteigerChargesParams {
        ComputeGasteigerChargesParams {
            n_iter: 12,
            throw_on_param_failure: false,
        }
    }
}

/// Computes Gasteiger charges, storing each atom's in its
/// `_GasteigerCharge` property and its implicit hydrogens' total in
/// `_GasteigerHCharge`. Mirrors `rdPartialCharges.ComputeGasteigerCharges`.
/// `Error::Runtime` when the molecule's property cache is not up to date.
pub fn compute_gasteiger_charges(mol: &mut Mol) -> Result<()> {
    compute_gasteiger_charges_with_params(mol, &ComputeGasteigerChargesParams::default())
}

/// Mirrors `rdPartialCharges.ComputeGasteigerCharges(mol, nIter,
/// throwOnParamFailure)`.
pub fn compute_gasteiger_charges_with_params(
    mol: &mut Mol,
    params: &ComputeGasteigerChargesParams,
) -> Result<()> {
    Ok(ffi::compute_gasteiger_charges(
        mol.pin_mut(),
        params.n_iter,
        params.throw_on_param_failure,
    )?)
}
