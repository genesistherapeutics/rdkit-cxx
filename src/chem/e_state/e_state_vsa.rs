// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.EState.EState_VSA`, a Python module: hybrid
//! EState-VSA descriptors, which bin EState indices by the atoms' Labute
//! surface areas (`VSA_EState`) or the surface areas by EState index
//! (`EState_VSA`), as the MOE VSA descriptors of `MolSurf` do.
//!
//! The default bins were chosen on the PP3K solubility data set to give
//! about as many atoms per bin. Python's `force` arguments only steer its
//! caches, so they are not ported, and each `VSA_EState*` and
//! `EState_VSA*` function is its bin.

use crate::Result;
use crate::chem::Mol;
use crate::chem::e_state::e_state::e_state_indices;
use crate::chem::mol_surf::labute_helper;
use crate::chem::py_util::bisect_right;

const VSA_BINS: [f64; 9] = [4.78, 5.00, 5.410, 5.740, 6.00, 6.07, 6.45, 7.00, 11.0];
const ESTATE_BINS: [f64; 10] = [
    -0.390, 0.290, 0.717, 1.165, 1.540, 1.807, 2.05, 4.69, 9.17, 15.0,
];

/// The surface-area bin boundaries of the `VSA_EState` descriptors. Mirrors
/// `EState_VSA.vsaBins`.
pub fn vsa_bins() -> &'static [f64] {
    &VSA_BINS
}

/// The EState bin boundaries of the `EState_VSA` descriptors. Mirrors
/// `EState_VSA.estateBins`.
pub fn estate_bins() -> &'static [f64] {
    &ESTATE_BINS
}

/// The sums of the atoms' EState indices, binned by their Labute surface
/// areas, with the default bins. Mirrors `EState_VSA.VSA_EState_`.
pub fn vsa_e_state(mol: &Mol) -> Result<Vec<f64>> {
    vsa_e_state_with_bins(mol, &VSA_BINS)
}

/// With the given bin boundaries, which, as in Python, need not be sorted.
/// Mirrors `EState_VSA.VSA_EState_(mol, bins)`.
pub fn vsa_e_state_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    let prop_contribs = e_state_indices(mol)?;
    let vol_contribs = labute_helper(mol, true)?;
    let mut ans = vec![0.0; bins.len() + 1];
    for (prop, vol) in prop_contribs.iter().zip(&vol_contribs[1..]) {
        ans[bisect_right(bins, *vol)] += prop;
    }
    Ok(ans)
}

/// The sums of the atoms' Labute surface areas, binned by their EState
/// indices, with the default bins. Mirrors `EState_VSA.EState_VSA_`.
pub fn e_state_vsa(mol: &Mol) -> Result<Vec<f64>> {
    e_state_vsa_with_bins(mol, &ESTATE_BINS)
}

/// With the given bin boundaries, which need not be sorted. Mirrors
/// `EState_VSA.EState_VSA_(mol, bins)`.
pub fn e_state_vsa_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    let prop_contribs = e_state_indices(mol)?;
    let vol_contribs = labute_helper(mol, true)?;
    let mut ans = vec![0.0; bins.len() + 1];
    for (prop, vol) in prop_contribs.iter().zip(&vol_contribs[1..]) {
        ans[bisect_right(bins, *prop)] += vol;
    }
    Ok(ans)
}

/// Defines one function per bin of a descriptor vector.
macro_rules! bins {
    ($vector:ident, $what:literal, $bins:literal: $($name:ident = $index:literal, $python:literal;)*) => {
        $(
            #[doc = concat!(
                "The ", $what, " bin ", stringify!($index), " (from 0) of [`", $bins,
                "`]. Mirrors `EState_VSA.", $python, "`."
            )]
            pub fn $name(mol: &Mol) -> Result<f64> {
                Ok($vector(mol)?[$index])
            }
        )*
    };
}

bins! {
    vsa_e_state, "sum of the EState indices of the atoms whose surface area falls in", "vsa_bins":
    vsa_e_state1 = 0, "VSA_EState1";
    vsa_e_state2 = 1, "VSA_EState2";
    vsa_e_state3 = 2, "VSA_EState3";
    vsa_e_state4 = 3, "VSA_EState4";
    vsa_e_state5 = 4, "VSA_EState5";
    vsa_e_state6 = 5, "VSA_EState6";
    vsa_e_state7 = 6, "VSA_EState7";
    vsa_e_state8 = 7, "VSA_EState8";
    vsa_e_state9 = 8, "VSA_EState9";
    vsa_e_state10 = 9, "VSA_EState10";
}

bins! {
    e_state_vsa, "surface area of the atoms whose EState index falls in", "estate_bins":
    e_state_vsa1 = 0, "EState_VSA1";
    e_state_vsa2 = 1, "EState_VSA2";
    e_state_vsa3 = 2, "EState_VSA3";
    e_state_vsa4 = 3, "EState_VSA4";
    e_state_vsa5 = 4, "EState_VSA5";
    e_state_vsa6 = 5, "EState_VSA6";
    e_state_vsa7 = 6, "EState_VSA7";
    e_state_vsa8 = 7, "EState_VSA8";
    e_state_vsa9 = 8, "EState_VSA9";
    e_state_vsa10 = 9, "EState_VSA10";
    e_state_vsa11 = 10, "EState_VSA11";
}
