// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.MolSurf`, a Python module: Labute's approximate
//! surface area, the topological polar surface area, and the MOE-like VSA
//! descriptors, which bin each atom's share of the surface area by its
//! Crippen logP or MR contribution or its Gasteiger charge.
//!
//! The descriptors call `rdkit_cxx::chem::rd_mol_descriptors`, as Python's
//! do. Python's `force` arguments, which only steer its caches and RDKit's
//! (CONVENTIONS.md U8), are not ported: every result is computed afresh.
//! Each `SlogP_VSA*`, `SMR_VSA*`, and `PEOE_VSA*` function in Python takes
//! its bin as a defaulted second argument; here each function is its bin.

use crate::Result;
use crate::chem::rd_mol_descriptors as descriptors;
use crate::chem::rd_partial_charges::compute_gasteiger_charges;
use crate::chem::{BondType, Mol};
use crate::py::bisect_right;
use crate::py::py_sum;

/// How much each bond type shortens the sum of its atoms' radii in Python's
/// reference Labute calculation. Mirrors `MolSurf.bondScaleFacts`, which
/// RDKit 2024.09 keeps as a list of the aromatic, single, double, and triple
/// factors.
pub fn bond_scale_facts() -> Vec<(BondType, f64)> {
    vec![
        (BondType::Aromatic, 0.1),
        (BondType::Single, 0.0),
        (BondType::Double, 0.2),
        (BondType::Triple, 0.3),
    ]
}

const MR_BINS: [f64; 9] = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0];
const LOGP_BINS: [f64; 11] = [-0.4, -0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6];
const CHG_BINS: [f64; 13] = [
    -0.3, -0.25, -0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
];

/// The bin boundaries of the `SMR_VSA` descriptors. Mirrors `MolSurf.mrBins`.
pub fn mr_bins() -> &'static [f64] {
    &MR_BINS
}

/// The bin boundaries of the `SlogP_VSA` descriptors. Mirrors
/// `MolSurf.logpBins`.
pub fn logp_bins() -> &'static [f64] {
    &LOGP_BINS
}

/// The bin boundaries of the `PEOE_VSA` descriptors. Mirrors
/// `MolSurf.chgBins`.
pub fn chg_bins() -> &'static [f64] {
    &CHG_BINS
}

/// The implicit hydrogens' contribution to Labute's surface area (0 unless
/// `include_hs`), then each atom's. Python's private `_LabuteHelper`.
pub(crate) fn labute_helper(mol: &Mol, include_hs: bool) -> Result<Vec<f64>> {
    descriptors::calc_labute_asa_contribs(mol, include_hs)
}

/// Labute's approximate surface area, with hydrogens. Mirrors
/// `MolSurf.LabuteASA`, which calls `rdMolDescriptors.CalcLabuteASA`.
pub fn labute_asa(mol: &Mol) -> Result<f64> {
    descriptors::calc_labute_asa(mol)
}

/// Without implicit hydrogens unless `include_hs`. Mirrors
/// `MolSurf.LabuteASA(mol, includeHs)`.
pub fn labute_asa_with_include_hs(mol: &Mol, include_hs: bool) -> Result<f64> {
    descriptors::calc_labute_asa_with_include_hs(mol, include_hs)
}

/// Labute's approximate surface area as Python sums it, from the atoms'
/// contributions (P. Labute, J. Mol. Graph. Mod. 18, 464-477 (2000)).
/// Mirrors `MolSurf.pyLabuteASA`.
pub fn py_labute_asa(mol: &Mol) -> Result<f64> {
    py_labute_asa_with_include_hs(mol, true)
}

/// Without implicit hydrogens unless `include_hs`. Mirrors
/// `MolSurf.pyLabuteASA(mol, includeHs)`.
pub fn py_labute_asa_with_include_hs(mol: &Mol, include_hs: bool) -> Result<f64> {
    Ok(py_sum(labute_helper(mol, include_hs)?))
}

/// The topological polar surface area, from nitrogens and oxygens (P. Ertl,
/// B. Rohde, and P. Selzer, J. Med. Chem. 43, 3714-3717 (2000)). Mirrors
/// `MolSurf.TPSA`, which calls `rdMolDescriptors.CalcTPSA`.
pub fn tpsa(mol: &Mol) -> Result<f64> {
    descriptors::calc_tpsa(mol)
}

/// With sulfur and phosphorus contributions when `include_s_and_p`. Mirrors
/// `MolSurf.TPSA(mol, includeSandP=...)`.
pub fn tpsa_with_include_s_and_p(mol: &Mol, include_s_and_p: bool) -> Result<f64> {
    descriptors::calc_tpsa_with_include_s_and_p(mol, include_s_and_p)
}

/// Adds each atom's surface area (`vol_contribs`, after the hydrogens') to
/// the bin of its property, as Python's reference VSA functions do.
fn binned(props: &[f64], vol_contribs: &[f64], bins: &[f64]) -> Vec<f64> {
    let mut ans = vec![0.0; bins.len() + 1];
    for (prop, vol) in props.iter().zip(&vol_contribs[1..]) {
        ans[bisect_right(bins, *prop)] += vol;
    }
    ans
}

/// The `SMR_VSA` bins computed in Rust, as Python's reference
/// implementation computes them, with the default bins. Mirrors
/// `MolSurf.pySMR_VSA_`.
pub fn py_smr_vsa(mol: &Mol) -> Result<Vec<f64>> {
    py_smr_vsa_with_bins(mol, &MR_BINS)
}

/// With the given bin boundaries, which, as in Python, need not be sorted.
/// Mirrors `MolSurf.pySMR_VSA_(mol, bins)`.
pub fn py_smr_vsa_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    let props: Vec<f64> = descriptors::calc_crippen_contribs(mol)?
        .into_iter()
        .map(|(_, mr)| mr)
        .collect();
    Ok(binned(&props, &labute_helper(mol, true)?, bins))
}

/// The `SlogP_VSA` bins computed as Python's reference implementation
/// computes them, with the default bins. Mirrors `MolSurf.pySlogP_VSA_`.
pub fn py_slogp_vsa(mol: &Mol) -> Result<Vec<f64>> {
    py_slogp_vsa_with_bins(mol, &LOGP_BINS)
}

/// With the given bin boundaries, which need not be sorted. Mirrors
/// `MolSurf.pySlogP_VSA_(mol, bins)`.
pub fn py_slogp_vsa_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    let props: Vec<f64> = descriptors::calc_crippen_contribs(mol)?
        .into_iter()
        .map(|(logp, _)| logp)
        .collect();
    Ok(binned(&props, &labute_helper(mol, true)?, bins))
}

/// The `PEOE_VSA` bins computed as Python's reference implementation
/// computes them, with the default bins. Python computes the Gasteiger
/// charges on the molecule it is given; this computes them on a copy.
/// Mirrors `MolSurf.pyPEOE_VSA_`.
pub fn py_peoe_vsa(mol: &Mol) -> Result<Vec<f64>> {
    py_peoe_vsa_with_bins(mol, &CHG_BINS)
}

/// With the given bin boundaries, which need not be sorted. Mirrors
/// `MolSurf.pyPEOE_VSA_(mol, bins)`.
pub fn py_peoe_vsa_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    let props = gasteiger_charges(mol)?;
    Ok(binned(&props, &labute_helper(mol, true)?, bins))
}

/// Each atom's Gasteiger charge, computed on a copy of the molecule.
pub(crate) fn gasteiger_charges(mol: &Mol) -> Result<Vec<f64>> {
    let mut copy = mol.clone();
    compute_gasteiger_charges(&mut copy)?;
    copy.get_atoms()
        .map(|atom| atom.get_double_prop("_GasteigerCharge"))
        .collect()
}

/// Defines one function per bin of a VSA descriptor vector.
macro_rules! vsa_bins {
    ($vector:path, $what:literal, $bins:literal: $($name:ident = $index:literal, $python:literal;)*) => {
        $(
            #[doc = concat!(
                "The surface area of the atoms whose ", $what, " falls in bin ", stringify!($index),
                " (from 0) of [`", $bins, "`]. Mirrors `MolSurf.", $python, "`."
            )]
            pub fn $name(mol: &Mol) -> Result<f64> {
                Ok($vector(mol)?[$index])
            }
        )*
    };
}

vsa_bins! {
    descriptors::smr_vsa, "Crippen MR contribution", "mr_bins":
    smr_vsa1 = 0, "SMR_VSA1";
    smr_vsa2 = 1, "SMR_VSA2";
    smr_vsa3 = 2, "SMR_VSA3";
    smr_vsa4 = 3, "SMR_VSA4";
    smr_vsa5 = 4, "SMR_VSA5";
    smr_vsa6 = 5, "SMR_VSA6";
    smr_vsa7 = 6, "SMR_VSA7";
    smr_vsa8 = 7, "SMR_VSA8";
    smr_vsa9 = 8, "SMR_VSA9";
    smr_vsa10 = 9, "SMR_VSA10";
}

vsa_bins! {
    descriptors::slogp_vsa, "Crippen logP contribution", "logp_bins":
    slogp_vsa1 = 0, "SlogP_VSA1";
    slogp_vsa2 = 1, "SlogP_VSA2";
    slogp_vsa3 = 2, "SlogP_VSA3";
    slogp_vsa4 = 3, "SlogP_VSA4";
    slogp_vsa5 = 4, "SlogP_VSA5";
    slogp_vsa6 = 5, "SlogP_VSA6";
    slogp_vsa7 = 6, "SlogP_VSA7";
    slogp_vsa8 = 7, "SlogP_VSA8";
    slogp_vsa9 = 8, "SlogP_VSA9";
    slogp_vsa10 = 9, "SlogP_VSA10";
    slogp_vsa11 = 10, "SlogP_VSA11";
    slogp_vsa12 = 11, "SlogP_VSA12";
}

vsa_bins! {
    descriptors::peoe_vsa, "Gasteiger charge", "chg_bins":
    peoe_vsa1 = 0, "PEOE_VSA1";
    peoe_vsa2 = 1, "PEOE_VSA2";
    peoe_vsa3 = 2, "PEOE_VSA3";
    peoe_vsa4 = 3, "PEOE_VSA4";
    peoe_vsa5 = 4, "PEOE_VSA5";
    peoe_vsa6 = 5, "PEOE_VSA6";
    peoe_vsa7 = 6, "PEOE_VSA7";
    peoe_vsa8 = 7, "PEOE_VSA8";
    peoe_vsa9 = 8, "PEOE_VSA9";
    peoe_vsa10 = 9, "PEOE_VSA10";
    peoe_vsa11 = 10, "PEOE_VSA11";
    peoe_vsa12 = 11, "PEOE_VSA12";
    peoe_vsa13 = 12, "PEOE_VSA13";
    peoe_vsa14 = 13, "PEOE_VSA14";
}
