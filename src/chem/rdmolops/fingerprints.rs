// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdmolops's path-based fingerprints.

use std::collections::BTreeMap;

use crate::Result;
use crate::chem::{Mol, MolBundle};
use crate::data_structs::{ExplicitBitVect, ULongSparseIntVect};
use crate::error::expect_ok;
use crate::ffi::fingerprints::ffi;

/// What Python's output arguments receive (CONVENTIONS.md A9).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FingerprintOutputs {
    /// The bits each atom sets. Python's `atomBits`.
    pub atom_bits: Vec<Vec<u64>>,
    /// The paths, as bond indices, that set each bit. Python's `bitInfo`.
    pub bit_info: BTreeMap<u64, Vec<Vec<u32>>>,
    /// The paths each atom is in, added to the counts the parameters gave.
    /// Python's `atomCounts`.
    pub atom_counts: Vec<u32>,
}

impl FingerprintOutputs {
    fn ffi() -> ffi::FingerprintOutputs {
        ffi::FingerprintOutputs {
            atom_bits: Vec::new(),
            bit_info_bits: Vec::new(),
            bit_info_num_paths: Vec::new(),
            bit_info_paths: Vec::new(),
            atom_counts: Vec::new(),
        }
    }

    fn from_ffi(outputs: ffi::FingerprintOutputs) -> FingerprintOutputs {
        let mut paths = outputs.bit_info_paths.into_iter();
        let bit_info = outputs
            .bit_info_bits
            .into_iter()
            .zip(outputs.bit_info_num_paths)
            .map(|(bit, n)| {
                let bit_paths = paths.by_ref().take(n as usize).map(|p| p.indices).collect();
                (bit, bit_paths)
            })
            .collect();
        FingerprintOutputs {
            atom_bits: outputs.atom_bits.into_iter().map(|b| b.values).collect(),
            bit_info,
            atom_counts: outputs.atom_counts,
        }
    }
}

/// Parameters for `rdk_fingerprint_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct RDKFingerprintParams {
    /// The fewest bonds in a path.
    pub min_path: u32,
    /// The most bonds in a path.
    pub max_path: u32,
    /// The number of bits.
    pub fp_size: u32,
    /// The bits each path sets.
    pub n_bits_per_hash: u32,
    /// Includes hydrogens in the paths.
    pub use_hs: bool,
    /// Folds the fingerprint until its density reaches this, if above 0.
    pub tgt_density: f64,
    /// The fewest bits folding leaves.
    pub min_size: u32,
    /// Uses branched subgraphs, not only linear paths.
    pub branched_paths: bool,
    /// Uses bond orders in the hashes.
    pub use_bond_order: bool,
    /// An invariant for each atom, in place of RDKit's.
    pub atom_invariants: Option<Vec<u32>>,
    /// The atoms paths must start from; `None` for all.
    pub from_atoms: Option<Vec<u32>>,
}

impl Default for RDKFingerprintParams {
    fn default() -> RDKFingerprintParams {
        RDKFingerprintParams {
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            n_bits_per_hash: 2,
            use_hs: true,
            tgt_density: 0.0,
            min_size: 128,
            branched_paths: true,
            use_bond_order: true,
            atom_invariants: None,
            from_atoms: None,
        }
    }
}

impl RDKFingerprintParams {
    fn to_ffi(&self) -> ffi::RDKFingerprintData {
        ffi::RDKFingerprintData {
            min_path: self.min_path,
            max_path: self.max_path,
            fp_size: self.fp_size,
            n_bits_per_hash: self.n_bits_per_hash,
            use_hs: self.use_hs,
            tgt_density: self.tgt_density,
            min_size: self.min_size,
            branched_paths: self.branched_paths,
            use_bond_order: self.use_bond_order,
            has_atom_invariants: self.atom_invariants.is_some(),
            atom_invariants: self.atom_invariants.clone().unwrap_or_default(),
            has_from_atoms: self.from_atoms.is_some(),
            from_atoms: self.from_atoms.clone().unwrap_or_default(),
        }
    }
}

/// RDKit's path fingerprint, a Daylight-like hash of the molecule's paths
/// and subgraphs. Mirrors `rdmolops.RDKFingerprint`.
pub fn rdk_fingerprint(mol: &Mol) -> Result<ExplicitBitVect> {
    rdk_fingerprint_with_params(mol, &RDKFingerprintParams::default())
}

/// Mirrors `rdmolops.RDKFingerprint(mol, minPath, ..., fromAtoms)`.
pub fn rdk_fingerprint_with_params(
    mol: &Mol,
    params: &RDKFingerprintParams,
) -> Result<ExplicitBitVect> {
    let mut outputs = FingerprintOutputs::ffi();
    Ok(ExplicitBitVect::from_ffi(ffi::rdk_fingerprint(
        mol.as_ffi(),
        &params.to_ffi(),
        false,
        &mut outputs,
    )?))
}

/// The fingerprint and Python's `atomBits` and `bitInfo` (A9). Mirrors
/// `rdmolops.RDKFingerprint(mol, ..., atomBits, bitInfo)`.
pub fn rdk_fingerprint_with_outputs(
    mol: &Mol,
    params: &RDKFingerprintParams,
) -> Result<(ExplicitBitVect, FingerprintOutputs)> {
    let mut outputs = FingerprintOutputs::ffi();
    let fp = ffi::rdk_fingerprint(mol.as_ffi(), &params.to_ffi(), true, &mut outputs)?;
    Ok((
        ExplicitBitVect::from_ffi(fp),
        FingerprintOutputs::from_ffi(outputs),
    ))
}

/// Parameters for `unfolded_rdk_fingerprint_count_based_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UnfoldedRDKFingerprintCountBasedParams {
    /// The fewest bonds in a path.
    pub min_path: u32,
    /// The most bonds in a path.
    pub max_path: u32,
    /// Includes hydrogens in the paths.
    pub use_hs: bool,
    /// Uses branched subgraphs, not only linear paths.
    pub branched_paths: bool,
    /// Uses bond orders in the hashes.
    pub use_bond_order: bool,
    /// An invariant for each atom, in place of RDKit's.
    pub atom_invariants: Option<Vec<u32>>,
    /// The atoms paths must start from; `None` for all.
    pub from_atoms: Option<Vec<u32>>,
}

impl Default for UnfoldedRDKFingerprintCountBasedParams {
    fn default() -> UnfoldedRDKFingerprintCountBasedParams {
        UnfoldedRDKFingerprintCountBasedParams {
            min_path: 1,
            max_path: 7,
            use_hs: true,
            branched_paths: true,
            use_bond_order: true,
            atom_invariants: None,
            from_atoms: None,
        }
    }
}

impl UnfoldedRDKFingerprintCountBasedParams {
    fn to_ffi(&self) -> ffi::RDKFingerprintData {
        RDKFingerprintParams {
            min_path: self.min_path,
            max_path: self.max_path,
            use_hs: self.use_hs,
            branched_paths: self.branched_paths,
            use_bond_order: self.use_bond_order,
            atom_invariants: self.atom_invariants.clone(),
            from_atoms: self.from_atoms.clone(),
            ..Default::default()
        }
        .to_ffi()
    }
}

/// The RDKit fingerprint unfolded, with a count for each path hash.
/// Mirrors `rdmolops.UnfoldedRDKFingerprintCountBased`.
pub fn unfolded_rdk_fingerprint_count_based(mol: &Mol) -> Result<ULongSparseIntVect> {
    unfolded_rdk_fingerprint_count_based_with_params(mol, &Default::default())
}

/// Mirrors `rdmolops.UnfoldedRDKFingerprintCountBased(mol, minPath, ...,
/// fromAtoms)`.
pub fn unfolded_rdk_fingerprint_count_based_with_params(
    mol: &Mol,
    params: &UnfoldedRDKFingerprintCountBasedParams,
) -> Result<ULongSparseIntVect> {
    let mut outputs = FingerprintOutputs::ffi();
    Ok(ULongSparseIntVect::from_ffi(
        ffi::unfolded_rdk_fingerprint_count_based(
            mol.as_ffi(),
            &params.to_ffi(),
            false,
            &mut outputs,
        )?,
    ))
}

/// The fingerprint and Python's `atomBits` and `bitInfo` (A9). Mirrors
/// `rdmolops.UnfoldedRDKFingerprintCountBased(mol, ..., atomBits, bitInfo)`.
pub fn unfolded_rdk_fingerprint_count_based_with_outputs(
    mol: &Mol,
    params: &UnfoldedRDKFingerprintCountBasedParams,
) -> Result<(ULongSparseIntVect, FingerprintOutputs)> {
    let mut outputs = FingerprintOutputs::ffi();
    let fp = ffi::unfolded_rdk_fingerprint_count_based(
        mol.as_ffi(),
        &params.to_ffi(),
        true,
        &mut outputs,
    )?;
    Ok((
        ULongSparseIntVect::from_ffi(fp),
        FingerprintOutputs::from_ffi(outputs),
    ))
}

/// The layers `layered_fingerprint` uses for substructure screening.
/// Mirrors `rdmolops.LayeredFingerprint_substructLayers`.
pub fn layered_fingerprint_substruct_layers() -> u32 {
    expect_ok(ffi::layered_fingerprint_substruct_layers())
}

/// Parameters for `layered_fingerprint_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredFingerprintParams {
    /// The layers to include, as bits.
    pub layer_flags: u32,
    /// The fewest bonds in a path.
    pub min_path: u32,
    /// The most bonds in a path.
    pub max_path: u32,
    /// The number of bits.
    pub fp_size: u32,
    /// Counts to add each atom's paths to, one per atom or more; the sums
    /// come back in `FingerprintOutputs::atom_counts`.
    pub atom_counts: Option<Vec<u32>>,
    /// Sets only the bits on in this vector, which needs `fp_size` bits.
    pub set_only_bits: Option<ExplicitBitVect>,
    /// Uses branched subgraphs, not only linear paths.
    pub branched_paths: bool,
    /// The atoms paths must start from; `None` for all.
    pub from_atoms: Option<Vec<u32>>,
}

impl Default for LayeredFingerprintParams {
    fn default() -> LayeredFingerprintParams {
        LayeredFingerprintParams {
            layer_flags: 0xFFFFFFFF,
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            atom_counts: None,
            set_only_bits: None,
            branched_paths: true,
            from_atoms: None,
        }
    }
}

impl LayeredFingerprintParams {
    fn to_ffi(&self) -> ffi::LayeredFingerprintData {
        ffi::LayeredFingerprintData {
            layer_flags: self.layer_flags,
            min_path: self.min_path,
            max_path: self.max_path,
            fp_size: self.fp_size,
            has_atom_counts: self.atom_counts.is_some(),
            atom_counts: self.atom_counts.clone().unwrap_or_default(),
            has_set_only_bits: self.set_only_bits.is_some(),
            branched_paths: self.branched_paths,
            has_from_atoms: self.from_atoms.is_some(),
            from_atoms: self.from_atoms.clone().unwrap_or_default(),
        }
    }
}

fn only_bits(bits: &Option<ExplicitBitVect>) -> ExplicitBitVect {
    // A placeholder when there is none, which the wrapper ignores.
    bits.clone().unwrap_or_default()
}

/// A path fingerprint in layers of atom and bond features, for
/// substructure screening. Mirrors `rdmolops.LayeredFingerprint`.
pub fn layered_fingerprint(mol: &Mol) -> Result<ExplicitBitVect> {
    layered_fingerprint_with_params(mol, &LayeredFingerprintParams::default())
}

/// Mirrors `rdmolops.LayeredFingerprint(mol, layerFlags, ..., fromAtoms)`.
pub fn layered_fingerprint_with_params(
    mol: &Mol,
    params: &LayeredFingerprintParams,
) -> Result<ExplicitBitVect> {
    Ok(layered_fingerprint_with_outputs(mol, params)?.0)
}

/// The fingerprint and the updated `atomCounts` (A9). Mirrors
/// `rdmolops.LayeredFingerprint(mol, ..., atomCounts)`.
pub fn layered_fingerprint_with_outputs(
    mol: &Mol,
    params: &LayeredFingerprintParams,
) -> Result<(ExplicitBitVect, FingerprintOutputs)> {
    let mut outputs = FingerprintOutputs::ffi();
    let fp = ffi::layered_fingerprint(
        mol.as_ffi(),
        &params.to_ffi(),
        only_bits(&params.set_only_bits).raw(),
        &mut outputs,
    )?;
    Ok((
        ExplicitBitVect::from_ffi(fp),
        FingerprintOutputs::from_ffi(outputs),
    ))
}

/// Parameters for `pattern_fingerprint_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct PatternFingerprintParams {
    /// The number of bits.
    pub fp_size: u32,
    /// Counts to add each atom's patterns to, one per atom or more; the sums
    /// come back in `FingerprintOutputs::atom_counts`.
    pub atom_counts: Option<Vec<u32>>,
    /// Sets only the bits on in this vector, which needs `fp_size` bits.
    pub set_only_bits: Option<ExplicitBitVect>,
    /// Makes the fingerprint the same for tautomers.
    pub tautomer_fingerprints: bool,
}

impl Default for PatternFingerprintParams {
    fn default() -> PatternFingerprintParams {
        PatternFingerprintParams {
            fp_size: 2048,
            atom_counts: None,
            set_only_bits: None,
            tautomer_fingerprints: false,
        }
    }
}

impl PatternFingerprintParams {
    fn to_ffi(&self) -> ffi::PatternFingerprintData {
        ffi::PatternFingerprintData {
            fp_size: self.fp_size,
            has_atom_counts: self.atom_counts.is_some(),
            atom_counts: self.atom_counts.clone().unwrap_or_default(),
            has_set_only_bits: self.set_only_bits.is_some(),
            tautomer_fingerprints: self.tautomer_fingerprints,
        }
    }
}

/// A fingerprint of SMARTS-like patterns, for substructure screening: a
/// molecule's bits include every substructure's. Mirrors
/// `rdmolops.PatternFingerprint`.
pub fn pattern_fingerprint(mol: &Mol) -> Result<ExplicitBitVect> {
    pattern_fingerprint_with_params(mol, &PatternFingerprintParams::default())
}

/// Mirrors `rdmolops.PatternFingerprint(mol, fpSize, atomCounts,
/// setOnlyBits, tautomerFingerprints)`.
pub fn pattern_fingerprint_with_params(
    mol: &Mol,
    params: &PatternFingerprintParams,
) -> Result<ExplicitBitVect> {
    Ok(pattern_fingerprint_with_outputs(mol, params)?.0)
}

/// The fingerprint and the updated `atomCounts` (A9). Mirrors
/// `rdmolops.PatternFingerprint(mol, ..., atomCounts)`.
pub fn pattern_fingerprint_with_outputs(
    mol: &Mol,
    params: &PatternFingerprintParams,
) -> Result<(ExplicitBitVect, FingerprintOutputs)> {
    let mut outputs = FingerprintOutputs::ffi();
    let fp = ffi::pattern_fingerprint(
        mol.as_ffi(),
        &params.to_ffi(),
        only_bits(&params.set_only_bits).raw(),
        &mut outputs,
    )?;
    Ok((
        ExplicitBitVect::from_ffi(fp),
        FingerprintOutputs::from_ffi(outputs),
    ))
}

/// The bits every molecule in the bundle sets in its pattern fingerprint.
/// Mirrors `rdmolops.PatternFingerprint(MolBundle)` (A5); its parameters
/// ignore `atom_counts`, which Python's overload lacks.
pub fn pattern_fingerprint_mol_bundle(bundle: &MolBundle) -> Result<ExplicitBitVect> {
    pattern_fingerprint_mol_bundle_with_params(bundle, &PatternFingerprintParams::default())
}

/// Mirrors `rdmolops.PatternFingerprint(bundle, fpSize, setOnlyBits,
/// tautomerFingerprints)`.
pub fn pattern_fingerprint_mol_bundle_with_params(
    bundle: &MolBundle,
    params: &PatternFingerprintParams,
) -> Result<ExplicitBitVect> {
    let set_only_bits = only_bits(&params.set_only_bits);
    let fp = bundle.with_ffi(|b| {
        ffi::pattern_fingerprint_mol_bundle(b, &params.to_ffi(), set_only_bits.raw())
    })?;
    Ok(ExplicitBitVect::from_ffi(fp))
}
