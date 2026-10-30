// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Reaction sanitization.

use crate::Result;
use crate::chem::rd_chem_reactions::ChemicalReaction;
use crate::chem::rdmolops::{AdjustQueryParameters, SanitizeFlags as MolSanitizeFlags};
use crate::error::expect_ok;
use crate::ffi::chem_reactions::ffi;

/// Reaction sanitization operations, as bits. Mirrors
/// `rdChemReactions.SanitizeFlags` (CONVENTIONS.md F1); the constants keep
/// Python's spelling.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct SanitizeFlags(u32);

macro_rules! sanitize_flags {
    ($($name:ident),* $(,)?) => {
        impl SanitizeFlags {
            $(
                #[doc = concat!("Mirrors `SanitizeFlags.", stringify!($name), "`.")]
                pub const $name: SanitizeFlags = SanitizeFlags(ffi::SanitizeRxnFlags::$name.repr);
            )*

            const NAMES: &[(&str, SanitizeFlags)] = &[$((stringify!($name), SanitizeFlags::$name)),*];
        }
    };
}

sanitize_flags!(
    SANITIZE_NONE,
    SANITIZE_RGROUP_NAMES,
    SANITIZE_ATOM_MAPS,
    SANITIZE_ADJUST_REACTANTS,
    SANITIZE_MERGEHS,
    SANITIZE_ALL,
);

crate::flags::flag_ops!(SanitizeFlags);

/// Parameters for `sanitize_rxn_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SanitizeRxnParams {
    /// The operations to run.
    pub sanitize_ops: SanitizeFlags,
    /// How `SANITIZE_ADJUST_REACTANTS` adjusts the reactant templates.
    /// Python's `None` is `get_default_adjust_params()`.
    pub params: AdjustQueryParameters,
    /// Returns the operation that failed instead of an error.
    pub catch_errors: bool,
}

impl Default for SanitizeRxnParams {
    fn default() -> SanitizeRxnParams {
        SanitizeRxnParams {
            sanitize_ops: SanitizeFlags::SANITIZE_ALL,
            params: get_default_adjust_params(),
            catch_errors: false,
        }
    }
}

/// Fixes the R groups and atom maps of a reaction drawn in a sketcher,
/// aromatizes its reactant templates, and merges their hydrogens, returning
/// `SANITIZE_NONE`. Atom maps that would overflow an `i32` are
/// `Error::Value`. Mirrors `rdChemReactions.SanitizeRxn`.
pub fn sanitize_rxn(rxn: &mut ChemicalReaction) -> Result<SanitizeFlags> {
    sanitize_rxn_with_params(rxn, &SanitizeRxnParams::default())
}

/// With `catch_errors`, returns the operation that failed. Mirrors
/// `rdChemReactions.SanitizeRxn(rxn, sanitizeOps, params, catchErrors)`.
pub fn sanitize_rxn_with_params(
    rxn: &mut ChemicalReaction,
    params: &SanitizeRxnParams,
) -> Result<SanitizeFlags> {
    Ok(SanitizeFlags(ffi::sanitize_rxn(
        rxn.pin_mut(),
        params.sanitize_ops.bits(),
        &params.params.to_ffi(),
        params.catch_errors,
    )?))
}

/// Sanitizes each template as a molecule. Mirrors
/// `rdChemReactions.SanitizeRxnAsMols`.
pub fn sanitize_rxn_as_mols(rxn: &mut ChemicalReaction) -> Result<()> {
    sanitize_rxn_as_mols_with_sanitize_ops(rxn, MolSanitizeFlags::SANITIZE_ALL)
}

/// Mirrors `rdChemReactions.SanitizeRxnAsMols(rxn, sanitizeOps)`.
pub fn sanitize_rxn_as_mols_with_sanitize_ops(
    rxn: &mut ChemicalReaction,
    sanitize_ops: MolSanitizeFlags,
) -> Result<()> {
    Ok(ffi::sanitize_rxn_as_mols(
        rxn.pin_mut(),
        sanitize_ops.bits(),
    )?)
}

/// The default adjustments for reactant templates: aromatize, nothing more.
/// Mirrors `rdChemReactions.GetDefaultAdjustParams`.
pub fn get_default_adjust_params() -> AdjustQueryParameters {
    AdjustQueryParameters::from_ffi(expect_ok(ffi::get_default_adjust_params()))
}

/// ChemDraw-style adjustments, logging that they are deprecated. Mirrors
/// `rdChemReactions.GetChemDrawRxnAdjustParams`.
pub fn get_chem_draw_rxn_adjust_params() -> AdjustQueryParameters {
    AdjustQueryParameters::from_ffi(expect_ok(ffi::get_chem_draw_rxn_adjust_params()))
}

/// Adjustments that let reactant templates match only at their R groups.
/// Mirrors `rdChemReactions.MatchOnlyAtRgroupsAdjustParams`.
pub fn match_only_at_rgroups_adjust_params() -> AdjustQueryParameters {
    AdjustQueryParameters::from_ffi(expect_ok(ffi::match_only_at_rgroups_adjust_params()))
}
