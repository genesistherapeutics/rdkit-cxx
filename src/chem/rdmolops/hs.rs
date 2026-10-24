// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Removing and merging hydrogens.

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::molops::ffi;

/// Which hydrogens `remove_hs_with_params` removes. Mirrors
/// `rdmolops.RemoveHsParameters`, plus `RemoveHs`'s `sanitize` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct RemoveHsParameters {
    /// Removes hydrogens with no bonds.
    pub remove_degree_zero: bool,
    /// Removes hydrogens with two or more bonds.
    pub remove_higher_degrees: bool,
    /// Removes hydrogens bonded only to hydrogens.
    pub remove_only_h_neighbors: bool,
    /// Removes hydrogens with isotopes.
    pub remove_isotopes: bool,
    /// Removes hydrogens with isotopes, recording them on their neighbors (`_isotopicHs`) for `add_hs`.
    pub remove_and_track_isotopes: bool,
    /// Removes hydrogens bonded to dummy atoms.
    pub remove_dummy_neighbors: bool,
    /// Removes hydrogens that define a double bond's stereo.
    pub remove_defining_bond_stereo: bool,
    /// Removes hydrogens with wedged bonds.
    pub remove_with_wedged_bond: bool,
    /// Removes hydrogens with queries.
    pub remove_with_query: bool,
    /// Removes hydrogens with atom map numbers.
    pub remove_mapped: bool,
    /// Removes hydrogens in substance groups.
    pub remove_in_s_groups: bool,
    /// Warns about some kinds of hydrogen it removes.
    pub show_warnings: bool,
    /// Removes hydrogens that are not implicit; deprecated, and `!implicitOnly` in Python's other overload.
    pub remove_nonimplicit: bool,
    /// Deprecated: counts removed hydrogens as explicit.
    pub update_explicit_count: bool,
    /// Removes hydrogens with charge -1.
    pub remove_hydrides: bool,
    /// Removes hydrogens bonded to atoms with non-tetrahedral stereo.
    pub remove_nontetrahedral_neighbors: bool,
    /// Sanitizes the result.
    pub sanitize: bool,
}

impl Default for RemoveHsParameters {
    fn default() -> RemoveHsParameters {
        let data = expect_ok(ffi::remove_hs_parameters_default());
        RemoveHsParameters {
            remove_degree_zero: data.remove_degree_zero,
            remove_higher_degrees: data.remove_higher_degrees,
            remove_only_h_neighbors: data.remove_only_h_neighbors,
            remove_isotopes: data.remove_isotopes,
            remove_and_track_isotopes: data.remove_and_track_isotopes,
            remove_dummy_neighbors: data.remove_dummy_neighbors,
            remove_defining_bond_stereo: data.remove_defining_bond_stereo,
            remove_with_wedged_bond: data.remove_with_wedged_bond,
            remove_with_query: data.remove_with_query,
            remove_mapped: data.remove_mapped,
            remove_in_s_groups: data.remove_in_s_groups,
            show_warnings: data.show_warnings,
            remove_nonimplicit: data.remove_nonimplicit,
            update_explicit_count: data.update_explicit_count,
            remove_hydrides: data.remove_hydrides,
            remove_nontetrahedral_neighbors: data.remove_nontetrahedral_neighbors,
            sanitize: data.sanitize,
        }
    }
}

impl RemoveHsParameters {
    fn to_ffi(&self) -> ffi::RemoveHsParametersData {
        ffi::RemoveHsParametersData {
            remove_degree_zero: self.remove_degree_zero,
            remove_higher_degrees: self.remove_higher_degrees,
            remove_only_h_neighbors: self.remove_only_h_neighbors,
            remove_isotopes: self.remove_isotopes,
            remove_and_track_isotopes: self.remove_and_track_isotopes,
            remove_dummy_neighbors: self.remove_dummy_neighbors,
            remove_defining_bond_stereo: self.remove_defining_bond_stereo,
            remove_with_wedged_bond: self.remove_with_wedged_bond,
            remove_with_query: self.remove_with_query,
            remove_mapped: self.remove_mapped,
            remove_in_s_groups: self.remove_in_s_groups,
            show_warnings: self.show_warnings,
            remove_nonimplicit: self.remove_nonimplicit,
            update_explicit_count: self.update_explicit_count,
            remove_hydrides: self.remove_hydrides,
            remove_nontetrahedral_neighbors: self.remove_nontetrahedral_neighbors,
            sanitize: self.sanitize,
        }
    }
}

/// A copy without the hydrogens that `add_hs` would put back. Mirrors
/// `rdmolops.RemoveHs`.
pub fn remove_hs(mol: &Mol) -> Result<Mol> {
    remove_hs_with_params(mol, &RemoveHsParameters::default())
}

/// Mirrors `rdmolops.RemoveHs(mol, params, sanitize)`, and `RemoveHs(mol,
/// implicitOnly, updateExplicitCount, sanitize)`, whose `implicitOnly` is
/// `!remove_nonimplicit`.
pub fn remove_hs_with_params(mol: &Mol, params: &RemoveHsParameters) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::remove_hs(
        mol.as_ffi(),
        &params.to_ffi(),
    )?))
}

/// A copy without any hydrogens, sanitized. Mirrors `rdmolops.RemoveAllHs`.
pub fn remove_all_hs(mol: &Mol) -> Result<Mol> {
    remove_all_hs_with_sanitize(mol, true)
}

/// Mirrors `rdmolops.RemoveAllHs(mol, sanitize)`.
pub fn remove_all_hs_with_sanitize(mol: &Mol, sanitize: bool) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::remove_all_hs(mol.as_ffi(), sanitize)?))
}

/// Parameters for `merge_query_hs_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MergeQueryHsParams {
    /// Merges only hydrogens without atom map numbers.
    pub merge_unmapped_only: bool,
    /// Merges hydrogens with isotopes too.
    pub merge_isotopes: bool,
}

/// A copy of a query with its explicit hydrogens merged into their
/// neighbors' queries, as `[C;!H0]` for `[C][H]`. Mirrors
/// `rdmolops.MergeQueryHs`.
pub fn merge_query_hs(mol: &Mol) -> Result<Mol> {
    merge_query_hs_with_params(mol, &MergeQueryHsParams::default())
}

/// Mirrors `rdmolops.MergeQueryHs(mol, mergeUnmappedOnly, mergeIsotopes)`.
pub fn merge_query_hs_with_params(mol: &Mol, params: &MergeQueryHsParams) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::merge_query_hs(
        mol.as_ffi(),
        params.merge_unmapped_only,
        params.merge_isotopes,
    )?))
}

/// Whether the query has hydrogen atoms, and whether any of them can't be
/// merged by `merge_query_hs`. Mirrors `rdmolops.HasQueryHs`.
pub fn has_query_hs(mol: &Mol) -> (bool, bool) {
    let pair = expect_ok(ffi::has_query_hs(mol.as_ffi()));
    (pair[0], pair[1])
}
