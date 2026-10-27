// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Result;
use crate::chem::rdchem::{Mol, MolBundle};
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::substruct::ffi;

/// Parameters for substructure matching. Mirrors
/// `rdchem.SubstructMatchParameters`, without its callback setter.
#[derive(Clone, Debug, PartialEq)]
pub struct SubstructMatchParameters {
    /// Requires matching chirality.
    pub use_chirality: bool,
    /// Uses enhanced stereochemistry (with `use_chirality`).
    pub use_enhanced_stereo: bool,
    /// Lets aromatic and conjugated bonds match each other.
    pub aromatic_matches_conjugated: bool,
    /// Matches generic atoms in the query.
    pub use_generic_matchers: bool,
    /// Matches queries against queries.
    pub use_query_query_matches: bool,
    /// Allows recursive queries.
    pub recursion_possible: bool,
    /// Drops matches with the same atoms in another order.
    pub uniquify: bool,
    /// The most matches to return.
    pub max_matches: u32,
    /// The most matches to consider for recursive queries.
    pub max_recursive_matches: u32,
    /// Threads to use; 0 or less uses every core.
    pub num_threads: i32,
    /// Atom properties whose values must match.
    pub atom_properties: Vec<String>,
    /// Bond properties whose values must match.
    pub bond_properties: Vec<String>,
    /// Lets specified stereo in the query match unspecified stereo.
    pub specified_stereo_query_matches_unspecified: bool,
}

impl Default for SubstructMatchParameters {
    fn default() -> SubstructMatchParameters {
        SubstructMatchParameters::from_ffi(expect_ok(ffi::substruct_match_parameters_default()))
    }
}

impl SubstructMatchParameters {
    pub(crate) fn from_ffi(data: ffi::SubstructMatchParametersData) -> SubstructMatchParameters {
        SubstructMatchParameters {
            use_chirality: data.use_chirality,
            use_enhanced_stereo: data.use_enhanced_stereo,
            aromatic_matches_conjugated: data.aromatic_matches_conjugated,
            use_generic_matchers: data.use_generic_matchers,
            use_query_query_matches: data.use_query_query_matches,
            recursion_possible: data.recursion_possible,
            uniquify: data.uniquify,
            max_matches: data.max_matches,
            max_recursive_matches: data.max_recursive_matches,
            num_threads: data.num_threads,
            atom_properties: data.atom_properties,
            bond_properties: data.bond_properties,
            specified_stereo_query_matches_unspecified: data
                .specified_stereo_query_matches_unspecified,
        }
    }

    pub(crate) fn to_ffi(&self) -> ffi::SubstructMatchParametersData {
        ffi::SubstructMatchParametersData {
            use_chirality: self.use_chirality,
            use_enhanced_stereo: self.use_enhanced_stereo,
            aromatic_matches_conjugated: self.aromatic_matches_conjugated,
            use_generic_matchers: self.use_generic_matchers,
            use_query_query_matches: self.use_query_query_matches,
            recursion_possible: self.recursion_possible,
            uniquify: self.uniquify,
            max_matches: self.max_matches,
            max_recursive_matches: self.max_recursive_matches,
            num_threads: self.num_threads,
            atom_properties: self.atom_properties.clone(),
            bond_properties: self.bond_properties.clone(),
            specified_stereo_query_matches_unspecified: self
                .specified_stereo_query_matches_unspecified,
        }
    }
}

impl Mol {
    /// Whether `query` matches part of the molecule. Ring queries need the
    /// target's rings found. Mirrors `Mol.HasSubstructMatch`.
    pub fn has_substruct_match(&self, query: &Mol) -> Result<bool> {
        self.has_substruct_match_with_params(query, &SubstructMatchParameters::default())
    }

    /// Mirrors `Mol.HasSubstructMatch(query, params)` and its keyword form.
    pub fn has_substruct_match_with_params(
        &self,
        query: &Mol,
        params: &SubstructMatchParameters,
    ) -> Result<bool> {
        Ok(ffi::mol_has_substruct_match(
            self.as_ffi(),
            query.as_ffi(),
            &params.to_ffi(),
        )?)
    }

    /// The target atom matching each query atom in the first match, or empty.
    /// Mirrors `Mol.GetSubstructMatch`.
    pub fn get_substruct_match(&self, query: &Mol) -> Result<Vec<u32>> {
        self.get_substruct_match_with_params(query, &SubstructMatchParameters::default())
    }

    /// Mirrors `Mol.GetSubstructMatch(query, params)` and its keyword form.
    pub fn get_substruct_match_with_params(
        &self,
        query: &Mol,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<u32>> {
        Ok(ffi::mol_get_substruct_match(
            self.as_ffi(),
            query.as_ffi(),
            &params.to_ffi(),
        )?)
    }

    /// Every match, as in `get_substruct_match`. Mirrors
    /// `Mol.GetSubstructMatches`.
    pub fn get_substruct_matches(&self, query: &Mol) -> Result<Vec<Vec<u32>>> {
        self.get_substruct_matches_with_params(query, &SubstructMatchParameters::default())
    }

    /// Mirrors `Mol.GetSubstructMatches(query, params)` and its keyword form.
    pub fn get_substruct_matches_with_params(
        &self,
        query: &Mol,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<Vec<u32>>> {
        Ok(IndexList::into_nested(ffi::mol_get_substruct_matches(
            self.as_ffi(),
            query.as_ffi(),
            &params.to_ffi(),
        )?))
    }
}

impl Mol {
    /// Whether any molecule in `query` matches part of the molecule. Mirrors
    /// `Mol.HasSubstructMatch(query: MolBundle)` (A5).
    pub fn has_substruct_match_mol_bundle(&self, query: &MolBundle) -> Result<bool> {
        self.has_substruct_match_mol_bundle_with_params(query, &SubstructMatchParameters::default())
    }

    /// Mirrors `Mol.HasSubstructMatch(query: MolBundle, params)` and its
    /// keyword form.
    pub fn has_substruct_match_mol_bundle_with_params(
        &self,
        query: &MolBundle,
        params: &SubstructMatchParameters,
    ) -> Result<bool> {
        let params = params.to_ffi();
        Ok(query.with_ffi(|query| {
            ffi::mol_has_substruct_match_mol_bundle(self.as_ffi(), query, &params)
        })?)
    }

    /// The first match of the first molecule in `query` that matches. Mirrors
    /// `Mol.GetSubstructMatch(query: MolBundle)` (A5).
    pub fn get_substruct_match_mol_bundle(&self, query: &MolBundle) -> Result<Vec<u32>> {
        self.get_substruct_match_mol_bundle_with_params(query, &SubstructMatchParameters::default())
    }

    /// Mirrors `Mol.GetSubstructMatch(query: MolBundle, params)` and its
    /// keyword form.
    pub fn get_substruct_match_mol_bundle_with_params(
        &self,
        query: &MolBundle,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<u32>> {
        let params = params.to_ffi();
        Ok(query.with_ffi(|query| {
            ffi::mol_get_substruct_match_mol_bundle(self.as_ffi(), query, &params)
        })?)
    }

    /// Every match of the first molecule in `query` that matches. Mirrors
    /// `Mol.GetSubstructMatches(query: MolBundle)` (A5).
    pub fn get_substruct_matches_mol_bundle(&self, query: &MolBundle) -> Result<Vec<Vec<u32>>> {
        self.get_substruct_matches_mol_bundle_with_params(
            query,
            &SubstructMatchParameters::default(),
        )
    }

    /// Mirrors `Mol.GetSubstructMatches(query: MolBundle, params)` and its
    /// keyword form.
    pub fn get_substruct_matches_mol_bundle_with_params(
        &self,
        query: &MolBundle,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<Vec<u32>>> {
        let params = params.to_ffi();
        Ok(IndexList::into_nested(query.with_ffi(|query| {
            ffi::mol_get_substruct_matches_mol_bundle(self.as_ffi(), query, &params)
        })?))
    }
}
