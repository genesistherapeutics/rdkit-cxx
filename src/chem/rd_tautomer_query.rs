// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdTautomerQuery`: substructure queries that match any
//! tautomer of the query molecule.
//!
//! A query is a template, in which the atoms and bonds tautomerism changes
//! match more widely, and the tautomers a match of the template must agree
//! with at those atoms and bonds. The wrapper checks what RDKit trusts
//! (CONVENTIONS.md U1): a tautomer transform file's bond and charge lists,
//! and a pickle's tautomers and indices. Its matcher is RDKit's, except that
//! the matches a recursive SMARTS in the template makes of its own pass
//! unchecked, where RDKit's reads past their end.

use std::fmt;
use std::path::Path;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::SubstructMatchParameters;
use crate::chem::rdmolops::mol_list_into_vec;
use crate::data_structs::ExplicitBitVect;
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::tautomer_query::ffi;
use crate::util::path_bytes;

const DEFAULT_FINGERPRINT_SIZE: u32 = 2048;

/// A substructure query that matches any tautomer of a molecule. Mirrors
/// `rdTautomerQuery.TautomerQuery`.
///
/// `TautomerQuery` is `Send` but not `Sync`: its pattern fingerprint finds the
/// template's rings through a read-only reference, as with `Mol`
/// (CONVENTIONS.md S1).
pub struct TautomerQuery {
    /// Never null.
    ptr: UniquePtr<ffi::TautomerQuery>,
}

impl TautomerQuery {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::TautomerQuery>) -> TautomerQuery {
        assert!(!ptr.is_null(), "RDKit returned a null tautomer query");
        TautomerQuery { ptr }
    }

    pub(crate) fn as_ffi(&self) -> &ffi::TautomerQuery {
        self.ptr.as_ref().unwrap()
    }

    /// A query for the tautomers of `mol` that RDKit's default transforms
    /// make. Mirrors `TautomerQuery(mol)`.
    pub fn new(mol: &Mol) -> Result<TautomerQuery> {
        Ok(TautomerQuery::from_ffi(ffi::tautomer_query_from_mol(
            mol.as_ffi(),
            b"",
        )?))
    }

    /// A query for the tautomers the transforms in a file make; an empty path
    /// uses RDKit's defaults. `Error::Os` for a file that cannot be opened,
    /// `Error::Value` for a transform whose SMARTS does not parse, or whose
    /// bond types or charges (when given) are fewer than its bonds or atoms,
    /// which RDKit would read past. Mirrors
    /// `TautomerQuery(mol, tautomerTransformFile)`.
    pub fn new_with_tautomer_transform_file<P: AsRef<Path>>(
        mol: &Mol,
        tautomer_transform_file: P,
    ) -> Result<TautomerQuery> {
        Ok(TautomerQuery::from_ffi(ffi::tautomer_query_from_mol(
            mol.as_ffi(),
            path_bytes(tautomer_transform_file.as_ref()),
        )?))
    }

    /// A query from a pickle (`to_binary`). `Error::Value` for one whose
    /// tautomers lack the template's atoms or whose modified atoms or bonds
    /// are out of range, which RDKit's matcher would read past; RDKit's
    /// errors for one it cannot read. Mirrors `TautomerQuery(pickle)`.
    ///
    /// # Safety
    ///
    /// `pickle` must be a query pickle RDKit wrote (this crate's `to_binary`
    /// or RDKit's), unmodified: it holds molecule pickles, which RDKit reads
    /// as `Mol::from_pkl_string` does.
    pub unsafe fn from_pickle(pickle: &[u8]) -> Result<TautomerQuery> {
        Ok(TautomerQuery::from_ffi(ffi::tautomer_query_from_pickle(
            pickle,
        )?))
    }

    /// Whether the query matches part of `target`. Mirrors
    /// `TautomerQuery.IsSubstructOf`.
    pub fn is_substruct_of(&self, target: &Mol) -> Result<bool> {
        self.is_substruct_of_with_params(target, &SubstructMatchParameters::default())
    }

    /// Mirrors `TautomerQuery.IsSubstructOf(target, params)` and its keyword
    /// form.
    pub fn is_substruct_of_with_params(
        &self,
        target: &Mol,
        params: &SubstructMatchParameters,
    ) -> Result<bool> {
        Ok(ffi::tautomer_query_is_substruct_of(
            self.as_ffi(),
            target.as_ffi(),
            &params.to_ffi(),
        )?)
    }

    /// The target atom matching each template atom in the first match, or
    /// empty. Mirrors `TautomerQuery.GetSubstructMatch`.
    pub fn get_substruct_match(&self, target: &Mol) -> Result<Vec<u32>> {
        self.get_substruct_match_with_params(target, &SubstructMatchParameters::default())
    }

    /// Mirrors `TautomerQuery.GetSubstructMatch(target, params)` and its
    /// keyword form.
    pub fn get_substruct_match_with_params(
        &self,
        target: &Mol,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<u32>> {
        Ok(ffi::tautomer_query_get_substruct_match(
            self.as_ffi(),
            target.as_ffi(),
            &params.to_ffi(),
        )?)
    }

    /// Every match, as in `get_substruct_match`. Mirrors
    /// `TautomerQuery.GetSubstructMatches`.
    pub fn get_substruct_matches(&self, target: &Mol) -> Result<Vec<Vec<u32>>> {
        self.get_substruct_matches_with_params(target, &SubstructMatchParameters::default())
    }

    /// Mirrors `TautomerQuery.GetSubstructMatches(target, params)` and its
    /// keyword form.
    pub fn get_substruct_matches_with_params(
        &self,
        target: &Mol,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<Vec<u32>>> {
        Ok(IndexList::into_nested(
            ffi::tautomer_query_get_substruct_matches(
                self.as_ffi(),
                target.as_ffi(),
                &params.to_ffi(),
            )?,
        ))
    }

    /// Every match with a copy of the tautomer it matched. Mirrors
    /// `TautomerQuery.GetSubstructMatchesWithTautomers`.
    pub fn get_substruct_matches_with_tautomers(
        &self,
        target: &Mol,
    ) -> Result<Vec<(Vec<u32>, Mol)>> {
        self.get_substruct_matches_with_tautomers_with_params(
            target,
            &SubstructMatchParameters::default(),
        )
    }

    /// Mirrors `TautomerQuery.GetSubstructMatchesWithTautomers(target,
    /// params)` and its keyword form.
    pub fn get_substruct_matches_with_tautomers_with_params(
        &self,
        target: &Mol,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<(Vec<u32>, Mol)>> {
        let mut matches = Vec::new();
        let tautomers = ffi::tautomer_query_get_substruct_matches_with_tautomers(
            self.as_ffi(),
            target.as_ffi(),
            &params.to_ffi(),
            &mut matches,
        )?;
        Ok(IndexList::into_nested(matches)
            .into_iter()
            .zip(mol_list_into_vec(tautomers))
            .collect())
    }

    /// The template's pattern fingerprint, for screening targets (with
    /// `pattern_fingerprint_tautomer_target`). Mirrors
    /// `TautomerQuery.PatternFingerprintTemplate`.
    pub fn pattern_fingerprint_template(&self) -> Result<ExplicitBitVect> {
        self.pattern_fingerprint_template_with_fingerprint_size(DEFAULT_FINGERPRINT_SIZE)
    }

    /// `Error::Value` for a size above `i32::MAX` (CONVENTIONS.md U3), and
    /// `Error::Runtime` for 0. Mirrors
    /// `TautomerQuery.PatternFingerprintTemplate(fingerprintSize)`.
    pub fn pattern_fingerprint_template_with_fingerprint_size(
        &self,
        fingerprint_size: u32,
    ) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(
            ffi::tautomer_query_pattern_fingerprint_template(self.as_ffi(), fingerprint_size)?,
        ))
    }

    /// A copy of the template the query matches first. Python returns the
    /// query's own molecule. Mirrors `TautomerQuery.GetTemplateMolecule`.
    pub fn get_template_molecule(&self) -> Mol {
        Mol::from_ffi(expect_ok(ffi::tautomer_query_get_template_molecule(
            self.as_ffi(),
        )))
    }

    /// The atoms tautomerism changes. Mirrors
    /// `TautomerQuery.GetModifiedAtoms`.
    pub fn get_modified_atoms(&self) -> Vec<u32> {
        expect_ok(ffi::tautomer_query_get_modified_atoms(self.as_ffi()))
    }

    /// The bonds tautomerism changes. Mirrors
    /// `TautomerQuery.GetModifiedBonds`.
    pub fn get_modified_bonds(&self) -> Vec<u32> {
        expect_ok(ffi::tautomer_query_get_modified_bonds(self.as_ffi()))
    }

    /// Copies of the tautomers. Mirrors `TautomerQuery.GetTautomers`.
    pub fn get_tautomers(&self) -> Vec<Mol> {
        mol_list_into_vec(expect_ok(ffi::tautomer_query_get_tautomers(self.as_ffi())))
    }

    /// A pickle of the query, a Boost.Serialization text archive.
    /// `Error::Runtime` if RDKit was built without Boost.Serialization
    /// (`tautomer_query_can_serialize`). Mirrors `TautomerQuery.ToBinary`.
    /// `Error::Value` for a template or tautomer a pickle cannot hold
    /// (`Mol::to_binary`).
    pub fn to_binary(&self) -> Result<Vec<u8>> {
        Ok(ffi::tautomer_query_to_binary(self.as_ffi())?)
    }
}

/// A copy, with its own template and tautomers.
impl Clone for TautomerQuery {
    fn clone(&self) -> TautomerQuery {
        TautomerQuery::from_ffi(expect_ok(ffi::tautomer_query_clone(self.as_ffi())))
    }
}

impl fmt::Debug for TautomerQuery {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TautomerQuery")
            .field("modified_atoms", &self.get_modified_atoms())
            .field("modified_bonds", &self.get_modified_bonds())
            .finish_non_exhaustive()
    }
}

/// A target's pattern fingerprint, for screening it against a query's
/// (`TautomerQuery::pattern_fingerprint_template`). Mirrors
/// `rdTautomerQuery.PatternFingerprintTautomerTarget`.
pub fn pattern_fingerprint_tautomer_target(target: &Mol) -> Result<ExplicitBitVect> {
    pattern_fingerprint_tautomer_target_with_fingerprint_size(target, DEFAULT_FINGERPRINT_SIZE)
}

/// `Error::Value` for a size above `i32::MAX` (CONVENTIONS.md U3), and
/// `Error::Runtime` for 0. Mirrors
/// `rdTautomerQuery.PatternFingerprintTautomerTarget(target,
/// fingerprintSize)`.
pub fn pattern_fingerprint_tautomer_target_with_fingerprint_size(
    target: &Mol,
    fingerprint_size: u32,
) -> Result<ExplicitBitVect> {
    Ok(ExplicitBitVect::from_ffi(
        ffi::pattern_fingerprint_tautomer_target(target.as_ffi(), fingerprint_size)?,
    ))
}

/// Whether RDKit was built with Boost.Serialization, which pickles need.
/// Mirrors `rdTautomerQuery.TautomerQueryCanSerialize`.
pub fn tautomer_query_can_serialize() -> bool {
    expect_ok(ffi::tautomer_query_can_serialize())
}
