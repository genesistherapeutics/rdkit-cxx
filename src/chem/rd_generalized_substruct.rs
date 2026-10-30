// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdGeneralizedSubstruct`: substructure searches with
//! extended query molecules, which enumerate link nodes, repeat units, and
//! variable attachment points, and match tautomers.
//!
//! The wrapper checks what RDKit trusts (CONVENTIONS.md U1): the link nodes
//! of a molecule to enumerate, and a pickle's or JSON's molecules, tautomer
//! queries, and base64. Tautomer queries match as `rd_tautomer_query`'s do.

use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::SubstructMatchParameters;
use crate::chem::rdmolops::AdjustQueryParameters;
use crate::data_structs::ExplicitBitVect;
use crate::error::{expect_ok, to_string};
use crate::ffi::generalized_substruct::ffi;
use crate::ffi::geometry::ffi::IndexList;

const DEFAULT_FINGERPRINT_SIZE: u32 = 2048;

/// A query that holds a molecule, a bundle of molecules, a tautomer query,
/// or a bundle of tautomer queries, as `create_extended_query_mol` makes
/// them. Mirrors `rdGeneralizedSubstruct.ExtendedQueryMol`.
///
/// `ExtendedQueryMol` is `Send` but not `Sync`: its pattern fingerprint finds
/// its molecules' rings through a read-only reference, as with `Mol`
/// (CONVENTIONS.md S1).
pub struct ExtendedQueryMol {
    /// Never null.
    ptr: UniquePtr<ffi::ExtendedQueryMol>,
}

impl ExtendedQueryMol {
    fn from_ffi(ptr: UniquePtr<ffi::ExtendedQueryMol>) -> ExtendedQueryMol {
        assert!(!ptr.is_null(), "RDKit returned a null extended query");
        ExtendedQueryMol { ptr }
    }

    pub(crate) fn as_ffi(&self) -> &ffi::ExtendedQueryMol {
        self.ptr.as_ref().unwrap()
    }

    /// A query from a binary pickle (`to_binary`) held in a string, as
    /// `from_data` reads bytes. Mirrors `ExtendedQueryMol(text)`.
    ///
    /// # Safety
    ///
    /// As for `from_data`.
    pub unsafe fn new(text: &str) -> Result<ExtendedQueryMol> {
        // SAFETY: the caller's contract is `from_data`'s.
        unsafe { ExtendedQueryMol::from_data(text.as_bytes()) }
    }

    /// A query from JSON (`to_json`) when `is_json`, else from a binary
    /// pickle. JSON whose molecules do not parse is `Error::Parse`, and
    /// base64 RDKit's decoder would write past its buffer for, `Error::Value`;
    /// malformed JSON is RDKit's `Error::Value`. Mirrors
    /// `ExtendedQueryMol(text, isJSON)`.
    ///
    /// # Safety
    ///
    /// `text` must be JSON (`to_json`) or a pickle (`to_binary`) RDKit
    /// wrote, unmodified. The JSON holds each molecule as a base64 molecule
    /// pickle, which RDKit reads as `Mol::from_pkl_string` does.
    pub unsafe fn new_with_is_json(text: &str, is_json: bool) -> Result<ExtendedQueryMol> {
        if is_json {
            Ok(ExtendedQueryMol::from_ffi(
                ffi::extended_query_mol_from_json(text)?,
            ))
        } else {
            // SAFETY: the caller's contract is `new`'s.
            unsafe { ExtendedQueryMol::new(text) }
        }
    }

    /// A query from a binary pickle (`to_binary`). `Error::Value` for one
    /// with an empty bundle, a tautomer query whose indices its matcher would
    /// read past, or a molecule with a monomer type its class lacks
    /// (CONVENTIONS.md U10); RDKit's errors for one it cannot read. A Rust
    /// addition: Python reads a pickle into a query with `InitFromBinary`.
    ///
    /// # Safety
    ///
    /// `data` must be a pickle RDKit wrote (this crate's `to_binary` or
    /// RDKit's), unmodified: it holds molecule pickles, which RDKit reads as
    /// `Mol::from_pkl_string` does.
    pub unsafe fn from_data(data: &[u8]) -> Result<ExtendedQueryMol> {
        Ok(ExtendedQueryMol::from_ffi(
            ffi::extended_query_mol_from_binary(data)?,
        ))
    }

    /// Replaces the query with one from a binary pickle, as `from_data`
    /// reads; on an error it is unchanged. Python's takes the pickle as a
    /// string. Mirrors `ExtendedQueryMol.InitFromBinary`.
    ///
    /// # Safety
    ///
    /// As for `from_data`.
    pub unsafe fn init_from_binary(&mut self, pkl: &[u8]) -> Result<()> {
        // SAFETY: the caller's contract is `from_data`'s.
        *self = unsafe { ExtendedQueryMol::from_data(pkl) }?;
        Ok(())
    }

    /// Replaces the query with one from JSON, as `new_with_is_json` reads;
    /// on an error it is unchanged. Mirrors `ExtendedQueryMol.InitFromJSON`.
    ///
    /// # Safety
    ///
    /// As for `new_with_is_json`.
    pub unsafe fn init_from_json(&mut self, text: &str) -> Result<()> {
        // SAFETY: the caller's contract is `new_with_is_json`'s.
        *self = unsafe { ExtendedQueryMol::new_with_is_json(text, true) }?;
        Ok(())
    }

    /// A binary pickle of the query. Bundles and tautomer queries need RDKit
    /// built with Boost.Serialization (`Error::Runtime` otherwise). RDKit
    /// writes a tautomer bundle's size in 16 bits, so one of more than 65535
    /// queries is `Error::Value`. Mirrors `ExtendedQueryMol.ToBinary`.
    /// `Error::Value` for a molecule a pickle cannot hold (`Mol::to_binary`).
    pub fn to_binary(&self) -> Result<Vec<u8>> {
        Ok(ffi::extended_query_mol_to_binary(self.as_ffi())?)
    }

    /// The query as JSON: each molecule's pickle in base64, with its SMILES
    /// or SMARTS. Mirrors `ExtendedQueryMol.ToJSON`. `Error::Value` for a
    /// molecule a pickle cannot hold (`Mol::to_binary`).
    pub fn to_json(&self) -> Result<String> {
        to_string(ffi::extended_query_mol_to_json(self.as_ffi())?)
    }

    /// The query's pattern fingerprint, for screening targets (with
    /// `pattern_fingerprint_target`); a bundle's is the bits all its members
    /// set. Mirrors `ExtendedQueryMol.PatternFingerprintQuery`.
    pub fn pattern_fingerprint_query(&self) -> Result<ExplicitBitVect> {
        self.pattern_fingerprint_query_with_fingerprint_size(DEFAULT_FINGERPRINT_SIZE)
    }

    /// `Error::Value` for a size above `i32::MAX` (CONVENTIONS.md U3), and
    /// `Error::Runtime` for 0. Mirrors
    /// `ExtendedQueryMol.PatternFingerprintQuery(fingerprintSize)`.
    pub fn pattern_fingerprint_query_with_fingerprint_size(
        &self,
        fingerprint_size: u32,
    ) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(
            ffi::extended_query_mol_pattern_fingerprint_query(self.as_ffi(), fingerprint_size)?,
        ))
    }
}

/// A copy with its own molecules.
impl Clone for ExtendedQueryMol {
    fn clone(&self) -> ExtendedQueryMol {
        ExtendedQueryMol::from_ffi(expect_ok(ffi::extended_query_mol_clone(self.as_ffi())))
    }
}

impl fmt::Debug for ExtendedQueryMol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ExtendedQueryMol").finish_non_exhaustive()
    }
}

/// Parameters for `create_extended_query_mol_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CreateExtendedQueryMolParams {
    /// Enumerates link nodes, repeat units, and variable attachment points.
    pub do_enumeration: bool,
    /// Makes tautomer queries.
    pub do_tautomers: bool,
    /// Adjusts each result's query properties.
    pub adjust_query_properties: bool,
    /// How to adjust them.
    pub adjust_query_parameters: AdjustQueryParameters,
}

impl Default for CreateExtendedQueryMolParams {
    fn default() -> CreateExtendedQueryMolParams {
        CreateExtendedQueryMolParams {
            do_enumeration: true,
            do_tautomers: true,
            adjust_query_properties: false,
            adjust_query_parameters: AdjustQueryParameters::default(),
        }
    }
}

/// A query from `mol`: its link nodes, repeat units, and variable attachment
/// points enumerated, then each result made a tautomer query, as each step
/// applies. Mirrors `rdGeneralizedSubstruct.CreateExtendedQueryMol`.
pub fn create_extended_query_mol(mol: &Mol) -> Result<ExtendedQueryMol> {
    create_extended_query_mol_with_params(mol, &CreateExtendedQueryMolParams::default())
}

/// `Error::Value` for a link node RDKit's enumerator would read past an
/// empty product list for: one on atoms past the 65th (whose isotope labels
/// wrap), with a bond that is not single, with its two bonds to one atom, or
/// on the atom of another node. Mirrors
/// `rdGeneralizedSubstruct.CreateExtendedQueryMol(mol, doEnumeration,
/// doTautomers, adjustQueryProperties, adjustQueryParameters)`.
pub fn create_extended_query_mol_with_params(
    mol: &Mol,
    params: &CreateExtendedQueryMolParams,
) -> Result<ExtendedQueryMol> {
    Ok(ExtendedQueryMol::from_ffi(ffi::create_extended_query_mol(
        mol.as_ffi(),
        params.do_enumeration,
        params.do_tautomers,
        params.adjust_query_properties,
        &params.adjust_query_parameters.to_ffi(),
    )?))
}

/// Whether `query` matches part of `mol`. Mirrors
/// `rdGeneralizedSubstruct.MolHasSubstructMatch`.
pub fn mol_has_substruct_match(mol: &Mol, query: &ExtendedQueryMol) -> Result<bool> {
    mol_has_substruct_match_with_params(mol, query, &SubstructMatchParameters::default())
}

/// Mirrors `rdGeneralizedSubstruct.MolHasSubstructMatch(mol, query,
/// params)`.
pub fn mol_has_substruct_match_with_params(
    mol: &Mol,
    query: &ExtendedQueryMol,
    params: &SubstructMatchParameters,
) -> Result<bool> {
    Ok(ffi::mol_has_substruct_match_extended(
        mol.as_ffi(),
        query.as_ffi(),
        &params.to_ffi(),
    )?)
}

/// The atom of `mol` matching each query atom in the first match, or empty.
/// Mirrors `rdGeneralizedSubstruct.MolGetSubstructMatch`.
pub fn mol_get_substruct_match(mol: &Mol, query: &ExtendedQueryMol) -> Result<Vec<u32>> {
    mol_get_substruct_match_with_params(mol, query, &SubstructMatchParameters::default())
}

/// Mirrors `rdGeneralizedSubstruct.MolGetSubstructMatch(mol, query,
/// params)`.
pub fn mol_get_substruct_match_with_params(
    mol: &Mol,
    query: &ExtendedQueryMol,
    params: &SubstructMatchParameters,
) -> Result<Vec<u32>> {
    Ok(ffi::mol_get_substruct_match_extended(
        mol.as_ffi(),
        query.as_ffi(),
        &params.to_ffi(),
    )?)
}

/// Every match, as in `mol_get_substruct_match`: of the bundle's first
/// member that matches, for a bundle. Mirrors
/// `rdGeneralizedSubstruct.MolGetSubstructMatches`.
pub fn mol_get_substruct_matches(mol: &Mol, query: &ExtendedQueryMol) -> Result<Vec<Vec<u32>>> {
    mol_get_substruct_matches_with_params(mol, query, &SubstructMatchParameters::default())
}

/// Mirrors `rdGeneralizedSubstruct.MolGetSubstructMatches(mol, query,
/// params)`.
pub fn mol_get_substruct_matches_with_params(
    mol: &Mol,
    query: &ExtendedQueryMol,
    params: &SubstructMatchParameters,
) -> Result<Vec<Vec<u32>>> {
    Ok(IndexList::into_nested(
        ffi::mol_get_substruct_matches_extended(mol.as_ffi(), query.as_ffi(), &params.to_ffi())?,
    ))
}

/// A target's pattern fingerprint, for screening it against a query's
/// (`ExtendedQueryMol::pattern_fingerprint_query`). Mirrors
/// `rdGeneralizedSubstruct.PatternFingerprintTarget`.
pub fn pattern_fingerprint_target(target: &Mol) -> Result<ExplicitBitVect> {
    pattern_fingerprint_target_with_fingerprint_size(target, DEFAULT_FINGERPRINT_SIZE)
}

/// `Error::Value` for a size above `i32::MAX` (CONVENTIONS.md U3), and
/// `Error::Runtime` for 0. Mirrors
/// `rdGeneralizedSubstruct.PatternFingerprintTarget(target,
/// fingerprintSize)`.
pub fn pattern_fingerprint_target_with_fingerprint_size(
    target: &Mol,
    fingerprint_size: u32,
) -> Result<ExplicitBitVect> {
    Ok(ExplicitBitVect::from_ffi(ffi::pattern_fingerprint_target(
        target.as_ffi(),
        fingerprint_size,
    )?))
}
