// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdSLNParse`: molecules and queries from Sybyl line
//! notation.
//!
//! The parser's macro expansion runs `std::regex`, whose recursion
//! overflows the stack when a `{` is followed by tens of thousands of
//! characters (Python crashes too); the process then aborts, as on any
//! stack overflow.

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::ffi::mol::ffi::RWMol;
use crate::ffi::sln_parse::ffi;
use crate::{Error, Result};

fn parsed(ptr: UniquePtr<RWMol>, sln: &str) -> Result<Mol> {
    if ptr.is_null() {
        return Err(Error::Parse {
            message: format!("could not parse SLN '{sln}'"),
        });
    }
    Ok(Mol::from_ffi(ptr))
}

/// Parameters for `mol_from_sln_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolFromSLNParams {
    /// Sanitizes the molecule (removing its explicit hydrogens).
    pub sanitize: bool,
    /// Prints the parser's trace to stderr.
    pub debug_parser: bool,
}

impl Default for MolFromSLNParams {
    fn default() -> MolFromSLNParams {
        MolFromSLNParams {
            sanitize: true,
            debug_parser: false,
        }
    }
}

/// A molecule from SLN with no query information. Mirrors
/// `rdSLNParse.MolFromSLN`. `Error::Parse` for SLN that does not parse
/// (Python's `None`, E2).
pub fn mol_from_sln(sln: &str) -> Result<Mol> {
    mol_from_sln_with_params(sln, &MolFromSLNParams::default())
}

/// Mirrors `rdSLNParse.MolFromSLN(SLN, sanitize, debugParser)`.
pub fn mol_from_sln_with_params(sln: &str, params: &MolFromSLNParams) -> Result<Mol> {
    parsed(
        ffi::mol_from_sln(sln, params.sanitize, params.debug_parser)?,
        sln,
    )
}

/// Parameters for `mol_from_query_sln_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolFromQuerySLNParams {
    /// Merges explicit hydrogens into their heavy atoms' queries.
    pub merge_hs: bool,
    /// Prints the parser's trace to stderr.
    pub debug_parser: bool,
}

impl Default for MolFromQuerySLNParams {
    fn default() -> MolFromQuerySLNParams {
        MolFromQuerySLNParams {
            merge_hs: true,
            debug_parser: false,
        }
    }
}

/// A substructure query from SLN. Mirrors `rdSLNParse.MolFromQuerySLN`.
pub fn mol_from_query_sln(sln: &str) -> Result<Mol> {
    mol_from_query_sln_with_params(sln, &MolFromQuerySLNParams::default())
}

/// Mirrors `rdSLNParse.MolFromQuerySLN(SLN, mergeHs, debugParser)`.
pub fn mol_from_query_sln_with_params(sln: &str, params: &MolFromQuerySLNParams) -> Result<Mol> {
    parsed(
        ffi::mol_from_query_sln(sln, params.merge_hs, params.debug_parser)?,
        sln,
    )
}
