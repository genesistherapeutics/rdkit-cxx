// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.MolKey.MolKey`, a Python module: molecule keys,
//! hashes of an InChI and a stereo category, for molecules that the Avalon
//! structure checker has checked.
//!
//! The functions that check structures (`CheckCTAB`, `GetInchiForCTAB`,
//! `GetKeyForCTAB`, `initStruchk`) need the Avalon structure checker, which
//! is not bound (see `pyAvalonTools.InitializeCheckMol` in the manifest), so
//! they are not ported; the module's codes, results, and helpers are.
//! Python's exceptions are `Error::Value`.

use crate::avalon::py_avalon_tools::StruChkResult;

/// The version of the keys. Mirrors `MolKey.MOL_KEY_VERSION`.
pub fn mol_key_version() -> &'static str {
    "1"
}

/// The structure checker's error bits, and the module's own, by name, in
/// Python's order. Mirrors `MolKey.ERROR_DICT`.
pub fn error_dict() -> Vec<(&'static str, u32)> {
    vec![
        ("BAD_MOLECULE", 1),
        ("ALIAS_CONVERSION_FAILED", 2),
        ("TRANSFORMED", 4),
        ("FRAGMENTS_FOUND", 8),
        ("EITHER_WARNING", 16),
        ("STEREO_ERROR", 32),
        ("DUBIOUS_STEREO_REMOVED", 64),
        ("ATOM_CLASH", 128),
        ("ATOM_CHECK_FAILED", 256),
        ("SIZE_CHECK_FAILED", 512),
        ("RECHARGED", 1024),
        ("STEREO_FORCED_BAD", 2048),
        ("STEREO_TRANSFORMED", 4096),
        ("TEMPLATE_TRANSFORMED", 8192),
        ("INCHI_COMPUTATION_ERROR", 65536),
        ("RDKIT_CONVERSION_ERROR", 131072),
        ("INCHI_READWRITE_ERROR", 262144),
        ("NULL_MOL", 524288),
    ]
}

/// InChI failed. Mirrors `MolKey.INCHI_COMPUTATION_ERROR`.
pub fn inchi_computation_error() -> u32 {
    65536
}
/// RDKit could not read the structure. Mirrors
/// `MolKey.RDKIT_CONVERSION_ERROR`.
pub fn rdkit_conversion_error() -> u32 {
    131072
}
/// InChI warned or failed. Mirrors `MolKey.INCHI_READWRITE_ERROR`.
pub fn inchi_readwrite_error() -> u32 {
    262144
}
/// The structure is empty. Mirrors `MolKey.NULL_MOL`.
pub fn null_mol() -> u32 {
    524288
}
/// The position of the chiral flag in a mol block's counts line. Mirrors
/// `MolKey.CHIRAL_POS`.
pub fn chiral_pos() -> usize {
    12
}

/// The bits that make a molecule bad: the structure checker's `bad_set` and
/// the module's errors. Mirrors `MolKey.BAD_SET`.
pub fn bad_set() -> u32 {
    StruChkResult::bad_set.bits()
        | inchi_computation_error()
        | rdkit_conversion_error()
        | inchi_readwrite_error()
        | null_mol()
}

/// The result for an empty structure: `null_mol()` and no structure. Mirrors
/// `MolKey.T_NULL_MOL`.
pub fn t_null_mol() -> (u32, &'static str) {
    (null_mol(), "")
}

/// The pattern of an InChI's stereo layers, uncompiled: the crate has no
/// regular expressions. The module's private racemate helper uses it.
/// Mirrors `MolKey.GET_STEREO_RE`.
pub fn get_stereo_re() -> &'static str {
    r"^InChI=1S(.*?)/(t.*?)/m\d/s1(.*$)"
}

/// The pattern of SMILES for no structure, uncompiled and case-insensitive.
/// Mirrors `MolKey.NULL_SMILES_RE`.
pub fn null_smiles_re() -> &'static str {
    r"^\s*$|^\s*NO_STRUCTURE\s*$"
}

/// The pattern of a mol block's counts line with no atoms. Mirrors
/// `MolKey.PATTERN_NULL_MOL`.
pub fn pattern_null_mol() -> &'static str {
    r"^([\s0]+[1-9]+[\s]+V[\w]*)"
}

/// The stereo categories, by code, in Python's order. Mirrors
/// `MolKey.stereo_code_dict`.
pub fn stereo_code_dict() -> Vec<(&'static str, u32)> {
    vec![
        ("DEFAULT", 0),
        ("S_ACHIR", 1),
        ("S_ABS", 2),
        ("S_REL", 3),
        ("S_PART", 4),
        ("S_UNKN", 5),
        ("S_ABS_ACHIR", 6),
        ("R_ONE", 11),
        ("R_REL", 12),
        ("R_OTHER", 13),
        ("MX_ENANT", 21),
        ("MX_DIAST", 22),
        ("MX_SP2", 31),
        ("MX_DIAST_ABS", 32),
        ("MX_DIAST_REL", 33),
        ("OTHER", 100),
        ("UNDEFINED", 200),
    ]
}

/// The names of the error bits set in `err`, in `error_dict`'s order.
/// Mirrors `MolKey.ErrorBitsToText`.
pub fn error_bits_to_text(err: u32) -> Vec<&'static str> {
    error_dict()
        .into_iter()
        .filter(|&(_, v)| err & v > 0)
        .map(|(k, _)| k)
        .collect()
}

/// An InChI from a checked structure. Mirrors the `InchiResult` named
/// tuple `MolKey.GetInchiForCTAB` returns.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct InchiResult {
    /// The error bits.
    pub error: u32,
    /// The InChI; `None` when the structure is bad.
    pub inchi: Option<String>,
    /// The structure as the checker fixed it.
    pub fixed_ctab: String,
}

/// A molecule key and how it was made. Mirrors the `MolKeyResult` named
/// tuple `MolKey.GetKeyForCTAB` returns.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct MolKeyResult {
    /// The key: the version, `|`, and the Base64 of an MD5 hash.
    pub mol_key: Option<String>,
    /// The error bits.
    pub error: u32,
    /// The InChI.
    pub inchi: String,
    /// The structure as the checker fixed it.
    pub fixed_ctab: String,
    /// The stereo category.
    pub stereo_code: Option<String>,
    /// More stereo information.
    pub stereo_comment: Option<String>,
}
