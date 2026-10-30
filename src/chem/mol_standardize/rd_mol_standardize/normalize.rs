// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Normalizing functional groups.

use std::path::Path;

use cxx::UniquePtr;

use super::CleanupParameters;
use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::standardize_normalize::ffi;
use crate::util::path_bytes;

/// Applies normalization transformations, each until it no longer changes
/// the molecule, starting the series again after any change. Mirrors
/// `rdMolStandardize.Normalizer`.
pub struct Normalizer {
    /// Never null.
    ptr: UniquePtr<ffi::Normalizer>,
}

impl std::fmt::Debug for Normalizer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Normalizer").finish_non_exhaustive()
    }
}

/// RDKit's built-in transformations and at most 200 restarts. Mirrors
/// `Normalizer()`.
impl Default for Normalizer {
    fn default() -> Normalizer {
        Normalizer::from_ffi(expect_ok(ffi::normalizer_default()))
    }
}

impl Normalizer {
    fn from_ffi(ptr: UniquePtr<ffi::Normalizer>) -> Normalizer {
        assert!(!ptr.is_null(), "RDKit returned a null Normalizer");
        Normalizer { ptr }
    }

    fn raw(&self) -> &ffi::Normalizer {
        self.ptr.as_ref().unwrap()
    }

    /// The transformations of a file (a name and reaction SMARTS,
    /// tab-separated, one per line; empty for RDKit's built-in ones), with
    /// at most `max_restarts` restarts. The file is read on each call and
    /// RDKit parses the text read (Python's RDKit reads a file once per
    /// process and keeps it, under its name). `Error::Os` for a file that
    /// cannot be read, and `Error::Parse` or `Error::Value` for reaction
    /// SMARTS that do not parse. Mirrors `Normalizer(normalizeFilename,
    /// maxRestarts)`.
    pub fn new<P: AsRef<Path>>(normalize_filename: P, max_restarts: u32) -> Result<Normalizer> {
        Ok(Normalizer::from_ffi(ffi::normalizer_new(
            path_bytes(normalize_filename.as_ref()),
            max_restarts,
        )?))
    }

    /// A normalized copy of the molecule: each fragment is normalized, and
    /// of a transformation's products the first by SMILES is kept. Mirrors
    /// `Normalizer.normalize`.
    pub fn normalize(&self, mol: &Mol) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::normalizer_normalize(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// Normalizes the molecule, as one, in place. `Error::Value` unless
    /// every transformation has one reactant and one product with no more
    /// atoms. Mirrors `Normalizer.normalizeInPlace`.
    pub fn normalize_in_place(&self, mol: &mut Mol) -> Result<()> {
        Ok(ffi::normalizer_normalize_in_place(
            self.raw(),
            mol.pin_mut(),
        )?)
    }
}

/// A normalizer from transformations, one per line: a name and reaction
/// SMARTS, tab-separated (lines starting `//` are comments), with
/// `params.max_restarts` restarts. `Error::Parse` or `Error::Value` for
/// reaction SMARTS that do not parse. Mirrors
/// `rdMolStandardize.NormalizerFromData`.
pub fn normalizer_from_data(param_data: &str, params: &CleanupParameters) -> Result<Normalizer> {
    Ok(Normalizer::from_ffi(ffi::normalizer_from_data(
        param_data,
        &params.to_ffi(),
    )?))
}

/// The normalizer the parameters describe: their normalization data, or
/// else their normalizations file, which is `Error::Value` if it holds no
/// transformations (see `CleanupParameters`). Mirrors
/// `rdMolStandardize.NormalizerFromParams`.
pub fn normalizer_from_params(params: &CleanupParameters) -> Result<Normalizer> {
    Ok(Normalizer::from_ffi(ffi::normalizer_from_params(
        &params.to_ffi(),
    )?))
}
