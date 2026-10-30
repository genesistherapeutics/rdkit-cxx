// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.inchi`, a Python module over `rdinchi` that
//! `rdkit.Chem` re-exports, as `rdkit_cxx::chem` re-exports this one: InChI
//! conversions that log InChI's warnings and errors to RDKit's log, as
//! Python's `RDLogger` does, and can treat them as errors.
//!
//! Python takes a `logging` level for InChI's messages on success; here it
//! is an [`InchiLogLevel`]. Python's `InchiReadWriteError` is
//! `Error::Value` with InChI's message: Python's exception also carries the
//! molecule or InChI made despite the warning, which `rdinchi` returns.

use crate::chem::Mol;
use crate::chem::rdinchi::{self, InchiToMolParams};
use crate::rd_base::log_message;
use crate::{Error, Result};

/// Whether the crate was built with InChI support: always. Mirrors
/// `inchi.INCHI_AVAILABLE`.
pub fn inchi_available() -> bool {
    true
}

/// The level at which InChI's messages are logged on success, as Python's
/// `logging` levels choose `RDLogger`'s functions
/// (`inchi.logLevelToLogFunctionLookup`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum InchiLogLevel {
    /// `logging.DEBUG`: `rdApp.debug`, as `DEBUG: ...`.
    Debug,
    /// `logging.INFO`: `rdApp.info`, as `INFO: ...`.
    Info,
    /// `logging.WARNING`: `rdApp.warning`, as `WARNING: ...`.
    Warning,
    /// `logging.ERROR`: `rdApp.error`, as `ERROR: ...`.
    Error,
    /// `logging.CRITICAL`: `rdApp.error`, as `CRITICAL: ...`.
    Critical,
}

impl InchiLogLevel {
    /// Logs as `RDLogger.logger()`'s method for the level does.
    fn log(self, msg: &str) {
        let (dest, prefix) = match self {
            InchiLogLevel::Debug => ("rdApp.debug", "DEBUG"),
            InchiLogLevel::Info => ("rdApp.info", "INFO"),
            InchiLogLevel::Warning => ("rdApp.warning", "WARNING"),
            InchiLogLevel::Error => ("rdApp.error", "ERROR"),
            InchiLogLevel::Critical => ("rdApp.error", "CRITICAL"),
        };
        log_message(dest, &format!("{prefix}: {msg}\n"));
    }
}

/// Logs InChI's message as Python's functions do.
fn log_result(retcode: i32, message: &str, log_level: Option<InchiLogLevel>) {
    if let Some(level) = log_level
        && retcode == 0
    {
        level.log(message);
    }
    if retcode != 0 {
        if retcode == 1 {
            InchiLogLevel::Warning.log(message);
        } else {
            InchiLogLevel::Error.log(message);
        }
    }
}

fn read_write_error(message: String) -> Error {
    Error::Value { message }
}

/// Parameters for `mol_from_inchi_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MolFromInchiParams {
    /// Sanitizes the molecule.
    pub sanitize: bool,
    /// Removes hydrogens, when `sanitize` is set.
    pub remove_hs: bool,
    /// Also logs InChI's message on success, at this level.
    pub log_level: Option<InchiLogLevel>,
    /// Makes a warning from InChI `Error::Value`, with its message.
    pub treat_warning_as_error: bool,
}

impl Default for MolFromInchiParams {
    fn default() -> MolFromInchiParams {
        MolFromInchiParams {
            sanitize: true,
            remove_hs: true,
            log_level: None,
            treat_warning_as_error: false,
        }
    }
}

/// The molecule an InChI describes; `None` when InChI cannot read it, or
/// when RDKit cannot sanitize it (Python's `ValueError`, which it logs).
/// InChI's warnings and errors go to RDKit's log. Mirrors
/// `inchi.MolFromInchi(inchi)`.
pub fn mol_from_inchi(inchi: &str) -> Result<Option<Mol>> {
    mol_from_inchi_with_params(inchi, &MolFromInchiParams::default())
}

/// As [`mol_from_inchi`]. Mirrors `inchi.MolFromInchi(inchi, sanitize,
/// removeHs, logLevel, treatWarningAsError)`.
pub fn mol_from_inchi_with_params(inchi: &str, params: &MolFromInchiParams) -> Result<Option<Mol>> {
    let result = match rdinchi::inchi_to_mol_with_params(
        inchi,
        &InchiToMolParams {
            sanitize: params.sanitize,
            remove_hs: params.remove_hs,
        },
    ) {
        Ok(result) => result,
        Err(
            e @ (Error::Value { .. }
            | Error::AtomValence { .. }
            | Error::AtomKekulize { .. }
            | Error::Kekulize { .. }
            | Error::MolSanitize { .. }),
        ) => {
            InchiLogLevel::Error.log(&e.to_string());
            return Ok(None);
        }
        Err(e) => return Err(e),
    };
    log_result(result.retcode, &result.message, params.log_level);
    if params.treat_warning_as_error && result.retcode != 0 {
        return Err(read_write_error(result.message));
    }
    Ok(result.mol)
}

/// Parameters for `mol_to_inchi_with_params` and
/// `mol_to_inchi_and_aux_info_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MolToInchiParams {
    /// InChI's options, as `rdinchi::mol_to_inchi_with_options` takes them.
    pub options: String,
    /// Also logs InChI's message on success, at this level.
    pub log_level: Option<InchiLogLevel>,
    /// Makes a warning from InChI `Error::Value`, with its message.
    pub treat_warning_as_error: bool,
}

/// The molecule's standard InChI and its AuxInfo. InChI's warnings and
/// errors go to RDKit's log; the molecules `rdinchi` refuses are its
/// errors. Mirrors `inchi.MolToInchiAndAuxInfo(mol)`.
pub fn mol_to_inchi_and_aux_info(mol: &Mol) -> Result<(String, String)> {
    mol_to_inchi_and_aux_info_with_params(mol, &MolToInchiParams::default())
}

/// As [`mol_to_inchi_and_aux_info`]. Mirrors
/// `inchi.MolToInchiAndAuxInfo(mol, options, logLevel,
/// treatWarningAsError)`.
pub fn mol_to_inchi_and_aux_info_with_params(
    mol: &Mol,
    params: &MolToInchiParams,
) -> Result<(String, String)> {
    let result = rdinchi::mol_to_inchi_with_options(mol, &params.options)?;
    log_result(result.retcode, &result.message, params.log_level);
    if params.treat_warning_as_error && result.retcode != 0 {
        return Err(read_write_error(result.message));
    }
    Ok((result.inchi, result.aux_info))
}

/// The molecule's standard InChI, without AuxInfo. Mirrors
/// `inchi.MolToInchi(mol)`.
pub fn mol_to_inchi(mol: &Mol) -> Result<String> {
    mol_to_inchi_with_params(mol, &MolToInchiParams::default())
}

/// As [`mol_to_inchi`]; `/AuxNone` is added to options without it. Mirrors
/// `inchi.MolToInchi(mol, options, logLevel, treatWarningAsError)`.
pub fn mol_to_inchi_with_params(mol: &Mol, params: &MolToInchiParams) -> Result<String> {
    let mut options = params.options.clone();
    if !options.contains("AuxNone") {
        if options.is_empty() {
            options.push_str("/AuxNone");
        } else {
            options.push_str(" /AuxNone");
        }
    }
    let params = MolToInchiParams {
        options,
        ..params.clone()
    };
    Ok(mol_to_inchi_and_aux_info_with_params(mol, &params)?.0)
}

/// The InChI key of an InChI; `None` when InChI rejects the string.
/// Mirrors `inchi.InchiToInchiKey`.
pub fn inchi_to_inchi_key(inchi: &str) -> Option<String> {
    let key = rdinchi::inchi_to_inchi_key(inchi);
    (!key.is_empty()).then_some(key)
}

/// The molecule's standard InChI key; empty when InChI fails. Mirrors
/// `inchi.MolToInchiKey(mol)`.
pub fn mol_to_inchi_key(mol: &Mol) -> Result<String> {
    mol_to_inchi_key_with_options(mol, "")
}

/// As [`mol_to_inchi_key`], with InChI's options. Mirrors
/// `inchi.MolToInchiKey(mol, options)`.
pub fn mol_to_inchi_key_with_options(mol: &Mol, options: &str) -> Result<String> {
    rdinchi::mol_to_inchi_key_with_options(mol, options)
}
