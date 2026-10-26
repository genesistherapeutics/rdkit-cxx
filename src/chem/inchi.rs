// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.inchi`, a pure-Python module over `rdinchi` that
//! `rdkit.Chem` re-exports, as `rdkit_cxx::chem` re-exports this one: InChI
//! conversions that log InChI's warnings and errors to RDKit's log, as
//! Python's `RDLogger` does, and can treat them as errors.
//!
//! Python takes a `logging` level for InChI's messages on success; here it
//! is an [`InchiLogLevel`]. Python's `InchiReadWriteError` is
//! `Error::Value` with InChI's message: Python's exception also carries the
//! molecule or InChI made despite the warning, which `rdinchi` returns.

use crate::chem::rdinchi::{self, InchiToMolParams};
use crate::chem::{Conformer, Mol, renumber_atoms};
use crate::geometry::rd_geometry::Point3D;
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

/// Parameters for `mol_from_inchi_with_params` and
/// `mol_from_inchi_and_aux_info_with_params`, with Python's defaults (A2).
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

/// The coordinates of an AuxInfo's `/rC:` layer, in the input's atom order,
/// and whether any is 3D; `None` when there are none, all are zero, or one
/// does not parse.
fn parse_auxinfo_coordinates(auxinfo: &str) -> Option<(Vec<Point3D>, bool)> {
    let layer = layer(auxinfo, "/rC:")?;
    let mut entries: Vec<&str> = layer.split(';').collect();
    while entries.last().is_some_and(|e| e.trim().is_empty()) {
        entries.pop();
    }
    let mut coords = Vec::new();
    for entry in entries {
        let entry = entry.trim();
        if entry.is_empty() {
            return None;
        }
        let parts: Vec<&str> = entry.split(',').collect();
        let value = |i: usize| -> Option<f64> {
            match parts.get(i) {
                Some(p) if !p.is_empty() => p.trim().parse().ok(),
                _ => Some(0.0),
            }
        };
        coords.push(Point3D::new(value(0)?, value(1)?, value(2)?));
    }
    if coords.is_empty()
        || coords
            .iter()
            .all(|p| p.x == 0.0 && p.y == 0.0 && p.z == 0.0)
    {
        return None;
    }
    let is_3d = coords.iter().any(|p| p.z != 0.0);
    Some((coords, is_3d))
}

/// The text of the first `tag` layer (`re.search(tag + '([^/]+)')`).
fn layer<'a>(auxinfo: &'a str, tag: &str) -> Option<&'a str> {
    let mut rest = auxinfo;
    while let Some(i) = rest.find(tag) {
        let after = &rest[i + tag.len()..];
        let end = after.find('/').unwrap_or(after.len());
        if end > 0 {
            return Some(&after[..end]);
        }
        rest = &rest[i + 1..];
    }
    None
}

/// The AuxInfo's `/N:` layer: for each InChI atom, the input's atom index
/// (from 0, so a 0 in the layer is -1).
fn parse_auxinfo_atom_order(auxinfo: &str) -> Option<Vec<i64>> {
    layer(auxinfo, "/N:")?
        .replace(';', ",")
        .split(',')
        .map(|t| t.trim().parse::<i64>().ok().map(|i| i - 1))
        .collect()
}

fn attach_conformer(mut mol: Mol, coords: &Option<(Vec<Point3D>, bool)>) -> Result<Mol> {
    if let Some((coords, is_3d)) = coords
        && coords.len() == mol.get_num_atoms() as usize
    {
        let mut conf = Conformer::new(mol.get_num_atoms());
        conf.set_3d(*is_3d);
        for (i, p) in coords.iter().enumerate() {
            conf.set_atom_position(i as u32, *p)?;
        }
        mol.add_conformer_with_assign_id(&conf, true)?;
    }
    Ok(mol)
}

/// The inverse of `atom_order`, for `renumber_atoms`; `None` for an index
/// past `size`. A negative index counts from the end, as a Python list
/// index does, and one before the start is `Error::Index`, as Python's
/// `IndexError`.
fn build_inverse_permutation(atom_order: &[i64], size: usize) -> Result<Option<Vec<u32>>> {
    let mut new_order = vec![0u32; size];
    for (inchi_idx, &orig_idx) in atom_order.iter().enumerate() {
        if orig_idx >= size as i64 {
            return Ok(None);
        }
        let pos = if orig_idx < 0 {
            orig_idx + size as i64
        } else {
            orig_idx
        };
        if pos < 0 {
            return Err(Error::Index {
                index: orig_idx.unsigned_abs() as u32,
                len: size as u32,
            });
        }
        new_order[pos as usize] = inchi_idx as u32;
    }
    Ok(Some(new_order))
}

/// The molecule an InChI describes, with its atoms in the order, and with
/// the coordinates, of the AuxInfo the InChI came with; as
/// `mol_from_inchi` where the AuxInfo gives no order that fits. Mirrors
/// `inchi.MolFromInchiAndAuxInfo(inchi, auxinfo)`.
pub fn mol_from_inchi_and_aux_info(inchi: &str, auxinfo: &str) -> Result<Option<Mol>> {
    mol_from_inchi_and_aux_info_with_params(inchi, auxinfo, &MolFromInchiParams::default())
}

/// As [`mol_from_inchi_and_aux_info`]. Mirrors
/// `inchi.MolFromInchiAndAuxInfo(inchi, auxinfo, sanitize, removeHs,
/// logLevel, treatWarningAsError)`.
pub fn mol_from_inchi_and_aux_info_with_params(
    inchi: &str,
    auxinfo: &str,
    params: &MolFromInchiParams,
) -> Result<Option<Mol>> {
    let Some(mol) = mol_from_inchi_with_params(inchi, params)? else {
        return Ok(None);
    };
    let coords = parse_auxinfo_coordinates(auxinfo);
    let Some(atom_order) = parse_auxinfo_atom_order(auxinfo) else {
        return attach_conformer(mol, &coords).map(Some);
    };
    let num_mol_atoms = mol.get_num_atoms() as usize;
    if params.remove_hs || atom_order.len() == num_mol_atoms {
        if atom_order.len() != num_mol_atoms {
            return Ok(Some(mol));
        }
        let Some(new_order) = build_inverse_permutation(&atom_order, num_mol_atoms)? else {
            return Ok(Some(mol));
        };
        return attach_conformer(renumber_atoms(&mol, &new_order)?, &coords).map(Some);
    }
    if atom_order.len() > num_mol_atoms {
        return Ok(Some(mol));
    }
    // The layer numbers the heavy atoms; the hydrogens follow them.
    let num_heavy = atom_order.len();
    let Some(mut new_order) = build_inverse_permutation(&atom_order, num_heavy)? else {
        return Ok(Some(mol));
    };
    new_order.resize(num_mol_atoms, 0);
    let mut h_slot = num_heavy;
    for i in 0..num_mol_atoms {
        if !atom_order.contains(&(i as i64)) {
            if h_slot >= num_mol_atoms {
                return Ok(Some(mol));
            }
            new_order[h_slot] = i as u32;
            h_slot += 1;
        }
    }
    attach_conformer(renumber_atoms(&mol, &new_order)?, &coords).map(Some)
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
