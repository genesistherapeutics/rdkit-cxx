// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Molecules and metadata in PNG text chunks.
//!
//! RDKit walks a PNG's chunks without checking their lengths, so the wrappers
//! check that the bytes are whole chunks up to `IEND`, and that each text
//! chunk holds its keyword, before RDKit reads them; anything else is
//! `Error::Parse`. The file forms read the file once and check those bytes.

use std::collections::BTreeMap;
use std::path::Path;

use cxx::UniquePtr;

use crate::chem::{
    Mol, SmilesParserParams,
};
use crate::error::expect_ok;
use crate::ffi::molops::ffi as molops;
use crate::ffi::png::ffi;
use crate::util::path_bytes;
use crate::{Error, Result};

/// What `mol_metadata_to_png_string_with_params` records, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolMetadataToPNGParams {
    /// Records the molecule's pickle, with the process-wide pickle
    /// properties (`get_default_pickle_properties`); `Error::Value` for a
    /// molecule a pickle cannot hold (`Mol::to_binary`).
    pub include_pkl: bool,
    /// Records its CXSMILES.
    pub include_smiles: bool,
    /// Records its mol block.
    pub include_mol: bool,
}

impl Default for MolMetadataToPNGParams {
    fn default() -> MolMetadataToPNGParams {
        MolMetadataToPNGParams {
            include_pkl: true,
            include_smiles: true,
            include_mol: false,
        }
    }
}

impl MolMetadataToPNGParams {
    fn to_ffi(&self) -> ffi::PNGMetadataParamsData {
        ffi::PNGMetadataParamsData {
            include_pkl: self.include_pkl,
            include_smiles: self.include_smiles,
            include_mol: self.include_mol,
        }
    }
}

/// The molecule recorded in a PNG: from the first pickle, SMILES, or mol
/// block in its metadata. Mirrors `rdmolfiles.MolFromPNGString`, which returns
/// `None` where this returns an error.
///
/// # Safety
///
/// A molecule pickle in `png`'s metadata, which RDKit reads first, must be
/// one RDKit wrote, unmodified, as for `Mol::from_pkl_string`.
pub unsafe fn mol_from_png_string(png: &[u8]) -> Result<Mol> {
    // SAFETY: the caller's contract is `mol_from_png_string_with_params`'s.
    unsafe { mol_from_png_string_with_params(png, &SmilesParserParams::default()) }
}

/// `params` parse a recorded SMILES, and their `sanitize` and `remove_hs` a
/// mol block. Mirrors `rdmolfiles.MolFromPNGString(png, params)`.
///
/// # Safety
///
/// As for `mol_from_png_string`.
pub unsafe fn mol_from_png_string_with_params(
    png: &[u8],
    params: &SmilesParserParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_png_string(
        png,
        &params.to_ffi(),
    )?))
}

/// As `mol_from_png_string`, for a PNG file; one that cannot be opened is
/// `Error::Os`. Mirrors `rdmolfiles.MolFromPNGFile`.
///
/// # Safety
///
/// As for `mol_from_png_string`, of the file's data.
pub unsafe fn mol_from_png_file(filename: impl AsRef<Path>) -> Result<Mol> {
    // SAFETY: the caller's contract is `mol_from_png_file_with_params`'s.
    unsafe { mol_from_png_file_with_params(filename, &SmilesParserParams::default()) }
}

/// Mirrors `rdmolfiles.MolFromPNGFile(filename, params)`.
///
/// # Safety
///
/// As for `mol_from_png_string`, of the file's data.
pub unsafe fn mol_from_png_file_with_params(
    filename: impl AsRef<Path>,
    params: &SmilesParserParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_png_file(
        path_bytes(filename.as_ref()),
        &params.to_ffi(),
    )?))
}

/// The tag RDKit gives pickles in PNG metadata. Mirrors `PNGData::pklTag`,
/// the default `tag` of `MolsFromPNGString`.
const PKL_TAG: &str = "rdkitPKL";

/// Parameters for `mols_from_png_string_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolsFromPNGStringParams {
    /// The prefix of the metadata keys to read: pickles (`rdkitPKL`), SMILES
    /// (`SMILES`), or mol blocks (`MOL`).
    pub tag: String,
    /// As `mol_from_png_string_with_params`'s.
    pub params: SmilesParserParams,
}

impl Default for MolsFromPNGStringParams {
    fn default() -> MolsFromPNGStringParams {
        MolsFromPNGStringParams {
            tag: PKL_TAG.into(),
            params: SmilesParserParams::default(),
        }
    }
}

/// Parameters for `mols_from_png_file_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolsFromPNGFileParams {
    /// As `MolsFromPNGStringParams::tag`.
    pub tag: String,
    /// As `mol_from_png_string_with_params`'s.
    pub params: SmilesParserParams,
}

impl Default for MolsFromPNGFileParams {
    fn default() -> MolsFromPNGFileParams {
        MolsFromPNGFileParams {
            tag: PKL_TAG.into(),
            params: SmilesParserParams::default(),
        }
    }
}

fn into_options(mut list: UniquePtr<molops::MolList>) -> Vec<Option<Mol>> {
    let size = expect_ok(molops::mol_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    (0..size)
        .map(|i| {
            let mol = expect_ok(molops::mol_list_take(list.pin_mut(), i));
            (!mol.is_null()).then(|| Mol::from_ffi(mol))
        })
        .collect()
}

/// A molecule for each pickle in a PNG's metadata, in file order; `None`
/// where Python's tuple holds `None` (a SMILES RDKit cannot parse). Mirrors
/// `rdmolfiles.MolsFromPNGString`.
///
/// # Safety
///
/// Each pickle in `png`'s metadata must be one RDKit wrote, unmodified, as
/// for `Mol::from_pkl_string`.
pub unsafe fn mols_from_png_string(png: &[u8]) -> Result<Vec<Option<Mol>>> {
    // SAFETY: the caller's contract is `mols_from_png_string_with_params`'s.
    unsafe { mols_from_png_string_with_params(png, &MolsFromPNGStringParams::default()) }
}

/// Mirrors `rdmolfiles.MolsFromPNGString(png, tag, params)`.
///
/// # Safety
///
/// When `params.tag` selects pickles, each must be one RDKit wrote,
/// unmodified, as for `Mol::from_pkl_string`.
pub unsafe fn mols_from_png_string_with_params(
    png: &[u8],
    params: &MolsFromPNGStringParams,
) -> Result<Vec<Option<Mol>>> {
    Ok(into_options(ffi::mols_from_png_string(
        png,
        &params.tag,
        &params.params.to_ffi(),
    )?))
}

/// As `mols_from_png_string`, for a PNG file. Python returns an empty tuple
/// where this returns an error. Mirrors `rdmolfiles.MolsFromPNGFile`.
///
/// # Safety
///
/// As for `mols_from_png_string`, of the file's data.
pub unsafe fn mols_from_png_file(filename: impl AsRef<Path>) -> Result<Vec<Option<Mol>>> {
    // SAFETY: the caller's contract is `mols_from_png_file_with_params`'s.
    unsafe { mols_from_png_file_with_params(filename, &MolsFromPNGFileParams::default()) }
}

/// Mirrors `rdmolfiles.MolsFromPNGFile(filename, tag, params)`.
///
/// # Safety
///
/// As for `mols_from_png_string_with_params`, of the file's data.
pub unsafe fn mols_from_png_file_with_params(
    filename: impl AsRef<Path>,
    params: &MolsFromPNGFileParams,
) -> Result<Vec<Option<Mol>>> {
    Ok(into_options(ffi::mols_from_png_file(
        path_bytes(filename.as_ref()),
        &params.tag,
        &params.params.to_ffi(),
    )?))
}

/// The PNG `png` with a molecule's pickle and CXSMILES added to its metadata.
/// Mirrors `rdmolfiles.MolMetadataToPNGString`.
pub fn mol_metadata_to_png_string(mol: &Mol, png: &[u8]) -> Result<Vec<u8>> {
    mol_metadata_to_png_string_with_params(mol, png, &MolMetadataToPNGParams::default())
}

/// Mirrors `rdmolfiles.MolMetadataToPNGString(mol, png, includePkl,
/// includeSmiles, includeMol)`.
pub fn mol_metadata_to_png_string_with_params(
    mol: &Mol,
    png: &[u8],
    params: &MolMetadataToPNGParams,
) -> Result<Vec<u8>> {
    Ok(ffi::mol_metadata_to_png_string(
        mol.as_ffi(),
        png,
        &params.to_ffi(),
    )?)
}

/// The PNG in a file with a molecule's metadata added; the file is not
/// changed. Mirrors `rdmolfiles.MolMetadataToPNGFile`.
pub fn mol_metadata_to_png_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<Vec<u8>> {
    mol_metadata_to_png_file_with_params(mol, filename, &MolMetadataToPNGParams::default())
}

/// Mirrors `rdmolfiles.MolMetadataToPNGFile(mol, filename, includePkl,
/// includeSmiles, includeMol)`.
pub fn mol_metadata_to_png_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolMetadataToPNGParams,
) -> Result<Vec<u8>> {
    Ok(ffi::mol_metadata_to_png_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        &params.to_ffi(),
    )?)
}

fn parallel(metadata: &BTreeMap<String, String>) -> (Vec<String>, Vec<String>) {
    (
        metadata.keys().cloned().collect(),
        metadata.values().cloned().collect(),
    )
}

/// The PNG `png` with compressed text chunks added, one per entry, in key
/// order. Mirrors `rdmolfiles.AddMetadataToPNGString`.
pub fn add_metadata_to_png_string(
    metadata: &BTreeMap<String, String>,
    png: &[u8],
) -> Result<Vec<u8>> {
    let (keys, values) = parallel(metadata);
    Ok(ffi::add_metadata_to_png_string(&keys, &values, png)?)
}

/// The PNG in a file with metadata added; the file is not changed. Mirrors
/// `rdmolfiles.AddMetadataToPNGFile`.
pub fn add_metadata_to_png_file(
    metadata: &BTreeMap<String, String>,
    filename: impl AsRef<Path>,
) -> Result<Vec<u8>> {
    let (keys, values) = parallel(metadata);
    Ok(ffi::add_metadata_to_png_file(
        &keys,
        &values,
        path_bytes(filename.as_ref()),
    )?)
}

fn into_list(metadata: UniquePtr<ffi::PNGMetadata>) -> Result<Vec<(String, Vec<u8>)>> {
    let metadata = metadata.as_ref().expect("RDKit returned null metadata");
    let size = expect_ok(ffi::png_metadata_size(metadata));
    (0..size)
        .map(|i| {
            let key = expect_ok(ffi::png_metadata_key(metadata, i))
                .to_str()
                .map_err(|err| Error::Value {
                    message: format!("invalid UTF-8 in a PNG metadata key: {err}"),
                })?
                .to_owned();
            Ok((key, expect_ok(ffi::png_metadata_value(metadata, i))))
        })
        .collect()
}

/// A PNG's text metadata, keyed by name; a later value replaces an earlier
/// one under the same key, as in Python's dict. A key that is not UTF-8 is
/// `Error::Value`. Mirrors `rdmolfiles.MetadataFromPNGString`.
pub fn metadata_from_png_string(png: &[u8]) -> Result<BTreeMap<String, Vec<u8>>> {
    Ok(metadata_from_png_string_as_list(png)?.into_iter().collect())
}

/// Every metadata entry, in file order. Mirrors
/// `rdmolfiles.MetadataFromPNGString(png, asList=True)`.
pub fn metadata_from_png_string_as_list(png: &[u8]) -> Result<Vec<(String, Vec<u8>)>> {
    into_list(ffi::metadata_from_png_string(png)?)
}

/// As `metadata_from_png_string`, for a PNG file. Mirrors
/// `rdmolfiles.MetadataFromPNGFile`.
pub fn metadata_from_png_file(filename: impl AsRef<Path>) -> Result<BTreeMap<String, Vec<u8>>> {
    Ok(metadata_from_png_file_as_list(filename)?
        .into_iter()
        .collect())
}

/// Mirrors `rdmolfiles.MetadataFromPNGFile(filename, asList=True)`.
pub fn metadata_from_png_file_as_list(
    filename: impl AsRef<Path>,
) -> Result<Vec<(String, Vec<u8>)>> {
    into_list(ffi::metadata_from_png_file(path_bytes(filename.as_ref()))?)
}
