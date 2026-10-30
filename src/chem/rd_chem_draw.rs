// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdChemDraw`: writing ChemDraw documents.
//!
//! RDKit reads CDXML and CDX through the ChemDraw SDK, which this crate does
//! not bind (see `MolsFromCDXML` in `inventory/manifest.toml`), so only the
//! writer is here.

use crate::Result;
use crate::chem::Mol;
use crate::error::Error;
use crate::ffi::molfiles::ffi;

/// A ChemDraw format. Mirrors `rdChemDraw.CDXFormat`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum CDXFormat {
    /// Binary CDX. Mirrors `CDXFormat.CDX`.
    Cdx,
    /// CDXML text. Mirrors `CDXFormat.CDXML`.
    Cdxml,
}

/// A CDXML document for a molecule. RDKit writes coordinates as 32-bit fixed
/// point, scaling 2D ones to ChemDraw's bond length, so coordinates that do
/// not fit, or a 2D molecule without bonds of nonzero length, are
/// `Error::Value`. Mirrors `rdChemDraw.MolToChemDrawBlock`.
pub fn mol_to_chem_draw_block(mol: &Mol) -> Result<String> {
    let bytes = mol_to_chem_draw_block_format(mol, CDXFormat::Cdxml)?;
    String::from_utf8(bytes).map_err(|err| Error::Value {
        message: format!("invalid UTF-8 in a string from RDKit: {err}"),
    })
}

/// A ChemDraw document in `format`, as bytes: binary for `Cdx`, text for
/// `Cdxml`. Mirrors `rdChemDraw.MolToChemDrawBlock(mol, format)`, which
/// returns `str` for either and so cannot decode a binary CDX document.
pub fn mol_to_chem_draw_block_format(mol: &Mol, format: CDXFormat) -> Result<Vec<u8>> {
    Ok(ffi::mol_to_chem_draw_block(
        mol.as_ffi(),
        format == CDXFormat::Cdx,
    )?)
}
