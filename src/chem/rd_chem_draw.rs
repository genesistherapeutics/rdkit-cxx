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
use crate::chem::rdmolfiles::{CDXMLFormat, mol_to_cdxml_block, mol_to_cdxml_block_format};

/// A ChemDraw format. Mirrors `rdChemDraw.CDXFormat`.
///
/// RDKit's `MolToCDXMLBlock` writes through `MolToChemDrawBlock`, mapping
/// `rdmolfiles.CDXMLFormat` onto this enum, so the writer here goes the other
/// way, through `rdmolfiles`' checked writer.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum CDXFormat {
    /// Binary CDX. Mirrors `CDXFormat.CDX`.
    Cdx,
    /// CDXML text. Mirrors `CDXFormat.CDXML`.
    Cdxml,
}

impl From<CDXFormat> for CDXMLFormat {
    fn from(format: CDXFormat) -> CDXMLFormat {
        match format {
            CDXFormat::Cdx => CDXMLFormat::Cdx,
            CDXFormat::Cdxml => CDXMLFormat::Cdxml,
        }
    }
}

/// A CDXML document for a molecule, as `rdmolfiles::mol_to_cdxml_block`
/// writes it: coordinates that do not fit ChemDraw's 32-bit fixed point, or
/// a 2D molecule without bonds of nonzero length, are `Error::Value`.
/// Mirrors `rdChemDraw.MolToChemDrawBlock`.
pub fn mol_to_chem_draw_block(mol: &Mol) -> Result<String> {
    mol_to_cdxml_block(mol)
}

/// A ChemDraw document in `format`, as bytes: binary for `Cdx`, text for
/// `Cdxml`. Mirrors `rdChemDraw.MolToChemDrawBlock(mol, format)`, which
/// returns `str` for either and so cannot decode a binary CDX document.
pub fn mol_to_chem_draw_block_format(mol: &Mol, format: CDXFormat) -> Result<Vec<u8>> {
    mol_to_cdxml_block_format(mol, format.into())
}
