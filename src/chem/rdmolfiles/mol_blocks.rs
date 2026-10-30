// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! V2000 and V3000 mol blocks, SCSR mol blocks, and SD property lists.

use std::path::Path;

use crate::Result;
use crate::chem::{Mol, MolWriterParams};
use crate::error::to_string;
use crate::ffi::molfiles::ffi;
use crate::util::path_bytes;

/// A V2000 mol block. A molecule with more than 999 atoms, bonds, or
/// substance groups is `Error::Value`. Mirrors `rdmolfiles.MolToV2KMolBlock`.
pub fn mol_to_v2k_mol_block(mol: &Mol) -> Result<String> {
    mol_to_v2k_mol_block_with_params(mol, &MolWriterParams::default())
}

/// Mirrors `rdmolfiles.MolToV2KMolBlock(mol, params, confId)`.
pub fn mol_to_v2k_mol_block_with_params(mol: &Mol, params: &MolWriterParams) -> Result<String> {
    to_string(ffi::mol_to_v2k_mol_block(mol.as_ffi(), &params.to_ffi())?)
}

/// A V3000 mol block. Mirrors `rdmolfiles.MolToV3KMolBlock`.
pub fn mol_to_v3k_mol_block(mol: &Mol) -> Result<String> {
    mol_to_v3k_mol_block_with_params(mol, &MolWriterParams::default())
}

/// `force_v3000` is ignored. Mirrors `rdmolfiles.MolToV3KMolBlock(mol, params,
/// confId)` and its keyword form.
pub fn mol_to_v3k_mol_block_with_params(mol: &Mol, params: &MolWriterParams) -> Result<String> {
    to_string(ffi::mol_to_v3k_mol_block(mol.as_ffi(), &params.to_ffi())?)
}

/// Writes a V3000 mol file. Mirrors `rdmolfiles.MolToV3KMolFile`.
pub fn mol_to_v3k_mol_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<()> {
    mol_to_v3k_mol_file_with_params(mol, filename, &MolWriterParams::default())
}

/// `force_v3000` is ignored. Mirrors `rdmolfiles.MolToV3KMolFile(mol,
/// filename, params, confId)` and its keyword form.
pub fn mol_to_v3k_mol_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolWriterParams,
) -> Result<()> {
    Ok(ffi::mol_to_v3k_mol_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        &params.to_ffi(),
    )?)
}

/// Parameters for the `create_*_property_list_with_params` functions, with
/// Python's defaults (A2). They all take the same optional arguments.
#[derive(Clone, Debug, PartialEq)]
pub struct CreatePropertyListParams {
    /// Written for atoms without the property; empty for `n/a`
    /// without a `[marker]` prefix.
    pub missing_value_marker: String,
    /// Wraps the list before lines this long. With no marker, it must exceed
    /// the first value's length plus one (`Error::Value`).
    pub line_size: u32,
}

impl Default for CreatePropertyListParams {
    fn default() -> CreatePropertyListParams {
        CreatePropertyListParams {
            missing_value_marker: String::new(),
            line_size: 190,
        }
    }
}

const INT: u8 = 0;
const DOUBLE: u8 = 1;
const BOOL: u8 = 2;
const STRING: u8 = 3;

fn create_property_list(
    mol: &mut Mol,
    kind: u8,
    prop_name: &str,
    params: &CreatePropertyListParams,
) -> Result<()> {
    Ok(ffi::create_property_list(
        mol.pin_mut(),
        kind,
        prop_name,
        &params.missing_value_marker,
        params.line_size,
    )?)
}

macro_rules! property_lists {
    ($($name:ident, $with:ident, $kind:expr, $python:literal;)*) => {
        $(
            #[doc = concat!(
                "Stores the atoms' values of a property as one molecule \
                 property that SD writers write as a property list. A value of \
                 another type is `Error::Runtime`. Mirrors `rdmolfiles.",
                $python,
                "`."
            )]
            pub fn $name(mol: &mut Mol, prop_name: &str) -> Result<()> {
                $with(mol, prop_name, &CreatePropertyListParams::default())
            }

            #[doc = concat!(
                "Mirrors `rdmolfiles.",
                $python,
                "(mol, propName, missingValueMarker, lineSize)`."
            )]
            pub fn $with(
                mol: &mut Mol,
                prop_name: &str,
                params: &CreatePropertyListParams,
            ) -> Result<()> {
                create_property_list(mol, $kind, prop_name, params)
            }
        )*
    };
}

property_lists! {
    create_atom_int_property_list, create_atom_int_property_list_with_params, INT,
        "CreateAtomIntPropertyList";
    create_atom_double_property_list, create_atom_double_property_list_with_params, DOUBLE,
        "CreateAtomDoublePropertyList";
    create_atom_bool_property_list, create_atom_bool_property_list_with_params, BOOL,
        "CreateAtomBoolPropertyList";
    create_atom_string_property_list, create_atom_string_property_list_with_params, STRING,
        "CreateAtomStringPropertyList";
}
