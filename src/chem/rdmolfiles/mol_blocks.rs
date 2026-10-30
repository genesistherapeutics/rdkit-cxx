// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! V2000 and V3000 mol blocks, SCSR mol blocks, and SD property lists.

use std::path::Path;

use crate::Result;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::{Mol, MolWriterParams};
use crate::error::{expect_ok, to_string};
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

rdkit_enum! {
    /// How SCSR templates are named in the expanded molecule. Mirrors
    /// `rdmolfiles.SCSRTemplateNames`.
    SCSRTemplateNames: ffi::SCSRTemplateNames {
        UseFirstName = UseFirstName,
        UseSecondName = UseSecondName,
        AsEntered = AsEntered,
    }
}

rdkit_enum! {
    /// Which hydrogen bonds between SCSR bases to keep. Mirrors
    /// `rdmolfiles.SCSRBaseHbondOptions`.
    SCSRBaseHbondOptions: ffi::SCSRBaseHbondOptions {
        Ignore = Ignore,
        UseSapAll = UseSapAll,
        UseSapOne = UseSapOne,
        Auto = Auto,
    }
}

/// How an SCSR mol block becomes a molecule. Mirrors
/// `rdmolfiles.MolFromSCSRParams`, plus `MolFromSCSRBlock`'s `sanitize` and
/// `removeHs` (CONVENTIONS.md A3).
#[derive(Clone, Debug, PartialEq)]
pub struct MolFromSCSRParams {
    /// Keeps leaving groups on atoms that are not bonded out of the template.
    pub include_leaving_groups: bool,
    /// Which template names to use.
    pub scsr_template_names: SCSRTemplateNames,
    /// How hydrogen bonds between bases are kept.
    pub scsr_base_hbond_options: SCSRBaseHbondOptions,
    /// Sanitizes the molecule.
    pub sanitize: bool,
    /// Removes hydrogens (only with `sanitize`).
    pub remove_hs: bool,
}

impl Default for MolFromSCSRParams {
    fn default() -> MolFromSCSRParams {
        let data = expect_ok(ffi::mol_from_scsr_params_default());
        MolFromSCSRParams {
            include_leaving_groups: data.include_leaving_groups,
            scsr_template_names: SCSRTemplateNames::from_ffi(data.scsr_template_names),
            scsr_base_hbond_options: SCSRBaseHbondOptions::from_ffi(data.scsr_base_hbond_options),
            sanitize: true,
            remove_hs: true,
        }
    }
}

impl MolFromSCSRParams {
    fn to_ffi(&self) -> ffi::MolFromSCSRParamsData {
        ffi::MolFromSCSRParamsData {
            include_leaving_groups: self.include_leaving_groups,
            scsr_template_names: self.scsr_template_names.to_ffi(),
            scsr_base_hbond_options: self.scsr_base_hbond_options.to_ffi(),
        }
    }
}

/// Expands an SCSR mol block, whose templates RDKit parses leniently, into a
/// full molecule. Mirrors `rdmolfiles.MolFromSCSRBlock`, which returns `None`
/// where this returns an error.
pub fn mol_from_scsr_block(mol_block: &str) -> Result<Mol> {
    mol_from_scsr_block_with_params(mol_block, &MolFromSCSRParams::default())
}

/// Mirrors `rdmolfiles.MolFromSCSRBlock(molBlock, sanitize, removeHs,
/// molFromSCSRParams)`.
pub fn mol_from_scsr_block_with_params(mol_block: &str, params: &MolFromSCSRParams) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_scsr_block(
        mol_block,
        params.sanitize,
        params.remove_hs,
        &params.to_ffi(),
    )?))
}

/// Reads an SCSR mol file. Mirrors `rdmolfiles.MolFromSCSRFile`.
pub fn mol_from_scsr_file(filename: impl AsRef<Path>) -> Result<Mol> {
    mol_from_scsr_file_with_params(filename, &MolFromSCSRParams::default())
}

/// Mirrors `rdmolfiles.MolFromSCSRFile(filename, sanitize, removeHs,
/// molFromSCSRParams)`.
pub fn mol_from_scsr_file_with_params(
    filename: impl AsRef<Path>,
    params: &MolFromSCSRParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_scsr_file(
        path_bytes(filename.as_ref()),
        params.sanitize,
        params.remove_hs,
        &params.to_ffi(),
    )?))
}

/// Parameters for the `create_*_property_list_with_params` functions, with
/// Python's defaults (A2). They all take the same optional arguments.
#[derive(Clone, Debug, PartialEq)]
pub struct CreatePropertyListParams {
    /// Written for atoms or bonds without the property; empty for `n/a`
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
