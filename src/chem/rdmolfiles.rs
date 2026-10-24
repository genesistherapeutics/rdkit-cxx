// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdmolfiles`.

use std::collections::BTreeMap;
use std::path::Path;

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::molfiles::ffi;

/// Parameters for SMILES parsing. Mirrors `rdmolfiles.SmilesParserParams`,
/// plus `MolFromSmiles`'s `replacements` keyword (CONVENTIONS.md A3).
#[derive(Clone, Debug, PartialEq)]
pub struct SmilesParserParams {
    /// The amount of parser debugging output.
    pub debug_parse: i32,
    /// Parses the molecule's name after the SMILES.
    pub parse_name: bool,
    /// Parses CXSMILES extensions.
    pub allow_cxsmiles: bool,
    /// Fails the parse when CXSMILES parsing fails.
    pub strict_cxsmiles: bool,
    /// Sanitizes the molecule.
    pub sanitize: bool,
    /// Removes hydrogens after parsing.
    pub remove_hs: bool,
    /// Text substitutions applied to the SMILES before parsing, repeated until
    /// it stops changing.
    pub replacements: BTreeMap<String, String>,
}

impl Default for SmilesParserParams {
    fn default() -> SmilesParserParams {
        let data = expect_ok(ffi::smiles_parser_params_default());
        SmilesParserParams {
            debug_parse: data.debug_parse,
            parse_name: data.parse_name,
            allow_cxsmiles: data.allow_cxsmiles,
            strict_cxsmiles: data.strict_cxsmiles,
            sanitize: data.sanitize,
            remove_hs: data.remove_hs,
            replacements: BTreeMap::new(),
        }
    }
}

impl SmilesParserParams {
    fn to_ffi(&self) -> ffi::SmilesParserParamsData {
        ffi::SmilesParserParamsData {
            debug_parse: self.debug_parse,
            parse_name: self.parse_name,
            allow_cxsmiles: self.allow_cxsmiles,
            strict_cxsmiles: self.strict_cxsmiles,
            sanitize: self.sanitize,
            remove_hs: self.remove_hs,
            replacement_keys: self.replacements.keys().cloned().collect(),
            replacement_values: self.replacements.values().cloned().collect(),
        }
    }
}

/// Parameters for SMARTS parsing. Mirrors `rdmolfiles.SmartsParserParams`,
/// plus `MolFromSmarts`'s `replacements` keyword (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct SmartsParserParams {
    /// The amount of parser debugging output.
    pub debug_parse: i32,
    /// Parses the molecule's name after the SMARTS.
    pub parse_name: bool,
    /// Parses CXSMILES extensions.
    pub allow_cxsmiles: bool,
    /// Fails the parse when CXSMILES parsing fails.
    pub strict_cxsmiles: bool,
    /// Merges hydrogen atoms into their neighbors' queries.
    pub merge_hs: bool,
    /// Text substitutions applied to the SMARTS before parsing.
    pub replacements: BTreeMap<String, String>,
}

impl Default for SmartsParserParams {
    fn default() -> SmartsParserParams {
        let data = expect_ok(ffi::smarts_parser_params_default());
        SmartsParserParams {
            debug_parse: data.debug_parse,
            parse_name: data.parse_name,
            allow_cxsmiles: data.allow_cxsmiles,
            strict_cxsmiles: data.strict_cxsmiles,
            merge_hs: data.merge_hs,
            replacements: BTreeMap::new(),
        }
    }
}

impl SmartsParserParams {
    fn to_ffi(&self) -> ffi::SmartsParserParamsData {
        ffi::SmartsParserParamsData {
            debug_parse: self.debug_parse,
            parse_name: self.parse_name,
            allow_cxsmiles: self.allow_cxsmiles,
            strict_cxsmiles: self.strict_cxsmiles,
            merge_hs: self.merge_hs,
            replacement_keys: self.replacements.keys().cloned().collect(),
            replacement_values: self.replacements.values().cloned().collect(),
        }
    }
}

/// Parameters for SMILES and SMARTS writing. Mirrors
/// `rdmolfiles.SmilesWriteParams`.
#[derive(Clone, Debug, PartialEq)]
pub struct SmilesWriteParams {
    /// Includes stereochemistry and isotopes.
    pub do_isomeric_smiles: bool,
    /// Writes the Kekulé form.
    pub do_kekule: bool,
    /// Writes the canonical form.
    pub canonical: bool,
    /// Removes stereo markers that canonicalization finds redundant.
    pub clean_stereo: bool,
    /// Writes every bond symbol.
    pub all_bonds_explicit: bool,
    /// Writes every hydrogen in brackets.
    pub all_hs_explicit: bool,
    /// Randomizes the atom order.
    pub do_random: bool,
    /// The atom to start at, or -1 for any.
    pub rooted_at_atom: i32,
    /// Writes dative bonds.
    pub include_dative_bonds: bool,
    /// Ignores atom map numbers when canonicalizing.
    pub ignore_atom_map_numbers: bool,
}

impl Default for SmilesWriteParams {
    fn default() -> SmilesWriteParams {
        let data = expect_ok(ffi::smiles_write_params_default());
        SmilesWriteParams {
            do_isomeric_smiles: data.do_isomeric_smiles,
            do_kekule: data.do_kekule,
            canonical: data.canonical,
            clean_stereo: data.clean_stereo,
            all_bonds_explicit: data.all_bonds_explicit,
            all_hs_explicit: data.all_hs_explicit,
            do_random: data.do_random,
            rooted_at_atom: data.rooted_at_atom,
            include_dative_bonds: data.include_dative_bonds,
            ignore_atom_map_numbers: data.ignore_atom_map_numbers,
        }
    }
}

impl SmilesWriteParams {
    pub(crate) fn to_ffi(&self) -> ffi::SmilesWriteParamsData {
        ffi::SmilesWriteParamsData {
            do_isomeric_smiles: self.do_isomeric_smiles,
            do_kekule: self.do_kekule,
            canonical: self.canonical,
            clean_stereo: self.clean_stereo,
            all_bonds_explicit: self.all_bonds_explicit,
            all_hs_explicit: self.all_hs_explicit,
            do_random: self.do_random,
            rooted_at_atom: self.rooted_at_atom,
            include_dative_bonds: self.include_dative_bonds,
            ignore_atom_map_numbers: self.ignore_atom_map_numbers,
        }
    }
}

/// Parameters for writing mol blocks and files. Mirrors
/// `rdmolfiles.MolWriterParams`, plus the writers' `confId` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct MolWriterParams {
    /// Includes stereochemistry.
    pub include_stereo: bool,
    /// Kekulizes the molecule before writing.
    pub kekulize: bool,
    /// Writes V3000 even when V2000 would do.
    pub force_v3000: bool,
    /// Coordinate precision (V3000 only).
    pub precision: u32,
    /// The conformer to write, or -1 for the default.
    pub conf_id: i32,
}

impl Default for MolWriterParams {
    fn default() -> MolWriterParams {
        let data = expect_ok(ffi::mol_writer_params_default());
        MolWriterParams {
            include_stereo: data.include_stereo,
            kekulize: data.kekulize,
            force_v3000: data.force_v3000,
            precision: data.precision,
            conf_id: data.conf_id,
        }
    }
}

impl MolWriterParams {
    fn to_ffi(&self) -> ffi::MolWriterParamsData {
        ffi::MolWriterParamsData {
            include_stereo: self.include_stereo,
            kekulize: self.kekulize,
            force_v3000: self.force_v3000,
            precision: self.precision,
            conf_id: self.conf_id,
        }
    }
}

/// Parameters for `mol_from_mol_block_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolFromMolBlockParams {
    /// Sanitizes the molecule.
    pub sanitize: bool,
    /// Removes hydrogens (only with `sanitize`).
    pub remove_hs: bool,
    /// Rejects malformed blocks rather than reading what it can.
    pub strict_parsing: bool,
}

impl Default for MolFromMolBlockParams {
    fn default() -> MolFromMolBlockParams {
        MolFromMolBlockParams {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
        }
    }
}

/// Parameters for `mol_from_mol_file_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolFromMolFileParams {
    /// Sanitizes the molecule.
    pub sanitize: bool,
    /// Removes hydrogens (only with `sanitize`).
    pub remove_hs: bool,
    /// Rejects malformed files rather than reading what it can.
    pub strict_parsing: bool,
}

impl Default for MolFromMolFileParams {
    fn default() -> MolFromMolFileParams {
        MolFromMolFileParams {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
        }
    }
}

fn path_bytes(path: &Path) -> &[u8] {
    path.as_os_str().as_encoded_bytes()
}

/// Parses a SMILES string, sanitizing the result. Mirrors
/// `rdmolfiles.MolFromSmiles`, which returns `None` where this returns
/// `Error::Parse` or a sanitization error (E2).
pub fn mol_from_smiles(smiles: &str) -> Result<Mol> {
    mol_from_smiles_with_params(smiles, &SmilesParserParams::default())
}

/// Mirrors `rdmolfiles.MolFromSmiles(SMILES, params)` and
/// `MolFromSmiles(SMILES, sanitize, replacements)`.
pub fn mol_from_smiles_with_params(smiles: &str, params: &SmilesParserParams) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_smiles(
        smiles,
        &params.to_ffi(),
    )?))
}

/// The canonical SMILES for a molecule. Mirrors `rdmolfiles.MolToSmiles`.
pub fn mol_to_smiles(mol: &Mol) -> Result<String> {
    mol_to_smiles_with_params(mol, &SmilesWriteParams::default())
}

/// Mirrors `rdmolfiles.MolToSmiles(mol, params)` and its keyword form.
pub fn mol_to_smiles_with_params(mol: &Mol, params: &SmilesWriteParams) -> Result<String> {
    to_string(ffi::mol_to_smiles(mol.as_ffi(), &params.to_ffi())?)
}

/// Parses a SMARTS query. Mirrors `rdmolfiles.MolFromSmarts`.
pub fn mol_from_smarts(smarts: &str) -> Result<Mol> {
    mol_from_smarts_with_params(smarts, &SmartsParserParams::default())
}

/// Mirrors `rdmolfiles.MolFromSmarts(SMARTS, params)` and
/// `MolFromSmarts(SMARTS, mergeHs, replacements)`.
pub fn mol_from_smarts_with_params(smarts: &str, params: &SmartsParserParams) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_smarts(
        smarts,
        &params.to_ffi(),
    )?))
}

/// A SMARTS for a molecule. Mirrors `rdmolfiles.MolToSmarts`.
pub fn mol_to_smarts(mol: &Mol) -> Result<String> {
    mol_to_smarts_with_params(mol, &SmilesWriteParams::default())
}

/// Mirrors `rdmolfiles.MolToSmarts(mol, params)` and its keyword form.
pub fn mol_to_smarts_with_params(mol: &Mol, params: &SmilesWriteParams) -> Result<String> {
    to_string(ffi::mol_to_smarts(mol.as_ffi(), &params.to_ffi())?)
}

/// Parses a mol block. Mirrors `rdmolfiles.MolFromMolBlock`.
pub fn mol_from_mol_block(mol_block: &str) -> Result<Mol> {
    mol_from_mol_block_with_params(mol_block, &MolFromMolBlockParams::default())
}

/// Mirrors `rdmolfiles.MolFromMolBlock(molBlock, sanitize, removeHs,
/// strictParsing)`.
pub fn mol_from_mol_block_with_params(
    mol_block: &str,
    params: &MolFromMolBlockParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_mol_block(
        mol_block,
        params.sanitize,
        params.remove_hs,
        params.strict_parsing,
    )?))
}

/// Reads a mol file. A file that cannot be opened is `Error::Os`. Mirrors
/// `rdmolfiles.MolFromMolFile`.
pub fn mol_from_mol_file(mol_file_name: impl AsRef<Path>) -> Result<Mol> {
    mol_from_mol_file_with_params(mol_file_name, &MolFromMolFileParams::default())
}

/// Mirrors `rdmolfiles.MolFromMolFile(molFileName, sanitize, removeHs,
/// strictParsing)`.
pub fn mol_from_mol_file_with_params(
    mol_file_name: impl AsRef<Path>,
    params: &MolFromMolFileParams,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::mol_from_mol_file(
        path_bytes(mol_file_name.as_ref()),
        params.sanitize,
        params.remove_hs,
        params.strict_parsing,
    )?))
}

/// A mol block for a molecule. Mirrors `rdmolfiles.MolToMolBlock`.
pub fn mol_to_mol_block(mol: &Mol) -> Result<String> {
    mol_to_mol_block_with_params(mol, &MolWriterParams::default())
}

/// Mirrors `rdmolfiles.MolToMolBlock(mol, params, confId)` and its keyword
/// form.
pub fn mol_to_mol_block_with_params(mol: &Mol, params: &MolWriterParams) -> Result<String> {
    to_string(ffi::mol_to_mol_block(mol.as_ffi(), &params.to_ffi())?)
}

/// Writes a mol file. Mirrors `rdmolfiles.MolToMolFile`.
pub fn mol_to_mol_file(mol: &Mol, filename: impl AsRef<Path>) -> Result<()> {
    mol_to_mol_file_with_params(mol, filename, &MolWriterParams::default())
}

/// Mirrors `rdmolfiles.MolToMolFile(mol, filename, params, confId)` and its
/// keyword form.
pub fn mol_to_mol_file_with_params(
    mol: &Mol,
    filename: impl AsRef<Path>,
    params: &MolWriterParams,
) -> Result<()> {
    Ok(ffi::mol_to_mol_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        &params.to_ffi(),
    )?)
}
