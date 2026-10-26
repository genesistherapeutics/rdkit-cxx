// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Reading and writing reactions: SMARTS and SMILES, rxn blocks, Marvin,
//! CDXML, molecules with reaction roles, and PNG metadata.

use std::collections::BTreeMap;
use std::path::Path;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rd_chem_reactions::ChemicalReaction;
use crate::chem::{Mol, SmilesWriteParams};
use crate::error::{expect_ok, to_string};
use crate::ffi::chem_reactions::ffi;
use crate::util::path_bytes;

/// Parameters for `reaction_from_smarts_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReactionFromSmartsParams {
    /// Text substitutions applied to each component before parsing.
    pub replacements: BTreeMap<String, String>,
    /// Parses the components as SMILES rather than SMARTS.
    pub use_smiles: bool,
}

/// Parameters for `reaction_to_cx_smarts_with_params`: Python's
/// `SmilesWriteParams` overload and its `flags` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionToCXSmartsParams {
    /// How to write each template.
    pub params: SmilesWriteParams,
    /// Which CXSMILES extensions to write: `CXSmilesFields` bits, all by
    /// default.
    pub flags: u32,
}

impl Default for ReactionToCXSmartsParams {
    fn default() -> ReactionToCXSmartsParams {
        ReactionToCXSmartsParams {
            params: SmilesWriteParams::default(),
            flags: expect_ok(ffi::cx_smiles_fields_all()),
        }
    }
}

/// Parameters for `reaction_to_cx_smiles_with_params`: those of the SMARTS
/// form. A Rust addition (A2).
pub type ReactionToCXSmilesParams = ReactionToCXSmartsParams;

/// Parameters for `reaction_from_rxn_block_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionFromRxnBlockParams {
    /// Sanitizes the templates.
    pub sanitize: bool,
    /// Removes hydrogens (only with `sanitize`).
    pub remove_hs: bool,
    /// Rejects malformed blocks rather than reading what it can.
    pub strict_parsing: bool,
}

impl Default for ReactionFromRxnBlockParams {
    fn default() -> ReactionFromRxnBlockParams {
        ReactionFromRxnBlockParams {
            sanitize: false,
            remove_hs: false,
            strict_parsing: true,
        }
    }
}

/// Parameters for `reaction_from_rxn_file_with_params`: those of the block
/// form. A Rust addition (A2).
pub type ReactionFromRxnFileParams = ReactionFromRxnBlockParams;

/// Parameters for `reaction_from_mrv_block_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReactionFromMrvBlockParams {
    /// Sanitizes the templates.
    pub sanitize: bool,
    /// Removes hydrogens.
    pub remove_hs: bool,
}

/// Parameters for `reaction_from_mrv_file_with_params`: those of the block
/// form. A Rust addition (A2).
pub type ReactionFromMrvFileParams = ReactionFromMrvBlockParams;

/// Parameters for `reactions_from_cdxml_block_with_params`, with Python's
/// defaults (A2).
pub type ReactionsFromCDXMLBlockParams = ReactionFromMrvBlockParams;

/// Parameters for `reactions_from_cdxml_file_with_params`, with Python's
/// defaults (A2).
pub type ReactionsFromCDXMLFileParams = ReactionFromMrvBlockParams;

/// Parameters for `reaction_to_rxn_block_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReactionToRxnBlockParams {
    /// Writes the agents in their own block rather than with the reactants.
    pub separate_agents: bool,
    /// Writes V3000.
    pub force_v3000: bool,
}

/// Parameters for `reaction_metadata_to_png_string_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionMetadataToPNGStringParams {
    /// Adds the reaction's pickle.
    pub include_pkl: bool,
    /// Adds its reaction SMILES.
    pub include_smiles: bool,
    /// Adds its reaction SMARTS.
    pub include_smarts: bool,
    /// Adds its V2000 rxn block, which depicts each template without a
    /// conformer, as `reaction_to_rxn_block` does.
    pub include_rxn: bool,
}

impl Default for ReactionMetadataToPNGStringParams {
    fn default() -> ReactionMetadataToPNGStringParams {
        ReactionMetadataToPNGStringParams {
            include_pkl: true,
            include_smiles: true,
            include_smarts: false,
            include_rxn: false,
        }
    }
}

/// Parameters for `reaction_metadata_to_png_file_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionMetadataToPNGFileParams {
    /// Adds the reaction's pickle.
    pub include_pkl: bool,
    /// Adds its reaction SMILES.
    pub include_smiles: bool,
    /// Adds its reaction SMARTS.
    pub include_smarts: bool,
    /// Adds its V2000 rxn block (Python's `includeMol`).
    pub include_mol: bool,
}

impl Default for ReactionMetadataToPNGFileParams {
    fn default() -> ReactionMetadataToPNGFileParams {
        ReactionMetadataToPNGFileParams {
            include_pkl: true,
            include_smiles: true,
            include_smarts: false,
            include_mol: false,
        }
    }
}

fn reaction(
    ptr: std::result::Result<UniquePtr<ffi::ChemicalReaction>, cxx::Exception>,
) -> Result<ChemicalReaction> {
    Ok(ChemicalReaction::from_ffi(ptr?))
}

fn reactions(
    list: std::result::Result<UniquePtr<ffi::ReactionList>, cxx::Exception>,
) -> Result<Vec<ChemicalReaction>> {
    let mut list = list?;
    let size = expect_ok(ffi::reaction_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    Ok((0..size)
        .map(|i| ChemicalReaction::from_ffi(expect_ok(ffi::reaction_list_take(list.pin_mut(), i))))
        .collect())
}

/// Parses reaction SMARTS: `reactants>agents>products`, with atom maps
/// linking reactant and product atoms. Unparsable SMARTS is `Error::Parse`.
/// Mirrors `rdChemReactions.ReactionFromSmarts`.
pub fn reaction_from_smarts(smarts: &str) -> Result<ChemicalReaction> {
    reaction_from_smarts_with_params(smarts, &ReactionFromSmartsParams::default())
}

/// Mirrors `rdChemReactions.ReactionFromSmarts(SMARTS, replacements,
/// useSmiles)`.
pub fn reaction_from_smarts_with_params(
    smarts: &str,
    params: &ReactionFromSmartsParams,
) -> Result<ChemicalReaction> {
    let keys = params.replacements.keys().cloned().collect();
    let values = params.replacements.values().cloned().collect();
    reaction(ffi::reaction_from_smarts(
        smarts,
        &keys,
        &values,
        params.use_smiles,
    ))
}

/// Parses reaction SMILES. Mirrors `rdChemReactions.ReactionFromSmiles`.
pub fn reaction_from_smiles(smiles: &str) -> Result<ChemicalReaction> {
    reaction_from_smiles_with_replacements(smiles, &BTreeMap::new())
}

/// Mirrors `rdChemReactions.ReactionFromSmiles(SMILES, replacements)`.
pub fn reaction_from_smiles_with_replacements(
    smiles: &str,
    replacements: &BTreeMap<String, String>,
) -> Result<ChemicalReaction> {
    reaction_from_smarts_with_params(
        smiles,
        &ReactionFromSmartsParams {
            replacements: replacements.clone(),
            use_smiles: true,
        },
    )
}

fn non_canonical() -> SmilesWriteParams {
    SmilesWriteParams {
        canonical: false,
        ..SmilesWriteParams::default()
    }
}

/// Reaction SMARTS for a reaction, its templates in order. Mirrors
/// `rdChemReactions.ReactionToSmarts`.
pub fn reaction_to_smarts(reaction: &ChemicalReaction) -> Result<String> {
    reaction_to_smarts_with_params(reaction, &non_canonical())
}

/// Mirrors `rdChemReactions.ReactionToSmarts(reaction, params)`.
pub fn reaction_to_smarts_with_params(
    reaction: &ChemicalReaction,
    params: &SmilesWriteParams,
) -> Result<String> {
    to_string(ffi::reaction_to_smarts(
        reaction.as_ffi(),
        &params.to_ffi(),
    )?)
}

/// Canonical reaction SMILES for a reaction. Mirrors
/// `rdChemReactions.ReactionToSmiles`.
pub fn reaction_to_smiles(reaction: &ChemicalReaction) -> Result<String> {
    reaction_to_smiles_with_params(reaction, &SmilesWriteParams::default())
}

/// Mirrors `rdChemReactions.ReactionToSmiles(reaction, params)` and
/// `ReactionToSmiles(reaction, canonical)`.
pub fn reaction_to_smiles_with_params(
    reaction: &ChemicalReaction,
    params: &SmilesWriteParams,
) -> Result<String> {
    to_string(ffi::reaction_to_smiles(
        reaction.as_ffi(),
        &params.to_ffi(),
    )?)
}

/// Reaction SMARTS with CXSMILES extensions. Mirrors
/// `rdChemReactions.ReactionToCXSmarts`.
pub fn reaction_to_cx_smarts(reaction: &ChemicalReaction) -> Result<String> {
    reaction_to_cx_smarts_with_params(
        reaction,
        &ReactionToCXSmartsParams {
            params: non_canonical(),
            ..Default::default()
        },
    )
}

/// Mirrors `rdChemReactions.ReactionToCXSmarts(reaction, params, flags)`.
pub fn reaction_to_cx_smarts_with_params(
    reaction: &ChemicalReaction,
    params: &ReactionToCXSmartsParams,
) -> Result<String> {
    to_string(ffi::reaction_to_cx_smarts(
        reaction.as_ffi(),
        &params.params.to_ffi(),
        params.flags,
    )?)
}

/// Canonical reaction SMILES with CXSMILES extensions. Mirrors
/// `rdChemReactions.ReactionToCXSmiles`.
pub fn reaction_to_cx_smiles(reaction: &ChemicalReaction) -> Result<String> {
    reaction_to_cx_smiles_with_params(reaction, &ReactionToCXSmilesParams::default())
}

/// Mirrors `rdChemReactions.ReactionToCXSmiles(reaction, params, flags)` and
/// `ReactionToCXSmiles(reaction, canonical)`.
pub fn reaction_to_cx_smiles_with_params(
    reaction: &ChemicalReaction,
    params: &ReactionToCXSmilesParams,
) -> Result<String> {
    to_string(ffi::reaction_to_cx_smiles(
        reaction.as_ffi(),
        &params.params.to_ffi(),
        params.flags,
    )?)
}

/// Parses an MDL rxn block (V2000 or V3000). Its templates are mol blocks,
/// checked as `mol_from_mol_block` checks one: past those limits it is
/// `Error::Value` (CONVENTIONS.md U12). Mirrors
/// `rdChemReactions.ReactionFromRxnBlock`.
pub fn reaction_from_rxn_block(rxnblock: &str) -> Result<ChemicalReaction> {
    reaction_from_rxn_block_with_params(rxnblock, &ReactionFromRxnBlockParams::default())
}

/// Mirrors `rdChemReactions.ReactionFromRxnBlock(rxnblock, sanitize,
/// removeHs, strictParsing)`.
pub fn reaction_from_rxn_block_with_params(
    rxnblock: &str,
    params: &ReactionFromRxnBlockParams,
) -> Result<ChemicalReaction> {
    reaction(ffi::reaction_from_rxn_block(
        rxnblock,
        params.sanitize,
        params.remove_hs,
        params.strict_parsing,
    ))
}

/// Reads an MDL rxn file. A file that cannot be opened is `Error::Os`, and
/// an empty one `Error::Parse`. Mirrors `rdChemReactions.ReactionFromRxnFile`,
/// which returns `None` for both.
pub fn reaction_from_rxn_file(filename: impl AsRef<Path>) -> Result<ChemicalReaction> {
    reaction_from_rxn_file_with_params(filename, &ReactionFromRxnFileParams::default())
}

/// Mirrors `rdChemReactions.ReactionFromRxnFile(filename, sanitize, removeHs,
/// strictParsing)`.
pub fn reaction_from_rxn_file_with_params(
    filename: impl AsRef<Path>,
    params: &ReactionFromRxnFileParams,
) -> Result<ChemicalReaction> {
    reaction(ffi::reaction_from_rxn_file(
        path_bytes(filename.as_ref()),
        params.sanitize,
        params.remove_hs,
        params.strict_parsing,
    ))
}

/// Parses a Marvin (MRV) reaction. Python returns `None` for any failure;
/// this returns RDKit's error: `Error::Parse` for Marvin data that is not a
/// reaction RDKit can read (and for a bond without two atoms, which RDKit
/// reads past), `Error::Runtime` for malformed XML. Mirrors
/// `rdChemReactions.ReactionFromMrvBlock`.
pub fn reaction_from_mrv_block(rxnblock: &str) -> Result<ChemicalReaction> {
    reaction_from_mrv_block_with_params(rxnblock, &ReactionFromMrvBlockParams::default())
}

/// Mirrors `rdChemReactions.ReactionFromMrvBlock(rxnblock, sanitize,
/// removeHs)`.
pub fn reaction_from_mrv_block_with_params(
    rxnblock: &str,
    params: &ReactionFromMrvBlockParams,
) -> Result<ChemicalReaction> {
    reaction(ffi::reaction_from_mrv_block(
        rxnblock,
        params.sanitize,
        params.remove_hs,
    ))
}

/// Reads a Marvin reaction file. A file that cannot be opened is
/// `Error::Os`. Mirrors `rdChemReactions.ReactionFromMrvFile`.
pub fn reaction_from_mrv_file(filename: impl AsRef<Path>) -> Result<ChemicalReaction> {
    reaction_from_mrv_file_with_params(filename, &ReactionFromMrvFileParams::default())
}

/// Mirrors `rdChemReactions.ReactionFromMrvFile(filename, sanitize,
/// removeHs)`.
pub fn reaction_from_mrv_file_with_params(
    filename: impl AsRef<Path>,
    params: &ReactionFromMrvFileParams,
) -> Result<ChemicalReaction> {
    reaction(ffi::reaction_from_mrv_file(
        path_bytes(filename.as_ref()),
        params.sanitize,
        params.remove_hs,
    ))
}

/// Whether Marvin data holds a reaction. Malformed XML is `Error::Runtime`.
/// Mirrors `rdChemReactions.MrvBlockIsReaction`.
pub fn mrv_block_is_reaction(mrv_data: &str) -> Result<bool> {
    Ok(ffi::mrv_block_is_reaction(mrv_data)?)
}

/// Whether a Marvin file holds a reaction. Mirrors
/// `rdChemReactions.MrvFileIsReaction`.
pub fn mrv_file_is_reaction(filename: impl AsRef<Path>) -> Result<bool> {
    Ok(ffi::mrv_file_is_reaction(path_bytes(filename.as_ref()))?)
}

/// The reactions in CDXML data, one per reaction step. Python returns an
/// empty tuple for any failure; this returns RDKit's error. Mirrors
/// `rdChemReactions.ReactionsFromCDXMLBlock`.
pub fn reactions_from_cdxml_block(rxnblock: &str) -> Result<Vec<ChemicalReaction>> {
    reactions_from_cdxml_block_with_params(rxnblock, &ReactionsFromCDXMLBlockParams::default())
}

/// Mirrors `rdChemReactions.ReactionsFromCDXMLBlock(rxnblock, sanitize,
/// removeHs)`.
pub fn reactions_from_cdxml_block_with_params(
    rxnblock: &str,
    params: &ReactionsFromCDXMLBlockParams,
) -> Result<Vec<ChemicalReaction>> {
    reactions(ffi::reactions_from_cdxml_block(
        rxnblock,
        params.sanitize,
        params.remove_hs,
    ))
}

/// The reactions in a CDXML file. A file that cannot be opened is
/// `Error::Os`. Mirrors `rdChemReactions.ReactionsFromCDXMLFile`.
pub fn reactions_from_cdxml_file(filename: impl AsRef<Path>) -> Result<Vec<ChemicalReaction>> {
    reactions_from_cdxml_file_with_params(filename, &ReactionsFromCDXMLFileParams::default())
}

/// Mirrors `rdChemReactions.ReactionsFromCDXMLFile(filename, sanitize,
/// removeHs)`.
pub fn reactions_from_cdxml_file_with_params(
    filename: impl AsRef<Path>,
    params: &ReactionsFromCDXMLFileParams,
) -> Result<Vec<ChemicalReaction>> {
    reactions(ffi::reactions_from_cdxml_file(
        path_bytes(filename.as_ref()),
        params.sanitize,
        params.remove_hs,
    ))
}

/// An MDL rxn block for a reaction, V2000 with the agents among the
/// reactants. The V2000 writer depicts each template without a conformer,
/// and this returns `Error::Value` for one the depictor would mishandle (see
/// `chem::rd_depictor`). Mirrors `rdChemReactions.ReactionToRxnBlock`.
pub fn reaction_to_rxn_block(reaction: &ChemicalReaction) -> Result<String> {
    reaction_to_rxn_block_with_params(reaction, &ReactionToRxnBlockParams::default())
}

/// Mirrors `rdChemReactions.ReactionToRxnBlock(reaction, separateAgents,
/// forceV3000)`.
pub fn reaction_to_rxn_block_with_params(
    reaction: &ChemicalReaction,
    params: &ReactionToRxnBlockParams,
) -> Result<String> {
    to_string(ffi::reaction_to_rxn_block(
        reaction.as_ffi(),
        params.separate_agents,
        params.force_v3000,
    )?)
}

/// A V3000 rxn block for a reaction. Mirrors
/// `rdChemReactions.ReactionToV3KRxnBlock`.
pub fn reaction_to_v3k_rxn_block(reaction: &ChemicalReaction) -> Result<String> {
    reaction_to_v3k_rxn_block_with_separate_agents(reaction, false)
}

/// Mirrors `rdChemReactions.ReactionToV3KRxnBlock(reaction, separateAgents)`.
pub fn reaction_to_v3k_rxn_block_with_separate_agents(
    reaction: &ChemicalReaction,
    separate_agents: bool,
) -> Result<String> {
    to_string(ffi::reaction_to_v3k_rxn_block(
        reaction.as_ffi(),
        separate_agents,
    )?)
}

/// A Marvin (MRV) document for a reaction. Mirrors
/// `rdChemReactions.ReactionToMrvBlock`.
pub fn reaction_to_mrv_block(reaction: &ChemicalReaction) -> Result<String> {
    reaction_to_mrv_block_with_pretty_print(reaction, false)
}

/// Indents the XML when `pretty_print`. Mirrors
/// `rdChemReactions.ReactionToMrvBlock(reaction, prettyPrint)`.
pub fn reaction_to_mrv_block_with_pretty_print(
    reaction: &ChemicalReaction,
    pretty_print: bool,
) -> Result<String> {
    to_string(ffi::reaction_to_mrv_block(reaction.as_ffi(), pretty_print)?)
}

/// Writes a Marvin file. A file that cannot be opened is `Error::Os`.
/// Mirrors `rdChemReactions.ReactionToMrvFile`.
pub fn reaction_to_mrv_file(reaction: &ChemicalReaction, filename: impl AsRef<Path>) -> Result<()> {
    reaction_to_mrv_file_with_pretty_print(reaction, filename, false)
}

/// Mirrors `rdChemReactions.ReactionToMrvFile(reaction, filename,
/// prettyPrint)`.
pub fn reaction_to_mrv_file_with_pretty_print(
    reaction: &ChemicalReaction,
    filename: impl AsRef<Path>,
    pretty_print: bool,
) -> Result<()> {
    Ok(ffi::reaction_to_mrv_file(
        reaction.as_ffi(),
        path_bytes(filename.as_ref()),
        pretty_print,
    )?)
}

/// A reaction from a molecule whose fragments' atoms carry the `molRxnRole`
/// property (1 reactant, 2 product, 3 agent); other fragments are left out.
/// Mirrors `rdChemReactions.ReactionFromMolecule`.
pub fn reaction_from_molecule(mol: &Mol) -> Result<ChemicalReaction> {
    reaction(ffi::reaction_from_molecule(mol.as_ffi()))
}

/// One molecule holding a reaction's templates, their atoms marked with
/// `molRxnRole`. It sets that property on the reaction's templates too, as
/// RDKit does. Mirrors `rdChemReactions.ReactionToMolecule`.
pub fn reaction_to_molecule(reaction: &ChemicalReaction) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::reaction_to_molecule(reaction.as_ffi())?))
}

/// A reaction from PNG data's metadata: its first reaction pickle, SMILES,
/// SMARTS, or rxn block. Data without one is `Error::Parse`. Mirrors
/// `rdChemReactions.ReactionFromPNGString`.
pub fn reaction_from_png_string(data: &[u8]) -> Result<ChemicalReaction> {
    reaction(ffi::reaction_from_png_string(data))
}

/// Mirrors `rdChemReactions.ReactionFromPNGFile`.
pub fn reaction_from_png_file(fname: impl AsRef<Path>) -> Result<ChemicalReaction> {
    reaction(ffi::reaction_from_png_file(path_bytes(fname.as_ref())))
}

/// PNG data with the reaction added to its metadata, as a pickle and
/// reaction SMILES. PNG data without a complete IEND chunk is
/// `Error::Parse`. Mirrors `rdChemReactions.ReactionMetadataToPNGString`.
pub fn reaction_metadata_to_png_string(mol: &ChemicalReaction, pngdata: &[u8]) -> Result<Vec<u8>> {
    reaction_metadata_to_png_string_with_params(mol, pngdata, &Default::default())
}

/// Mirrors `rdChemReactions.ReactionMetadataToPNGString(mol, pngdata,
/// includePkl, includeSmiles, includeSmarts, includeRxn)`.
pub fn reaction_metadata_to_png_string_with_params(
    mol: &ChemicalReaction,
    pngdata: &[u8],
    params: &ReactionMetadataToPNGStringParams,
) -> Result<Vec<u8>> {
    Ok(ffi::reaction_metadata_to_png_string(
        mol.as_ffi(),
        pngdata,
        params.include_pkl,
        params.include_smiles,
        params.include_smarts,
        params.include_rxn,
    )?)
}

/// A PNG file's data with the reaction added to its metadata; the file is
/// not changed. Mirrors `rdChemReactions.ReactionMetadataToPNGFile`.
pub fn reaction_metadata_to_png_file(
    mol: &ChemicalReaction,
    filename: impl AsRef<Path>,
) -> Result<Vec<u8>> {
    reaction_metadata_to_png_file_with_params(mol, filename, &Default::default())
}

/// Mirrors `rdChemReactions.ReactionMetadataToPNGFile(mol, filename,
/// includePkl, includeSmiles, includeSmarts, includeMol)`.
pub fn reaction_metadata_to_png_file_with_params(
    mol: &ChemicalReaction,
    filename: impl AsRef<Path>,
    params: &ReactionMetadataToPNGFileParams,
) -> Result<Vec<u8>> {
    Ok(ffi::reaction_metadata_to_png_file(
        mol.as_ffi(),
        path_bytes(filename.as_ref()),
        params.include_pkl,
        params.include_smiles,
        params.include_smarts,
        params.include_mol,
    )?)
}
