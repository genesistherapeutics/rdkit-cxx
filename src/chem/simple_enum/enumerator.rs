// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.SimpleEnum.Enumerator`, a Python module: reaction
//! enumeration with reactant templates whose atoms are labeled with
//! functional groups, as in rxn files' atom values.

use std::path::PathBuf;

use crate::chem::all_chem::enumerate_library_from_reaction;
use crate::chem::rd_chem_reactions::{self, ChemicalReaction, PreprocessReactionResult};
use crate::chem::{
    CXSmilesFields, Mol, MolToCXSmilesParams, mol_to_cx_smiles_with_params,
    parse_mol_query_def_file,
};
use crate::{Error, Result};

/// Parameters for `preprocess_reaction_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct PreprocessReactionParams {
    /// A file of functional group queries, as `ParseMolQueryDefFile` reads
    /// them; `None` for RDKit's functional group hierarchy.
    pub func_group_filename: Option<PathBuf>,
    /// The atom property that holds the labels.
    pub prop_name: String,
}

impl Default for PreprocessReactionParams {
    fn default() -> PreprocessReactionParams {
        PreprocessReactionParams {
            func_group_filename: None,
            prop_name: "molFileValue".to_owned(),
        }
    }
}

/// Validates the reaction and, if it has no errors, replaces the reactant
/// atoms labeled in their `molFileValue` property (a functional group's
/// name, or names separated by commas) with the group's query. Returns the
/// numbers of warnings and errors, of reactant and product templates, and
/// each reactant's labels. A label with no group is `Error::Key`, as in
/// Python. Mirrors `Enumerator.PreprocessReaction(reaction)`.
pub fn preprocess_reaction(reaction: &mut ChemicalReaction) -> Result<PreprocessReactionResult> {
    preprocess_reaction_with_params(reaction, &PreprocessReactionParams::default())
}

/// As [`preprocess_reaction`], with the groups of a file. A file that cannot
/// be read is `Error::Os`, as Python's `IOError`. Mirrors
/// `Enumerator.PreprocessReaction(reaction, funcGroupFilename, propName)`.
pub fn preprocess_reaction_with_params(
    reaction: &mut ChemicalReaction,
    params: &PreprocessReactionParams,
) -> Result<PreprocessReactionResult> {
    let queries = match &params.func_group_filename {
        Some(file) => parse_mol_query_def_file(file).map_err(|_| Error::Os {
            message: format!("cannot open {}", file.display()),
        })?,
        None => Default::default(),
    };
    rd_chem_reactions::preprocess_reaction_with_params(
        reaction,
        &rd_chem_reactions::PreprocessReactionParams {
            queries,
            prop_name: params.prop_name.clone(),
        },
    )
}

/// The products of the reaction for each combination of building blocks,
/// one from each list, after `preprocess_reaction` with its defaults, made
/// as they are asked for. A reaction with validation errors is
/// `Error::Value` ("bad reaction"), and so is a number of lists other than
/// its reactants, as in Python. Mirrors
/// `Enumerator.EnumerateReaction(reaction, bbLists)`.
pub fn enumerate_reaction<'a>(
    reaction: &'a mut ChemicalReaction,
    bb_lists: &'a [Vec<Mol>],
) -> Result<Box<dyn Iterator<Item = Result<Vec<Mol>>> + 'a>> {
    enumerate_reaction_with_unique_products_only(reaction, bb_lists, false)
}

/// As [`enumerate_reaction`]; `unique_products_only` yields each set of
/// products, by its sorted CXSMILES with enhanced stereo, once, and skips
/// empty ones. Mirrors `Enumerator.EnumerateReaction(reaction, bbLists,
/// uniqueProductsOnly)`.
pub fn enumerate_reaction_with_unique_products_only<'a>(
    reaction: &'a mut ChemicalReaction,
    bb_lists: &'a [Vec<Mol>],
    unique_products_only: bool,
) -> Result<Box<dyn Iterator<Item = Result<Vec<Mol>>> + 'a>> {
    // Python preprocesses with the defaults, whatever its funcGroupFilename
    // and propName arguments are.
    let PreprocessReactionResult {
        n_error, n_reacts, ..
    } = rd_chem_reactions::preprocess_reaction(reaction)?;
    if n_error != 0 {
        return Err(Error::Value {
            message: "bad reaction".to_owned(),
        });
    }
    if bb_lists.len() != n_reacts as usize {
        return Err(Error::Value {
            message: format!(
                "{n_reacts} reactants in reaction, {} bb lists supplied",
                bb_lists.len()
            ),
        });
    }
    let ps = enumerate_library_from_reaction(reaction, bb_lists)?;
    if !unique_products_only {
        return Ok(Box::new(ps));
    }
    let params = MolToCXSmilesParams {
        flags: CXSmilesFields::CX_ENHANCEDSTEREO,
        ..Default::default()
    };
    let mut seen: Vec<String> = Vec::new();
    Ok(Box::new(ps.filter_map(move |entry| {
        let entry = match entry {
            Ok(entry) => entry,
            Err(e) => return Some(Err(e)),
        };
        if entry.is_empty() {
            return None;
        }
        let smis: Result<Vec<String>> = entry
            .iter()
            .map(|x| mol_to_cx_smiles_with_params(x, &params))
            .collect();
        let mut smis = match smis {
            Ok(smis) => smis,
            Err(e) => return Some(Err(e)),
        };
        smis.sort();
        let smi = smis.join(".");
        if seen.contains(&smi) {
            return None;
        }
        seen.push(smi);
        Some(Ok(entry))
    })))
}
