// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdChemReactions`: chemical reactions, reading and
//! writing them, running them on molecules, and their fingerprints.
//!
//! ```no_run
//! use rdkit_cxx::chem::mol_from_smiles;
//! use rdkit_cxx::chem::mol_to_smiles;
//! use rdkit_cxx::chem::rd_chem_reactions::reaction_from_smarts;
//!
//! let mut rxn = reaction_from_smarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]")?;
//! let reactants = [mol_from_smiles("C(=O)O")?, mol_from_smiles("CNC")?];
//! let products = rxn.run_reactants(&reactants)?;
//! assert_eq!(mol_to_smiles(&products[0][0])?, "CN(C)C=O");
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```
//!
//! RDKit's `ChemicalReactionException` and `ChemicalReactionParserException`
//! (Python's `ValueError`s) are `Error::Value` and `Error::Parse`
//! (CONVENTIONS.md E4).

mod enumerate;
mod fingerprints;
mod io;
mod reaction;
mod sanitize;

use std::collections::BTreeMap;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::ffi::chem_reactions::ffi;

pub use enumerate::*;
pub use fingerprints::*;
pub use io::*;
pub use reaction::{
    AddRecursiveQueriesToReactionParams, ChemicalReaction, ReactantLabels,
    RemoveUnmappedProductTemplatesParams, RemoveUnmappedReactantTemplatesParams, TemplateRef,
    Templates,
};
pub use sanitize::*;

use reaction::{AGENT, PRODUCT, REACTANT, group_labels};

/// Whether a molecule counts as an agent of a reaction: whether fewer than
/// `agent_threshold` of its heavy atoms are mapped. Mirrors
/// `rdChemReactions.IsReactionTemplateMoleculeAgent`.
pub fn is_reaction_template_molecule_agent(molecule: &Mol, agent_threshold: f64) -> bool {
    crate::error::expect_ok(ffi::is_reaction_template_molecule_agent(
        molecule.as_ffi(),
        agent_threshold,
    ))
}

/// Whether any reactant or product template has a mapped atom. Mirrors
/// `rdChemReactions.HasReactionAtomMapping`.
pub fn has_reaction_atom_mapping(rxn: &ChemicalReaction) -> bool {
    crate::error::expect_ok(ffi::has_reaction_atom_mapping(rxn.as_ffi()))
}

/// Whether the query reaction's reactant and product templates are
/// substructures of the reaction's. Mirrors
/// `rdChemReactions.HasReactionSubstructMatch`.
pub fn has_reaction_substruct_match(
    reaction: &ChemicalReaction,
    query_reaction: &ChemicalReaction,
) -> Result<bool> {
    has_reaction_substruct_match_with_include_agents(reaction, query_reaction, false)
}

/// Compares the agents too when `include_agents`. Mirrors
/// `rdChemReactions.HasReactionSubstructMatch(reaction, queryReaction,
/// includeAgents)`.
pub fn has_reaction_substruct_match_with_include_agents(
    reaction: &ChemicalReaction,
    query_reaction: &ChemicalReaction,
    include_agents: bool,
) -> Result<bool> {
    Ok(ffi::has_reaction_substruct_match(
        reaction.as_ffi(),
        query_reaction.as_ffi(),
        include_agents,
    )?)
}

fn has_template_substruct_match(
    reaction: &ChemicalReaction,
    query_reaction: &ChemicalReaction,
    kind: u8,
) -> Result<bool> {
    Ok(ffi::has_template_substruct_match(
        reaction.as_ffi(),
        query_reaction.as_ffi(),
        kind,
    )?)
}

/// Whether one of the query reaction's agent templates is a substructure of
/// one of the reaction's. Mirrors
/// `rdChemReactions.HasAgentTemplateSubstructMatch`.
pub fn has_agent_template_substruct_match(
    reaction: &ChemicalReaction,
    query_reaction: &ChemicalReaction,
) -> Result<bool> {
    has_template_substruct_match(reaction, query_reaction, AGENT)
}

/// Whether one of the query reaction's product templates is a substructure
/// of one of the reaction's. Mirrors
/// `rdChemReactions.HasProductTemplateSubstructMatch`.
pub fn has_product_template_substruct_match(
    reaction: &ChemicalReaction,
    query_reaction: &ChemicalReaction,
) -> Result<bool> {
    has_template_substruct_match(reaction, query_reaction, PRODUCT)
}

/// Whether one of the query reaction's reactant templates is a substructure
/// of one of the reaction's. Mirrors
/// `rdChemReactions.HasReactantTemplateSubstructMatch`.
pub fn has_reactant_template_substruct_match(
    reaction: &ChemicalReaction,
    query_reaction: &ChemicalReaction,
) -> Result<bool> {
    has_template_substruct_match(reaction, query_reaction, REACTANT)
}

/// Sets the `molInversionFlag` of the product template atoms from the
/// stereochemistry of the mapped atoms, as the parsers do. Mirrors
/// `rdChemReactions.UpdateProductsStereochemistry`.
pub fn update_products_stereochemistry(reaction: &mut ChemicalReaction) -> Result<()> {
    Ok(ffi::update_products_stereochemistry(reaction.pin_mut())?)
}

/// A reaction product cut down to the atoms the reaction added, with dummy
/// atoms where they were attached. Mirrors
/// `rdChemReactions.ReduceProductToSideChains`.
pub fn reduce_product_to_side_chains(product: &Mol) -> Result<Mol> {
    reduce_product_to_side_chains_with_add_dummy_atoms(product, true)
}

/// Adds no dummy atoms unless `add_dummy_atoms`. Mirrors
/// `rdChemReactions.ReduceProductToSideChains(product, addDummyAtoms)`.
pub fn reduce_product_to_side_chains_with_add_dummy_atoms(
    product: &Mol,
    add_dummy_atoms: bool,
) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::reduce_product_to_side_chains(
        product.as_ffi(),
        add_dummy_atoms,
    )?))
}

/// Removes the atom maps from every template. Mirrors
/// `rdChemReactions.RemoveMappingNumbersFromReactions`.
pub fn remove_mapping_numbers_from_reactions(reaction: &mut ChemicalReaction) {
    crate::error::expect_ok(ffi::remove_mapping_numbers_from_reactions(
        reaction.pin_mut(),
    ))
}

/// Parameters for `preprocess_reaction_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug)]
pub struct PreprocessReactionParams {
    /// The queries, by the names the atoms' labels use. Empty uses RDKit's
    /// functional group hierarchy.
    pub queries: BTreeMap<String, Mol>,
    /// The atom property holding the labels.
    pub prop_name: String,
}

impl Default for PreprocessReactionParams {
    fn default() -> PreprocessReactionParams {
        PreprocessReactionParams {
            queries: BTreeMap::new(),
            prop_name: "molFileValue".to_owned(),
        }
    }
}

/// What `preprocess_reaction` returns: Python's `nWarn, nError, nReacts,
/// nProds, reactantLabels`.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct PreprocessReactionResult {
    /// The number of warnings from validation.
    pub n_warn: u32,
    /// The number of errors from validation.
    pub n_error: u32,
    /// The number of reactant templates.
    pub n_reacts: u32,
    /// The number of product templates.
    pub n_prods: u32,
    /// Each reactant's labels; none when validation failed.
    pub reactant_labels: ReactantLabels,
}

/// Validates the reaction and, if it has no errors, replaces the reactant
/// atoms labeled in their `molFileValue` property with the functional group
/// queries the labels name. A label with no query is `Error::Key`. Mirrors
/// `rdChemReactions.PreprocessReaction`.
pub fn preprocess_reaction(reaction: &mut ChemicalReaction) -> Result<PreprocessReactionResult> {
    preprocess_reaction_with_params(reaction, &PreprocessReactionParams::default())
}

/// Mirrors `rdChemReactions.PreprocessReaction(reaction, queries,
/// propName)`.
pub fn preprocess_reaction_with_params(
    reaction: &mut ChemicalReaction,
    params: &PreprocessReactionParams,
) -> Result<PreprocessReactionResult> {
    let data = with_borrowed_mols_shared(
        params
            .queries
            .iter()
            .map(|(name, mol)| (name.as_str(), mol)),
        |list| ffi::preprocess_reaction(reaction.pin_mut(), list, &params.prop_name),
    )?;
    Ok(PreprocessReactionResult {
        n_warn: data.num_warnings,
        n_error: data.num_errors,
        n_reacts: data.num_reactants,
        n_prods: data.num_products,
        reactant_labels: group_labels(data.labels, data.num_label_lists)?,
    })
}

/// Parameters for `compute_2d_coords_for_reaction_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct Compute2DCoordsForReactionParams {
    /// The space between the reaction's components.
    pub spacing: f64,
    /// Computes conjugation and hybridization of the templates first.
    pub update_props: bool,
    /// Orients each template canonically.
    pub canon_orient: bool,
    /// The rotatable bonds flipped at random per sample.
    pub n_flips_per_sample: u32,
    /// The number of random samples.
    pub n_sample: u32,
    /// The sampling's random seed.
    pub sample_seed: i32,
    /// Permutes the bonds around atoms with four neighbors while sampling.
    pub permute_deg4_nodes: bool,
    /// The bond length, or RDKit's default when not positive.
    pub bond_length: f64,
}

impl Default for Compute2DCoordsForReactionParams {
    fn default() -> Compute2DCoordsForReactionParams {
        Compute2DCoordsForReactionParams {
            spacing: 1.0,
            update_props: true,
            canon_orient: true,
            n_flips_per_sample: 0,
            n_sample: 0,
            sample_seed: 0,
            permute_deg4_nodes: false,
            bond_length: -1.0,
        }
    }
}

/// Lays out the reactant and product templates in 2D, left to right.
/// Templates over 46340 atoms are `Error::Value`, since RDKit computes their
/// distance matrices in `i32`, as are templates the depictor would mishandle
/// and a bond length past its coordinate bound (see `chem::rd_depictor`).
/// Mirrors `rdChemReactions.Compute2DCoordsForReaction`.
pub fn compute_2d_coords_for_reaction(reaction: &mut ChemicalReaction) -> Result<()> {
    compute_2d_coords_for_reaction_with_params(reaction, &Default::default())
}

/// Mirrors `rdChemReactions.Compute2DCoordsForReaction(reaction, spacing,
/// updateProps, canonOrient, nFlipsPerSample, nSample, sampleSeed,
/// permuteDeg4Nodes, bondLength)`.
pub fn compute_2d_coords_for_reaction_with_params(
    reaction: &mut ChemicalReaction,
    params: &Compute2DCoordsForReactionParams,
) -> Result<()> {
    Ok(ffi::compute_2d_coords_for_reaction(
        reaction.pin_mut(),
        params.spacing,
        params.update_props,
        params.canon_orient,
        params.n_flips_per_sample,
        params.n_sample,
        params.sample_seed,
        params.permute_deg4_nodes,
        params.bond_length,
    )?)
}
