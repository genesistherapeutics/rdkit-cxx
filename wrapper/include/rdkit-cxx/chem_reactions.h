/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdChemReactions: chemical reactions, their parsers and writers, running
// them, and their fingerprints, sanitization, and depiction.
//
// Every template of a reaction Rust holds is an RWMol of its own: the
// wrappers copy the molecules Rust adds, and normalize_templates replaces
// any template an RDKit parser made as a plain ROMol, so the template views
// Rust lends are RWMols (CONVENTIONS.md W1) and no two reactions share a
// template.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionFingerprints.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using ChemicalReaction = RDKit::ChemicalReaction;
using FingerprintType = RDKit::FingerprintType;
using RDProps = RDKit::RDProps;
// A C enum, which cxx checks against the bridge's values.
using SanitizeRxnFlags = RDKit::RxnOps::SanitizeRxnFlags;

// Shared structs, defined by the bridges (src/ffi/chem_reactions.rs and the
// bridges it borrows them from).
struct ReactionFingerprintParamsData;
struct ReactantLabel;
struct PreprocessReactionData;
struct SmilesWriteParamsData;
struct SubstructMatchParametersData;
struct AdjustQueryParametersData;

// Product sets from running a reaction: one list of products per match, each
// with one product per product template.
class ProductSets {
 public:
  std::vector<std::vector<std::unique_ptr<RWMol>>> sets;
};

// Makes every template of `rxn` an RWMol that no other reaction shares, and
// refuses null templates (the V2000 rxn parser adds a null agent for an
// agent block cut short). Called on every reaction RDKit makes before Rust
// sees it.
void normalize_templates(ChemicalReaction &rxn);
// check_pickled_mol (mol_checks.h) on every template of a reaction read
// from a pickle, whose templates are molecule pickles (CONVENTIONS.md U10).
void check_reaction_monomer_types(const ChemicalReaction &rxn);

// Copies RDKit's product sets for Rust, clearing their matrix caches (U6).
std::unique_ptr<ProductSets> to_product_sets(
    const std::vector<RDKit::MOL_SPTR_VECT> &products);

std::size_t product_sets_size(const ProductSets &sets);
std::size_t product_sets_set_size(const ProductSets &sets, std::size_t i);
std::unique_ptr<RWMol> product_sets_take(ProductSets &sets, std::size_t i,
                                         std::size_t j);

// The ChemicalReaction class. Templates are selected by `kind`: 0 for
// reactants, 1 for products, 2 for agents.
std::unique_ptr<ChemicalReaction> chemical_reaction_new();
std::unique_ptr<ChemicalReaction> chemical_reaction_clone(
    const ChemicalReaction &rxn);
std::unique_ptr<ChemicalReaction> chemical_reaction_from_bin_str(
    rust::Slice<const std::uint8_t> bin_str);
rust::Vec<std::uint8_t> chemical_reaction_to_binary(const ChemicalReaction &rxn,
                                                    bool has_flags,
                                                    unsigned flags);
const RDProps &chemical_reaction_as_props(const ChemicalReaction &rxn);
unsigned chemical_reaction_get_num_templates(const ChemicalReaction &rxn,
                                             std::uint8_t kind);
// Checks `which` (E3).
const RWMol &chemical_reaction_get_template(const ChemicalReaction &rxn,
                                            std::uint8_t kind, unsigned which);
unsigned chemical_reaction_add_template(ChemicalReaction &rxn,
                                        std::uint8_t kind, const RWMol &mol);
// The removed templates, copied.
std::unique_ptr<MolList> chemical_reaction_remove_unmapped_templates(
    ChemicalReaction &rxn, std::uint8_t kind, double threshold_unmapped_atoms,
    bool move_to_agent_templates);
std::unique_ptr<MolList> chemical_reaction_remove_agent_templates(
    ChemicalReaction &rxn);
void chemical_reaction_initialize(ChemicalReaction &rxn, bool silent);
bool chemical_reaction_is_initialized(const ChemicalReaction &rxn);
// (numWarnings, numErrors).
rust::Vec<std::uint32_t> chemical_reaction_validate(ChemicalReaction &rxn,
                                                    bool silent);
// The runs initialize the reaction first, as Python's do, and run on copies
// of the reactants, which RDKit writes.
std::unique_ptr<ProductSets> chemical_reaction_run_reactants(
    ChemicalReaction &rxn, BorrowedMols &reactants, unsigned max_products);
std::unique_ptr<ProductSets> chemical_reaction_run_reactant(
    ChemicalReaction &rxn, const RWMol &reactant, unsigned reaction_idx);
bool chemical_reaction_run_reactant_in_place(ChemicalReaction &rxn,
                                             RWMol &reactant,
                                             bool remove_unmatched_atoms);
bool chemical_reaction_is_molecule_template(const ChemicalReaction &rxn,
                                            std::uint8_t kind,
                                            const RWMol &mol);
rust::Vec<IndexList> chemical_reaction_get_reacting_atoms(
    const ChemicalReaction &rxn, bool mapped_atoms_only);
// The queries are named by the list's names. Returns the labels when
// `get_labels`.
rust::Vec<ReactantLabel> chemical_reaction_add_recursive_queries(
    ChemicalReaction &rxn, BorrowedMols &queries, rust::Str prop_name,
    bool get_labels);
SubstructMatchParametersData chemical_reaction_get_substruct_params(
    const ChemicalReaction &rxn);
void chemical_reaction_set_substruct_params(
    ChemicalReaction &rxn, const SubstructMatchParametersData &params);

// Parsers and writers, as in Code/GraphMol/ChemReactions/Wrap.
std::unique_ptr<ChemicalReaction> reaction_from_smarts(
    rust::Str smarts, const rust::Vec<rust::String> &replacement_keys,
    const rust::Vec<rust::String> &replacement_values, bool use_smiles);
std::unique_ptr<std::string> reaction_to_smarts(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params);
std::unique_ptr<std::string> reaction_to_smiles(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params);
std::unique_ptr<std::string> reaction_to_cx_smarts(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params,
    std::uint32_t flags);
std::unique_ptr<std::string> reaction_to_cx_smiles(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params,
    std::uint32_t flags);
std::uint32_t cx_smiles_fields_all();
std::unique_ptr<ChemicalReaction> reaction_from_rxn_block(rust::Str rxn_block,
                                                          bool sanitize,
                                                          bool remove_hs,
                                                          bool strict_parsing);
std::unique_ptr<ChemicalReaction> reaction_from_rxn_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing);
std::unique_ptr<ChemicalReaction> reaction_from_mrv_block(rust::Str mrv,
                                                          bool sanitize,
                                                          bool remove_hs);
std::unique_ptr<ChemicalReaction> reaction_from_mrv_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs);
bool mrv_block_is_reaction(rust::Str mrv);
bool mrv_file_is_reaction(rust::Slice<const std::uint8_t> path);
std::unique_ptr<std::string> reaction_to_rxn_block(const ChemicalReaction &rxn,
                                                   bool separate_agents,
                                                   bool force_v3000);
std::unique_ptr<std::string> reaction_to_v3k_rxn_block(
    const ChemicalReaction &rxn, bool separate_agents);
std::unique_ptr<std::string> reaction_to_mrv_block(const ChemicalReaction &rxn,
                                                   bool pretty_print);
void reaction_to_mrv_file(const ChemicalReaction &rxn,
                          rust::Slice<const std::uint8_t> path,
                          bool pretty_print);
std::unique_ptr<ChemicalReaction> reaction_from_molecule(const RWMol &mol);
std::unique_ptr<RWMol> reaction_to_molecule(const ChemicalReaction &rxn);

// PNG metadata. The wrappers read files themselves and check the PNG chunks
// RDKit reads without checking (check_png in chem_reactions.cpp).
std::unique_ptr<ChemicalReaction> reaction_from_png_string(
    rust::Slice<const std::uint8_t> data);
std::unique_ptr<ChemicalReaction> reaction_from_png_file(
    rust::Slice<const std::uint8_t> path);
rust::Vec<std::uint8_t> reaction_metadata_to_png_string(
    const ChemicalReaction &rxn, rust::Slice<const std::uint8_t> png,
    bool include_pkl, bool include_smiles, bool include_smarts,
    bool include_rxn);
rust::Vec<std::uint8_t> reaction_metadata_to_png_file(
    const ChemicalReaction &rxn, rust::Slice<const std::uint8_t> path,
    bool include_pkl, bool include_smiles, bool include_smarts,
    bool include_rxn);

// Depiction, under lock_depiction() (U4).
void compute_2d_coords_for_reaction(ChemicalReaction &rxn, double spacing,
                                    bool update_props, bool canon_orient,
                                    unsigned n_flips_per_sample,
                                    unsigned n_sample, int sample_seed,
                                    bool permute_deg4_nodes,
                                    double bond_length);

// Fingerprints.
ReactionFingerprintParamsData reaction_fingerprint_params_default();
ReactionFingerprintParamsData default_difference_fp_params();
ReactionFingerprintParamsData default_structural_fp_params();
std::unique_ptr<SparseIntVect> create_difference_fingerprint_for_reaction(
    const ChemicalReaction &rxn, const ReactionFingerprintParamsData &params);
std::unique_ptr<BitVect> create_structural_fingerprint_for_reaction(
    const ChemicalReaction &rxn, const ReactionFingerprintParamsData &params);

// ReactionUtils.h and friends.
bool is_reaction_template_molecule_agent(const RWMol &mol,
                                         double agent_threshold);
bool has_reaction_atom_mapping(const ChemicalReaction &rxn);
bool has_reaction_substruct_match(const ChemicalReaction &rxn,
                                  const ChemicalReaction &query_rxn,
                                  bool include_agents);
bool has_template_substruct_match(const ChemicalReaction &rxn,
                                  const ChemicalReaction &query_rxn,
                                  std::uint8_t kind);
void update_products_stereochemistry(ChemicalReaction &rxn);
std::unique_ptr<RWMol> reduce_product_to_side_chains(const RWMol &product,
                                                     bool add_dummy_atoms);
void remove_mapping_numbers_from_reactions(ChemicalReaction &rxn);
// An empty query list uses RDKit's functional group hierarchy, as Python's
// does.
PreprocessReactionData preprocess_reaction(ChemicalReaction &rxn,
                                           BorrowedMols &queries,
                                           rust::Str prop_name);

// Sanitization (SanitizeRxn.h).
std::uint32_t sanitize_rxn(ChemicalReaction &rxn, std::uint32_t sanitize_ops,
                           const AdjustQueryParametersData &params,
                           bool catch_errors);
void sanitize_rxn_as_mols(ChemicalReaction &rxn, std::uint32_t sanitize_ops);
AdjustQueryParametersData get_default_adjust_params();
AdjustQueryParametersData get_chem_draw_rxn_adjust_params();
AdjustQueryParametersData match_only_at_rgroups_adjust_params();

}  // namespace rdkit_cxx
