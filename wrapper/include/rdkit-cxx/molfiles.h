/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolfiles's parsers and writers: SMILES, SMARTS, mol blocks, and the other
// formats. PNG metadata is in rdkit-cxx/png.h, suppliers in
// rdkit-cxx/suppliers.h, and writer classes in rdkit-cxx/writers.h.

#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/FileParsers/CDXMLParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/chemdraw.h>

#include "rdkit-cxx/atom.h"
#include "rdkit-cxx/bond.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// RDKit's enums, by the names the bridge gives them.
using CXSmilesFields = RDKit::SmilesWrite::CXSmilesFields;
using RestoreBondDirOption = RDKit::RestoreBondDirOption;
using SCSRTemplateNames = RDKit::v2::FileParsers::SCSRTemplateNames;
using SCSRBaseHbondOptions = RDKit::v2::FileParsers::SCSRBaseHbondOptions;
using CDXMLFormat = RDKit::v2::CDXMLParser::CDXMLFormat;
using CDXFormat = RDKit::v2::CDXFormat;

// Shared structs, defined by the bridge (src/ffi/molfiles.rs).
struct SmilesParserParamsData;
struct SmartsParserParamsData;
struct SmilesWriteParamsData;
struct MolWriterParamsData;
struct MolFromSCSRParamsData;

// Each returns the C++ class's defaults (CONVENTIONS.md A3).
SmilesParserParamsData smiles_parser_params_default();
SmartsParserParamsData smarts_parser_params_default();
SmilesWriteParamsData smiles_write_params_default();
MolWriterParamsData mol_writer_params_default();
MolFromSCSRParamsData mol_from_scsr_params_default();

// Conversions the other rdmolfiles bridges share.
RDKit::SmilesWriteParams write_params(const SmilesWriteParamsData &data);
RDKit::MolWriterParams writer_params(const MolWriterParamsData &data);
// Builds a v1 SmilesParserParams whose replacements point into `repl`, which
// must outlive it; throws for a debugParse the parsers refuse.
RDKit::SmilesParserParams smiles_parser_params(
    const SmilesParserParamsData &data,
    std::map<std::string, std::string> &repl);

std::unique_ptr<RWMol> mol_from_smiles(rust::Str smiles,
                                       const SmilesParserParamsData &params);
std::unique_ptr<std::string> mol_to_smiles(const RWMol &mol,
                                           const SmilesWriteParamsData &params);
std::unique_ptr<RWMol> mol_from_smarts(rust::Str smarts,
                                       const SmartsParserParamsData &params);
std::unique_ptr<std::string> mol_to_smarts(const RWMol &mol,
                                           const SmilesWriteParamsData &params);
// These write the parsers' debug globals (the Rust functions are unsafe).
std::unique_ptr<Atom> atom_from_smiles(rust::Str smiles);
std::unique_ptr<Atom> atom_from_smarts(rust::Str smarts);
std::unique_ptr<Bond> bond_from_smiles(rust::Str smiles);
std::unique_ptr<Bond> bond_from_smarts(rust::Str smarts);
std::unique_ptr<std::string> mol_to_cx_smiles(
    const RWMol &mol, const SmilesWriteParamsData &params, std::uint32_t flags,
    RestoreBondDirOption restore_bond_dirs);
std::unique_ptr<std::string> mol_to_cx_smarts(const RWMol &mol,
                                              bool isomeric_smiles);
// The fragment functions check their indices (E3) and symbol counts.
std::unique_ptr<std::string> mol_fragment_to_smiles(
    const RWMol &mol, const SmilesWriteParamsData &params,
    rust::Slice<const std::uint32_t> atoms_to_use, bool has_bonds,
    rust::Slice<const std::uint32_t> bonds_to_use, bool has_atom_symbols,
    const rust::Vec<rust::String> &atom_symbols, bool has_bond_symbols,
    const rust::Vec<rust::String> &bond_symbols, bool cx);
std::unique_ptr<std::string> mol_fragment_to_smarts(
    const RWMol &mol, rust::Slice<const std::uint32_t> atoms_to_use,
    bool has_bonds, rust::Slice<const std::uint32_t> bonds_to_use,
    bool isomeric_smarts, bool cx);
std::unique_ptr<std::vector<std::string>> mol_to_random_smiles_vect(
    const RWMol &mol, std::uint32_t num_smiles, std::uint32_t random_seed,
    bool isomeric_smiles, bool kekule_smiles, bool all_bonds_explicit,
    bool all_hs_explicit);
rust::Vec<std::uint32_t> canonical_rank_atoms(const RWMol &mol, bool break_ties,
                                              bool include_chirality,
                                              bool include_isotopes,
                                              bool include_atom_maps,
                                              bool include_chiral_presence);
rust::Vec<std::int32_t> canonical_rank_atoms_in_fragment(
    const RWMol &mol, rust::Slice<const std::uint32_t> atoms_to_use,
    bool has_bonds, rust::Slice<const std::uint32_t> bonds_to_use,
    bool has_atom_symbols, const rust::Vec<rust::String> &atom_symbols,
    bool break_ties, bool include_chirality, bool include_isotopes,
    bool include_atom_maps, bool include_chiral_presence);
void canonicalize_enhanced_stereo(RWMol &mol);
void create_property_list(RWMol &mol, std::uint8_t kind, rust::Str prop_name,
                          rust::Str missing_value_marker,
                          std::uint32_t line_size);

std::unique_ptr<RWMol> mol_from_mol_block(rust::Str mol_block, bool sanitize,
                                          bool remove_hs, bool strict_parsing);
std::unique_ptr<RWMol> mol_from_mol_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs,
                                         bool strict_parsing);
std::unique_ptr<std::string> mol_to_mol_block(
    const RWMol &mol, const MolWriterParamsData &params);
void mol_to_mol_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     const MolWriterParamsData &params);
std::unique_ptr<std::string> mol_to_v2k_mol_block(
    const RWMol &mol, const MolWriterParamsData &params);
std::unique_ptr<std::string> mol_to_v3k_mol_block(
    const RWMol &mol, const MolWriterParamsData &params);
void mol_to_v3k_mol_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                         const MolWriterParamsData &params);
std::unique_ptr<RWMol> mol_from_scsr_block(rust::Str mol_block, bool sanitize,
                                           bool remove_hs,
                                           const MolFromSCSRParamsData &params);
std::unique_ptr<RWMol> mol_from_scsr_file(rust::Slice<const std::uint8_t> path,
                                          bool sanitize, bool remove_hs,
                                          const MolFromSCSRParamsData &params);

std::unique_ptr<RWMol> mol_from_tpl_block(rust::Str tpl_block, bool sanitize,
                                          bool skip_first_conf);
std::unique_ptr<RWMol> mol_from_tpl_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool skip_first_conf);
std::unique_ptr<std::string> mol_to_tpl_block(const RWMol &mol,
                                              rust::Str partial_charge_prop,
                                              bool write_first_conf_twice);
void mol_to_tpl_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     rust::Str partial_charge_prop,
                     bool write_first_conf_twice);
std::unique_ptr<RWMol> mol_from_mol2_block(rust::Str mol2_block, bool sanitize,
                                           bool remove_hs,
                                           bool cleanup_substructures);
std::unique_ptr<RWMol> mol_from_mol2_file(rust::Slice<const std::uint8_t> path,
                                          bool sanitize, bool remove_hs,
                                          bool cleanup_substructures);
std::unique_ptr<RWMol> mol_from_mrv_block(rust::Str mrv_block, bool sanitize,
                                          bool remove_hs);
std::unique_ptr<RWMol> mol_from_mrv_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs);
std::unique_ptr<std::string> mol_to_mrv_block(const RWMol &mol,
                                              bool include_stereo,
                                              std::int32_t conf_id,
                                              bool kekulize, bool pretty_print);
void mol_to_mrv_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     bool include_stereo, std::int32_t conf_id, bool kekulize,
                     bool pretty_print);
std::unique_ptr<RWMol> mol_from_xyz_block(rust::Str xyz_block);
std::unique_ptr<RWMol> mol_from_xyz_file(rust::Slice<const std::uint8_t> path);
std::unique_ptr<std::string> mol_to_xyz_block(const RWMol &mol,
                                              std::int32_t conf_id,
                                              std::uint32_t precision);
void mol_to_xyz_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     std::int32_t conf_id, std::uint32_t precision);
std::unique_ptr<RWMol> mol_from_pdb_block(rust::Str mol_block, bool sanitize,
                                          bool remove_hs, std::uint32_t flavor,
                                          bool proximity_bonding);
std::unique_ptr<RWMol> mol_from_pdb_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs,
                                         std::uint32_t flavor,
                                         bool proximity_bonding);
std::unique_ptr<std::string> mol_to_pdb_block(const RWMol &mol,
                                              std::int32_t conf_id,
                                              std::uint32_t flavor);
void mol_to_pdb_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     std::int32_t conf_id, std::uint32_t flavor);
std::unique_ptr<std::string> mol_to_cml_block(const RWMol &mol,
                                              std::int32_t conf_id,
                                              bool kekulize);
void mol_to_cml_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     std::int32_t conf_id, bool kekulize);
// Writes the SMILES parser's debug global (the Rust function is unsafe).
std::unique_ptr<RWMol> mol_from_rdkit_svg(rust::Str svg, bool sanitize,
                                          bool remove_hs);
std::unique_ptr<RWMol> mol_from_sequence(rust::Str text, bool sanitize,
                                         std::int32_t flavor);
std::unique_ptr<RWMol> mol_from_fasta(rust::Str text, bool sanitize,
                                      std::int32_t flavor);
std::unique_ptr<RWMol> mol_from_helm(rust::Str text, bool sanitize);
std::unique_ptr<std::string> mol_to_sequence(const RWMol &mol);
std::unique_ptr<std::string> mol_to_fasta(const RWMol &mol);
std::unique_ptr<std::string> mol_to_helm(const RWMol &mol);
bool has_chem_draw_cdx_support();
// Checks the coordinates RDKit converts to fixed point.
rust::Vec<std::uint8_t> mol_to_chem_draw_block(const RWMol &mol,
                                               CDXFormat format);

}  // namespace rdkit_cxx
