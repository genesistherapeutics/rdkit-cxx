/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdfragcatalog: catalogs of molecular fragments and fingerprints from them.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/FragCatalog/FragCatGenerator.h>
#include <GraphMol/FragCatalog/FragCatParams.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using FragCatParams = RDKit::FragCatParams;
using FragCatalog = RDKit::FragCatalog;

// Reads the functional groups from a file (A8).
std::unique_ptr<FragCatParams> frag_cat_params_new(
    std::uint32_t l_len, std::uint32_t u_len,
    rust::Slice<const std::uint8_t> fgroup_filename, double tol);
std::unique_ptr<std::string> frag_cat_params_get_type_string(
    const FragCatParams &params);
std::uint32_t frag_cat_params_get_upper_frag_length(
    const FragCatParams &params);
std::uint32_t frag_cat_params_get_lower_frag_length(
    const FragCatParams &params);
double frag_cat_params_get_tolerance(const FragCatParams &params);
std::uint32_t frag_cat_params_get_num_func_groups(const FragCatParams &params);
// A copy of a functional group; checks its index (E3).
std::unique_ptr<RWMol> frag_cat_params_get_func_group(
    const FragCatParams &params, std::uint32_t fid);
rust::Vec<std::uint8_t> frag_cat_params_serialize(const FragCatParams &params);

std::unique_ptr<FragCatalog> frag_catalog_new();
std::unique_ptr<FragCatalog> frag_catalog_from_params(
    const FragCatParams &params);
// Pickles are checked first (catalog_checks.h).
std::unique_ptr<FragCatalog> frag_catalog_from_pickle(
    rust::Slice<const std::uint8_t> pickle);
std::uint32_t frag_catalog_get_num_entries(const FragCatalog &catalog);
std::uint32_t frag_catalog_get_fp_length(const FragCatalog &catalog);
// A copy of the catalog's parameters, or null when it has none.
std::unique_ptr<FragCatParams> frag_catalog_get_catalog_params(
    const FragCatalog &catalog);
rust::Vec<std::uint8_t> frag_catalog_serialize(const FragCatalog &catalog);
// The lookups check their indices (E3); the bit lookups also refuse a bit no
// entry has.
std::unique_ptr<std::string> frag_catalog_get_bit_description(
    const FragCatalog &catalog, std::uint32_t idx);
std::uint32_t frag_catalog_get_bit_order(const FragCatalog &catalog,
                                         std::uint32_t idx);
rust::Vec<std::int32_t> frag_catalog_get_bit_func_group_ids(
    const FragCatalog &catalog, std::uint32_t idx);
// -1 when no entry has the bit.
std::int64_t frag_catalog_get_bit_entry_id(const FragCatalog &catalog,
                                           std::uint32_t idx);
rust::Vec<double> frag_catalog_get_bit_discrims(const FragCatalog &catalog,
                                                std::uint32_t idx);
std::int32_t frag_catalog_get_entry_bit_id(const FragCatalog &catalog,
                                           std::uint32_t idx);
std::unique_ptr<std::string> frag_catalog_get_entry_description(
    const FragCatalog &catalog, std::uint32_t idx);
std::uint32_t frag_catalog_get_entry_order(const FragCatalog &catalog,
                                           std::uint32_t idx);
rust::Vec<std::int32_t> frag_catalog_get_entry_func_group_ids(
    const FragCatalog &catalog, std::uint32_t idx);
rust::Vec<std::uint32_t> frag_catalog_get_entry_down_ids(
    const FragCatalog &catalog, std::uint32_t idx);

std::uint32_t frag_cat_generator_add_frags_from_mol(const RWMol &mol,
                                                    FragCatalog &fcat);
std::unique_ptr<BitVect> frag_fp_generator_get_fp_for_mol(
    const RWMol &mol, const FragCatalog &fcat);

}  // namespace rdkit_cxx
