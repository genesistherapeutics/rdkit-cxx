/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolCatalog: hierarchical catalogs of molecules.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/MolCatalog/MolCatalog.h>
#include <GraphMol/MolCatalog/MolCatalogEntry.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using MolCatalog = RDKit::MolCatalog;
using MolCatalogEntry = RDKit::MolCatalogEntry;

std::unique_ptr<MolCatalog> create_mol_catalog();
// Pickles are checked first (catalog_checks.h).
std::unique_ptr<MolCatalog> mol_catalog_from_pickle(
    rust::Slice<const std::uint8_t> pickle);
std::uint32_t mol_catalog_get_num_entries(const MolCatalog &catalog);
std::uint32_t mol_catalog_get_fp_length(const MolCatalog &catalog);
rust::Vec<std::uint8_t> mol_catalog_serialize(const MolCatalog &catalog);
// The lookups check their indices (E3).
std::unique_ptr<std::string> mol_catalog_get_bit_description(
    const MolCatalog &catalog, std::uint32_t idx);
// -1 when no entry has the bit.
std::int64_t mol_catalog_get_bit_entry_id(const MolCatalog &catalog,
                                          std::uint32_t idx);
std::int32_t mol_catalog_get_entry_bit_id(const MolCatalog &catalog,
                                          std::uint32_t idx);
std::unique_ptr<std::string> mol_catalog_get_entry_description(
    const MolCatalog &catalog, std::uint32_t idx);
rust::Vec<std::uint32_t> mol_catalog_get_entry_down_ids(
    const MolCatalog &catalog, std::uint32_t idx);
std::uint32_t mol_catalog_add_entry(MolCatalog &catalog,
                                    const MolCatalogEntry &entry);
void mol_catalog_add_edge(MolCatalog &catalog, std::uint32_t id1,
                          std::uint32_t id2);

std::unique_ptr<MolCatalogEntry> mol_catalog_entry_new();
std::unique_ptr<MolCatalogEntry> mol_catalog_entry_from_pickle(
    rust::Slice<const std::uint8_t> pickle);
std::unique_ptr<std::string> mol_catalog_entry_get_description(
    const MolCatalogEntry &entry);
void mol_catalog_entry_set_description(MolCatalogEntry &entry, rust::Str val);
// A copy of the entry's molecule, or null when it has none.
std::unique_ptr<RWMol> mol_catalog_entry_get_mol(const MolCatalogEntry &entry);
void mol_catalog_entry_set_mol(MolCatalogEntry &entry, const RWMol &mol);
std::uint32_t mol_catalog_entry_get_order(const MolCatalogEntry &entry);
void mol_catalog_entry_set_order(MolCatalogEntry &entry, std::uint32_t order);
rust::Vec<std::uint8_t> mol_catalog_entry_serialize(
    const MolCatalogEntry &entry);

}  // namespace rdkit_cxx
