// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_catalog.h"

#include <climits>

#include <GraphMol/MolCatalog/MolCatalogParams.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/catalog_checks.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

// A MolCatalogEntry's pickle, as MolCatalogEntry::initFromStream reads it:
// a molecule, a bit id, an order, and a description.
void check_entry(std::istream &ss, std::size_t size) {
  check_entry_mol(ss);
  read_int32(ss);  // the bit id
  read_int32(ss);  // the order
  check_entry_text(ss, size);
}

// MolCatalogEntry::toStream writes a description's length as an int32, which
// a length of INT_MAX or more overflows, and the readers refuse INT_MAX.
void check_description_length(const MolCatalogEntry &entry) {
  if (entry.getDescription().size() >= static_cast<std::size_t>(INT_MAX)) {
    throw ValueErrorException(
        "a catalog entry's description is too long to pickle");
  }
}

const MolCatalogEntry &entry_with_idx(const MolCatalog &catalog,
                                      std::uint32_t idx) {
  check_index(idx, catalog.getNumEntries());
  return *catalog.getEntryWithIdx(idx);
}

std::unique_ptr<std::string> copy(const std::string &text) {
  return std::make_unique<std::string>(text);
}

}  // namespace

std::unique_ptr<MolCatalog> create_mol_catalog() {
  RDKit::MolCatalogParams params;
  return std::make_unique<MolCatalog>(&params);
}

std::unique_ptr<MolCatalog> mol_catalog_from_pickle(
    rust::Slice<const std::uint8_t> pickle) {
  auto text = to_std(pickle);
  check_catalog_pickle<RDKit::MolCatalogParams>(text, check_entry);
  return std::make_unique<MolCatalog>(text);
}

std::uint32_t mol_catalog_get_num_entries(const MolCatalog &catalog) {
  return catalog.getNumEntries();
}

std::uint32_t mol_catalog_get_fp_length(const MolCatalog &catalog) {
  return catalog.getFPLength();
}

rust::Vec<std::uint8_t> mol_catalog_serialize(const MolCatalog &catalog) {
  for (unsigned i = 0; i < catalog.getNumEntries(); ++i) {
    check_description_length(*catalog.getEntryWithIdx(i));
  }
  return to_bytes(catalog.Serialize());
}

// The bit lookups check the bit against the fingerprint length, as RDKit
// does, but RDKit dereferences the entry it finds, and a pickle can give a
// catalog bits no entry has.
std::unique_ptr<std::string> mol_catalog_get_bit_description(
    const MolCatalog &catalog, std::uint32_t idx) {
  check_index(idx, catalog.getFPLength());
  const auto *entry = catalog.getEntryWithBitId(idx);
  if (!entry) {
    throw ValueErrorException("no catalog entry has bit " +
                              std::to_string(idx));
  }
  return copy(entry->getDescription());
}

std::int64_t mol_catalog_get_bit_entry_id(const MolCatalog &catalog,
                                          std::uint32_t idx) {
  check_index(idx, catalog.getFPLength());
  return catalog.getIdOfEntryWithBitId(idx);
}

std::int32_t mol_catalog_get_entry_bit_id(const MolCatalog &catalog,
                                          std::uint32_t idx) {
  return entry_with_idx(catalog, idx).getBitId();
}

std::unique_ptr<std::string> mol_catalog_get_entry_description(
    const MolCatalog &catalog, std::uint32_t idx) {
  return copy(entry_with_idx(catalog, idx).getDescription());
}

// getDownEntryList indexes the catalog's graph with idx unchecked.
rust::Vec<std::uint32_t> mol_catalog_get_entry_down_ids(
    const MolCatalog &catalog, std::uint32_t idx) {
  check_index(idx, catalog.getNumEntries());
  rust::Vec<std::uint32_t> ids;
  for (int id : catalog.getDownEntryList(idx)) {
    ids.push_back(static_cast<std::uint32_t>(id));
  }
  return ids;
}

std::uint32_t mol_catalog_add_entry(MolCatalog &catalog,
                                    const MolCatalogEntry &entry) {
  // The catalog owns the copy once its graph holds it.
  return catalog.addEntry(new MolCatalogEntry(entry));
}

void mol_catalog_add_edge(MolCatalog &catalog, std::uint32_t id1,
                          std::uint32_t id2) {
  check_index(id1, catalog.getNumEntries());
  check_index(id2, catalog.getNumEntries());
  catalog.addEdge(id1, id2);
}

std::unique_ptr<MolCatalogEntry> mol_catalog_entry_new() {
  return std::make_unique<MolCatalogEntry>();
}

std::unique_ptr<MolCatalogEntry> mol_catalog_entry_from_pickle(
    rust::Slice<const std::uint8_t> pickle) {
  auto text = to_std(pickle);
  auto ss = pickle_stream(text);
  check_entry(ss, text.size());
  return std::make_unique<MolCatalogEntry>(text);
}

std::unique_ptr<std::string> mol_catalog_entry_get_description(
    const MolCatalogEntry &entry) {
  return copy(entry.getDescription());
}

void mol_catalog_entry_set_description(MolCatalogEntry &entry, rust::Str val) {
  entry.setDescription(std::string(val));
}

std::unique_ptr<RWMol> mol_catalog_entry_get_mol(const MolCatalogEntry &entry) {
  const auto *mol = entry.getMol();
  if (!mol) {
    return nullptr;
  }
  auto result = std::make_unique<RWMol>(*mol);
  clear_matrix_caches(*result);
  return result;
}

void mol_catalog_entry_set_mol(MolCatalogEntry &entry, const RWMol &mol) {
  entry.setMol(new RDKit::ROMol(mol));  // the entry owns it
}

std::uint32_t mol_catalog_entry_get_order(const MolCatalogEntry &entry) {
  return entry.getOrder();
}

void mol_catalog_entry_set_order(MolCatalogEntry &entry, std::uint32_t order) {
  entry.setOrder(order);
}

rust::Vec<std::uint8_t> mol_catalog_entry_serialize(
    const MolCatalogEntry &entry) {
  check_description_length(entry);
  return to_bytes(entry.Serialize());
}

}  // namespace rdkit_cxx
