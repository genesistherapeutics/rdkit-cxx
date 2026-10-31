// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/frag_catalog.h"

#include <climits>
#include <mutex>
#include <tuple>

#include <GraphMol/FragCatalog/FragCatalogEntry.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/utils.h>

#include "rdkit-cxx/catalog_checks.h"
#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

using RDKit::FragCatalogEntry;

// A FragCatalogEntry's pickle, as FragCatalogEntry::initFromStream reads it:
// a molecule, a bit id, a description, an order, and a map from atoms to
// functional groups, each a count of int32s.
void check_entry(std::istream &ss, std::size_t size) {
  check_entry_mol(ss);
  read_int32(ss);  // the bit id
  check_entry_text(ss, size);
  read_int32(ss);  // the order
  auto num_atoms = read_int32(ss);
  for (std::int32_t i = 0; i < num_atoms; ++i) {
    read_int32(ss);  // the atom
    auto num_groups = read_int32(ss);
    for (std::int32_t j = 0; j < num_groups; ++j) {
      read_int32(ss);
    }
  }
}

const FragCatalogEntry &entry_with_idx(const FragCatalog &catalog,
                                       std::uint32_t idx) {
  check_index(idx, catalog.getNumEntries());
  return *catalog.getEntryWithIdx(idx);
}

// RDKit's bit lookups check the bit against the fingerprint length but
// dereference the entry they find, and a pickle can give a catalog bits no
// entry has.
const FragCatalogEntry &entry_with_bit(const FragCatalog &catalog,
                                       std::uint32_t idx) {
  check_index(idx, catalog.getFPLength());
  const auto *entry = catalog.getEntryWithBitId(idx);
  if (!entry) {
    throw ValueErrorException("no catalog entry has bit " +
                              std::to_string(idx));
  }
  return *entry;
}

// The functional groups' ids, in the order of their atoms, as
// GetEntryFuncGroupIds in Wrap/FragCatalog.cpp lists them.
rust::Vec<std::int32_t> func_group_ids(const FragCatalogEntry &entry) {
  rust::Vec<std::int32_t> ids;
  for (const auto &groups : entry.getFuncGroupMap()) {
    for (int id : groups.second) {
      ids.push_back(id);
    }
  }
  return ids;
}

// RDKit's generators dereference the catalog's parameters, which a
// FragCatalog() lacks.
const FragCatParams &params_of(const FragCatalog &catalog) {
  const auto *params = catalog.getCatalogParams();
  if (!params) {
    throw ValueErrorException("the fragment catalog has no parameters");
  }
  return *params;
}

// For each path of two or more bonds, RDKit's generators divide the path's
// product of primes by firstThousandPrimes[bond], with the bond's index in
// the molecule left when its functional groups are removed, which reads past
// the table at index NUM_PRIMES_AVAIL (computeIntVectPrimesProduct takes the
// index modulo the table's size; these divisions do not).
void check_num_bonds(const RWMol &mol, const FragCatParams &params) {
  if (params.getUpperFragLength() >= 2 &&
      mol.getNumBonds() > static_cast<unsigned>(RDKit::NUM_PRIMES_AVAIL)) {
    throw ValueErrorException(
        "fragments of two or more bonds need a molecule of at most " +
        std::to_string(RDKit::NUM_PRIMES_AVAIL) + " bonds");
  }
}

std::unique_ptr<std::string> copy(const std::string &text) {
  return std::make_unique<std::string>(text);
}

}  // namespace

std::unique_ptr<FragCatParams> frag_cat_params_new(
    std::uint32_t l_len, std::uint32_t u_len,
    rust::Slice<const std::uint8_t> fgroup_filename, double tol) {
  return std::make_unique<FragCatParams>(l_len, u_len, to_path(fgroup_filename),
                                         tol);
}

std::unique_ptr<std::string> frag_cat_params_get_type_string(
    const FragCatParams &params) {
  return copy(params.getTypeStr());
}

std::uint32_t frag_cat_params_get_upper_frag_length(
    const FragCatParams &params) {
  return params.getUpperFragLength();
}

std::uint32_t frag_cat_params_get_lower_frag_length(
    const FragCatParams &params) {
  return params.getLowerFragLength();
}

double frag_cat_params_get_tolerance(const FragCatParams &params) {
  return params.getTolerance();
}

std::uint32_t frag_cat_params_get_num_func_groups(const FragCatParams &params) {
  return params.getNumFuncGroups();
}

std::unique_ptr<RWMol> frag_cat_params_get_func_group(
    const FragCatParams &params, std::uint32_t fid) {
  check_index(fid, params.getNumFuncGroups());
  auto result = std::make_unique<RWMol>(*params.getFuncGroup(fid));
  clear_matrix_caches(*result);
  return result;
}

rust::Vec<std::uint8_t> frag_cat_params_serialize(const FragCatParams &params) {
  return to_bytes(params.Serialize());
}

std::unique_ptr<FragCatalog> frag_catalog_new() {
  return std::make_unique<FragCatalog>();
}

std::unique_ptr<FragCatalog> frag_catalog_from_params(
    const FragCatParams &params) {
  return std::make_unique<FragCatalog>(&params);
}

std::unique_ptr<FragCatalog> frag_catalog_from_pickle(
    rust::Slice<const std::uint8_t> pickle) {
  auto text = to_std(pickle);
  check_mol_pickle_versions(text);
  check_catalog_pickle<FragCatParams>(text, check_entry);
  return std::make_unique<FragCatalog>(text);
}

std::uint32_t frag_catalog_get_num_entries(const FragCatalog &catalog) {
  return catalog.getNumEntries();
}

std::uint32_t frag_catalog_get_fp_length(const FragCatalog &catalog) {
  return catalog.getFPLength();
}

std::unique_ptr<FragCatParams> frag_catalog_get_catalog_params(
    const FragCatalog &catalog) {
  const auto *params = catalog.getCatalogParams();
  if (!params) {
    return nullptr;
  }
  return std::make_unique<FragCatParams>(*params);
}

// FragCatalogEntry::toStream writes a description's length as an int32,
// which a length of INT_MAX or more overflows, and the readers refuse
// INT_MAX. RDKit checks that the catalog has parameters.
rust::Vec<std::uint8_t> frag_catalog_serialize(const FragCatalog &catalog) {
  for (unsigned i = 0; i < catalog.getNumEntries(); ++i) {
    if (catalog.getEntryWithIdx(i)->getDescription().size() >=
        static_cast<std::size_t>(INT_MAX)) {
      throw ValueErrorException(
          "a catalog entry's description is too long to pickle");
    }
  }
  return to_bytes(catalog.Serialize());
}

std::unique_ptr<std::string> frag_catalog_get_bit_description(
    const FragCatalog &catalog, std::uint32_t idx) {
  return copy(entry_with_bit(catalog, idx).getDescription());
}

std::uint32_t frag_catalog_get_bit_order(const FragCatalog &catalog,
                                         std::uint32_t idx) {
  return entry_with_bit(catalog, idx).getOrder();
}

rust::Vec<std::int32_t> frag_catalog_get_bit_func_group_ids(
    const FragCatalog &catalog, std::uint32_t idx) {
  return func_group_ids(entry_with_bit(catalog, idx));
}

std::int64_t frag_catalog_get_bit_entry_id(const FragCatalog &catalog,
                                           std::uint32_t idx) {
  check_index(idx, catalog.getFPLength());
  return catalog.getIdOfEntryWithBitId(idx);
}

rust::Vec<double> frag_catalog_get_bit_discrims(const FragCatalog &catalog,
                                                std::uint32_t idx) {
  auto discrims = entry_with_bit(catalog, idx).getDiscrims();
  rust::Vec<double> values;
  values.push_back(std::get<0>(discrims));
  values.push_back(std::get<1>(discrims));
  values.push_back(std::get<2>(discrims));
  return values;
}

std::int32_t frag_catalog_get_entry_bit_id(const FragCatalog &catalog,
                                           std::uint32_t idx) {
  return entry_with_idx(catalog, idx).getBitId();
}

std::unique_ptr<std::string> frag_catalog_get_entry_description(
    const FragCatalog &catalog, std::uint32_t idx) {
  return copy(entry_with_idx(catalog, idx).getDescription());
}

std::uint32_t frag_catalog_get_entry_order(const FragCatalog &catalog,
                                           std::uint32_t idx) {
  return entry_with_idx(catalog, idx).getOrder();
}

rust::Vec<std::int32_t> frag_catalog_get_entry_func_group_ids(
    const FragCatalog &catalog, std::uint32_t idx) {
  return func_group_ids(entry_with_idx(catalog, idx));
}

// getDownEntryList indexes the catalog's graph with idx unchecked.
rust::Vec<std::uint32_t> frag_catalog_get_entry_down_ids(
    const FragCatalog &catalog, std::uint32_t idx) {
  check_index(idx, catalog.getNumEntries());
  rust::Vec<std::uint32_t> ids;
  for (int id : catalog.getDownEntryList(idx)) {
    ids.push_back(static_cast<std::uint32_t>(id));
  }
  return ids;
}

// The generators match the functional groups, and every substructure match
// constructs a ControlCHandler (control_c_lock.h, U4). The catalog's entries
// hold parts of the molecule, which serializing it pickles; the parameters
// are written as text.
std::uint32_t frag_cat_generator_add_frags_from_mol(const RWMol &mol,
                                                    FragCatalog &fcat) {
  const auto &params = params_of(fcat);
  check_num_bonds(mol, params);
  check_picklable_parts(mol, params.getUpperFragLength());
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return RDKit::FragCatGenerator().addFragsFromMol(mol, &fcat);
}

std::unique_ptr<BitVect> frag_fp_generator_get_fp_for_mol(
    const RWMol &mol, const FragCatalog &fcat) {
  check_num_bonds(mol, params_of(fcat));
  // The fingerprint has a bit per bit id, and its int arithmetic needs at
  // most INT_MAX (U3).
  auto num_bits = fcat.getFPLength();
  if (num_bits > static_cast<unsigned>(INT_MAX)) {
    throw ValueErrorException(
        "a fragment fingerprint has at most 2^31 - 1 bits");
  }
  // RDKit allocates the fingerprint before these throw, and leaks it: the
  // catalog needs entries of one bond (an invariant RDKit checks), and a
  // pickle can give an entry a bit past the fingerprint's end, which setBit
  // refuses.
  fcat.getEntriesOfOrder(1);
  for (unsigned i = 0; i < fcat.getNumEntries(); ++i) {
    auto bit = fcat.getEntryWithIdx(i)->getBitId();
    if (bit >= 0 && static_cast<unsigned>(bit) >= num_bits) {
      throw ValueErrorException("a catalog entry has bit " +
                                std::to_string(bit) + " of " +
                                std::to_string(num_bits));
    }
  }
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::unique_ptr<ExplicitBitVect> fp(
      RDKit::FragFPGenerator().getFPForMol(mol, fcat));
  return wrap_bit_vect(std::move(fp));
}

}  // namespace rdkit_cxx
