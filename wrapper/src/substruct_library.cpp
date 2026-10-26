// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/substruct_library.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/substruct.rs.h"
#include "rdkit-cxx/ffi/substruct_library.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/substruct.h"
#include "rdkit-cxx/tautomer_checks.h"

namespace rdkit_cxx {

namespace {

using RDKit::CachedMolHolder;
using RDKit::CachedSmilesMolHolder;
using RDKit::CachedTrustedSmilesMolHolder;
using RDKit::KeyFromPropHolder;
using RDKit::PatternHolder;
using RDKit::TautomerPatternHolder;

// RDKit's searches count matches in int (internalGetMatches' counter, and
// molbundleGetMatches' running total), one per molecule, so no holder grows
// past INT_MAX entries.
void check_room(std::uint32_t size) {
  if (size >= static_cast<std::uint32_t>(INT_MAX)) {
    throw ValueErrorException("a holder holds at most INT_MAX entries");
  }
}

// The SMILES holders keep a molecule as its SMILES, whose writer removes
// atoms from copies of a molecule of several fragments (check_removals).
void check_smiles_holder_removals(const RDKit::MolHolderBase &holder,
                                  const RDKit::ROMol &mol) {
  if (dynamic_cast<const CachedSmilesMolHolder *>(&holder) ||
      dynamic_cast<const CachedTrustedSmilesMolHolder *>(&holder)) {
    check_removals(mol);
  }
}

// PatternFingerprintMol makes ExplicitBitVects of the holder's size, which
// count their bits in int (U3).
void check_num_bits(std::uint32_t num_bits) {
  if (num_bits > static_cast<std::uint32_t>(INT_MAX)) {
    throw ValueErrorException("a pattern holder has at most INT_MAX bits");
  }
}

const ExplicitBitVect &explicit_bit_vect(const BitVect &v) {
  const auto *bv = std::get_if<ExplicitBitVect>(&v.v);
  if (!bv) {
    throw ValueErrorException("pattern fingerprints are ExplicitBitVects");
  }
  return *bv;
}

PatternHolder &pattern_holder(const FpHolderHandle &holder) {
  // Every FPHolderBase Rust can make, or read from a pickle, is one.
  auto *res = dynamic_cast<PatternHolder *>(holder.ptr.get());
  if (!res) {
    throw ValueErrorException("the fingerprint holder is not a PatternHolder");
  }
  return *res;
}

// FPHolderBase::passesFilter calls AllProbeBitsMatch, which reads the
// query's blocks from the stored fingerprint (dynamic_bitset::is_subset_of
// checks the sizes only in debug builds), past the end of a shorter one.
// Every fingerprint a holder keeps has the holder's size (checked where it
// comes from outside RDKit), and so does every query RDKit makes for it.
void check_size(const ExplicitBitVect &v, std::uint32_t num_bits) {
  if (v.getNumBits() != num_bits) {
    throw ValueErrorException("a fingerprint for this holder has " +
                              std::to_string(num_bits) + " bits, not " +
                              std::to_string(v.getNumBits()));
  }
}

std::unique_ptr<RWMol> copy_out(const RDKit::ROMol &mol) {
  auto copy = std::make_unique<RWMol>(mol);
  clear_matrix_caches(*copy);
  return copy;
}

rust::Vec<std::uint32_t> to_vec(const std::vector<unsigned int> &values) {
  rust::Vec<std::uint32_t> res;
  res.reserve(values.size());
  for (auto value : values) {
    res.push_back(value);
  }
  return res;
}

// The indices a search reads, as RDKit's internalGetMatches computes them.
struct Range {
  unsigned int start;
  unsigned int end;
};

Range search_range(const SubstructLibrary &view, bool ranged,
                   std::uint32_t start_idx, std::uint32_t end_idx) {
  if (!ranged) {
    // As the Python overloads without a range; RDKit's preconditions
    // refuse an empty library.
    return Range{0, view.size()};
  }
  // RDKit's preconditions, as index errors (E3).
  check_index(start_idx, view.size());
  const auto &order = view.getSearchOrder();
  if (!order.empty()) {
    check_index(start_idx, static_cast<unsigned>(order.size()));
  }
  return Range{start_idx, end_idx};
}

}  // namespace

std::unique_ptr<MolHolderHandle> mol_holder_new(std::uint8_t kind) {
  switch (kind) {
    case 0:
      return std::make_unique<MolHolderHandle>(
          boost::make_shared<RDKit::MolHolder>());
    case 1:
      return std::make_unique<MolHolderHandle>(
          boost::make_shared<CachedMolHolder>());
    case 2:
      return std::make_unique<MolHolderHandle>(
          boost::make_shared<CachedSmilesMolHolder>());
    case 3:
      return std::make_unique<MolHolderHandle>(
          boost::make_shared<CachedTrustedSmilesMolHolder>());
    default:
      throw ValueErrorException("unknown molecule holder kind");
  }
}

std::uint8_t mol_holder_kind(const MolHolderHandle &holder) {
  const auto *base = holder.ptr.get();
  if (dynamic_cast<const RDKit::MolHolder *>(base)) {
    return 0;
  }
  if (dynamic_cast<const CachedMolHolder *>(base)) {
    return 1;
  }
  if (dynamic_cast<const CachedSmilesMolHolder *>(base)) {
    return 2;
  }
  if (dynamic_cast<const CachedTrustedSmilesMolHolder *>(base)) {
    return 3;
  }
  throw ValueErrorException("unknown molecule holder kind");
}

std::uint32_t mol_holder_len(const MolHolderHandle &holder) {
  return holder.ptr->size();
}

// CachedSmilesMolHolder writes SMILES, and CachedMolHolder pickles; neither
// matches substructures.
std::uint32_t mol_holder_add_mol(MolHolderHandle &holder, const RWMol &mol) {
  check_room(holder.ptr->size());
  check_smiles_holder_removals(*holder.ptr, mol);
  return holder.ptr->addMol(mol);
}

std::uint32_t mol_holder_add_binary(MolHolderHandle &holder,
                                    rust::Slice<const std::uint8_t> pickle) {
  auto *cached = dynamic_cast<CachedMolHolder *>(holder.ptr.get());
  if (!cached) {
    throw ValueErrorException("the holder is not a CachedMolHolder");
  }
  check_room(cached->size());
  return cached->addBinary(std::string(
      reinterpret_cast<const char *>(pickle.data()), pickle.size()));
}

std::uint32_t mol_holder_add_smiles(MolHolderHandle &holder, rust::Str smiles) {
  check_room(holder.ptr->size());
  std::string text(smiles);
  if (auto *h = dynamic_cast<CachedSmilesMolHolder *>(holder.ptr.get())) {
    return h->addSmiles(text);
  }
  if (auto *h =
          dynamic_cast<CachedTrustedSmilesMolHolder *>(holder.ptr.get())) {
    return h->addSmiles(text);
  }
  throw ValueErrorException("the holder is not a SMILES holder");
}

// The cached holders parse their molecules here (molecule pickles and
// SMILES, which RDKit's parsers check), so a pickle's monomer infos are
// checked as any molecule pickle's are (U10).
std::unique_ptr<RWMol> mol_holder_get_mol(const MolHolderHandle &holder,
                                          std::uint32_t idx) {
  check_index(idx, holder.ptr->size());
  auto mol = holder.ptr->getMol(idx);
  if (!mol) {
    throw RDKit::SmilesParse::SmilesParseException(
        "the molecule's SMILES does not parse");
  }
  check_monomer_types(*mol);
  return copy_out(*mol);
}

std::unique_ptr<FpHolderHandle> fp_holder_new(std::uint8_t kind,
                                              std::uint32_t num_bits) {
  check_num_bits(num_bits);
  switch (kind) {
    case 0:
      return std::make_unique<FpHolderHandle>(
          boost::make_shared<PatternHolder>(num_bits));
    case 1:
      return std::make_unique<FpHolderHandle>(
          boost::make_shared<TautomerPatternHolder>(num_bits));
    default:
      throw ValueErrorException("unknown fingerprint holder kind");
  }
}

std::uint8_t fp_holder_kind(const FpHolderHandle &holder) {
  return dynamic_cast<const TautomerPatternHolder *>(holder.ptr.get()) ? 1 : 0;
}

std::uint32_t fp_holder_num_bits(const FpHolderHandle &holder) {
  return pattern_holder(holder).getNumBits();
}

std::uint32_t fp_holder_len(const FpHolderHandle &holder) {
  return holder.ptr->size();
}

// Pattern fingerprints match their patterns (control_c_lock.h).
std::uint32_t fp_holder_add_mol(FpHolderHandle &holder, const RWMol &mol) {
  check_room(holder.ptr->size());
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return holder.ptr->addMol(mol);
}

std::uint32_t fp_holder_add_fingerprint(FpHolderHandle &holder,
                                        const BitVect &v) {
  const auto &bv = explicit_bit_vect(v);
  check_size(bv, fp_holder_num_bits(holder));
  check_room(holder.ptr->size());
  return holder.ptr->addFingerprint(bv);
}

std::unique_ptr<BitVect> fp_holder_get_fingerprint(const FpHolderHandle &holder,
                                                   std::uint32_t idx) {
  check_index(idx, holder.ptr->size());
  return std::make_unique<BitVect>(holder.ptr->getFingerprint(idx));
}

bool fp_holder_passes_filter(const FpHolderHandle &holder, std::uint32_t idx,
                             const BitVect &query) {
  check_index(idx, holder.ptr->size());
  const auto &bv = explicit_bit_vect(query);
  check_size(bv, holder.ptr->getFingerprint(idx).getNumBits());
  return holder.ptr->passesFilter(idx, bv);
}

std::unique_ptr<BitVect> fp_holder_make_fingerprint(
    const FpHolderHandle &holder, const RWMol &mol) {
  std::unique_ptr<ExplicitBitVect> fp;
  {
    std::lock_guard<std::mutex> control_c(control_c_mutex());
    fp.reset(holder.ptr->makeFingerprint(mol));
  }
  return std::make_unique<BitVect>(std::move(*fp));
}

std::unique_ptr<KeyHolderHandle> key_holder_new(rust::Str propname) {
  return std::make_unique<KeyHolderHandle>(
      boost::make_shared<KeyFromPropHolder>(std::string(propname)));
}

std::uint32_t key_holder_len(const KeyHolderHandle &holder) {
  return holder.ptr->size();
}

std::uint32_t key_holder_add_mol(KeyHolderHandle &holder, const RWMol &mol) {
  check_room(holder.ptr->size());
  return holder.ptr->addMol(mol);
}

std::uint32_t key_holder_add_key(KeyHolderHandle &holder, rust::Str key) {
  check_room(holder.ptr->size());
  return holder.ptr->addKey(std::string(key));
}

std::unique_ptr<std::string> key_holder_get_key(const KeyHolderHandle &holder,
                                                std::uint32_t idx) {
  check_index(idx, holder.ptr->size());
  return std::make_unique<std::string>(holder.ptr->getKey(idx));
}

// RDKit's getKeys reads each index with at(), whose out_of_range would be a
// RuntimeError; the indices are checked as E3 says.
std::unique_ptr<std::vector<std::string>> key_holder_get_keys(
    const KeyHolderHandle &holder, rust::Slice<const std::uint32_t> indices) {
  for (auto idx : indices) {
    check_index(idx, holder.ptr->size());
  }
  return std::make_unique<std::vector<std::string>>(holder.ptr->getKeys(
      std::vector<unsigned int>(indices.begin(), indices.end())));
}

std::unique_ptr<std::string> key_holder_get_prop_name(
    const KeyHolderHandle &holder) {
  const auto *keys = dynamic_cast<const KeyFromPropHolder *>(holder.ptr.get());
  if (!keys) {
    throw ValueErrorException("the key holder is not a KeyFromPropHolder");
  }
  return std::make_unique<std::string>(keys->getPropName());
}

std::unique_ptr<SubstructLibrary> substruct_library_view(
    const MolHolderHandle &mols) {
  return std::make_unique<SubstructLibrary>(mols.ptr);
}

void substruct_library_view_set_fp_holder(SubstructLibrary &view,
                                          const FpHolderHandle &fps) {
  view.getFpHolder() = fps.ptr;
  view.resetHolders();
}

void substruct_library_view_set_key_holder(SubstructLibrary &view,
                                           const KeyHolderHandle &keys) {
  view.getKeyHolder() = keys.ptr;
}

// RDKit indexes its molecules, and the bitset of molecules found, by each
// entry; setSearchOrder checks them with its own IndexErrorException, so
// they are checked here as E3 says. The molecule holder may have been
// swapped since the order was set.
void substruct_library_view_set_search_order(
    SubstructLibrary &view, rust::Slice<const std::uint32_t> order) {
  auto size = view.size();
  for (auto idx : order) {
    check_index(idx, size);
  }
  view.setSearchOrder(std::vector<unsigned int>(order.begin(), order.end()));
}

// Pattern fingerprints match their patterns (control_c_lock.h). As in
// RDKit, a holder that fails leaves the ones before it with the molecule.
std::uint32_t substruct_library_add_mol(SubstructLibrary &view,
                                        const RWMol &mol) {
  check_room(view.getMolHolder()->size());
  check_smiles_holder_removals(*view.getMolHolder(), mol);
  if (view.getFpHolder()) {
    check_room(view.getFpHolder()->size());
  }
  if (view.getKeyHolder()) {
    check_room(view.getKeyHolder()->size());
  }
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return view.addMol(mol);
}

namespace {

// One thread (substruct_library.h); RDKit's matches and the queries'
// pattern fingerprints construct ControlCHandlers (control_c_lock.h).
constexpr int kOneThread = 1;

std::vector<unsigned int> mol_matches(const SubstructLibrary &view,
                                      const RWMol &query, Range range,
                                      const RDKit::SubstructMatchParameters &p,
                                      int max_results) {
  const RDKit::ROMol &q = query;
  return view.getMatches(q, range.start, range.end, p, kOneThread, max_results);
}

// As molbundleGetMatches in SubstructLibrary.cpp, which searches for each
// member in turn, skipping the molecules earlier members matched, after
// subtracting its running total of matches from maxResults in int: the
// subtractions add up to as much as (members - 1) times the molecules
// matched, which overflows past 2^31. Here the arithmetic is in 64 bits,
// and each member's search is RDKit's, on a library whose search order
// leaves out the molecules earlier members matched, as RDKit's bitset of
// found molecules does.
std::vector<unsigned int> bundle_matches(
    const SubstructLibrary &view, const MolBundle &query, Range range,
    const RDKit::SubstructMatchParameters &p, std::int64_t max_results) {
  std::vector<unsigned int> res;
  const auto &members = query.getMols();
  if (members.empty()) {
    return res;
  }
  // RDKit's preconditions, which each member's search checks before it
  // returns for a maxResults of 0.
  view.getMatches(*members[0], range.start, range.end, p, kOneThread, 0);
  const auto &order = view.getSearchOrder();
  unsigned int end = std::min(range.end, view.size());
  if (!order.empty()) {
    end = std::min(end, static_cast<unsigned int>(order.size()));
  }
  std::vector<bool> found(view.size());
  for (const auto &member : members) {
    max_results -= static_cast<std::int64_t>(res.size());
    if (max_results == 0) {
      continue;
    }
    std::vector<unsigned int> remaining;
    for (unsigned int pos = range.start; pos < end; ++pos) {
      unsigned int idx = order.empty() ? pos : order[pos];
      if (!found[idx]) {
        remaining.push_back(idx);
      }
    }
    if (remaining.empty()) {
      continue;
    }
    SubstructLibrary rest(view.getMolHolder(), view.getFpHolder());
    rest.setSearchOrder(remaining);
    int limit = max_results < 0 ? -1 : static_cast<int>(max_results);
    for (auto idx : rest.getMatches(*member, 0,
                                    static_cast<unsigned int>(remaining.size()),
                                    p, kOneThread, limit)) {
      found[idx] = true;
      res.push_back(idx);
    }
  }
  return res;
}

std::vector<unsigned int> find(const SubstructLibrary &view, const RWMol &query,
                               Range range,
                               const RDKit::SubstructMatchParameters &p,
                               int max_results) {
  return mol_matches(view, query, range, p, max_results);
}

std::vector<unsigned int> find(const SubstructLibrary &view,
                               const MolBundle &query, Range range,
                               const RDKit::SubstructMatchParameters &p,
                               int max_results) {
  return bundle_matches(view, query, range, p, max_results);
}

std::vector<unsigned int> find(const SubstructLibrary &view,
                               const TautomerQuery &query, Range range,
                               const RDKit::SubstructMatchParameters &p,
                               int max_results) {
  check_rdkit_tautomer_match(query);
  return view.getMatches(query, range.start, range.end, p, kOneThread,
                         max_results);
}

std::vector<unsigned int> find(const SubstructLibrary &view,
                               const ExtendedQueryMol &query, Range range,
                               const RDKit::SubstructMatchParameters &p,
                               int max_results) {
  using XQ = ExtendedQueryMol;
  if (const auto *tautomers = std::get_if<XQ::TautomerQuery_T>(&query.xqmol)) {
    check_rdkit_tautomer_match(**tautomers);
  }
  if (const auto *bundle = std::get_if<XQ::TautomerBundle_T>(&query.xqmol)) {
    for (const auto &tautomers : **bundle) {
      check_rdkit_tautomer_match(*tautomers);
    }
  }
  return view.getMatches(query, range.start, range.end, p, kOneThread,
                         max_results);
}

template <class Query>
std::vector<unsigned int> search(const SubstructLibrary &view,
                                 const Query &query, bool ranged,
                                 std::uint32_t start_idx, std::uint32_t end_idx,
                                 const SubstructMatchParametersData &data,
                                 int max_results) {
  auto range = search_range(view, ranged, start_idx, end_idx);
  auto params = substruct_params_from_data(data);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return find(view, query, range, params, max_results);
}

template <class Query>
rust::Vec<std::uint32_t> get_matches(const SubstructLibrary &view,
                                     const Query &query, bool ranged,
                                     std::uint32_t start_idx,
                                     std::uint32_t end_idx,
                                     const SubstructMatchParametersData &data,
                                     std::int32_t max_results) {
  return to_vec(
      search(view, query, ranged, start_idx, end_idx, data, max_results));
}

// countMatches is getMatches' search with every match counted.
template <class Query>
std::uint32_t count_matches(const SubstructLibrary &view, const Query &query,
                            bool ranged, std::uint32_t start_idx,
                            std::uint32_t end_idx,
                            const SubstructMatchParametersData &data) {
  return static_cast<std::uint32_t>(
      search(view, query, ranged, start_idx, end_idx, data, -1).size());
}

// hasMatch is getMatches with one result.
template <class Query>
bool has_match(const SubstructLibrary &view, const Query &query, bool ranged,
               std::uint32_t start_idx, std::uint32_t end_idx,
               const SubstructMatchParametersData &data) {
  return !search(view, query, ranged, start_idx, end_idx, data, 1).empty();
}

}  // namespace

rust::Vec<std::uint32_t> substruct_library_get_matches(
    const SubstructLibrary &view, const RWMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results) {
  return get_matches(view, query, ranged, start_idx, end_idx, params,
                     max_results);
}

rust::Vec<std::uint32_t> substruct_library_get_matches_mol_bundle(
    const SubstructLibrary &view, const MolBundle &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results) {
  return get_matches(view, query, ranged, start_idx, end_idx, params,
                     max_results);
}

std::uint32_t substruct_library_count_matches(
    const SubstructLibrary &view, const RWMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params) {
  return count_matches(view, query, ranged, start_idx, end_idx, params);
}

std::uint32_t substruct_library_count_matches_mol_bundle(
    const SubstructLibrary &view, const MolBundle &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params) {
  return count_matches(view, query, ranged, start_idx, end_idx, params);
}

bool substruct_library_has_match(const SubstructLibrary &view,
                                 const RWMol &query, bool ranged,
                                 std::uint32_t start_idx, std::uint32_t end_idx,
                                 const SubstructMatchParametersData &params) {
  return has_match(view, query, ranged, start_idx, end_idx, params);
}

bool substruct_library_has_match_mol_bundle(
    const SubstructLibrary &view, const MolBundle &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params) {
  return has_match(view, query, ranged, start_idx, end_idx, params);
}

rust::Vec<std::uint32_t> substruct_library_get_matches_tautomer_query(
    const SubstructLibrary &view, const TautomerQuery &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results) {
  return get_matches(view, query, ranged, start_idx, end_idx, params,
                     max_results);
}

std::uint32_t substruct_library_count_matches_tautomer_query(
    const SubstructLibrary &view, const TautomerQuery &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params) {
  return count_matches(view, query, ranged, start_idx, end_idx, params);
}

bool substruct_library_has_match_tautomer_query(
    const SubstructLibrary &view, const TautomerQuery &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params) {
  return has_match(view, query, ranged, start_idx, end_idx, params);
}

rust::Vec<std::uint32_t> substruct_library_get_matches_extended_query_mol(
    const SubstructLibrary &view, const ExtendedQueryMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results) {
  return get_matches(view, query, ranged, start_idx, end_idx, params,
                     max_results);
}

std::uint32_t substruct_library_count_matches_extended_query_mol(
    const SubstructLibrary &view, const ExtendedQueryMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params) {
  return count_matches(view, query, ranged, start_idx, end_idx, params);
}

bool substruct_library_has_match_extended_query_mol(
    const SubstructLibrary &view, const ExtendedQueryMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params) {
  return has_match(view, query, ranged, start_idx, end_idx, params);
}

// RDKit's addPatterns resizes the holder's fingerprints to the library's
// size and fills them in place, on std::async threads: when reading a
// molecule throws (a stored SMILES that does not sanitize), the holder keeps
// null fingerprints that its getters and filters dereference (Python
// segfaults reading one). Here the fingerprints are made first, on the
// calling thread as RDKit's fillPatterns makes them (an empty molecule's for
// a molecule that does not parse), and replace the holder's only once all
// are made; RDKit's resize leaks the ones it replaces, which are freed here.
void add_patterns(const SubstructLibrary &view, FpHolderHandle &patterns) {
  const auto &mols = view.getMolecules();
  std::vector<std::unique_ptr<ExplicitBitVect>> fps;
  {
    std::lock_guard<std::mutex> control_c(control_c_mutex());
    for (unsigned int idx = 0; idx < mols.size(); ++idx) {
      auto mol = mols.getMol(idx);
      fps.emplace_back(
          patterns.ptr->makeFingerprint(mol ? *mol : RDKit::ROMol()));
    }
  }
  auto &held = patterns.ptr->getFingerprints();
  for (auto *fp : held) {
    delete fp;
  }
  held.clear();
  for (auto &fp : fps) {
    held.push_back(fp.release());
  }
}

bool substruct_library_can_serialize() {
  return RDKit::SubstructLibraryCanSerialize();
}

rust::Vec<std::uint8_t> substruct_library_serialize(
    const SubstructLibrary &view) {
  if (!RDKit::SubstructLibraryCanSerialize()) {
    throw ValueErrorException("RDKit was built without Boost serialization");
  }
  auto text = view.Serialize();
  rust::Vec<std::uint8_t> res;
  res.reserve(text.size());
  for (char c : text) {
    res.push_back(static_cast<std::uint8_t>(c));
  }
  return res;
}

// What RDKit-written archives may still get wrong, and what the wrappers'
// checks above assume of every library.
std::unique_ptr<SubstructLibrary> substruct_library_from_pickle(
    rust::Slice<const std::uint8_t> data) {
  if (!RDKit::SubstructLibraryCanSerialize()) {
    throw ValueErrorException("RDKit was built without Boost serialization");
  }
  auto lib = std::make_unique<SubstructLibrary>(
      std::string(reinterpret_cast<const char *>(data.data()), data.size()));
  const auto &mols = lib->getMolHolder();
  if (!mols) {
    throw ValueErrorException("the library has no molecule holder");
  }
  check_room(mols->size());
  // MolHolder's molecules are molecule pickles (U10); the cached holders'
  // are parsed, and checked, when read.
  if (auto *held = dynamic_cast<RDKit::MolHolder *>(mols.get())) {
    for (const auto &mol : held->getMols()) {
      if (!mol) {
        throw ValueErrorException("the library holds a null molecule");
      }
      check_monomer_types(*mol);
    }
  }
  for (auto idx : lib->getSearchOrder()) {
    check_index(idx, mols->size());
  }
  if (const auto &fps = lib->getFpHolder()) {
    FpHolderHandle handle(fps);
    auto num_bits = pattern_holder(handle).getNumBits();
    check_num_bits(num_bits);
    check_room(fps->size());
    for (const auto *fp : fps->getFingerprints()) {
      if (!fp) {
        throw ValueErrorException("the library holds a null fingerprint");
      }
      check_size(*fp, num_bits);
    }
  }
  if (const auto &keys = lib->getKeyHolder()) {
    check_room(keys->size());
  }
  return lib;
}

std::unique_ptr<MolHolderHandle> substruct_library_mol_holder(
    const SubstructLibrary &lib) {
  return std::make_unique<MolHolderHandle>(lib.getMolHolder());
}

std::unique_ptr<FpHolderHandle> substruct_library_fp_holder(
    const SubstructLibrary &lib) {
  if (!lib.getFpHolder()) {
    return nullptr;
  }
  return std::make_unique<FpHolderHandle>(lib.getFpHolder());
}

std::unique_ptr<KeyHolderHandle> substruct_library_key_holder(
    const SubstructLibrary &lib) {
  if (!lib.getKeyHolder()) {
    return nullptr;
  }
  return std::make_unique<KeyHolderHandle>(lib.getKeyHolder());
}

rust::Vec<std::uint32_t> substruct_library_search_order(
    const SubstructLibrary &lib) {
  return to_vec(lib.getSearchOrder());
}

}  // namespace rdkit_cxx
