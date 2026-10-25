/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdSubstructLibrary: libraries of molecules searched by substructure.
//
// Each holder is owned by one Rust value, through the shared_ptr RDKit's
// library takes. A Rust library owns its holders and its search order, and
// each call builds an RDKit library (a view) that shares them for the call,
// so a holder Rust swaps or edits between calls is the one the next call
// reads.
//
// RDKit's searches split their work across std::async threads when asked
// for more than one: an exception in one (a stored SMILES that does not
// sanitize, a fingerprint index past the fingerprint holder's end) unwinds
// internalGetMatches, which destroys the threads' result vectors before the
// futures that wait for the other threads, so those threads write freed
// memory (Python aborts with heap corruption). Each thread would also
// construct a ControlCHandler, which control_c_mutex cannot order. The
// wrappers therefore always ask RDKit for one thread, which searches on the
// calling thread under control_c_mutex.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/SubstructLibrary/SubstructLibrary.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/generalized_substruct.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_bundle.h"
#include "rdkit-cxx/tautomer_query.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using SubstructLibrary = RDKit::SubstructLibrary;

// Shared struct, defined by the substruct bridge (src/ffi/substruct.rs).
struct SubstructMatchParametersData;

class MolHolderHandle {
 public:
  explicit MolHolderHandle(boost::shared_ptr<RDKit::MolHolderBase> p)
      : ptr(std::move(p)) {}
  boost::shared_ptr<RDKit::MolHolderBase> ptr;
};

class FpHolderHandle {
 public:
  explicit FpHolderHandle(boost::shared_ptr<RDKit::FPHolderBase> p)
      : ptr(std::move(p)) {}
  boost::shared_ptr<RDKit::FPHolderBase> ptr;
};

class KeyHolderHandle {
 public:
  explicit KeyHolderHandle(boost::shared_ptr<RDKit::KeyHolderBase> p)
      : ptr(std::move(p)) {}
  boost::shared_ptr<RDKit::KeyHolderBase> ptr;
};

// Molecule holders by kind: 0 MolHolder, 1 CachedMolHolder,
// 2 CachedSmilesMolHolder, 3 CachedTrustedSmilesMolHolder.
std::unique_ptr<MolHolderHandle> mol_holder_new(std::uint8_t kind);
std::uint8_t mol_holder_kind(const MolHolderHandle &holder);
std::uint32_t mol_holder_len(const MolHolderHandle &holder);
std::uint32_t mol_holder_add_mol(MolHolderHandle &holder, const RWMol &mol);
// CachedMolHolder's AddBinary and the SMILES holders' AddSmiles store the
// input unchecked, as RDKit does; RDKit parses it when it reads the
// molecule.
std::uint32_t mol_holder_add_binary(MolHolderHandle &holder,
                                    rust::Slice<const std::uint8_t> pickle);
std::uint32_t mol_holder_add_smiles(MolHolderHandle &holder, rust::Str smiles);
// A copy of the molecule; a stored SMILES that does not parse is
// SmilesParseException (Python returns None).
std::unique_ptr<RWMol> mol_holder_get_mol(const MolHolderHandle &holder,
                                          std::uint32_t idx);

// Fingerprint holders by kind: 0 PatternHolder, 1 TautomerPatternHolder.
std::unique_ptr<FpHolderHandle> fp_holder_new(std::uint8_t kind,
                                              std::uint32_t num_bits);
std::uint8_t fp_holder_kind(const FpHolderHandle &holder);
std::uint32_t fp_holder_num_bits(const FpHolderHandle &holder);
std::uint32_t fp_holder_len(const FpHolderHandle &holder);
std::uint32_t fp_holder_add_mol(FpHolderHandle &holder, const RWMol &mol);
std::uint32_t fp_holder_add_fingerprint(FpHolderHandle &holder,
                                        const BitVect &v);
std::unique_ptr<BitVect> fp_holder_get_fingerprint(const FpHolderHandle &holder,
                                                   std::uint32_t idx);
bool fp_holder_passes_filter(const FpHolderHandle &holder, std::uint32_t idx,
                             const BitVect &query);
std::unique_ptr<BitVect> fp_holder_make_fingerprint(
    const FpHolderHandle &holder, const RWMol &mol);

// KeyFromPropHolder, the only key holder.
std::unique_ptr<KeyHolderHandle> key_holder_new(rust::Str propname);
std::uint32_t key_holder_len(const KeyHolderHandle &holder);
std::uint32_t key_holder_add_mol(KeyHolderHandle &holder, const RWMol &mol);
std::uint32_t key_holder_add_key(KeyHolderHandle &holder, rust::Str key);
std::unique_ptr<std::string> key_holder_get_key(const KeyHolderHandle &holder,
                                                std::uint32_t idx);
std::unique_ptr<std::vector<std::string>> key_holder_get_keys(
    const KeyHolderHandle &holder, rust::Slice<const std::uint32_t> indices);
std::unique_ptr<std::string> key_holder_get_prop_name(
    const KeyHolderHandle &holder);

// A library sharing the holders, for one call.
std::unique_ptr<SubstructLibrary> substruct_library_view(
    const MolHolderHandle &mols);
void substruct_library_view_set_fp_holder(SubstructLibrary &view,
                                          const FpHolderHandle &fps);
void substruct_library_view_set_key_holder(SubstructLibrary &view,
                                           const KeyHolderHandle &keys);
void substruct_library_view_set_search_order(
    SubstructLibrary &view, rust::Slice<const std::uint32_t> order);

std::uint32_t substruct_library_add_mol(SubstructLibrary &view,
                                        const RWMol &mol);

// The searches, for a molecule, a bundle, a tautomer query, or an extended
// query molecule. With ranged false, the whole library (start_idx and
// end_idx are ignored). RDKit always runs them on one thread (above), and
// matches tautomer queries with its own TautomerQuery::substructOf, so
// their templates are checked (tautomer_checks.h).
rust::Vec<std::uint32_t> substruct_library_get_matches(
    const SubstructLibrary &view, const RWMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results);
rust::Vec<std::uint32_t> substruct_library_get_matches_mol_bundle(
    const SubstructLibrary &view, const MolBundle &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results);
std::uint32_t substruct_library_count_matches(
    const SubstructLibrary &view, const RWMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params);
std::uint32_t substruct_library_count_matches_mol_bundle(
    const SubstructLibrary &view, const MolBundle &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params);
bool substruct_library_has_match(const SubstructLibrary &view,
                                 const RWMol &query, bool ranged,
                                 std::uint32_t start_idx, std::uint32_t end_idx,
                                 const SubstructMatchParametersData &params);
bool substruct_library_has_match_mol_bundle(
    const SubstructLibrary &view, const MolBundle &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params);

rust::Vec<std::uint32_t> substruct_library_get_matches_tautomer_query(
    const SubstructLibrary &view, const TautomerQuery &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results);
std::uint32_t substruct_library_count_matches_tautomer_query(
    const SubstructLibrary &view, const TautomerQuery &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params);
bool substruct_library_has_match_tautomer_query(
    const SubstructLibrary &view, const TautomerQuery &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params);
rust::Vec<std::uint32_t> substruct_library_get_matches_extended_query_mol(
    const SubstructLibrary &view, const ExtendedQueryMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params, std::int32_t max_results);
std::uint32_t substruct_library_count_matches_extended_query_mol(
    const SubstructLibrary &view, const ExtendedQueryMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params);
bool substruct_library_has_match_extended_query_mol(
    const SubstructLibrary &view, const ExtendedQueryMol &query, bool ranged,
    std::uint32_t start_idx, std::uint32_t end_idx,
    const SubstructMatchParametersData &params);

// Python's AddPatterns: fingerprints of the view's molecules, which replace
// the holder's.
void add_patterns(const SubstructLibrary &view, FpHolderHandle &patterns);

bool substruct_library_can_serialize();
rust::Vec<std::uint8_t> substruct_library_serialize(
    const SubstructLibrary &view);
// A library from Serialize's bytes, checked; its parts, which Rust takes.
std::unique_ptr<SubstructLibrary> substruct_library_from_pickle(
    rust::Slice<const std::uint8_t> data);
std::unique_ptr<MolHolderHandle> substruct_library_mol_holder(
    const SubstructLibrary &lib);
// Null when the library has none.
std::unique_ptr<FpHolderHandle> substruct_library_fp_holder(
    const SubstructLibrary &lib);
std::unique_ptr<KeyHolderHandle> substruct_library_key_holder(
    const SubstructLibrary &lib);
rust::Vec<std::uint32_t> substruct_library_search_order(
    const SubstructLibrary &lib);

}  // namespace rdkit_cxx
