// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/substruct.h"

#include <mutex>

#include <GraphMol/MolBundle.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/substruct.rs.h"

namespace rdkit_cxx {

namespace {

RDKit::SubstructMatchParameters params_from(
    const SubstructMatchParametersData &data) {
  RDKit::SubstructMatchParameters params;
  params.useChirality = data.use_chirality;
  params.useEnhancedStereo = data.use_enhanced_stereo;
  params.aromaticMatchesConjugated = data.aromatic_matches_conjugated;
  params.useGenericMatchers = data.use_generic_matchers;
  params.useQueryQueryMatches = data.use_query_query_matches;
  params.recursionPossible = data.recursion_possible;
  params.uniquify = data.uniquify;
  params.maxMatches = data.max_matches;
  params.maxRecursiveMatches = data.max_recursive_matches;
  params.numThreads = data.num_threads;
  for (const auto &prop : data.atom_properties) {
    params.atomProperties.emplace_back(prop);
  }
  for (const auto &prop : data.bond_properties) {
    params.bondProperties.emplace_back(prop);
  }
  return params;
}

// As convertMatch in substructmethods.h: the target atom for each query atom.
rust::Vec<std::uint32_t> convert_match(const RDKit::MatchVectType &match) {
  std::vector<std::uint32_t> res(match.size());
  for (const auto &[query_idx, mol_idx] : match) {
    res[query_idx] = mol_idx;
  }
  rust::Vec<std::uint32_t> out;
  for (auto idx : res) {
    out.push_back(idx);
  }
  return out;
}

// As helpHasSubstructMatch, helpGetSubstructMatch, and
// helpGetSubstructMatches, for a molecule or bundle and a molecule or bundle.
// Every match constructs a ControlCHandler (control_c_lock.h).
template <class Target, class Query>
bool has_match(const Target &target, const Query &query,
               const SubstructMatchParametersData &data) {
  auto params = params_from(data);
  params.maxMatches = 1;
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return !RDKit::SubstructMatch(target, query, params).empty();
}

template <class Target, class Query>
rust::Vec<std::uint32_t> get_match(const Target &target, const Query &query,
                                   const SubstructMatchParametersData &data) {
  auto params = params_from(data);
  params.maxMatches = 1;
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  auto matches = RDKit::SubstructMatch(target, query, params);
  return matches.empty() ? rust::Vec<std::uint32_t>()
                         : convert_match(matches[0]);
}

template <class Target, class Query>
rust::Vec<IndexList> get_matches(const Target &target, const Query &query,
                                 const SubstructMatchParametersData &data) {
  rust::Vec<IndexList> result;
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  for (const auto &match :
       RDKit::SubstructMatch(target, query, params_from(data))) {
    result.push_back(IndexList{convert_match(match)});
  }
  return result;
}

}  // namespace

SubstructMatchParametersData substruct_params_to_data(
    const RDKit::SubstructMatchParameters &p) {
  rust::Vec<rust::String> atom_properties;
  for (const auto &prop : p.atomProperties) {
    atom_properties.push_back(prop);
  }
  rust::Vec<rust::String> bond_properties;
  for (const auto &prop : p.bondProperties) {
    bond_properties.push_back(prop);
  }
  return SubstructMatchParametersData{p.useChirality,
                                      p.useEnhancedStereo,
                                      p.aromaticMatchesConjugated,
                                      p.useGenericMatchers,
                                      p.useQueryQueryMatches,
                                      p.recursionPossible,
                                      p.uniquify,
                                      p.maxMatches,
                                      p.maxRecursiveMatches,
                                      p.numThreads,
                                      std::move(atom_properties),
                                      std::move(bond_properties)};
}

RDKit::SubstructMatchParameters substruct_params_from_data(
    const SubstructMatchParametersData &data) {
  return params_from(data);
}

SubstructMatchParametersData substruct_match_parameters_default() {
  return substruct_params_to_data(RDKit::SubstructMatchParameters());
}

bool mol_has_substruct_match(const RWMol &target, const RWMol &query,
                             const SubstructMatchParametersData &params) {
  return has_match(target, query, params);
}

rust::Vec<std::uint32_t> mol_get_substruct_match(
    const RWMol &target, const RWMol &query,
    const SubstructMatchParametersData &params) {
  return get_match(target, query, params);
}

rust::Vec<IndexList> mol_get_substruct_matches(
    const RWMol &target, const RWMol &query,
    const SubstructMatchParametersData &params) {
  return get_matches(target, query, params);
}

bool mol_has_substruct_match_mol_bundle(
    const RWMol &target, const MolBundle &query,
    const SubstructMatchParametersData &params) {
  return has_match(target, query, params);
}

rust::Vec<std::uint32_t> mol_get_substruct_match_mol_bundle(
    const RWMol &target, const MolBundle &query,
    const SubstructMatchParametersData &params) {
  return get_match(target, query, params);
}

rust::Vec<IndexList> mol_get_substruct_matches_mol_bundle(
    const RWMol &target, const MolBundle &query,
    const SubstructMatchParametersData &params) {
  return get_matches(target, query, params);
}

bool mol_bundle_has_substruct_match(
    const MolBundle &target, const RWMol &query,
    const SubstructMatchParametersData &params) {
  return has_match(target, query, params);
}

rust::Vec<std::uint32_t> mol_bundle_get_substruct_match(
    const MolBundle &target, const RWMol &query,
    const SubstructMatchParametersData &params) {
  return get_match(target, query, params);
}

rust::Vec<IndexList> mol_bundle_get_substruct_matches(
    const MolBundle &target, const RWMol &query,
    const SubstructMatchParametersData &params) {
  return get_matches(target, query, params);
}

bool mol_bundle_has_substruct_match_mol_bundle(
    const MolBundle &target, const MolBundle &query,
    const SubstructMatchParametersData &params) {
  return has_match(target, query, params);
}

rust::Vec<std::uint32_t> mol_bundle_get_substruct_match_mol_bundle(
    const MolBundle &target, const MolBundle &query,
    const SubstructMatchParametersData &params) {
  return get_match(target, query, params);
}

rust::Vec<IndexList> mol_bundle_get_substruct_matches_mol_bundle(
    const MolBundle &target, const MolBundle &query,
    const SubstructMatchParametersData &params) {
  return get_matches(target, query, params);
}

rust::Vec<std::uint32_t> resonance_mol_supplier_get_substruct_match(
    ResonanceMolSupplier &suppl, const RWMol &query, bool use_chirality,
    bool use_query_query_matches) {
  RDKit::SubstructMatchParameters params;
  params.useChirality = use_chirality;
  params.useQueryQueryMatches = use_query_query_matches;
  params.maxMatches = 1;
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  auto matches = RDKit::SubstructMatch(suppl, query, params);
  return matches.empty() ? rust::Vec<std::uint32_t>()
                         : convert_match(matches[0]);
}

rust::Vec<IndexList> resonance_mol_supplier_get_substruct_matches(
    ResonanceMolSupplier &suppl, const RWMol &query, bool uniquify,
    bool use_chirality, bool use_query_query_matches, std::uint32_t max_matches,
    int num_threads) {
  RDKit::SubstructMatchParameters params;
  params.uniquify = uniquify;
  params.useChirality = use_chirality;
  params.useQueryQueryMatches = use_query_query_matches;
  params.maxMatches = max_matches;
  // RDKit's threads would match outside control_c_mutex, which they cannot
  // take, so the match runs on this thread.
  (void)num_threads;
  params.numThreads = 1;
  rust::Vec<IndexList> result;
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  for (const auto &match : RDKit::SubstructMatch(suppl, query, params)) {
    result.push_back(IndexList{convert_match(match)});
  }
  return result;
}

}  // namespace rdkit_cxx
