// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/substruct.h"

#include <GraphMol/Substruct/SubstructMatch.h>

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
  params.aromaticMatchesSingleOrDouble = data.aromatic_matches_single_or_double;
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
  params.specifiedStereoQueryMatchesUnspecified =
      data.specified_stereo_query_matches_unspecified;
  params.extraAtomCheckOverridesDefaultCheck =
      data.extra_atom_check_overrides_default_check;
  params.extraBondCheckOverridesDefaultCheck =
      data.extra_bond_check_overrides_default_check;
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

}  // namespace

SubstructMatchParametersData substruct_match_parameters_default() {
  RDKit::SubstructMatchParameters p;
  rust::Vec<rust::String> atom_properties;
  for (const auto &prop : p.atomProperties) {
    atom_properties.push_back(prop);
  }
  rust::Vec<rust::String> bond_properties;
  for (const auto &prop : p.bondProperties) {
    bond_properties.push_back(prop);
  }
  return SubstructMatchParametersData{
      p.useChirality,
      p.useEnhancedStereo,
      p.aromaticMatchesConjugated,
      p.aromaticMatchesSingleOrDouble,
      p.useGenericMatchers,
      p.useQueryQueryMatches,
      p.recursionPossible,
      p.uniquify,
      p.maxMatches,
      p.maxRecursiveMatches,
      p.numThreads,
      std::move(atom_properties),
      std::move(bond_properties),
      p.specifiedStereoQueryMatchesUnspecified,
      p.extraAtomCheckOverridesDefaultCheck,
      p.extraBondCheckOverridesDefaultCheck};
}

bool mol_has_substruct_match(const RWMol &mol, const RWMol &query,
                             const SubstructMatchParametersData &data) {
  auto params = params_from(data);
  params.maxMatches = 1;
  return !RDKit::SubstructMatch(mol, query, params).empty();
}

rust::Vec<std::uint32_t> mol_get_substruct_match(
    const RWMol &mol, const RWMol &query,
    const SubstructMatchParametersData &data) {
  auto params = params_from(data);
  params.maxMatches = 1;
  auto matches = RDKit::SubstructMatch(mol, query, params);
  return matches.empty() ? rust::Vec<std::uint32_t>()
                         : convert_match(matches[0]);
}

rust::Vec<IndexList> mol_get_substruct_matches(
    const RWMol &mol, const RWMol &query,
    const SubstructMatchParametersData &data) {
  rust::Vec<IndexList> result;
  for (const auto &match : RDKit::SubstructMatch(mol, query, params_from(data))) {
    result.push_back(IndexList{convert_match(match)});
  }
  return result;
}

}  // namespace rdkit_cxx
