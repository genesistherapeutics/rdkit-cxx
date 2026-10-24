// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/graph.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/graph.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

rust::Vec<double> copy_matrix(const double *matrix, unsigned n) {
  rust::Vec<double> values;
  std::size_t size = static_cast<std::size_t>(n) * n;
  values.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    values.push_back(matrix[i]);
  }
  return values;
}

template <class Paths>
void append_paths(const Paths &paths, rust::Vec<IndexList> &lists) {
  for (const auto &path : paths) {
    rust::Vec<std::uint32_t> indices;
    indices.reserve(path.size());
    for (auto idx : path) {
      indices.push_back(static_cast<std::uint32_t>(idx));
    }
    lists.push_back(IndexList{std::move(indices)});
  }
}

template <class Paths>
rust::Vec<IndexList> to_lists(const Paths &paths) {
  rust::Vec<IndexList> lists;
  append_paths(paths, lists);
  return lists;
}

std::vector<unsigned> to_vector(rust::Slice<const std::uint32_t> indices) {
  return std::vector<unsigned>(indices.begin(), indices.end());
}

// As getResidue and getChainId in Code/GraphMol/nbWrap/MolOps.cpp, but with
// dynamic_cast: Python static_casts any monomer info whose type is
// PDBRESIDUE, and a plain AtomMonomerInfo can carry that type.
const RDKit::AtomPDBResidueInfo *residue_info(const RDKit::Atom *atom) {
  return dynamic_cast<const RDKit::AtomPDBResidueInfo *>(
      atom->getMonomerInfo());
}

std::string get_residue(const RDKit::ROMol &, const RDKit::Atom *atom) {
  const auto *info = residue_info(atom);
  return info ? info->getResidueName() : "";
}

std::string get_chain_id(const RDKit::ROMol &, const RDKit::Atom *atom) {
  const auto *info = residue_info(atom);
  return info ? info->getChainId() : "";
}

// getMolFragsWithQuery leaves a skipped atom's group as T(), so with a white
// list that keeps the "" group and skips a neighbor, it bonds to atom -1 and
// fails a range check. The groups are filtered here instead, which gives
// each kept group the same molecule.
std::unique_ptr<MolList> split_mol(const RWMol &mol,
                                   std::string (*query)(const RDKit::ROMol &,
                                                        const RDKit::Atom *),
                                   bool has_white_list,
                                   rust::Slice<const rust::String> white_list,
                                   bool negate_list) {
  std::vector<std::string> names;
  for (const auto &name : white_list) {
    names.emplace_back(name);
  }
  auto frags = RDKit::MolOps::getMolFragsWithQuery(mol, query, false);
  auto result = std::make_unique<MolList>();
  for (const auto &[name, frag] : frags) {
    bool listed = std::find(names.begin(), names.end(), name) != names.end();
    if (has_white_list && listed == negate_list) {
      continue;
    }
    result->push_named(name, std::make_unique<RWMol>(*frag));
  }
  return result;
}

void set_options(RDKit::SubsetOptions &options, const SubsetOptionsData &data) {
  options.sanitize = data.sanitize;
  options.clearComputedProps = data.clear_computed_props;
  options.copyAsQuery = data.copy_as_query;
  options.copyCoordinates = data.copy_coordinates;
  options.conformerIdx = static_cast<unsigned>(data.conformer_idx);
  options.method = data.method;
}

void set_info(SubsetInfoData &data, const RDKit::SubsetInfo &info) {
  data.atom_mapping.clear();
  for (const auto &[from, to] : info.atomMapping) {
    data.atom_mapping.push_back(IndexPair{from, to});
  }
  data.bond_mapping.clear();
  for (const auto &[from, to] : info.bondMapping) {
    data.bond_mapping.push_back(IndexPair{from, to});
  }
}

}  // namespace

// The matrices are copied out before the guard clears RDKit's cache (U6).
rust::Vec<double> get_distance_matrix(const RWMol &mol, bool use_bo,
                                      bool use_atom_wts) {
  check_matrix_size(mol);
  MatrixCacheGuard guard(mol);
  const double *matrix =
      RDKit::MolOps::getDistanceMat(mol, use_bo, use_atom_wts, true);
  return copy_matrix(matrix, mol.getNumAtoms());
}

rust::Vec<double> get_3d_distance_matrix(const RWMol &mol, std::int32_t conf_id,
                                         bool use_atom_wts) {
  check_matrix_size(mol);
  MatrixCacheGuard guard(mol);
  const double *matrix =
      RDKit::MolOps::get3DDistanceMat(mol, conf_id, use_atom_wts, true);
  return copy_matrix(matrix, mol.getNumAtoms());
}

// RDKit passes emptyVal to memset, so only its low byte reaches the matrix's
// doubles; here every unbonded pair, the diagonal included, is empty_val.
rust::Vec<double> get_adjacency_matrix(const RWMol &mol, bool use_bo,
                                       std::int32_t empty_val) {
  check_matrix_size(mol);
  MatrixCacheGuard guard(mol);
  const double *matrix =
      RDKit::MolOps::getAdjacencyMatrix(mol, use_bo, 0, true);
  auto values = copy_matrix(matrix, mol.getNumAtoms());
  std::size_t n = mol.getNumAtoms();
  std::vector<bool> bonded(n * n);
  for (const auto *bond : mol.bonds()) {
    std::size_t begin = bond->getBeginAtomIdx();
    std::size_t end = bond->getEndAtomIdx();
    bonded[begin * n + end] = true;
    bonded[end * n + begin] = true;
  }
  for (std::size_t i = 0; i < n * n; ++i) {
    if (!bonded[i]) {
      values[i] = empty_val;
    }
  }
  return values;
}

// A path visits each atom once, except that it may close a ring with its
// last step, so no path is longer than numAtoms + 1. RDKit's loops run to
// the requested length regardless, and add one to it for bond paths, so a
// longer request returns no paths here without calling RDKit.
rust::Vec<IndexList> find_all_paths_of_length_n(const RWMol &mol,
                                                std::uint32_t length,
                                                bool use_bonds, bool use_hs,
                                                std::int32_t rooted_at_atom,
                                                bool only_shortest_paths) {
  check_matrix_size(mol);
  if (length > static_cast<std::uint64_t>(mol.getNumAtoms()) + 1) {
    return {};
  }
  // onlyShortestPaths reads the cached distance matrix (U6).
  MatrixCacheGuard guard(mol);
  return to_lists(RDKit::findAllPathsOfLengthN(
      mol, length, use_bonds, use_hs, rooted_at_atom, only_shortest_paths));
}

// A subgraph of n bonds has n distinct bonds, so a longer request returns
// none here, without RDKit's walk to that length.
rust::Vec<IndexList> find_all_subgraphs_of_length_n(
    const RWMol &mol, std::uint32_t length, bool use_hs,
    std::int32_t rooted_at_atom) {
  if (length > mol.getNumBonds()) {
    return {};
  }
  return to_lists(
      RDKit::findAllSubgraphsOfLengthN(mol, length, use_hs, rooted_at_atom));
}

// As findAllSubgraphsOfLengthsMtoNHelper in Code/GraphMol/nbWrap/MolOps.cpp,
// flattened: counts[i] is the number of subgraphs of length min + i. Lengths
// past the number of bonds have none; RDKit would loop to max, which never
// ends for UINT_MAX.
rust::Vec<IndexList> find_all_subgraphs_of_length_m_to_n(
    const RWMol &mol, std::uint32_t min, std::uint32_t max, bool use_hs,
    std::int32_t rooted_at_atom, rust::Vec<std::uint32_t> &counts) {
  if (min > max) {
    throw ValueErrorException("lowerLen > upperLen");
  }
  counts.clear();
  rust::Vec<IndexList> lists;
  auto upper = std::min(max, mol.getNumBonds());
  if (min > upper) {
    return lists;
  }
  auto paths = RDKit::findAllSubgraphsOfLengthsMtoN(mol, min, upper, use_hs,
                                                    rooted_at_atom);
  for (auto length = min; length <= upper; ++length) {
    const auto &group = paths[static_cast<int>(length)];
    counts.push_back(static_cast<std::uint32_t>(group.size()));
    append_paths(group, lists);
  }
  return lists;
}

rust::Vec<IndexList> find_unique_subgraphs_of_length_n(
    const RWMol &mol, std::uint32_t length, bool use_hs, bool use_bo,
    std::int32_t rooted_at_atom) {
  if (length > mol.getNumBonds()) {
    return {};
  }
  return to_lists(RDKit::findUniqueSubgraphsOfLengthN(mol, length, use_hs,
                                                      use_bo, rooted_at_atom));
}

rust::Vec<std::uint32_t> find_atom_environment_of_radius_n(
    const RWMol &mol, std::uint32_t radius, std::uint32_t rooted_at_atom,
    bool use_hs, bool enforce_size, rust::Vec<IndexPair> &atom_map) {
  std::unordered_map<unsigned, unsigned> map;
  auto path = RDKit::findAtomEnvironmentOfRadiusN(mol, radius, rooted_at_atom,
                                                  use_hs, enforce_size, &map);
  std::map<unsigned, unsigned> sorted(map.begin(), map.end());
  atom_map.clear();
  for (const auto &[atom, distance] : sorted) {
    atom_map.push_back(IndexPair{atom, distance});
  }
  rust::Vec<std::uint32_t> result;
  for (auto idx : path) {
    result.push_back(static_cast<std::uint32_t>(idx));
  }
  return result;
}

// copyMolSubset skips bond indices past the end, as Python does.
std::unique_ptr<RWMol> path_to_submol(const RWMol &mol,
                                      rust::Slice<const std::uint32_t> path,
                                      bool use_query,
                                      rust::Vec<IndexPair> &atom_map) {
  RDKit::PATH_TYPE bonds;
  for (auto idx : path) {
    // PATH_TYPE holds ints; copyMolSubset converts them back to unsigned.
    bonds.push_back(static_cast<int>(idx));
  }
  std::map<int, int> mapping;
  std::unique_ptr<RDKit::ROMol> submol(
      RDKit::Subgraphs::pathToSubmol(mol, bonds, use_query, mapping));
  atom_map.clear();
  for (const auto &[from, to] : mapping) {
    atom_map.push_back(IndexPair{static_cast<std::uint32_t>(from),
                                 static_cast<std::uint32_t>(to)});
  }
  return std::make_unique<RWMol>(*submol);
}

std::unique_ptr<MolList> split_mol_by_pdb_residues(
    const RWMol &mol, bool has_white_list,
    rust::Slice<const rust::String> white_list, bool negate_list) {
  return split_mol(mol, get_residue, has_white_list, white_list, negate_list);
}

std::unique_ptr<MolList> split_mol_by_pdb_chain_id(
    const RWMol &mol, bool has_white_list,
    rust::Slice<const rust::String> white_list, bool negate_list) {
  return split_mol(mol, get_chain_id, has_white_list, white_list, negate_list);
}

// As renumberAtomsHelper in Code/GraphMol/nbWrap/MolOps.cpp; renumberAtoms
// checks the values.
std::unique_ptr<RWMol> renumber_atoms(
    const RWMol &mol, rust::Slice<const std::uint32_t> new_order) {
  if (new_order.size() != mol.getNumAtoms()) {
    throw ValueErrorException(
        "newOrder must have one entry for each atom in the molecule");
  }
  std::unique_ptr<RDKit::ROMol> result(
      RDKit::MolOps::renumberAtoms(mol, to_vector(new_order)));
  return std::make_unique<RWMol>(*result);
}

SubsetOptionsData subset_options_default() {
  RDKit::SubsetOptions options;
  return SubsetOptionsData{options.sanitize,
                           options.clearComputedProps,
                           options.copyAsQuery,
                           options.copyCoordinates,
                           static_cast<std::int32_t>(options.conformerIdx),
                           options.method};
}

// The path form skips indices past the end, as Python does.
std::unique_ptr<RWMol> copy_mol_subset(const RWMol &mol,
                                       rust::Slice<const std::uint32_t> path,
                                       const SubsetOptionsData &data,
                                       SubsetInfoData &info_data) {
  RDKit::SubsetOptions options;
  set_options(options, data);
  RDKit::SubsetInfo info;
  auto result = RDKit::copyMolSubset(mol, to_vector(path), info, options);
  set_info(info_data, info);
  return result;
}

// copyMolSubset throws RDKit's IndexErrorException, which has no length,
// for indices past the end; the wrapper checks them first (E3).
std::unique_ptr<RWMol> copy_mol_subset_atoms_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> atom_indices,
    rust::Slice<const std::uint32_t> bond_indices,
    const SubsetOptionsData &data, SubsetInfoData &info_data) {
  for (auto idx : atom_indices) {
    check_index(idx, mol.getNumAtoms());
  }
  for (auto idx : bond_indices) {
    check_index(idx, mol.getNumBonds());
  }
  RDKit::SubsetOptions options;
  set_options(options, data);
  RDKit::SubsetInfo info;
  auto result = RDKit::copyMolSubset(mol, to_vector(atom_indices),
                                     to_vector(bond_indices), info, options);
  set_info(info_data, info);
  return result;
}

}  // namespace rdkit_cxx
