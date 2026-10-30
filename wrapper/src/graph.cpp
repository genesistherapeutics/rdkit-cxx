// Copyright (c) 2026 Genesis Molecular AI, Inc.
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

// As getResidue and getChainId in Code/GraphMol/Wrap/MolOps.cpp, which
// read only info whose type is PDBRESIDUE, but with dynamic_cast where
// Python static_casts (only an AtomPDBResidueInfo can have the type,
// CONVENTIONS.md U10).
const RDKit::AtomPDBResidueInfo *residue_info(const RDKit::Atom *atom) {
  const auto *info = atom->getMonomerInfo();
  if (!info || info->getMonomerType() != RDKit::AtomMonomerInfo::PDBRESIDUE) {
    return nullptr;
  }
  return dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info);
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

// As findAllSubgraphsOfLengthsMtoNHelper in Code/GraphMol/Wrap/MolOps.cpp,
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

// RDKit 2024.03's pathToSubmol fetches each bond with getBondWithIdx,
// whose range check would throw with the half-built molecule unowned, so the
// indices are checked first; Python raises from that check.
std::unique_ptr<RWMol> path_to_submol(const RWMol &mol,
                                      rust::Slice<const std::uint32_t> path,
                                      bool use_query,
                                      rust::Vec<IndexPair> &atom_map) {
  RDKit::PATH_TYPE bonds;
  for (auto idx : path) {
    check_index(idx, mol.getNumBonds());
    // PATH_TYPE holds ints; pathToSubmol converts them back to unsigned.
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

// As renumberAtomsHelper in Code/GraphMol/Wrap/MolOps.cpp, which rejects a
// short newOrder and values past the atoms, and (pythonObjectToVect in
// Code/RDBoost/Wrap.h treating an empty list as None) an empty one even for a
// molecule without atoms. Here newOrder needs one entry per atom, each an
// atom's index, and each atom once: RDKit 2024.03's renumberAtoms refuses
// only values past the atom count, so a value equal to it writes past its
// reverse order, and it accepts repeated atoms.
std::unique_ptr<RWMol> renumber_atoms(
    const RWMol &mol, rust::Slice<const std::uint32_t> new_order) {
  if (new_order.size() != mol.getNumAtoms()) {
    throw ValueErrorException(
        "newOrder must have one entry for each atom in the molecule");
  }
  std::vector<bool> seen(new_order.size(), false);
  for (auto idx : new_order) {
    if (idx >= new_order.size()) {
      throw ValueErrorException("newOrder has an index past the atoms");
    }
    if (seen[idx]) {
      throw ValueErrorException("newOrder names an atom more than once");
    }
    seen[idx] = true;
  }
  std::unique_ptr<RDKit::ROMol> result(
      RDKit::MolOps::renumberAtoms(mol, to_vector(new_order)));
  return std::make_unique<RWMol>(*result);
}

}  // namespace rdkit_cxx
