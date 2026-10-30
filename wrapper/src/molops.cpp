// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/molops.h"

#include <string>
#include <vector>

#include <GraphMol/MolOps.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/molops.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

// As sanitizeMol in Code/GraphMol/Wrap/MolOps.cpp: with catch_errors, a
// failure returns the operation that failed instead of throwing.
std::uint32_t sanitize_mol(RWMol &mol, std::uint32_t sanitize_ops,
                           bool catch_errors) {
  unsigned int operation_that_failed = 0;
  if (catch_errors) {
    try {
      RDKit::MolOps::sanitizeMol(mol, operation_that_failed, sanitize_ops);
    } catch (...) {
    }
  } else {
    RDKit::MolOps::sanitizeMol(mol, operation_that_failed, sanitize_ops);
  }
  return operation_that_failed;
}

void kekulize(RWMol &mol, bool clear_aromatic_flags, bool canonical) {
  RDKit::MolOps::Kekulize(mol, clear_aromatic_flags, canonical);
}

AddHsParametersData add_hs_parameters_default() {
  RDKit::MolOps::AddHsParameters params;
  return AddHsParametersData{params.explicitOnly,
                             params.addCoords,
                             params.addResidueInfo,
                             params.skipQueries,
                             false,
                             {}};
}

// As addHs2 in Code/GraphMol/Wrap/MolOps.cpp: a copy gains the hydrogens.
// Python (pythonObjectToVect in Code/RDBoost/Wrap.h) treats an empty
// onlyOnAtoms list as None, adding hydrogens to every atom; here an empty list
// adds none.
std::unique_ptr<RWMol> add_hs(const RWMol &mol,
                              const AddHsParametersData &data) {
  RDKit::MolOps::AddHsParameters params;
  params.explicitOnly = data.explicit_only;
  params.addCoords = data.add_coords;
  params.addResidueInfo = data.add_residue_info;
  params.skipQueries = data.skip_queries;
  std::vector<unsigned> only_on;
  if (data.has_only_on_atoms) {
    for (auto idx : data.only_on_atoms) {
      check_index(idx, mol.getNumAtoms());
      only_on.push_back(idx);
    }
  }
  auto res = std::make_unique<RWMol>(mol);
  RDKit::MolOps::addHs(*res, params,
                       data.has_only_on_atoms ? &only_on : nullptr);
  return res;
}

rust::Vec<IndexList> get_mol_frags(const RWMol &mol) {
  std::vector<std::vector<int>> frags;
  RDKit::MolOps::getMolFrags(mol, frags);
  rust::Vec<IndexList> result;
  for (const auto &frag : frags) {
    IndexList list;
    for (auto idx : frag) {
      list.indices.push_back(static_cast<std::uint32_t>(idx));
    }
    result.push_back(std::move(list));
  }
  return result;
}

// Each fragment is copied into an RWMol (W1, W6). RDKit makes some
// fragments by removing the other atoms from a copy (getTheFrags), which
// trusts what check_removals checks.
std::unique_ptr<MolList> get_mol_frags_as_mols(const RWMol &mol,
                                               bool sanitize_frags) {
  check_removals(mol);
  std::vector<std::unique_ptr<RDKit::ROMol>> frags;
  RDKit::MolOps::getMolFrags(mol, frags, sanitize_frags);
  auto list = std::make_unique<MolList>();
  for (const auto &frag : frags) {
    list->push(std::make_unique<RWMol>(*frag));
  }
  return list;
}

RemoveHsParametersData remove_hs_parameters_default() {
  RDKit::MolOps::RemoveHsParameters params;
  return RemoveHsParametersData{params.removeDegreeZero,
                                params.removeHigherDegrees,
                                params.removeOnlyHNeighbors,
                                params.removeIsotopes,
                                params.removeAndTrackIsotopes,
                                params.removeDummyNeighbors,
                                params.removeDefiningBondStereo,
                                params.removeWithWedgedBond,
                                params.removeWithQuery,
                                params.removeMapped,
                                params.removeInSGroups,
                                params.showWarnings,
                                params.removeNonimplicit,
                                params.updateExplicitCount,
                                params.removeHydrides,
                                params.removeNontetrahedralNeighbors,
                                true};
}

std::unique_ptr<RWMol> remove_hs(const RWMol &mol,
                                 const RemoveHsParametersData &data) {
  RDKit::MolOps::RemoveHsParameters params;
  params.removeDegreeZero = data.remove_degree_zero;
  params.removeHigherDegrees = data.remove_higher_degrees;
  params.removeOnlyHNeighbors = data.remove_only_h_neighbors;
  params.removeIsotopes = data.remove_isotopes;
  params.removeAndTrackIsotopes = data.remove_and_track_isotopes;
  params.removeDummyNeighbors = data.remove_dummy_neighbors;
  params.removeDefiningBondStereo = data.remove_defining_bond_stereo;
  params.removeWithWedgedBond = data.remove_with_wedged_bond;
  params.removeWithQuery = data.remove_with_query;
  params.removeMapped = data.remove_mapped;
  params.removeInSGroups = data.remove_in_s_groups;
  params.showWarnings = data.show_warnings;
  params.removeNonimplicit = data.remove_nonimplicit;
  params.updateExplicitCount = data.update_explicit_count;
  params.removeHydrides = data.remove_hydrides;
  params.removeNontetrahedralNeighbors = data.remove_nontetrahedral_neighbors;
  // Removing atoms trusts what check_removals checks.
  check_removals(mol);
  auto res = std::make_unique<RWMol>(mol);
  RDKit::MolOps::removeHs(*res, params, data.sanitize);
  return res;
}

std::unique_ptr<RWMol> remove_all_hs(const RWMol &mol, bool sanitize) {
  // Removing atoms trusts what check_removals checks.
  check_removals(mol);
  auto res = std::make_unique<RWMol>(mol);
  RDKit::MolOps::removeAllHs(*res, sanitize);
  return res;
}

std::unique_ptr<RWMol> merge_query_hs(const RWMol &mol,
                                      bool merge_unmapped_only,
                                      bool merge_isotopes) {
  // Removing atoms trusts what check_removals checks.
  check_removals(mol);
  auto res = std::make_unique<RWMol>(mol);
  RDKit::MolOps::mergeQueryHs(*res, merge_unmapped_only, merge_isotopes);
  return res;
}

rust::Vec<bool> has_query_hs(const RWMol &mol) {
  auto [has, unmergeable] = RDKit::MolOps::hasQueryHs(mol);
  rust::Vec<bool> result;
  result.push_back(has);
  result.push_back(unmergeable);
  return result;
}

std::size_t mol_list_size(const MolList &list) { return list.size(); }

const std::string &mol_list_name(const MolList &list, std::size_t i) {
  return list.name(i);
}

std::unique_ptr<RWMol> mol_list_take(MolList &list, std::size_t i) {
  return list.take(i);
}

std::unique_ptr<BorrowedMols> borrowed_mols_new() {
  return std::make_unique<BorrowedMols>();
}

void borrowed_mols_push(BorrowedMols &list, rust::Str name, RWMol &mol) {
  list.push(std::string(name), mol);
}

void borrowed_mols_push_shared(BorrowedMols &list, rust::Str name,
                               const RWMol &mol) {
  list.push(std::string(name), const_cast<RWMol &>(mol));
}

}  // namespace rdkit_cxx
