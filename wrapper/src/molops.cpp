// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/molops.h"

#include <string>
#include <vector>

#include <GraphMol/MolOps.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/molops.rs.h"

namespace rdkit_cxx {

// As sanitizeMol in Code/GraphMol/nbWrap/MolOps.cpp: with catch_errors, a
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

// As addHs2 in Code/GraphMol/nbWrap/MolOps.cpp: a copy gains the hydrogens.
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

// Each fragment is copied into an RWMol (W1, W6).
std::unique_ptr<MolList> get_mol_frags_as_mols(const RWMol &mol,
                                               bool sanitize_frags) {
  std::vector<std::unique_ptr<RDKit::ROMol>> frags;
  RDKit::MolOps::getMolFrags(mol, frags, sanitize_frags);
  auto list = std::make_unique<MolList>();
  for (const auto &frag : frags) {
    list->push(std::make_unique<RWMol>(*frag));
  }
  return list;
}

std::size_t mol_list_size(const MolList &list) { return list.size(); }

std::unique_ptr<RWMol> mol_list_take(MolList &list, std::size_t i) {
  return list.take(i);
}

}  // namespace rdkit_cxx
