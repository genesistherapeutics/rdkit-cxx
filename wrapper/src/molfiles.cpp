// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#include "rdkit-cxx/molfiles.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace rdkit_cxx {

std::unique_ptr<RWMol> mol_from_smiles(rust::Str smiles) {
  std::string text(smiles);
  auto mol = RDKit::v2::SmilesParse::MolFromSmiles(text);
  if (!mol) {
    // Python returns None here (E2).
    throw RDKit::SmilesParse::SmilesParseException("could not parse SMILES: " +
                                                   text);
  }
  return mol;
}

std::unique_ptr<std::string> mol_to_smiles(const RWMol &mol) {
  return std::make_unique<std::string>(
      RDKit::MolToSmiles(mol, RDKit::SmilesWriteParams()));
}

}  // namespace rdkit_cxx
