// Copyright (c) 2026, Steven Kearnes
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

std::unique_ptr<RWMol> mol_from_smiles(rust::Str smiles);
std::unique_ptr<std::string> mol_to_smiles(const RWMol &mol);

}  // namespace rdkit_cxx
