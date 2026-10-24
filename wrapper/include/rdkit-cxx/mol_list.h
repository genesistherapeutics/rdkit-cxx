/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "rdkit-cxx/mol.h"

namespace rdkit_cxx {

// Molecules a function returns, handed to Rust one at a time (CONVENTIONS.md
// W6).
class MolList {
 public:
  void push(std::unique_ptr<RWMol> mol) { d_mols.push_back(std::move(mol)); }
  std::size_t size() const { return d_mols.size(); }
  // Moves out the molecule at i; each is taken once.
  std::unique_ptr<RWMol> take(std::size_t i) { return std::move(d_mols.at(i)); }

 private:
  std::vector<std::unique_ptr<RWMol>> d_mols;
};

}  // namespace rdkit_cxx
