/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/ROMol.h>

#include "rdkit-cxx/mol.h"

namespace rdkit_cxx {

// Molecules a function returns, handed to Rust one at a time (CONVENTIONS.md
// W6). A function returning a map pushes every molecule with its name.
class MolList {
 public:
  void push(std::unique_ptr<RWMol> mol) { d_mols.push_back(std::move(mol)); }
  void push_named(std::string name, std::unique_ptr<RWMol> mol) {
    d_names.push_back(std::move(name));
    push(std::move(mol));
  }
  std::size_t size() const { return d_mols.size(); }
  const std::string &name(std::size_t i) const { return d_names.at(i); }
  // Moves out the molecule at i; each is taken once.
  std::unique_ptr<RWMol> take(std::size_t i) { return std::move(d_mols.at(i)); }

 private:
  std::vector<std::string> d_names;
  std::vector<std::unique_ptr<RWMol>> d_mols;
};

// Molecules Rust lends to one call, as the shared_ptrs RDKit takes, with
// no-op deleters (W8). Rust builds the list, makes the call, and drops the
// list before its borrows of the molecules end.
class BorrowedMols {
 public:
  void push(std::string name, RDKit::ROMol &mol) {
    d_names.push_back(std::move(name));
    d_mols.emplace_back(&mol, [](RDKit::ROMol *) {});
  }
  const std::vector<std::string> &names() const { return d_names; }
  std::vector<boost::shared_ptr<RDKit::ROMol>> &mols() { return d_mols; }

 private:
  std::vector<std::string> d_names;
  std::vector<boost::shared_ptr<RDKit::ROMol>> d_mols;
};

}  // namespace rdkit_cxx
