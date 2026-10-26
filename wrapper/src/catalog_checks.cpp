// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/catalog_checks.h"

#include <climits>

#include <GraphMol/MolPickler.h>
#include <GraphMol/ROMol.h>

#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

std::stringstream pickle_stream(const std::string &text) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  ss.write(text.data(), static_cast<std::streamsize>(text.size()));
  return ss;
}

std::int32_t read_int32(std::istream &ss) {
  std::int32_t value;
  RDKit::streamRead(ss, value);
  return value;
}

void check_entry_mol(std::istream &ss) {
  RDKit::ROMol mol;
  RDKit::MolPickler::molFromPickle(ss, mol);
  check_pickled_mol(mol);
}

void check_entry_text(std::istream &ss, std::size_t size) {
  auto length = read_int32(ss);
  auto pos = static_cast<std::size_t>(ss.tellg());
  // RDKit allocates length + 1 bytes, reads length, and writes a NUL at
  // length: a negative length writes before its buffer, INT_MAX overflows
  // the sum, and bytes past the end are left as they were allocated.
  if (length < 0 || length == INT_MAX ||
      static_cast<std::size_t>(length) > size - pos) {
    throw ValueErrorException(
        "a catalog pickle's description has a bad length: " +
        std::to_string(length));
  }
  ss.seekg(length, std::ios_base::cur);
}

}  // namespace rdkit_cxx
