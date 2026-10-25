// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/file_checks.h"

#include <fstream>
#include <iterator>
#include <sstream>

#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

std::string read_file(const std::string &name) {
  std::ifstream in(name, std::ios_base::binary);
  if (!in || in.bad()) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  std::ostringstream bytes;
  bytes << in.rdbuf();
  if (in.bad()) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  return std::move(bytes).str();
}

std::string read_supplier_file(const std::string &name) {
  auto bytes = read_file(name);
  if (bytes.empty()) {
    throw RDKit::BadFileException("Invalid input file " + name);
  }
  return bytes;
}

void check_precision(unsigned precision) {
  if (precision > kMaxPrecision) {
    throw ValueErrorException("precision must be at most " +
                              std::to_string(kMaxPrecision));
  }
}

void check_debug_parse(int debug_parse) {
  if (debug_parse != 0) {
    throw ValueErrorException(
        "debugParse must be 0: RDKit stores it in a global that every SMILES "
        "and SMARTS parse reads without a lock");
  }
}

void check_pdb_residue_infos(const RDKit::ROMol &mol, bool every_info) {
  for (const auto atom : mol.atoms()) {
    const auto *info = atom->getMonomerInfo();
    if (!info || dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info)) {
      continue;
    }
    if (every_info ||
        info->getMonomerType() == RDKit::AtomMonomerInfo::PDBRESIDUE) {
      throw ValueErrorException(
          "atom " + std::to_string(atom->getIdx()) +
          " has monomer info that is not an AtomPDBResidueInfo" +
          (every_info ? "" : " but is marked PDBRESIDUE"));
    }
  }
}

}  // namespace rdkit_cxx
