// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/conformer_parser.h"

#include <fstream>
#include <string>
#include <vector>

#include <GraphMol/Conformer.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

// RDKit installs Contrib/ConformerParser's library but not its header, so
// this is its readAmberTrajectory and addConformersFromList
// (ConformerParser.cpp), with checks: RDKit divides by three times the atom
// count, which is zero for a molecule with no atoms, and asserts that the
// file holds numConfs conformers. Python clears the conformers before
// reading the file; this reads it first, so a bad file leaves them alone.
rust::Vec<std::uint32_t> add_conformers_from_amber_trajectory(
    RWMol &mol, rust::Slice<const std::uint8_t> path, std::int32_t num_confs,
    bool clear_confs) {
  unsigned int num_atoms = mol.getNumAtoms();
  if (num_atoms == 0) {
    throw ValueErrorException("the molecule has no atoms");
  }
  auto file_name = to_path(path);
  std::ifstream in(file_name.c_str());
  if (!in || in.bad()) {
    throw RDKit::BadFileException("Bad input file " + file_name);
  }
  std::string title;
  std::getline(in, title);
  std::vector<double> coords;
  while (true) {
    double c;
    if (!(in >> c)) {
      if (!in.eof()) {
        throw ValueErrorException("Error while reading file");
      }
      break;
    }
    coords.push_back(c);
  }
  std::size_t per_conf = 3 * static_cast<std::size_t>(num_atoms);
  if (coords.size() % per_conf != 0) {
    throw ValueErrorException("Wrong number of coordinates");
  }
  std::size_t available = coords.size() / per_conf;
  std::size_t wanted =
      num_confs < 0 ? available : static_cast<std::size_t>(num_confs);
  if (wanted > available) {
    throw ValueErrorException(
        "numConfs is greater than the number of "
        "conformations in the file");
  }
  if (clear_confs) {
    mol.clearConformers();
  }
  rust::Vec<std::uint32_t> ids;
  for (std::size_t i = 0; i < wanted; ++i) {
    auto conf = std::make_unique<RDKit::Conformer>(num_atoms);
    for (unsigned int atom = 0; atom < num_atoms; ++atom) {
      const double *xyz = &coords[i * per_conf + 3 * atom];
      conf->setAtomPos(atom, RDGeom::Point3D(xyz[0], xyz[1], xyz[2]));
    }
    ids.push_back(mol.addConformer(conf.release(), true));
  }
  return ids;
}

}  // namespace rdkit_cxx
