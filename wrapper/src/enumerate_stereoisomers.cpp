// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/enumerate_stereoisomers.h"

#include <mutex>

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/MolOps.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/embedder.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/enumerate_stereoisomers.rs.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

namespace ES = RDKit::EnumerateStereoisomers;

ES::StereoEnumerationOptions from_data(
    const StereoEnumerationOptionsData &data) {
  ES::StereoEnumerationOptions options;
  options.tryEmbedding = data.try_embedding;
  options.onlyUnassigned = data.only_unassigned;
  options.onlyStereoGroups = data.only_stereo_groups;
  options.unique = data.unique;
  options.maxIsomers = data.max_isomers;
  options.randomSeed = data.random_seed;
  return options;
}

// The parameters of DGeomHelpers::EmbedMolecule(mol), which embeddable()
// (EnumerateStereoisomers.cpp) calls: a default EmbedParameters with the
// overload's default arguments, which leave its checks' fields as they are.
RDKit::DGeomHelpers::EmbedParameters embed_molecule_params() {
  RDKit::DGeomHelpers::EmbedParameters params;
  params.randomSeed = -1;
  params.useRandomCoords = false;
  return params;
}

}  // namespace

StereoEnumerationOptionsData stereo_enumeration_options_default() {
  ES::StereoEnumerationOptions options;
  StereoEnumerationOptionsData data{};
  data.try_embedding = options.tryEmbedding;
  data.only_unassigned = options.onlyUnassigned;
  data.only_stereo_groups = options.onlyStereoGroups;
  data.unique = options.unique;
  data.max_isomers = options.maxIsomers;
  data.random_seed = options.randomSeed;
  return data;
}

// The constructor finds the potential stereo elements and wedges bonds from
// the first conformer, if any; it draws no random numbers (the generator is
// the enumerator's own) and matches no substructures.
std::unique_ptr<StereoisomerEnumerator> stereoisomer_enumerator_new(
    const RWMol &mol, const StereoEnumerationOptionsData &data, bool verbose) {
  auto options = from_data(data);
  // next() writes each isomer, a copy of the molecule, as CXSMILES, which
  // removes atoms from copies of a molecule of several fragments
  // (check_removals); the enumerator keeps its own copy.
  check_removals(mol);
  auto res = std::make_unique<StereoisomerEnumerator>();
  res->enumerator =
      std::make_unique<ES::StereoisomerEnumerator>(mol, options, verbose);
  // Without stereo elements to flip, next() returns the molecule itself,
  // and does not embed it.
  res->embeds =
      options.tryEmbedding && res->enumerator->getStereoisomerCount() > 1;
  if (res->embeds) {
    res->unchecked = std::make_unique<RDKit::ROMol>(mol);
  }
  return res;
}

std::uint64_t stereoisomer_enumerator_get_stereoisomer_count(
    const StereoisomerEnumerator &enumerator) {
  return enumerator.enumerator->getStereoisomerCount();
}

// With tryEmbedding, next() embeds each isomer, with hydrogens added, as
// EmbedMolecule(mol): without a seed, so from RDKit's random generator
// (random_mutex), and with a ControlCHandler (control_c_mutex), taken in
// U4's order for the whole call. The embedder's checks are made once, on the
// molecule with hydrogens, whose atoms and fragments every isomer's share.
// Otherwise next() flips stereo elements, perceives stereo, and writes
// CXSMILES, none of which shares state with other threads.
std::unique_ptr<RWMol> stereoisomer_enumerator_next(
    StereoisomerEnumerator &enumerator) {
  std::unique_lock<std::mutex> control_c;
  std::unique_lock<std::mutex> random;
  if (enumerator.embeds) {
    if (enumerator.unchecked) {
      std::unique_ptr<RDKit::ROMol> with_hs(
          RDKit::MolOps::addHs(*enumerator.unchecked));
      check_embed_mol(*with_hs, embed_molecule_params());
      enumerator.unchecked.reset();
    }
    control_c = std::unique_lock<std::mutex>(control_c_mutex());
    random = std::unique_lock<std::mutex>(random_mutex());
  }
  auto isomer = enumerator.enumerator->next();
  if (!isomer) {
    return nullptr;
  }
  auto res = std::make_unique<RWMol>(*isomer);
  clear_matrix_caches(*res);
  return res;
}

}  // namespace rdkit_cxx
