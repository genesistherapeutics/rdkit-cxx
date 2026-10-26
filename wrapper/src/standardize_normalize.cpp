// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize_normalize.h"

#include <mutex>
#include <sstream>

#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/standardize.rs.h"
#include "rdkit-cxx/ffi/standardize_normalize.rs.h"
#include "rdkit-cxx/standardize.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace MS = RDKit::MolStandardize;

namespace {

// The objects' methods do not change them, but are not const.
template <class T>
T &mut(const T &object) {
  return const_cast<T &>(object);
}

MS::MetalDisconnectorOptions to_options(
    const MetalDisconnectorOptionsData &data) {
  MS::MetalDisconnectorOptions options;
  options.splitGrignards = data.split_grignards;
  options.splitAromaticC = data.split_aromatic_c;
  options.adjustCharges = data.adjust_charges;
  options.removeHapticDummies = data.remove_haptic_dummies;
  return options;
}

std::unique_ptr<std::string> smarts(const RDKit::ROMol *query) {
  if (!query) {
    throw ValueErrorException("the disconnector has no query");
  }
  return std::make_unique<std::string>(RDKit::MolToSmarts(*query));
}

}  // namespace

std::unique_ptr<Normalizer> normalizer_default() {
  return std::make_unique<Normalizer>();
}

std::unique_ptr<Normalizer> normalizer_new(
    rust::Slice<const std::uint8_t> normalize_filename,
    std::uint32_t max_restarts) {
  return std::make_unique<Normalizer>(to_path(normalize_filename),
                                      max_restarts);
}

// As normalizerFromDataAndParams in nbWrap/Normalize.cpp.
std::unique_ptr<Normalizer> normalizer_from_data(
    rust::Str param_data, const CleanupParametersData &params) {
  std::istringstream in{std::string(param_data)};
  return std::make_unique<Normalizer>(
      in, to_cleanup_parameters(params).maxRestarts);
}

std::unique_ptr<Normalizer> normalizer_from_params(
    const CleanupParametersData &params) {
  return std::unique_ptr<Normalizer>(
      MS::normalizerFromParams(to_cleanup_parameters(params)));
}

// Normalizing runs reactions, and disconnecting metals matches patterns
// (control_c_lock.h).
std::unique_ptr<RWMol> normalizer_normalize(const Normalizer &normalizer,
                                            const RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return adopt_mol(mut(normalizer).normalize(mol));
}

void normalizer_normalize_in_place(const Normalizer &normalizer, RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  mut(normalizer).normalizeInPlace(mol);
}

MetalDisconnectorOptionsData metal_disconnector_options_default() {
  MS::MetalDisconnectorOptions options;
  return MetalDisconnectorOptionsData{
      options.splitGrignards, options.splitAromaticC, options.adjustCharges,
      options.removeHapticDummies};
}

std::unique_ptr<MetalDisconnector> metal_disconnector_new(
    const MetalDisconnectorOptionsData &options) {
  return std::make_unique<MetalDisconnector>(to_options(options));
}

std::unique_ptr<std::string> metal_disconnector_get_metal_nof(
    const MetalDisconnector &disconnector) {
  return smarts(mut(disconnector).getMetalNof());
}

std::unique_ptr<std::string> metal_disconnector_get_metal_non(
    const MetalDisconnector &disconnector) {
  return smarts(mut(disconnector).getMetalNon());
}

void check_metal_query(const RDKit::ROMol &query) {
  auto atoms = query.getNumAtoms();
  if (atoms == 0) {
    return;
  }
  if (atoms != 2 || !query.getBondBetweenAtoms(0, 1)) {
    throw ValueErrorException(
        "a metal disconnector's query needs two atoms joined by a bond (or "
        "none): the metal and the atom bonded to it");
  }
}

void metal_disconnector_set_metal_nof(MetalDisconnector &disconnector,
                                      const RWMol &mol) {
  check_metal_query(mol);
  disconnector.setMetalNof(mol);
}

void metal_disconnector_set_metal_non(MetalDisconnector &disconnector,
                                      const RWMol &mol) {
  check_metal_query(mol);
  disconnector.setMetalNon(mol);
}

std::unique_ptr<RWMol> metal_disconnector_disconnect(
    const MetalDisconnector &disconnector, const RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return adopt_mol(mut(disconnector).disconnect(mol));
}

void metal_disconnector_disconnect_in_place(
    const MetalDisconnector &disconnector, RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  mut(disconnector).disconnectInPlace(mol);
}

std::unique_ptr<RWMol> disconnect_organometallics(const RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return adopt_mol(
      MS::disconnectOrganometallics(static_cast<const RDKit::ROMol &>(mol)));
}

std::unique_ptr<RWMol> disconnect_organometallics_with_options(
    const RWMol &mol, const MetalDisconnectorOptionsData &options) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return adopt_mol(MS::disconnectOrganometallics(
      static_cast<const RDKit::ROMol &>(mol), to_options(options)));
}

void disconnect_organometallics_in_place(RWMol &mol) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  MS::disconnectOrganometallicsInPlace(mol);
}

void disconnect_organometallics_in_place_with_options(
    RWMol &mol, const MetalDisconnectorOptionsData &options) {
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  MS::disconnectOrganometallicsInPlace(mol, to_options(options));
}

}  // namespace rdkit_cxx
