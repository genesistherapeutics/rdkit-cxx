// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/deprotect.h"

#include <mutex>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

std::unique_ptr<std::string> copy(const std::string &text) {
  return std::make_unique<std::string>(text);
}

std::unique_ptr<RWMol> to_rwmol(std::unique_ptr<RDKit::ROMol> mol) {
  auto result = std::make_unique<RWMol>(*mol);
  clear_matrix_caches(*result);
  return result;
}

}  // namespace

std::unique_ptr<DeprotectData> deprotect_data_new(rust::Str deprotection_class,
                                                  rust::Str reaction_smarts,
                                                  rust::Str abbreviation,
                                                  rust::Str full_name) {
  // DeprotectData parses the reaction SMARTS with its CXSMILES extension
  // (file_checks.h).
  check_reaction_cx(std::string(reaction_smarts), true);
  return std::make_unique<DeprotectData>(
      std::string(deprotection_class), std::string(reaction_smarts),
      std::string(abbreviation), std::string(full_name));
}

std::unique_ptr<DeprotectData> deprotect_data_clone(const DeprotectData &data) {
  return std::make_unique<DeprotectData>(data);
}

std::unique_ptr<std::string> deprotect_data_deprotection_class(
    const DeprotectData &data) {
  return copy(data.deprotection_class);
}

std::unique_ptr<std::string> deprotect_data_full_name(
    const DeprotectData &data) {
  return copy(data.full_name);
}

std::unique_ptr<std::string> deprotect_data_abbreviation(
    const DeprotectData &data) {
  return copy(data.abbreviation);
}

std::unique_ptr<std::string> deprotect_data_reaction_smarts(
    const DeprotectData &data) {
  return copy(data.reaction_smarts);
}

std::unique_ptr<std::string> deprotect_data_example(const DeprotectData &data) {
  return copy(data.example);
}

bool deprotect_data_is_valid(const DeprotectData &data) {
  return data.isValid();
}

std::size_t num_default_deprotections() {
  return RDKit::Deprotect::getDeprotections().size();
}

std::unique_ptr<DeprotectData> default_deprotection(std::size_t i) {
  return std::make_unique<DeprotectData>(
      RDKit::Deprotect::getDeprotections().at(i));
}

std::unique_ptr<DeprotectList> deprotect_list_new() {
  return std::make_unique<DeprotectList>();
}

void deprotect_list_push(DeprotectList &list, const DeprotectData &data) {
  list.items.push_back(data);
}

// Deprotection runs reactions, which match substructures (control_c_lock.h).
std::unique_ptr<RWMol> deprotect_with(const RWMol &mol,
                                      const DeprotectList &list) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_rwmol(RDKit::Deprotect::deprotect(mol, list.items));
}

std::unique_ptr<RWMol> deprotect_default(const RWMol &mol) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_rwmol(
      RDKit::Deprotect::deprotect(mol, RDKit::Deprotect::getDeprotections()));
}

// In place, the reactions remove atoms from the molecule (runReactant),
// which trusts what check_removals checks, checked before the first edit.
bool deprotect_in_place_with(RWMol &mol, const DeprotectList &list) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return RDKit::Deprotect::deprotectInPlace(mol, list.items);
}

bool deprotect_in_place_default(RWMol &mol) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return RDKit::Deprotect::deprotectInPlace(
      mol, RDKit::Deprotect::getDeprotections());
}

}  // namespace rdkit_cxx
