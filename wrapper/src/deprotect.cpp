// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/deprotect.h"

#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/MolOps.h>
#include <RDGeneral/RDLog.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"

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

// RDKit's deprotections run each reaction on the molecule, or on the
// product of the last, and sanitize the products they keep, up to 100
// times; each product's rings are new. These are Deprotect::deprotect and
// deprotectInPlace (Deprotect.cpp), with each run's products, and each
// molecule before it is sanitized, checked first (relevant_cycles.h).
constexpr unsigned int kMaxDeprotections = 100;

std::unique_ptr<RDKit::ROMol> checked_deprotect(
    const RDKit::ROMol &mol,
    const std::vector<RDKit::Deprotect::DeprotectData> &deprotections) {
  using RDKit::MolOps::SanitizeFlags;
  std::vector<std::string> deprotections_used;
  bool something_happened = true;
  auto m = boost::make_shared<RDKit::ROMol>(mol);
  while (something_happened) {
    something_happened = false;
    for (const auto &deprotect : deprotections) {
      if (!deprotect.isValid()) {
        continue;
      }
      check_products_relevant_cycles(*deprotect.rxn, {m}, 0, true, 0);
      for (auto &prods : deprotect.rxn->runReactant(m, 0)) {
        m = prods[0];
        try {
          RDLog::LogStateSetter blocker;
          RDKit::MolOps::sanitizeMol(*dynamic_cast<RDKit::RWMol *>(m.get()));
        } catch (RDKit::MolSanitizeException &) {
          continue;
        }
        deprotections_used.push_back(deprotect.abbreviation);
        if (deprotections_used.size() >= kMaxDeprotections) {
          BOOST_LOG(rdErrorLog)
              << "Too many deprotections, halting..." << std::endl;
        } else {
          something_happened = true;
        }
        break;
      }
    }
  }
  m->setProp("DEPROTECTIONS", deprotections_used);
  m->setProp<int>("DEPROTECTION_COUNT", deprotections_used.size());
  return std::make_unique<RDKit::ROMol>(*m);
}

// The in-place runs change the molecule before it is checked, so they run on
// a copy, which replaces the molecule once they have all run.
bool checked_deprotect_in_place(
    RWMol &mol,
    const std::vector<RDKit::Deprotect::DeprotectData> &deprotections) {
  using RDKit::MolOps::SanitizeFlags;
  RWMol work(mol);
  std::vector<std::string> deprotections_used;
  bool modified = false;
  bool something_happened = true;
  while (something_happened) {
    something_happened = false;
    for (const auto &deprotect : deprotections) {
      if (!deprotect.isValid()) {
        continue;
      }
      bool changes = deprotect.rxn->runReactant(work);
      if (changes) {
        check_sanitize_relevant_cycles(work, SanitizeFlags::SANITIZE_ALL);
        try {
          RDLog::LogStateSetter blocker;
          RDKit::MolOps::sanitizeMol(work);
        } catch (RDKit::MolSanitizeException &) {
          continue;
        }
        modified = true;
        deprotections_used.push_back(deprotect.abbreviation);
        if (deprotections_used.size() >= kMaxDeprotections) {
          BOOST_LOG(rdErrorLog)
              << "Too many deprotections, halting..." << std::endl;
        } else {
          something_happened = true;
        }
      }
    }
  }
  work.setProp("DEPROTECTIONS", deprotections_used);
  work.setProp<int>("DEPROTECTION_COUNT", deprotections_used.size());
  mol = std::move(work);
  return modified;
}

}  // namespace

std::unique_ptr<DeprotectData> deprotect_data_new(rust::Str deprotection_class,
                                                  rust::Str reaction_smarts,
                                                  rust::Str abbreviation,
                                                  rust::Str full_name) {
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
  return to_rwmol(checked_deprotect(mol, list.items));
}

std::unique_ptr<RWMol> deprotect_default(const RWMol &mol) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_rwmol(checked_deprotect(mol, RDKit::Deprotect::getDeprotections()));
}

// In place, the reactions remove atoms from the molecule (runReactant),
// which trusts what check_removals checks, checked before the first edit.
bool deprotect_in_place_with(RWMol &mol, const DeprotectList &list) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return checked_deprotect_in_place(mol, list.items);
}

bool deprotect_in_place_default(RWMol &mol) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return checked_deprotect_in_place(mol, RDKit::Deprotect::getDeprotections());
}

}  // namespace rdkit_cxx
