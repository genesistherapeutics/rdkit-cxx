// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize_charge.h"

#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/standardize_charge.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/standardize.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace MS = RDKit::MolStandardize;

namespace {

std::string text(const rust::Vec<std::uint8_t> &bytes) {
  return to_std(rust::Slice<const std::uint8_t>(bytes.data(), bytes.size()));
}

// Reionizer::reionizeInPlace parses each correction's SMARTS and matches
// the result without checking it, so SMARTS that do not parse would be
// dereferenced as null.
std::vector<MS::ChargeCorrection> checked(
    rust::Slice<const ChargeCorrectionData> ccs) {
  std::vector<MS::ChargeCorrection> result;
  for (const auto &cc : ccs) {
    std::string smarts = text(cc.smarts);
    std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
    if (!query) {
      throw ValueErrorException("charge correction SMARTS do not parse: " +
                                smarts);
    }
    result.emplace_back(text(cc.name), smarts, cc.charge);
  }
  return result;
}

// The objects' methods do not change them, but are not const.
template <class T>
T &mut(const T &object) {
  return const_cast<T &>(object);
}

}  // namespace

rust::Vec<ChargeCorrectionData> charge_corrections() {
  rust::Vec<ChargeCorrectionData> result;
  for (const auto &cc : MS::CHARGE_CORRECTIONS) {
    result.push_back(ChargeCorrectionData{to_bytes(cc.Name),
                                          to_bytes(cc.Smarts), cc.Charge});
  }
  return result;
}

std::unique_ptr<Reionizer> reionizer_default() {
  return std::make_unique<Reionizer>();
}

std::unique_ptr<Reionizer> reionizer_new(
    rust::Slice<const std::uint8_t> acidbase_file,
    rust::Slice<const ChargeCorrectionData> ccs, bool default_ccs) {
  auto name = to_path(acidbase_file);
  if (default_ccs) {
    return std::make_unique<Reionizer>(name);
  }
  return std::make_unique<Reionizer>(name, checked(ccs));
}

std::unique_ptr<Reionizer> reionizer_from_data(
    rust::Str param_data, rust::Slice<const ChargeCorrectionData> ccs) {
  auto corrections = checked(ccs);
  std::istringstream in{std::string(param_data)};
  return std::make_unique<Reionizer>(in, corrections);
}

std::unique_ptr<RWMol> reionizer_reionize(const Reionizer &reionizer,
                                          const RWMol &mol) {
  check_reionizer_atoms(mol);
  // The reionizer and uncharger match substructures (control_c_lock.h).
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return adopt_mol(mut(reionizer).reionize(mol));
}

void reionizer_reionize_in_place(const Reionizer &reionizer, RWMol &mol) {
  check_reionizer_atoms(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  mut(reionizer).reionizeInPlace(mol);
}

std::unique_ptr<Uncharger> uncharger_new(bool canonical_order, bool force,
                                         bool protonation_only) {
  return std::make_unique<Uncharger>(canonical_order, force, protonation_only);
}

// The uncharger sums the molecule's formal charges in int after it
// neutralizes some, which only shrinks their magnitudes (U9).
std::unique_ptr<RWMol> uncharger_uncharge(const Uncharger &uncharger,
                                          const RWMol &mol) {
  check_formal_charge_sum(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  return adopt_mol(mut(uncharger).uncharge(mol));
}

void uncharger_uncharge_in_place(const Uncharger &uncharger, RWMol &mol) {
  check_formal_charge_sum(mol);
  std::lock_guard<std::mutex> control_c(control_c_mutex());
  mut(uncharger).unchargeInPlace(mol);
}

}  // namespace rdkit_cxx
