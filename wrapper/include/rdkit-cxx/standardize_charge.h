/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolStandardize's ChargeCorrection, Reionizer, and Uncharger (Charge.h).

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolStandardize/Charge.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// A shared struct, defined by the bridge (src/ffi/standardize_charge.rs).
struct ChargeCorrectionData;

using Reionizer = RDKit::MolStandardize::Reionizer;
using Uncharger = RDKit::MolStandardize::Uncharger;

rust::Vec<ChargeCorrectionData> charge_corrections();

// The constructors check the charge corrections' SMARTS (Charge.cpp parses
// them on every reionization and dereferences what SmartsToMol returns,
// null for SMARTS that do not parse).
std::unique_ptr<Reionizer> reionizer_default();
std::unique_ptr<Reionizer> reionizer_new(
    rust::Slice<const std::uint8_t> acidbase_file,
    rust::Slice<const ChargeCorrectionData> ccs, bool default_ccs);
std::unique_ptr<Reionizer> reionizer_from_data(
    rust::Str param_data, rust::Slice<const ChargeCorrectionData> ccs);
std::unique_ptr<RWMol> reionizer_reionize(const Reionizer &reionizer,
                                          const RWMol &mol);
void reionizer_reionize_in_place(const Reionizer &reionizer, RWMol &mol);

std::unique_ptr<Uncharger> uncharger_new(bool canonical_order, bool force,
                                         bool protonation_only);
std::unique_ptr<RWMol> uncharger_uncharge(const Uncharger &uncharger,
                                          const RWMol &mol);
void uncharger_uncharge_in_place(const Uncharger &uncharger, RWMol &mol);

}  // namespace rdkit_cxx
