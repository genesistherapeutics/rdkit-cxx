/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdReducedGraphs: extended reduced graphs and ErG fingerprints.

#pragma once

#include <cstdint>
#include <memory>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

std::unique_ptr<RWMol> generate_mol_extended_reduced_graph(const RWMol &mol);
rust::Vec<double> generate_erg_fingerprint_for_reduced_graph(
    const RWMol &mol, double fuzz_increment, std::int32_t min_path,
    std::int32_t max_path);
rust::Vec<double> get_erg_fingerprint(const RWMol &mol, double fuzz_increment,
                                      std::int32_t min_path,
                                      std::int32_t max_path);

}  // namespace rdkit_cxx
