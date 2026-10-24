/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit's periodic table, a process-wide singleton. The functions take an
// atomic number; Rust resolves element symbols with
// periodic_table_get_atomic_number first.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "rust/cxx.h"

namespace rdkit_cxx {

unsigned periodic_table_get_atomic_number(rust::Str element_symbol);
unsigned periodic_table_get_max_atomic_number();
std::unique_ptr<std::string> periodic_table_get_element_symbol(
    unsigned atomic_number);
std::unique_ptr<std::string> periodic_table_get_element_name(
    unsigned atomic_number);
double periodic_table_get_atomic_weight(unsigned atomic_number);
double periodic_table_get_rb0(unsigned atomic_number);
double periodic_table_get_rcovalent(unsigned atomic_number);
double periodic_table_get_rvdw(unsigned atomic_number);
int periodic_table_get_default_valence(unsigned atomic_number);
rust::Vec<std::int32_t> periodic_table_get_valence_list(unsigned atomic_number);
int periodic_table_get_n_outer_elecs(unsigned atomic_number);
int periodic_table_get_most_common_isotope(unsigned atomic_number);
double periodic_table_get_most_common_isotope_mass(unsigned atomic_number);
unsigned periodic_table_get_row(unsigned atomic_number);
double periodic_table_get_mass_for_isotope(unsigned atomic_number,
                                           unsigned isotope);
double periodic_table_get_abundance_for_isotope(unsigned atomic_number,
                                                unsigned isotope);

}  // namespace rdkit_cxx
