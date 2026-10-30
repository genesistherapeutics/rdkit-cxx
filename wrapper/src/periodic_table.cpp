/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/periodic_table.h"

#include <GraphMol/PeriodicTable.h>

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

namespace {

const RDKit::PeriodicTable &table() {
  return *RDKit::PeriodicTable::getTable();
}

// RDKit 2024.03's table holds atomic numbers 0 to 118 (atomic_data.cpp),
// and has no getter for its size.
constexpr unsigned kMaxAtomicNumber = 118;

// RDKit checks atomic numbers with PRECONDITION; this reports them as
// IndexError, like other indices (CONVENTIONS.md E3).
unsigned checked(unsigned atomic_number) {
  check_index(atomic_number, kMaxAtomicNumber + 1);
  return atomic_number;
}

}  // namespace

unsigned periodic_table_get_atomic_number(rust::Str element_symbol) {
  return table().getAtomicNumber(std::string(element_symbol));
}

std::unique_ptr<std::string> periodic_table_get_element_symbol(
    unsigned atomic_number) {
  return std::make_unique<std::string>(
      table().getElementSymbol(checked(atomic_number)));
}

std::unique_ptr<std::string> periodic_table_get_element_name(
    unsigned atomic_number) {
  return std::make_unique<std::string>(
      table().getElementName(checked(atomic_number)));
}

double periodic_table_get_atomic_weight(unsigned atomic_number) {
  return table().getAtomicWeight(checked(atomic_number));
}

double periodic_table_get_rb0(unsigned atomic_number) {
  return table().getRb0(checked(atomic_number));
}

double periodic_table_get_rcovalent(unsigned atomic_number) {
  return table().getRcovalent(checked(atomic_number));
}

double periodic_table_get_rvdw(unsigned atomic_number) {
  return table().getRvdw(checked(atomic_number));
}

int periodic_table_get_default_valence(unsigned atomic_number) {
  return table().getDefaultValence(checked(atomic_number));
}

rust::Vec<std::int32_t> periodic_table_get_valence_list(
    unsigned atomic_number) {
  rust::Vec<std::int32_t> result;
  for (auto v : table().getValenceList(checked(atomic_number))) {
    result.push_back(v);
  }
  return result;
}

int periodic_table_get_n_outer_elecs(unsigned atomic_number) {
  return table().getNouterElecs(checked(atomic_number));
}

int periodic_table_get_most_common_isotope(unsigned atomic_number) {
  return table().getMostCommonIsotope(checked(atomic_number));
}

double periodic_table_get_most_common_isotope_mass(unsigned atomic_number) {
  return table().getMostCommonIsotopeMass(checked(atomic_number));
}

double periodic_table_get_mass_for_isotope(unsigned atomic_number,
                                           unsigned isotope) {
  return table().getMassForIsotope(checked(atomic_number), isotope);
}

double periodic_table_get_abundance_for_isotope(unsigned atomic_number,
                                                unsigned isotope) {
  return table().getAbundanceForIsotope(checked(atomic_number), isotope);
}

}  // namespace rdkit_cxx
