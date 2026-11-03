/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/resonance.h"

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"
#include "rdkit-cxx/threads.h"

namespace rdkit_cxx {

namespace {

// RDKit returns new ROMols; Rust molecules own RWMols (W1). Resonance
// structures are copies of the supplier's molecule, which carries the distance
// matrix it computed, so the copies drop it (U6).
std::unique_ptr<RWMol> to_rwmol(RDKit::ROMol *mol) {
  std::unique_ptr<RDKit::ROMol> owned(mol);
  if (!owned) {
    return nullptr;
  }
  auto result = std::make_unique<RWMol>(*owned);
  clear_matrix_caches(*result);
  return result;
}

}  // namespace

// Each conjugated group sums the formal charges of the conjugated atoms in
// int (ConjElectrons in Resonance.cpp), from the supplier's own copy of the
// molecule (U9).
std::unique_ptr<ResonanceMolSupplier> resonance_mol_supplier_new(
    const RWMol &mol, std::uint32_t flags, std::uint32_t max_structs) {
  check_formal_charge_sum(mol);
  // Each resonance structure, which has the molecule's graph with other
  // bond orders, is sanitized for its radicals and hydrogens
  // (ResonanceUtils::sanitizeMol in Resonance.cpp), afresh.
  check_sanitize_relevant_cycles(mol, RDKit::MolOps::SANITIZE_FINDRADICALS |
                                          RDKit::MolOps::SANITIZE_ADJUSTHS);
  // The supplier copies the molecule, so the const_cast does not write it.
  return std::make_unique<ResonanceMolSupplier>(const_cast<RWMol &>(mol), flags,
                                                max_structs);
}

std::uint32_t resonance_mol_supplier_get_num_conj_grps(
    const ResonanceMolSupplier &suppl) {
  return suppl.getNumConjGrps();
}

int resonance_mol_supplier_get_bond_conj_grp_idx(
    const ResonanceMolSupplier &suppl, std::uint32_t bi) {
  check_index(bi, suppl.mol().getNumBonds());
  return suppl.getBondConjGrpIdx(bi);
}

int resonance_mol_supplier_get_atom_conj_grp_idx(
    const ResonanceMolSupplier &suppl, std::uint32_t ai) {
  check_index(ai, suppl.mol().getNumAtoms());
  return suppl.getAtomConjGrpIdx(ai);
}

// setNumThreads counts the threads with getNumThreadsToUse, which negates a
// negative count in int, undefined for INT_MIN; it gets the positive count
// num_threads_to_use computes, at most the hardware's or num_threads.
void resonance_mol_supplier_set_num_threads(ResonanceMolSupplier &suppl,
                                            int num_threads) {
  suppl.setNumThreads(static_cast<int>(num_threads_to_use(num_threads)));
}

bool resonance_mol_supplier_was_canceled(const ResonanceMolSupplier &suppl) {
  return suppl.wasCanceled();
}

void resonance_mol_supplier_enumerate(ResonanceMolSupplier &suppl) {
  suppl.enumerate();
}

bool resonance_mol_supplier_get_is_enumerated(
    const ResonanceMolSupplier &suppl) {
  // getIsEnumerated only reads a flag but is not declared const.
  return const_cast<ResonanceMolSupplier &>(suppl).getIsEnumerated();
}

std::uint32_t resonance_mol_supplier_length(ResonanceMolSupplier &suppl) {
  return suppl.length();
}

void resonance_mol_supplier_reset(ResonanceMolSupplier &suppl) {
  suppl.reset();
}

bool resonance_mol_supplier_at_end(ResonanceMolSupplier &suppl) {
  return suppl.atEnd();
}

std::unique_ptr<RWMol> resonance_mol_supplier_next(
    ResonanceMolSupplier &suppl) {
  return to_rwmol(suppl.next());
}

std::unique_ptr<RWMol> resonance_mol_supplier_get_item(
    ResonanceMolSupplier &suppl, std::uint32_t idx) {
  check_index(idx, suppl.length());
  return to_rwmol(suppl[idx]);
}

}  // namespace rdkit_cxx
