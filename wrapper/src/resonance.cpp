/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/resonance.h"

#include "rdkit-cxx/errors.h"

namespace rdkit_cxx {

namespace {

// RDKit returns new ROMols; Rust molecules own RWMols (W1).
std::unique_ptr<RWMol> to_rwmol(RDKit::ROMol *mol) {
  std::unique_ptr<RDKit::ROMol> owned(mol);
  return owned ? std::make_unique<RWMol>(*owned) : nullptr;
}

}  // namespace

std::unique_ptr<ResonanceMolSupplier> resonance_mol_supplier_new(
    const RWMol &mol, std::uint32_t flags, std::uint32_t max_structs) {
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

void resonance_mol_supplier_set_num_threads(ResonanceMolSupplier &suppl,
                                            int num_threads) {
  suppl.setNumThreads(num_threads);
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
