/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Resonance structure enumeration. The substructure searches are in
// rdkit-cxx/substruct.h.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/Resonance.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using ResonanceMolSupplier = RDKit::ResonanceMolSupplier;
using ResonanceFlags = RDKit::ResonanceMolSupplier::ResonanceFlags;

std::unique_ptr<ResonanceMolSupplier> resonance_mol_supplier_new(
    const RWMol &mol, std::uint32_t flags, std::uint32_t max_structs);
std::uint32_t resonance_mol_supplier_get_num_conj_grps(
    const ResonanceMolSupplier &suppl);
// The atom or bond's conjugated group, or -1; indices are checked against the
// molecule (CONVENTIONS.md E3).
int resonance_mol_supplier_get_bond_conj_grp_idx(
    const ResonanceMolSupplier &suppl, std::uint32_t bi);
int resonance_mol_supplier_get_atom_conj_grp_idx(
    const ResonanceMolSupplier &suppl, std::uint32_t ai);
void resonance_mol_supplier_set_num_threads(ResonanceMolSupplier &suppl,
                                            int num_threads);
bool resonance_mol_supplier_was_canceled(const ResonanceMolSupplier &suppl);
void resonance_mol_supplier_enumerate(ResonanceMolSupplier &suppl);
bool resonance_mol_supplier_get_is_enumerated(
    const ResonanceMolSupplier &suppl);
// These enumerate the structures first, as RDKit's do.
std::uint32_t resonance_mol_supplier_length(ResonanceMolSupplier &suppl);
void resonance_mol_supplier_reset(ResonanceMolSupplier &suppl);
bool resonance_mol_supplier_at_end(ResonanceMolSupplier &suppl);
// The next structure, or null at the end.
std::unique_ptr<RWMol> resonance_mol_supplier_next(ResonanceMolSupplier &suppl);
std::unique_ptr<RWMol> resonance_mol_supplier_get_item(
    ResonanceMolSupplier &suppl, std::uint32_t idx);

}  // namespace rdkit_cxx
