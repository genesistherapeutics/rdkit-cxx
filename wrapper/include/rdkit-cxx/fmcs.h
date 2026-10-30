/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdFMCS: maximum common substructures, with the callbacks Rust provides
// for comparing atoms and bonds, tracking progress, and checking matches.

#pragma once

#include <cstddef>
#include <cstdint>

#include <GraphMol/FMCS/FMCS.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using AtomComparator = RDKit::AtomComparator;
using BondComparator = RDKit::BondComparator;

// Shared structs, defined by the bridge (src/ffi/fmcs.rs).
struct MCSAtomCompareParametersData;
struct MCSBondCompareParametersData;
struct MCSParametersData;
struct MCSResultData;
// The Rust side of one FindMCS call: the molecules and the closures the
// trampolines in fmcs.cpp call (src/ffi/fmcs.rs).
struct McsCallbacks;

MCSAtomCompareParametersData mcs_atom_compare_parameters_default();
MCSBondCompareParametersData mcs_bond_compare_parameters_default();
MCSParametersData mcs_parameters_default();

MCSResultData find_mcs(BorrowedMols &mols, const MCSParametersData &params,
                       const McsCallbacks &callbacks);

// The predefined checks Python's MCSAtomCompare and MCSBondCompare offer
// their subclasses.
bool mcs_atom_compare_check_atom_ring_match(
    const MCSAtomCompareParametersData &params, const RWMol &mol1,
    std::uint32_t atom1, const RWMol &mol2, std::uint32_t atom2);
bool mcs_atom_compare_check_atom_charge(
    const MCSAtomCompareParametersData &params, const RWMol &mol1,
    std::uint32_t atom1, const RWMol &mol2, std::uint32_t atom2);
bool mcs_atom_compare_check_atom_chirality(
    const MCSAtomCompareParametersData &params, const RWMol &mol1,
    std::uint32_t atom1, const RWMol &mol2, std::uint32_t atom2);
bool mcs_bond_compare_check_bond_stereo(
    const MCSBondCompareParametersData &params, const RWMol &mol1,
    std::uint32_t bond1, const RWMol &mol2, std::uint32_t bond2);
bool mcs_bond_compare_check_bond_ring_match(
    const MCSBondCompareParametersData &params, const RWMol &mol1,
    std::uint32_t bond1, const RWMol &mol2, std::uint32_t bond2);

}  // namespace rdkit_cxx
