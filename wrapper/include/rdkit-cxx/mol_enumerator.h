/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolEnumerator: the molecules a molecule's position variation bonds, link
// nodes, and repeat units stand for.

#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <GraphMol/MolEnumerator/MolEnumerator.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared types, defined by the bridge (src/ffi/mol_enumerator.rs).
struct MolEnumeratorParamsData;
enum class EnumeratorType : std::uint8_t;

// Throws ValueErrorException for what MolEnumerator::enumerate would trust
// in an enumeration with these operations: more variation points than its
// recursion can take, link nodes its reactions read past an empty vector on,
// and the removals check_removals checks (mol_enumerator.cpp says why).
// Every wrapper whose RDKit call enumerates calls it first, with the
// operations RDKit runs; MolEnumerator::enumerate(mol) runs all three.
void check_enumerate(const RDKit::ROMol &mol, bool position_variation,
                     bool repeat_unit, bool link_node);

// Throws ValueErrorException if MolEnumerator::enumerate(mol, ops) would
// symmetrize the rings of a link-node product with more than INT_MAX
// relevant cycles (relevant_cycles.h; mol_enumerator.cpp says when). Called
// under control_c_mutex, as the enumeration is.
void check_enumerate_relevant_cycles(
    const RDKit::ROMol &mol,
    const std::vector<RDKit::MolEnumerator::MolEnumeratorParams> &ops);
// The operations MolEnumerator::enumerate(mol, max_per_operation) runs.
std::vector<RDKit::MolEnumerator::MolEnumeratorParams> default_enumerate_ops(
    std::size_t max_per_operation);

// A default-constructed MolEnumeratorParams, which has no operation.
MolEnumeratorParamsData mol_enumerator_params_default();
// As MolEnumerator::enumerate, with the bundle's molecules copied out.
std::unique_ptr<MolList> mol_enumerator_enumerate(
    const RWMol &mol, std::uint32_t max_per_operation);
std::unique_ptr<MolList> mol_enumerator_enumerate_with_params(
    const RWMol &mol, const MolEnumeratorParamsData &params);

}  // namespace rdkit_cxx
