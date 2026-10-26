/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolEnumerator: the molecules a molecule's position variation bonds, link
// nodes, and repeat units stand for.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolEnumerator/MolEnumerator.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared types, defined by the bridge (src/ffi/mol_enumerator.rs).
struct MolEnumeratorParamsData;
enum class EnumeratorType : std::uint8_t;

// A default-constructed MolEnumeratorParams, which has no operation.
MolEnumeratorParamsData mol_enumerator_params_default();
// As MolEnumerator::enumerate, with the bundle's molecules copied out.
std::unique_ptr<MolList> mol_enumerator_enumerate(
    const RWMol &mol, std::uint32_t max_per_operation);
std::unique_ptr<MolList> mol_enumerator_enumerate_with_params(
    const RWMol &mol, const MolEnumeratorParamsData &params);

}  // namespace rdkit_cxx
