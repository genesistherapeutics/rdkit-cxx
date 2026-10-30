/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Molecule bundles. The Rust MolBundle owns its molecules; each call builds a
// C++ bundle that borrows them for the length of the call
// (src/chem/rdchem/mol_bundle.rs).

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/MolBundle.h>

#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using MolBundle = RDKit::MolBundle;

std::unique_ptr<MolBundle> mol_bundle_new();
// Adds `mol` without taking ownership: the bundle must not outlive it.
void mol_bundle_add_borrowed(MolBundle &bundle, const RWMol &mol);
// Throws as FixedMolSizeMolBundle::addMol does when `nmol` cannot join a
// bundle whose first molecule is `first`.
void fixed_mol_size_mol_bundle_check(const RWMol &first, const RWMol &nmol);
// As MolBundle(pklString) and BundleToBinary in
// Code/GraphMol/Wrap/MolBundle.cpp, except that without serialization support
// BundleToBinary raises RuntimeError there and ValueError here.
std::unique_ptr<MolList> mol_bundle_from_pkl_string(
    rust::Slice<const std::uint8_t> pkl);
rust::Vec<std::uint8_t> mol_bundle_to_binary(const MolBundle &bundle);
bool mol_bundle_can_serialize();

}  // namespace rdkit_cxx
