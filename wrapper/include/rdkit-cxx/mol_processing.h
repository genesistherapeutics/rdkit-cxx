/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolProcessing: fingerprints for the molecules in a file.
//
// RDKit's getFingerprintsForMolsInFile (MolProcessing.cpp) opens the file
// itself (GeneralMolSupplier::getSupplier), so its input cannot be checked as
// it is read (CONVENTIONS.md U12), it reads MAE files with the
// MaeMolSupplier the crate does not bind, and its threads would match
// substructures under no lock (U4). These wrappers do what it does, with
// RDKit's own suppliers and generators, over a checked stream.

#pragma once

#include <cstdint>
#include <memory>

#include "rdkit-cxx/fingerprint_generator.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared struct, defined by the bridge (src/ffi/mol_processing.rs).
struct SupplierOptionsData;

// A default-constructed GeneralMolSupplier::SupplierOptions.
SupplierOptionsData supplier_options_default();

// As getFingerprintsForMolsInFile with no generator: Morgan fingerprints of
// radius 3 from a 32-bit generator, as Python's first overload makes.
std::unique_ptr<FingerprintList> get_fingerprints_for_mols_in_file(
    rust::Slice<const std::uint8_t> path, const SupplierOptionsData &options);
// As getFingerprintsForMolsInFile with a (64-bit) generator.
std::unique_ptr<FingerprintList>
get_fingerprints_for_mols_in_file_with_generator(
    rust::Slice<const std::uint8_t> path, const FpGenerator &generator,
    const SupplierOptionsData &options);

}  // namespace rdkit_cxx
