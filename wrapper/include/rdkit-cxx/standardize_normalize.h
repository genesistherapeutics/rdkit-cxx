/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolStandardize's Normalizer (Normalize.h) and MetalDisconnector
// (Metal.h).

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/MolStandardize/Metal.h>
#include <GraphMol/MolStandardize/Normalize.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/standardize.rs and
// src/ffi/standardize_normalize.rs).
struct CleanupParametersData;
struct MetalDisconnectorOptionsData;

using Normalizer = RDKit::MolStandardize::Normalizer;
using MetalDisconnector = RDKit::MolStandardize::MetalDisconnector;

std::unique_ptr<Normalizer> normalizer_default();
std::unique_ptr<Normalizer> normalizer_new(
    rust::Slice<const std::uint8_t> normalize_filename,
    std::uint32_t max_restarts);
std::unique_ptr<Normalizer> normalizer_from_data(
    rust::Str param_data, const CleanupParametersData &params);
std::unique_ptr<Normalizer> normalizer_from_params(
    const CleanupParametersData &params);
std::unique_ptr<RWMol> normalizer_normalize(const Normalizer &normalizer,
                                            const RWMol &mol);
void normalizer_normalize_in_place(const Normalizer &normalizer, RWMol &mol);

MetalDisconnectorOptionsData metal_disconnector_options_default();
std::unique_ptr<MetalDisconnector> metal_disconnector_new(
    const MetalDisconnectorOptionsData &options);
std::unique_ptr<std::string> metal_disconnector_get_metal_nof(
    const MetalDisconnector &disconnector);
std::unique_ptr<std::string> metal_disconnector_get_metal_non(
    const MetalDisconnector &disconnector);
// The setters check the query's shape (see check_metal_query).
void metal_disconnector_set_metal_nof(MetalDisconnector &disconnector,
                                      const RWMol &mol);
void metal_disconnector_set_metal_non(MetalDisconnector &disconnector,
                                      const RWMol &mol);
std::unique_ptr<RWMol> metal_disconnector_disconnect(
    const MetalDisconnector &disconnector, const RWMol &mol);
void metal_disconnector_disconnect_in_place(
    const MetalDisconnector &disconnector, RWMol &mol);

// DisconnectOrganometallics without options uses the C++ defaults of
// disconnectOrganometallics, not MetalDisconnectorOptions'.
std::unique_ptr<RWMol> disconnect_organometallics(const RWMol &mol);
std::unique_ptr<RWMol> disconnect_organometallics_with_options(
    const RWMol &mol, const MetalDisconnectorOptionsData &options);
void disconnect_organometallics_in_place(RWMol &mol);
void disconnect_organometallics_in_place_with_options(
    RWMol &mol, const MetalDisconnectorOptionsData &options);

// MetalDisconnector::disconnect takes each match of its queries for a
// metal (the match's first atom) and the atom bonded to it (the second): it
// reads match[1] and dereferences the bond between them, which with fewer
// than two query atoms is read past the match's end, and with atoms that
// are not bonded, or a pair matched twice (a third atom matched two ways,
// after the first match removed the bond), is null. Throws
// ValueErrorException unless the query has no atoms, or two joined by a
// bond: SubstructMatch then finds each bonded pair once.
void check_metal_query(const RDKit::ROMol &query);

}  // namespace rdkit_cxx
