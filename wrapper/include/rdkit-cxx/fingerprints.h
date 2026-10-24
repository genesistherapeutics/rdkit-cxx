/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdmolops's path-based fingerprints: RDKFingerprint, its unfolded count
// form, LayeredFingerprint, and PatternFingerprint.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/Fingerprints/Fingerprints.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_bundle.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridge (src/ffi/fingerprints.rs).
struct RDKFingerprintData;
struct LayeredFingerprintData;
struct PatternFingerprintData;
struct FingerprintOutputs;

// The outputs (Python's atomBits and bitInfo) are filled only when asked
// for, as Python fills them only when given.
std::unique_ptr<BitVect> rdk_fingerprint(const RWMol &mol,
                                         const RDKFingerprintData &params,
                                         bool with_outputs,
                                         FingerprintOutputs &outputs);
std::unique_ptr<SparseIntVect> unfolded_rdk_fingerprint_count_based(
    const RWMol &mol, const RDKFingerprintData &params, bool with_outputs,
    FingerprintOutputs &outputs);
std::unique_ptr<BitVect> layered_fingerprint(
    const RWMol &mol, const LayeredFingerprintData &params,
    const BitVect &set_only_bits, FingerprintOutputs &outputs);
std::unique_ptr<BitVect> pattern_fingerprint(
    const RWMol &mol, const PatternFingerprintData &params,
    const BitVect &set_only_bits, FingerprintOutputs &outputs);
std::unique_ptr<BitVect> pattern_fingerprint_mol_bundle(
    const MolBundle &bundle, const PatternFingerprintData &params,
    const BitVect &set_only_bits);
std::uint32_t layered_fingerprint_substruct_layers();

}  // namespace rdkit_cxx
