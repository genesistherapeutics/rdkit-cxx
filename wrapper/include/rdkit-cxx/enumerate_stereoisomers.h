/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdEnumerateStereoisomers: the stereoisomers of a molecule, one at a time.

#pragma once

#include <cstdint>
#include <memory>

#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared type, defined by the bridge (src/ffi/enumerate_stereoisomers.rs).
struct StereoEnumerationOptionsData;

// RDKit's enumerator, which owns a copy of the molecule (its flippers point
// into it) and its own random generator. When its isomers are embedded, the
// wrapper also keeps a copy of the molecule until it has checked what the
// embedder trusts.
class StereoisomerEnumerator {
 public:
  std::unique_ptr<RDKit::EnumerateStereoisomers::StereoisomerEnumerator>
      enumerator;
  bool embeds = false;
  std::unique_ptr<RDKit::ROMol> unchecked;
};

// A default-constructed StereoEnumerationOptions.
StereoEnumerationOptionsData stereo_enumeration_options_default();
std::unique_ptr<StereoisomerEnumerator> stereoisomer_enumerator_new(
    const RWMol &mol, const StereoEnumerationOptionsData &options,
    bool verbose);
std::uint64_t stereoisomer_enumerator_get_stereoisomer_count(
    const StereoisomerEnumerator &enumerator);
// The next isomer, or null after the last.
std::unique_ptr<RWMol> stereoisomer_enumerator_next(
    StereoisomerEnumerator &enumerator);

}  // namespace rdkit_cxx
