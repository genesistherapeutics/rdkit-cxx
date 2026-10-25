/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMHFPFingerprint: MinHash (MHFP) and SECFP fingerprints.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>

#include <GraphMol/Fingerprints/MHFP.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using MHFPEncoder = RDKit::MHFPFingerprints::MHFPEncoder;

// The shingling options every encoder method takes.
struct ShinglingOptions;

std::unique_ptr<MHFPEncoder> mhfp_encoder_new(std::uint32_t n_permutations,
                                              std::uint32_t seed);
rust::Vec<std::uint32_t> mhfp_from_string_array(
    const MHFPEncoder &encoder, rust::Slice<const rust::Str> vec);
rust::Vec<std::uint32_t> mhfp_from_array(const MHFPEncoder &encoder,
                                         rust::Slice<const std::uint32_t> vec);
rust::Vec<rust::String> mhfp_shingling_from_smiles(
    const MHFPEncoder &encoder, rust::Str smiles,
    const ShinglingOptions &options);
rust::Vec<rust::String> mhfp_shingling_from_mol(
    const MHFPEncoder &encoder, const RWMol &mol,
    const ShinglingOptions &options);
rust::Vec<std::uint32_t> mhfp_encode_smiles(const MHFPEncoder &encoder,
                                            rust::Str smiles,
                                            const ShinglingOptions &options);
rust::Vec<std::uint32_t> mhfp_encode_mol(const MHFPEncoder &encoder,
                                         const RWMol &mol,
                                         const ShinglingOptions &options);
std::unique_ptr<BitVect> mhfp_encode_secfp_smiles(
    const MHFPEncoder &encoder, rust::Str smiles,
    const ShinglingOptions &options, std::size_t length);
std::unique_ptr<BitVect> mhfp_encode_secfp_mol(const MHFPEncoder &encoder,
                                               const RWMol &mol,
                                               const ShinglingOptions &options,
                                               std::size_t length);
double mhfp_distance(rust::Slice<const std::uint32_t> a,
                     rust::Slice<const std::uint32_t> b);

}  // namespace rdkit_cxx
