/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdFingerprintGenerator: Morgan, atom-pair, RDKit, and topological-torsion
// fingerprint generators and their invariant generators.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <GraphMol/Fingerprints/FingerprintGenerator.h>

#include "rdkit-cxx/data_structs.h"
#include "rdkit-cxx/mol.h"
#include "rdkit-cxx/mol_list.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using FpGenerator = RDKit::FingerprintGenerator<std::uint64_t>;
using FPType = RDKit::FPType;
using AtomInvGen = RDKit::AtomInvariantsGenerator;
using BondInvGen = RDKit::BondInvariantsGenerator;

// Fingerprints a bulk call returns, handed to Rust one at a time.
class FingerprintList {
 public:
  std::vector<std::unique_ptr<BitVect>> bit_vects;
  std::vector<std::unique_ptr<SparseIntVect>> int_vects;
};

// Shared types, defined by the bridge (src/ffi/fingerprint_generator.rs).
enum class FingerprintOutput : std::uint8_t;
struct FingerprintOptionsData;
struct MorganGeneratorData;
struct AtomPairGeneratorData;
struct RDKitFPGeneratorData;
struct TopologicalTorsionGeneratorData;
struct FingerprintArgsData;
struct AdditionalOutputData;

std::unique_ptr<FpGenerator> get_morgan_generator(
    const MorganGeneratorData &data, std::unique_ptr<AtomInvGen> atom_inv_gen,
    std::unique_ptr<BondInvGen> bond_inv_gen);
std::unique_ptr<FpGenerator> get_atom_pair_generator(
    const AtomPairGeneratorData &data,
    std::unique_ptr<AtomInvGen> atom_inv_gen);
std::unique_ptr<FpGenerator> get_rdkit_fp_generator(
    const RDKitFPGeneratorData &data, std::unique_ptr<AtomInvGen> atom_inv_gen);
std::unique_ptr<FpGenerator> get_topological_torsion_generator(
    const TopologicalTorsionGeneratorData &data,
    std::unique_ptr<AtomInvGen> atom_inv_gen);
std::unique_ptr<FpGenerator> fingerprint_generator_from_json(rust::Str json);

std::unique_ptr<AtomInvGen> get_morgan_atom_inv_gen(
    bool include_ring_membership);
std::unique_ptr<AtomInvGen> get_morgan_feature_atom_inv_gen(
    bool has_patterns, BorrowedMols &patterns);
std::unique_ptr<BondInvGen> get_morgan_bond_inv_gen(bool use_bond_types,
                                                    bool use_chirality);
std::unique_ptr<AtomInvGen> get_atom_pair_atom_inv_gen(bool include_chirality);
std::unique_ptr<AtomInvGen> get_rdkit_atom_inv_gen();
std::unique_ptr<AtomInvGen> atom_inv_gen_clone(const AtomInvGen &gen);
std::unique_ptr<BondInvGen> bond_inv_gen_clone(const BondInvGen &gen);
std::unique_ptr<std::string> atom_inv_gen_info_string(const AtomInvGen &gen);
std::unique_ptr<std::string> bond_inv_gen_info_string(const BondInvGen &gen);

FingerprintOptionsData fp_generator_get_options(const FpGenerator &gen);
void fp_generator_set_options(FpGenerator &gen,
                              const FingerprintOptionsData &data);
std::unique_ptr<std::string> fp_generator_info_string(const FpGenerator &gen);
std::unique_ptr<std::string> fp_generator_to_json(const FpGenerator &gen);

// One fingerprint of the given kind; bit vectors come back in `bits`, count
// vectors in `counts`, and the NumPy forms in `values`.
std::unique_ptr<BitVect> fp_generator_bits(const FpGenerator &gen,
                                           const RWMol &mol,
                                           FingerprintOutput output,
                                           const FingerprintArgsData &args,
                                           AdditionalOutputData &outputs);
std::unique_ptr<SparseIntVect> fp_generator_counts(
    const FpGenerator &gen, const RWMol &mol, FingerprintOutput output,
    const FingerprintArgsData &args, AdditionalOutputData &outputs);
rust::Vec<std::uint32_t> fp_generator_values(const FpGenerator &gen,
                                             const RWMol &mol,
                                             FingerprintOutput output,
                                             const FingerprintArgsData &args,
                                             AdditionalOutputData &outputs);

std::unique_ptr<FingerprintList> fp_generator_bulk(const FpGenerator &gen,
                                                   BorrowedMols &mols,
                                                   FingerprintOutput output,
                                                   std::int32_t num_threads);
// As getFPBulk and its siblings: a default generator of the type, one
// molecule at a time.
std::unique_ptr<FingerprintList> fp_type_bulk(BorrowedMols &mols,
                                              FPType fp_type,
                                              FingerprintOutput output);
std::size_t fingerprint_list_size(const FingerprintList &list);
std::unique_ptr<BitVect> fingerprint_list_take_bits(FingerprintList &list,
                                                    std::size_t i);
std::unique_ptr<SparseIntVect> fingerprint_list_take_counts(
    FingerprintList &list, std::size_t i);

}  // namespace rdkit_cxx
