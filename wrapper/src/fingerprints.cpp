/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/fingerprints.h"

#include <limits>
#include <map>
#include <vector>

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/fingerprints.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

// As pythonObjectToVect(fromAtoms, mol.getNumAtoms()) in the Python
// wrappers: the atoms paths start from.
std::unique_ptr<std::vector<std::uint32_t>> from_atoms(
    const RWMol &mol, bool has, const rust::Vec<std::uint32_t> &atoms) {
  if (!has) {
    return nullptr;
  }
  auto result = std::make_unique<std::vector<std::uint32_t>>();
  for (auto idx : atoms) {
    check_index(idx, mol.getNumAtoms());
    result->push_back(idx);
  }
  return result;
}

// RDKit's own preconditions check the sizes of these.
template <class T>
std::unique_ptr<std::vector<T>> optional_vector(
    bool has, const rust::Vec<std::uint32_t> &values) {
  if (!has) {
    return nullptr;
  }
  return std::make_unique<std::vector<T>>(values.begin(), values.end());
}

// Fingerprints are ExplicitBitVects, which hold at most INT_MAX bits
// (CONVENTIONS.md U3).
void check_fp_size(std::uint32_t fp_size) {
  if (fp_size > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
    throw ValueErrorException("bit vectors hold at most 2^31 - 1 bits");
  }
}

ExplicitBitVect *only_bits(bool has, const BitVect &bits) {
  // RDKit only reads it.
  return has ? const_cast<ExplicitBitVect *>(&as_explicit(bits)) : nullptr;
}

template <class Bit>
void fill_outputs(const std::vector<std::vector<Bit>> &atom_bits,
                  const std::map<Bit, std::vector<std::vector<int>>> &bit_info,
                  FingerprintOutputs &outputs) {
  for (const auto &bits : atom_bits) {
    U64List list;
    for (auto bit : bits) {
      list.values.push_back(bit);
    }
    outputs.atom_bits.push_back(std::move(list));
  }
  for (const auto &[bit, paths] : bit_info) {
    outputs.bit_info_bits.push_back(bit);
    outputs.bit_info_num_paths.push_back(paths.size());
    for (const auto &path : paths) {
      IndexList list;
      for (auto idx : path) {
        list.indices.push_back(static_cast<std::uint32_t>(idx));
      }
      outputs.bit_info_paths.push_back(std::move(list));
    }
  }
}

}  // namespace

std::unique_ptr<BitVect> rdk_fingerprint(const RWMol &mol,
                                         const RDKFingerprintData &params,
                                         bool with_outputs,
                                         FingerprintOutputs &outputs) {
  check_fp_size(params.fp_size);
  auto invariants = optional_vector<std::uint32_t>(params.has_atom_invariants,
                                                   params.atom_invariants);
  auto from = from_atoms(mol, params.has_from_atoms, params.from_atoms);
  std::vector<std::vector<std::uint32_t>> atom_bits(mol.getNumAtoms());
  std::map<std::uint32_t, std::vector<std::vector<int>>> bit_info;
  std::unique_ptr<ExplicitBitVect> fp(RDKit::RDKFingerprintMol(
      mol, params.min_path, params.max_path, params.fp_size,
      params.n_bits_per_hash, params.use_hs, params.tgt_density,
      params.min_size, params.branched_paths, params.use_bond_order,
      invariants.get(), from.get(), with_outputs ? &atom_bits : nullptr,
      with_outputs ? &bit_info : nullptr));
  if (with_outputs) {
    fill_outputs(atom_bits, bit_info, outputs);
  }
  return wrap_bit_vect(std::move(fp));
}

std::unique_ptr<SparseIntVect> unfolded_rdk_fingerprint_count_based(
    const RWMol &mol, const RDKFingerprintData &params, bool with_outputs,
    FingerprintOutputs &outputs) {
  auto invariants = optional_vector<std::uint32_t>(params.has_atom_invariants,
                                                   params.atom_invariants);
  auto from = from_atoms(mol, params.has_from_atoms, params.from_atoms);
  std::vector<std::vector<std::uint64_t>> atom_bits(mol.getNumAtoms());
  std::map<std::uint64_t, std::vector<std::vector<int>>> bit_info;
  std::unique_ptr<RDKit::SparseIntVect<std::uint64_t>> fp(
      RDKit::getUnfoldedRDKFingerprintMol(
          mol, params.min_path, params.max_path, params.use_hs,
          params.branched_paths, params.use_bond_order, invariants.get(),
          from.get(), with_outputs ? &atom_bits : nullptr,
          with_outputs ? &bit_info : nullptr));
  if (with_outputs) {
    fill_outputs(atom_bits, bit_info, outputs);
  }
  return wrap_sparse_int_vect(std::move(fp));
}

std::unique_ptr<BitVect> layered_fingerprint(
    const RWMol &mol, const LayeredFingerprintData &params,
    const BitVect &set_only_bits, FingerprintOutputs &outputs) {
  check_fp_size(params.fp_size);
  auto counts =
      optional_vector<unsigned int>(params.has_atom_counts, params.atom_counts);
  auto from = from_atoms(mol, params.has_from_atoms, params.from_atoms);
  std::unique_ptr<ExplicitBitVect> fp(RDKit::LayeredFingerprintMol(
      mol, params.layer_flags, params.min_path, params.max_path, params.fp_size,
      counts.get(), only_bits(params.has_set_only_bits, set_only_bits),
      params.branched_paths, from.get()));
  if (counts) {
    for (auto count : *counts) {
      outputs.atom_counts.push_back(count);
    }
  }
  return wrap_bit_vect(std::move(fp));
}

std::unique_ptr<BitVect> pattern_fingerprint(
    const RWMol &mol, const PatternFingerprintData &params,
    const BitVect &set_only_bits, FingerprintOutputs &outputs) {
  check_fp_size(params.fp_size);
  auto counts =
      optional_vector<unsigned int>(params.has_atom_counts, params.atom_counts);
  std::unique_ptr<ExplicitBitVect> fp(RDKit::PatternFingerprintMol(
      mol, params.fp_size, counts.get(),
      only_bits(params.has_set_only_bits, set_only_bits),
      params.tautomer_fingerprints));
  if (counts) {
    for (auto count : *counts) {
      outputs.atom_counts.push_back(count);
    }
  }
  return wrap_bit_vect(std::move(fp));
}

std::unique_ptr<BitVect> pattern_fingerprint_mol_bundle(
    const MolBundle &bundle, const PatternFingerprintData &params,
    const BitVect &set_only_bits) {
  check_fp_size(params.fp_size);
  std::unique_ptr<ExplicitBitVect> fp(RDKit::PatternFingerprintMol(
      bundle, params.fp_size,
      only_bits(params.has_set_only_bits, set_only_bits),
      params.tautomer_fingerprints));
  return wrap_bit_vect(std::move(fp));
}

std::uint32_t layered_fingerprint_substruct_layers() {
  return RDKit::substructLayers;
}

}  // namespace rdkit_cxx
