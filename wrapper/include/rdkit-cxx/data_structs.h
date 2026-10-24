/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Bit vectors, sparse integer vectors, and value vectors. Types that differ
// only by a template argument or share an interface sit behind one variant
// class each, so one set of wrappers serves every Rust type that mirrors them
// (src/data_structs/c_data_structs/).

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <variant>

#include <DataStructs/BitVects.h>
#include <DataStructs/DiscreteValueVect.h>
#include <DataStructs/RealValueVect.h>
#include <DataStructs/SparseIntVect.h>

#include "rust/cxx.h"

namespace rdkit_cxx {

// An ExplicitBitVect or a SparseBitVect.
class BitVect {
 public:
  explicit BitVect(ExplicitBitVect bv) : v(std::move(bv)) {}
  explicit BitVect(SparseBitVect bv) : v(std::move(bv)) {}
  std::variant<ExplicitBitVect, SparseBitVect> v;
};

// A SparseIntVect by index type: int32, int64, uint32, uint64, in the order
// of SparseIntVectKind in src/ffi/data_structs.rs. Indices cross as uint64;
// none is negative.
class SparseIntVect {
 public:
  template <class IndexType>
  explicit SparseIntVect(RDKit::SparseIntVect<IndexType> siv)
      : v(std::move(siv)) {}
  std::variant<
      RDKit::SparseIntVect<std::int32_t>, RDKit::SparseIntVect<std::int64_t>,
      RDKit::SparseIntVect<std::uint32_t>, RDKit::SparseIntVect<std::uint64_t>>
      v;
};

using DiscreteValueVect = RDKit::DiscreteValueVect;
using DiscreteValueType = RDKit::DiscreteValueVect::DiscreteValueType;
using RealValueVect = RDKit::RealValueVect;

// Shared structs, defined by the bridge (src/ffi/data_structs.rs).
struct SparseIntVectElements;

// For other areas' wrappers: the RDKit objects behind the Rust types, and
// wrappers for the objects RDKit returns. A mismatched kind is a bug in the
// bindings and throws std::bad_variant_access.
const ExplicitBitVect &as_explicit(const BitVect &bv);
std::unique_ptr<BitVect> wrap_bit_vect(std::unique_ptr<ExplicitBitVect> bv);
std::unique_ptr<BitVect> wrap_bit_vect(std::unique_ptr<SparseBitVect> bv);
template <class IndexType>
const RDKit::SparseIntVect<IndexType> &as_sparse_int_vect(
    const SparseIntVect &siv) {
  return std::get<RDKit::SparseIntVect<IndexType>>(siv.v);
}
template <class IndexType>
std::unique_ptr<SparseIntVect> wrap_sparse_int_vect(
    std::unique_ptr<RDKit::SparseIntVect<IndexType>> siv) {
  return std::make_unique<SparseIntVect>(std::move(*siv));
}

// Bit vectors. Indices are checked (CONVENTIONS.md E3), and so are the sizes
// for &, |, and ^, which RDKit leaves to boost::dynamic_bitset's assert.
std::unique_ptr<BitVect> bit_vect_new_explicit(std::uint32_t size,
                                               bool bits_set);
std::unique_ptr<BitVect> bit_vect_new_sparse(std::uint32_t size);
std::unique_ptr<BitVect> bit_vect_from_pkl(bool sparse,
                                           rust::Slice<const std::uint8_t> pkl);
std::unique_ptr<BitVect> bit_vect_clone(const BitVect &bv);
bool bit_vect_set_bit(BitVect &bv, std::uint32_t which);
bool bit_vect_un_set_bit(BitVect &bv, std::uint32_t which);
void bit_vect_set_bits_from_list(BitVect &bv,
                                 rust::Slice<const std::uint32_t> on_bit_list);
void bit_vect_un_set_bits_from_list(
    BitVect &bv, rust::Slice<const std::uint32_t> off_bit_list);
bool bit_vect_get_bit(const BitVect &bv, std::uint32_t which);
std::uint32_t bit_vect_get_num_bits(const BitVect &bv);
std::uint32_t bit_vect_get_num_on_bits(const BitVect &bv);
std::uint32_t bit_vect_get_num_off_bits(const BitVect &bv);
rust::Vec<std::uint32_t> bit_vect_get_on_bits(const BitVect &bv);
rust::Vec<std::uint8_t> bit_vect_to_binary(const BitVect &bv);
void bit_vect_from_base64(BitVect &bv, rust::Str in_d);
std::unique_ptr<std::string> bit_vect_to_base64(const BitVect &bv);
rust::Vec<bool> bit_vect_to_list(const BitVect &bv);
// `op` is '&', '|', '^', or '+' (explicit only).
std::unique_ptr<BitVect> bit_vect_binary_op(const BitVect &a, const BitVect &b,
                                            std::uint8_t op);
void bit_vect_append(BitVect &bv, const BitVect &other);
std::unique_ptr<BitVect> bit_vect_invert(const BitVect &bv);
bool bit_vect_eq(const BitVect &a, const BitVect &b);

// As the metrics in DataStructs/BitOps.h, by the Metric numbers in
// src/data_structs/c_data_structs/similarity.rs; `a` and `b` are Tversky's.
double bit_vect_similarity(const BitVect &a, const BitVect &b,
                           std::uint8_t metric, double alpha, double beta,
                           bool return_distance);
double bit_vect_similarity_pkl(const BitVect &a,
                               rust::Slice<const std::uint8_t> pkl,
                               std::uint8_t metric, double alpha, double beta,
                               bool return_distance);
rust::Vec<double> bit_vect_proj_similarity(const BitVect &a, const BitVect &b,
                                           bool on_bits);
std::int32_t bit_vect_num_bits_in_common(const BitVect &a, const BitVect &b);
rust::Vec<std::uint32_t> bit_vect_bits_in_common(const BitVect &a,
                                                 const BitVect &b,
                                                 bool on_bits);
std::unique_ptr<BitVect> bit_vect_fold(const BitVect &bv,
                                       std::uint32_t fold_factor);
bool bit_vect_all_probe_bits_match(const BitVect &probe, const BitVect &ref);
bool bit_vect_all_probe_bits_match_pkl(const BitVect &probe,
                                       rust::Slice<const std::uint8_t> ref);
std::unique_ptr<std::string> bit_vect_to_text(const BitVect &bv);
std::unique_ptr<std::string> bit_vect_to_fps_text(const BitVect &bv);
rust::Vec<std::uint8_t> bit_vect_to_binary_text(const BitVect &bv);
void bit_vect_init_from_daylight_string(BitVect &bv, rust::Str s);
std::unique_ptr<BitVect> convert_to_explicit(const BitVect &sbv);
// As createFromBitString, createFromFPSText, and createFromBinaryText in
// Code/DataStructs/nbWrap/wrap_Utils.cpp.
std::unique_ptr<BitVect> create_from_bit_string(rust::Str bits);
std::unique_ptr<BitVect> create_from_fps_text(rust::Str fps);
std::unique_ptr<BitVect> create_from_binary_text(
    rust::Slice<const std::uint8_t> fps);

// Sparse integer vectors. Values are ints for every index type; arithmetic
// that would overflow an int, or divide by zero, throws ValueErrorException.
std::unique_ptr<SparseIntVect> sparse_int_vect_new(std::uint8_t kind,
                                                   std::uint64_t length);
std::unique_ptr<SparseIntVect> sparse_int_vect_from_pkl(
    std::uint8_t kind, rust::Slice<const std::uint8_t> pkl);
std::unique_ptr<SparseIntVect> sparse_int_vect_clone(const SparseIntVect &siv);
std::int32_t sparse_int_vect_get_val(const SparseIntVect &siv,
                                     std::uint64_t idx);
void sparse_int_vect_set_val(SparseIntVect &siv, std::uint64_t idx,
                             std::int32_t val);
std::uint64_t sparse_int_vect_get_length(const SparseIntVect &siv);
std::int32_t sparse_int_vect_get_total_val(const SparseIntVect &siv,
                                           bool use_abs);
rust::Vec<std::uint8_t> sparse_int_vect_to_binary(const SparseIntVect &siv);
void sparse_int_vect_update_from_sequence(SparseIntVect &siv,
                                          rust::Slice<const std::uint64_t> seq);
SparseIntVectElements sparse_int_vect_get_nonzero_elements(
    const SparseIntVect &siv);
rust::Vec<std::int32_t> sparse_int_vect_to_list(const SparseIntVect &siv);
// `op` is '&', '|', '+', or '-'.
std::unique_ptr<SparseIntVect> sparse_int_vect_binary_op(const SparseIntVect &a,
                                                         const SparseIntVect &b,
                                                         std::uint8_t op);
// `op` is '+' or '-'.
void sparse_int_vect_assign_op(SparseIntVect &a, const SparseIntVect &b,
                               std::uint8_t op);
// `op` is '+', '-', '*', or '/'.
void sparse_int_vect_scalar_op(SparseIntVect &siv, std::uint8_t op,
                               std::int32_t val);
bool sparse_int_vect_eq(const SparseIntVect &a, const SparseIntVect &b);
// `metric` is 0 for Dice, 1 for Tanimoto, 2 for Tversky.
double sparse_int_vect_similarity(const SparseIntVect &a,
                                  const SparseIntVect &b, std::uint8_t metric,
                                  double alpha, double beta,
                                  bool return_distance, double bounds);

// Discrete and real value vectors.
std::unique_ptr<DiscreteValueVect> discrete_value_vect_new(
    DiscreteValueType val_type, std::uint32_t length);
std::unique_ptr<DiscreteValueVect> discrete_value_vect_from_pkl(
    rust::Slice<const std::uint8_t> pkl);
std::unique_ptr<DiscreteValueVect> discrete_value_vect_clone(
    const DiscreteValueVect &dvv);
std::uint32_t discrete_value_vect_get_val(const DiscreteValueVect &dvv,
                                          std::uint32_t i);
void discrete_value_vect_set_val(DiscreteValueVect &dvv, std::uint32_t i,
                                 std::uint32_t val);
std::uint32_t discrete_value_vect_get_length(const DiscreteValueVect &dvv);
DiscreteValueType discrete_value_vect_get_value_type(
    const DiscreteValueVect &dvv);
std::uint32_t discrete_value_vect_get_total_val(const DiscreteValueVect &dvv);
rust::Vec<std::uint8_t> discrete_value_vect_to_binary(
    const DiscreteValueVect &dvv);
// `op` is '&', '|', '+', or '-'.
std::unique_ptr<DiscreteValueVect> discrete_value_vect_binary_op(
    const DiscreteValueVect &a, const DiscreteValueVect &b, std::uint8_t op);
void discrete_value_vect_assign_op(DiscreteValueVect &a,
                                   const DiscreteValueVect &b, std::uint8_t op);
std::uint32_t discrete_value_vect_compute_l1_norm(const DiscreteValueVect &a,
                                                  const DiscreteValueVect &b);

std::unique_ptr<RealValueVect> real_value_vect_new(std::uint32_t length);
std::unique_ptr<RealValueVect> real_value_vect_from_pkl(
    rust::Slice<const std::uint8_t> pkl);
std::unique_ptr<RealValueVect> real_value_vect_clone(const RealValueVect &rvv);
double real_value_vect_get_val(const RealValueVect &rvv, std::uint32_t i);
void real_value_vect_set_val(RealValueVect &rvv, std::uint32_t i, double val);
std::uint32_t real_value_vect_get_length(const RealValueVect &rvv);
double real_value_vect_get_total_val(const RealValueVect &rvv);
rust::Vec<std::uint8_t> real_value_vect_to_binary(const RealValueVect &rvv);
std::unique_ptr<RealValueVect> real_value_vect_binary_op(const RealValueVect &a,
                                                         const RealValueVect &b,
                                                         std::uint8_t op);
void real_value_vect_assign_op(RealValueVect &a, const RealValueVect &b,
                               std::uint8_t op);
double real_value_vect_compute_l1_norm(const RealValueVect &a,
                                       const RealValueVect &b);

}  // namespace rdkit_cxx
