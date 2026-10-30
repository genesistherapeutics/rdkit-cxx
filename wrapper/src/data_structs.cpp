/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/data_structs.h"

#include <algorithm>
#include <limits>
#include <type_traits>

#include <DataStructs/BitOps.h>
#include <DataStructs/BitVectUtils.h>
#include <DataStructs/base64.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/data_structs.rs.h"
#include "rdkit-cxx/pickle_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

// Bit indices in ascending order. A SparseBitVect keeps its bits in a
// std::set<int>, so bits past INT_MAX (from fingerprint generators) come
// first until sorted as the unsigned indices they are.
template <class T>
rust::Vec<std::uint32_t> to_indices(const std::vector<T> &v) {
  std::vector<std::uint32_t> sorted;
  for (auto i : v) {
    sorted.push_back(static_cast<std::uint32_t>(i));
  }
  std::sort(sorted.begin(), sorted.end());
  rust::Vec<std::uint32_t> result;
  for (auto i : sorted) {
    result.push_back(i);
  }
  return result;
}

// An index error for 64-bit indices, which may not fit IndexError's signed
// index.
void check_index_64(std::uint64_t idx, std::uint64_t len) {
  if (idx >= len) {
    constexpr std::uint64_t most = std::numeric_limits<std::int64_t>::max();
    throw IndexError(static_cast<std::int64_t>(std::min(idx, most)), len);
  }
}

// Calls f with both operands' RDKit objects, which must be the same kind.
template <class V, class F>
auto visit2(const V &a, const V &b, F f) {
  if (a.v.index() != b.v.index()) {
    throw ValueErrorException("vectors of different types");
  }
  return std::visit(
      [&](const auto &x) {
        using T = std::decay_t<decltype(x)>;
        return f(x, std::get<T>(b.v));
      },
      a.v);
}

// RDKit indexes bits with int in places (getOnBits, FoldFingerprint, and
// SparseBitVect's storage), so a vector holds at most INT_MAX bits.
void check_num_bits(std::uint64_t size) {
  if (size > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw ValueErrorException("bit vectors hold at most 2^31 - 1 bits");
  }
}

// Fingerprint generators make SparseBitVects of 2^32 - 1 bits (U3), whose
// bits past INT_MAX RDKit stores as negative ints. Pickling one subtracts
// those ints in int, which can overflow, and converting one makes an
// ExplicitBitVect too large for its int arithmetic, so both refuse them.
void check_int_sized(const BitVect &bv) {
  check_num_bits(bit_vect_get_num_bits(bv));
}

void check_same_size(const BitVect &a, const BitVect &b) {
  if (bit_vect_get_num_bits(a) != bit_vect_get_num_bits(b)) {
    throw ValueErrorException("BitVects must be same length");
  }
}

template <class T>
double metric_value(const T &a, const T &b, std::uint8_t metric, double alpha,
                    double beta) {
  switch (metric) {
    case 0:
      return TanimotoSimilarity(a, b);
    case 1:
      return CosineSimilarity(a, b);
    case 2:
      return KulczynskiSimilarity(a, b);
    case 3:
      return DiceSimilarity(a, b);
    case 4:
      return SokalSimilarity(a, b);
    case 5:
      return McConnaugheySimilarity(a, b);
    case 6:
      return AsymmetricSimilarity(a, b);
    case 7:
      return BraunBlanquetSimilarity(a, b);
    case 8:
      return RusselSimilarity(a, b);
    case 9:
      return RogotGoldbergSimilarity(a, b);
    case 10:
      return OnBitSimilarity(a, b);
    case 11:
      return AllBitSimilarity(a, b);
    case 12:
      return TverskySimilarity(a, b, alpha, beta);
  }
  throw ValueErrorException("unknown similarity metric");
}

std::uint8_t kind_of(const SparseIntVect &siv) {
  return static_cast<std::uint8_t>(siv.v.index());
}

template <class F>
auto with_kind(std::uint8_t kind, F f) {
  switch (kind) {
    case 0:
      return f(std::int32_t{});
    case 1:
      return f(std::int64_t{});
    case 2:
      return f(std::uint32_t{});
    case 3:
      return f(std::uint64_t{});
  }
  throw ValueErrorException("unknown SparseIntVect kind");
}

int checked(long long value) {
  if (value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw ValueErrorException("SparseIntVect value overflow");
  }
  return static_cast<int>(value);
}

// Checks that adding (or subtracting) b to a overflows no value.
template <class IndexType>
void check_add(const RDKit::SparseIntVect<IndexType> &a,
               const RDKit::SparseIntVect<IndexType> &b, int sign) {
  if (a.getLength() != b.getLength()) {
    throw ValueErrorException("SparseIntVect size mismatch");
  }
  for (const auto &[idx, val] : b.getNonzeroElements()) {
    checked(static_cast<long long>(a.getVal(idx)) +
            sign * static_cast<long long>(val));
  }
}

template <class V>
void check_lengths(const V &a, const V &b) {
  if (a.getLength() != b.getLength()) {
    throw ValueErrorException("length mismatch");
  }
}

}  // namespace

const ExplicitBitVect &as_explicit(const BitVect &bv) {
  return std::get<ExplicitBitVect>(bv.v);
}

std::unique_ptr<BitVect> wrap_bit_vect(std::unique_ptr<ExplicitBitVect> bv) {
  return std::make_unique<BitVect>(std::move(*bv));
}

std::unique_ptr<BitVect> wrap_bit_vect(std::unique_ptr<SparseBitVect> bv) {
  return std::make_unique<BitVect>(std::move(*bv));
}

std::unique_ptr<BitVect> bit_vect_new_explicit(std::uint32_t size,
                                               bool bits_set) {
  check_num_bits(size);
  return std::make_unique<BitVect>(ExplicitBitVect(size, bits_set));
}

std::unique_ptr<BitVect> bit_vect_new_sparse(std::uint32_t size) {
  check_num_bits(size);
  return std::make_unique<BitVect>(SparseBitVect(size));
}

std::unique_ptr<BitVect> bit_vect_from_pkl(
    bool sparse, rust::Slice<const std::uint8_t> pkl) {
  auto text = to_std(pkl);
  check_bit_vect_pickle(text);
  if (sparse) {
    return std::make_unique<BitVect>(SparseBitVect(text));
  }
  return std::make_unique<BitVect>(ExplicitBitVect(text));
}

std::unique_ptr<BitVect> bit_vect_clone(const BitVect &bv) {
  return std::make_unique<BitVect>(bv);
}

bool bit_vect_set_bit(BitVect &bv, std::uint32_t which) {
  check_index(which, bit_vect_get_num_bits(bv));
  return std::visit([&](auto &x) { return x.setBit(which); }, bv.v);
}

bool bit_vect_un_set_bit(BitVect &bv, std::uint32_t which) {
  check_index(which, bit_vect_get_num_bits(bv));
  return std::visit([&](auto &x) { return x.unsetBit(which); }, bv.v);
}

void bit_vect_set_bits_from_list(BitVect &bv,
                                 rust::Slice<const std::uint32_t> on_bit_list) {
  for (auto which : on_bit_list) {
    check_index(which, bit_vect_get_num_bits(bv));
  }
  for (auto which : on_bit_list) {
    bit_vect_set_bit(bv, which);
  }
}

void bit_vect_un_set_bits_from_list(
    BitVect &bv, rust::Slice<const std::uint32_t> off_bit_list) {
  for (auto which : off_bit_list) {
    check_index(which, bit_vect_get_num_bits(bv));
  }
  for (auto which : off_bit_list) {
    bit_vect_un_set_bit(bv, which);
  }
}

bool bit_vect_get_bit(const BitVect &bv, std::uint32_t which) {
  check_index(which, bit_vect_get_num_bits(bv));
  return std::visit([&](const auto &x) { return x.getBit(which); }, bv.v);
}

std::uint32_t bit_vect_get_num_bits(const BitVect &bv) {
  return std::visit([](const auto &x) { return x.getNumBits(); }, bv.v);
}

std::uint32_t bit_vect_get_num_on_bits(const BitVect &bv) {
  return std::visit([](const auto &x) { return x.getNumOnBits(); }, bv.v);
}

std::uint32_t bit_vect_get_num_off_bits(const BitVect &bv) {
  return std::visit([](const auto &x) { return x.getNumOffBits(); }, bv.v);
}

rust::Vec<std::uint32_t> bit_vect_get_on_bits(const BitVect &bv) {
  IntVect res;
  std::visit([&](const auto &x) { x.getOnBits(res); }, bv.v);
  return to_indices(res);
}

rust::Vec<std::uint8_t> bit_vect_to_binary(const BitVect &bv) {
  check_int_sized(bv);
  return to_bytes(std::visit([](const auto &x) { return x.toString(); }, bv.v));
}

void bit_vect_from_base64(BitVect &bv, rust::Str in_d) {
  // Base64Decode writes past its buffer for some lengths. It skips '=', so
  // padding makes the buffer long enough without changing the result; it
  // stops at a NUL, so the padding goes there. (It also reads one
  // uninitialized table entry, for byte 0xFF, which UTF-8 never has.)
  std::string text(in_d);
  text = text.substr(0, text.find('\0')) + "====";
  unsigned int len = 0;
  std::unique_ptr<char[]> decoded(Base64Decode(text.c_str(), &len));
  std::string pkl(decoded.get(), len);
  check_bit_vect_pickle(pkl);
  // RDKit's in-place parse leaves the vector broken when it fails, so parse
  // a new vector and replace this one.
  std::visit(
      [&](auto &x) {
        using T = std::decay_t<decltype(x)>;
        x = T(pkl);
      },
      bv.v);
}

std::unique_ptr<std::string> bit_vect_to_base64(const BitVect &bv) {
  check_int_sized(bv);
  auto tmp = std::visit([](const auto &x) { return x.toString(); }, bv.v);
  std::unique_ptr<const char[]> txt(Base64Encode(tmp.c_str(), tmp.length()));
  return std::make_unique<std::string>(txt.get());
}

rust::Vec<bool> bit_vect_to_list(const BitVect &bv) {
  rust::Vec<bool> result;
  auto n = bit_vect_get_num_bits(bv);
  result.reserve(n);
  for (std::uint32_t i = 0; i < n; ++i) {
    result.push_back(false);
  }
  for (auto i : bit_vect_get_on_bits(bv)) {
    result[i] = true;
  }
  return result;
}

std::unique_ptr<BitVect> bit_vect_binary_op(const BitVect &a, const BitVect &b,
                                            std::uint8_t op) {
  if (op != '+') {
    check_same_size(a, b);
  }
  return visit2(a, b, [&](const auto &x, const auto &y) {
    using T = std::decay_t<decltype(x)>;
    switch (op) {
      case '&':
        return std::make_unique<BitVect>(T(x & y));
      case '|':
        return std::make_unique<BitVect>(T(x | y));
      case '^':
        return std::make_unique<BitVect>(T(x ^ y));
      case '+':
        if constexpr (std::is_same_v<T, ExplicitBitVect>) {
          check_num_bits(std::uint64_t{x.getNumBits()} + y.getNumBits());
          return std::make_unique<BitVect>(T(x + y));
        }
    }
    throw ValueErrorException("unsupported bit vector operation");
  });
}

void bit_vect_append(BitVect &bv, const BitVect &other) {
  check_num_bits(std::uint64_t{bit_vect_get_num_bits(bv)} +
                 bit_vect_get_num_bits(other));
  std::get<ExplicitBitVect>(bv.v) += std::get<ExplicitBitVect>(other.v);
}

std::unique_ptr<BitVect> bit_vect_invert(const BitVect &bv) {
  return std::visit(
      [](const auto &x) {
        using T = std::decay_t<decltype(x)>;
        return std::make_unique<BitVect>(T(~x));
      },
      bv.v);
}

bool bit_vect_eq(const BitVect &a, const BitVect &b) {
  return a.v.index() == b.v.index() &&
         visit2(a, b, [](const auto &x, const auto &y) { return x == y; });
}

double bit_vect_similarity(const BitVect &a, const BitVect &b,
                           std::uint8_t metric, double alpha, double beta,
                           bool return_distance) {
  auto sim = visit2(a, b, [&](const auto &x, const auto &y) {
    return metric_value(x, y, metric, alpha, beta);
  });
  return return_distance ? 1.0 - sim : sim;
}

double bit_vect_similarity_pkl(const BitVect &a,
                               rust::Slice<const std::uint8_t> pkl,
                               std::uint8_t metric, double alpha, double beta,
                               bool return_distance) {
  auto text = to_std(pkl);
  check_bit_vect_pickle(text);
  auto sim = std::visit(
      [&](const auto &x) {
        using T = std::decay_t<decltype(x)>;
        T other(text);
        return metric_value(x, other, metric, alpha, beta);
      },
      a.v);
  return return_distance ? 1.0 - sim : sim;
}

rust::Vec<double> bit_vect_proj_similarity(const BitVect &a, const BitVect &b,
                                           bool on_bits) {
  auto res = visit2(a, b, [&](const auto &x, const auto &y) {
    return on_bits ? OnBitProjSimilarity(x, y) : OffBitProjSimilarity(x, y);
  });
  rust::Vec<double> result;
  for (auto v : res) {
    result.push_back(v);
  }
  return result;
}

std::int32_t bit_vect_num_bits_in_common(const BitVect &a, const BitVect &b) {
  // RDKit's ExplicitBitVect overload has no size check.
  check_same_size(a, b);
  return visit2(
      a, b, [](const auto &x, const auto &y) { return NumBitsInCommon(x, y); });
}

rust::Vec<std::uint32_t> bit_vect_bits_in_common(const BitVect &a,
                                                 const BitVect &b,
                                                 bool on_bits) {
  return to_indices(visit2(a, b, [&](const auto &x, const auto &y) {
    return on_bits ? OnBitsInCommon(x, y) : OffBitsInCommon(x, y);
  }));
}

std::unique_ptr<BitVect> bit_vect_fold(const BitVect &bv,
                                       std::uint32_t fold_factor) {
  return std::visit(
      [&](const auto &x) {
        using T = std::decay_t<decltype(x)>;
        return wrap_bit_vect(
            std::unique_ptr<T>(FoldFingerprint(x, fold_factor)));
      },
      bv.v);
}

bool bit_vect_all_probe_bits_match(const BitVect &probe, const BitVect &ref) {
  // RDKit's ExplicitBitVect overload has no size check.
  check_same_size(probe, ref);
  return visit2(probe, ref, [](const auto &x, const auto &y) {
    return AllProbeBitsMatch(x, y);
  });
}

// RDKit's pickle overload reads the pickle through raw pointers with no
// bounds checks, so this reads the pickle into a vector and compares bits
// as that overload does: whether each of the probe's on bits is on in the
// reference, whatever their sizes.
bool bit_vect_all_probe_bits_match_pkl(const BitVect &probe,
                                       rust::Slice<const std::uint8_t> ref) {
  auto text = to_std(ref);
  check_bit_vect_pickle(text);
  return std::visit(
      [&](const auto &x) {
        using T = std::decay_t<decltype(x)>;
        T other(text);
        std::vector<int> on;
        x.getOnBits(on);
        return std::all_of(on.begin(), on.end(), [&](int i) {
          return static_cast<unsigned>(i) < other.getNumBits() &&
                 other.getBit(static_cast<unsigned>(i));
        });
      },
      probe.v);
}

std::unique_ptr<std::string> bit_vect_to_text(const BitVect &bv) {
  return std::make_unique<std::string>(
      std::visit([](const auto &x) { return BitVectToText(x); }, bv.v));
}

std::unique_ptr<std::string> bit_vect_to_fps_text(const BitVect &bv) {
  return std::make_unique<std::string>(
      std::visit([](const auto &x) { return BitVectToFPSText(x); }, bv.v));
}

rust::Vec<std::uint8_t> bit_vect_to_binary_text(const BitVect &bv) {
  return to_bytes(
      std::visit([](const auto &x) { return BitVectToBinaryText(x); }, bv.v));
}

void bit_vect_init_from_daylight_string(BitVect &bv, rust::Str s) {
  // FromDaylightString reads past the string unless, after an optional
  // newline, it is groups of four characters and a pad character.
  std::string text(s);
  auto length = text.size();
  if (length > 0 && text[length - 1] == '\n') {
    --length;
  }
  if (length == 0 || (length - 1) % 4 != 0 ||
      (text[length - 1] != '1' && text[length - 1] != '2' &&
       text[length - 1] != '3')) {
    throw ValueErrorException("bad daylight fingerprint string");
  }
  std::visit([&](auto &x) { FromDaylightString(x, text); }, bv.v);
}

std::unique_ptr<BitVect> convert_to_explicit(const BitVect &sbv) {
  check_int_sized(sbv);
  return wrap_bit_vect(std::unique_ptr<ExplicitBitVect>(
      convertToExplicit(&std::get<SparseBitVect>(sbv.v))));
}

std::unique_ptr<BitVect> create_from_bit_string(rust::Str bits) {
  std::string text(bits);
  check_num_bits(text.size());
  ExplicitBitVect bv(text.size());
  FromBitString(bv, text);
  return std::make_unique<BitVect>(std::move(bv));
}

std::unique_ptr<BitVect> create_from_fps_text(rust::Str fps) {
  std::string text(fps);
  if (text.size() % 2) {
    throw ValueErrorException(
        "input string must have an even number of characters");
  }
  check_num_bits(std::uint64_t{text.size()} * 4);
  ExplicitBitVect bv(text.size() * 4);
  UpdateBitVectFromFPSText(bv, text);
  return std::make_unique<BitVect>(std::move(bv));
}

std::unique_ptr<BitVect> create_from_binary_text(
    rust::Slice<const std::uint8_t> fps) {
  check_num_bits(std::uint64_t{fps.size()} * 8);
  ExplicitBitVect bv(fps.size() * 8);
  UpdateBitVectFromBinaryText(bv, to_std(fps));
  return std::make_unique<BitVect>(std::move(bv));
}

std::unique_ptr<SparseIntVect> sparse_int_vect_new(std::uint8_t kind,
                                                   std::uint64_t length) {
  return with_kind(kind, [&](auto index) {
    using IndexType = decltype(index);
    if (length >
        static_cast<std::uint64_t>(std::numeric_limits<IndexType>::max())) {
      throw ValueErrorException("SparseIntVect length out of range");
    }
    return std::make_unique<SparseIntVect>(
        RDKit::SparseIntVect<IndexType>(static_cast<IndexType>(length)));
  });
}

std::unique_ptr<SparseIntVect> sparse_int_vect_from_pkl(
    std::uint8_t kind, rust::Slice<const std::uint8_t> pkl) {
  return with_kind(kind, [&](auto index) {
    using IndexType = decltype(index);
    RDKit::SparseIntVect<IndexType> siv(to_std(pkl));
    check_sparse_int_vect(siv);
    return std::make_unique<SparseIntVect>(std::move(siv));
  });
}

std::unique_ptr<SparseIntVect> sparse_int_vect_clone(const SparseIntVect &siv) {
  return std::make_unique<SparseIntVect>(siv);
}

std::int32_t sparse_int_vect_get_val(const SparseIntVect &siv,
                                     std::uint64_t idx) {
  check_index_64(idx, sparse_int_vect_get_length(siv));
  return std::visit(
      [&](const auto &x) {
        using IndexType =
            typename std::decay_t<decltype(x)>::StorageType::key_type;
        return x.getVal(static_cast<IndexType>(idx));
      },
      siv.v);
}

void sparse_int_vect_set_val(SparseIntVect &siv, std::uint64_t idx,
                             std::int32_t val) {
  check_index_64(idx, sparse_int_vect_get_length(siv));
  std::visit(
      [&](auto &x) {
        using IndexType =
            typename std::decay_t<decltype(x)>::StorageType::key_type;
        x.setVal(static_cast<IndexType>(idx), val);
      },
      siv.v);
}

std::uint64_t sparse_int_vect_get_length(const SparseIntVect &siv) {
  return std::visit(
      [](const auto &x) { return static_cast<std::uint64_t>(x.getLength()); },
      siv.v);
}

std::int32_t sparse_int_vect_get_total_val(const SparseIntVect &siv,
                                           bool use_abs) {
  return std::visit(
      [&](const auto &x) {
        check_sparse_int_vect_sums(x);
        return x.getTotalVal(use_abs);
      },
      siv.v);
}

rust::Vec<std::uint8_t> sparse_int_vect_to_binary(const SparseIntVect &siv) {
  return to_bytes(
      std::visit([](const auto &x) { return x.toString(); }, siv.v));
}

void sparse_int_vect_update_from_sequence(
    SparseIntVect &siv, rust::Slice<const std::uint64_t> seq) {
  auto length = sparse_int_vect_get_length(siv);
  for (auto idx : seq) {
    check_index_64(idx, length);
  }
  for (auto idx : seq) {
    sparse_int_vect_set_val(
        siv, idx,
        checked(static_cast<long long>(sparse_int_vect_get_val(siv, idx)) + 1));
  }
}

SparseIntVectElements sparse_int_vect_get_nonzero_elements(
    const SparseIntVect &siv) {
  SparseIntVectElements elements;
  std::visit(
      [&](const auto &x) {
        for (const auto &[idx, val] : x.getNonzeroElements()) {
          elements.indices.push_back(static_cast<std::uint64_t>(idx));
          elements.values.push_back(val);
        }
      },
      siv.v);
  return elements;
}

rust::Vec<std::int32_t> sparse_int_vect_to_list(const SparseIntVect &siv) {
  rust::Vec<std::int32_t> result;
  auto length = sparse_int_vect_get_length(siv);
  result.reserve(length);
  for (std::uint64_t i = 0; i < length; ++i) {
    result.push_back(0);
  }
  auto elements = sparse_int_vect_get_nonzero_elements(siv);
  for (std::size_t i = 0; i < elements.indices.size(); ++i) {
    result[elements.indices[i]] = elements.values[i];
  }
  return result;
}

std::unique_ptr<SparseIntVect> sparse_int_vect_binary_op(const SparseIntVect &a,
                                                         const SparseIntVect &b,
                                                         std::uint8_t op) {
  return visit2(a, b, [&](const auto &x, const auto &y) {
    switch (op) {
      case '&':
        return std::make_unique<SparseIntVect>(x & y);
      case '|':
        return std::make_unique<SparseIntVect>(x | y);
      case '+':
        check_add(x, y, 1);
        return std::make_unique<SparseIntVect>(x + y);
      case '-':
        check_add(x, y, -1);
        return std::make_unique<SparseIntVect>(x - y);
    }
    throw ValueErrorException("unsupported SparseIntVect operation");
  });
}

void sparse_int_vect_assign_op(SparseIntVect &a, const SparseIntVect &b,
                               std::uint8_t op) {
  if (a.v.index() != b.v.index()) {
    throw ValueErrorException("vectors of different types");
  }
  std::visit(
      [&](auto &x) {
        using T = std::decay_t<decltype(x)>;
        const auto &y = std::get<T>(b.v);
        check_add(x, y, op == '+' ? 1 : -1);
        if (op == '+') {
          x += y;
        } else {
          x -= y;
        }
      },
      a.v);
}

void sparse_int_vect_scalar_op(SparseIntVect &siv, std::uint8_t op,
                               std::int32_t val) {
  if (op == '/' && val == 0) {
    throw ValueErrorException("SparseIntVect division by zero");
  }
  std::visit(
      [&](auto &x) {
        // RDKit applies these to the nonzero values only.
        for (const auto &[idx, v] : x.getNonzeroElements()) {
          long long a = v;
          long long b = val;
          switch (op) {
            case '+':
              checked(a + b);
              break;
            case '-':
              checked(a - b);
              break;
            case '*':
              checked(a * b);
              break;
            case '/':
              checked(a / b);
              break;
          }
        }
        switch (op) {
          case '+':
            x += val;
            break;
          case '-':
            x -= val;
            break;
          case '*':
            x *= val;
            break;
          case '/':
            x /= val;
            break;
        }
      },
      siv.v);
}

bool sparse_int_vect_eq(const SparseIntVect &a, const SparseIntVect &b) {
  return a.v.index() == b.v.index() &&
         visit2(a, b, [](const auto &x, const auto &y) { return x == y; });
}

double sparse_int_vect_similarity(const SparseIntVect &a,
                                  const SparseIntVect &b, std::uint8_t metric,
                                  double alpha, double beta,
                                  bool return_distance, double bounds) {
  return visit2(a, b, [&](const auto &x, const auto &y) {
    check_sparse_int_vect_sums(x);
    check_sparse_int_vect_sums(y);
    switch (metric) {
      case 0:
        return DiceSimilarity(x, y, return_distance, bounds);
      case 1:
        return TanimotoSimilarity(x, y, return_distance, bounds);
      case 2:
        return TverskySimilarity(x, y, alpha, beta, return_distance, bounds);
    }
    throw ValueErrorException("unknown similarity metric");
  });
}

std::unique_ptr<DiscreteValueVect> discrete_value_vect_new(
    DiscreteValueType val_type, std::uint32_t length) {
  // RDKit's count of storage ints wraps for lengths near UINT_MAX.
  if (length > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
    throw ValueErrorException(
        "DiscreteValueVects hold at most 2^31 - 1 values");
  }
  return std::make_unique<DiscreteValueVect>(val_type, length);
}

std::unique_ptr<DiscreteValueVect> discrete_value_vect_from_pkl(
    rust::Slice<const std::uint8_t> pkl) {
  auto text = to_std(pkl);
  check_discrete_value_vect_pickle(text);
  return std::make_unique<DiscreteValueVect>(text);
}

std::unique_ptr<DiscreteValueVect> discrete_value_vect_clone(
    const DiscreteValueVect &dvv) {
  return std::make_unique<DiscreteValueVect>(dvv);
}

std::uint32_t discrete_value_vect_get_val(const DiscreteValueVect &dvv,
                                          std::uint32_t i) {
  check_index(i, dvv.getLength());
  return dvv.getVal(i);
}

void discrete_value_vect_set_val(DiscreteValueVect &dvv, std::uint32_t i,
                                 std::uint32_t val) {
  check_index(i, dvv.getLength());
  dvv.setVal(i, val);
}

std::uint32_t discrete_value_vect_get_length(const DiscreteValueVect &dvv) {
  return dvv.getLength();
}

DiscreteValueType discrete_value_vect_get_value_type(
    const DiscreteValueVect &dvv) {
  return dvv.getValueType();
}

std::uint32_t discrete_value_vect_get_total_val(const DiscreteValueVect &dvv) {
  return dvv.getTotalVal();
}

rust::Vec<std::uint8_t> discrete_value_vect_to_binary(
    const DiscreteValueVect &dvv) {
  return to_bytes(dvv.toString());
}

std::unique_ptr<DiscreteValueVect> discrete_value_vect_binary_op(
    const DiscreteValueVect &a, const DiscreteValueVect &b, std::uint8_t op) {
  check_lengths(a, b);
  switch (op) {
    case '&':
      return std::make_unique<DiscreteValueVect>(a & b);
    case '|':
      return std::make_unique<DiscreteValueVect>(a | b);
    case '+':
      return std::make_unique<DiscreteValueVect>(a + b);
    case '-':
      return std::make_unique<DiscreteValueVect>(a - b);
  }
  throw ValueErrorException("unsupported DiscreteValueVect operation");
}

void discrete_value_vect_assign_op(DiscreteValueVect &a,
                                   const DiscreteValueVect &b,
                                   std::uint8_t op) {
  check_lengths(a, b);
  if (op == '+') {
    a += b;
  } else {
    a -= b;
  }
}

std::uint32_t discrete_value_vect_compute_l1_norm(const DiscreteValueVect &a,
                                                  const DiscreteValueVect &b) {
  return RDKit::computeL1Norm(a, b);
}

}  // namespace rdkit_cxx
