// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/info_theory.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <variant>
#include <vector>

#include <ML/InfoTheory/InfoGainFuncs.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

constexpr std::uint64_t kIntMax = std::numeric_limits<int>::max();

// ChiSquare and InfoEntropyGain (InfoGainFuncs.h) truncate the matrix's
// row or column sums to int and add those in int: undefined for a sum past
// int's range, or NaN, which Python's then turns into garbage.
void check_int_sums(const std::vector<double> &sums) {
  std::int64_t total = 0;
  for (auto sum : sums) {
    // Truncating to int is defined for (INT_MIN - 1, INT_MAX + 1).
    if (!(sum > -2147483649.0 && sum < 2147483648.0)) {
      throw ValueErrorException(
          "the matrix's sums must be within int's range: RDKit adds them in "
          "int");
    }
    total += static_cast<std::int64_t>(sum);
    if (total > static_cast<std::int64_t>(kIntMax) ||
        total < std::numeric_limits<int>::min()) {
      throw ValueErrorException(
          "the matrix's sums must be within int's range: RDKit adds them in "
          "int");
    }
  }
}

// The row sums ChiSquare computes, in its order.
std::vector<double> row_sums(rust::Slice<const double> values,
                             std::uint64_t rows, std::uint64_t cols) {
  std::vector<double> sums(rows, 0.0);
  for (std::uint64_t i = 0; i < rows; ++i) {
    for (std::uint64_t j = 0; j < cols; ++j) {
      sums[i] += values[i * cols + j];
    }
  }
  return sums;
}

// The column sums InfoEntropyGain computes, in its order.
std::vector<double> col_sums(rust::Slice<const double> values,
                             std::uint64_t rows, std::uint64_t cols) {
  std::vector<double> sums(cols, 0.0);
  for (std::uint64_t i = 0; i < cols; ++i) {
    for (std::uint64_t j = 0; j < rows; ++j) {
      sums[i] += values[j * cols + i];
    }
  }
  return sums;
}

void check_values(rust::Slice<const double> values, std::uint64_t rows,
                  std::uint64_t cols) {
  if (values.size() != rows * cols) {
    throw ValueErrorException("the matrix must have rows * cols values");
  }
}

// The ranker numbers bits in int (getTopN's pairs, whose bit past INT_MAX
// would index the counts as negative), and InfoEntropyGain adds the
// classes' counts, each up to 65535, in int, which more than 32768 classes
// can overflow.
constexpr std::uint32_t kMaxClasses = 32768;

}  // namespace

// As infoEntropy in Code/ML/InfoTheory/Wrap/rdInfoTheory.cpp. InfoEntropy
// loops over the values in int.
double info_entropy(rust::Slice<const double> res_arr) {
  CHECK_INVARIANT(res_arr.size() > 0, "");
  if (res_arr.size() > kIntMax) {
    throw ValueErrorException("at most 2^31 - 1 values: RDKit counts in int");
  }
  // InfoEntropy reads but does not change the values.
  return RDInfoTheory::InfoEntropy(const_cast<double *>(res_arr.data()),
                                   static_cast<long>(res_arr.size()));
}

// As infoGain. InfoEntropyGain runs InfoEntropy, which loops in int, over
// each row.
double info_gain(rust::Slice<const double> values, std::uint64_t rows,
                 std::uint64_t cols) {
  check_values(values, rows, cols);
  if (cols > kIntMax) {
    throw ValueErrorException(
        "at most 2^31 - 1 columns: RDKit counts them in int");
  }
  check_int_sums(col_sums(values, rows, cols));
  // InfoEntropyGain reads but does not change the matrix.
  return RDInfoTheory::InfoEntropyGain(const_cast<double *>(values.data()),
                                       static_cast<long>(rows),
                                       static_cast<long>(cols));
}

// As chiSquare. ChiSquare indexes the matrix in int.
double chi_square(rust::Slice<const double> values, std::uint64_t rows,
                  std::uint64_t cols) {
  check_values(values, rows, cols);
  if (rows > kIntMax || cols > kIntMax || rows * cols > kIntMax) {
    throw ValueErrorException(
        "at most 2^31 - 1 values: RDKit indexes them in int");
  }
  check_int_sums(row_sums(values, rows, cols));
  // ChiSquare reads but does not change the matrix.
  return RDInfoTheory::ChiSquare(const_cast<double *>(values.data()),
                                 static_cast<long>(rows),
                                 static_cast<long>(cols));
}

std::unique_ptr<InfoBitRanker> info_bit_ranker_new(std::uint32_t n_bits,
                                                   std::uint32_t n_classes,
                                                   InfoType info_type) {
  if (n_bits > kIntMax) {
    throw ValueErrorException(
        "an InfoBitRanker ranks at most 2^31 - 1 bits: RDKit numbers them "
        "in int");
  }
  if (n_classes > kMaxClasses) {
    throw ValueErrorException(
        "an InfoBitRanker takes at most 32768 classes: RDKit adds their "
        "counts in int");
  }
  return std::make_unique<InfoBitRanker>(n_bits, n_classes, info_type);
}

std::uint32_t info_bit_ranker_get_num_classes(const InfoBitRanker &ranker) {
  return ranker.ranker.getNumClasses();
}

// As AccumulateVotes in Code/ML/InfoTheory/Wrap/InfoBitRanker.cpp; RDKit
// checks the label and the vector's size.
void info_bit_ranker_accumulate_votes(InfoBitRanker &ranker,
                                      const BitVect &bit_vect,
                                      std::uint32_t label) {
  std::visit([&](const auto &bv) { ranker.ranker.accumulateVotes(bv, label); },
             bit_vect.v);
}

// RDKit's setBiasList keeps the list before it checks it, so a list it
// refuses stays, and BiasCheckBit then indexes the classes' fractions with
// the list's ids (a segfault in Python for an id past the classes). Its
// checks, made here first, throw as RDKit's do.
void info_bit_ranker_set_bias_list(InfoBitRanker &ranker,
                                   rust::Slice<const std::uint32_t> classes) {
  unsigned int n_classes = ranker.ranker.getNumClasses();
  RDKit::INT_VECT list;
  for (auto c : classes) {
    list.push_back(static_cast<int>(c));
  }
  URANGE_CHECK(list.size(), n_classes);
  RDKit::INT_VECT sorted(list);
  std::sort(sorted.begin(), sorted.end());
  CHECK_INVARIANT(std::unique(sorted.begin(), sorted.end()) == sorted.end(),
                  "There are duplicates in the class bias list");
  for (auto c : sorted) {
    URANGE_CHECK(static_cast<unsigned int>(c), n_classes);
  }
  ranker.ranker.setBiasList(list);
}

// RDKit's setMaskBits replaces the mask before it sets the bits, which
// ExplicitBitVect::setBit refuses past the vector's end, so a refused list
// would leave a partial mask. The indices are checked first (E3).
void info_bit_ranker_set_mask_bits(InfoBitRanker &ranker,
                                   rust::Slice<const std::uint32_t> bits) {
  RDKit::INT_VECT list;
  for (auto bit : bits) {
    check_index(bit, ranker.n_bits);
    list.push_back(static_cast<int>(bit));
  }
  ranker.ranker.setMaskBits(list);
}

// As getTopNbits. With num 0, RDKit's getTopN compares each bit's score
// with the top of its empty priority queue (a segfault in Python); and it
// sizes and indexes its num * (2 + nClasses) results in int and unsigned
// int, which wrap.
rust::Vec<double> info_bit_ranker_get_top_n(InfoBitRanker &ranker,
                                            std::uint32_t num) {
  if (num == 0) {
    throw ValueErrorException(
        "num must be at least 1: RDKit reads its empty ranking");
  }
  std::uint64_t ncols = 2 + std::uint64_t{ranker.ranker.getNumClasses()};
  if (num * ncols > kIntMax) {
    throw ValueErrorException(
        "num * (2 + nClasses) must be at most 2^31 - 1: RDKit indexes the "
        "results in int");
  }
  const double *top = ranker.ranker.getTopN(num);
  rust::Vec<double> result;
  result.reserve(num * ncols);
  for (std::uint64_t i = 0; i < num * ncols; ++i) {
    result.push_back(top[i]);
  }
  return result;
}

// RDKit's writeTopBitsToFile opens the file itself and ignores a failed
// write; the wrapper opens it, has RDKit write to the stream, and checks
// the stream (U12).
void info_bit_ranker_write_top_bits_to_file(
    const InfoBitRanker &ranker, rust::Slice<const std::uint8_t> path) {
  auto name = to_path(path);
  std::ofstream file(name);
  if (!file || file.bad()) {
    throw RDKit::BadFileException("Bad output file " + name);
  }
  ranker.ranker.writeTopBitsToStream(&file);
  file.close();
  if (file.fail()) {
    throw RDKit::BadFileException("could not write the output");
  }
}

std::unique_ptr<BitCorrMatGenerator> bit_corr_mat_generator_new() {
  return std::make_unique<BitCorrMatGenerator>();
}

// RDKit's setBitIdList sizes the matrix as nd * (nd - 1) / 2 in int, which
// overflows for more than 46341 bits (bad_array_new_length in Python, or a
// matrix too small for its writes).
void bit_corr_mat_generator_set_bit_list(
    BitCorrMatGenerator &generator, rust::Slice<const std::uint32_t> bits) {
  if (bits.size() > 46341) {
    throw ValueErrorException(
        "at most 46341 bits: RDKit sizes their matrix in int");
  }
  RDKit::INT_VECT list;
  for (auto bit : bits) {
    list.push_back(static_cast<int>(bit));
  }
  generator.setBitIdList(list);
}

// As CollectVotes in Code/ML/InfoTheory/Wrap/BitCorrMatGenerator.cpp.
// RDKit reads each listed bit of the vector, and a bit past its end throws
// IndexErrorException; the bits are checked first (E3).
void bit_corr_mat_generator_collect_votes(BitCorrMatGenerator &generator,
                                          const BitVect &bit_vect) {
  std::visit(
      [&](const auto &bv) {
        for (auto bit : generator.getCorrBitList()) {
          check_index(static_cast<unsigned int>(bit), bv.getNumBits());
        }
        generator.collectVotes(bv);
      },
      bit_vect.v);
}

// As getCorrMatrix. RDKit has no matrix before a bit list is set.
rust::Vec<double> bit_corr_mat_generator_get_corr_matrix(
    const BitCorrMatGenerator &generator) {
  std::size_t n = generator.getCorrBitList().size();
  std::size_t len = n > 1 ? n * (n - 1) / 2 : 0;
  rust::Vec<double> result;
  if (len == 0) {
    return result;
  }
  // getCorrMat only returns the matrix, but is not const.
  const double *matrix =
      const_cast<BitCorrMatGenerator &>(generator).getCorrMat();
  result.reserve(len);
  for (std::size_t i = 0; i < len; ++i) {
    result.push_back(matrix[i]);
  }
  return result;
}

}  // namespace rdkit_cxx
