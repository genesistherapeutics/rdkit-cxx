/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdInfoTheory: entropy, information gain, and chi-squared statistics of
// count matrices, a ranker of bits by them, and a bit correlation matrix
// (Code/ML/InfoTheory, whose headers are in wrapper/vendor).

#pragma once

#include <cstdint>
#include <memory>

#include <ML/InfoTheory/CorrMatGenerator.h>
#include <ML/InfoTheory/InfoBitRanker.h>

#include "rdkit-cxx/data_structs.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using InfoType = RDInfoTheory::InfoBitRanker::InfoType;
using BitCorrMatGenerator = RDInfoTheory::BitCorrMatGenerator;

// RDKit's ranker, with its size, which it keeps private. RDKit's class
// deletes its arrays and copies none, so neither is copied.
class InfoBitRanker {
 public:
  InfoBitRanker(unsigned int n_bits, unsigned int n_classes, InfoType type)
      : ranker(n_bits, n_classes, type), n_bits(n_bits) {}
  InfoBitRanker(const InfoBitRanker &) = delete;
  InfoBitRanker &operator=(const InfoBitRanker &) = delete;
  RDInfoTheory::InfoBitRanker ranker;
  unsigned int n_bits;
};

// The matrices are rows by cols values, by rows.
double info_entropy(rust::Slice<const double> res_arr);
double info_gain(rust::Slice<const double> values, std::uint64_t rows,
                 std::uint64_t cols);
double chi_square(rust::Slice<const double> values, std::uint64_t rows,
                  std::uint64_t cols);

std::unique_ptr<InfoBitRanker> info_bit_ranker_new(std::uint32_t n_bits,
                                                   std::uint32_t n_classes,
                                                   InfoType info_type);
std::uint32_t info_bit_ranker_get_num_classes(const InfoBitRanker &ranker);
void info_bit_ranker_accumulate_votes(InfoBitRanker &ranker,
                                      const BitVect &bit_vect,
                                      std::uint32_t label);
void info_bit_ranker_set_bias_list(InfoBitRanker &ranker,
                                   rust::Slice<const std::uint32_t> classes);
void info_bit_ranker_set_mask_bits(InfoBitRanker &ranker,
                                   rust::Slice<const std::uint32_t> bits);
// num rows of 2 + nClasses values: the bit, its score, and its count in
// each class.
rust::Vec<double> info_bit_ranker_get_top_n(InfoBitRanker &ranker,
                                            std::uint32_t num);
void info_bit_ranker_write_top_bits_to_file(
    const InfoBitRanker &ranker, rust::Slice<const std::uint8_t> path);

std::unique_ptr<BitCorrMatGenerator> bit_corr_mat_generator_new();
void bit_corr_mat_generator_set_bit_list(BitCorrMatGenerator &generator,
                                         rust::Slice<const std::uint32_t> bits);
void bit_corr_mat_generator_collect_votes(BitCorrMatGenerator &generator,
                                          const BitVect &bit_vect);
rust::Vec<double> bit_corr_mat_generator_get_corr_matrix(
    const BitCorrMatGenerator &generator);

}  // namespace rdkit_cxx
