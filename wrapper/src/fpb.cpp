/*
 * Copyright (c) 2026 Steven Kearnes All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/fpb.h"

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/fpb.rs.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

// RDKit reads nBits/8 bytes from the pointer it is given.
const std::uint8_t *checked_fp(rust::Slice<const std::uint8_t> bytes,
                               unsigned n_bits) {
  if (bytes.size() != n_bits / 8) {
    throw ValueErrorException(
        "the fingerprint has " + std::to_string(bytes.size()) +
        " bytes; the file's have " + std::to_string(n_bits / 8));
  }
  return bytes.data();
}

rust::Vec<FPBNeighbor> to_neighbors(
    const std::vector<std::pair<double, unsigned int>> &nbrs) {
  rust::Vec<FPBNeighbor> result;
  for (const auto &[sim, idx] : nbrs) {
    result.push_back(FPBNeighbor{sim, idx});
  }
  return result;
}

rust::Vec<MultiFPBNeighbor> to_neighbors(
    const std::vector<MultiFPBReader::ResultTuple> &nbrs) {
  rust::Vec<MultiFPBNeighbor> result;
  for (const auto &nbr : nbrs) {
    result.push_back(
        MultiFPBNeighbor{std::get<0>(nbr), std::get<1>(nbr), std::get<2>(nbr)});
  }
  return result;
}

}  // namespace

std::unique_ptr<FPBReader> fpb_reader_new(rust::Slice<const std::uint8_t> path,
                                          bool lazy) {
  return std::make_unique<FPBReader>(to_std(path), lazy);
}

void fpb_reader_init(FPBReader &reader) { reader.init(); }

std::uint32_t fpb_reader_length(const FPBReader &reader) {
  return reader.length();
}

std::uint32_t fpb_reader_get_num_bits(const FPBReader &reader) {
  return reader.nBits();
}

std::unique_ptr<BitVect> fpb_reader_get_fp(const FPBReader &reader,
                                           std::uint32_t idx) {
  check_index(idx, reader.length());
  auto fp = reader.getFP(idx);
  if (!fp) {
    throw ValueErrorException("null fingerprint pointer");
  }
  return std::make_unique<BitVect>(*fp);
}

rust::Vec<std::uint8_t> fpb_reader_get_bytes(const FPBReader &reader,
                                             std::uint32_t which) {
  check_index(which, reader.length());
  auto bv = reader.getBytes(which);
  rust::Vec<std::uint8_t> result;
  for (unsigned i = 0; i < reader.nBits() / 8; ++i) {
    result.push_back(bv[i]);
  }
  return result;
}

std::unique_ptr<std::string> fpb_reader_get_id(const FPBReader &reader,
                                               std::uint32_t idx) {
  check_index(idx, reader.length());
  return std::make_unique<std::string>(reader.getId(idx));
}

double fpb_reader_get_tanimoto(const FPBReader &reader, std::uint32_t which,
                               rust::Slice<const std::uint8_t> bytes) {
  check_index(which, reader.length());
  return reader.getTanimoto(which, checked_fp(bytes, reader.nBits()));
}

double fpb_reader_get_tversky(const FPBReader &reader, std::uint32_t which,
                              rust::Slice<const std::uint8_t> bytes, double ca,
                              double cb) {
  check_index(which, reader.length());
  return reader.getTversky(which, checked_fp(bytes, reader.nBits()), ca, cb);
}

rust::Vec<FPBNeighbor> fpb_reader_get_tanimoto_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv,
    double threshold) {
  return to_neighbors(
      reader.getTanimotoNeighbors(checked_fp(bv, reader.nBits()), threshold));
}

rust::Vec<FPBNeighbor> fpb_reader_get_tversky_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv, double ca,
    double cb, double threshold) {
  return to_neighbors(reader.getTverskyNeighbors(checked_fp(bv, reader.nBits()),
                                                 ca, cb, threshold));
}

rust::Vec<std::uint32_t> fpb_reader_get_containing_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv) {
  return to_vec(reader.getContainingNeighbors(checked_fp(bv, reader.nBits())));
}

std::unique_ptr<MultiFPBReader> multi_fpb_reader_new(bool init_on_search) {
  return std::make_unique<MultiFPBReader>(init_on_search);
}

std::uint32_t multi_fpb_reader_add_reader(MultiFPBReader &multi,
                                          FPBReader &reader) {
  return multi.addReader(&reader);
}

void multi_fpb_reader_init(MultiFPBReader &multi) { multi.init(); }

std::uint32_t multi_fpb_reader_get_num_bits(const MultiFPBReader &multi) {
  return multi.nBits();
}

rust::Vec<MultiFPBNeighbor> multi_fpb_reader_get_tanimoto_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv,
    double threshold, std::uint32_t num_threads) {
  return to_neighbors(multi.getTanimotoNeighbors(checked_fp(bv, multi.nBits()),
                                                 threshold, num_threads));
}

rust::Vec<MultiFPBNeighbor> multi_fpb_reader_get_tversky_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv, double ca,
    double cb, double threshold, std::uint32_t num_threads) {
  return to_neighbors(multi.getTverskyNeighbors(
      checked_fp(bv, multi.nBits()), ca, cb, threshold, num_threads));
}

rust::Vec<MultiFPBContaining> multi_fpb_reader_get_containing_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv,
    std::uint32_t num_threads) {
  rust::Vec<MultiFPBContaining> result;
  for (const auto &[idx, reader] : multi.getContainingNeighbors(
           checked_fp(bv, multi.nBits()), num_threads)) {
    result.push_back(MultiFPBContaining{idx, reader});
  }
  return result;
}

}  // namespace rdkit_cxx
