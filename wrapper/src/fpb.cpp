/*
 * Copyright (c) 2026 Steven Kearnes. All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/fpb.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <fstream>
#include <future>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <DataStructs/BitOps.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDThreads.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/fpb.rs.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

constexpr std::uint64_t kMaxU32 = std::numeric_limits<std::uint32_t>::max();

[[noreturn]] void malformed(const std::string &what) {
  throw ValueErrorException("malformed FPB file: " + what);
}

template <class T>
T read_le(const std::string &data, std::uint64_t offset) {
  unsigned char bytes[sizeof(T)];
  std::memcpy(bytes, data.data() + offset, sizeof(T));
  if constexpr (std::endian::native == std::endian::big) {
    std::reverse(bytes, bytes + sizeof(T));
  }
  T value;
  std::memcpy(&value, bytes, sizeof(T));
  return value;
}

// Checks the chunks RDKit's FPBReader::init reads (FPBReader.cpp), in the
// same order, against the invariants its later reads rely on but never
// check: that each chunk fits in the file, and that the arena's and the id
// table's sizes and offsets fit in their chunks and in the 32-bit arithmetic
// RDKit does with them.
FpbLayout check_fpb(const std::string &data) {
  constexpr char kMagic[] = "FPB1\r\n\0\0";
  if (data.size() < 8 || std::memcmp(data.data(), kMagic, 8) != 0) {
    throw RDKit::BadFileException("Invalid FPB magic");
  }
  FpbLayout layout;
  int num_arenas = 0, num_ids = 0;
  std::uint64_t pos = 8;
  while (true) {
    if (data.size() - pos < 12) {
      malformed("no FEND chunk");
    }
    auto size = read_le<std::uint64_t>(data, pos);
    // RDKit reads the tag as a C string.
    const char *tag_bytes = data.data() + pos + 8;
    std::string tag(tag_bytes, strnlen(tag_bytes, 4));
    pos += 12;
    if (size > data.size() - pos) {
      malformed("the " + tag + " chunk runs past the end of the file");
    }
    if (tag == "FEND") {
      break;
    } else if (tag == "POPC") {
      // RDKit checks the size itself, and bounds the scans the counts give.
      layout.num_pop_counts += size / 4;
    } else if (tag == "AREN") {
      if (++num_arenas > 1) {
        malformed("more than one AREN chunk");
      }
      if (size < 9) {
        malformed("a short AREN chunk");
      }
      layout.num_bytes = read_le<std::uint32_t>(data, pos);
      layout.storage_size = read_le<std::uint32_t>(data, pos + 4);
      std::uint8_t spacer = data[pos + 8];
      if (layout.storage_size == 0) {
        malformed("a zero fingerprint storage size");
      }
      if (layout.num_bytes > layout.storage_size) {
        malformed("fingerprints larger than their storage");
      }
      // ExplicitBitVects hold at most INT_MAX bits, and the lazy search
      // buffer is 1000 fingerprints, sized in 32 bits.
      if (layout.num_bytes > std::numeric_limits<int>::max() / 8 ||
          layout.storage_size > kMaxU32 / 1000) {
        malformed("fingerprints too large");
      }
      if (size - 9 < spacer || (size - 9 - spacer) % layout.storage_size) {
        malformed("an AREN chunk that is not a whole number of fingerprints");
      }
      // RDKit finds fingerprints and skips the arena with 32-bit offsets.
      if (size - 9 - spacer > kMaxU32) {
        malformed("4 GiB of fingerprints or more");
      }
    } else if (tag == "FPID") {
      if (++num_ids > 1) {
        malformed("more than one FPID chunk");
      }
      if (size < 8) {
        malformed("a short FPID chunk");
      }
      std::uint64_t num4 = read_le<std::uint32_t>(data, pos);
      std::uint64_t num8 = read_le<std::uint32_t>(data, pos + 4);
      // RDKit misindexes 8-byte offsets ("not yet tested" in its source),
      // and sizes the table of 4-byte ones in 32 bits.
      if (num8 != 0) {
        malformed("8-byte id offsets, which RDKit misreads");
      }
      if ((num4 + 1) * 4 > kMaxU32 || (num4 + 1) * 4 > size - 8) {
        malformed("an id table larger than its chunk");
      }
      std::uint64_t table = pos + size - (num4 + 1) * 4;
      std::uint64_t prev = 0;
      for (std::uint64_t i = 0; i <= num4; ++i) {
        std::uint64_t offset = read_le<std::uint32_t>(data, table + 4 * i);
        if (offset < prev || offset > size) {
          malformed("id offsets out of order or past their chunk");
        }
        prev = offset;
      }
    }
    pos += size;
  }
  return layout;
}

double popcount(const FPBReader &reader,
                const std::vector<std::uint8_t> &query) {
  return CalcBitmapPopcount(query.data(), reader.layout().storage_size);
}

// The popcount bounds RDKit's similarity searches compute (equations 24 and
// 25 of Swamidass and Baldi), or false when they would index past the
// popcount table. Searching without the screen then gives the same results.
bool tanimoto_screen(const FPBReader &reader,
                     const std::vector<std::uint8_t> &query, double threshold) {
  double probe = popcount(reader, query);
  // RDKit's bound for thresholds this small is a byte count, not a bit
  // count, so it would miss neighbors; scan everything instead.
  if (!(threshold > 1e-6)) {
    return false;
  }
  return reader.screen_in_range(std::floor(threshold * probe),
                                std::ceil(probe / threshold));
}

bool tversky_screen(const FPBReader &reader,
                    const std::vector<std::uint8_t> &query, double ca,
                    double cb, double threshold) {
  double probe = popcount(reader, query);
  if (!(threshold * cb > 1e-6)) {
    return false;
  }
  return reader.screen_in_range(
      std::floor((threshold * probe * ca) / (1. - threshold + threshold * ca)),
      std::ceil(probe * (1 - threshold + threshold * cb) / (threshold * cb)));
}

std::vector<std::pair<double, unsigned int>> tanimoto_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv,
    double threshold) {
  auto query = reader.query(bv);
  return reader.rdkit().getTanimotoNeighbors(
      query.data(), threshold, tanimoto_screen(reader, query, threshold));
}

std::vector<std::pair<double, unsigned int>> tversky_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv, double ca,
    double cb, double threshold) {
  auto query = reader.query(bv);
  return reader.rdkit().getTverskyNeighbors(
      query.data(), ca, cb, threshold,
      tversky_screen(reader, query, ca, cb, threshold));
}

// The containment search always screens, with the query's popcount as the
// index, which the padding keeps at most the number of bits.
std::vector<unsigned int> containing_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv) {
  auto query = reader.query(bv);
  return reader.rdkit().getContainingNeighbors(query.data());
}

rust::Vec<FPBNeighbor> to_neighbors(
    const std::vector<std::pair<double, unsigned int>> &nbrs) {
  rust::Vec<FPBNeighbor> result;
  for (const auto &[sim, idx] : nbrs) {
    result.push_back(FPBNeighbor{sim, idx});
  }
  return result;
}

// As generic_nbr_helper and get_containing_nbrs in MultiFPBReader.cpp: reader
// i is searched by thread i % num_threads, each reader by one thread.
template <class R, class F>
std::vector<R> search_readers(const MultiFPBReader &multi,
                              std::uint32_t num_threads, F search) {
  multi.check_searchable();
  const auto &readers = multi.readers();
  std::vector<R> results(readers.size());
  auto work = [&](unsigned int tid, unsigned int n) {
    for (std::size_t i = tid; i < readers.size(); i += n) {
      if (multi.init_on_search()) {
        readers[i]->init();
      }
      results[i] = search(*readers[i]);
    }
  };
  auto n = RDKit::getNumThreadsToUse(static_cast<int>(
      std::min<std::uint32_t>(num_threads, std::numeric_limits<int>::max())));
  if (n == 1) {
    work(0, 1);
  } else {
    std::vector<std::future<void>> tasks;
    for (unsigned int tid = 0; tid < n && tid < readers.size(); ++tid) {
      tasks.push_back(std::async(std::launch::async, work, tid, n));
    }
    for (auto &task : tasks) {
      task.get();
    }
  }
  return results;
}

// Sorted as MultiFPBReader.cpp's tplSorter: by similarity, then reader, then
// index.
rust::Vec<MultiFPBNeighbor> merge_neighbors(
    const std::vector<std::vector<std::pair<double, unsigned int>>> &results) {
  std::vector<std::tuple<double, unsigned int, unsigned int>> all;
  for (unsigned int r = 0; r < results.size(); ++r) {
    for (const auto &[sim, idx] : results[r]) {
      all.emplace_back(sim, idx, r);
    }
  }
  std::sort(all.begin(), all.end(), [](const auto &a, const auto &b) {
    if (std::get<0>(a) != std::get<0>(b)) {
      return std::get<0>(a) > std::get<0>(b);
    }
    if (std::get<2>(a) != std::get<2>(b)) {
      return std::get<2>(a) < std::get<2>(b);
    }
    return std::get<1>(a) < std::get<1>(b);
  });
  rust::Vec<MultiFPBNeighbor> result;
  for (const auto &[sim, idx, r] : all) {
    result.push_back(MultiFPBNeighbor{sim, idx, r});
  }
  return result;
}

}  // namespace

FPBReader::FPBReader(std::string data, bool lazy)
    : d_layout(check_fpb(data)), d_lazy(lazy) {
  d_stream = new std::istringstream(std::move(data));
  d_reader = std::make_unique<RDKit::FPBReader>(d_stream, true, lazy);
}

void FPBReader::init() {
  d_reader->init();
  // Readers that aren't lazy never read their stream again.
  if (!d_lazy) {
    d_stream->str(std::string());
  }
}

std::vector<std::uint8_t> FPBReader::query(
    rust::Slice<const std::uint8_t> bytes) const {
  auto num_bytes = rdkit().nBits() / 8;  // throws unless initialized
  if (bytes.size() != num_bytes) {
    throw ValueErrorException(
        "the fingerprint has " + std::to_string(bytes.size()) +
        " bytes; the file's have " + std::to_string(num_bytes));
  }
  std::vector<std::uint8_t> query(d_layout.storage_size, 0);
  std::copy(bytes.begin(), bytes.end(), query.begin());
  return query;
}

bool FPBReader::screen_in_range(double min_count, double max_count) const {
  // RDKit screens only with a count for each popcount, 0 to nBits, plus one.
  double num_bits = rdkit().nBits();
  if (d_layout.num_pop_counts != num_bits + 2) {
    return true;  // RDKit won't screen
  }
  return min_count >= 0 && min_count <= num_bits + 1 && max_count >= 0 &&
         max_count <= num_bits;
}

std::uint32_t MultiFPBReader::add_reader(FPBReader &reader) {
  d_readers.push_back(&reader);
  return d_readers.size();
}

// As MultiFPBReader::init.
void MultiFPBReader::init() {
  unsigned int num_bits = 0;
  for (auto *reader : d_readers) {
    reader->init();
    if (!num_bits) {
      num_bits = reader->rdkit().nBits();
    } else if (reader->rdkit().nBits() != num_bits) {
      throw ValueErrorException("bit lengths of child readers don't match");
    }
  }
  d_init = true;
}

std::uint32_t MultiFPBReader::num_bits() const {
  if (d_readers.empty()) {
    throw std::runtime_error("no readers");
  }
  if (!d_init) {
    throw std::runtime_error("not initialized");
  }
  return d_readers[0]->rdkit().nBits();
}

void MultiFPBReader::check_searchable() const {
  if (!d_init && !d_init_on_search) {
    throw std::runtime_error("not initialized");
  }
}

std::unique_ptr<FPBReader> fpb_reader_new(rust::Slice<const std::uint8_t> path,
                                          bool lazy) {
  // As FPBReader::_initFromFilename, which also opens the file eagerly.
  auto name = to_std(path);
  std::ifstream file(name, std::ios_base::binary);
  if (!file || file.bad()) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  std::string data((std::istreambuf_iterator<char>(file)),
                   std::istreambuf_iterator<char>());
  if (file.bad()) {
    throw RDKit::BadFileException("Bad input file " + name);
  }
  return std::make_unique<FPBReader>(std::move(data), lazy);
}

void fpb_reader_init(FPBReader &reader) { reader.init(); }

std::uint32_t fpb_reader_length(const FPBReader &reader) {
  return reader.rdkit().length();
}

std::uint32_t fpb_reader_get_num_bits(const FPBReader &reader) {
  return reader.rdkit().nBits();
}

std::unique_ptr<BitVect> fpb_reader_get_fp(const FPBReader &reader,
                                           std::uint32_t idx) {
  check_index(idx, reader.rdkit().length());
  auto fp = reader.rdkit().getFP(idx);
  if (!fp) {
    throw ValueErrorException("null fingerprint pointer");
  }
  return std::make_unique<BitVect>(*fp);
}

// getBytes copies storage_size bytes, at least num_bytes (check_fpb).
rust::Vec<std::uint8_t> fpb_reader_get_bytes(const FPBReader &reader,
                                             std::uint32_t which) {
  check_index(which, reader.rdkit().length());
  auto bv = reader.rdkit().getBytes(which);
  rust::Vec<std::uint8_t> result;
  for (unsigned i = 0; i < reader.layout().num_bytes; ++i) {
    result.push_back(bv[i]);
  }
  return result;
}

std::unique_ptr<std::string> fpb_reader_get_id(const FPBReader &reader,
                                               std::uint32_t idx) {
  check_index(idx, reader.rdkit().length());
  return std::make_unique<std::string>(reader.rdkit().getId(idx));
}

double fpb_reader_get_tanimoto(const FPBReader &reader, std::uint32_t which,
                               rust::Slice<const std::uint8_t> bytes) {
  check_index(which, reader.rdkit().length());
  auto query = reader.query(bytes);
  return reader.rdkit().getTanimoto(which, query.data());
}

double fpb_reader_get_tversky(const FPBReader &reader, std::uint32_t which,
                              rust::Slice<const std::uint8_t> bytes, double ca,
                              double cb) {
  check_index(which, reader.rdkit().length());
  auto query = reader.query(bytes);
  return reader.rdkit().getTversky(which, query.data(), ca, cb);
}

rust::Vec<FPBNeighbor> fpb_reader_get_tanimoto_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv,
    double threshold) {
  return to_neighbors(tanimoto_neighbors(reader, bv, threshold));
}

rust::Vec<FPBNeighbor> fpb_reader_get_tversky_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv, double ca,
    double cb, double threshold) {
  return to_neighbors(tversky_neighbors(reader, bv, ca, cb, threshold));
}

rust::Vec<std::uint32_t> fpb_reader_get_containing_neighbors(
    const FPBReader &reader, rust::Slice<const std::uint8_t> bv) {
  return to_vec(containing_neighbors(reader, bv));
}

std::unique_ptr<MultiFPBReader> multi_fpb_reader_new(bool init_on_search) {
  return std::make_unique<MultiFPBReader>(init_on_search);
}

std::uint32_t multi_fpb_reader_add_reader(MultiFPBReader &multi,
                                          FPBReader &reader) {
  return multi.add_reader(reader);
}

void multi_fpb_reader_init(MultiFPBReader &multi) { multi.init(); }

std::uint32_t multi_fpb_reader_get_num_bits(const MultiFPBReader &multi) {
  return multi.num_bits();
}

rust::Vec<MultiFPBNeighbor> multi_fpb_reader_get_tanimoto_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv,
    double threshold, std::uint32_t num_threads) {
  using Neighbors = std::vector<std::pair<double, unsigned int>>;
  return merge_neighbors(search_readers<Neighbors>(
      multi, num_threads, [&](const FPBReader &reader) {
        return tanimoto_neighbors(reader, bv, threshold);
      }));
}

rust::Vec<MultiFPBNeighbor> multi_fpb_reader_get_tversky_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv, double ca,
    double cb, double threshold, std::uint32_t num_threads) {
  using Neighbors = std::vector<std::pair<double, unsigned int>>;
  return merge_neighbors(search_readers<Neighbors>(
      multi, num_threads, [&](const FPBReader &reader) {
        return tversky_neighbors(reader, bv, ca, cb, threshold);
      }));
}

// Sorted as MultiFPBReader.cpp's pairSorter: by index, then reader.
rust::Vec<MultiFPBContaining> multi_fpb_reader_get_containing_neighbors(
    const MultiFPBReader &multi, rust::Slice<const std::uint8_t> bv,
    std::uint32_t num_threads) {
  auto results = search_readers<std::vector<unsigned int>>(
      multi, num_threads, [&](const FPBReader &reader) {
        return containing_neighbors(reader, bv);
      });
  std::vector<std::pair<unsigned int, unsigned int>> all;
  for (unsigned int r = 0; r < results.size(); ++r) {
    for (auto idx : results[r]) {
      all.emplace_back(idx, r);
    }
  }
  std::sort(all.begin(), all.end());
  rust::Vec<MultiFPBContaining> result;
  for (const auto &[idx, r] : all) {
    result.push_back(MultiFPBContaining{idx, r});
  }
  return result;
}

}  // namespace rdkit_cxx
