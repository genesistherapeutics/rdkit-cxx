/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/pickle_checks.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <optional>
#include <string>

#include <DataStructs/DiscreteValueVect.h>

namespace rdkit_cxx {

namespace {

constexpr long long kMaxInt = std::numeric_limits<int>::max();

// The little-endian value at `offset`, as RDKit's streamRead reads it, or
// nothing when the bytes run out first (and streamRead throws).
template <class T>
std::optional<T> read(std::string_view pkl, std::size_t offset) {
  if (offset > pkl.size() || pkl.size() - offset < sizeof(T)) {
    return std::nullopt;
  }
  unsigned char bytes[sizeof(T)];
  std::memcpy(bytes, pkl.data() + offset, sizeof(T));
  if constexpr (std::endian::native == std::endian::big) {
    std::reverse(bytes, bytes + sizeof(T));
  }
  T value;
  std::memcpy(&value, bytes, sizeof(T));
  return value;
}

// RDKit negates the leading version, which overflows for INT_MIN.
void check_version(std::string_view pkl, const char *name) {
  auto version = read<std::int32_t>(pkl, 0);
  if (version && *version == std::numeric_limits<std::int32_t>::min()) {
    throw ValueErrorException(std::string("bad version in ") + name +
                              " pickle");
  }
}

// RDKit's grids cast their sizes to int and their extents in grid units to
// unsigned int when copied.
void check_grid(std::uint32_t nx, std::uint32_t ny, std::uint32_t nz,
                double spacing, double ox, double oy, double oz,
                std::uint64_t storage_length, const char *name) {
  auto bad = [&](const char *what) {
    throw ValueErrorException(std::string(what) + " in " + name + " pickle");
  };
  if (!(std::isfinite(spacing) && spacing > 0.0)) {
    bad("bad spacing");
  }
  if (!(std::isfinite(ox) && std::isfinite(oy) && std::isfinite(oz))) {
    bad("bad offset");
  }
  if (nx == 0 || ny == 0 || nz == 0) {
    bad("an empty axis");
  }
  for (auto n : {nx, ny, nz}) {
    if (!std::isfinite(n * spacing)) {
      bad("bad extent");
    }
  }
  std::uint64_t points = std::uint64_t{nx} * ny;
  if (points > static_cast<std::uint64_t>(kMaxInt) ||
      points * nz > static_cast<std::uint64_t>(kMaxInt)) {
    bad("too many points");
  }
  if (points * nz != storage_length) {
    bad("a storage size mismatch");
  }
}

}  // namespace

void check_bit_vect_pickle(std::string_view pkl) {
  check_version(pkl, "BitVect");
  // RDKit takes a negative size as a huge unsigned one.
  auto version = read<std::int32_t>(pkl, 0);
  auto size = read<std::int32_t>(pkl, 4);
  if (version && *version < 0 && size && *size < 0) {
    throw ValueErrorException("negative size in BitVect pickle");
  }
}

void check_discrete_value_vect_pickle(std::string_view pkl) {
  check_version(pkl, "DiscreteValueVect");
  // Header: version, type, bits per value, mask, length, number of ints;
  // then the ints. RDKit trusts every field, so check them all.
  auto version = read<std::int32_t>(pkl, 0);
  auto type = read<std::uint32_t>(pkl, 4);
  auto bits = read<std::uint32_t>(pkl, 8);
  auto mask = read<std::uint32_t>(pkl, 12);
  auto length = read<std::uint32_t>(pkl, 16);
  auto num_ints = read<std::uint32_t>(pkl, 20);
  auto bad = [](const char *what) {
    throw ValueErrorException(std::string(what) +
                              " in DiscreteValueVect pickle");
  };
  // RDKit checks none of its reads: it compares a version it may not have
  // read, and divides by the bits per value before reading the rest.
  if (!version) {
    bad("a truncated header");
  }
  if (*version != -1) {
    // RDKit rejects any other version.
    return;
  }
  if (!type || !bits || !mask || !length || !num_ints) {
    bad("a truncated header");
  }
  if (*type > RDKit::DiscreteValueVect::SIXTEENBITVALUE ||
      *bits != (1u << *type) || *mask != (1ull << *bits) - 1) {
    bad("a bad value type");
  }
  if (*length > kMaxInt) {
    bad("a bad length");
  }
  std::uint64_t vals_per_int = 32 / *bits;
  if (*num_ints != (*length + vals_per_int - 1) / vals_per_int) {
    bad("a bad data size");
  }
  if (pkl.size() - 24 < std::uint64_t{*num_ints} * 4) {
    bad("truncated data");
  }
}

void check_real_value_vect_pickle(std::string_view pkl) {
  check_version(pkl, "RealValueVect");
  // Header: version, length; then the values. RDKit allocates the length it
  // reads before reading them, so a short pickle could ask for 32 GB.
  auto version = read<std::int32_t>(pkl, 0);
  auto length = read<std::uint32_t>(pkl, 4);
  if (!version || *version != -1 || !length) {
    return;
  }
  if (*length > kMaxInt) {
    throw ValueErrorException("a bad length in RealValueVect pickle");
  }
  if ((pkl.size() - 8) / sizeof(double) < *length) {
    throw ValueErrorException("truncated data in RealValueVect pickle");
  }
}

void check_uniform_grid3d_pickle(std::string_view pkl) {
  check_version(pkl, "UniformGrid3D");
  // Header: version, nx, ny, nz, spacing, offset x, y, z, storage size; then
  // the storage, a DiscreteValueVect pickle.
  auto version = read<std::int32_t>(pkl, 0);
  auto nx = read<std::uint32_t>(pkl, 4);
  auto ny = read<std::uint32_t>(pkl, 8);
  auto nz = read<std::uint32_t>(pkl, 12);
  auto spacing = read<double>(pkl, 16);
  auto ox = read<double>(pkl, 24);
  auto oy = read<double>(pkl, 32);
  auto oz = read<double>(pkl, 40);
  auto storage_size = read<std::uint32_t>(pkl, 48);
  if (!version || *version != -1 || !storage_size) {
    return;
  }
  // RDKit reads the storage without checking for a short read.
  if (pkl.size() - 52 < *storage_size) {
    throw ValueErrorException("truncated UniformGrid3D pickle");
  }
  auto storage = pkl.substr(52, *storage_size);
  check_discrete_value_vect_pickle(storage);
  auto length = read<std::uint32_t>(storage, 16);
  if (!length) {
    throw ValueErrorException("truncated UniformGrid3D pickle");
  }
  check_grid(*nx, *ny, *nz, *spacing, *ox, *oy, *oz, *length, "UniformGrid3D");
}

void check_uniform_real_value_grid3d_pickle(std::string_view pkl) {
  check_version(pkl, "UniformRealValueGrid3D");
  auto version = read<std::int32_t>(pkl, 0);
  auto nx = read<std::uint32_t>(pkl, 4);
  auto ny = read<std::uint32_t>(pkl, 8);
  auto nz = read<std::uint32_t>(pkl, 12);
  auto spacing = read<double>(pkl, 16);
  auto ox = read<double>(pkl, 24);
  auto oy = read<double>(pkl, 32);
  auto oz = read<double>(pkl, 40);
  auto storage_size = read<std::uint32_t>(pkl, 48);
  if (!version || *version != -1 || !storage_size) {
    return;
  }
  // RDKit allocates the storage size it reads before reading the storage.
  if (pkl.size() - 52 < *storage_size) {
    throw ValueErrorException("truncated UniformRealValueGrid3D pickle");
  }
  auto storage = pkl.substr(52, *storage_size);
  check_real_value_vect_pickle(storage);
  auto length = read<std::uint32_t>(storage, 4);
  if (!length) {
    throw ValueErrorException("truncated UniformRealValueGrid3D pickle");
  }
  check_grid(*nx, *ny, *nz, *spacing, *ox, *oy, *oz, *length,
             "UniformRealValueGrid3D");
}

}  // namespace rdkit_cxx
