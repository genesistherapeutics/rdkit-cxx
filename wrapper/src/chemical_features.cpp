// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/chemical_features.h"

#include <cstring>

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

Point3D to_point(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

RDGeom::Point3D from_point(const Point3D &p) {
  return RDGeom::Point3D(p.x, p.y, p.z);
}

std::uint32_t read_u32(const std::string &pickle, std::size_t at) {
  std::uint32_t value = 0;
  for (int i = 3; i >= 0; --i) {
    value = (value << 8) | static_cast<unsigned char>(pickle[at + i]);
  }
  return value;
}

[[noreturn]] void bad_pickle(const std::string &why) {
  throw ValueErrorException("bad FreeChemicalFeature pickle: " + why);
}

// FreeChemicalFeature::initFromString (ChemicalFeatures/
// FreeChemicalFeature.cpp) reads, little-endian, a version (0x10, or 0x20
// with an int id after it), then the family and the type, each a uint32
// length and that many bytes, then three doubles. It throws a const char *
// for another version, which no std::exception handler catches. It copies
// each string into a buffer of the length the pickle gives, reads into it
// without checking how much it read, and takes the string up to the first
// NUL: past the end of the buffer when the pickle has none in those bytes
// (or the length is 0), and through bytes it never read when the pickle
// ends early. So each string must be there whole and hold a NUL (CONVENTIONS.md
// U2).
void check_pickle(const std::string &pickle) {
  std::size_t at = 0;
  auto need = [&](std::size_t n, const char *what) {
    if (pickle.size() - at < n) {
      bad_pickle(std::string("it ends before the ") + what);
    }
  };
  need(4, "version");
  auto version = read_u32(pickle, at);
  at += 4;
  if (version == 0x20) {
    need(4, "id");
    at += 4;
  } else if (version != 0x10) {
    bad_pickle("unknown version " + std::to_string(version));
  }
  for (const char *what : {"family", "type"}) {
    need(4, what);
    auto len = read_u32(pickle, at);
    at += 4;
    need(len, what);
    if (!std::memchr(pickle.data() + at, '\0', len)) {
      bad_pickle(std::string("the ") + what + " is not NUL-terminated");
    }
    at += len;
  }
  need(3 * sizeof(double), "position");
}

}  // namespace

std::unique_ptr<FreeChemicalFeature> free_chemical_feature_new() {
  return std::make_unique<FreeChemicalFeature>();
}

std::unique_ptr<FreeChemicalFeature> free_chemical_feature_from_parts(
    rust::Str family, rust::Str type, const Point3D &loc, std::int32_t id) {
  return std::make_unique<FreeChemicalFeature>(
      std::string(family), std::string(type), from_point(loc), id);
}

std::unique_ptr<FreeChemicalFeature> free_chemical_feature_from_pickle(
    rust::Slice<const std::uint8_t> pickle) {
  auto text = to_std(pickle);
  check_pickle(text);
  return std::make_unique<FreeChemicalFeature>(text);
}

std::unique_ptr<FreeChemicalFeature> free_chemical_feature_clone(
    const FreeChemicalFeature &feat) {
  return std::make_unique<FreeChemicalFeature>(feat);
}

rust::Vec<std::uint8_t> free_chemical_feature_to_pickle(
    const FreeChemicalFeature &feat) {
  return to_bytes(feat.toString());
}

std::int32_t free_chemical_feature_get_id(const FreeChemicalFeature &feat) {
  return feat.getId();
}

std::unique_ptr<std::string> free_chemical_feature_get_family(
    const FreeChemicalFeature &feat) {
  return std::make_unique<std::string>(feat.getFamily());
}

std::unique_ptr<std::string> free_chemical_feature_get_type(
    const FreeChemicalFeature &feat) {
  return std::make_unique<std::string>(feat.getType());
}

Point3D free_chemical_feature_get_pos(const FreeChemicalFeature &feat) {
  return to_point(feat.getPos());
}

void free_chemical_feature_set_id(FreeChemicalFeature &feat, std::int32_t id) {
  feat.setId(id);
}

void free_chemical_feature_set_family(FreeChemicalFeature &feat,
                                      rust::Str family) {
  feat.setFamily(std::string(family));
}

void free_chemical_feature_set_type(FreeChemicalFeature &feat, rust::Str type) {
  feat.setType(std::string(type));
}

void free_chemical_feature_set_pos(FreeChemicalFeature &feat,
                                   const Point3D &loc) {
  feat.setPos(from_point(loc));
}

}  // namespace rdkit_cxx
