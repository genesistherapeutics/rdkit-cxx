// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mhfp.h"

#include <climits>
#include <string>
#include <vector>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/mhfp.rs.h"
#include "rdkit-cxx/file_checks.h"

namespace rdkit_cxx {

namespace {

// RDKit declares the encoder's methods non-const, but they only read its
// permutations, so a shared encoder can call them.
MHFPEncoder &encoder_mut(const MHFPEncoder &encoder) {
  return const_cast<MHFPEncoder &>(encoder);
}

// CreateShingling(smiles) asserts that the SMILES parsed; a parse error here
// is a ParseError instead (E2).
std::unique_ptr<RWMol> parse(rust::Str smiles) {
  check_smiles_cx(std::string(smiles));
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(std::string(smiles)));
  if (!mol) {
    throw RDKit::SmilesParse::SmilesParseException("could not parse SMILES '" +
                                                   std::string(smiles) + "'");
  }
  return std::make_unique<RWMol>(*mol);
}

std::vector<std::string> shingling(const MHFPEncoder &encoder,
                                   const RDKit::ROMol &mol,
                                   const ShinglingOptions &options) {
  return encoder_mut(encoder).CreateShingling(
      mol, options.radius, options.rings, options.isomeric, options.kekulize,
      options.min_radius);
}

rust::Vec<rust::String> to_strings(const std::vector<std::string> &values) {
  rust::Vec<rust::String> result;
  result.reserve(values.size());
  for (const auto &value : values) {
    // Shingles are SMILES, which RDKit writes in ASCII.
    result.push_back(rust::String(value));
  }
  return result;
}

rust::Vec<std::uint32_t> to_vec(const std::vector<std::uint32_t> &values) {
  rust::Vec<std::uint32_t> result;
  result.reserve(values.size());
  for (auto value : values) {
    result.push_back(value);
  }
  return result;
}

// EncodeSECFP folds into length bits with length taken modulo 2^32: zero
// divides by zero, and the bit vector holds at most INT_MAX bits (U3).
void check_length(std::size_t length) {
  if (length == 0 || length > static_cast<std::size_t>(INT_MAX)) {
    throw ValueErrorException("length must be between 1 and 2^31 - 1");
  }
}

std::unique_ptr<BitVect> secfp(const MHFPEncoder &encoder, RDKit::ROMol &mol,
                               const ShinglingOptions &options,
                               std::size_t length) {
  check_length(length);
  return wrap_bit_vect(
      std::make_unique<ExplicitBitVect>(encoder_mut(encoder).EncodeSECFP(
          mol, options.radius, options.rings, options.isomeric,
          options.kekulize, options.min_radius, length)));
}

}  // namespace

std::unique_ptr<MHFPEncoder> mhfp_encoder_new(std::uint32_t n_permutations,
                                              std::uint32_t seed) {
  return std::make_unique<MHFPEncoder>(n_permutations, seed);
}

rust::Vec<std::uint32_t> mhfp_from_string_array(
    const MHFPEncoder &encoder, rust::Slice<const rust::Str> vec) {
  std::vector<std::string> strings;
  strings.reserve(vec.size());
  for (auto s : vec) {
    strings.emplace_back(s);
  }
  return to_vec(encoder_mut(encoder).FromStringArray(strings));
}

rust::Vec<std::uint32_t> mhfp_from_array(const MHFPEncoder &encoder,
                                         rust::Slice<const std::uint32_t> vec) {
  return to_vec(encoder_mut(encoder).FromArray(
      std::vector<std::uint32_t>(vec.begin(), vec.end())));
}

rust::Vec<rust::String> mhfp_shingling_from_smiles(
    const MHFPEncoder &encoder, rust::Str smiles,
    const ShinglingOptions &options) {
  return to_strings(shingling(encoder, *parse(smiles), options));
}

rust::Vec<rust::String> mhfp_shingling_from_mol(
    const MHFPEncoder &encoder, const RWMol &mol,
    const ShinglingOptions &options) {
  return to_strings(shingling(encoder, mol, options));
}

rust::Vec<std::uint32_t> mhfp_encode_smiles(const MHFPEncoder &encoder,
                                            rust::Str smiles,
                                            const ShinglingOptions &options) {
  return to_vec(encoder_mut(encoder).FromStringArray(
      shingling(encoder, *parse(smiles), options)));
}

rust::Vec<std::uint32_t> mhfp_encode_mol(const MHFPEncoder &encoder,
                                         const RWMol &mol,
                                         const ShinglingOptions &options) {
  return to_vec(
      encoder_mut(encoder).FromStringArray(shingling(encoder, mol, options)));
}

std::unique_ptr<BitVect> mhfp_encode_secfp_smiles(
    const MHFPEncoder &encoder, rust::Str smiles,
    const ShinglingOptions &options, std::size_t length) {
  return secfp(encoder, *parse(smiles), options, length);
}

std::unique_ptr<BitVect> mhfp_encode_secfp_mol(const MHFPEncoder &encoder,
                                               const RWMol &mol,
                                               const ShinglingOptions &options,
                                               std::size_t length) {
  // As Python's EncodeSECFPMol, which encodes a copy (EncodeSECFP takes a
  // non-const molecule).
  RDKit::ROMol copy(mol);
  return secfp(encoder, copy, options, length);
}

double mhfp_distance(rust::Slice<const std::uint32_t> a,
                     rust::Slice<const std::uint32_t> b) {
  // Distance reads b at every index of a.
  if (a.size() != b.size()) {
    throw ValueErrorException("the MinHashes must be the same length");
  }
  return MHFPEncoder::Distance(std::vector<std::uint32_t>(a.begin(), a.end()),
                               std::vector<std::uint32_t>(b.begin(), b.end()));
}

}  // namespace rdkit_cxx
