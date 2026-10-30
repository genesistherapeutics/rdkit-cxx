// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/png.h"

#include <cstring>
#include <map>
#include <mutex>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string/predicate.hpp>

#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/molfiles.rs.h"
#include "rdkit-cxx/ffi/png.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/molfiles.h"
#include "rdkit-cxx/streams.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace FP = RDKit::v2::FileParsers;

std::uint32_t big_endian(const char *bytes) {
  auto b = reinterpret_cast<const unsigned char *>(bytes);
  return (std::uint32_t{b[0]} << 24) | (std::uint32_t{b[1]} << 16) |
         (std::uint32_t{b[2]} << 8) | std::uint32_t{b[3]};
}

// RDKit's PNG readers (PNGStreamToMetadata and addMetadataToPNGStream in
// PNGParser.cpp) walk the chunks without checking their lengths:
// addMetadataToPNGStream reads the length and type past the end of a file
// without an IEND chunk as uninitialized values and sizes a buffer by them,
// and both read a text chunk's keyword up to a NUL past the chunk. This
// checks that the bytes are whole chunks up to IEND and that each text chunk
// RDKit reads holds its keyword's NUL.
void check_png(const std::string &png) {
  static const char header[8] = {'\x89', 'P',  'N',    'G',
                                 '\r',   '\n', '\x1a', '\n'};
  if (png.size() < 8 || std::memcmp(png.data(), header, 8) != 0) {
    throw RDKit::FileParseException("PNG header not recognized");
  }
  std::size_t pos = 8;
  for (;;) {
    if (png.size() - pos < 8) {
      throw RDKit::FileParseException("did not find IEND block in PNG");
    }
    std::uint64_t len = big_endian(png.data() + pos);
    const char *type = png.data() + pos + 4;
    pos += 8;
    if (png.size() - pos < len + 4) {
      throw RDKit::FileParseException(
          "a PNG chunk runs past the end of the data");
    }
    if (std::memcmp(type, "IEND", 4) == 0) {
      return;
    }
    // As RDKit tells text chunks (tEXt and zTXt) apart.
    bool text = ((type[0] == 't' && type[1] == 'E') ||
                 (type[0] == 'z' && type[1] == 'T')) &&
                type[2] == 'X' && type[3] == 't';
    if (len > 0 && text && !std::memchr(png.data() + pos, '\0', len)) {
      throw RDKit::FileParseException("a PNG text chunk has no keyword");
    }
    pos += len + 4;
  }
}

// As Mol::from_pkl_string, a pickle's monomer types are checked (U10).
std::unique_ptr<RWMol> from_pickle(const std::string &pkl) {
  auto mol = std::make_unique<RWMol>();
  RDKit::MolPickler::molFromPickle(pkl, mol.get());
  check_pickled_mol(*mol);
  return mol;
}

// As PNGStreamToMols in PNGParser.cpp: a molecule for one metadata value,
// or null for a SMILES RDKit cannot parse. The mol block is checked (U12),
// and so is a SMILES's CXSMILES extension (check_smiles_cx).
std::unique_ptr<RWMol> from_metadata(const std::string &key,
                                     const std::string &value,
                                     const RDKit::SmilesParserParams &params) {
  if (boost::starts_with(key, RDKit::PNGData::pklTag)) {
    return from_pickle(value);
  }
  if (boost::starts_with(key, RDKit::PNGData::smilesTag)) {
    check_smiles_cx(value, params.allowCXSMILES, params.strictCXSMILES);
    return std::unique_ptr<RWMol>(RDKit::SmilesToMol(value, params));
  }
  check_input_text(value);
  FP::MolFileParserParams mol_params;
  mol_params.sanitize = params.sanitize;
  mol_params.removeHs = params.removeHs;
  return FP::MolFromMolBlock(value, mol_params);
}

bool has_mol_tag(const std::string &key) {
  return boost::starts_with(key, RDKit::PNGData::pklTag) ||
         boost::starts_with(key, RDKit::PNGData::smilesTag) ||
         boost::starts_with(key, RDKit::PNGData::molTag);
}

// As PNGStreamToMol in PNGParser.cpp: the first record RDKit can use.
std::unique_ptr<RWMol> png_to_mol(const std::string &png,
                                  const SmilesParserParamsData &data) {
  check_png(png);
  std::map<std::string, std::string> repl;
  auto params = smiles_parser_params(data, repl);
  for (const auto &[key, value] : RDKit::PNGStringToMetadata(png)) {
    if (has_mol_tag(key)) {
      auto mol = from_metadata(key, value, params);
      if (!mol) {
        throw RDKit::FileParseException("could not parse the PNG's " + key);
      }
      return mol;
    }
  }
  throw RDKit::FileParseException("No suitable metadata found.");
}

std::unique_ptr<MolList> png_to_mols(const std::string &png,
                                     const std::string &tag,
                                     const SmilesParserParamsData &data) {
  check_png(png);
  std::map<std::string, std::string> repl;
  auto params = smiles_parser_params(data, repl);
  auto result = std::make_unique<MolList>();
  for (const auto &[key, value] : RDKit::PNGStringToMetadata(png)) {
    if (boost::starts_with(key, tag) && has_mol_tag(key)) {
      result->push(from_metadata(key, value, params));
    }
  }
  return result;
}

// RDKit 2024.03 pickles with the default properties and writes the CXSMILES
// and mol block with the writers' defaults (addMolToPNGStream).
rust::Vec<std::uint8_t> mol_to_png(const RWMol &mol, const std::string &png,
                                   const PNGMetadataParamsData &data) {
  check_png(png);
  if (data.include_pkl) {
    check_picklable(mol);
  }
  if (data.include_smiles) {
    // The SMILES writer removes atoms from copies (mol_to_smiles).
    check_removals(mol);
  }
  // The mol block computes 2D coordinates for a molecule without a
  // conformer, with compute2DCoords, which trusts the molecule
  // (depiction_checks.h).
  std::unique_lock<std::mutex> depicting;
  if (data.include_mol && mol.getNumConformers() == 0) {
    depicting = lock_depiction();
    check_compute_2d_coords_mol(mol, false);
  }
  return to_bytes(RDKit::addMolToPNGString(
      mol, png, data.include_pkl, data.include_smiles, data.include_mol));
}

rust::Vec<std::uint8_t> add_metadata(const rust::Vec<rust::String> &keys,
                                     const rust::Vec<rust::String> &values,
                                     const std::string &png) {
  check_png(png);
  std::vector<std::pair<std::string, std::string>> metadata;
  for (std::size_t i = 0; i < keys.size() && i < values.size(); ++i) {
    metadata.emplace_back(std::string(keys[i]), std::string(values[i]));
  }
  return to_bytes(RDKit::addMetadataToPNGString(png, metadata));
}

std::unique_ptr<PNGMetadata> metadata(const std::string &png) {
  check_png(png);
  auto result = std::make_unique<PNGMetadata>();
  result->items = RDKit::PNGStringToMetadata(png);
  return result;
}

}  // namespace

std::unique_ptr<RWMol> mol_from_png_string(
    rust::Slice<const std::uint8_t> png, const SmilesParserParamsData &params) {
  return png_to_mol(to_std(png), params);
}

// Each file form reads the file once and checks those bytes (U2).
std::unique_ptr<RWMol> mol_from_png_file(rust::Slice<const std::uint8_t> path,
                                         const SmilesParserParamsData &params) {
  return png_to_mol(read_file(to_path(path)), params);
}

std::unique_ptr<MolList> mols_from_png_string(
    rust::Slice<const std::uint8_t> png, rust::Str tag,
    const SmilesParserParamsData &params) {
  return png_to_mols(to_std(png), std::string(tag), params);
}

std::unique_ptr<MolList> mols_from_png_file(
    rust::Slice<const std::uint8_t> path, rust::Str tag,
    const SmilesParserParamsData &params) {
  return png_to_mols(read_file(to_path(path)), std::string(tag), params);
}

rust::Vec<std::uint8_t> mol_metadata_to_png_string(
    const RWMol &mol, rust::Slice<const std::uint8_t> png,
    const PNGMetadataParamsData &params) {
  return mol_to_png(mol, to_std(png), params);
}

rust::Vec<std::uint8_t> mol_metadata_to_png_file(
    const RWMol &mol, rust::Slice<const std::uint8_t> path,
    const PNGMetadataParamsData &params) {
  return mol_to_png(mol, read_file(to_path(path)), params);
}

rust::Vec<std::uint8_t> add_metadata_to_png_string(
    const rust::Vec<rust::String> &keys, const rust::Vec<rust::String> &values,
    rust::Slice<const std::uint8_t> png) {
  return add_metadata(keys, values, to_std(png));
}

rust::Vec<std::uint8_t> add_metadata_to_png_file(
    const rust::Vec<rust::String> &keys, const rust::Vec<rust::String> &values,
    rust::Slice<const std::uint8_t> path) {
  return add_metadata(keys, values, read_file(to_path(path)));
}

std::unique_ptr<PNGMetadata> metadata_from_png_string(
    rust::Slice<const std::uint8_t> png) {
  return metadata(to_std(png));
}

std::unique_ptr<PNGMetadata> metadata_from_png_file(
    rust::Slice<const std::uint8_t> path) {
  return metadata(read_file(to_path(path)));
}

std::size_t png_metadata_size(const PNGMetadata &metadata) {
  return metadata.items.size();
}

const std::string &png_metadata_key(const PNGMetadata &metadata,
                                    std::size_t i) {
  return metadata.items.at(i).first;
}

rust::Vec<std::uint8_t> png_metadata_value(const PNGMetadata &metadata,
                                           std::size_t i) {
  return to_bytes(metadata.items.at(i).second);
}

}  // namespace rdkit_cxx
