// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/molfiles.h"

#include <map>
#include <mutex>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/FileParseException.h>

#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/ffi/molfiles.rs.h"
#include "rdkit-cxx/random_lock.h"

namespace rdkit_cxx {

namespace {

std::string to_std(rust::Slice<const std::uint8_t> bytes) {
  return std::string(reinterpret_cast<const char *>(bytes.data()),
                     bytes.size());
}

std::map<std::string, std::string> replacements(
    const rust::Vec<rust::String> &keys,
    const rust::Vec<rust::String> &values) {
  std::map<std::string, std::string> result;
  for (std::size_t i = 0; i < keys.size() && i < values.size(); ++i) {
    result[std::string(keys[i])] = std::string(values[i]);
  }
  return result;
}

RDKit::SmilesWriteParams write_params(const SmilesWriteParamsData &data) {
  RDKit::SmilesWriteParams params;
  params.doIsomericSmiles = data.do_isomeric_smiles;
  params.doKekule = data.do_kekule;
  params.canonical = data.canonical;
  params.cleanStereo = data.clean_stereo;
  params.allBondsExplicit = data.all_bonds_explicit;
  params.allHsExplicit = data.all_hs_explicit;
  params.doRandom = data.do_random;
  params.rootedAtAtom = data.rooted_at_atom;
  params.includeDativeBonds = data.include_dative_bonds;
  params.ignoreAtomMapNumbers = data.ignore_atom_map_numbers;
  return params;
}

RDKit::MolWriterParams writer_params(const MolWriterParamsData &data) {
  RDKit::MolWriterParams params;
  params.includeStereo = data.include_stereo;
  params.kekulize = data.kekulize;
  params.forceV3000 = data.force_v3000;
  params.precision = data.precision;
  return params;
}

// Python returns None where these return null (E2).
std::unique_ptr<RWMol> parsed(RDKit::RWMol *mol, const std::string &what,
                              const std::string &text) {
  if (!mol) {
    throw RDKit::SmilesParse::SmilesParseException("could not parse " + what +
                                                   ": " + text);
  }
  return std::unique_ptr<RWMol>(mol);
}

}  // namespace

SmilesParserParamsData smiles_parser_params_default() {
  RDKit::SmilesParserParams params;
  return SmilesParserParamsData{params.debugParse,
                                params.parseName,
                                params.allowCXSMILES,
                                params.strictCXSMILES,
                                params.sanitize,
                                params.removeHs,
                                {},
                                {}};
}

SmartsParserParamsData smarts_parser_params_default() {
  RDKit::SmartsParserParams params;
  return SmartsParserParamsData{params.debugParse,
                                params.parseName,
                                params.allowCXSMILES,
                                params.strictCXSMILES,
                                params.mergeHs,
                                {},
                                {}};
}

SmilesWriteParamsData smiles_write_params_default() {
  RDKit::SmilesWriteParams params;
  return SmilesWriteParamsData{
      params.doIsomericSmiles,   params.doKekule,
      params.canonical,          params.cleanStereo,
      params.allBondsExplicit,   params.allHsExplicit,
      params.doRandom,           params.rootedAtAtom,
      params.includeDativeBonds, params.ignoreAtomMapNumbers};
}

MolWriterParamsData mol_writer_params_default() {
  RDKit::MolWriterParams params;
  return MolWriterParamsData{params.includeStereo, params.kekulize,
                             params.forceV3000, params.precision, -1};
}

// As MolFromSmilesHelper in Code/GraphMol/nbWrap/rdmolfiles.cpp.
std::unique_ptr<RWMol> mol_from_smiles(rust::Str smiles,
                                       const SmilesParserParamsData &data) {
  std::string text(smiles);
  auto repl = replacements(data.replacement_keys, data.replacement_values);
  RDKit::SmilesParserParams params;
  params.debugParse = data.debug_parse;
  params.parseName = data.parse_name;
  params.allowCXSMILES = data.allow_cxsmiles;
  params.strictCXSMILES = data.strict_cxsmiles;
  params.sanitize = data.sanitize;
  params.removeHs = data.remove_hs;
  params.replacements = &repl;
  return parsed(RDKit::SmilesToMol(text, params), "SMILES", text);
}

std::unique_ptr<std::string> mol_to_smiles(const RWMol &mol,
                                           const SmilesWriteParamsData &data) {
  // Random SMILES draw from RDKit's global generator (Canon.cpp,
  // SmilesWrite.cpp). MolToSmarts ignores doRandom.
  std::unique_lock<std::mutex> lock;
  if (data.do_random) {
    lock = lock_random_generator();
  }
  return std::make_unique<std::string>(
      RDKit::MolToSmiles(mol, write_params(data)));
}

// As MolFromSmartsHelper in Code/GraphMol/nbWrap/rdmolfiles.cpp.
std::unique_ptr<RWMol> mol_from_smarts(rust::Str smarts,
                                       const SmartsParserParamsData &data) {
  std::string text(smarts);
  auto repl = replacements(data.replacement_keys, data.replacement_values);
  RDKit::SmartsParserParams params;
  params.debugParse = data.debug_parse;
  params.parseName = data.parse_name;
  params.allowCXSMILES = data.allow_cxsmiles;
  params.strictCXSMILES = data.strict_cxsmiles;
  params.mergeHs = data.merge_hs;
  params.replacements = &repl;
  return parsed(RDKit::SmartsToMol(text, params), "SMARTS", text);
}

std::unique_ptr<std::string> mol_to_smarts(const RWMol &mol,
                                           const SmilesWriteParamsData &data) {
  return std::make_unique<std::string>(
      RDKit::MolToSmarts(mol, write_params(data)));
}

std::unique_ptr<RWMol> mol_from_mol_block(rust::Str mol_block, bool sanitize,
                                          bool remove_hs, bool strict_parsing) {
  RDKit::v2::FileParsers::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.strictParsing = strict_parsing;
  auto mol =
      RDKit::v2::FileParsers::MolFromMolBlock(std::string(mol_block), params);
  if (!mol) {
    throw RDKit::FileParseException("could not parse the mol block");
  }
  return mol;
}

std::unique_ptr<RWMol> mol_from_mol_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs,
                                         bool strict_parsing) {
  RDKit::v2::FileParsers::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.strictParsing = strict_parsing;
  auto name = to_std(path);
  auto mol = RDKit::v2::FileParsers::MolFromMolFile(name, params);
  if (!mol) {
    throw RDKit::FileParseException("could not parse the mol file " + name);
  }
  return mol;
}

namespace {

// The mol-block writers compute 2D coordinates for stereo when the molecule
// has no conformer (prepareMol in MolFileWriter.cpp).
std::unique_lock<std::mutex> lock_if_depicting(
    const RWMol &mol, const MolWriterParamsData &data) {
  if (data.include_stereo && mol.getNumConformers() == 0) {
    return lock_depiction();
  }
  return {};
}

}  // namespace

std::unique_ptr<std::string> mol_to_mol_block(const RWMol &mol,
                                              const MolWriterParamsData &data) {
  auto lock = lock_if_depicting(mol, data);
  return std::make_unique<std::string>(
      RDKit::MolToMolBlock(mol, writer_params(data), data.conf_id));
}

void mol_to_mol_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     const MolWriterParamsData &data) {
  auto lock = lock_if_depicting(mol, data);
  RDKit::MolToMolFile(mol, to_std(path), writer_params(data), data.conf_id);
}

}  // namespace rdkit_cxx
