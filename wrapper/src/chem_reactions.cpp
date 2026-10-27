// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/chem_reactions.h"

#include <climits>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/FilterCatalog/FunctionalGroupHierarchy.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>
#include <GraphMol/MarvinParse/MarvinParser.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "rdkit-cxx/adjust_query.h"
#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/adjust_query.rs.h"
#include "rdkit-cxx/ffi/chem_reactions.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/molfiles.rs.h"
#include "rdkit-cxx/ffi/substruct.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/molfiles.h"
#include "rdkit-cxx/streams.h"
#include "rdkit-cxx/substruct.h"

namespace rdkit_cxx {

namespace {

using RDKit::MOL_SPTR_VECT;
using RDKit::ROMOL_SPTR;

// The reaction's templates, as check_removals checks a molecule: the SMILES
// writer writes a template of several components a component at a time,
// which it makes by removing the other atoms from copies (getMolFrags), and
// SANITIZE_MERGEHS removes the reactants' query hydrogens.
void check_template_removals(const ChemicalReaction &rxn) {
  for (const auto *templates :
       {&rxn.getReactants(), &rxn.getAgents(), &rxn.getProducts()}) {
    for (const auto &templ : *templates) {
      check_removals(*templ);
    }
  }
}

std::string to_std(rust::Slice<const std::uint8_t> bytes) {
  return std::string(reinterpret_cast<const char *>(bytes.data()),
                     bytes.size());
}

rust::Vec<std::uint8_t> to_vec(const std::string &text) {
  rust::Vec<std::uint8_t> res;
  res.reserve(text.size());
  for (char c : text) {
    res.push_back(static_cast<std::uint8_t>(c));
  }
  return res;
}

std::unique_ptr<std::string> own(std::string text) {
  return std::make_unique<std::string>(std::move(text));
}

// The bytes of a file, read once so RDKit parses exactly what was checked
// (U2).
std::string read_path(rust::Slice<const std::uint8_t> path) {
  return read_file(to_std(path));
}

const MOL_SPTR_VECT &templates(const ChemicalReaction &rxn, std::uint8_t kind) {
  switch (kind) {
    case 0:
      return rxn.getReactants();
    case 1:
      return rxn.getProducts();
    case 2:
      return rxn.getAgents();
    default:
      throw std::logic_error("unknown reaction template kind");
  }
}

void normalize(MOL_SPTR_VECT::iterator begin, MOL_SPTR_VECT::iterator end,
               const char *what) {
  for (auto it = begin; it != end; ++it) {
    // Without this, RDKit dereferences the null template in every later
    // call that reads the templates.
    if (!*it) {
      throw RDKit::ChemicalReactionParserException(
          std::string("the reaction has a null ") + what + " template");
    }
    if (!dynamic_cast<RWMol *>(it->get()) || it->use_count() > 1) {
      auto *mol = new RWMol(**it);
      clear_matrix_caches(*mol);
      *it = ROMOL_SPTR(mol);
    }
  }
}

std::unique_ptr<MolList> copies(const MOL_SPTR_VECT &mols) {
  auto res = std::make_unique<MolList>();
  for (const auto &mol : mols) {
    auto copy = std::make_unique<RWMol>(*mol);
    clear_matrix_caches(*copy);
    res->push(std::move(copy));
  }
  return res;
}

std::unique_ptr<ChemicalReaction> finish(std::unique_ptr<ChemicalReaction> rxn,
                                         const std::string &what) {
  if (!rxn) {
    throw RDKit::ChemicalReactionParserException("could not parse " + what);
  }
  normalize_templates(*rxn);
  return rxn;
}

std::unique_ptr<RWMol> to_rwmol(std::unique_ptr<RDKit::ROMol> mol) {
  std::unique_ptr<RWMol> res;
  if (auto *rw = dynamic_cast<RWMol *>(mol.get())) {
    mol.release();
    res.reset(rw);
  } else {
    res = std::make_unique<RWMol>(*mol);
  }
  clear_matrix_caches(*res);
  return res;
}

// Python's runs initialize the reaction first; a reaction that fails
// validation stays uninitialized, and RDKit refuses to run it.
void initialize_if_needed(ChemicalReaction &rxn) {
  if (!rxn.isInitialized()) {
    rxn.initReactantMatchers();
  }
}

std::map<std::string, ROMOL_SPTR> query_map(BorrowedMols &queries) {
  std::map<std::string, ROMOL_SPTR> res;
  for (std::size_t i = 0; i < queries.mols().size(); ++i) {
    res[queries.names()[i]] = ROMOL_SPTR(new RDKit::ROMol(*queries.mols()[i]));
  }
  return res;
}

rust::Vec<ReactantLabel> to_labels(
    const std::vector<std::vector<std::pair<unsigned int, std::string>>>
        &labels) {
  rust::Vec<ReactantLabel> res;
  for (std::size_t r = 0; r < labels.size(); ++r) {
    for (const auto &[idx, label] : labels[r]) {
      res.push_back(ReactantLabel{static_cast<std::uint32_t>(r), idx,
                                  std::make_unique<std::string>(label)});
    }
  }
  return res;
}

bool any_template_without_conformer(const ChemicalReaction &rxn) {
  for (std::uint8_t kind = 0; kind < 3; ++kind) {
    for (const auto &mol : templates(rxn, kind)) {
      if (!mol->getNumConformers()) {
        return true;
      }
    }
  }
  return false;
}

// The V2000 rxn writer writes each template with MolToMolBlock and
// includeStereo, which computes 2D coordinates for a molecule without a
// conformer (prepareMol in MolFileWriter.cpp) with compute2DCoords, so it
// holds the lock and checks those templates under it (U4); compute2DCoords
// trusts the molecule (depiction_checks.h). The V3000 writer writes the
// templates' CTABs without depicting them.
std::unique_lock<std::mutex> lock_if_depicting(const ChemicalReaction &rxn,
                                               bool force_v3000) {
  if (force_v3000 || !any_template_without_conformer(rxn)) {
    return {};
  }
  auto lock = lock_depiction();
  for (std::uint8_t kind = 0; kind < 3; ++kind) {
    for (const auto &mol : templates(rxn, kind)) {
      if (!mol->getNumConformers()) {
        check_compute_2d_coords_mol(*mol, false);
      }
    }
  }
  return lock;
}

// PNG files: the 8-byte signature, then chunks of a big-endian length, a
// 4-byte tag, the data, and a CRC.
constexpr unsigned char kPngHeader[8] = {137, 80, 78, 71, 13, 10, 26, 10};

std::uint32_t read_be32(const std::string &data, std::size_t pos) {
  std::uint32_t res = 0;
  for (std::size_t i = 0; i < 4; ++i) {
    res = (res << 8) | static_cast<unsigned char>(data[pos + i]);
  }
  return res;
}

// RDKit's checkPNGHeader reads the signature a byte at a time without
// checking that the reads succeed, so a short input that matches so far is
// compared through uninitialized bytes.
void check_png_header(const std::string &data) {
  if (data.size() < sizeof(kPngHeader) ||
      std::memcmp(data.data(), kPngHeader, sizeof(kPngHeader)) != 0) {
    throw RDKit::FileParseException("PNG header not recognized");
  }
}

// Walks the chunks as PNGStreamToMetadata (FileParsers/PNGParser.cpp) does.
// For a tEXt or zTXt chunk whose key and its NUL fill the chunk, RDKit takes
// &front() of the empty value string. Inputs RDKit rejects itself, with an
// exception, are left to it.
void check_png_metadata(const std::string &data) {
  check_png_header(data);
  std::size_t pos = sizeof(kPngHeader);
  while (data.size() - pos >= 8) {
    auto len = read_be32(data, pos);
    const char *tag = data.data() + pos + 4;
    std::size_t begin = pos + 8;
    if (std::memcmp(tag, "IEND", 4) == 0) {
      return;
    }
    bool text = len > 0 &&
                ((tag[0] == 't' && tag[1] == 'E') ||
                 (tag[0] == 'z' && tag[1] == 'T')) &&
                tag[2] == 'X' && tag[3] == 't';
    if (text) {
      auto nul = data.find('\0', begin);
      if (nul == std::string::npos) {
        // getline reads the key to the end of the input; the value's read
        // then fails, and RDKit throws, unless the value is empty.
        nul = data.size();
        if (nul == begin) {
          return;
        }
      }
      if (nul - begin + 1 == len) {
        throw RDKit::FileParseException("empty text chunk in PNG");
      }
      if (nul == data.size()) {
        return;
      }
    }
    // RDKit skips blockLen + 4 bytes, summed in uint32.
    pos = begin + static_cast<std::uint32_t>(len + 4);
    if (pos > data.size()) {
      return;
    }
  }
}

// addMetadataToPNGStream reads every chunk's length and tag, and the final
// CRC, without checking that they were read, so a PNG that ends before its
// IEND chunk is complete gives it uninitialized values.
void check_png_for_metadata(const std::string &data) {
  check_png_header(data);
  std::size_t pos = sizeof(kPngHeader);
  while (true) {
    if (data.size() - pos < 8) {
      break;
    }
    auto len = read_be32(data, pos);
    if (std::memcmp(data.data() + pos + 4, "IEND", 4) == 0) {
      if (data.size() - pos >= 12) {
        return;
      }
      break;
    }
    // RDKit reads blockLen + 4 bytes, summed in uint32.
    std::size_t next = pos + 8 + static_cast<std::uint32_t>(len + 4);
    if (next > data.size()) {
      break;
    }
    pos = next;
  }
  throw RDKit::FileParseException("did not find IEND block in PNG");
}

// As ReactionFromPNGStream in ChemReactions/PNGParser.cpp: the first record
// RDKit can use, with its pickle's templates checked (U10) and its rxn block
// checked before RDKit reads it (U12).
std::unique_ptr<ChemicalReaction> reaction_from_png(const std::string &data) {
  check_png_metadata(data);
  namespace RP = RDKit::v2::ReactionParser;
  for (const auto &[key, value] : RDKit::PNGStringToMetadata(data)) {
    std::unique_ptr<ChemicalReaction> rxn;
    if (boost::starts_with(key, RDKit::PNGData::rxnPklTag)) {
      rxn = std::make_unique<ChemicalReaction>(value);
      check_reaction_monomer_types(*rxn);
    } else if (boost::starts_with(key, RDKit::PNGData::rxnSmilesTag)) {
      check_cx_reaction(value);
      rxn = RP::ReactionFromSmiles(value);
    } else if (boost::starts_with(key, RDKit::PNGData::rxnSmartsTag)) {
      check_cx_reaction(value);
      rxn = RP::ReactionFromSmarts(value);
    } else if (boost::starts_with(key, RDKit::PNGData::rxnRxnTag)) {
      check_input_text(value);
      rxn = RP::ReactionFromRxnBlock(value);
    } else {
      continue;
    }
    return finish(std::move(rxn), "the PNG metadata");
  }
  throw RDKit::FileParseException("No suitable metadata found.");
}

rust::Vec<std::uint8_t> add_metadata(const ChemicalReaction &rxn,
                                     const std::string &png, bool include_pkl,
                                     bool include_smiles, bool include_smarts,
                                     bool include_rxn) {
  check_png_for_metadata(png);
  if (include_pkl) {
    check_picklable_templates(rxn);
  }
  if (include_smiles) {
    check_template_removals(rxn);
  }
  // The rxn block is the V2000 writer's.
  std::unique_lock<std::mutex> lock;
  if (include_rxn) {
    lock = lock_if_depicting(rxn, false);
  }
  return to_vec(RDKit::addChemicalReactionToPNGString(
      rxn, png, include_pkl, include_smiles, include_smarts, include_rxn));
}

RDKit::ReactionFingerprintParams fp_params_from(
    const ReactionFingerprintParamsData &data) {
  return RDKit::ReactionFingerprintParams(
      data.include_agents, data.bit_ratio_agents, data.non_agent_weight,
      data.agent_weight, data.fp_size, data.fp_type);
}

ReactionFingerprintParamsData fp_params_to_data(
    const RDKit::ReactionFingerprintParams &params) {
  return ReactionFingerprintParamsData{
      params.includeAgents, params.bitRatioAgents, params.nonAgentWeight,
      params.agentWeight,   params.fpSize,         params.fpType};
}

// A template's count fingerprint, as generateFingerprint in
// ReactionFingerprints.cpp makes it; null for the types it refuses.
std::unique_ptr<RDKit::SparseIntVect<std::uint32_t>> count_fingerprint(
    RDKit::ROMol &mol, unsigned fp_size, FingerprintType type) {
  mol.updatePropertyCache(false);
  switch (type) {
    case FingerprintType::AtomPairFP: {
      RDKit::AtomPair::AtomPairArguments args;
      args.d_fpSize = fp_size;
      std::unique_ptr<RDKit::FingerprintGenerator<std::uint32_t>> gen{
          RDKit::AtomPair::getAtomPairGenerator<std::uint32_t>(args)};
      return std::unique_ptr<RDKit::SparseIntVect<std::uint32_t>>(
          gen->getCountFingerprint(mol));
    }
    case FingerprintType::TopologicalTorsionFP: {
      RDKit::TopologicalTorsion::TopologicalTorsionArguments args;
      args.d_fpSize = fp_size;
      std::unique_ptr<RDKit::FingerprintGenerator<std::uint64_t>> gen{
          RDKit::TopologicalTorsion::getTopologicalTorsionGenerator<
              std::uint64_t>(args)};
      return std::unique_ptr<RDKit::SparseIntVect<std::uint32_t>>(
          gen->getCountFingerprint(mol));
    }
    case FingerprintType::MorganFP: {
      if (!mol.getRingInfo()->isInitialized()) {
        mol.updatePropertyCache();
        RDKit::MolOps::findSSSR(mol);
      }
      RDKit::MorganFingerprint::MorganArguments args;
      args.d_radius = 2;
      args.d_fpSize = fp_size;
      std::unique_ptr<RDKit::FingerprintGenerator<std::uint32_t>> gen{
          RDKit::MorganFingerprint::getMorganGenerator<std::uint32_t>(args)};
      return std::unique_ptr<RDKit::SparseIntVect<std::uint32_t>>(
          gen->getCountFingerprint(mol));
    }
    default:
      return nullptr;
  }
}

// Each bit's count summed over one kind of template, as RDKit sums them.
std::map<std::uint32_t, long long> count_sums(
    const ChemicalReaction &rxn, std::uint8_t kind,
    const RDKit::ReactionFingerprintParams &params) {
  std::map<std::uint32_t, long long> sums;
  for (const auto &mol : templates(rxn, kind)) {
    auto fp = count_fingerprint(*mol, params.fpSize, params.fpType);
    if (!fp) {
      continue;
    }
    for (const auto &[idx, val] : fp->getNonzeroElements()) {
      sums[idx] += val;
    }
  }
  return sums;
}

void check_int(long long value, const char *what) {
  if (value > INT_MAX || value < INT_MIN) {
    throw ValueErrorException(std::string("the difference fingerprint's ") +
                              what + " overflow RDKit's int counts");
  }
}

// DifferenceFingerprintChemReaction sums each kind of template's counts, then
// multiplies them by the weights and adds and subtracts them, all in int.
// This computes the counts it will, and refuses any value that overflows.
void check_difference_fingerprint(
    const ChemicalReaction &rxn,
    const RDKit::ReactionFingerprintParams &params) {
  auto type = static_cast<int>(params.fpType);
  if (params.fpSize == 0 || type <= 0 || type >= 4) {
    return;  // RDKit refuses these itself.
  }
  if (params.includeAgents && params.nonAgentWeight > INT_MAX) {
    // RDKit passes the weight to SparseIntVect::operator*=(int).
    throw ValueErrorException("nonAgentWeight is larger than INT_MAX");
  }
  auto reactants = count_sums(rxn, 0, params);
  auto products = count_sums(rxn, 1, params);
  std::map<std::uint32_t, long long> agents;
  if (params.includeAgents) {
    agents = count_sums(rxn, 2, params);
  }
  std::map<std::uint32_t, bool> bits;
  for (const auto *sums : {&reactants, &products, &agents}) {
    for (const auto &[idx, val] : *sums) {
      check_int(val, "summed counts");
      bits[idx] = true;
    }
  }
  for (const auto &[idx, unused] : bits) {
    long long r = reactants.count(idx) ? reactants[idx] : 0;
    long long p = products.count(idx) ? products[idx] : 0;
    if (!params.includeAgents) {
      check_int(p - r, "differences");
      continue;
    }
    long long a = agents.count(idx) ? agents[idx] : 0;
    long long w = params.nonAgentWeight;
    check_int(r * w, "weighted counts");
    check_int(p * w, "weighted counts");
    check_int(a * params.agentWeight, "weighted counts");
    check_int(p * w - r * w, "differences");
    check_int(p * w - r * w + a * params.agentWeight, "differences");
  }
}

// StructuralFingerprintChemReaction sizes its parts from fpSize and
// bitRatioAgents. It casts ceil(fpSize * bitRatioAgents) to int, which is
// undefined for a value outside int, and a bit vector over INT_MAX bits
// overflows RDKit's int arithmetic (U3).
void check_structural_fingerprint(
    const RDKit::ReactionFingerprintParams &params) {
  if (params.fpSize == 0) {
    return;  // RDKit refuses this itself.
  }
  unsigned final_size = params.fpSize / 2;
  if (params.includeAgents) {
    unsigned agent_size = params.fpSize / 3;
    if (params.bitRatioAgents < 1.0) {
      double scaled =
          std::ceil(static_cast<double>(params.fpSize) * params.bitRatioAgents);
      if (!(scaled >= INT_MIN && scaled <= INT_MAX)) {
        throw ValueErrorException(
            "fpSize * bitRatioAgents does not fit in an int");
      }
      agent_size = static_cast<unsigned>(static_cast<int>(scaled));
    }
    unsigned scaling = !(agent_size % 2) ? agent_size : agent_size - 1;
    final_size = (params.fpSize - scaling) / 2;
  }
  unsigned agent_size = params.fpSize - 2 * final_size;
  std::uint64_t total = 2 * std::uint64_t{final_size} + agent_size;
  if (total > INT_MAX) {
    throw ValueErrorException(
        "the structural fingerprint would have more than INT_MAX bits");
  }
}

// RDKit's fingerprints of the templates cache distance matrices on them
// (atom pairs), and compute them in int (U6).
std::vector<std::unique_ptr<MatrixCacheGuard>> guard_templates(
    const ChemicalReaction &rxn, bool check_size) {
  std::vector<std::unique_ptr<MatrixCacheGuard>> guards;
  for (std::uint8_t kind = 0; kind < 3; ++kind) {
    for (const auto &mol : templates(rxn, kind)) {
      if (check_size) {
        check_matrix_size(*mol);
      }
      guards.push_back(std::make_unique<MatrixCacheGuard>(*mol));
    }
  }
  return guards;
}

// fixAtomMaps (SanitizeRxn.cpp) gives a reactant R group without an atom map
// the map int(R label) + the largest atom map, summed in int.
void check_atom_map_room(const ChemicalReaction &rxn) {
  long long max_map = 0;
  for (std::uint8_t kind = 0; kind < 3; ++kind) {
    for (const auto &mol : templates(rxn, kind)) {
      for (const auto atom : mol->atoms()) {
        int map = 0;
        if (atom->getPropIfPresent<int>(
                RDKit::common_properties::molAtomMapNumber, map) &&
            map > max_map) {
          max_map = map;
        }
      }
    }
  }
  for (const auto &mol : rxn.getReactants()) {
    for (const auto atom : mol->atoms()) {
      unsigned int rlabel = 0;
      int map = 0;
      atom->getPropIfPresent(RDKit::common_properties::_MolFileRLabel, rlabel);
      atom->getPropIfPresent(RDKit::common_properties::molAtomMapNumber, map);
      if (rlabel && !map &&
          static_cast<int>(rlabel) + max_map >
              static_cast<long long>(INT_MAX)) {
        throw ValueErrorException(
            "an R label plus the largest atom map overflows an int");
      }
    }
  }
}

}  // namespace

void normalize_templates(ChemicalReaction &rxn) {
  normalize(rxn.beginReactantTemplates(), rxn.endReactantTemplates(),
            "reactant");
  normalize(rxn.beginProductTemplates(), rxn.endProductTemplates(), "product");
  normalize(rxn.beginAgentTemplates(), rxn.endAgentTemplates(), "agent");
}

std::unique_ptr<ProductSets> to_product_sets(
    const std::vector<MOL_SPTR_VECT> &products) {
  auto res = std::make_unique<ProductSets>();
  for (const auto &set : products) {
    std::vector<std::unique_ptr<RWMol>> mols;
    for (const auto &mol : set) {
      if (!mol) {
        throw std::logic_error("RDKit returned a null product");
      }
      auto copy = std::make_unique<RWMol>(*mol);
      clear_matrix_caches(*copy);
      mols.push_back(std::move(copy));
    }
    res->sets.push_back(std::move(mols));
  }
  return res;
}

std::size_t product_sets_size(const ProductSets &sets) {
  return sets.sets.size();
}

std::size_t product_sets_set_size(const ProductSets &sets, std::size_t i) {
  return sets.sets.at(i).size();
}

std::unique_ptr<RWMol> product_sets_take(ProductSets &sets, std::size_t i,
                                         std::size_t j) {
  return std::move(sets.sets.at(i).at(j));
}

std::size_t reaction_list_size(const ReactionList &list) {
  return list.items.size();
}

std::unique_ptr<ChemicalReaction> reaction_list_take(ReactionList &list,
                                                     std::size_t i) {
  return std::move(list.items.at(i));
}

std::unique_ptr<ChemicalReaction> chemical_reaction_new() {
  return std::make_unique<ChemicalReaction>();
}

std::unique_ptr<ChemicalReaction> chemical_reaction_clone(
    const ChemicalReaction &rxn) {
  // The copy constructor copies every template into a new RWMol.
  return std::make_unique<ChemicalReaction>(rxn);
}

void check_reaction_monomer_types(const ChemicalReaction &rxn) {
  for (std::uint8_t kind = 0; kind < 3; ++kind) {
    for (const auto &mol : templates(rxn, kind)) {
      if (mol) {
        check_pickled_mol(*mol);
      }
    }
  }
}

std::unique_ptr<ChemicalReaction> chemical_reaction_from_bin_str(
    rust::Slice<const std::uint8_t> bin_str) {
  // The templates are molecule pickles, which RDKit fuzzes (U2), checked as
  // those are (U10); the rest is counts and tags it checks, the
  // substructure parameters' JSON, and properties read as molecule pickles
  // read them.
  auto rxn = std::make_unique<ChemicalReaction>();
  RDKit::ReactionPickler::reactionFromPickle(to_std(bin_str), rxn.get());
  check_reaction_monomer_types(*rxn);
  normalize_templates(*rxn);
  return rxn;
}

rust::Vec<std::uint8_t> chemical_reaction_to_binary(const ChemicalReaction &rxn,
                                                    bool has_flags,
                                                    unsigned flags) {
  check_picklable_templates(rxn);
  std::string res;
  RDKit::ReactionPickler::pickleReaction(
      rxn, res,
      has_flags ? flags : RDKit::MolPickler::getDefaultPickleProperties());
  return to_vec(res);
}

const RDProps &chemical_reaction_as_props(const ChemicalReaction &rxn) {
  return rxn;
}

unsigned chemical_reaction_get_num_templates(const ChemicalReaction &rxn,
                                             std::uint8_t kind) {
  return templates(rxn, kind).size();
}

const RWMol &chemical_reaction_get_template(const ChemicalReaction &rxn,
                                            std::uint8_t kind, unsigned which) {
  const auto &mols = templates(rxn, kind);
  check_index(which, mols.size());
  // Every template is an RWMol (normalize_templates).
  const auto *mol = dynamic_cast<const RWMol *>(mols[which].get());
  if (!mol) {
    throw std::logic_error("a reaction template is not an RWMol");
  }
  return *mol;
}

unsigned chemical_reaction_add_template(ChemicalReaction &rxn,
                                        std::uint8_t kind, const RWMol &mol) {
  // Python's binding makes the reaction share Python's molecule; this one
  // adds a copy, so no reaction shares a molecule Rust owns.
  auto *copy = new RWMol(mol);
  clear_matrix_caches(*copy);
  ROMOL_SPTR tmpl(copy);
  switch (kind) {
    case 0:
      return rxn.addReactantTemplate(tmpl);
    case 1:
      return rxn.addProductTemplate(tmpl);
    case 2:
      return rxn.addAgentTemplate(tmpl);
    default:
      throw std::logic_error("unknown reaction template kind");
  }
}

std::unique_ptr<MolList> chemical_reaction_remove_unmapped_templates(
    ChemicalReaction &rxn, std::uint8_t kind, double threshold_unmapped_atoms,
    bool move_to_agent_templates) {
  MOL_SPTR_VECT removed;
  if (kind == 0) {
    rxn.removeUnmappedReactantTemplates(threshold_unmapped_atoms,
                                        move_to_agent_templates, &removed);
  } else {
    rxn.removeUnmappedProductTemplates(threshold_unmapped_atoms,
                                       move_to_agent_templates, &removed);
  }
  // Copies, since the agent templates may hold the same molecules.
  return copies(removed);
}

std::unique_ptr<MolList> chemical_reaction_remove_agent_templates(
    ChemicalReaction &rxn) {
  MOL_SPTR_VECT removed;
  rxn.removeAgentTemplates(&removed);
  return copies(removed);
}

void chemical_reaction_initialize(ChemicalReaction &rxn, bool silent) {
  rxn.initReactantMatchers(silent);
}

bool chemical_reaction_is_initialized(const ChemicalReaction &rxn) {
  return rxn.isInitialized();
}

rust::Vec<std::uint32_t> chemical_reaction_validate(ChemicalReaction &rxn,
                                                    bool silent) {
  unsigned int num_warnings = 0;
  unsigned int num_errors = 0;
  rxn.validate(num_warnings, num_errors, silent);
  rust::Vec<std::uint32_t> res;
  res.push_back(num_warnings);
  res.push_back(num_errors);
  return res;
}

// Running a reaction, testing a molecule against its templates, comparing
// reactions, preprocessing, and the fingerprints (which can be pattern
// fingerprints) match substructures (control_c_lock.h).
std::unique_ptr<ProductSets> chemical_reaction_run_reactants(
    ChemicalReaction &rxn, BorrowedMols &reactants, unsigned max_products) {
  initialize_if_needed(rxn);
  // run_Reactants clears every reactant's atom bookmarks (ReactionRunner.cpp),
  // though Python passes them as const, so it runs on copies.
  MOL_SPTR_VECT copies;
  for (const auto &mol : reactants.mols()) {
    copies.emplace_back(new RWMol(*mol));
  }
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_product_sets(rxn.runReactants(copies, max_products));
}

std::unique_ptr<ProductSets> chemical_reaction_run_reactant(
    ChemicalReaction &rxn, const RWMol &reactant, unsigned reaction_idx) {
  initialize_if_needed(rxn);
  // RDKit checks reactionIdx with a PRECONDITION (Python's RuntimeError),
  // after it has cleared the reactant's bookmarks.
  check_index(reaction_idx, rxn.getNumReactantTemplates());
  ROMOL_SPTR copy(new RWMol(reactant));
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_product_sets(rxn.runReactant(copy, reaction_idx));
}

bool chemical_reaction_run_reactant_in_place(ChemicalReaction &rxn,
                                             RWMol &reactant,
                                             bool remove_unmatched_atoms) {
  initialize_if_needed(rxn);
  // Python's binding static_casts its ROMol to RWMol, undefined for a Mol;
  // a Rust molecule is an RWMol (W1). Its edits keep cached matrices (U6),
  // and its removals trust what check_removals checks, checked before the
  // first edit.
  check_removals(reactant);
  MatrixCacheGuard guard(reactant);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  try {
    return rxn.runReactant(reactant, remove_unmatched_atoms);
  } catch (...) {
    // run_Reactant edits in a batch, which an exception leaves open.
    reactant.rollbackBatchEdit();
    throw;
  }
}

bool chemical_reaction_is_molecule_template(const ChemicalReaction &rxn,
                                            std::uint8_t kind,
                                            const RWMol &mol) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  switch (kind) {
    case 0:
      return RDKit::isMoleculeReactantOfReaction(rxn, mol);
    case 1:
      return RDKit::isMoleculeProductOfReaction(rxn, mol);
    case 2:
      return RDKit::isMoleculeAgentOfReaction(rxn, mol);
    default:
      throw std::logic_error("unknown reaction template kind");
  }
}

rust::Vec<IndexList> chemical_reaction_get_reacting_atoms(
    const ChemicalReaction &rxn, bool mapped_atoms_only) {
  return to_index_lists(RDKit::getReactingAtoms(rxn, mapped_atoms_only));
}

rust::Vec<ReactantLabel> chemical_reaction_add_recursive_queries(
    ChemicalReaction &rxn, BorrowedMols &queries, rust::Str prop_name,
    bool get_labels) {
  auto qs = query_map(queries);
  std::string prop(prop_name);
  if (!get_labels) {
    RDKit::addRecursiveQueriesToReaction(rxn, qs, prop);
    return {};
  }
  std::vector<std::vector<std::pair<unsigned int, std::string>>> labels;
  RDKit::addRecursiveQueriesToReaction(rxn, qs, prop, &labels);
  return to_labels(labels);
}

SubstructMatchParametersData chemical_reaction_get_substruct_params(
    const ChemicalReaction &rxn) {
  return substruct_params_to_data(rxn.getSubstructParams());
}

void chemical_reaction_set_substruct_params(
    ChemicalReaction &rxn, const SubstructMatchParametersData &params) {
  rxn.getSubstructParams() = substruct_params_from_data(params);
}

std::unique_ptr<ChemicalReaction> reaction_from_smarts(
    rust::Str smarts, const rust::Vec<rust::String> &replacement_keys,
    const rust::Vec<rust::String> &replacement_values, bool use_smiles) {
  std::map<std::string, std::string> replacements;
  for (std::size_t i = 0;
       i < replacement_keys.size() && i < replacement_values.size(); ++i) {
    replacements[std::string(replacement_keys[i])] =
        std::string(replacement_values[i]);
  }
  std::string text(smarts);
  check_cx_reaction(text);
  return finish(
      std::unique_ptr<ChemicalReaction>(
          RDKit::RxnSmartsToChemicalReaction(text, &replacements, use_smiles)),
      text);
}

namespace {

// The reaction writers write each template with RDKit's SMILES or SMARTS
// writer, whose doRandom advances RDKit's global generator (U4).
std::unique_lock<std::mutex> random_lock(const SmilesWriteParamsData &params) {
  return params.do_random ? std::unique_lock<std::mutex>(random_mutex())
                          : std::unique_lock<std::mutex>();
}

}  // namespace

std::unique_ptr<std::string> reaction_to_smarts(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params) {
  auto lock = random_lock(params);
  return own(RDKit::ChemicalReactionToRxnSmarts(rxn, write_params(params)));
}

std::unique_ptr<std::string> reaction_to_smiles(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params) {
  check_template_removals(rxn);
  auto lock = random_lock(params);
  return own(RDKit::ChemicalReactionToRxnSmiles(rxn, write_params(params)));
}

std::unique_ptr<std::string> reaction_to_cx_smarts(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params,
    std::uint32_t flags) {
  auto lock = random_lock(params);
  return own(
      RDKit::ChemicalReactionToCXRxnSmarts(rxn, write_params(params), flags));
}

std::unique_ptr<std::string> reaction_to_cx_smiles(
    const ChemicalReaction &rxn, const SmilesWriteParamsData &params,
    std::uint32_t flags) {
  check_template_removals(rxn);
  auto lock = random_lock(params);
  return own(
      RDKit::ChemicalReactionToCXRxnSmiles(rxn, write_params(params), flags));
}

std::uint32_t cx_smiles_fields_all() {
  return RDKit::SmilesWrite::CXSmilesFields::CX_ALL;
}

std::unique_ptr<ChemicalReaction> reaction_from_rxn_block(rust::Str rxn_block,
                                                          bool sanitize,
                                                          bool remove_hs,
                                                          bool strict_parsing) {
  std::string text(rxn_block);
  // Its templates are mol blocks, read by RDKit's line-oriented reader (U12).
  check_input_text(text);
  RDKit::v2::FileParsers::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.strictParsing = strict_parsing;
  return finish(RDKit::v2::ReactionParser::ReactionFromRxnBlock(text, params),
                "the rxn block");
}

std::unique_ptr<ChemicalReaction> reaction_from_rxn_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs,
    bool strict_parsing) {
  // Read once and checked, as reaction_from_rxn_block checks a block (U2,
  // U12); RDKit's ReactionFromRxnFile parses an opened file as it parses a
  // block.
  auto text = read_path(path);
  check_input_text(text);
  RDKit::v2::FileParsers::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.strictParsing = strict_parsing;
  return finish(RDKit::v2::ReactionParser::ReactionFromRxnBlock(text, params),
                "the rxn file " + to_std(path));
}

std::unique_ptr<ChemicalReaction> reaction_from_mrv_block(rust::Str mrv,
                                                          bool sanitize,
                                                          bool remove_hs) {
  std::string text(mrv);
  check_mrv(text);
  RDKit::v2::MarvinParser::MrvParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  return finish(RDKit::v2::MarvinParser::ReactionFromMrvBlock(text, params),
                "the Marvin block");
}

std::unique_ptr<ChemicalReaction> reaction_from_mrv_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs) {
  // Read here so the parser reads exactly what was checked.
  auto text = read_path(path);
  check_mrv(text);
  RDKit::v2::MarvinParser::MrvParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  return finish(RDKit::v2::MarvinParser::ReactionFromMrvBlock(text, params),
                "the Marvin file " + to_std(path));
}

bool mrv_block_is_reaction(rust::Str mrv) {
  return RDKit::v2::MarvinParser::MrvBlockIsReaction(std::string(mrv));
}

bool mrv_file_is_reaction(rust::Slice<const std::uint8_t> path) {
  return RDKit::v2::MarvinParser::MrvFileIsReaction(to_std(path));
}

namespace {

bool is_xml_space(char c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

// The ChemDraw library reads a reaction step's atom maps
// (ReactionStepAtomMap, ...Manual, ...Auto) as ints until the end of the
// value, and never reaches it once a token fails to read as one: it appends
// to the map until memory runs out. Every token must be an int, then, and
// an entity, which would decode to text this cannot see, is refused too.
void check_cdxml_atom_maps(const std::string &text) {
  static const std::string kName = "ReactionStepAtomMap";
  for (auto at = text.find(kName); at != std::string::npos;
       at = text.find(kName, at + 1)) {
    if (at == 0 || !is_xml_space(text[at - 1])) {
      continue;
    }
    auto pos = at + kName.size();
    for (const char *suffix : {"Manual", "Auto"}) {
      if (text.compare(pos, std::strlen(suffix), suffix) == 0) {
        pos += std::strlen(suffix);
        break;
      }
    }
    while (pos < text.size() && is_xml_space(text[pos])) {
      ++pos;
    }
    if (pos >= text.size() || text[pos] != '=') {
      continue;
    }
    ++pos;
    while (pos < text.size() && is_xml_space(text[pos])) {
      ++pos;
    }
    if (pos >= text.size() || (text[pos] != '"' && text[pos] != '\'')) {
      continue;
    }
    auto end = text.find(text[pos], pos + 1);
    if (end == std::string::npos) {
      continue;
    }
    std::istringstream tokens(text.substr(pos + 1, end - pos - 1));
    std::string token;
    while (tokens >> token) {
      auto digits =
          token.find_first_not_of("+-") == 1 ? token.substr(1) : token;
      bool ok = token.find('&') == std::string::npos && !digits.empty() &&
                digits.find_first_not_of("0123456789") == std::string::npos;
      if (ok) {
        try {
          std::size_t used = 0;
          long long value = std::stoll(token, &used);
          ok = used == token.size() && value >= INT_MIN && value <= INT_MAX;
        } catch (const std::exception &) {
          ok = false;
        }
      }
      if (!ok) {
        throw RDKit::FileParseException("a reaction step's atom map holds \"" +
                                        token + "\", which is not an int");
      }
    }
  }
}

std::unique_ptr<ReactionList> to_reaction_list(
    std::vector<std::unique_ptr<ChemicalReaction>> rxns) {
  auto res = std::make_unique<ReactionList>();
  for (auto &rxn : rxns) {
    // The CDXML parser shares templates between the reactions it returns;
    // normalize_templates gives each its own.
    res->items.push_back(finish(std::move(rxn), "the CDXML"));
  }
  return res;
}

}  // namespace

std::unique_ptr<ReactionList> reactions_from_cdxml_block(rust::Str cdxml,
                                                         bool sanitize,
                                                         bool remove_hs) {
  std::string text(cdxml);
  check_cdxml_atom_maps(text);
  return to_reaction_list(
      RDKit::CDXMLToChemicalReactions(text, sanitize, remove_hs));
}

std::unique_ptr<ReactionList> reactions_from_cdxml_file(
    rust::Slice<const std::uint8_t> path, bool sanitize, bool remove_hs) {
  // RDKit returns no reactions for a file it cannot open.
  auto text = read_path(path);
  check_cdxml_atom_maps(text);
  return to_reaction_list(
      RDKit::CDXMLToChemicalReactions(text, sanitize, remove_hs));
}

std::unique_ptr<std::string> reaction_to_rxn_block(const ChemicalReaction &rxn,
                                                   bool separate_agents,
                                                   bool force_v3000) {
  auto lock = lock_if_depicting(rxn, force_v3000);
  return own(
      RDKit::ChemicalReactionToRxnBlock(rxn, separate_agents, force_v3000));
}

std::unique_ptr<std::string> reaction_to_v3k_rxn_block(
    const ChemicalReaction &rxn, bool separate_agents) {
  return own(RDKit::ChemicalReactionToV3KRxnBlock(rxn, separate_agents));
}

// The Marvin reaction writer copies each template and writes the
// coordinates it has, without computing any (MarvinWriter.cpp).
std::unique_ptr<std::string> reaction_to_mrv_block(const ChemicalReaction &rxn,
                                                   bool pretty_print) {
  return own(RDKit::ChemicalReactionToMrvBlock(rxn, pretty_print));
}

void reaction_to_mrv_file(const ChemicalReaction &rxn,
                          rust::Slice<const std::uint8_t> path,
                          bool pretty_print) {
  RDKit::ChemicalReactionToMrvFile(rxn, to_std(path), pretty_print);
}

// RxnMolToChemicalReaction splits the molecule into its fragments
// (getMolFrags), removing atoms from copies of it.
std::unique_ptr<ChemicalReaction> reaction_from_molecule(const RWMol &mol) {
  check_removals(mol);
  return finish(
      std::unique_ptr<ChemicalReaction>(RDKit::RxnMolToChemicalReaction(mol)),
      "the molecule");
}

std::unique_ptr<RWMol> reaction_to_molecule(const ChemicalReaction &rxn) {
  return to_rwmol(
      std::unique_ptr<RDKit::ROMol>(RDKit::ChemicalReactionToRxnMol(rxn)));
}

std::unique_ptr<ChemicalReaction> reaction_from_png_string(
    rust::Slice<const std::uint8_t> data) {
  return reaction_from_png(to_std(data));
}

std::unique_ptr<ChemicalReaction> reaction_from_png_file(
    rust::Slice<const std::uint8_t> path) {
  return reaction_from_png(read_path(path));
}

rust::Vec<std::uint8_t> reaction_metadata_to_png_string(
    const ChemicalReaction &rxn, rust::Slice<const std::uint8_t> png,
    bool include_pkl, bool include_smiles, bool include_smarts,
    bool include_rxn) {
  return add_metadata(rxn, to_std(png), include_pkl, include_smiles,
                      include_smarts, include_rxn);
}

rust::Vec<std::uint8_t> reaction_metadata_to_png_file(
    const ChemicalReaction &rxn, rust::Slice<const std::uint8_t> path,
    bool include_pkl, bool include_smiles, bool include_smarts,
    bool include_rxn) {
  // RDKit opens the file without checking that it opened.
  return add_metadata(rxn, read_path(path), include_pkl, include_smiles,
                      include_smarts, include_rxn);
}

void compute_2d_coords_for_reaction(ChemicalReaction &rxn, double spacing,
                                    bool update_props, bool canon_orient,
                                    unsigned n_flips_per_sample,
                                    unsigned n_sample, int sample_seed,
                                    bool permute_deg4_nodes,
                                    double bond_length) {
  check_bond_length(bond_length);
  auto lock = lock_depiction();
  // Random bond flips (nSample and nFlipsPerSample) draw from RDKit's
  // generator (U4).
  std::lock_guard<std::mutex> random_lock(random_mutex());
  // compute2DCoordsForReaction runs compute2DCoords on each reactant and
  // product template (ReactionDepict.cpp), which trusts the molecule
  // (depiction_checks.h); its updateProps changes nothing the checks read.
  // RDDepict computes the templates' distance matrices, in int, and caches
  // them (U6). Agents are not depicted.
  std::vector<std::unique_ptr<MatrixCacheGuard>> guards;
  for (std::uint8_t kind = 0; kind < 2; ++kind) {
    for (const auto &mol : templates(rxn, kind)) {
      check_matrix_size(*mol);
      check_compute_2d_coords_mol(*mol, false);
      guards.push_back(std::make_unique<MatrixCacheGuard>(*mol));
    }
  }
  BondLengthOverride setter(bond_length);
  RDDepict::compute2DCoordsForReaction(rxn, spacing, update_props, canon_orient,
                                       n_flips_per_sample, n_sample,
                                       sample_seed, permute_deg4_nodes);
}

ReactionFingerprintParamsData reaction_fingerprint_params_default() {
  return fp_params_to_data(RDKit::ReactionFingerprintParams());
}

ReactionFingerprintParamsData default_difference_fp_params() {
  return fp_params_to_data(RDKit::DefaultDifferenceFPParams);
}

ReactionFingerprintParamsData default_structural_fp_params() {
  return fp_params_to_data(RDKit::DefaultStructuralFPParams);
}

std::unique_ptr<SparseIntVect> create_difference_fingerprint_for_reaction(
    const ChemicalReaction &rxn, const ReactionFingerprintParamsData &data) {
  auto params = fp_params_from(data);
  auto guards =
      guard_templates(rxn, params.fpType == FingerprintType::AtomPairFP);
  check_difference_fingerprint(rxn, params);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::unique_ptr<RDKit::SparseIntVect<std::uint32_t>> fp(
      RDKit::DifferenceFingerprintChemReaction(rxn, params));
  return std::make_unique<SparseIntVect>(std::move(*fp));
}

std::unique_ptr<BitVect> create_structural_fingerprint_for_reaction(
    const ChemicalReaction &rxn, const ReactionFingerprintParamsData &data) {
  auto params = fp_params_from(data);
  check_structural_fingerprint(params);
  auto guards =
      guard_templates(rxn, params.fpType == FingerprintType::AtomPairFP);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::unique_ptr<ExplicitBitVect> fp(
      RDKit::StructuralFingerprintChemReaction(rxn, params));
  return std::make_unique<BitVect>(std::move(*fp));
}

bool is_reaction_template_molecule_agent(const RWMol &mol,
                                         double agent_threshold) {
  return RDKit::isReactionTemplateMoleculeAgent(mol, agent_threshold);
}

bool has_reaction_atom_mapping(const ChemicalReaction &rxn) {
  return RDKit::hasReactionAtomMapping(rxn);
}

bool has_reaction_substruct_match(const ChemicalReaction &rxn,
                                  const ChemicalReaction &query_rxn,
                                  bool include_agents) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return RDKit::hasReactionSubstructMatch(rxn, query_rxn, include_agents);
}

bool has_template_substruct_match(const ChemicalReaction &rxn,
                                  const ChemicalReaction &query_rxn,
                                  std::uint8_t kind) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  switch (kind) {
    case 0:
      return RDKit::hasReactantTemplateSubstructMatch(rxn, query_rxn);
    case 1:
      return RDKit::hasProductTemplateSubstructMatch(rxn, query_rxn);
    case 2:
      return RDKit::hasAgentTemplateSubstructMatch(rxn, query_rxn);
    default:
      throw std::logic_error("unknown reaction template kind");
  }
}

void update_products_stereochemistry(ChemicalReaction &rxn) {
  RDKit::updateProductsStereochem(&rxn);
}

std::unique_ptr<RWMol> reduce_product_to_side_chains(const RWMol &product,
                                                     bool add_dummy_atoms) {
  // RDKit copies the product before it edits, and removes atoms from the
  // copy (check_removals).
  check_removals(product);
  ROMOL_SPTR borrowed(const_cast<RWMol *>(&product), [](RDKit::ROMol *) {});
  return to_rwmol(std::unique_ptr<RDKit::ROMol>(
      RDKit::reduceProductToSideChains(borrowed, add_dummy_atoms)));
}

void remove_mapping_numbers_from_reactions(ChemicalReaction &rxn) {
  RDKit::removeMappingNumbersFromReactions(rxn);
}

PreprocessReactionData preprocess_reaction(ChemicalReaction &rxn,
                                           BorrowedMols &queries,
                                           rust::Str prop_name) {
  // As PreprocessReaction in Wrap/rdChemReactions.cpp.
  std::map<std::string, ROMOL_SPTR> qs;
  if (queries.mols().empty()) {
    const bool normalized = true;
    qs = RDKit::GetFlattenedFunctionalGroupHierarchy(normalized);
  } else {
    qs = query_map(queries);
  }
  unsigned int num_reactants = rxn.getNumReactantTemplates();
  unsigned int num_products = rxn.getNumProductTemplates();
  unsigned int num_warnings = 0;
  unsigned int num_errors = 0;
  rxn.validate(num_warnings, num_errors);
  std::vector<std::vector<std::pair<unsigned int, std::string>>> labels;
  if (!num_errors) {
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    RDKit::preprocessReaction(rxn, num_warnings, num_errors, labels, qs,
                              std::string(prop_name));
  }
  return PreprocessReactionData{num_warnings,
                                num_errors,
                                num_reactants,
                                num_products,
                                static_cast<std::uint32_t>(labels.size()),
                                to_labels(labels)};
}

namespace {

// As RxnOps::adjustTemplates, which runs adjustQueryProperties on each
// template: that replaces bonds with RWMol::replaceBond, which would leave a
// template's stereo groups pointing at the bonds it deletes (replace_bond),
// so the groups are guarded.
void adjust_templates(const RDKit::MOL_SPTR_VECT &templates,
                      const RDKit::MolOps::AdjustQueryParameters &params) {
  std::vector<std::unique_ptr<StereoGroupBondGuard>> guards;
  for (const auto &templ : templates) {
    guards.push_back(std::make_unique<StereoGroupBondGuard>(*templ));
  }
  RDKit::RxnOps::adjustTemplates(templates, params);
}

}  // namespace

std::uint32_t sanitize_rxn(ChemicalReaction &rxn, std::uint32_t sanitize_ops,
                           const AdjustQueryParametersData &data,
                           bool catch_errors) {
  namespace RxnOps = RDKit::RxnOps;
  auto params = adjust_query_params_from_data(data);
  if (sanitize_ops & RxnOps::SANITIZE_MERGEHS) {
    // Checked before the first step edits the reaction.
    check_template_removals(rxn);
  }
  unsigned int failed = RxnOps::SANITIZE_NONE;
  try {
    // RxnOps::sanitizeRxn's steps, with fixAtomMaps checked first.
    if (sanitize_ops & RxnOps::SANITIZE_RGROUP_NAMES) {
      failed = RxnOps::SANITIZE_RGROUP_NAMES;
      RxnOps::fixRGroups(rxn);
    }
    if (sanitize_ops & RxnOps::SANITIZE_ATOM_MAPS) {
      failed = RxnOps::SANITIZE_ATOM_MAPS;
      check_atom_map_room(rxn);
      RxnOps::fixAtomMaps(rxn);
    }
    if (sanitize_ops & RxnOps::SANITIZE_ADJUST_REACTANTS) {
      failed = RxnOps::SANITIZE_ADJUST_REACTANTS;
      adjust_templates(rxn.getReactants(), params);
    }
    if (sanitize_ops & RxnOps::SANITIZE_ADJUST_PRODUCTS) {
      failed = RxnOps::SANITIZE_ADJUST_PRODUCTS;
      adjust_templates(rxn.getProducts(), params);
    }
    if (sanitize_ops & RxnOps::SANITIZE_MERGEHS) {
      failed = RxnOps::SANITIZE_MERGEHS;
      RxnOps::fixHs(rxn);
    }
    failed = RxnOps::SANITIZE_NONE;
  } catch (...) {
    // As sanitizeReaction in Wrap/rdChemReactions.cpp.
    if (!catch_errors) {
      throw;
    }
  }
  return failed;
}

void sanitize_rxn_as_mols(ChemicalReaction &rxn, std::uint32_t sanitize_ops) {
  // RDKit dereferences dynamic_cast<RWMol *> of each template, which
  // normalize_templates makes valid.
  RDKit::RxnOps::sanitizeRxnAsMols(rxn, sanitize_ops);
}

AdjustQueryParametersData get_default_adjust_params() {
  return adjust_query_params_to_data(RDKit::RxnOps::DefaultRxnAdjustParams());
}

AdjustQueryParametersData get_chem_draw_rxn_adjust_params() {
  return adjust_query_params_to_data(RDKit::RxnOps::ChemDrawRxnAdjustParams());
}

AdjustQueryParametersData match_only_at_rgroups_adjust_params() {
  return adjust_query_params_to_data(
      RDKit::RxnOps::MatchOnlyAtRgroupsAdjustParams());
}

}  // namespace rdkit_cxx
