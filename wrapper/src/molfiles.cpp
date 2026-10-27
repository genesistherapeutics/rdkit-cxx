// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/molfiles.h"

#include <climits>
#include <cmath>
#include <map>
#include <mutex>
#include <optional>
#include <sstream>

#include <GraphMol/Canon.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/MarvinParse/MarvinParser.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/new_canon.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/StreamOps.h>
#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/molfiles.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/streams.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace FP = RDKit::v2::FileParsers;

std::map<std::string, std::string> replacements(
    const rust::Vec<rust::String> &keys,
    const rust::Vec<rust::String> &values) {
  std::map<std::string, std::string> result;
  for (std::size_t i = 0; i < keys.size() && i < values.size(); ++i) {
    result[std::string(keys[i])] = std::string(values[i]);
  }
  return result;
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

std::unique_ptr<RWMol> parsed(std::unique_ptr<RWMol> mol,
                              const std::string &what) {
  if (!mol) {
    throw RDKit::FileParseException("could not parse the " + what);
  }
  return mol;
}

std::unique_ptr<std::string> str(std::string s) {
  return std::make_unique<std::string>(std::move(s));
}

std::vector<int> atom_indices(const RWMol &mol,
                              rust::Slice<const std::uint32_t> atoms) {
  std::vector<int> result;
  result.reserve(atoms.size());
  for (auto idx : atoms) {
    check_index(idx, mol.getNumAtoms());
    result.push_back(static_cast<int>(idx));
  }
  return result;
}

std::vector<int> bond_indices(const RWMol &mol,
                              rust::Slice<const std::uint32_t> bonds) {
  std::vector<int> result;
  result.reserve(bonds.size());
  for (auto idx : bonds) {
    check_index(idx, mol.getNumBonds());
    result.push_back(static_cast<int>(idx));
  }
  return result;
}

std::vector<std::string> strings(const rust::Vec<rust::String> &values) {
  std::vector<std::string> result;
  result.reserve(values.size());
  for (const auto &value : values) {
    result.emplace_back(value);
  }
  return result;
}

// The randomized writers draw the root atom (SmilesWrite.cpp) and the
// traversal order (Canon.cpp) from RDKit's global random number generator,
// which has no lock, so two threads writing random SMILES race on it (U4).
std::unique_lock<std::mutex> random_lock(bool random) {
  return random ? std::unique_lock<std::mutex>(random_mutex())
                : std::unique_lock<std::mutex>();
}

// A copy of RDKit's walk through a TPL block in MolFromTPLDataStream
// (TplFileParser.cpp), which indexes the second token of the counts line and
// of the CONFS line without checking that the line has one: a counts line or
// a CONFS line with a single token reads past the end of the token vector.
void check_tpl(const std::string &text) {
  std::istringstream in(text);
  auto line = [&in]() { return RDKit::getLine(in); };
  auto tokens = [](const std::string &s) {
    std::vector<std::string> result;
    boost::split(result, s, boost::is_any_of(" \t"), boost::token_compress_on);
    return result;
  };
  auto too_few = []() {
    throw RDKit::FileParseException(
        "a TPL counts or CONFS line needs two fields");
  };
  // The format, comment, and name lines; RDKit returns null at EOF.
  std::string tempStr;
  for (int i = 0; i < 3; ++i) {
    tempStr = line();
    if (in.eof()) {
      return;
    }
  }
  if (tempStr.size() >= 4 && tempStr.substr(0, 4) == "NAME") {
    tempStr = line();
    if (in.eof()) {
      return;
    }
  }
  if (tempStr.size() >= 4 && tempStr.substr(0, 4) == "PROP") {
    tempStr = line();
    if (in.eof()) {
      return;
    }
  }
  auto counts = tokens(tempStr);
  if (counts.size() < 2) {
    too_few();
  }
  unsigned int nAtoms, nBonds;
  try {
    nAtoms =
        RDKit::FileParserUtils::stripSpacesAndCast<unsigned int>(counts[0]);
    nBonds =
        RDKit::FileParserUtils::stripSpacesAndCast<unsigned int>(counts[1]);
  } catch (const boost::bad_lexical_cast &) {
    return;  // RDKit throws the same cast error
  }
  // RDKit reads one line per atom and bond, stopping with an error at EOF.
  // It sizes a conformer for nAtoms first, so missing atom lines are refused
  // here, with RDKit's message, before it allocates for a count the block
  // cannot hold.
  for (std::uint64_t i = 0; i < std::uint64_t{nAtoms} + nBonds; ++i) {
    line();
    if (in.eof()) {
      if (i < nAtoms) {
        throw RDKit::FileParseException("EOF hit while reading atoms.");
      }
      return;
    }
  }
  tempStr = line();
  if (in.eof()) {
    return;
  }
  if (tempStr.size() >= 5 && tempStr.substr(0, 5) == "CONFS" &&
      tokens(tempStr).size() < 2) {
    too_few();
  }
}

// RDKit's SVG reader (ptreeToMol in SVGParser.cpp) subtracts 1 from each
// bond's atom indices in int, which overflows for INT_MIN.
void check_svg(const std::string &svg) {
  namespace pt = boost::property_tree;
  std::stringstream in(svg);
  pt::ptree tree;
  pt::read_xml(in, tree);
  pt::ptree empty;
  for (const auto &[key, mol] : tree.get_child("svg.metadata", empty)) {
    if (key != "rdkit:mol") {
      continue;
    }
    for (const auto &[bkey, bond] : mol) {
      if (bkey != "rdkit:bond") {
        continue;
      }
      for (const char *attr :
           {"<xmlattr>.begin-atom-idx", "<xmlattr>.end-atom-idx"}) {
        auto idx = bond.get_optional<int>(attr);
        if (idx && *idx == INT_MIN) {
          throw ValueErrorException("a bond atom index of INT_MIN");
        }
      }
    }
  }
}

// RDKit's proximity bonding (ConnectTheDots_Large in ProximityBonds.cpp)
// casts each coordinate divided by 5.45 to int and sums the casts with
// multipliers up to 571 in int, which overflows for a coordinate beyond about
// 1.7e7; PDB fields reach 1e8. RDKit bonds the first conformer of the
// molecule it read, which this reads first without proximity bonding.
constexpr double kMaxProximityCoordinate = 1e7;

void check_pdb_coordinates(const RDKit::ROMol &probe) {
  if (!probe.getNumConformers()) {
    return;
  }
  for (const auto &pos : probe.getConformer().getPositions()) {
    for (double c : {pos.x, pos.y, pos.z}) {
      if (!(std::fabs(c) <= kMaxProximityCoordinate)) {
        throw ValueErrorException(
            "PDB coordinates must be finite and at most 1e7 in magnitude for "
            "proximity bonding");
      }
    }
  }
}

FP::PDBParserParams pdb_params(bool sanitize, bool remove_hs,
                               std::uint32_t flavor, bool proximity_bonding) {
  FP::PDBParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.flavor = flavor;
  params.proximityBonding = proximity_bonding;
  return params;
}

// RDKit's property lists (getPropertyList in FileParserUtils.h) drop the
// last character of the line so far when a value fills the line; with no
// missing-value marker the first line is empty, and dropping from it is
// undefined. That happens only for the first value.
template <typename T, typename Range>
void check_first_list_value(Range items, const std::string &prop,
                            const std::string &marker, unsigned line_size) {
  if (!marker.empty()) {
    return;
  }
  for (const auto item : items) {
    std::string value = "n/a";
    T tVal;
    if (item->getPropIfPresent(prop, tVal)) {
      value = boost::lexical_cast<std::string>(tVal);
    }
    if (value.length() + 1 >= line_size) {
      throw ValueErrorException(
          "lineSize must be longer than the first value when there is no "
          "missing-value marker");
    }
    return;
  }
}

}  // namespace

RDKit::SmilesWriteParams write_params(const SmilesWriteParamsData &data) {
  RDKit::SmilesWriteParams params;
  params.doIsomericSmiles = data.do_isomeric_smiles;
  params.doKekule = data.do_kekule;
  params.canonical = data.canonical;
  params.allBondsExplicit = data.all_bonds_explicit;
  params.allHsExplicit = data.all_hs_explicit;
  params.doRandom = data.do_random;
  params.rootedAtAtom = data.rooted_at_atom;
  params.includeDativeBonds = data.include_dative_bonds;
  params.ignoreAtomMapNumbers = data.ignore_atom_map_numbers;
  return params;
}

RDKit::MolWriterParams writer_params(const MolWriterParamsData &data) {
  check_precision(data.precision);
  RDKit::MolWriterParams params;
  params.includeStereo = data.include_stereo;
  params.kekulize = data.kekulize;
  params.forceV3000 = data.force_v3000;
  params.precision = data.precision;
  return params;
}

RDKit::SmilesParserParams smiles_parser_params(
    const SmilesParserParamsData &data,
    std::map<std::string, std::string> &repl) {
  check_debug_parse(data.debug_parse);
  repl = replacements(data.replacement_keys, data.replacement_values);
  RDKit::SmilesParserParams params;
  params.debugParse = data.debug_parse;
  params.parseName = data.parse_name;
  params.allowCXSMILES = data.allow_cxsmiles;
  params.strictCXSMILES = data.strict_cxsmiles;
  params.sanitize = data.sanitize;
  params.removeHs = data.remove_hs;
  params.replacements = &repl;
  return params;
}

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
      params.doIsomericSmiles,     params.doKekule,
      params.canonical,            params.allBondsExplicit,
      params.allHsExplicit,        params.doRandom,
      params.rootedAtAtom,         params.includeDativeBonds,
      params.ignoreAtomMapNumbers};
}

MolWriterParamsData mol_writer_params_default() {
  RDKit::MolWriterParams params;
  return MolWriterParamsData{params.includeStereo, params.kekulize,
                             params.forceV3000, params.precision, -1};
}

// As MolFromSmilesHelper in Code/GraphMol/Wrap/rdmolfiles.cpp.
std::unique_ptr<RWMol> mol_from_smiles(rust::Str smiles,
                                       const SmilesParserParamsData &data) {
  std::string text(smiles);
  std::map<std::string, std::string> repl;
  auto params = smiles_parser_params(data, repl);
  return parsed(RDKit::SmilesToMol(text, params), "SMILES", text);
}

// MolToSmiles writes a molecule of several fragments a fragment at a time,
// making some by removing the other atoms from a copy (getMolFrags), which
// trusts what check_removals checks.
std::unique_ptr<std::string> mol_to_smiles(const RWMol &mol,
                                           const SmilesWriteParamsData &data) {
  check_removals(mol);
  auto lock = random_lock(data.do_random);
  return str(RDKit::MolToSmiles(mol, write_params(data)));
}

// As MolFromSmartsHelper in Code/GraphMol/Wrap/rdmolfiles.cpp.
std::unique_ptr<RWMol> mol_from_smarts(rust::Str smarts,
                                       const SmartsParserParamsData &data) {
  check_debug_parse(data.debug_parse);
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

// MolToSmarts ignores doRandom (SmartsWrite.cpp passes false to Canon), so it
// needs no random_mutex.
std::unique_ptr<std::string> mol_to_smarts(const RWMol &mol,
                                           const SmilesWriteParamsData &data) {
  return str(RDKit::MolToSmarts(mol, write_params(data)));
}

std::unique_ptr<Atom> atom_from_smiles(rust::Str smiles) {
  std::string text(smiles);
  auto atom = RDKit::v2::SmilesParse::AtomFromSmiles(text);
  if (!atom) {
    throw RDKit::SmilesParse::SmilesParseException(
        "could not parse SMILES atom: " + text);
  }
  return atom;
}

std::unique_ptr<Atom> atom_from_smarts(rust::Str smarts) {
  std::string text(smarts);
  auto atom = RDKit::v2::SmilesParse::AtomFromSmarts(text);
  if (!atom) {
    throw RDKit::SmilesParse::SmilesParseException(
        "could not parse SMARTS atom: " + text);
  }
  return atom;
}

std::unique_ptr<Bond> bond_from_smiles(rust::Str smiles) {
  std::string text(smiles);
  auto bond = RDKit::v2::SmilesParse::BondFromSmiles(text);
  if (!bond) {
    throw RDKit::SmilesParse::SmilesParseException(
        "could not parse SMILES bond: " + text);
  }
  return bond;
}

std::unique_ptr<Bond> bond_from_smarts(rust::Str smarts) {
  std::string text(smarts);
  auto bond = RDKit::v2::SmilesParse::BondFromSmarts(text);
  if (!bond) {
    throw RDKit::SmilesParse::SmilesParseException(
        "could not parse SMARTS bond: " + text);
  }
  return bond;
}

std::unique_ptr<std::string> mol_to_cx_smiles(
    const RWMol &mol, const SmilesWriteParamsData &data, std::uint32_t flags,
    RestoreBondDirOption restore_bond_dirs) {
  check_removals(mol);
  auto lock = random_lock(data.do_random);
  return str(
      RDKit::MolToCXSmiles(mol, write_params(data), flags, restore_bond_dirs));
}

std::unique_ptr<std::string> mol_to_cx_smarts(const RWMol &mol,
                                              bool isomeric_smiles) {
  return str(RDKit::MolToCXSmarts(mol, isomeric_smiles));
}

// As MolFragmentToSmilesHelper1 in Code/GraphMol/Wrap/rdmolfiles.cpp, with
// the indices checked: RDKit sets them in bitsets without a range check.
// Python (pythonObjectToVect in Code/RDBoost/Wrap.h) treats empty bondsToUse,
// atomSymbols, and bondSymbols lists as None; here an empty bondsToUse is an
// empty set, and empty symbol lists must still match the atoms or bonds.
std::unique_ptr<std::string> mol_fragment_to_smiles(
    const RWMol &mol, const SmilesWriteParamsData &data,
    rust::Slice<const std::uint32_t> atoms_to_use, bool has_bonds,
    rust::Slice<const std::uint32_t> bonds_to_use, bool has_atom_symbols,
    const rust::Vec<rust::String> &atom_symbols, bool has_bond_symbols,
    const rust::Vec<rust::String> &bond_symbols, bool cx) {
  if (atoms_to_use.empty()) {
    throw ValueErrorException("atomsToUse must not be empty");
  }
  auto atoms = atom_indices(mol, atoms_to_use);
  std::optional<std::vector<int>> bonds;
  if (has_bonds) {
    bonds = bond_indices(mol, bonds_to_use);
  }
  std::optional<std::vector<std::string>> asymbols, bsymbols;
  if (has_atom_symbols) {
    if (atom_symbols.size() != mol.getNumAtoms()) {
      throw ValueErrorException(
          "length of atom symbol list != number of atoms");
    }
    asymbols = strings(atom_symbols);
  }
  if (has_bond_symbols) {
    if (bond_symbols.size() != mol.getNumBonds()) {
      throw ValueErrorException(
          "length of bond symbol list != number of bonds");
    }
    bsymbols = strings(bond_symbols);
  }
  auto params = write_params(data);
  // Without bonds, a root atom has RDKit split the molecule into fragments.
  check_removals(mol);
  auto lock = random_lock(data.do_random);
  const auto *b = bonds ? &*bonds : nullptr;
  const auto *as = asymbols ? &*asymbols : nullptr;
  const auto *bs = bsymbols ? &*bsymbols : nullptr;
  return str(cx ? RDKit::MolFragmentToCXSmiles(mol, params, atoms, b, as, bs)
                : RDKit::MolFragmentToSmiles(mol, params, atoms, b, as, bs));
}

// As molFragmentToSmarts in Code/GraphMol/Wrap/rdmolfiles.cpp, with the
// indices checked. Python (pythonObjectToVect in Code/RDBoost/Wrap.h) treats
// an empty bondsToUse list as None, using every bond; here it uses none.
std::unique_ptr<std::string> mol_fragment_to_smarts(
    const RWMol &mol, rust::Slice<const std::uint32_t> atoms_to_use,
    bool has_bonds, rust::Slice<const std::uint32_t> bonds_to_use,
    bool isomeric_smarts, bool cx) {
  if (atoms_to_use.empty()) {
    throw ValueErrorException("atomsToUse argument must be non-empty");
  }
  auto atoms = atom_indices(mol, atoms_to_use);
  std::optional<std::vector<int>> bonds;
  if (has_bonds) {
    bonds = bond_indices(mol, bonds_to_use);
  }
  const auto *b = bonds ? &*bonds : nullptr;
  return str(cx ? RDKit::MolFragmentToCXSmarts(mol, atoms, b, isomeric_smarts)
                : RDKit::MolFragmentToSmarts(mol, atoms, b, isomeric_smarts));
}

std::unique_ptr<std::vector<std::string>> mol_to_random_smiles_vect(
    const RWMol &mol, std::uint32_t num_smiles, std::uint32_t random_seed,
    bool isomeric_smiles, bool kekule_smiles, bool all_bonds_explicit,
    bool all_hs_explicit) {
  check_removals(mol);
  // RDKit seeds the global generator and draws from it for each SMILES.
  std::lock_guard<std::mutex> lock(random_mutex());
  return std::make_unique<std::vector<std::string>>(
      RDKit::MolToRandomSmilesVect(mol, num_smiles, random_seed,
                                   isomeric_smiles, kekule_smiles,
                                   all_bonds_explicit, all_hs_explicit));
}

// As CanonicalRankAtoms in Code/GraphMol/Wrap/rdmolfiles.cpp.
rust::Vec<std::uint32_t> canonical_rank_atoms(const RWMol &mol, bool break_ties,
                                              bool include_chirality,
                                              bool include_isotopes,
                                              bool include_atom_maps,
                                              bool include_chiral_presence) {
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  const bool includeStereoGroups = true;
  const bool useNonStereoRanks = false;
  RDKit::Canon::rankMolAtoms(mol, ranks, break_ties, include_chirality,
                             include_isotopes, include_atom_maps,
                             include_chiral_presence, includeStereoGroups,
                             useNonStereoRanks);
  return to_vec(ranks);
}

// As CanonicalRankAtomsInFragment in Code/GraphMol/Wrap/rdmolfiles.cpp,
// which sets the indices in bitsets; they are checked here. Boost.Python
// passes includeChiralPresence where Canon::rankFragmentAtoms takes
// includeAtomMaps and drops includeAtomMaps (rdmolfiles.cpp:509-512,
// new_canon.h:917); here each goes where RDKit's signature puts it. Python
// also treats an empty atomSymbols list as None.
rust::Vec<std::int32_t> canonical_rank_atoms_in_fragment(
    const RWMol &mol, rust::Slice<const std::uint32_t> atoms_to_use,
    bool has_bonds, rust::Slice<const std::uint32_t> bonds_to_use,
    bool has_atom_symbols, const rust::Vec<rust::String> &atom_symbols,
    bool break_ties, bool include_chirality, bool include_isotopes,
    bool include_atom_maps, bool include_chiral_presence) {
  if (atoms_to_use.empty()) {
    throw ValueErrorException("atomsToUse must not be empty");
  }
  auto avect = atom_indices(mol, atoms_to_use);
  std::vector<int> bvect;
  if (has_bonds) {
    bvect = bond_indices(mol, bonds_to_use);
  }
  std::optional<std::vector<std::string>> asymbols;
  if (has_atom_symbols) {
    if (atom_symbols.size() != mol.getNumAtoms()) {
      throw ValueErrorException(
          "length of atom symbol list != number of atoms");
    }
    asymbols = strings(atom_symbols);
  }
  boost::dynamic_bitset<> atoms(mol.getNumAtoms());
  for (auto idx : avect) {
    atoms[idx] = true;
  }
  boost::dynamic_bitset<> bonds(mol.getNumBonds());
  for (auto idx : bvect) {
    bonds[idx] = true;
  }
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  RDKit::Canon::rankFragmentAtoms(mol, ranks, atoms, bonds,
                                  asymbols ? &*asymbols : nullptr, break_ties,
                                  include_chirality, include_isotopes,
                                  include_atom_maps, include_chiral_presence);
  rust::Vec<std::int32_t> result;
  result.reserve(mol.getNumAtoms());
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    result.push_back(atoms[i] ? static_cast<std::int32_t>(ranks[i]) : -1);
  }
  return result;
}

void canonicalize_enhanced_stereo(RWMol &mol) {
  RDKit::Canon::canonicalizeEnhancedStereo(mol);
}

// RDKit 2024.03 reads integer properties as int64_t.
void create_property_list(RWMol &mol, std::uint8_t kind, rust::Str prop_name,
                          rust::Str missing_value_marker,
                          std::uint32_t line_size) {
  namespace FPU = RDKit::FileParserUtils;
  std::string prop(prop_name), marker(missing_value_marker);
  switch (kind) {
    case 0:
      check_first_list_value<std::int64_t>(mol.atoms(), prop, marker,
                                           line_size);
      FPU::createAtomIntPropertyList(mol, prop, marker, line_size);
      break;
    case 1:
      check_first_list_value<double>(mol.atoms(), prop, marker, line_size);
      FPU::createAtomDoublePropertyList(mol, prop, marker, line_size);
      break;
    case 2:
      check_first_list_value<bool>(mol.atoms(), prop, marker, line_size);
      FPU::createAtomBoolPropertyList(mol, prop, marker, line_size);
      break;
    default:
      check_first_list_value<std::string>(mol.atoms(), prop, marker,
                                          line_size);
      FPU::createAtomStringPropertyList(mol, prop, marker, line_size);
      break;
  }
}

// The mol-block readers check their input as U12 describes (streams.h); a file
// is read through a checked stream, so RDKit reads only checked bytes.
std::unique_ptr<RWMol> mol_from_mol_block(rust::Str mol_block, bool sanitize,
                                          bool remove_hs, bool strict_parsing) {
  FP::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.strictParsing = strict_parsing;
  std::string text(mol_block);
  check_input_text(text);
  return parsed(FP::MolFromMolBlock(text, params), "mol block");
}

// As MolFromMolFile in MolFileParser.cpp, through a checked stream.
std::unique_ptr<RWMol> mol_from_mol_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs,
                                         bool strict_parsing) {
  FP::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.strictParsing = strict_parsing;
  auto name = to_path(path);
  auto in = open_checked_file(name, false);
  unsigned int line = 0;
  std::unique_ptr<RWMol> mol;
  try {
    mol = FP::MolFromMolDataStream(*in, line, params);
  } catch (...) {
    in->check();  // a refused input explains RDKit's failure
    throw;
  }
  in->check();
  return parsed(std::move(mol), "mol file " + name);
}

namespace {

// The mol-block and Marvin writers compute 2D coordinates for stereo when
// the molecule has no conformer (prepareMol in MolFileWriter.cpp,
// MolToMrvBlock in MarvinWriter.cpp), with compute2DCoords, which trusts the
// molecule (depiction_checks.h).
std::unique_lock<std::mutex> lock_if_depicting(const RWMol &mol,
                                               bool include_stereo) {
  if (include_stereo && mol.getNumConformers() == 0) {
    auto lock = lock_depiction();
    check_compute_2d_coords_mol(mol, false);
    return lock;
  }
  return {};
}

std::unique_lock<std::mutex> lock_if_depicting(
    const RWMol &mol, const MolWriterParamsData &data) {
  return lock_if_depicting(mol, data.include_stereo);
}

}  // namespace

std::unique_ptr<std::string> mol_to_mol_block(const RWMol &mol,
                                              const MolWriterParamsData &data) {
  auto lock = lock_if_depicting(mol, data);
  return str(RDKit::MolToMolBlock(mol, writer_params(data), data.conf_id));
}

void mol_to_mol_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     const MolWriterParamsData &data) {
  auto lock = lock_if_depicting(mol, data);
  RDKit::MolToMolFile(mol, to_path(path), writer_params(data), data.conf_id);
}

std::unique_ptr<std::string> mol_to_v2k_mol_block(
    const RWMol &mol, const MolWriterParamsData &data) {
  auto lock = lock_if_depicting(mol, data);
  return str(RDKit::MolToV2KMolBlock(mol, writer_params(data), data.conf_id));
}

std::unique_ptr<std::string> mol_to_v3k_mol_block(
    const RWMol &mol, const MolWriterParamsData &data) {
  auto lock = lock_if_depicting(mol, data);
  return str(RDKit::MolToV3KMolBlock(mol, writer_params(data), data.conf_id));
}

void mol_to_v3k_mol_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                         const MolWriterParamsData &data) {
  auto lock = lock_if_depicting(mol, data);
  RDKit::MolToV3KMolFile(mol, to_path(path), writer_params(data), data.conf_id);
}


namespace {

std::unique_ptr<RWMol> tpl_from_text(const std::string &text, bool sanitize,
                                     bool skip_first_conf,
                                     const std::string &what) {
  check_tpl(text);
  std::istringstream in(text);
  unsigned int line = 0;
  FP::TPLParserParams params;
  params.sanitize = sanitize;
  params.skipFirstConf = skip_first_conf;
  return parsed(FP::MolFromTPLDataStream(in, line, params), what);
}

}  // namespace

// As MolFromTPLBlock in Code/GraphMol/Wrap/rdmolfiles.cpp.
std::unique_ptr<RWMol> mol_from_tpl_block(rust::Str tpl_block, bool sanitize,
                                          bool skip_first_conf) {
  return tpl_from_text(std::string(tpl_block), sanitize, skip_first_conf,
                       "TPL block");
}

// The file is read once and checked; RDKit parses those bytes (U2).
std::unique_ptr<RWMol> mol_from_tpl_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool skip_first_conf) {
  auto name = to_path(path);
  return tpl_from_text(read_file(name), sanitize, skip_first_conf,
                       "TPL file " + name);
}

std::unique_ptr<std::string> mol_to_tpl_block(const RWMol &mol,
                                              rust::Str partial_charge_prop,
                                              bool write_first_conf_twice) {
  return str(RDKit::MolToTPLText(mol, std::string(partial_charge_prop),
                                 write_first_conf_twice));
}

void mol_to_tpl_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     rust::Str partial_charge_prop,
                     bool write_first_conf_twice) {
  RDKit::MolToTPLFile(mol, to_path(path), std::string(partial_charge_prop),
                      write_first_conf_twice);
}

namespace {

FP::Mol2ParserParams mol2_params(bool sanitize, bool remove_hs,
                                 bool cleanup_substructures) {
  FP::Mol2ParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  params.variant = FP::Mol2Type::CORINA;
  params.cleanupSubstructures = cleanup_substructures;
  return params;
}

}  // namespace

std::unique_ptr<RWMol> mol_from_mol2_block(rust::Str mol2_block, bool sanitize,
                                           bool remove_hs,
                                           bool cleanup_substructures) {
  return parsed(FP::MolFromMol2Block(
                    std::string(mol2_block),
                    mol2_params(sanitize, remove_hs, cleanup_substructures)),
                "Mol2 block");
}

std::unique_ptr<RWMol> mol_from_mol2_file(rust::Slice<const std::uint8_t> path,
                                          bool sanitize, bool remove_hs,
                                          bool cleanup_substructures) {
  auto name = to_path(path);
  return parsed(FP::MolFromMol2File(name, mol2_params(sanitize, remove_hs,
                                                      cleanup_substructures)),
                "Mol2 file " + name);
}

namespace {

std::unique_ptr<RWMol> mrv_from_text(const std::string &text, bool sanitize,
                                     bool remove_hs, const std::string &what) {
  check_mrv(text);
  std::istringstream in(text);
  RDKit::v2::MarvinParser::MrvParserParams params;
  params.sanitize = sanitize;
  params.removeHs = remove_hs;
  return parsed(RDKit::v2::MarvinParser::MolFromMrvDataStream(in, params),
                what);
}

}  // namespace

std::unique_ptr<RWMol> mol_from_mrv_block(rust::Str mrv_block, bool sanitize,
                                          bool remove_hs) {
  return mrv_from_text(std::string(mrv_block), sanitize, remove_hs,
                       "Marvin block");
}

// The file is read once and checked; RDKit parses those bytes (U2).
std::unique_ptr<RWMol> mol_from_mrv_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs) {
  auto name = to_path(path);
  return mrv_from_text(read_file(name), sanitize, remove_hs,
                       "Marvin file " + name);
}

std::unique_ptr<std::string> mol_to_mrv_block(const RWMol &mol,
                                              bool include_stereo,
                                              std::int32_t conf_id,
                                              bool kekulize,
                                              bool pretty_print) {
  auto lock = lock_if_depicting(mol, include_stereo);
  return str(RDKit::MolToMrvBlock(mol, include_stereo, conf_id, kekulize,
                                  pretty_print));
}

void mol_to_mrv_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     bool include_stereo, std::int32_t conf_id, bool kekulize,
                     bool pretty_print) {
  auto lock = lock_if_depicting(mol, include_stereo);
  RDKit::MolToMrvFile(mol, to_path(path), include_stereo, conf_id, kekulize,
                      pretty_print);
}

// As MolFromXYZBlock in Code/GraphMol/Wrap/rdmolfiles.cpp, which reads an
// empty block as an empty molecule.
std::unique_ptr<RWMol> mol_from_xyz_block(rust::Str xyz_block) {
  std::istringstream in{std::string(xyz_block)};
  return parsed(FP::MolFromXYZDataStream(in), "XYZ block");
}

std::unique_ptr<RWMol> mol_from_xyz_file(rust::Slice<const std::uint8_t> path) {
  auto name = to_path(path);
  return parsed(FP::MolFromXYZFile(name), "XYZ file " + name);
}

std::unique_ptr<std::string> mol_to_xyz_block(const RWMol &mol,
                                              std::int32_t conf_id,
                                              std::uint32_t precision) {
  check_precision(precision);
  return str(RDKit::MolToXYZBlock(mol, conf_id, precision));
}

void mol_to_xyz_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     std::int32_t conf_id, std::uint32_t precision) {
  check_precision(precision);
  RDKit::MolToXYZFile(mol, to_path(path), conf_id, precision);
}

namespace {

// As Python's MolFromPDBBlock and MolFromPDBFile, through
// MolFromPDBDataStream, which reads up to the END record.
std::unique_ptr<RWMol> pdb_from_text(const std::string &text, bool sanitize,
                                     bool remove_hs, std::uint32_t flavor,
                                     bool proximity_bonding,
                                     const std::string &what) {
  if (proximity_bonding) {
    std::istringstream in(text);
    auto probe =
        FP::MolFromPDBDataStream(in, pdb_params(false, false, flavor, false));
    if (probe) {
      check_pdb_coordinates(*probe);
    }
  }
  std::istringstream in(text);
  return parsed(
      FP::MolFromPDBDataStream(
          in, pdb_params(sanitize, remove_hs, flavor, proximity_bonding)),
      what);
}

}  // namespace

std::unique_ptr<RWMol> mol_from_pdb_block(rust::Str mol_block, bool sanitize,
                                          bool remove_hs, std::uint32_t flavor,
                                          bool proximity_bonding) {
  return pdb_from_text(std::string(mol_block), sanitize, remove_hs, flavor,
                       proximity_bonding, "PDB block");
}

// The file is read once and checked; RDKit parses those bytes (U2).
std::unique_ptr<RWMol> mol_from_pdb_file(rust::Slice<const std::uint8_t> path,
                                         bool sanitize, bool remove_hs,
                                         std::uint32_t flavor,
                                         bool proximity_bonding) {
  auto name = to_path(path);
  return pdb_from_text(read_file(name), sanitize, remove_hs, flavor,
                       proximity_bonding, "PDB file " + name);
}

std::unique_ptr<std::string> mol_to_pdb_block(const RWMol &mol,
                                              std::int32_t conf_id,
                                              std::uint32_t flavor) {
  return str(RDKit::MolToPDBBlock(mol, conf_id, flavor));
}

void mol_to_pdb_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     std::int32_t conf_id, std::uint32_t flavor) {
  RDKit::MolToPDBFile(mol, to_path(path), conf_id, flavor);
}

// RDKit's CML writer (molToPTree in CMLWriter.cpp) sums the formal charges
// in int (U9).
std::unique_ptr<std::string> mol_to_cml_block(const RWMol &mol,
                                              std::int32_t conf_id,
                                              bool kekulize) {
  check_formal_charge_sum(mol);
  return str(RDKit::MolToCMLBlock(mol, conf_id, kekulize));
}

void mol_to_cml_file(const RWMol &mol, rust::Slice<const std::uint8_t> path,
                     std::int32_t conf_id, bool kekulize) {
  check_formal_charge_sum(mol);
  RDKit::MolToCMLFile(mol, to_path(path), conf_id, kekulize);
}

// As MolFromSVG in Code/GraphMol/Wrap/rdmolfiles.cpp.
std::unique_ptr<RWMol> mol_from_rdkit_svg(rust::Str svg, bool sanitize,
                                          bool remove_hs) {
  std::string text(svg);
  check_svg(text);
  return parsed(
      std::unique_ptr<RWMol>(RDKit::RDKitSVGToMol(text, sanitize, remove_hs)),
      "RDKit SVG");
}

std::unique_ptr<RWMol> mol_from_sequence(rust::Str text, bool sanitize,
                                         std::int32_t flavor) {
  return parsed(std::unique_ptr<RWMol>(
                    RDKit::SequenceToMol(std::string(text), sanitize, flavor)),
                "sequence");
}

// RDKit's FASTAToMol (SequenceParsers.cpp) sets the title on SequenceToMol's
// result without checking it, dereferencing null for a titled sequence it
// cannot parse. This is FASTAToMol with that check.
std::unique_ptr<RWMol> mol_from_fasta(rust::Str text, bool sanitize,
                                      std::int32_t flavor) {
  std::string seq(text);
  const char *ptr = seq.c_str();
  std::string title;
  if (ptr[0] == '>') {
    ptr++;
    while (*ptr && *ptr != '\n' && *ptr != '\r') {
      title += *ptr++;
    }
  }
  auto mol = parsed(
      std::unique_ptr<RWMol>(RDKit::SequenceToMol(ptr, sanitize, flavor)),
      "FASTA");
  if (!title.empty()) {
    mol->setProp(RDKit::common_properties::_Name, title);
  }
  return mol;
}

std::unique_ptr<RWMol> mol_from_helm(rust::Str text, bool sanitize) {
  return parsed(
      std::unique_ptr<RWMol>(RDKit::HELMToMol(std::string(text), sanitize)),
      "HELM");
}

std::unique_ptr<std::string> mol_to_sequence(const RWMol &mol) {
  return str(RDKit::MolToSequence(mol));
}

std::unique_ptr<std::string> mol_to_fasta(const RWMol &mol) {
  return str(RDKit::MolToFASTA(mol));
}

std::unique_ptr<std::string> mol_to_helm(const RWMol &mol) {
  return str(RDKit::MolToHELM(mol));
}

}  // namespace rdkit_cxx
