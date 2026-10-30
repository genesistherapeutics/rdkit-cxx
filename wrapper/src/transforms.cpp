// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/transforms.h"

#include <algorithm>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <RDGeneral/Exceptions.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/transforms.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

// RDKit returns new ROMols; Rust molecules own RWMols (W1). Null stays null
// for the functions that return None. None of these calls computes a matrix
// cache, but the copy drops any the input carried (U6).
std::unique_ptr<RWMol> to_rwmol(RDKit::ROMol *mol) {
  std::unique_ptr<RDKit::ROMol> owned(mol);
  if (!owned) {
    return nullptr;
  }
  auto result = std::make_unique<RWMol>(*owned);
  clear_matrix_caches(*result);
  return result;
}

std::unique_ptr<RWMol> to_rwmol(std::unique_ptr<RDKit::ROMol> mol) {
  return to_rwmol(mol.release());
}

std::vector<RDKit::MatchVectType> to_matches(
    rust::Slice<const IndexList> matches) {
  // As seqOfSeqsToMatchVectTypeVect in Code/GraphMol/Wrap/MolOps.cpp.
  if (matches.empty()) {
    throw ValueErrorException("matches must not be empty");
  }
  std::vector<RDKit::MatchVectType> result;
  for (const auto &match : matches) {
    if (match.indices.empty()) {
      throw ValueErrorException("tuples in matches must not be empty");
    }
    RDKit::MatchVectType pairs;
    for (std::size_t i = 0; i < match.indices.size(); ++i) {
      // ScoreMatchesByDegreeOfCoreSubstitution range-checks both indices
      // through getAtomWithIdx; an index past INT_MAX becomes negative here
      // and fails that check.
      pairs.emplace_back(static_cast<int>(i),
                         static_cast<int>(match.indices[i]));
    }
    result.push_back(std::move(pairs));
  }
  return result;
}

std::vector<std::pair<unsigned, unsigned>> to_label_pairs(
    rust::Slice<const IndexPair> labels) {
  std::vector<std::pair<unsigned, unsigned>> result;
  for (const auto &label : labels) {
    result.emplace_back(label.first, label.second);
  }
  return result;
}

// As fragmentOnBondsHelper and fragmentOnSomeBondsHelper in
// Code/GraphMol/Wrap/MolOps.cpp, which check the bond indices before
// fragmentOnBonds calls getBondWithIdx, and the bondTypes size, which RDKit
// only checks with PRECONDITION. Python treats an empty bondTypes list as None
// and leaves dummyLabels' size unchecked; here an empty list is an empty list,
// whose size is checked, and dummyLabels must cover the bonds too.
std::vector<unsigned> check_fragment_args(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    bool has_dummy_labels, rust::Slice<const IndexPair> dummy_labels,
    bool has_bond_types, rust::Slice<const BondType> bond_types) {
  if (bond_indices.empty()) {
    throw ValueErrorException("empty bond indices");
  }
  std::vector<unsigned> indices;
  for (auto idx : bond_indices) {
    check_index(idx, mol.getNumBonds());
    indices.push_back(idx);
  }
  // RDKit 2024.09's fragmentOnBonds finds every bond before it removes any,
  // and reads a repeated one after removing it; RDKit 2026.03.5 refuses them
  // (rdkit/rdkit#9430).
  std::vector<unsigned> sorted(indices);
  std::sort(sorted.begin(), sorted.end());
  if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
    throw ValueErrorException("bondIndices contains duplicates");
  }
  if (has_dummy_labels && dummy_labels.size() < bond_indices.size()) {
    throw ValueErrorException("dummyLabels shorter than bondIndices");
  }
  if (has_bond_types && bond_types.size() != bond_indices.size()) {
    throw ValueErrorException("bondTypes shorter than bondIndices");
  }
  return indices;
}

RDKit::MolzipParams to_params(const MolzipParamsData &data) {
  RDKit::MolzipParams params;
  params.label = data.label;
  params.atomSymbols.clear();
  for (const auto &symbol : data.atom_symbols) {
    params.atomSymbols.emplace_back(symbol);
  }
  params.enforceValenceRules = data.enforce_valence_rules;
  params.generateCoordinates = data.generate_coordinates;
  return params;
}

// Whether the atom carries a label, as get_label in MolZip.cpp.
bool has_zip_label(const RDKit::Atom *atom, const RDKit::MolzipParams &params) {
  switch (params.label) {
    case RDKit::MolzipLabel::AtomMapNumber:
      return atom->getAtomicNum() == 0 && atom->getAtomMapNum() != 0;
    case RDKit::MolzipLabel::Isotope:
      return atom->getAtomicNum() == 0 && atom->getIsotope() != 0;
    case RDKit::MolzipLabel::FragmentOnBonds:
      return atom->getAtomicNum() == 0;
    case RDKit::MolzipLabel::AtomType:
      return std::find(params.atomSymbols.begin(), params.atomSymbols.end(),
                       atom->getSymbol()) != params.atomSymbols.end();
    default:
      return false;
  }
}

// molzip takes the atom bonded to each labeled atom from get_other_atom,
// which is null unless the labeled atom has exactly one neighbor. RDKit
// dereferences that null pointer with the FragmentOnBonds label, and when
// zipping a list, which tracks attachment indices, so such an atom is
// refused there. Zipping a list also dereferences a null dummy for every
// linker ([*:1][*:2], two labeled atoms bonded to each other), so a list may
// not hold one.
void check_zip_labels(const RDKit::ROMol &mol,
                      const RDKit::MolzipParams &params, bool zipping_list) {
  if (params.label != RDKit::MolzipLabel::FragmentOnBonds && !zipping_list) {
    return;
  }
  for (const auto *atom : mol.atoms()) {
    if (!has_zip_label(atom, params)) {
      continue;
    }
    if (atom->getDegree() != 1) {
      throw ValueErrorException("molzip: labeled atom " +
                                std::to_string(atom->getIdx()) +
                                " must have exactly one neighbor");
    }
    if (zipping_list && params.label != RDKit::MolzipLabel::FragmentOnBonds) {
      for (const auto *nbr : mol.atomNeighbors(atom)) {
        if (has_zip_label(nbr, params)) {
          throw ValueErrorException(
              "molzip: RDKit cannot zip a list that holds a linker (atoms " +
              std::to_string(atom->getIdx()) + " and " +
              std::to_string(nbr->getIdx()) + ")");
        }
      }
    }
  }
}

// fragmentOnBRICSBonds fills its pattern globals on first use, with no lock,
// and then matches every call against those shared query molecules (U4).
// It comes before control_c_mutex, so a closure that holds that one must
// not wait for it.
std::mutex &brics_mutex() {
  check_not_in_control_c_callback();
  static std::mutex mutex;
  return mutex;
}

}  // namespace

// deleteSubstructs, replaceSubstructs, replaceSidechains, replaceCore with a
// query, and fragmentOnBRICSBonds match substructures (control_c_lock.h).
// These, MurckoDecompose, and the fragmenters remove atoms or bonds from a
// copy of the molecule, which trusts what check_removals checks.
std::unique_ptr<RWMol> delete_substructs(const RWMol &mol, const RWMol &query,
                                         bool only_frags, bool use_chirality) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_rwmol(
      RDKit::deleteSubstructs(mol, query, only_frags, use_chirality));
}

// As replaceSubstructures in Code/GraphMol/Wrap/MolOps.cpp.
std::unique_ptr<MolList> replace_substructs(
    const RWMol &mol, const RWMol &query, const RWMol &replacement,
    bool replace_all, std::uint32_t replacement_connection_point,
    bool use_chirality) {
  // The copy gains the replacement's substance groups before RDKit removes
  // the match.
  check_removals({&mol, &replacement});
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  auto results =
      RDKit::replaceSubstructs(mol, query, replacement, replace_all,
                               replacement_connection_point, use_chirality);
  auto list = std::make_unique<MolList>();
  for (const auto &result : results) {
    auto copy = std::make_unique<RWMol>(*result);
    clear_matrix_caches(*copy);
    list->push(std::move(copy));
  }
  return list;
}

std::unique_ptr<RWMol> replace_sidechains(const RWMol &mol,
                                          const RWMol &core_query,
                                          bool use_chirality) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_rwmol(RDKit::replaceSidechains(mol, core_query, use_chirality));
}

// As replaceCoreHelper in Code/GraphMol/Wrap/MolOps.cpp, with the match as
// (core atom, molecule atom) pairs; replaceCore checks both ranges. An index
// past INT_MAX becomes negative here and fails that check.
std::unique_ptr<RWMol> replace_core_with_matches(
    const RWMol &mol, const RWMol &core, rust::Slice<const IndexPair> matches,
    bool replace_dummies, bool label_by_index, bool require_dummy_match) {
  RDKit::MatchVectType match;
  for (const auto &pair : matches) {
    match.emplace_back(static_cast<int>(pair.first),
                       static_cast<int>(pair.second));
  }
  // RDKit 2024.09 reads the match of each core atom that shares its molecule
  // atom with another core atom's neighbor without checking that the match
  // has one; RDKit 2026.03 throws this (rdkit/rdkit#9450). Indices out of
  // range are left to RDKit, which refuses them first.
  std::vector<unsigned> mapped(mol.getNumAtoms(), 0);
  bool in_range = true;
  for (const auto &[core_idx, mol_idx] : match) {
    if (core_idx < 0 || core_idx >= static_cast<int>(core.getNumAtoms()) ||
        mol_idx < 0 || mol_idx >= static_cast<int>(mol.getNumAtoms())) {
      in_range = false;
      break;
    }
    ++mapped[mol_idx];
  }
  for (const auto &[core_idx, mol_idx] : match) {
    if (!in_range || mapped[mol_idx] < 2) {
      continue;
    }
    const auto *core_atom = core.getAtomWithIdx(core_idx);
    if (core_atom->getDegree() != 1) {
      continue;  // RDKit's invariant throws
    }
    auto neighbor = static_cast<int>(
        core[*core.getAtomNeighbors(core_atom).first]->getIdx());
    if (std::none_of(match.begin(), match.end(),
                     [neighbor](const std::pair<int, int> &p) {
                       return p.first == neighbor;
                     })) {
      throw ValueErrorException(
          "Supplied MatchVect is missing the neighbor of a multiply mapped "
          "core atom");
    }
  }
  check_removals(mol);
  return to_rwmol(RDKit::replaceCore(mol, core, match, replace_dummies,
                                     label_by_index, require_dummy_match));
}

std::unique_ptr<RWMol> replace_core(const RWMol &mol, const RWMol &core_query,
                                    bool replace_dummies, bool label_by_index,
                                    bool require_dummy_match,
                                    bool use_chirality) {
  check_removals(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_rwmol(RDKit::replaceCore(mol, core_query, replace_dummies,
                                     label_by_index, require_dummy_match,
                                     use_chirality));
}

// MurckoDecompose computes the distance matrix on its const input, with
// force, and reads its paths (U6).
std::unique_ptr<RWMol> murcko_decompose(const RWMol &mol) {
  check_removals(mol);
  MatrixCacheGuard guard(mol);
  return to_rwmol(RDKit::MurckoDecompose(mol));
}

std::unique_ptr<RWMol> combine_mols(const RWMol &mol1, const RWMol &mol2,
                                    const Point3D &offset) {
  return to_rwmol(RDKit::combineMols(
      mol1, mol2, RDGeom::Point3D(offset.x, offset.y, offset.z)));
}

// As getMostSubstitutedCoreMatchHelper in Code/GraphMol/Wrap/MolOps.cpp.
// RDKit returns a reference into the matches, which live until the copy.
rust::Vec<std::uint32_t> get_most_substituted_core_match(
    const RWMol &mol, const RWMol &core, rust::Slice<const IndexList> matches) {
  auto match_vects = to_matches(matches);
  const auto &best = RDKit::getMostSubstitutedCoreMatch(mol, core, match_vects);
  rust::Vec<std::uint32_t> result;
  for (const auto &pair : best) {
    result.push_back(static_cast<std::uint32_t>(pair.second));
  }
  return result;
}

rust::Vec<IndexList> sort_matches_by_degree_of_core_substitution(
    const RWMol &mol, const RWMol &core, rust::Slice<const IndexList> matches) {
  auto match_vects = to_matches(matches);
  rust::Vec<IndexList> result;
  for (const auto &match :
       RDKit::sortMatchesByDegreeOfCoreSubstitution(mol, core, match_vects)) {
    rust::Vec<std::uint32_t> indices;
    for (const auto &pair : match) {
      indices.push_back(static_cast<std::uint32_t>(pair.second));
    }
    result.push_back(IndexList{std::move(indices)});
  }
  return result;
}

// As addRecursiveQuery in Code/GraphMol/Wrap/MolOps.cpp.
void add_recursive_query(RWMol &mol, const RWMol &query, std::uint32_t atom_idx,
                         bool preserve_existing_query) {
  if (atom_idx >= mol.getNumAtoms()) {
    throw ValueErrorException("atom index exceeds mol.GetNumAtoms()");
  }
  auto q =
      std::make_unique<RDKit::RecursiveStructureQuery>(new RDKit::ROMol(query));
  auto *atom = mol.getAtomWithIdx(atom_idx);
  if (!atom->hasQuery()) {
    RDKit::QueryAtom query_atom(*atom);
    mol.replaceAtom(atom_idx, &query_atom);
    atom = mol.getAtomWithIdx(atom_idx);
  }
  if (!preserve_existing_query) {
    atom->setQuery(q.release());
  } else {
    atom->expandQuery(q.release(), Queries::COMPOSITE_AND);
  }
}

// addRecursiveQueries copies each query molecule into its query.
void add_recursive_queries(RWMol &mol, BorrowedMols &queries,
                           rust::Str prop_name) {
  std::map<std::string, RDKit::ROMOL_SPTR> map;
  for (std::size_t i = 0; i < queries.names().size(); ++i) {
    map[queries.names()[i]] = queries.mols()[i];
  }
  check_query_props(mol, std::string(prop_name));
  RDKit::addRecursiveQueries(mol, map, std::string(prop_name));
}

namespace {

// As parseQueryDefFile (ChemTransforms.cpp) reads each line's SMARTS field,
// which it parses with SmartsToMol, reading CX extensions (file_checks.h).
void check_query_defs(const std::string &text, const std::string &delim,
                      const std::string &comment_text,
                      std::uint32_t smarts_column) {
  std::istringstream lines(text);
  std::string line;
  boost::char_separator<char> sep(delim.c_str());
  while (std::getline(lines, line)) {
    if (line.empty() || line.find(comment_text) == 0) {
      continue;
    }
    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
    std::uint32_t pos = 0;
    for (auto token = tokens.begin(); token != tokens.end(); ++token, ++pos) {
      if (pos == smarts_column) {
        check_smiles_cx(boost::trim_copy_if(*token, boost::is_any_of(" \t")),
                        true, true);
      }
    }
  }
}

std::unique_ptr<MolList> query_defs(const std::string &text, bool standardize,
                                    const std::string &delim,
                                    const std::string &comment_text,
                                    std::uint32_t name_column,
                                    std::uint32_t smarts_column) {
  check_query_defs(text, delim, comment_text, smarts_column);
  std::map<std::string, RDKit::ROMOL_SPTR> defs;
  RDKit::parseQueryDefText(text, defs, standardize, delim, comment_text,
                           name_column, smarts_column);
  auto list = std::make_unique<MolList>();
  for (const auto &[name, mol] : defs) {
    list->push_named(name, std::make_unique<RWMol>(*mol));
  }
  return list;
}

}  // namespace

std::unique_ptr<MolList> parse_mol_query_def_file(
    rust::Slice<const std::uint8_t> path, bool standardize, rust::Str delimiter,
    rust::Str comment, std::uint32_t name_column, std::uint32_t smarts_column) {
  // Read once, checked, and parsed from the text, rather than letting
  // parseQueryDefFile open the file again.
  return query_defs(read_file(to_path(path)), standardize,
                    std::string(delimiter), std::string(comment), name_column,
                    smarts_column);
}

std::unique_ptr<MolList> parse_mol_query_def_text(
    rust::Str text, bool standardize, rust::Str delimiter, rust::Str comment,
    std::uint32_t name_column, std::uint32_t smarts_column) {
  return query_defs(std::string(text), standardize, std::string(delimiter),
                    std::string(comment), name_column, smarts_column);
}

std::unique_ptr<RWMol> fragment_on_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    bool add_dummies, bool has_dummy_labels,
    rust::Slice<const IndexPair> dummy_labels, bool has_bond_types,
    rust::Slice<const BondType> bond_types,
    rust::Vec<std::uint32_t> &cuts_per_atom) {
  auto indices = check_fragment_args(mol, bond_indices, has_dummy_labels,
                                     dummy_labels, has_bond_types, bond_types);
  check_removals(mol);
  auto labels = to_label_pairs(dummy_labels);
  std::vector<BondType> types(bond_types.begin(), bond_types.end());
  std::vector<unsigned> cuts(mol.getNumAtoms());
  auto *result = RDKit::MolFragmenter::fragmentOnBonds(
      mol, indices, add_dummies, has_dummy_labels ? &labels : nullptr,
      has_bond_types ? &types : nullptr, &cuts);
  cuts_per_atom = to_vec(cuts);
  return to_rwmol(result);
}

// As fragmentOnSomeBondsHelper in Code/GraphMol/Wrap/MolOps.cpp, with the
// arguments checked as check_fragment_args says.
std::unique_ptr<MolList> fragment_on_some_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    std::uint32_t num_to_break, bool add_dummies, bool has_dummy_labels,
    rust::Slice<const IndexPair> dummy_labels, bool has_bond_types,
    rust::Slice<const BondType> bond_types, bool return_cuts_per_atom,
    rust::Vec<IndexList> &cuts_per_atom) {
  auto indices = check_fragment_args(mol, bond_indices, has_dummy_labels,
                                     dummy_labels, has_bond_types, bond_types);
  check_removals(mol);
  auto labels = to_label_pairs(dummy_labels);
  std::vector<BondType> types(bond_types.begin(), bond_types.end());
  std::vector<std::vector<unsigned>> cuts;
  std::vector<RDKit::ROMOL_SPTR> frags;
  auto list = std::make_unique<MolList>();
  cuts_per_atom.clear();
  // RDKit 2024.09 shifts 1 by num_to_break, undefined from 64; with more
  // bonds to break than it is given, it makes no fragments, as RDKit 2026.03
  // returns before it shifts (rdkit/rdkit#9451). It refuses more than 63
  // bonds first.
  if (indices.size() <= 63 && num_to_break > indices.size()) {
    return list;
  }
  RDKit::MolFragmenter::fragmentOnSomeBonds(
      mol, indices, frags, num_to_break, add_dummies,
      has_dummy_labels ? &labels : nullptr, has_bond_types ? &types : nullptr,
      return_cuts_per_atom ? &cuts : nullptr);
  for (const auto &frag : frags) {
    auto copy = std::make_unique<RWMol>(*frag);
    clear_matrix_caches(*copy);
    list->push(std::move(copy));
  }
  for (const auto &cut : cuts) {
    cuts_per_atom.push_back(IndexList{to_vec(cut)});
  }
  return list;
}

std::unique_ptr<RWMol> fragment_on_brics_bonds(const RWMol &mol) {
  check_removals(mol);
  std::lock_guard<std::mutex> lock(brics_mutex());
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_rwmol(RDKit::MolFragmenter::fragmentOnBRICSBonds(mol));
}

MolzipParamsData molzip_params_default() {
  RDKit::MolzipParams params;
  rust::Vec<rust::String> symbols;
  for (const auto &symbol : params.atomSymbols) {
    symbols.push_back(rust::String(symbol));
  }
  return MolzipParamsData{params.label, std::move(symbols),
                          params.enforceValenceRules,
                          params.generateCoordinates};
}

namespace {

// With generateCoordinates, molzip computes 2D coordinates.
std::unique_lock<std::mutex> lock_if_depicting(
    const RDKit::MolzipParams &params) {
  if (params.generateCoordinates) {
    return lock_depiction();
  }
  return {};
}

// With generateCoordinates, molzip of a list (MolFragmenter.cpp) depicts the
// zipped molecule with compute2DCoords, which trusts it (depiction_checks.h).
// Zipping is deterministic, so this zips once without coordinates and checks
// the result first. The caller holds lock_depiction().
//
// RDKit 2024.09 then gives each of `mols` (as RDKit orders them) the zipped
// atoms' coordinates, finding each atom in the zipped molecule without
// checking that it is there. It is not for a molecule RDKit leaves out: all
// but the first with FragmentOnBonds, and one after the first whose SMILES
// repeats an earlier one's. Nor is there a zipped molecule when RDKit cannot
// zip, which it reads regardless. So those are refused, and false says
// that RDKit returns no molecule, without the call that would read it.
template <typename Zip>
bool check_zipped_depiction(const RDKit::MolzipParams &params, Zip zip,
                            const std::vector<RDKit::ROMOL_SPTR> &mols) {
  if (!params.generateCoordinates) {
    return true;
  }
  if (params.label == RDKit::MolzipLabel::FragmentOnBonds && mols.size() > 1) {
    throw ValueErrorException(
        "generateCoordinates with FragmentOnBonds zips only the first "
        "molecule, and RDKit 2024.09 reads the others' atoms past the zipped "
        "molecule's");
  }
  std::set<std::string> seen;
  for (std::size_t i = 1; i < mols.size(); ++i) {
    if (!seen.insert(RDKit::MolToSmiles(*mols[i])).second) {
      throw ValueErrorException(
          "generateCoordinates with two molecules of the same SMILES after "
          "the first: RDKit zips one, and RDKit 2024.09 reads the other's "
          "atoms past the zipped molecule's");
    }
  }
  auto probe = params;
  probe.generateCoordinates = false;
  auto zipped = zip(probe);
  if (!zipped) {
    return false;
  }
  if (zipped->getNumAtoms() > 0) {
    check_compute_2d_coords_mol(*zipped, false);
  }
  return true;
}

}  // namespace

// molzip returns null when it cannot zip. It combines the molecules and then
// removes their dummy atoms (check_removals); a list's fragments after the
// first are written as SMILES first, which removes atoms from copies of a
// fragment with several components.
std::unique_ptr<RWMol> molzip(const RWMol &a, const RWMol &b,
                              const MolzipParamsData &data) {
  auto params = to_params(data);
  check_zip_labels(a, params, false);
  check_zip_labels(b, params, false);
  check_removals({&a, &b});
  auto lock = lock_if_depicting(params);
  return to_rwmol(RDKit::molzip(a, b, params));
}

std::unique_ptr<RWMol> molzip_single(const RWMol &a,
                                     const MolzipParamsData &data) {
  auto params = to_params(data);
  check_zip_labels(a, params, false);
  check_removals(a);
  auto lock = lock_if_depicting(params);
  return to_rwmol(RDKit::molzip(a, params));
}

// As molzipHelper in Code/GraphMol/Wrap/MolOps.cpp. With
// generateCoordinates, RDKit gives each input molecule a conformer, which is
// why Rust lends them mutably.
std::unique_ptr<RWMol> molzip_fragments(BorrowedMols &mols,
                                        const MolzipParamsData &data) {
  if (mols.mols().empty()) {
    return nullptr;
  }
  auto params = to_params(data);
  // With FragmentOnBonds, RDKit zips only the first molecule.
  for (const auto &mol : mols.mols()) {
    check_zip_labels(*mol, params, true);
    if (params.label == RDKit::MolzipLabel::FragmentOnBonds) {
      break;
    }
  }
  std::vector<const RDKit::ROMol *> parts;
  for (const auto &mol : mols.mols()) {
    parts.push_back(mol.get());
    if (params.label == RDKit::MolzipLabel::FragmentOnBonds) {
      break;
    }
  }
  check_removals(parts);
  auto lock = lock_if_depicting(params);
  if (!check_zipped_depiction(
          params,
          [&](const RDKit::MolzipParams &probe) {
            return RDKit::molzip(mols.mols(), probe);
          },
          mols.mols())) {
    return nullptr;
  }
  return to_rwmol(RDKit::molzip(mols.mols(), params));
}

// As rgroupRowZipHelper in Code/GraphMol/Wrap/MolOps.cpp; RDKit requires a
// "Core" entry.
std::unique_ptr<RWMol> molzip_row(BorrowedMols &row,
                                  const MolzipParamsData &data) {
  auto params = to_params(data);
  std::map<std::string, RDKit::ROMOL_SPTR> map;
  for (std::size_t i = 0; i < row.names().size(); ++i) {
    map[row.names()[i]] = row.mols()[i];
  }
  if (map.find("Core") == map.end()) {
    throw ValueErrorException("RGroup has no Core, cannot molzip");
  }
  // RDKit zips the core first, and with FragmentOnBonds only the core.
  for (const auto &[name, mol] : map) {
    if (params.label != RDKit::MolzipLabel::FragmentOnBonds || name == "Core") {
      check_zip_labels(*mol, params, true);
    }
  }
  std::vector<const RDKit::ROMol *> parts;
  for (const auto &[name, mol] : map) {
    if (params.label != RDKit::MolzipLabel::FragmentOnBonds || name == "Core") {
      parts.push_back(mol.get());
    }
  }
  check_removals(parts);
  auto lock = lock_if_depicting(params);
  // RDKit zips them as a list, the core first.
  std::vector<RDKit::ROMOL_SPTR> ordered{map["Core"]};
  for (const auto &[name, mol] : map) {
    if (name != "Core") {
      ordered.push_back(mol);
    }
  }
  if (!check_zipped_depiction(
          params,
          [&](const RDKit::MolzipParams &probe) {
            return RDKit::molzip(map, probe);
          },
          ordered)) {
    return nullptr;
  }
  return to_rwmol(RDKit::molzip(map, params));
}

}  // namespace rdkit_cxx
