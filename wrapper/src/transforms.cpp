// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/transforms.h"

#include <algorithm>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/StreamOps.h>
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
// whose size is checked, and dummyLabels must cover the bonds too. RDKit
// 2024.03's fragmentOnBonds (MolFragmenter.cpp) removes a repeated bond
// twice, the second time through the pointer it has deleted, and
// fragmentOnSomeBonds passes it such a pair; 2026.03 refuses repeats with
// this message.
std::vector<unsigned> check_fragment_args(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    bool has_dummy_labels, rust::Slice<const IndexPair> dummy_labels,
    bool has_bond_types, rust::Slice<const BondType> bond_types) {
  if (bond_indices.empty()) {
    throw ValueErrorException("empty bond indices");
  }
  std::vector<unsigned> indices;
  std::unordered_set<unsigned> seen;
  for (auto idx : bond_indices) {
    check_index(idx, mol.getNumBonds());
    if (!seen.insert(idx).second) {
      throw ValueErrorException("bondIndices contains duplicates");
    }
    indices.push_back(idx);
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
// refused there.
//
// RDKit 2024.03's molzip (MolFragmenter.cpp) has no linker handling, which
// 2026.03 added (MolZip.cpp's isLinker and linkerBonds). In a linker
// ([*:1][*:2], two labeled atoms bonded to each other) each labeled atom is
// the atom another label bonds, which the batch edit deletes as a used
// dummy; restore_chirality then reads it (restore, getOwningMol) after it is
// freed. So no call may hold one, for any label: with FragmentOnBonds every
// dummy is labeled, so no dummy may be bonded to a dummy. (The atom a
// FragmentOnBonds isotope names, which RDKit bonds, can then be a dummy only
// if RDKit finds no dummy beside it naming the other atom back, and gives up
// before it deletes anything.) `mols` are the molecules RDKit combines.
void check_zip_labels(const std::vector<const RDKit::ROMol *> &mols,
                      const RDKit::MolzipParams &params, bool zipping_list) {
  const bool fob = params.label == RDKit::MolzipLabel::FragmentOnBonds;
  for (const auto *mol : mols) {
    for (const auto *atom : mol->atoms()) {
      if (!has_zip_label(atom, params)) {
        continue;
      }
      if ((fob || zipping_list) && atom->getDegree() != 1) {
        throw ValueErrorException("molzip: labeled atom " +
                                  std::to_string(atom->getIdx()) +
                                  " must have exactly one neighbor");
      }
      for (const auto *nbr : mol->atomNeighbors(atom)) {
        if (has_zip_label(nbr, params)) {
          throw ValueErrorException(
              "molzip: RDKit 2024.03 cannot zip a linker (labeled atoms " +
              std::to_string(atom->getIdx()) + " and " +
              std::to_string(nbr->getIdx()) + " bonded to each other)");
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

namespace {

// RDKit 2024.03's replaceCore (ChemTransforms.cpp), for each pair whose
// molecule atom is matched more than once, takes the core atom's one
// neighbor and dereferences the find_if for the pair that matches it,
// end() when none does (2026.03 checks). replaceCore first range-checks
// every pair, throwing for one out of range, and refuses a core atom without
// exactly one neighbor (CHECK_INVARIANT), but only once it reaches it, so
// every pair is checked here.
void check_core_match(const RWMol &mol, const RWMol &core,
                      const RDKit::MatchVectType &match) {
  const auto num_core = static_cast<int>(core.getNumAtoms());
  const auto num_mol = static_cast<int>(mol.getNumAtoms());
  std::map<int, unsigned> times_matched;
  std::unordered_set<int> core_atoms;
  for (const auto &[core_idx, mol_idx] : match) {
    if (core_idx < 0 || core_idx >= num_core || mol_idx < 0 ||
        mol_idx >= num_mol) {
      return;
    }
    ++times_matched[mol_idx];
    core_atoms.insert(core_idx);
  }
  for (const auto &[core_idx, mol_idx] : match) {
    const auto *core_atom = core.getAtomWithIdx(core_idx);
    if (times_matched[mol_idx] < 2 || core_atom->getDegree() != 1) {
      continue;
    }
    const auto *nbr = *core.atomNeighbors(core_atom).begin();
    if (!core_atoms.count(static_cast<int>(nbr->getIdx()))) {
      throw ValueErrorException(
          "molecule atom " + std::to_string(mol_idx) +
          " is matched more than once, and the neighbor " +
          std::to_string(nbr->getIdx()) + " of core atom " +
          std::to_string(core_idx) +
          " is not in the match, which RDKit 2024.03 then reads past the "
          "match's end for");
    }
  }
}

}  // namespace

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
  check_core_match(mol, core, match);
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

// addRecursiveQueries copies each query molecule into its query, and parses
// an atom's property value that names no query as SMARTS
// (check_query_props).
void add_recursive_queries(RWMol &mol, BorrowedMols &queries,
                           rust::Str prop_name) {
  std::map<std::string, RDKit::ROMOL_SPTR> map;
  for (std::size_t i = 0; i < queries.names().size(); ++i) {
    map[queries.names()[i]] = queries.mols()[i];
  }
  const std::string prop(prop_name);
  check_query_props(mol, prop);
  RDKit::addRecursiveQueries(mol, map, prop);
}

namespace {

// RDKit 2024.03's parseQueryDefFile (ChemTransforms.cpp) parses each line's
// SMARTS column with SmartsToMol and CXSMILES on, whose atomProp block reads
// past the end of the text (file_checks.h). This reads the lines as it does
// (getLine, the comment test, the tokenizer over `delimiter`, and the trim)
// and checks each SMARTS it would parse.
void check_query_defs(const std::string &text, const std::string &delimiter,
                      const std::string &comment, unsigned name_column,
                      unsigned smarts_column) {
  std::istringstream in(text);
  boost::char_separator<char> sep(delimiter.c_str());
  while (!in.eof() && !in.fail()) {
    std::string line = RDKit::getLine(in);
    if (line.empty() || line.find(comment) == 0) {
      continue;
    }
    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
    std::string name;
    std::string smarts;
    unsigned pos = 0;
    for (auto token = tokens.begin(); token != tokens.end(); ++token, ++pos) {
      if (pos == name_column) {
        name = *token;
      } else if (pos == smarts_column) {
        smarts = *token;
      }
    }
    boost::trim_if(name, boost::is_any_of(" \t"));
    boost::trim_if(smarts, boost::is_any_of(" \t"));
    if (!name.empty() && !smarts.empty()) {
      check_smiles_cx(smarts, true, true);
    }
  }
}

std::unique_ptr<MolList> query_defs(const std::string &text, bool standardize,
                                    const std::string &delimiter,
                                    const std::string &comment,
                                    std::uint32_t name_column,
                                    std::uint32_t smarts_column) {
  check_query_defs(text, delimiter, comment, name_column, smarts_column);
  std::map<std::string, RDKit::ROMOL_SPTR> defs;
  RDKit::parseQueryDefText(text, defs, standardize, delimiter, comment,
                           name_column, smarts_column);
  auto list = std::make_unique<MolList>();
  for (const auto &[name, mol] : defs) {
    list->push_named(name, std::make_unique<RWMol>(*mol));
  }
  return list;
}

}  // namespace

// The file is read once (as parseQueryDefFile reads it, a stream of its
// text) and checked; RDKit parses the text checked.
std::unique_ptr<MolList> parse_mol_query_def_file(
    rust::Slice<const std::uint8_t> path, bool standardize, rust::Str delimiter,
    rust::Str comment, std::uint32_t name_column, std::uint32_t smarts_column) {
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
// arguments checked as check_fragment_args says. RDKit 2024.03's
// fragmentOnSomeBonds (MolFragmenter.cpp) shifts 1 by num_to_break, which is
// undefined from 64 on and, as x86 masks it, fragments one bond 64 times; so
// this returns no fragments when there are not num_to_break bonds, as
// 2026.03 does, and 2024.03 does for counts below 64 (and for 0).
std::unique_ptr<MolList> fragment_on_some_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    std::uint32_t num_to_break, bool add_dummies, bool has_dummy_labels,
    rust::Slice<const IndexPair> dummy_labels, bool has_bond_types,
    rust::Slice<const BondType> bond_types, bool return_cuts_per_atom,
    rust::Vec<IndexList> &cuts_per_atom) {
  auto indices = check_fragment_args(mol, bond_indices, has_dummy_labels,
                                     dummy_labels, has_bond_types, bond_types);
  if (indices.size() > 63) {
    throw ValueErrorException("currently can only fragment on up to 63 bonds");
  }
  cuts_per_atom.clear();
  if (num_to_break == 0 || num_to_break > indices.size()) {
    return std::make_unique<MolList>();
  }
  check_removals(mol);
  auto labels = to_label_pairs(dummy_labels);
  std::vector<BondType> types(bond_types.begin(), bond_types.end());
  std::vector<std::vector<unsigned>> cuts;
  std::vector<RDKit::ROMOL_SPTR> frags;
  RDKit::MolFragmenter::fragmentOnSomeBonds(
      mol, indices, frags, num_to_break, add_dummies,
      has_dummy_labels ? &labels : nullptr, has_bond_types ? &types : nullptr,
      return_cuts_per_atom ? &cuts : nullptr);
  auto list = std::make_unique<MolList>();
  for (const auto &frag : frags) {
    auto copy = std::make_unique<RWMol>(*frag);
    clear_matrix_caches(*copy);
    list->push(std::move(copy));
  }
  cuts_per_atom.clear();
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

// With generateCoordinates, molzip of a list (RDKit 2024.03's
// MolFragmenter.cpp) depicts the zipped molecule with compute2DCoords,
// which trusts it (depiction_checks.h), and then dereferences it unchecked,
// null when molzip cannot zip (2026.03 checks). It then finds, for each atom
// of every input, the zipped atom it became and dereferences the find_if
// unchecked (2026.03 checks): end() for the atoms of the molecules it leaves
// out of the zip, which are every one after the first with FragmentOnBonds,
// and otherwise those after the first whose SMILES repeats an earlier one's
// (RDKit's ring dedup), so such lists are refused. Zipping is
// deterministic, so this zips once without coordinates and checks the
// result first; false when RDKit cannot zip, for the caller to return null
// without the call that would dereference it. `mols` are the inputs in
// RDKit's order. The caller holds lock_depiction().
template <typename Zip>
bool check_zipped_depiction(const std::vector<RDKit::ROMOL_SPTR> &mols,
                            const RDKit::MolzipParams &params, Zip zip) {
  if (!params.generateCoordinates) {
    return true;
  }
  if (params.label == RDKit::MolzipLabel::FragmentOnBonds) {
    if (mols.size() > 1) {
      throw ValueErrorException(
          "molzip: with generate_coordinates and the FragmentOnBonds label, "
          "RDKit 2024.03 reads past the zipped atoms for every molecule "
          "after the first");
    }
  } else {
    std::vector<std::string> seen;
    for (std::size_t i = 1; i < mols.size(); ++i) {
      auto smiles = RDKit::MolToSmiles(*mols[i]);
      if (std::find(seen.begin(), seen.end(), smiles) != seen.end()) {
        throw ValueErrorException(
            "molzip: with generate_coordinates, RDKit 2024.03 reads past the "
            "zipped atoms for molecule " +
            std::to_string(i) +
            ", which it leaves out as its SMILES repeats an earlier one's");
      }
      seen.push_back(std::move(smiles));
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
  check_zip_labels({&a, &b}, params, false);
  check_removals({&a, &b});
  auto lock = lock_if_depicting(params);
  return to_rwmol(RDKit::molzip(a, b, params));
}

std::unique_ptr<RWMol> molzip_single(const RWMol &a,
                                     const MolzipParamsData &data) {
  auto params = to_params(data);
  check_zip_labels({&a}, params, false);
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
  std::vector<const RDKit::ROMol *> parts;
  for (const auto &mol : mols.mols()) {
    parts.push_back(mol.get());
    if (params.label == RDKit::MolzipLabel::FragmentOnBonds) {
      break;
    }
  }
  check_zip_labels(parts, params, true);
  check_removals(parts);
  auto lock = lock_if_depicting(params);
  if (!check_zipped_depiction(mols.mols(), params,
                              [&](const RDKit::MolzipParams &probe) {
                                return RDKit::molzip(mols.mols(), probe);
                              })) {
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
  // RDKit zips the core first and then the others in the map's order, and
  // with FragmentOnBonds only the core.
  std::vector<RDKit::ROMOL_SPTR> ordered{map["Core"]};
  for (const auto &[name, mol] : map) {
    if (name != "Core") {
      ordered.push_back(mol);
    }
  }
  std::vector<const RDKit::ROMol *> parts;
  for (const auto &mol : ordered) {
    parts.push_back(mol.get());
    if (params.label == RDKit::MolzipLabel::FragmentOnBonds) {
      break;
    }
  }
  check_zip_labels(parts, params, true);
  check_removals(parts);
  auto lock = lock_if_depicting(params);
  if (!check_zipped_depiction(ordered, params,
                              [&](const RDKit::MolzipParams &probe) {
                                return RDKit::molzip(map, probe);
                              })) {
    return nullptr;
  }
  return to_rwmol(RDKit::molzip(map, params));
}

}  // namespace rdkit_cxx
