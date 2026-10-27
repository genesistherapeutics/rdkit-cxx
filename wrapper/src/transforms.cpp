// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/transforms.h"

#include <algorithm>
#include <map>
#include <mutex>
#include <optional>
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
#include <boost/make_shared.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/transforms.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace RDKit {
// MolFragmenter.cpp's molzip, which also records which atom took each removed
// dummy's place. RDKit 2025.03 exports it but does not declare it.
std::unique_ptr<ROMol> molzip(const ROMol &a, const ROMol &b,
                              const MolzipParams &params,
                              std::optional<std::map<int, int>> &mapping);
}  // namespace RDKit

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
// whose size is checked, and dummyLabels must cover the bonds too. A bond
// may be given once: fragmentOnBonds takes pointers to all of them first and
// then removes each, so a second copy is used after it is freed (Python
// crashes).
std::vector<unsigned> check_fragment_args(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    bool has_dummy_labels, rust::Slice<const IndexPair> dummy_labels,
    bool has_bond_types, rust::Slice<const BondType> bond_types) {
  if (bond_indices.empty()) {
    throw ValueErrorException("empty bond indices");
  }
  std::vector<unsigned> indices;
  std::vector<bool> seen(mol.getNumBonds());
  for (auto idx : bond_indices) {
    check_index(idx, mol.getNumBonds());
    if (seen[idx]) {
      throw ValueErrorException("bond index " + std::to_string(idx) +
                                " is repeated");
    }
    seen[idx] = true;
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

namespace {

// For a molecule atom that several core atoms match, replaceCore takes each
// such core atom's one neighbor (it checks the degree) and finds the pair
// that matches that neighbor with std::find_if, dereferencing the end when
// there is none (Python crashes). It does this, in the match's order, only
// once every index is in range.
void check_core_match(const RWMol &mol, const RWMol &core,
                      const RDKit::MatchVectType &match) {
  const int num_core = core.getNumAtoms();
  const int num_mol = mol.getNumAtoms();
  std::vector<int> times(num_mol);
  for (const auto &[core_idx, mol_idx] : match) {
    if (core_idx < 0 || core_idx >= num_core || mol_idx < 0 ||
        mol_idx >= num_mol) {
      return;  // replaceCore raises.
    }
    ++times[mol_idx];
  }
  for (const auto &[core_idx, mol_idx] : match) {
    if (times[mol_idx] < 2) {
      continue;
    }
    const auto *atom = core.getAtomWithIdx(core_idx);
    if (atom->getDegree() != 1) {
      return;  // replaceCore raises.
    }
    const int nbr = core[*core.getAtomNeighbors(atom).first]->getIdx();
    if (std::none_of(match.begin(), match.end(),
                     [nbr](const auto &pair) { return pair.first == nbr; })) {
      throw ValueErrorException(
          "core atom " + std::to_string(core_idx) + " matches molecule atom " +
          std::to_string(mol_idx) +
          ", which another core atom matches too, but its neighbor " +
          std::to_string(nbr) + " is not in the match");
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

// addRecursiveQueries copies each query molecule into its query.
void add_recursive_queries(RWMol &mol, BorrowedMols &queries,
                           rust::Str prop_name) {
  std::map<std::string, RDKit::ROMOL_SPTR> map;
  for (std::size_t i = 0; i < queries.names().size(); ++i) {
    map[queries.names()[i]] = queries.mols()[i];
  }
  RDKit::addRecursiveQueries(mol, map, std::string(prop_name));
}

std::unique_ptr<MolList> parse_mol_query_def_file(
    rust::Slice<const std::uint8_t> path, bool standardize, rust::Str delimiter,
    rust::Str comment, std::uint32_t name_column, std::uint32_t smarts_column) {
  std::map<std::string, RDKit::ROMOL_SPTR> defs;
  RDKit::parseQueryDefFile(to_path(path), defs, standardize,
                           std::string(delimiter), std::string(comment),
                           name_column, smarts_column);
  auto list = std::make_unique<MolList>();
  for (const auto &[name, mol] : defs) {
    list->push_named(name, std::make_unique<RWMol>(*mol));
  }
  return list;
}

std::unique_ptr<MolList> parse_mol_query_def_text(
    rust::Str text, bool standardize, rust::Str delimiter, rust::Str comment,
    std::uint32_t name_column, std::uint32_t smarts_column) {
  std::map<std::string, RDKit::ROMOL_SPTR> defs;
  RDKit::parseQueryDefText(std::string(text), defs, standardize,
                           std::string(delimiter), std::string(comment),
                           name_column, smarts_column);
  auto list = std::make_unique<MolList>();
  for (const auto &[name, mol] : defs) {
    list->push_named(name, std::make_unique<RWMol>(*mol));
  }
  return list;
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
// arguments checked as check_fragment_args says. RDKit starts from
// (1 << numToBreak) - 1 in 64 bits, which is undefined for 64 or more; with
// at most 63 bonds (RDKit's limit, which it checks first) and more to break
// than there are, it finds no fragments, so such a call returns none.
std::unique_ptr<MolList> fragment_on_some_bonds(
    const RWMol &mol, rust::Slice<const std::uint32_t> bond_indices,
    std::uint32_t num_to_break, bool add_dummies, bool has_dummy_labels,
    rust::Slice<const IndexPair> dummy_labels, bool has_bond_types,
    rust::Slice<const BondType> bond_types, bool return_cuts_per_atom,
    rust::Vec<IndexList> &cuts_per_atom) {
  auto indices = check_fragment_args(mol, bond_indices, has_dummy_labels,
                                     dummy_labels, has_bond_types, bond_types);
  check_removals(mol);
  cuts_per_atom.clear();
  if (indices.size() <= 63 && num_to_break > indices.size()) {
    return std::make_unique<MolList>();
  }
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

// With generateCoordinates, molzip of a list (MolFragmenter.cpp in RDKit
// 2025.03) numbers the inputs' atoms in a property, zips them, depicts the
// zipped molecule with compute2DCoords, which trusts it
// (depiction_checks.h), and then finds each input atom in it, or the atom
// that took a removed dummy's place, with std::find_if. It dereferences a
// null zipped molecule (when it cannot zip) and the end of the zipped atoms
// for an input atom it did not keep (a duplicate fragment it dropped, or any
// fragment past the first with the FragmentOnBonds label); RDKit 2026.03
// checks both. Zipping is deterministic, so this replays the zip on copies,
// without coordinates, and refuses those cases before RDKit runs. The caller
// holds lock_depiction().
void check_zip_list(const std::vector<RDKit::ROMOL_SPTR> &decomposition,
                    const RDKit::MolzipParams &params) {
  if (!params.generateCoordinates || decomposition.empty()) {
    return;
  }
  const std::string index_prop = "__zipIndex";
  std::vector<RDKit::ROMOL_SPTR> copies;
  int num_indices = 0;
  for (const auto &mol : decomposition) {
    auto copy = boost::make_shared<RDKit::ROMol>(*mol);
    for (auto atom : copy->atoms()) {
      atom->setProp(index_prop, num_indices++);
    }
    copies.push_back(copy);
  }
  std::vector<RDKit::ROMOL_SPTR> mols;
  if (params.label != RDKit::MolzipLabel::FragmentOnBonds &&
      copies.size() > 1) {
    std::vector<std::string> seen;
    for (std::size_t idx = 1; idx < copies.size(); ++idx) {
      auto smiles = RDKit::MolToSmiles(*copies[idx]);
      if (std::find(seen.begin(), seen.end(), smiles) == seen.end()) {
        mols.push_back(copies[idx]);
        seen.push_back(smiles);
      }
    }
  }
  auto combined = copies[0];
  for (const auto &mol : mols) {
    combined.reset(RDKit::combineMols(*combined, *mol));
  }
  auto probe = params;
  probe.generateCoordinates = false;
  std::optional<std::map<int, int>> mapping = std::map<int, int>();
  static const RDKit::ROMol none;
  auto zipped = RDKit::molzip(*combined, none, probe, mapping);
  if (!zipped) {
    throw ValueErrorException(
        "molzip cannot zip these molecules, and RDKit 2025.03 would then "
        "generate coordinates for no molecule");
  }
  if (zipped->getNumAtoms() == 0) {
    return;
  }
  check_compute_2d_coords_mol(*zipped, false);
  std::set<int> kept;
  for (const auto atom : zipped->atoms()) {
    int index;
    if (atom->getPropIfPresent(index_prop, index)) {
      kept.insert(index);
    }
  }
  for (int index = 0; index < num_indices; ++index) {
    auto found = mapping->find(index);
    if (!kept.count(found == mapping->end() ? index : found->second)) {
      throw ValueErrorException(
          "molzip with generateCoordinates cannot place every input atom in "
          "the zipped molecule, which RDKit 2025.03 does not check");
    }
  }
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
  check_zip_list(mols.mols(), params);
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
  // As RDKit's molzip of a row: the core, then the rest by name.
  std::vector<RDKit::ROMOL_SPTR> decomposition{map["Core"]};
  for (const auto &[name, mol] : map) {
    if (name != "Core") {
      decomposition.push_back(mol);
    }
  }
  check_zip_list(decomposition, params);
  return to_rwmol(RDKit::molzip(map, params));
}

}  // namespace rdkit_cxx
