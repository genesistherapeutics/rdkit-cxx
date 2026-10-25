// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/tautomer_query.h"

#include <limits>
#include <map>
#include <set>
#include <span>
#include <string>
#include <utility>

#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/MolStandardize/TautomerCatalog/TautomerCatalogParams.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <RDGeneral/Exceptions.h>
#include <boost/dynamic_bitset.hpp>

#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/substruct.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/substruct.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace MolStandardize = RDKit::MolStandardize;

// RDKit applies a transform's bond types and charges by position
// (TautomerEnumerator::enumerate in Tautomer.cpp), one per bond and one per
// matched atom, without checking that the transform file listed that many:
// a shorter list is read past its end.
void check_transforms(const MolStandardize::TautomerCatalogParams &params) {
  for (const auto &transform : params.getTransforms()) {
    const auto &mol = *transform.Mol;
    std::string name;
    mol.getPropIfPresent(RDKit::common_properties::_Name, name);
    if (!transform.BondTypes.empty() &&
        transform.BondTypes.size() < mol.getNumBonds()) {
      throw ValueErrorException("tautomer transform '" + name +
                                "' lists fewer bond types than its SMARTS "
                                "has bonds");
    }
    if (!transform.Charges.empty() &&
        transform.Charges.size() < mol.getNumAtoms()) {
      throw ValueErrorException("tautomer transform '" + name +
                                "' lists fewer charges than its SMARTS has "
                                "atoms");
    }
  }
}

// As TautomerQuery::fromMol in TautomerQuery.cpp, which reads the transform
// file itself; here the transforms are read once and checked before use.
std::unique_ptr<TautomerQuery> from_mol(
    const RDKit::ROMol &query, MolStandardize::TautomerCatalogParams &params) {
  MolStandardize::TautomerEnumerator enumerator(
      new MolStandardize::TautomerCatalog(&params));
  const auto res = enumerator.enumerate(query);

  std::vector<size_t> modified_atoms;
  modified_atoms.reserve(res.modifiedAtoms().count());
  for (size_t i = 0; i < query.getNumAtoms(); i++) {
    if (res.modifiedAtoms().test(i)) {
      modified_atoms.push_back(i);
    }
  }
  std::vector<size_t> modified_bonds;
  modified_bonds.reserve(res.modifiedBonds().count());
  for (size_t i = 0; i < query.getNumBonds(); i++) {
    if (res.modifiedBonds().test(i)) {
      modified_bonds.push_back(i);
    }
  }

  auto templ = std::make_unique<RDKit::RWMol>(query);
  for (auto idx : modified_atoms) {
    const auto *atom = query.getAtomWithIdx(idx);
    RDKit::QueryAtom query_atom(atom->getAtomicNum());
    // Forward original queries.
    if (atom->hasQuery()) {
      query_atom.setQuery(
          static_cast<const RDKit::QueryAtom *>(atom)->getQuery()->copy());
    }
    const bool update_label = false;
    const bool preserve_props = true;
    templ->replaceAtom(idx, &query_atom, update_label, preserve_props);
  }
  for (auto idx : modified_bonds) {
    RDKit::QueryBond query_bond;
    query_bond.setQuery(RDKit::makeSingleOrDoubleOrAromaticBondQuery());
    templ->replaceBond(idx, &query_bond, true);
  }

  // As fromMol: the enumerator's quick copies lack the molecule's
  // properties, which CXSMILES output (link nodes) reads.
  auto tautomers = res.tautomers();
  for (auto &tautomer : tautomers) {
    tautomer->updateProps(query);
  }
  return std::make_unique<TautomerQuery>(
      std::move(tautomers), static_cast<RDKit::ROMol *>(templ.release()),
      std::move(modified_atoms), std::move(modified_bonds));
}

// As TautomerQuery::matchTautomer. The query is checked
// (check_tautomer_query), and the match is a full match of the template, so
// its indices are in range; the target bond is checked all the same.
bool match_tautomer(const std::vector<size_t> &atoms,
                    const std::vector<size_t> &bonds, const RDKit::ROMol &mol,
                    const RDKit::ROMol &tautomer,
                    std::span<const unsigned int> match,
                    const RDKit::SubstructMatchParameters &params) {
  for (auto idx : atoms) {
    if (!RDKit::atomCompat(tautomer.getAtomWithIdx(idx),
                           mol.getAtomWithIdx(match[idx]), params)) {
      return false;
    }
  }
  for (auto idx : bonds) {
    const auto *query_bond = tautomer.getBondWithIdx(idx);
    const auto *target_bond =
        mol.getBondBetweenAtoms(match[query_bond->getBeginAtomIdx()],
                                match[query_bond->getEndAtomIdx()]);
    if (!target_bond || !RDKit::bondCompat(query_bond, target_bond, params)) {
      return false;
    }
  }
  return true;
}

// As removeTautomerDuplicates in TautomerQuery.cpp: keeps the first match of
// each set of target atoms, with its tautomer.
void remove_duplicates(std::vector<RDKit::MatchVectType> &matches,
                       std::vector<RDKit::ROMOL_SPTR> *matching_tautomers,
                       unsigned int num_atoms) {
  std::set<boost::dynamic_bitset<>> seen;
  std::vector<RDKit::MatchVectType> res;
  std::vector<RDKit::ROMOL_SPTR> res_tautomers;
  for (size_t i = 0; i < matches.size(); ++i) {
    boost::dynamic_bitset<> val(num_atoms);
    for (const auto &pair : matches[i]) {
      val.set(pair.second);
    }
    if (seen.insert(std::move(val)).second) {
      res.push_back(std::move(matches[i]));
      if (matching_tautomers) {
        res_tautomers.push_back((*matching_tautomers)[i]);
      }
    }
  }
  matches = std::move(res);
  if (matching_tautomers) {
    *matching_tautomers = std::move(res_tautomers);
  }
}

RDKit::SubstructMatchParameters params_with_max_matches_one(
    const SubstructMatchParametersData &data) {
  auto params = substruct_params_from_data(data);
  params.maxMatches = 1;
  return params;
}

// As convertMatch in nbWrap/substructmethods.h: the target atom for each
// query atom.
rust::Vec<std::uint32_t> convert_match(const RDKit::MatchVectType &match) {
  std::vector<std::uint32_t> res(match.size());
  for (const auto &[query_idx, mol_idx] : match) {
    res[query_idx] = mol_idx;
  }
  return to_vec(res);
}

std::unique_ptr<RWMol> copy_mol(const RDKit::ROMol &mol) {
  auto result = std::make_unique<RWMol>(mol);
  clear_matrix_caches(*result);
  return result;
}

rust::Vec<std::uint32_t> to_indices(const std::vector<size_t> &indices) {
  rust::Vec<std::uint32_t> res;
  res.reserve(indices.size());
  for (auto idx : indices) {
    res.push_back(static_cast<std::uint32_t>(idx));
  }
  return res;
}

}  // namespace

// The matcher reads the template's match at each modified atom and at the
// ends of each tautomer's modified bonds, so every tautomer must have the
// template's atoms and every modified bond; RDKit trusts a pickle's or JSON's
// indices.
void check_tautomer_query(const TautomerQuery &query) {
  const auto &templ = query.getTemplateMolecule();
  check_monomer_types(templ);
  const auto num_atoms = templ.getNumAtoms();
  const auto tautomers = query.getTautomers();
  for (const auto &tautomer : tautomers) {
    if (!tautomer) {
      throw ValueErrorException("a TautomerQuery has a null tautomer");
    }
    check_monomer_types(*tautomer);
    if (tautomer->getNumAtoms() != num_atoms) {
      throw ValueErrorException(
          "a tautomer and its query's template have different numbers of "
          "atoms");
    }
  }
  for (auto idx : query.getModifiedAtoms()) {
    if (idx >= num_atoms) {
      throw ValueErrorException(
          "a TautomerQuery's modified atom is out of range");
    }
  }
  for (auto idx : query.getModifiedBonds()) {
    for (const auto &tautomer : tautomers) {
      if (idx >= tautomer->getNumBonds()) {
        throw ValueErrorException(
            "a TautomerQuery's modified bond is out of range");
      }
    }
  }
}

// As TautomerQuery::substructOf, whose matcher runs as the template match's
// extraFinalCheck. RDKit's RecursiveMatcher hands that check each recursive
// SMARTS's own matches as well, which RDKit's matcher indexes by the
// template's atoms, past their end. This one passes a match of another size
// unchecked, as a recursive query without the check would, and pairs each
// final match with the tautomer its check accepted.
std::vector<RDKit::MatchVectType> tautomer_query_substruct_of(
    const TautomerQuery &query, const RDKit::ROMol &mol,
    const RDKit::SubstructMatchParameters &params,
    std::vector<RDKit::ROMOL_SPTR> *matching_tautomers) {
  if (matching_tautomers) {
    matching_tautomers->clear();
  }
  const auto &templ = query.getTemplateMolecule();
  const auto num_atoms = templ.getNumAtoms();
  const auto tautomers = query.getTautomers();
  const auto atoms = query.getModifiedAtoms();
  const auto bonds = query.getModifiedBonds();
  std::map<std::vector<unsigned int>, RDKit::ROMOL_SPTR> accepted;

  RDKit::SubstructMatchParameters template_params(params);
  // RDKit checks every mapping of the template to the target.
  template_params.uniquify = false;
  template_params.extraFinalCheck = [&](const RDKit::ROMol &target,
                                        std::span<const unsigned int> match) {
    if (match.size() != num_atoms) {
      return true;
    }
    for (const auto &tautomer : tautomers) {
      if (match_tautomer(atoms, bonds, target, *tautomer, match, params)) {
        if (matching_tautomers) {
          accepted.emplace(
              std::vector<unsigned int>(match.begin(), match.end()), tautomer);
        }
        return true;
      }
    }
    return false;
  };
  auto matches = RDKit::SubstructMatch(mol, templ, template_params);
  if (matching_tautomers) {
    for (const auto &match : matches) {
      std::vector<unsigned int> key(num_atoms);
      for (const auto &[query_idx, mol_idx] : match) {
        key[query_idx] = mol_idx;
      }
      matching_tautomers->push_back(accepted.at(key));
    }
  }
  if (params.uniquify && matches.size() > 1) {
    remove_duplicates(matches, matching_tautomers, mol.getNumAtoms());
  }
  return matches;
}

void check_fingerprint_size(std::uint32_t fingerprint_size) {
  if (fingerprint_size >
      static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
    throw ValueErrorException("bit vectors hold at most 2^31 - 1 bits");
  }
}

// RDKit's default constructor leaves the template null, which every method
// dereferences; this query's template is an empty molecule.
std::unique_ptr<TautomerQuery> tautomer_query_new_empty() {
  return std::make_unique<TautomerQuery>(
      std::vector<RDKit::ROMOL_SPTR>(), new RDKit::ROMol(),
      std::vector<size_t>(), std::vector<size_t>());
}

std::unique_ptr<TautomerQuery> tautomer_query_clone(
    const TautomerQuery &query) {
  return std::make_unique<TautomerQuery>(query);
}

std::unique_ptr<TautomerQuery> tautomer_query_from_mol(
    const RWMol &mol, rust::Slice<const std::uint8_t> tautomer_transform_file) {
  // An empty path reads RDKit's default transforms.
  MolStandardize::TautomerCatalogParams params(
      to_path(tautomer_transform_file));
  check_transforms(params);
  return from_mol(mol, params);
}

std::unique_ptr<TautomerQuery> tautomer_query_from_pickle(
    rust::Slice<const std::uint8_t> pickle) {
  auto query = std::make_unique<TautomerQuery>(to_std(pickle));
  check_tautomer_query(*query);
  return query;
}

rust::Vec<std::uint8_t> tautomer_query_to_binary(const TautomerQuery &query) {
  return to_bytes(query.serialize());
}

bool tautomer_query_is_substruct_of(
    const TautomerQuery &query, const RWMol &target,
    const SubstructMatchParametersData &params) {
  return !tautomer_query_substruct_of(
              query, target, params_with_max_matches_one(params), nullptr)
              .empty();
}

rust::Vec<std::uint32_t> tautomer_query_get_substruct_match(
    const TautomerQuery &query, const RWMol &target,
    const SubstructMatchParametersData &params) {
  auto matches = tautomer_query_substruct_of(
      query, target, params_with_max_matches_one(params), nullptr);
  return matches.empty() ? rust::Vec<std::uint32_t>()
                         : convert_match(matches[0]);
}

rust::Vec<IndexList> tautomer_query_get_substruct_matches(
    const TautomerQuery &query, const RWMol &target,
    const SubstructMatchParametersData &params) {
  rust::Vec<IndexList> result;
  for (const auto &match : tautomer_query_substruct_of(
           query, target, substruct_params_from_data(params), nullptr)) {
    result.push_back(IndexList{convert_match(match)});
  }
  return result;
}

std::unique_ptr<MolList> tautomer_query_get_substruct_matches_with_tautomers(
    const TautomerQuery &query, const RWMol &target,
    const SubstructMatchParametersData &params, rust::Vec<IndexList> &matches) {
  std::vector<RDKit::ROMOL_SPTR> tautomers;
  auto res = tautomer_query_substruct_of(
      query, target, substruct_params_from_data(params), &tautomers);
  auto list = std::make_unique<MolList>();
  for (size_t i = 0; i < res.size(); ++i) {
    matches.push_back(IndexList{convert_match(res[i])});
    list->push(copy_mol(*tautomers.at(i)));
  }
  return list;
}

std::unique_ptr<BitVect> tautomer_query_pattern_fingerprint_template(
    const TautomerQuery &query, std::uint32_t fingerprint_size) {
  check_fingerprint_size(fingerprint_size);
  return wrap_bit_vect(std::unique_ptr<ExplicitBitVect>(
      query.patternFingerprintTemplate(fingerprint_size)));
}

std::unique_ptr<RWMol> tautomer_query_get_template_molecule(
    const TautomerQuery &query) {
  return copy_mol(query.getTemplateMolecule());
}

// Every index is below the template's number of atoms or bonds: fromMol's by
// construction, a pickle's by check_tautomer_query.
rust::Vec<std::uint32_t> tautomer_query_get_modified_atoms(
    const TautomerQuery &query) {
  return to_indices(query.getModifiedAtoms());
}

rust::Vec<std::uint32_t> tautomer_query_get_modified_bonds(
    const TautomerQuery &query) {
  return to_indices(query.getModifiedBonds());
}

std::unique_ptr<MolList> tautomer_query_get_tautomers(
    const TautomerQuery &query) {
  auto list = std::make_unique<MolList>();
  for (const auto &tautomer : query.getTautomers()) {
    list->push(copy_mol(*tautomer));
  }
  return list;
}

std::unique_ptr<BitVect> pattern_fingerprint_tautomer_target(
    const RWMol &target, std::uint32_t fingerprint_size) {
  check_fingerprint_size(fingerprint_size);
  return wrap_bit_vect(std::unique_ptr<ExplicitBitVect>(
      TautomerQuery::patternFingerprintTarget(target, fingerprint_size)));
}

bool tautomer_query_can_serialize() {
  return RDKit::TautomerQueryCanSerialize();
}

}  // namespace rdkit_cxx
