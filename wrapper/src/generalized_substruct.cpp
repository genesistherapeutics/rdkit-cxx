// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/generalized_substruct.h"

#include <cstdint>
#include <cstring>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <GraphMol/MolBundle.h>
#include <GraphMol/MolEnumerator/LinkNode.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/types.h>
#include <boost/make_shared.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "rdkit-cxx/adjust_query.h"
#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/ffi/adjust_query.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/substruct.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/mol_enumerator.h"
#include "rdkit-cxx/substruct.h"
#include "rdkit-cxx/tautomer_query.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace bpt = boost::property_tree;
using XQ = ExtendedQueryMol;

template <class... Fs>
struct overloaded : Fs... {
  using Fs::operator()...;
};

// RDKit's Base64Decode (DataStructs/base64.cpp) sizes its buffer as three
// quarters of the text's length, writes a NUL at the buffer's last byte
// (before the buffer for text under two characters), and writes up to two
// bytes past the buffer for some unpadded lengths; it reads an uninitialized
// table entry for byte 0xFF. It counts the text to its first NUL and decodes
// the characters of the alphabet, skipping the rest.
void check_base64(const std::string &text) {
  const std::size_t in_len = std::strlen(text.c_str());
  std::size_t valid = 0;
  for (std::size_t i = 0; i < in_len; ++i) {
    const auto c = static_cast<unsigned char>(text[i]);
    if (c == 0xFF) {
      throw ValueErrorException("byte 0xFF in a base64 pickle");
    }
    if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
        (c >= '0' && c <= '9') || c == '+' || c == '/') {
      ++valid;
    }
  }
  const std::size_t out_len = 3 * in_len / 4;
  const std::size_t rest = valid % 4;
  const std::size_t written =
      3 * (valid / 4) + (rest == 2 ? 1 : (rest == 3 ? 3 : 0));
  if (out_len == 0 || written > out_len) {
    throw ValueErrorException(
        "a base64 pickle too short for its characters, which RDKit's decoder "
        "would write past its buffer");
  }
}

// RDKit's pt_to_mol (GeneralizedSubstruct/TextIO.cpp) returns SmilesToMol's
// or SmartsToMol's null for text that does not parse, which the query then
// holds and dereferences. Parsing is deterministic, so text that parses here
// parses there.
void check_mol_element(const bpt::ptree &elem) {
  auto pkl = elem.get<std::string>("pkl", "");
  if (!pkl.empty()) {
    check_base64(pkl);
    return;
  }
  auto smiles = elem.get<std::string>("smiles", "");
  if (!smiles.empty()) {
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
      throw RDKit::SmilesParse::SmilesParseException(
          "could not parse SMILES '" + smiles + "'");
    }
    return;
  }
  // A missing "smarts" is RDKit's error.
  if (auto smarts = elem.get_optional<std::string>("smarts")) {
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmartsToMol(*smarts));
    if (!mol) {
      throw RDKit::SmilesParse::SmilesParseException(
          "could not parse SMARTS '" + *smarts + "'");
    }
  }
}

template <class F>
void check_children(const bpt::ptree &pt, const char *key, F check) {
  if (auto children = pt.get_child_optional(key)) {
    for (const auto &child : *children) {
      check(child.second);
    }
  }
}

void check_tautomer_query_element(const bpt::ptree &pt) {
  check_children(pt, "tautomers", check_mol_element);
  if (auto templ = pt.get_child_optional("template")) {
    check_mol_element(*templ);
  }
}

// Checks the molecules of an ExtendedQueryMol's JSON where
// ExtendedQueryMol::initFromJSON reads them; malformed JSON and missing keys
// are left for RDKit to report.
void check_json(const std::string &text) {
  bpt::ptree pt;
  std::istringstream ss(text);
  try {
    bpt::read_json(ss, pt);
  } catch (const bpt::ptree_error &) {
    return;
  }
  auto type = pt.get_optional<unsigned char>("xqm_type");
  if (!type) {
    return;
  }
  switch (*type) {
    case XQ::XQM_MOL:
      if (auto mol = pt.get_child_optional("mol")) {
        check_mol_element(*mol);
      }
      break;
    case XQ::XQM_MOLBUNDLE:
      check_children(pt, "mols", check_mol_element);
      break;
    case XQ::XQM_TAUTOMERQUERY:
      check_tautomer_query_element(pt);
      break;
    case XQ::XQM_TAUTOMERBUNDLE:
      check_children(pt, "tautomerQueries", check_tautomer_query_element);
      break;
    default:
      break;
  }
}

// A query read from a pickle or JSON: its molecules' monomer types are
// checked as molecule pickles' are (U10), its tautomer queries as
// tautomer_query.cpp checks them, and an empty bundle, whose pattern
// fingerprint RDKit returns as null, is refused.
void check_extended_query_mol(const XQ &query) {
  std::visit(
      overloaded{
          [](const XQ::RWMol_T &mol) {
            if (!mol) {
              throw ValueErrorException("an ExtendedQueryMol has no molecule");
            }
            check_pickled_mol(*mol);
          },
          [](const XQ::MolBundle_T &bundle) {
            if (!bundle || bundle->empty()) {
              throw ValueErrorException(
                  "an ExtendedQueryMol's bundle is empty");
            }
            for (const auto &mol : bundle->getMols()) {
              if (!mol) {
                throw ValueErrorException(
                    "an ExtendedQueryMol's bundle has a null molecule");
              }
              check_pickled_mol(*mol);
            }
          },
          [](const XQ::TautomerQuery_T &query) {
            if (!query) {
              throw ValueErrorException(
                  "an ExtendedQueryMol has no tautomer query");
            }
            check_tautomer_query(*query);
          },
          [](const XQ::TautomerBundle_T &queries) {
            if (!queries || queries->empty()) {
              throw ValueErrorException(
                  "an ExtendedQueryMol's tautomer bundle is empty");
            }
            for (const auto &query : *queries) {
              if (!query) {
                throw ValueErrorException(
                    "an ExtendedQueryMol's tautomer bundle has a null query");
              }
              check_tautomer_query(*query);
            }
          },
      },
      query.xqmol);
}

// As SubstructMatch(mol, ExtendedQueryMol, params) in XQMol.cpp, with the
// tautomer queries matched by tautomer_query.cpp's matcher, holding
// control_c_mutex (control_c_lock.h).
std::vector<RDKit::MatchVectType> substruct_match(
    const RDKit::ROMol &mol, const XQ &query,
    const RDKit::SubstructMatchParameters &params) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return std::visit(
      overloaded{
          [&](const XQ::RWMol_T &q) {
            return RDKit::SubstructMatch(mol, *q, params);
          },
          [&](const XQ::MolBundle_T &q) {
            return RDKit::SubstructMatch(mol, *q, params);
          },
          [&](const XQ::TautomerQuery_T &q) {
            return tautomer_query_substruct_of(*q, mol, params, nullptr);
          },
          [&](const XQ::TautomerBundle_T &q) {
            std::vector<RDKit::MatchVectType> res;
            for (const auto &tq : *q) {
              res = tautomer_query_substruct_of(*tq, mol, params, nullptr);
              if (!res.empty()) {
                break;
              }
            }
            return res;
          },
      },
      query.xqmol);
}

// As convertMatches in Code/GraphMol/Wrap/substructmethods.h: the target atom
// for each query atom.
rust::Vec<std::uint32_t> convert_match(const RDKit::MatchVectType &match) {
  std::vector<std::uint32_t> res(match.size());
  for (const auto &[query_idx, mol_idx] : match) {
    res[query_idx] = mol_idx;
  }
  return to_vec(res);
}

RDKit::SubstructMatchParameters first_match_only(
    const SubstructMatchParametersData &data) {
  auto params = substruct_params_from_data(data);
  params.maxMatches = 1;
  return params;
}

}  // namespace

// RDKit's LinkNodeOp (MolEnumerator/LinkNode.cpp) labels each link node's
// atoms with isotopes 1000 * (index + 1), which an atom stores in 16 bits,
// and replaces the node with a reaction whose SMARTS needs those isotopes,
// distinct ends, and single bonds. It takes the reaction's first product
// without checking there is one (ASSERT_INVARIANT, compiled out), so a node
// whose labels wrap, whose two bonds share their end, or whose bonds are not
// single reads past an empty vector. So does a second node on the same atom,
// which the first node's enumeration has replaced with a chain. The specs
// themselves are parsed as the enumerator parses them, which refuses the
// malformed ones.
void check_link_nodes(const RDKit::ROMol &mol) {
  if (!mol.hasProp(RDKit::common_properties::molFileLinkNodes)) {
    return;
  }
  std::vector<unsigned int> var_atoms;
  for (const auto &node : RDKit::MolEnumerator::utils::getMolLinkNodes(mol)) {
    const auto var_atom = node.bondAtoms.at(0).first;
    for (const auto &[begin, end] : node.bondAtoms) {
      if (begin > 64 || end > 64) {
        throw ValueErrorException(
            "a link node's atoms must be among the first 65, whose isotope "
            "labels the enumerator can store");
      }
      if (mol.getBondBetweenAtoms(begin, end)->getBondType() !=
          RDKit::Bond::SINGLE) {
        throw ValueErrorException("a link node's bonds must be single bonds");
      }
    }
    if (node.bondAtoms.at(0).second == node.bondAtoms.at(1).second) {
      throw ValueErrorException("a link node's two bonds must differ");
    }
    for (auto other : var_atoms) {
      if (other == var_atom) {
        throw ValueErrorException("two link nodes share an atom");
      }
    }
    var_atoms.push_back(var_atom);
  }
}

std::unique_ptr<ExtendedQueryMol> extended_query_mol_from_binary(
    rust::Slice<const std::uint8_t> pickle) {
  auto query = std::make_unique<XQ>(to_std(pickle), false);
  check_extended_query_mol(*query);
  return query;
}

std::unique_ptr<ExtendedQueryMol> extended_query_mol_from_json(rust::Str text) {
  std::string json(text);
  check_json(json);
  auto query = std::make_unique<XQ>(json, true);
  check_extended_query_mol(*query);
  return query;
}

// RDKit's copy shares a bundle's molecules with the original, which a
// pattern fingerprint writes through (it finds their rings), so a copy that
// moves to another thread copies them.
std::unique_ptr<ExtendedQueryMol> extended_query_mol_clone(
    const ExtendedQueryMol &query) {
  return std::visit(
      overloaded{
          [](const XQ::RWMol_T &q) {
            return std::make_unique<XQ>(std::make_unique<RWMol>(*q));
          },
          [](const XQ::MolBundle_T &q) {
            auto bundle = std::make_unique<RDKit::MolBundle>();
            for (const auto &mol : q->getMols()) {
              bundle->addMol(boost::make_shared<RDKit::ROMol>(*mol));
            }
            return std::make_unique<XQ>(std::move(bundle));
          },
          // TautomerQuery's copy constructor copies the molecules.
          [](const XQ::TautomerQuery_T &q) {
            return std::make_unique<XQ>(std::make_unique<TautomerQuery>(*q));
          },
          [](const XQ::TautomerBundle_T &q) {
            auto queries =
                std::make_unique<std::vector<std::unique_ptr<TautomerQuery>>>();
            for (const auto &tq : *q) {
              queries->push_back(std::make_unique<TautomerQuery>(*tq));
            }
            return std::make_unique<XQ>(std::move(queries));
          },
      },
      query.xqmol);
}

namespace {

// toBinary and toJSON (TextIO.cpp) pickle every molecule the query holds
// (check_picklable).
void check_picklable_query(const XQ &query) {
  std::visit(overloaded{
                 [](const XQ::RWMol_T &mol) { check_picklable(*mol); },
                 [](const XQ::MolBundle_T &bundle) {
                   for (const auto &mol : bundle->getMols()) {
                     check_picklable(*mol);
                   }
                 },
                 [](const XQ::TautomerQuery_T &q) {
                   check_picklable_tautomer_query(*q);
                 },
                 [](const XQ::TautomerBundle_T &queries) {
                   for (const auto &q : *queries) {
                     check_picklable_tautomer_query(*q);
                   }
                 },
             },
             query.xqmol);
}

}  // namespace

// ExtendedQueryMol::toBinary writes a tautomer bundle's size in 16 bits and
// then every query, so a larger bundle (JSON can hold one, and so can an
// enumeration) would read back as its first size mod 2^16 queries.
rust::Vec<std::uint8_t> extended_query_mol_to_binary(
    const ExtendedQueryMol &query) {
  if (const auto *bundle = std::get_if<XQ::TautomerBundle_T>(&query.xqmol);
      bundle && *bundle &&
      (*bundle)->size() > std::numeric_limits<std::uint16_t>::max()) {
    throw ValueErrorException(
        "RDKit's binary form holds at most 65535 tautomer queries");
  }
  check_picklable_query(query);
  return to_bytes(query.toBinary());
}

std::unique_ptr<std::string> extended_query_mol_to_json(
    const ExtendedQueryMol &query) {
  check_picklable_query(query);
  return std::make_unique<std::string>(query.toJSON());
}

std::unique_ptr<BitVect> extended_query_mol_pattern_fingerprint_query(
    const ExtendedQueryMol &query, std::uint32_t fingerprint_size) {
  check_fingerprint_size(fingerprint_size);
  // Pattern fingerprints match their patterns (control_c_lock.h).
  std::unique_lock<ControlCMutex> control_c(control_c_mutex());
  auto fp = query.patternFingerprintQuery(fingerprint_size);
  control_c.unlock();
  if (!fp) {
    throw ValueErrorException("an empty query has no pattern fingerprint");
  }
  return wrap_bit_vect(std::move(fp));
}

namespace {

// As MolOps::adjustQueryProperties(mol, params), which copies the molecule
// and replaces its bonds with RWMol::replaceBond: that leaves the copy's
// stereo groups pointing at the deleted bonds (replace_bond), so the copy's
// groups are guarded.
boost::shared_ptr<RDKit::ROMol> adjusted(
    const RDKit::ROMol &mol,
    const RDKit::MolOps::AdjustQueryParameters &params) {
  auto res = boost::make_shared<RDKit::RWMol>(mol);
  StereoGroupBondGuard guard(*res);
  RDKit::MolOps::adjustQueryProperties(*res, &params);
  return res;
}

// As GeneralizedSubstruct::createExtendedQueryMol in XQMol.cpp, with the
// query properties adjusted by `adjusted` and the tautomer queries made by
// make_tautomer_query, since RDKit's adjustQueryProperties and
// TautomerQuery::fromMol leave stereo groups pointing at bonds they delete.
XQ create_extended_query(const RDKit::ROMol &mol, bool do_enumeration,
                         bool do_tautomers, bool adjust_query_properties,
                         const RDKit::MolOps::AdjustQueryParameters &params) {
  RDKit::MolBundle bundle;
  if (do_enumeration) {
    bundle = RDKit::MolEnumerator::enumerate(mol);
  }
  if (bundle.empty()) {
    // Nothing enumerated: the molecule itself.
    bundle.addMol(boost::make_shared<RDKit::ROMol>(mol));
  }
  RDKit::MolBundle adjusted_bundle;
  for (const auto &bmol : bundle.getMols()) {
    adjusted_bundle.addMol(adjust_query_properties ? adjusted(*bmol, params)
                                                   : bmol);
  }
  if (adjusted_bundle.size() == 1) {
    const auto &lmol = *adjusted_bundle.getMol(0);
    if (do_tautomers) {
      auto tq = make_tautomer_query(lmol);
      if (tq->getTautomers().size() != 1) {
        return XQ(std::move(tq));
      }
    }
    // No tautomers, just one molecule.
    return XQ(std::make_unique<RWMol>(lmol));
  }
  if (do_tautomers) {
    bool had_tautomers = false;
    auto queries =
        std::make_unique<std::vector<std::unique_ptr<TautomerQuery>>>();
    for (const auto &bmol : adjusted_bundle.getMols()) {
      auto tq = make_tautomer_query(*bmol);
      had_tautomers = had_tautomers || tq->getTautomers().size() > 1;
      queries->push_back(std::move(tq));
    }
    if (had_tautomers) {
      return XQ(std::move(queries));
    }
  }
  return XQ(std::make_unique<RDKit::MolBundle>(adjusted_bundle));
}

}  // namespace

std::unique_ptr<ExtendedQueryMol> create_extended_query_mol(
    const RWMol &mol, bool do_enumeration, bool do_tautomers,
    bool adjust_query_properties,
    const AdjustQueryParametersData &adjust_query_parameters) {
  if (do_enumeration) {
    check_link_nodes(mol);
    // createExtendedQueryMol enumerates with MolEnumerator::enumerate(mol).
    check_enumerate(mol, true, true, true);
  }
  // Enumerating removes atoms from copies (position variation, repeat
  // units), and tautomer enumeration writes SMILES, which removes atoms from
  // copies of a molecule of several fragments (check_removals); adjusting
  // removes nothing.
  if (do_enumeration || do_tautomers) {
    check_removals(mol);
  }
  // Tautomer enumeration matches its transforms (control_c_lock.h).
  std::unique_lock<ControlCMutex> control_c(control_c_mutex());
  auto query = std::make_unique<XQ>(create_extended_query(
      mol, do_enumeration, do_tautomers, adjust_query_properties,
      adjust_query_params_from_data(adjust_query_parameters)));
  check_extended_query_mol(*query);
  return query;
}

bool mol_has_substruct_match_extended(
    const RWMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &params) {
  return !substruct_match(mol, query, first_match_only(params)).empty();
}

rust::Vec<std::uint32_t> mol_get_substruct_match_extended(
    const RWMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &params) {
  auto matches = substruct_match(mol, query, first_match_only(params));
  return matches.empty() ? rust::Vec<std::uint32_t>()
                         : convert_match(matches[0]);
}

rust::Vec<IndexList> mol_get_substruct_matches_extended(
    const RWMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParametersData &params) {
  rust::Vec<IndexList> result;
  for (const auto &match :
       substruct_match(mol, query, substruct_params_from_data(params))) {
    result.push_back(IndexList{convert_match(match)});
  }
  return result;
}

std::unique_ptr<BitVect> pattern_fingerprint_target(
    const RWMol &target, std::uint32_t fingerprint_size) {
  check_fingerprint_size(fingerprint_size);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return wrap_bit_vect(RDKit::GeneralizedSubstruct::patternFingerprintTargetMol(
      target, fingerprint_size));
}

}  // namespace rdkit_cxx
