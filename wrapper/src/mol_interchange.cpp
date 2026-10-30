// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_interchange.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/MolInterchange/details.h>
#include <GraphMol/MolPickler.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>

// As Parser.cpp: boost::json header-only, its sources compiled once, here.
#define BOOST_JSON_NO_LIB
#ifndef BOOST_CONTAINER_NO_LIB
#define BOOST_CONTAINER_NO_LIB
#endif
#include <boost/json.hpp>
#include <boost/json/src.hpp>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/mol_interchange.rs.h"

namespace rdkit_cxx {

namespace {

namespace bj = boost::json;
namespace MI = RDKit::MolInterchange;
using RDKit::MolPickler;

// Every query JSON makes compares the target's value with its own by
// subtracting in int (Queries::queryCmp), and the values they read are below
// 2^30 in magnitude, as for rdqueries (rdqueries.cpp), so query values and
// tolerances within 2^30 keep the difference inside int (CONVENTIONS.md U7).
constexpr std::int64_t kMaxQueryValue = std::int64_t{1} << 30;

[[noreturn]] void refuse(const std::string &why) {
  throw ValueErrorException("the JSON " + why);
}

// RDKit reads the JSON's strings with c_str(), so they end at a NUL.
std::string text(const bj::string &s) { return std::string(s.c_str()); }

const bj::value *member(const bj::value &v, std::string_view key) {
  const auto *obj = v.if_object();
  if (!obj) {
    return nullptr;
  }
  auto it = obj->find(key);
  return it == obj->end() ? nullptr : &it->value();
}

// The descriptions finalizeQueryFromDescription (QueryOps.cpp) gives a data
// function, for queries that compare an atom's or bond's value with theirs.
// RDKit makes each from EqualityQuery or a subclass of it, and its writers
// and matchers take a query of one of these descriptions for one
// (static_cast<ATOM_EQUALS_QUERY *> in SmartsWrite.cpp, MolFileWriter.cpp,
// QueryOps.cpp, QueryAtom.cpp, Reaction.cpp, ...).
bool is_atom_value_query(const std::string &descr) {
  static const std::set<std::string> descriptions{
      "AtomRingBondCount",
      "AtomHasRingBond",
      "AtomRingSize",
      "AtomMinRingSize",
      "AtomImplicitValence",
      "AtomTotalValence",
      "AtomAtomicNum",
      "AtomExplicitDegree",
      "AtomTotalDegree",
      "AtomHeavyAtomDegree",
      "AtomHCount",
      "AtomImplicitHCount",
      "AtomHasImplicitH",
      "AtomIsAromatic",
      "AtomIsAliphatic",
      "AtomUnsaturated",
      "AtomMass",
      "AtomIsotope",
      "AtomFormalCharge",
      "AtomNegativeFormalCharge",
      "AtomHybridization",
      "AtomInRing",
      "AtomInNRings",
      "AtomHasHeteroatomNeighbors",
      "AtomNumHeteroatomNeighbors",
      "AtomNonHydrogenDegree",
      "AtomHasAliphaticHeteroatomNeighbors",
      "AtomNumAliphaticHeteroatomNeighbors",
      "AtomType",
      "AtomNumRadicalElectrons"};
  return descriptions.count(descr) > 0;
}

bool is_bond_value_query(const std::string &descr) {
  static const std::set<std::string> descriptions{
      "BondRingSize",
      "BondMinRingSize",
      "BondOrder",
      "BondDir",
      "BondInRing",
      "BondInNRings",
      "SingleOrAromaticBond",
      "SingleOrDoubleBond",
      "DoubleOrAromaticBond",
      "SingleOrDoubleOrAromaticBond"};
  return descriptions.count(descr) > 0;
}

// EqualityQuery and its subclasses (Queries/*.h).
bool is_equality_tag(int tag) {
  return tag == MolPickler::QUERY_EQUALS || tag == MolPickler::QUERY_GREATER ||
         tag == MolPickler::QUERY_GREATEREQUAL ||
         tag == MolPickler::QUERY_LESS || tag == MolPickler::QUERY_LESSEQUAL;
}

bool is_boolean_tag(int tag) {
  return tag == MolPickler::QUERY_AND || tag == MolPickler::QUERY_OR ||
         tag == MolPickler::QUERY_XOR;
}

// Checks, before RDKit's JSON parser (MolInterchange/Parser.cpp) reads a
// document, what it trusts, following it read by read. Where RDKit would
// throw for what it reads (a value of the wrong type, a missing key), the
// check stops there and leaves the error to RDKit, which reads no further;
// what RDKit would read without a check is refused with ValueErrorException.
class DocumentChecker {
 public:
  // As DocToMols.
  explicit DocumentChecker(const bj::value &doc) {
    if (!doc.is_object()) {
      return;
    }
    if (const auto *header = member(doc, "commonchem")) {
      const auto *version = member(*header, "version");
      if (!version || !version->is_int64() ||
          version->as_int64() != MI::currentMolJSONVersion) {
        return;
      }
    } else if (const auto *header = member(doc, "rdkitjson")) {
      const auto *version = member(*header, "version");
      if (!version || !version->is_int64() ||
          static_cast<int>(version->as_int64()) > MI::currentRDKitJSONVersion ||
          static_cast<int>(version->as_int64()) < 10) {
        return;
      }
    } else {
      return;
    }
    if (const auto *defaults = member(doc, "defaults")) {
      if (!defaults->is_object()) {
        return;
      }
      const auto *atom = member(*defaults, "atom");
      if (atom && !atom->is_object()) {
        return;
      }
      d_bond_defaults = member(*defaults, "bond");
      if (d_bond_defaults && !d_bond_defaults->is_object()) {
        return;
      }
    }
    if (const auto *mols = member(doc, "molecules")) {
      if (const auto *list = mols->if_array()) {
        for (const auto &mol : *list) {
          check_mol(mol);
        }
      }
    }
  }

 private:
  // A bond's stereo, as getStringDefaultValue reads it, or nothing where it
  // throws.
  std::optional<std::string> bond_stereo(const bj::value &bond) const {
    const bj::value *stereo = member(bond, "stereo");
    if (!stereo && d_bond_defaults) {
      stereo = member(*d_bond_defaults, "stereo");
    }
    if (!stereo) {
      return std::string();
    }
    if (const auto *s = stereo->if_string()) {
      return text(*s);
    }
    return std::nullopt;
  }

  // As processMol, for a molecule and for a recursive query's.
  void check_mol(const bj::value &mol) {
    if (!mol.is_object()) {
      return;
    }
    if (const auto *name = member(mol, "name"); name && !name->is_string()) {
      return;
    }
    const auto *atoms = member(mol, "atoms");
    const auto *bonds = member(mol, "bonds");
    if (!atoms || !bonds || !atoms->is_array() || !bonds->is_array()) {
      return;
    }
    const std::size_t num_atoms = atoms->as_array().size();
    const std::size_t num_bonds = bonds->as_array().size();

    // readBond, then readBondStereo for every bond once one has stereo.
    std::vector<std::string> stereo;
    bool need_stereo_loop = false;
    for (const auto &bond : bonds->as_array()) {
      if (!bond.is_object()) {
        return;
      }
      auto s = bond_stereo(bond);
      if (!s) {
        return;
      }
      need_stereo_loop |= *s != "unspecified";
      stereo.push_back(std::move(*s));
    }
    if (need_stereo_loop) {
      for (std::size_t i = 0; i < num_bonds; ++i) {
        if (stereo[i] == "unspecified") {
          continue;
        }
        if (MI::stereoBondlookup.find(stereo[i]) ==
            MI::stereoBondlookup.end()) {
          return;
        }
        if (const auto *atoms_of =
                member(bonds->as_array()[i], "stereoAtoms")) {
          const auto *list = atoms_of->if_array();
          if (!list) {
            return;
          }
          // readBondStereo reads aids[0] and aids[1] without checking that
          // they exist.
          if (list->size() < 2) {
            refuse("gives a bond fewer than two stereo atoms");
          }
        }
      }
    }

    // readStereoGroups: the bonds the groups hold pointers to.
    std::set<unsigned> group_bonds;
    if (const auto *groups = member(mol, "stereoGroups")) {
      const auto *list = groups->if_array();
      if (!list) {
        return;
      }
      for (const auto &group : *list) {
        if (const auto *ids = member(group, "bonds")) {
          const auto *id_list = ids->if_array();
          if (!id_list) {
            return;
          }
          for (const auto &id : *id_list) {
            if (!id.is_int64()) {
              return;
            }
            group_bonds.insert(static_cast<unsigned>(id.as_int64()));
          }
        }
      }
    }

    // readSubstanceGroups. RDKit checks the atoms, bonds, and parent atoms
    // as it sets them, but not the brackets, the bonds of the cstates, or
    // the atoms of the attachment points, which the crate trusts as
    // SubstanceGroup's own setters keep them (substance_group.cpp).
    if (const auto *sgroups = member(mol, "substanceGroups")) {
      const auto *list = sgroups->if_array();
      if (!list) {
        return;
      }
      for (const auto &sgroup : *list) {
        if (!sgroup.is_object()) {
          return;
        }
        if (const auto *brackets = member(sgroup, "brackets")) {
          const auto *bracket_list = brackets->if_array();
          if (!bracket_list) {
            return;
          }
          for (const auto &bracket : *bracket_list) {
            const auto *points = bracket.if_array();
            if (!points) {
              return;
            }
            // RDKit writes each point into a Bracket, an array of three,
            // past its end for a fourth.
            if (points->size() > 3) {
              refuse(
                  "gives a substance group a bracket of more than three "
                  "points");
            }
          }
        }
        if (const auto *cstates = member(sgroup, "cstates")) {
          const auto *cstate_list = cstates->if_array();
          if (!cstate_list) {
            return;
          }
          for (const auto &cstate : *cstate_list) {
            const auto *bond = member(cstate, "bond");
            if (!bond || !bond->is_int64()) {
              return;
            }
            if (static_cast<unsigned>(bond->as_int64()) >= num_bonds) {
              refuse(
                  "gives a substance group a cstate for a bond the "
                  "molecule lacks");
            }
          }
        }
        if (const auto *points = member(sgroup, "attachPoints")) {
          const auto *point_list = points->if_array();
          if (!point_list) {
            return;
          }
          for (const auto &point : *point_list) {
            const auto *a_idx = member(point, "aIdx");
            if (!a_idx || !a_idx->is_int64()) {
              return;
            }
            if (static_cast<unsigned>(a_idx->as_int64()) >= num_atoms) {
              refuse(
                  "gives a substance group an attachment point at an "
                  "atom the molecule lacks");
            }
            if (const auto *lv_idx = member(point, "lvIdx")) {
              if (!lv_idx->is_int64()) {
                return;
              }
              // As RDKit stores it, in AttachPoint's int.
              auto lv =
                  static_cast<int>(static_cast<unsigned>(lv_idx->as_int64()));
              if (lv != -1 &&
                  (lv < 0 || static_cast<std::size_t>(lv) >= num_atoms)) {
                refuse(
                    "gives a substance group an attachment point that "
                    "leaves an atom the molecule lacks");
              }
            }
          }
        }
      }
    }

    // The extensions, in order.
    std::size_t replaced_bonds = 0;
    if (const auto *extensions = member(mol, "extensions")) {
      const auto *list = extensions->if_array();
      if (!list) {
        return;
      }
      for (const auto &extension : *list) {
        const auto *name_value = member(extension, "name");
        const auto *name = name_value ? name_value->if_string() : nullptr;
        if (!name) {
          return;
        }
        const auto *version = member(extension, "formatVersion");
        if (*name == "rdkitRepresentation") {
          if (!version || !version->is_int64()) {
            return;
          }
          if (version->as_int64() > MI::currentRDKitRepresentationVersion) {
            continue;
          }
          if (!check_representation(extension)) {
            return;
          }
        } else if (*name == "rdkitQueries") {
          if (!version || !version->is_int64()) {
            return;
          }
          if (version->as_int64() > MI::currentQueryRepresentationVersion) {
            continue;
          }
          std::size_t replaced = 0;
          if (!check_queries(extension, replaced)) {
            return;
          }
          replaced_bonds = std::max(replaced_bonds, replaced);
        }
      }
    }
    // readQueries replaces the first bonds, one per bond query
    // (RWMol::replaceBond), deleting each old bond, but replaceBond, unlike
    // replaceAtom, leaves stereo groups pointing at the bond it deletes.
    for (auto idx : group_bonds) {
      if (idx < replaced_bonds && idx < num_bonds) {
        refuse(
            "gives a bond in a stereo group a query; RDKit's replaceBond "
            "leaves the group pointing at the bond it deletes");
      }
    }
  }

  // As readRDKitRepresentation; false where RDKit throws.
  bool check_representation(const bj::value &rep) {
    if (const auto *rings = member(rep, "atomRings")) {
      const auto *list = rings->if_array();
      if (!list) {
        return false;
      }
      for (const auto &ring : *list) {
        const auto *atoms = ring.if_array();
        if (!atoms) {
          return false;
        }
        // RDKit reads an empty ring's element sz - 1, past the end.
        if (atoms->empty()) {
          refuse("gives an empty atom ring");
        }
      }
    }
    return true;
  }

  // As readQueries; counts the bond queries, which replace bonds from the
  // first; false where RDKit throws.
  bool check_queries(const bj::value &rep, std::size_t &replaced_bonds) {
    if (const auto *queries = member(rep, "atomQueries")) {
      const auto *list = queries->if_array();
      if (!list) {
        return false;
      }
      for (const auto &query : *list) {
        if (!query.is_object()) {
          return false;
        }
        if (member(query, "tag")) {
          check_atom_query(query);
        }
      }
    }
    if (const auto *queries = member(rep, "bondQueries")) {
      const auto *list = queries->if_array();
      if (!list) {
        return false;
      }
      for (const auto &query : *list) {
        if (!query.is_object()) {
          return false;
        }
        if (member(query, "tag")) {
          ++replaced_bonds;
          check_bond_query(query, false);
        }
      }
    }
    return true;
  }

  static std::optional<int> query_tag(const bj::value &query) {
    const auto *tag = member(query, "tag");
    if (!tag || !tag->is_int64()) {
      return std::nullopt;
    }
    // As RDKit reads it.
    return static_cast<int>(tag->as_int64());
  }

  static std::optional<std::string> query_descr(const bj::value &query) {
    const auto *descr = member(query, "descr");
    const auto *s = descr ? descr->if_string() : nullptr;
    if (!s) {
      return std::nullopt;
    }
    return text(*s);
  }

  // The values RDKit reads for the query's tag (readBaseQuery), each within
  // U7's bound when RDKit can read it.
  static void check_values(const bj::value &query, int tag) {
    auto bounded = [&](const char *key) {
      const auto *value = member(query, key);
      if (value && value->is_int64() &&
          (value->as_int64() > kMaxQueryValue ||
           value->as_int64() < -kMaxQueryValue)) {
        refuse(std::string("gives a query a ") + key +
               " beyond 2^30, where RDKit's int comparison can overflow");
      }
    };
    if (is_equality_tag(tag) || tag == MolPickler::QUERY_ATOMRING) {
      bounded("val");
      bounded("tol");
    } else if (tag == MolPickler::QUERY_RANGE) {
      bounded("lower");
      bounded("upper");
      bounded("tol");
    }
  }

  // finishQuery's children: two or more for AND, OR, and XOR, as RDKit's
  // writers read (SmartsWrite.cpp, MolFileWriter.cpp, QueryOps.cpp,
  // PatternFingerprints.cpp take the first two without checking), and none
  // for the queries RDKit makes without them.
  template <class F>
  static void check_children(const bj::value &query, int tag, F &&check) {
    const auto *children = member(query, "children");
    const bj::array *list = nullptr;
    if (children) {
      list = children->if_array();
      if (!list) {
        return;
      }
    }
    std::size_t n = list ? list->size() : 0;
    if (is_boolean_tag(tag) && n < 2) {
      refuse("gives an AND, OR, or XOR query fewer than two children");
    }
    if (!is_boolean_tag(tag) && n > 0) {
      refuse("gives children to a query that is not AND, OR, or XOR");
    }
    if (list) {
      for (const auto &child : *list) {
        // readQuery throws for a child that is not an object with a tag.
        if (!child.is_object() || !member(child, "tag")) {
          return;
        }
        check(child);
      }
    }
  }

  // As readQuery for atoms.
  void check_atom_query(const bj::value &query) {
    auto tag = query_tag(query);
    if (!tag) {
      return;
    }
    if (*tag == MolPickler::QUERY_RECURSIVE) {
      // RDKit reads the molecule before the description.
      const auto *mol = member(query, "subquery");
      if (!mol) {
        return;
      }
      check_mol(*mol);
      const auto *atoms = member(*mol, "atoms");
      // RecursiveStructureQuery's matcher reads the molecule's first atom.
      if (atoms && atoms->is_array() && atoms->as_array().empty()) {
        refuse("gives a recursive query a molecule without atoms");
      }
    }
    auto descr = query_descr(query);
    if (!descr) {
      return;
    }
    // finalizeQueryFromDescription strips a range query's prefix, which the
    // SMARTS parser adds for its range, less-than, and greater-than queries
    // (smarts.yy) and gives those classes.
    int prefixed = -1;
    std::string base = *descr;
    for (const auto &[prefix, prefix_tag] :
         {std::pair<std::string, int>{"range_", MolPickler::QUERY_RANGE},
          {"less_", MolPickler::QUERY_LESSEQUAL},
          {"greater_", MolPickler::QUERY_GREATEREQUAL}}) {
      if (descr->rfind(prefix, 0) == 0) {
        prefixed = prefix_tag;
        base = descr->substr(prefix.size());
        break;
      }
    }
    bool known = false;
    if (prefixed >= 0) {
      known = is_atom_value_query(base) && *tag == prefixed;
    } else if (is_atom_value_query(base)) {
      // AtomRingQuery is an EqualityQuery, which RDKit makes for
      // AtomInNRings and AtomRingBondCount (QueryOps.h, QueryOps.cpp).
      known = is_equality_tag(*tag) || *tag == MolPickler::QUERY_ATOMRING;
    } else if (base == "AtomNull") {
      known = *tag == MolPickler::QUERY_NULL;
    } else if (base == "RecursiveStructure") {
      known = *tag == MolPickler::QUERY_RECURSIVE;
    } else if (base == "AtomAnd") {
      known = *tag == MolPickler::QUERY_AND;
    } else if (base == "AtomOr") {
      known = *tag == MolPickler::QUERY_OR;
    } else if (base == "AtomXor") {
      known = *tag == MolPickler::QUERY_XOR;
    }
    // RDKit takes the query's class from its tag and how to read it from
    // its description; RDKit's code that reads queries takes the
    // description to imply the class (a RecursiveStructure query is cast to
    // RecursiveStructureQuery in SubstructMatch.cpp, an AtomAtomicNum query
    // to EqualityQuery in SmartsWrite.cpp), and without a data function for
    // its description a query dereferences null when it matches.
    if (!known) {
      refuse("gives an atom a query RDKit does not make: description '" +
             *descr + "' with tag " + std::to_string(*tag));
    }
    check_values(query, *tag);
    check_children(query, *tag,
                   [&](const bj::value &child) { check_atom_query(child); });
  }

  // As readQuery for bonds. A child of an AND or OR query is not a null
  // query: SmartsWrite.cpp casts each child that is not AND or OR to
  // EqualityQuery. XOR bond queries are not read: QueryBond.cpp's matcher
  // casts the other bond's query to EqualityQuery whatever its class.
  void check_bond_query(const bj::value &query, bool in_and_or) {
    auto tag = query_tag(query);
    auto descr = query_descr(query);
    if (!tag || !descr) {
      return;
    }
    bool known = false;
    if (is_bond_value_query(*descr)) {
      known = is_equality_tag(*tag);
    } else if (*descr == "BondNull") {
      known = *tag == MolPickler::QUERY_NULL && !in_and_or;
    } else if (*descr == "BondAnd") {
      known = *tag == MolPickler::QUERY_AND;
    } else if (*descr == "BondOr") {
      known = *tag == MolPickler::QUERY_OR;
    }
    if (!known) {
      refuse("gives a bond a query RDKit does not make: description '" +
             *descr + "' with tag " + std::to_string(*tag) +
             (in_and_or ? " in an AND or OR query" : ""));
    }
    check_values(query, *tag);
    check_children(query, *tag, [&](const bj::value &child) {
      check_bond_query(child, true);
    });
  }

  const bj::value *d_bond_defaults = nullptr;
};

MI::JSONParseParameters parse_params(const JSONParseParametersData &data) {
  MI::JSONParseParameters params;
  params.setAromaticBonds = data.set_aromatic_bonds;
  params.strictValenceCheck = data.strict_valence_check;
  params.parseProperties = data.parse_properties;
  params.parseConformers = data.parse_conformers;
  params.useHCounts = data.use_h_counts;
  return params;
}

MI::JSONWriteParameters write_params(const JSONWriteParametersData &data) {
  MI::JSONWriteParameters params;
  params.useRDKitExtensions = data.use_rdkit_extensions;
  return params;
}

}  // namespace

JSONParseParametersData json_parse_parameters_default() {
  MI::JSONParseParameters params;
  JSONParseParametersData data;
  data.set_aromatic_bonds = params.setAromaticBonds;
  data.strict_valence_check = params.strictValenceCheck;
  data.parse_properties = params.parseProperties;
  data.parse_conformers = params.parseConformers;
  data.use_h_counts = params.useHCounts;
  return data;
}

JSONWriteParametersData json_write_parameters_default() {
  MI::JSONWriteParameters params;
  JSONWriteParametersData data;
  data.use_rdkit_extensions = params.useRDKitExtensions;
  return data;
}

// RDKit 2026.09's readSubstanceGroups (MolInterchange/Parser.cpp) declares
// an AttachPoint, whose int lvIdx has no initializer, and sets lvIdx only
// when the point has the key; its writer leaves the key out for -1, so its
// own documents read back an uninitialized value. The checked document
// gets -1 there, which RDKit reads through unsigned into the int as -1.
// (rdkit-2024.09's wrapper does the same.)
namespace {
void fill_leaving_atoms(bj::value &doc) {
  if (!doc.is_object()) {
    return;
  }
  auto *mols = doc.as_object().if_contains("molecules");
  if (!mols || !mols->is_array()) {
    return;
  }
  for (auto &mol : mols->as_array()) {
    auto *sgroups = mol.is_object()
                        ? mol.as_object().if_contains("substanceGroups")
                        : nullptr;
    if (!sgroups || !sgroups->is_array()) {
      continue;
    }
    for (auto &sgroup : sgroups->as_array()) {
      auto *points = sgroup.is_object()
                         ? sgroup.as_object().if_contains("attachPoints")
                         : nullptr;
      if (!points || !points->is_array()) {
        continue;
      }
      for (auto &point : points->as_array()) {
        if (point.is_object() && !point.as_object().contains("lvIdx")) {
          point.as_object()["lvIdx"] = std::int64_t{-1};
        }
      }
    }
  }
}
}  // namespace

// RDKit parses the document it is given with boost::json's default options
// (JSONDataToMols), as this does. It is handed the checked document,
// serialized again, so it parses exactly what was checked (U2), whatever its
// boost::json makes of the original text.
std::unique_ptr<MolList> json_to_mols(rust::Str json_block,
                                      const JSONParseParametersData &params) {
  auto p = parse_params(params);
  std::vector<boost::shared_ptr<RDKit::ROMol>> mols;
  try {
    bj::value doc =
        bj::parse(std::string_view(json_block.data(), json_block.size()));
    DocumentChecker check(doc);
    fill_leaving_atoms(doc);
    mols = MI::JSONDataToMols(bj::serialize(doc), p);
  } catch (const boost::system::system_error &e) {
    // boost::json's errors, for text that is not JSON and for values RDKit
    // reads as another type, are parse errors here (E4).
    throw RDKit::FileParseException(e.what());
  }
  auto list = std::make_unique<MolList>();
  for (const auto &mol : mols) {
    list->push(std::make_unique<RWMol>(*mol));
  }
  return list;
}

std::unique_ptr<std::string> mol_to_json(
    const RWMol &mol, const JSONWriteParametersData &params) {
  return std::make_unique<std::string>(
      MI::MolToJSONData(mol, write_params(params)));
}

std::unique_ptr<std::string> mols_to_json(
    BorrowedMols &mols, const JSONWriteParametersData &params) {
  std::vector<const RDKit::ROMol *> list;
  for (const auto &mol : mols.mols()) {
    list.push_back(mol.get());
  }
  return std::make_unique<std::string>(
      MI::MolsToJSONData(list, write_params(params)));
}

}  // namespace rdkit_cxx
