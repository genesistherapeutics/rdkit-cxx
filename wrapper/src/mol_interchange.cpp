// Copyright (c) 2026 Genesis Molecular AI, Inc.
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
#include <tuple>
#include <utility>
#include <vector>

#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/MolInterchange/details.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>
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

// How RapidJSON, which RDKit 2025.03 parses JSON with (JSONDataToMols), types
// the value boost::json serializes: an integer as an integer (an int if it
// fits in 32 bits), and a double with an exponent, so as a double.
bool rj_int(const bj::value *v) {
  return v && v->is_int64() && v->as_int64() >= INT_MIN &&
         v->as_int64() <= INT_MAX;
}
// What RapidJSON's GetInt and GetUint read without checking the type: any
// integer of 32 bits, signed or not, read through its union as the other
// (RDKit writes an attachment point's absent lvIdx as -1 and reads it with
// GetUint).
bool rj_int32(const bj::value *v) {
  return v && v->is_int64() && v->as_int64() >= INT_MIN &&
         v->as_int64() <= static_cast<std::int64_t>(UINT_MAX);
}
bool rj_number(const bj::value *v) { return v && v->is_number(); }
bool rj_double(const bj::value *v) { return v && v->is_double(); }
bool rj_string(const bj::value *v) { return v && v->is_string(); }
bool rj_bool(const bj::value *v) { return v && v->is_bool(); }
bool rj_array(const bj::value *v) { return v && v->is_array(); }
bool rj_object(const bj::value *v) { return v && v->is_object(); }

// RapidJSON's accessors do not check the type of what they read (its
// asserts are compiled out), so a value of another type is read as the one
// RDKit expects: an int from a string's pointer, an array from an object's
// members, an element past an array's end. Such a document is malformed, a
// parse error (E4), as boost::json's type errors were for RDKit 2026.03.
void need(bool ok, const std::string &what) {
  if (!ok) {
    throw RDKit::FileParseException(
        "the JSON " + what + ", which RDKit 2025.03 reads without checking");
  }
}

// Checks, before RDKit's JSON parser (MolInterchange/Parser.cpp) reads a
// document, what it trusts, following it read by read. Each check returns
// false where RDKit would throw for what it reads (a value of the wrong
// type it checks, a missing key, an index out of range), and RDKit reads no
// further; what RDKit would read without a check is refused with
// ValueErrorException. Where RDKit's throw depends on more than the document
// (a stereo atom that is no neighbor, an invalid sanitization), the checks go
// on, so they may refuse a document RDKit would also reject.
class DocumentChecker {
 public:
  // As DocToMols.
  DocumentChecker(const bj::value &doc, const MI::JSONParseParameters &params)
      : d_params(params) {
    if (!doc.is_object()) {
      return;
    }
    if (const auto *header = member(doc, "commonchem")) {
      const auto *version = member(*header, "version");
      if (!version) {
        return;
      }
      need(rj_int32(version), "has a commonchem version that is not an int");
      if (version->as_int64() != MI::currentMolJSONVersion) {
        return;
      }
    } else if (const auto *header = member(doc, "rdkitjson")) {
      const auto *version = member(*header, "version");
      if (!version) {
        return;
      }
      need(rj_int32(version), "has an rdkitjson version that is not an int");
      if (version->as_int64() > MI::currentRDKitJSONVersion ||
          version->as_int64() < 10) {
        return;
      }
    } else {
      return;
    }
    // GetValueByPointer("/defaults/atom") and ".../bond".
    if (const auto *defaults = member(doc, "defaults")) {
      d_atom_defaults = member(*defaults, "atom");
      if (d_atom_defaults && !d_atom_defaults->is_object()) {
        return;
      }
      d_bond_defaults = member(*defaults, "bond");
      if (d_bond_defaults && !d_bond_defaults->is_object()) {
        return;
      }
    }
    if (const auto *mols = member(doc, "molecules")) {
      if (!mols->is_array()) {
        return;
      }
      for (const auto &mol : mols->as_array()) {
        if (!check_mol(mol)) {
          return;
        }
      }
    }
  }

 private:
  // getIntDefaultValue and the like: false where RDKit throws for the
  // value's type.
  template <class P>
  static bool typed(const bj::value &from, const bj::value *defaults,
                    const char *key, P &&is) {
    if (const auto *v = member(from, key)) {
      return is(v);
    }
    if (const auto *v = defaults ? member(*defaults, key) : nullptr) {
      return is(v);
    }
    return true;
  }
  static std::string string_default(const bj::value &from,
                                    const bj::value *defaults,
                                    const char *key) {
    const auto *v = member(from, key);
    if (!v && defaults) {
      v = member(*defaults, key);
    }
    return v ? text(v->as_string()) : std::string();
  }
  static std::int64_t int_default(const bj::value &from,
                                  const bj::value *defaults, const char *key) {
    const auto *v = member(from, key);
    if (!v && defaults) {
      v = member(*defaults, key);
    }
    return v ? v->as_int64() : 0;
  }

  // As processMol, for a molecule and for a recursive query's; false where
  // RDKit throws.
  bool check_mol(const bj::value &mol) {
    need(mol.is_object(), "gives a molecule that is not an object");
    if (const auto *name = member(mol, "name")) {
      need(rj_string(name), "gives a molecule a name that is not a string");
    }
    const auto *atoms = member(mol, "atoms");
    const auto *bonds = member(mol, "bonds");
    if (!atoms || !bonds) {
      return false;
    }
    need(atoms->is_array(), "gives a molecule atoms that are not an array");
    need(bonds->is_array(), "gives a molecule bonds that are not an array");
    const std::size_t num_atoms = atoms->as_array().size();
    const std::size_t num_bonds = bonds->as_array().size();

    // readAtom.
    for (const auto &atom : atoms->as_array()) {
      need(atom.is_object(), "gives an atom that is not an object");
      for (const char *key : {"z", "impHs", "chg", "nRad", "isotope"}) {
        if (std::string_view(key) == "impHs" && !d_params.useHCounts) {
          continue;
        }
        if (!typed(atom, d_atom_defaults, key, rj_int)) {
          return false;
        }
      }
      if (!typed(atom, d_atom_defaults, "stereo", rj_string) ||
          !MI::chilookup.count(
              string_default(atom, d_atom_defaults, "stereo"))) {
        return false;
      }
    }

    // readBond, then readBondStereo for every bond once one has stereo.
    std::set<std::pair<std::int64_t, std::int64_t>> pairs;
    std::vector<std::string> stereo;
    bool need_stereo_loop = false;
    for (const auto &bond : bonds->as_array()) {
      need(bond.is_object(), "gives a bond that is not an object");
      const auto *ends = member(bond, "atoms");
      need(rj_array(ends) && ends->as_array().size() >= 2 &&
               rj_int32(&ends->as_array()[0]) && rj_int32(&ends->as_array()[1]),
           "gives a bond without two int atoms");
      // RWMol::addBond checks the atoms, as unsigned.
      auto a = static_cast<std::int64_t>(
          static_cast<unsigned>(ends->as_array()[0].as_int64()));
      auto b = static_cast<std::int64_t>(
          static_cast<unsigned>(ends->as_array()[1].as_int64()));
      if (a >= static_cast<std::int64_t>(num_atoms) ||
          b >= static_cast<std::int64_t>(num_atoms) || a == b ||
          !pairs.insert({std::min(a, b), std::max(a, b)}).second) {
        return false;
      }
      if (!typed(bond, d_bond_defaults, "bo", rj_int) ||
          !MI::bolookup.count(static_cast<unsigned>(
              int_default(bond, d_bond_defaults, "bo")))) {
        return false;
      }
      if (!typed(bond, d_bond_defaults, "stereo", rj_string)) {
        return false;
      }
      stereo.push_back(string_default(bond, d_bond_defaults, "stereo"));
      need_stereo_loop |= stereo.back() != "unspecified";
    }
    if (need_stereo_loop) {
      for (std::size_t i = 0; i < num_bonds; ++i) {
        if (stereo[i] == "unspecified") {
          continue;
        }
        if (!MI::stereoBondlookup.count(stereo[i])) {
          return false;
        }
        if (const auto *ids = member(bonds->as_array()[i], "stereoAtoms")) {
          need(rj_array(ids), "gives a bond stereo atoms not in an array");
          // readBondStereo reads aids[0] and aids[1] without checking that
          // they exist.
          if (ids->as_array().size() < 2) {
            refuse("gives a bond fewer than two stereo atoms");
          }
          need(rj_int32(&ids->as_array()[0]) && rj_int32(&ids->as_array()[1]),
               "gives a bond stereo atoms that are not ints");
        } else if (stereo[i] == "cis" || stereo[i] == "trans") {
          return false;
        }
      }
    }

    // readStereoGroups: the bonds the groups hold pointers to.
    std::set<unsigned> group_bonds;
    if (const auto *groups = member(mol, "stereoGroups")) {
      need(groups->is_array(), "gives stereo groups that are not an array");
      for (const auto &group : groups->as_array()) {
        need(group.is_object(), "gives a stereo group that is not an object");
        const auto *type = member(group, "type");
        const auto *ids = member(group, "atoms");
        const auto *bond_ids = member(group, "bonds");
        if (!type || (!ids && !bond_ids)) {
          return false;
        }
        need(rj_string(type),
             "gives a stereo group a type that is not a string");
        auto found = MI::stereoGrouplookup.find(text(type->as_string()));
        if (found == MI::stereoGrouplookup.end()) {
          return false;
        }
        if (const auto *id = member(group, "id");
            id && found->second != RDKit::StereoGroupType::STEREO_ABSOLUTE) {
          need(rj_int32(id), "gives a stereo group an id that is not unsigned");
        }
        for (const auto &[list, count, what] :
             {std::tuple{ids, num_atoms, "atom"},
              std::tuple{bond_ids, num_bonds, "bond"}}) {
          if (!list) {
            continue;
          }
          need(list->is_array(), std::string("gives a stereo group ") + what +
                                     "s that are not an array");
          std::set<std::int64_t> seen;
          for (const auto &id : list->as_array()) {
            need(rj_int32(&id), std::string("gives a stereo group an ") + what +
                                    " that is not unsigned");
            if (static_cast<std::size_t>(id.as_int64()) >= count) {
              return false;
            }
            if (!seen.insert(id.as_int64()).second) {
              refuse(std::string("gives a stereo group an ") + what +
                     " twice, which RDKit 2025.03 keeps, leaving it to dangle "
                     "once removed");
            }
            if (list == bond_ids) {
              group_bonds.insert(static_cast<unsigned>(id.as_int64()));
            }
          }
        }
      }
    }

    // readSubstanceGroups.
    if (const auto *sgroups = member(mol, "substanceGroups")) {
      need(sgroups->is_array(), "gives substance groups that are not an array");
      for (const auto &sgroup : sgroups->as_array()) {
        if (!check_substance_group(sgroup, num_atoms, num_bonds)) {
          return false;
        }
      }
    }

    // readConformer, for each conformer.
    if (const auto *confs = member(mol, "conformers");
        confs && d_params.parseConformers) {
      need(confs->is_array(), "gives conformers that are not an array");
      for (const auto &conf : confs->as_array()) {
        need(conf.is_object(), "gives a conformer that is not an object");
        const auto *dim = member(conf, "dim");
        if (!dim) {
          return false;
        }
        need(rj_int32(dim), "gives a conformer a dim that is not an int");
        if (dim->as_int64() != 2 && dim->as_int64() != 3) {
          return false;
        }
        const auto *coords = member(conf, "coords");
        if (!coords) {
          return false;
        }
        need(coords->is_array(),
             "gives conformer coords that are not an array");
        for (const auto &point : coords->as_array()) {
          need(point.is_array(), "gives a coordinate that is not an array");
          if (static_cast<std::int64_t>(point.as_array().size()) !=
              dim->as_int64()) {
            return false;
          }
          for (const auto &x : point.as_array()) {
            need(rj_number(&x), "gives a coordinate that is not a number");
          }
        }
        if (coords->as_array().size() != num_atoms) {
          return false;
        }
      }
    }
    if (const auto *props = member(mol, "properties");
        props && d_params.parseProperties) {
      need(props->is_object(), "gives properties that are not an object");
    }

    // The extensions, in order.
    std::size_t replaced_bonds = 0;
    bool rings = false;
    if (const auto *extensions = member(mol, "extensions")) {
      need(extensions->is_array(), "gives extensions that are not an array");
      for (const auto &extension : extensions->as_array()) {
        need(extension.is_object(), "gives an extension that is not an object");
        const auto *name = member(extension, "name");
        if (!name) {
          return false;
        }
        need(rj_string(name), "gives an extension a name that is not a string");
        const auto which = text(name->as_string());
        if (which != "rdkitRepresentation" && which != "partialCharges" &&
            which != "rdkitQueries") {
          continue;
        }
        const auto *version = member(extension, "formatVersion");
        if (!version) {
          return false;
        }
        need(rj_int32(version),
             "gives an extension a formatVersion that is not an int");
        if (which == "rdkitRepresentation") {
          if (version->as_int64() > MI::currentRDKitRepresentationVersion) {
            continue;
          }
          if (!check_representation(extension, num_atoms, num_bonds, pairs,
                                    rings)) {
            return false;
          }
        } else if (which == "partialCharges") {
          if (version->as_int64() > MI::currentChargeRepresentationVersion) {
            continue;
          }
          if (const auto *values = member(extension, "values")) {
            need(values->is_array(),
                 "gives partial charges that are not an array");
            if (values->as_array().size() != num_atoms) {
              return false;
            }
            for (const auto &value : values->as_array()) {
              if (!rj_double(&value)) {
                return false;
              }
            }
          }
        } else {
          if (version->as_int64() > MI::currentQueryRepresentationVersion) {
            continue;
          }
          std::size_t replaced = 0;
          if (!check_queries(extension, num_atoms, num_bonds, replaced)) {
            return false;
          }
          replaced_bonds = std::max(replaced_bonds, replaced);
        }
      }
    }
    // readQueries replaces the first bonds, one per bond query
    // (RWMol::replaceBond), deleting each old bond, but replaceBond, unlike
    // replaceAtom, leaves stereo groups pointing at the bond it deletes.
    for (auto idx : group_bonds) {
      if (idx < replaced_bonds) {
        refuse(
            "gives a bond in a stereo group a query; RDKit's replaceBond "
            "leaves the group pointing at the bond it deletes");
      }
    }
    return true;
  }

  // As readSubstanceGroups for one group. RDKit checks the atoms, bonds,
  // and parent atoms as it sets them, but not the brackets, the bonds of the
  // cstates, or the atoms of the attachment points, which the crate trusts
  // as SubstanceGroup's own setters keep them (substance_group.cpp).
  bool check_substance_group(const bj::value &sgroup, std::size_t num_atoms,
                             std::size_t num_bonds) {
    need(sgroup.is_object(), "gives a substance group that is not an object");
    const auto *props = member(sgroup, "properties");
    if (!props) {
      return false;
    }
    need(props->is_object(),
         "gives substance group properties that are not an object");
    const auto *type = member(*props, "TYPE");
    if (!type) {
      return false;
    }
    need(rj_string(type),
         "gives a substance group a TYPE that is not a string");
    if (!RDKit::SubstanceGroupChecks::isValidType(text(type->as_string()))) {
      return false;
    }
    for (const auto &[key, valid] :
         {std::pair{"SUBTYPE", &RDKit::SubstanceGroupChecks::isValidSubType},
          std::pair{"CONNECT",
                    &RDKit::SubstanceGroupChecks::isValidConnectType}}) {
      // parseProperties keeps only ints, doubles, and strings; no number
      // is a valid value.
      if (const auto *v = member(*props, key)) {
        if (rj_int32(v) || rj_double(v) ||
            (rj_string(v) && !valid(text(v->as_string())))) {
          return false;
        }
      }
    }
    for (const auto &[key, count] :
         {std::pair{"atoms", num_atoms}, std::pair{"bonds", num_bonds},
          std::pair{"parentAtoms", num_atoms}}) {
      if (const auto *ids = member(sgroup, key)) {
        need(ids->is_array(), std::string("gives a substance group ") + key +
                                  " that are not an array");
        for (const auto &id : ids->as_array()) {
          need(rj_int32(&id), std::string("gives a substance group ") + key +
                                  " that are not unsigned");
          if (static_cast<std::size_t>(id.as_int64()) >= count) {
            return false;
          }
        }
      }
    }
    auto point = [](const bj::value &p, const char *what) {
      need(p.is_array(),
           std::string("gives a ") + what + " that is not an array");
      if (p.as_array().size() != 3) {
        return false;
      }
      for (const auto &x : p.as_array()) {
        need(rj_number(&x), std::string("gives a ") + what +
                                " coordinate that is not a number");
      }
      return true;
    };
    if (const auto *brackets = member(sgroup, "brackets")) {
      need(brackets->is_array(), "gives brackets that are not an array");
      for (const auto &bracket : brackets->as_array()) {
        need(bracket.is_array(), "gives a bracket that is not an array");
        std::size_t n = 0;
        for (const auto &p : bracket.as_array()) {
          if (!point(p, "bracket point")) {
            return false;
          }
          // RDKit writes each point into a Bracket, an array of three, past
          // its end for a fourth.
          if (++n > 3) {
            refuse(
                "gives a substance group a bracket of more than three "
                "points");
          }
        }
      }
    }
    if (const auto *cstates = member(sgroup, "cstates")) {
      need(cstates->is_array(), "gives cstates that are not an array");
      for (const auto &cstate : cstates->as_array()) {
        need(cstate.is_object(), "gives a cstate that is not an object");
        const auto *bond = member(cstate, "bond");
        need(rj_int32(bond), "gives a cstate a bond that is not unsigned");
        if (const auto *vector = member(cstate, "vector")) {
          if (!point(*vector, "cstate vector")) {
            return false;
          }
        }
        if (static_cast<std::size_t>(bond->as_int64()) >= num_bonds) {
          refuse(
              "gives a substance group a cstate for a bond the molecule "
              "lacks");
        }
      }
    }
    if (const auto *points = member(sgroup, "attachPoints")) {
      need(points->is_array(), "gives attachment points that are not an array");
      for (const auto &attach : points->as_array()) {
        need(attach.is_object(),
             "gives an attachment point that is not an object");
        const auto *a_idx = member(attach, "aIdx");
        need(rj_int32(a_idx),
             "gives an attachment point an aIdx that is not unsigned");
        const auto *lv_idx = member(attach, "lvIdx");
        if (lv_idx) {
          need(rj_int32(lv_idx),
               "gives an attachment point an lvIdx that is not unsigned");
        }
        if (const auto *id = member(attach, "id")) {
          need(rj_string(id),
               "gives an attachment point an id that is not a string");
        }
        if (static_cast<std::size_t>(a_idx->as_int64()) >= num_atoms) {
          refuse(
              "gives a substance group an attachment point at an atom the "
              "molecule lacks");
        }
        if (lv_idx) {
          // As RDKit stores it, in AttachPoint's int.
          auto lv = static_cast<int>(static_cast<unsigned>(lv_idx->as_int64()));
          if (lv != -1 &&
              (lv < 0 || static_cast<std::size_t>(lv) >= num_atoms)) {
            refuse(
                "gives a substance group an attachment point that leaves an "
                "atom the molecule lacks");
          }
        }
      }
    }
    return true;
  }

  // As readRDKitRepresentation; false where RDKit throws.
  static bool check_representation(
      const bj::value &rep, std::size_t num_atoms, std::size_t num_bonds,
      const std::set<std::pair<std::int64_t, std::int64_t>> &pairs,
      bool &rings) {
    auto index = [](const bj::value &v, std::size_t n) {
      return static_cast<unsigned>(v.as_int64()) < n;
    };
    for (const auto &[key, n] : {std::pair{"aromaticAtoms", num_atoms},
                                 std::pair{"aromaticBonds", num_bonds}}) {
      if (const auto *list = member(rep, key)) {
        need(list->is_array(),
             std::string("gives ") + key + " not in an array");
        for (const auto &v : list->as_array()) {
          if (!rj_int(&v) || !index(v, n)) {
            return false;
          }
        }
      }
    }
    if (const auto *ranks = member(rep, "cipRanks")) {
      need(ranks->is_array(), "gives CIP ranks not in an array");
      std::size_t i = 0;
      for (const auto &v : ranks->as_array()) {
        if (!rj_int(&v) || i++ >= num_atoms) {
          return false;
        }
      }
    }
    if (const auto *codes = member(rep, "cipCodes")) {
      need(codes->is_array(), "gives CIP codes not in an array");
      for (const auto &v : codes->as_array()) {
        if (!v.is_array()) {
          return false;
        }
        const auto &pair = v.as_array();
        need(pair.size() >= 2 && rj_int32(&pair[0]) && rj_string(&pair[1]),
             "gives a CIP code that is not an int and a string");
        if (!index(pair[0], num_atoms)) {
          return false;
        }
      }
    }
    if (const auto *ring_list = member(rep, "atomRings")) {
      need(ring_list->is_array(), "gives atom rings not in an array");
      // RDKit initializes the ring information once.
      if (rings) {
        return false;
      }
      rings = true;
      for (const auto &ring : ring_list->as_array()) {
        if (!ring.is_array()) {
          return false;
        }
        const auto &atoms = ring.as_array();
        // RDKit reads an empty ring's element sz - 1, past the end.
        if (atoms.empty()) {
          refuse("gives an empty atom ring");
        }
        for (std::size_t i = 0; i < atoms.size(); ++i) {
          const auto &next = atoms[(i + 1) % atoms.size()];
          need(rj_int32(&atoms[i]) && rj_int32(&next),
               "gives a ring atom that is not an int");
          if (!index(atoms[i], num_atoms) || !index(next, num_atoms)) {
            return false;
          }
          auto a = static_cast<std::int64_t>(atoms[i].as_int64());
          auto b = static_cast<std::int64_t>(next.as_int64());
          if (!pairs.count({std::min(a, b), std::max(a, b)})) {
            return false;
          }
        }
      }
    }
    return true;
  }

  // As readQueries; counts the bond queries, which replace bonds from the
  // first; false where RDKit throws.
  bool check_queries(const bj::value &rep, std::size_t num_atoms,
                     std::size_t num_bonds, std::size_t &replaced_bonds) {
    for (const auto &[key, n] : {std::pair{"atomQueries", num_atoms},
                                 std::pair{"bondQueries", num_bonds}}) {
      const auto *queries = member(rep, key);
      if (!queries) {
        continue;
      }
      need(queries->is_array(),
           std::string("gives ") + key + " not in an array");
      std::size_t idx = 0;
      for (const auto &query : queries->as_array()) {
        if (!query.is_object()) {
          return false;
        }
        if (!member(query, "tag")) {
          continue;
        }
        if (idx++ >= n) {
          return false;
        }
        bool ok = std::string_view(key) == "atomQueries"
                      ? check_atom_query(query)
                      : check_bond_query(query, false);
        if (std::string_view(key) == "bondQueries") {
          ++replaced_bonds;
        }
        if (!ok) {
          return false;
        }
      }
    }
    return true;
  }

  // readBaseQuery: its values, each an int RDKit reads unchecked and, where
  // RDKit can read it, within U7's bound; false where RDKit throws.
  static bool check_base_query(const bj::value &query, int tag) {
    if (!member(query, "descr")) {
      return false;
    }
    auto value = [&](const char *key, bool required) {
      const auto *v = member(query, key);
      if (!v && !required) {
        return;
      }
      need(rj_int32(v),
           std::string("gives a query a ") + key + " that is not an int");
      if (std::string_view(key) != "ends" &&
          (v->as_int64() > kMaxQueryValue || v->as_int64() < -kMaxQueryValue)) {
        refuse(std::string("gives a query a ") + key +
               " beyond 2^30, where RDKit's int comparison can overflow");
      }
    };
    if (is_equality_tag(tag)) {
      value("val", true);
      value("tol", false);
    } else if (tag == MolPickler::QUERY_RANGE) {
      value("lower", true);
      value("upper", true);
      value("tol", false);
      value("ends", false);
    } else if (tag == MolPickler::QUERY_SET) {
      if (const auto *set = member(query, "set")) {
        need(set->is_array(), "gives a query a set that is not an array");
        for (const auto &v : set->as_array()) {
          need(rj_int32(&v), "gives a query a set member that is not an int");
        }
      }
    } else if (!is_boolean_tag(tag) && tag != MolPickler::QUERY_NULL) {
      return false;  // an unknown tag
    }
    return true;
  }

  // finishQuery's own fields: the description, type, and negation it reads
  // unchecked.
  static std::string check_finish(const bj::value &query) {
    const auto *descr = member(query, "descr");
    need(rj_string(descr), "gives a query a description that is not a string");
    if (const auto *type = member(query, "type")) {
      need(rj_string(type), "gives a query a type that is not a string");
    }
    if (const auto *negated = member(query, "negated")) {
      need(rj_bool(negated), "gives a query a negation that is not a bool");
    }
    return text(descr->as_string());
  }

  // finishQuery's children: two or more for AND, OR, and XOR, as RDKit's
  // writers read (SmartsWrite.cpp, MolFileWriter.cpp, QueryOps.cpp,
  // PatternFingerprints.cpp take the first two without checking), and none
  // for the queries RDKit makes without them. readQuery reads each child's
  // tag as an object's member.
  template <class F>
  static bool check_children(const bj::value &query, int tag, F &&check) {
    const auto *children = member(query, "children");
    if (children) {
      need(children->is_array(), "gives query children not in an array");
    }
    std::size_t n = children ? children->as_array().size() : 0;
    if (is_boolean_tag(tag) && n < 2) {
      refuse("gives an AND, OR, or XOR query fewer than two children");
    }
    if (!is_boolean_tag(tag) && n > 0) {
      refuse("gives children to a query that is not AND, OR, or XOR");
    }
    if (children) {
      for (const auto &child : children->as_array()) {
        need(child.is_object(), "gives a query child that is not an object");
        if (!member(child, "tag") || !check(child)) {
          return false;
        }
      }
    }
    return true;
  }

  // As readQuery for atoms.
  bool check_atom_query(const bj::value &query) {
    const auto *tag_value = member(query, "tag");
    need(rj_int32(tag_value), "gives a query a tag that is not an int");
    const int tag = static_cast<int>(tag_value->as_int64());
    if (tag == MolPickler::QUERY_RECURSIVE) {
      const auto *mol = member(query, "subquery");
      if (!mol || !check_mol(*mol)) {
        return false;
      }
      // RecursiveStructureQuery's matcher reads the molecule's first atom.
      if (member(*mol, "atoms")->as_array().empty()) {
        refuse("gives a recursive query a molecule without atoms");
      }
    } else if (tag == MolPickler::QUERY_ATOMRING) {
      const auto *val = member(query, "val");
      need(rj_int32(val), "gives a ring query a val that is not an int");
      if (const auto *tol = member(query, "tol")) {
        need(rj_int32(tol), "gives a ring query a tol that is not an int");
      }
      for (const auto *v : {val, member(query, "tol")}) {
        if (v && (v->as_int64() > kMaxQueryValue ||
                  v->as_int64() < -kMaxQueryValue)) {
          refuse(
              "gives a query a value beyond 2^30, where RDKit's int "
              "comparison can overflow");
        }
      }
    } else if (!check_base_query(query, tag)) {
      return false;
    }
    const auto descr = check_finish(query);
    // finalizeQueryFromDescription strips a range query's prefix, which the
    // SMARTS parser adds for its range, less-than, and greater-than queries
    // (smarts.yy) and gives those classes.
    int prefixed = -1;
    std::string base = descr;
    for (const auto &[prefix, prefix_tag] :
         {std::pair<std::string, int>{"range_", MolPickler::QUERY_RANGE},
          {"less_", MolPickler::QUERY_LESSEQUAL},
          {"greater_", MolPickler::QUERY_GREATEREQUAL}}) {
      if (descr.rfind(prefix, 0) == 0) {
        prefixed = prefix_tag;
        base = descr.substr(prefix.size());
        break;
      }
    }
    bool known = false;
    if (prefixed >= 0) {
      // RDKit 2025.03 reads an AtomRingSize query's value through
      // static_cast<ATOM_EQUALS_QUERY *> (finalizeQueryFromDescription in
      // QueryOps.cpp), and RangeQuery is no EqualityQuery, as the less and
      // greater queries are (rdkit/rdkit#9172, fixed in 2026.03.1).
      known = is_atom_value_query(base) && tag == prefixed &&
              !(prefixed == MolPickler::QUERY_RANGE && base == "AtomRingSize");
    } else if (is_atom_value_query(base)) {
      // AtomRingQuery is an EqualityQuery, which RDKit makes for
      // AtomInNRings and AtomRingBondCount (QueryOps.h, QueryOps.cpp).
      known = is_equality_tag(tag) || tag == MolPickler::QUERY_ATOMRING;
    } else if (base == "AtomNull") {
      known = tag == MolPickler::QUERY_NULL;
    } else if (base == "RecursiveStructure") {
      known = tag == MolPickler::QUERY_RECURSIVE;
    } else if (base == "AtomAnd") {
      known = tag == MolPickler::QUERY_AND;
    } else if (base == "AtomOr") {
      known = tag == MolPickler::QUERY_OR;
    } else if (base == "AtomXor") {
      known = tag == MolPickler::QUERY_XOR;
    }
    // RDKit takes the query's class from its tag and how to read it from
    // its description; RDKit's code that reads queries takes the
    // description to imply the class (a RecursiveStructure query is cast to
    // RecursiveStructureQuery in SubstructMatch.cpp, an AtomAtomicNum query
    // to EqualityQuery in SmartsWrite.cpp), and without a data function for
    // its description a query dereferences null when it matches.
    if (!known) {
      refuse("gives an atom a query RDKit does not make: description '" +
             descr + "' with tag " + std::to_string(tag));
    }
    return check_children(query, tag, [&](const bj::value &child) {
      return check_atom_query(child);
    });
  }

  // As readQuery for bonds. A child of an AND or OR query is not a null
  // query: SmartsWrite.cpp casts each child that is not AND or OR to
  // EqualityQuery. XOR bond queries are not read: QueryBond.cpp's matcher
  // casts the other bond's query to EqualityQuery whatever its class.
  static bool check_bond_query(const bj::value &query, bool in_and_or) {
    const auto *tag_value = member(query, "tag");
    need(rj_int32(tag_value), "gives a query a tag that is not an int");
    const int tag = static_cast<int>(tag_value->as_int64());
    if (!check_base_query(query, tag)) {
      return false;
    }
    const auto descr = check_finish(query);
    bool known = false;
    if (is_bond_value_query(descr)) {
      known = is_equality_tag(tag);
    } else if (descr == "BondNull") {
      known = tag == MolPickler::QUERY_NULL && !in_and_or;
    } else if (descr == "BondAnd") {
      known = tag == MolPickler::QUERY_AND;
    } else if (descr == "BondOr") {
      known = tag == MolPickler::QUERY_OR;
    }
    if (!known) {
      refuse("gives a bond a query RDKit does not make: description '" + descr +
             "' with tag " + std::to_string(tag) +
             (in_and_or ? " in an AND or OR query" : ""));
    }
    return check_children(query, tag, [&](const bj::value &child) {
      return check_bond_query(child, true);
    });
  }

  const MI::JSONParseParameters &d_params;
  const bj::value *d_atom_defaults = nullptr;
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

// RDKit 2025.03's readSubstanceGroups (MolInterchange/Parser.cpp) leaves an
// attachment point's lvIdx uninitialized when the point has none, as RDKit's
// writer leaves out -1; RDKit 2025.09 makes it -1. The checked document gets
// 4294967295 there, which RDKit reads with GetUint as -1 into AttachPoint's
// int.
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
          point.as_object()["lvIdx"] = std::uint64_t{4294967295u};
        }
      }
    }
  }
}

// RDKit 2025.03 parses the document it is given with RapidJSON
// (JSONDataToMols), without checking for errors. This parses it with
// boost::json, which refuses what is not JSON (and nesting deeper than 32,
// which RapidJSON's recursive parser would follow down the stack), checks it,
// and hands RDKit the checked document, serialized again, so RDKit parses
// exactly what was checked (U2), with the types DocumentChecker expects.
std::unique_ptr<MolList> json_to_mols(rust::Str json_block,
                                      const JSONParseParametersData &params) {
  auto p = parse_params(params);
  std::vector<boost::shared_ptr<RDKit::ROMol>> mols;
  try {
    bj::value doc =
        bj::parse(std::string_view(json_block.data(), json_block.size()));
    DocumentChecker check(doc, p);
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
