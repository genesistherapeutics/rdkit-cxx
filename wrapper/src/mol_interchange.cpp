// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_interchange.h"

#include <climits>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/MolInterchange/details.h>
#include <GraphMol/MolPickler.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>

// RDKit 2024.09 parses with rapidjson, which it does not install, so the
// check parses with boost::json, header-only, its sources compiled once,
// here.
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

// RDKit reads conformer coordinates with rapidjson's GetFloat, a
// static_cast<float> of the double, undefined beyond float's range. Below
// FLT_MAX (3.40282e38) by far more than rapidjson's parse can differ from
// boost::json's.
constexpr double kMaxCoordinate = 3.4e38;

[[noreturn]] void refuse(const std::string &why) {
  throw ValueErrorException("the JSON " + why);
}

// RDKit's reader throws here, and reads nothing further; the check stops.
struct RDKitThrows {};
[[noreturn]] void rdkit_throws() { throw RDKitThrows{}; }

// RDKit reads the JSON's strings with GetString(), as C strings, so they
// end at a NUL.
std::string text(const bj::string &s) { return std::string(s.c_str()); }

const bj::value *member(const bj::object &obj, std::string_view key) {
  auto it = obj.find(key);
  return it == obj.end() ? nullptr : &it->value();
}

// What RDKit's rapidjson says of the value it parses where boost::json
// parsed `v` (json_to_mols hands RDKit boost::json's serialization).
// boost::json writes an int64 or a uint64 as its digits, which rapidjson
// types by value (IsInt within int, IsUint within unsigned), and a double
// always with an exponent, which rapidjson types as a double, not an int.
bool rj_int(const bj::value &v) {
  return v.is_int64() && v.get_int64() >= INT_MIN && v.get_int64() <= INT_MAX;
}

bool rj_uint(const bj::value &v) {
  return v.is_int64() && v.get_int64() >= 0 && v.get_int64() <= UINT_MAX;
}

// The value RDKit reads without checking that it is there or has the type
// it reads it as: rapidjson's operator[] gives a null value for a missing
// key, and its Get* and member lookups on a value of another type read the
// value's storage as that type (RAPIDJSON_ASSERT is compiled out).
const bj::value &present(const bj::value *v, const std::string &what) {
  if (!v) {
    refuse("lacks " + what + ", which RDKit reads without checking");
  }
  return *v;
}

const bj::object &object_of(const bj::value &v, const std::string &what) {
  if (!v.is_object()) {
    refuse("gives " + what + " that is not an object");
  }
  return v.get_object();
}

const bj::array &array_of(const bj::value &v, const std::string &what) {
  if (!v.is_array()) {
    refuse("gives " + what + " that is not an array");
  }
  return v.get_array();
}

int int_of(const bj::value &v, const std::string &what) {
  if (!rj_int(v)) {
    refuse("gives " + what + " that is not an int");
  }
  return static_cast<int>(v.get_int64());
}

unsigned uint_of(const bj::value &v, const std::string &what) {
  if (!rj_uint(v)) {
    refuse("gives " + what + " that is not an unsigned int");
  }
  return static_cast<unsigned>(v.get_int64());
}

std::string string_of(const bj::value &v, const std::string &what) {
  if (!v.is_string()) {
    refuse("gives " + what + " that is not a string");
  }
  return text(v.get_string());
}

void bool_of(const bj::value &v, const std::string &what) {
  if (!v.is_bool()) {
    refuse("gives " + what + " that is not a bool");
  }
}

// GetDouble.
double number_of(const bj::value &v, const std::string &what) {
  if (!v.is_number()) {
    refuse("gives " + what + " that is not a number");
  }
  return v.to_number<double>();
}

// The array RDKit reads elements 0 to n - 1 of without checking its size.
const bj::array &array_of_at_least(const bj::value &v, std::size_t n,
                                   const std::string &what) {
  const auto &list = array_of(v, what);
  if (list.size() < n) {
    refuse("gives " + what + " of fewer than " + std::to_string(n) + " values");
  }
  return list;
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

// Checks, before RDKit's JSON parser (MolInterchange/Parser.cpp, which
// reads with rapidjson) reads a document, what it trusts, following it read
// by read. RDKit 2024.09 checks few of the types it reads: rapidjson's
// operator[] gives a null value for a key that is missing, and its GetInt,
// GetUint, GetDouble, GetFloat, GetString, GetBool, GetArray, GetObject,
// HasMember, FindMember, and operator[] read a value of another type, or an
// array's element past its end, as if it were one. Each such read is
// checked here, as rapidjson types the value it reads (rj_int and rj_uint),
// and refused with ValueErrorException, as are the indices RDKit trusts.
// Where RDKit would throw for what it reads, the check stops there
// (RDKitThrows) and leaves the error to RDKit, which reads no further.
class DocumentChecker {
 public:
  // As DocToMols.
  DocumentChecker(const bj::value &doc, const MI::JSONParseParameters &params)
      : d_params(params) {
    try {
      check_doc(doc);
    } catch (const RDKitThrows &) {
    }
  }

 private:
  void check_doc(const bj::value &doc) {
    if (!doc.is_object()) {
      rdkit_throws();
    }
    const auto &root = doc.get_object();
    if (const auto *header = member(root, "commonchem")) {
      if (!header->is_object() || !member(header->get_object(), "version")) {
        rdkit_throws();
      }
      int version = int_of(*member(header->get_object(), "version"),
                           "a commonchem version");
      if (version != MI::currentMolJSONVersion) {
        rdkit_throws();
      }
    } else if (const auto *header = member(root, "rdkitjson")) {
      if (!header->is_object() || !member(header->get_object(), "version")) {
        rdkit_throws();
      }
      int version = int_of(*member(header->get_object(), "version"),
                           "an rdkitjson version");
      if (version > MI::currentRDKitJSONVersion || version < 10) {
        rdkit_throws();
      }
    } else {
      rdkit_throws();
    }
    // RDKit finds the defaults with a JSON pointer, which finds none in
    // defaults that are not an object, and an empty DefaultValueCache finds
    // none in its null value.
    if (const auto *defaults = member(root, "defaults");
        defaults && defaults->is_object()) {
      const auto *atom = member(defaults->get_object(), "atom");
      if (atom && !atom->is_object()) {
        rdkit_throws();
      }
      const auto *bond = member(defaults->get_object(), "bond");
      if (bond && !bond->is_object()) {
        rdkit_throws();
      }
      d_bond_defaults = bond ? &bond->get_object() : nullptr;
    }
    if (const auto *mols = member(root, "molecules")) {
      if (!mols->is_array()) {
        rdkit_throws();
      }
      for (const auto &mol : mols->get_array()) {
        check_mol(mol);
      }
    }
  }

  // A bond's stereo, as getStringDefaultValue reads it.
  std::string bond_stereo(const bj::object &bond) const {
    const bj::value *stereo = member(bond, "stereo");
    if (!stereo && d_bond_defaults) {
      stereo = member(*d_bond_defaults, "stereo");
    }
    if (!stereo) {
      return std::string();
    }
    if (!stereo->is_string()) {
      rdkit_throws();
    }
    return text(stereo->get_string());
  }

  // As processMol, for a molecule and for a recursive query's.
  void check_mol(const bj::value &mol_value) {
    const auto &mol = object_of(mol_value, "a molecule");
    if (const auto *name = member(mol, "name")) {
      string_of(*name, "a molecule's name");
    }
    const auto *atoms_value = member(mol, "atoms");
    const auto *bonds_value = member(mol, "bonds");
    if (!atoms_value || !bonds_value) {
      rdkit_throws();
    }
    // readAtom finds the atom's values in it.
    const auto &atoms = array_of(*atoms_value, "a molecule's atoms");
    for (const auto &atom : atoms) {
      object_of(atom, "an atom");
    }
    const std::size_t num_atoms = atoms.size();

    // readBond, then readBondStereo for every bond once one has stereo.
    const auto &bonds = array_of(*bonds_value, "a molecule's bonds");
    const std::size_t num_bonds = bonds.size();
    std::vector<std::string> stereo;
    bool need_stereo_loop = false;
    for (const auto &bond_value : bonds) {
      const auto &bond = object_of(bond_value, "a bond");
      // readBond reads aids[0] and aids[1] without checking that they
      // exist; RWMol::addBond checks the atoms they name.
      const auto &ends =
          array_of_at_least(present(member(bond, "atoms"), "a bond's atoms"), 2,
                            "a bond's atoms");
      int_of(ends[0], "a bond's atom");
      int_of(ends[1], "a bond's atom");
      auto s = bond_stereo(bond);
      need_stereo_loop |= s != "unspecified";
      stereo.push_back(std::move(s));
    }
    if (need_stereo_loop) {
      for (std::size_t i = 0; i < num_bonds; ++i) {
        if (stereo[i] == "unspecified") {
          continue;
        }
        if (MI::stereoBondlookup.find(stereo[i]) ==
            MI::stereoBondlookup.end()) {
          rdkit_throws();
        }
        const auto &bond = bonds[i].get_object();
        if (const auto *atoms_of = member(bond, "stereoAtoms")) {
          // readBondStereo reads aids[0] and aids[1] without checking that
          // they exist; setStereoAtoms checks the atoms they name.
          const auto &list = array_of(*atoms_of, "a bond's stereo atoms");
          if (list.size() < 2) {
            refuse("gives a bond fewer than two stereo atoms");
          }
          int_of(list[0], "a bond's stereo atom");
          int_of(list[1], "a bond's stereo atom");
        } else if (stereo[i] != "either") {
          rdkit_throws();
        }
      }
    }

    // readStereoGroups. getAtomWithIdx checks the atoms. RDKit 2024.09 reads
    // no bonds for a group, so the groups hold no bonds for readQueries'
    // replaceBond to leave them pointing at (as RDKit 2026.03's can).
    if (const auto *groups = member(mol, "stereoGroups")) {
      for (const auto &group_value :
           array_of(*groups, "a molecule's stereo groups")) {
        const auto &group = object_of(group_value, "a stereo group");
        const auto *type = member(group, "type");
        const auto *ids = member(group, "atoms");
        if (!type || !ids) {
          rdkit_throws();
        }
        auto typ = MI::stereoGrouplookup.find(
            string_of(*type, "a stereo group's type"));
        if (typ == MI::stereoGrouplookup.end()) {
          rdkit_throws();
        }
        if (const auto *id = member(group, "id");
            id && typ->second != RDKit::StereoGroupType::STEREO_ABSOLUTE) {
          uint_of(*id, "a stereo group's id");
        }
        // RDKit 2024.09 keeps an atom listed twice, which dangles once the
        // atom is removed (check_stereo_group_members in mol_checks.cpp).
        std::set<unsigned> seen;
        for (const auto &atom : array_of(*ids, "a stereo group's atoms")) {
          if (!seen.insert(uint_of(atom, "a stereo group's atom")).second) {
            refuse(
                "gives a stereo group an atom twice, which RDKit 2024.09 "
                "keeps, leaving it to dangle once removed");
          }
        }
      }
    }

    // readSubstanceGroups. RDKit checks the atoms, bonds, and parent atoms
    // as it sets them, but not the brackets, the bonds of the cstates, or
    // the atoms of the attachment points, which the crate trusts as
    // SubstanceGroup's own setters keep them (substance_group.cpp).
    if (const auto *sgroups = member(mol, "substanceGroups")) {
      for (const auto &sgroup_value :
           array_of(*sgroups, "a molecule's substance groups")) {
        check_substance_group(object_of(sgroup_value, "a substance group"),
                              num_atoms, num_bonds);
      }
    }

    if (d_params.parseConformers) {
      if (const auto *confs = member(mol, "conformers")) {
        for (const auto &conf : array_of(*confs, "a molecule's conformers")) {
          check_conformer(object_of(conf, "a conformer"));
        }
      }
    }

    // parseProperties checks the type of each value it reads.
    if (d_params.parseProperties) {
      if (const auto *props = member(mol, "properties")) {
        object_of(*props, "a molecule's properties");
      }
    }

    // The extensions, in order.
    if (const auto *extensions = member(mol, "extensions")) {
      for (const auto &extension_value :
           array_of(*extensions, "a molecule's extensions")) {
        const auto &extension = object_of(extension_value, "an extension");
        const auto *name_value = member(extension, "name");
        if (!name_value) {
          rdkit_throws();
        }
        // RDKit compares the name as a C string.
        auto name = string_of(*name_value, "an extension's name");
        int current;
        if (name == "rdkitRepresentation") {
          current = MI::currentRDKitRepresentationVersion;
        } else if (name == "partialCharges") {
          current = MI::currentChargeRepresentationVersion;
        } else if (name == "rdkitQueries") {
          current = MI::currentQueryRepresentationVersion;
        } else {
          continue;
        }
        const auto *version = member(extension, "formatVersion");
        if (!version) {
          rdkit_throws();
        }
        // RDKit ignores, with a warning, an extension of a later version.
        if (int_of(*version, "an extension's formatVersion") > current) {
          continue;
        }
        if (name == "rdkitRepresentation") {
          check_representation(extension);
        } else if (name == "partialCharges") {
          check_partial_charges(extension);
        } else {
          check_queries(extension, num_atoms, num_bonds);
        }
      }
    }
  }

  // As readSubstanceGroups, for one group.
  void check_substance_group(const bj::object &sgroup, std::size_t num_atoms,
                             std::size_t num_bonds) {
    const auto *props = member(sgroup, "properties");
    if (!props) {
      rdkit_throws();
    }
    const auto &properties =
        object_of(*props, "a substance group's properties");
    const auto *type = member(properties, "TYPE");
    if (!type) {
      rdkit_throws();
    }
    string_of(*type, "a substance group's TYPE");
    for (const char *key : {"atoms", "bonds", "parentAtoms"}) {
      if (const auto *ids = member(sgroup, key)) {
        for (const auto &id :
             array_of(*ids, std::string("a substance group's ") + key)) {
          uint_of(id, std::string("a substance group's ") + key);
        }
      }
    }
    if (const auto *brackets = member(sgroup, "brackets")) {
      for (const auto &bracket :
           array_of(*brackets, "a substance group's brackets")) {
        std::size_t idx = 0;
        for (const auto &point_value :
             array_of(bracket, "a substance group's bracket")) {
          const auto &point =
              array_of(point_value, "a substance group's bracket point");
          if (point.size() != 3) {
            rdkit_throws();
          }
          for (const auto &x : point) {
            number_of(x, "a bracket point's coordinate");
          }
          // RDKit writes each point into a Bracket, an array of three,
          // past its end for a fourth.
          if (idx++ >= 3) {
            refuse(
                "gives a substance group a bracket of more than three "
                "points");
          }
        }
      }
    }
    if (const auto *cstates = member(sgroup, "cstates")) {
      for (const auto &cstate_value :
           array_of(*cstates, "a substance group's cstates")) {
        const auto &cstate = object_of(cstate_value, "a cstate");
        if (uint_of(present(member(cstate, "bond"), "a cstate's bond"),
                    "a cstate's bond") >= num_bonds) {
          refuse(
              "gives a substance group a cstate for a bond the "
              "molecule lacks");
        }
        if (const auto *vector = member(cstate, "vector")) {
          const auto &point = array_of(*vector, "a cstate's vector");
          if (point.size() != 3) {
            rdkit_throws();
          }
          for (const auto &x : point) {
            number_of(x, "a cstate vector's coordinate");
          }
        }
      }
    }
    if (const auto *points = member(sgroup, "attachPoints")) {
      for (const auto &point_value :
           array_of(*points, "a substance group's attachment points")) {
        const auto &point = object_of(point_value, "an attachment point");
        if (uint_of(
                present(member(point, "aIdx"), "an attachment point's aIdx"),
                "an attachment point's aIdx") >= num_atoms) {
          refuse(
              "gives a substance group an attachment point at an "
              "atom the molecule lacks");
        }
        if (const auto *lv_idx = member(point, "lvIdx")) {
          // As RDKit stores it, in AttachPoint's int.
          auto lv =
              static_cast<int>(uint_of(*lv_idx, "an attachment point's lvIdx"));
          if (lv != -1 &&
              (lv < 0 || static_cast<std::size_t>(lv) >= num_atoms)) {
            refuse(
                "gives a substance group an attachment point that "
                "leaves an atom the molecule lacks");
          }
        }
        if (const auto *id = member(point, "id")) {
          string_of(*id, "an attachment point's id");
        }
      }
    }
  }

  // As readConformer. Coordinates for more atoms than the molecule has
  // grow the conformer, which ROMol::addConformer refuses.
  void check_conformer(const bj::object &conf) {
    const auto *dim_value = member(conf, "dim");
    if (!dim_value) {
      rdkit_throws();
    }
    int dim = int_of(*dim_value, "a conformer's dim");
    if (dim != 2 && dim != 3) {
      rdkit_throws();
    }
    const auto *coords = member(conf, "coords");
    if (!coords) {
      rdkit_throws();
    }
    for (const auto &point_value : array_of(*coords, "a conformer's coords")) {
      const auto &point = array_of(point_value, "a conformer's coordinate");
      if (point.size() != static_cast<std::size_t>(dim)) {
        rdkit_throws();
      }
      for (const auto &x : point) {
        double d = number_of(x, "a conformer's coordinate");
        if (!(d >= -kMaxCoordinate && d <= kMaxCoordinate)) {
          refuse("gives a conformer a coordinate beyond float's range");
        }
      }
    }
  }

  // As readRDKitRepresentation. getAtomWithIdx and getBondWithIdx check the
  // aromatic atoms and bonds, CIP ranks, and CIP codes' atoms, which RDKit
  // checks are ints, and getBondBetweenAtoms a ring's atoms.
  void check_representation(const bj::object &rep) {
    for (const char *key : {"aromaticAtoms", "aromaticBonds", "cipRanks"}) {
      if (const auto *list = member(rep, key)) {
        array_of(*list, std::string("an extension's ") + key);
      }
    }
    if (const auto *codes = member(rep, "cipCodes")) {
      for (const auto &code : array_of(*codes, "an extension's cipCodes")) {
        if (!code.is_array()) {
          rdkit_throws();
        }
        // RDKit reads val[0] and val[1] without checking that they exist.
        const auto &pair = array_of_at_least(code, 2, "a CIP code");
        int_of(pair[0], "a CIP code's atom");
        string_of(pair[1], "a CIP code");
      }
    }
    if (const auto *rings = member(rep, "atomRings")) {
      for (const auto &ring : array_of(*rings, "an extension's atomRings")) {
        if (!ring.is_array()) {
          rdkit_throws();
        }
        // RDKit reads an empty ring's element sz - 1, past the end.
        if (ring.get_array().empty()) {
          refuse("gives an empty atom ring");
        }
        for (const auto &atom : ring.get_array()) {
          int_of(atom, "an atom ring's atom");
        }
      }
    }
  }

  // As readPartialCharges, which checks the values' count and types.
  void check_partial_charges(const bj::object &rep) {
    if (const auto *values = member(rep, "values")) {
      array_of(*values, "an extension's partial charges");
    }
  }

  // As readQueries. It checks that there are no more queries than atoms or
  // bonds.
  void check_queries(const bj::object &rep, std::size_t num_atoms,
                     std::size_t num_bonds) {
    for (const char *key : {"atomQueries", "bondQueries"}) {
      const auto *queries = member(rep, key);
      if (!queries) {
        continue;
      }
      bool atoms = key == std::string("atomQueries");
      std::size_t idx = 0;
      for (const auto &query :
           array_of(*queries, std::string("an extension's ") + key)) {
        if (!query.is_object()) {
          rdkit_throws();
        }
        if (!member(query.get_object(), "tag")) {
          continue;
        }
        if (idx++ >= (atoms ? num_atoms : num_bonds)) {
          rdkit_throws();
        }
        if (atoms) {
          check_atom_query(query.get_object());
        } else {
          check_bond_query(query.get_object(), false);
        }
      }
    }
  }

  // readQuery's tag; RDKit throws for a query without one.
  static int query_tag(const bj::object &query) {
    const auto *tag = member(query, "tag");
    if (!tag) {
      rdkit_throws();
    }
    return int_of(*tag, "a query's tag");
  }

  static void optional_int(const bj::object &query, const char *key) {
    if (const auto *value = member(query, key)) {
      int_of(*value, std::string("a query's ") + key);
    }
  }

  static void required_int(const bj::object &query, const char *key) {
    int_of(present(member(query, key), std::string("a query's ") + key),
           std::string("a query's ") + key);
  }

  // As readBaseQuery, which checks that the query has a description.
  static void check_base_query(const bj::object &query, int tag) {
    if (!member(query, "descr")) {
      rdkit_throws();
    }
    switch (tag) {
      case MolPickler::QUERY_AND:
      case MolPickler::QUERY_OR:
      case MolPickler::QUERY_XOR:
      case MolPickler::QUERY_NULL:
        break;
      case MolPickler::QUERY_EQUALS:
      case MolPickler::QUERY_GREATER:
      case MolPickler::QUERY_GREATEREQUAL:
      case MolPickler::QUERY_LESS:
      case MolPickler::QUERY_LESSEQUAL:
        required_int(query, "val");
        optional_int(query, "tol");
        break;
      case MolPickler::QUERY_RANGE:
        required_int(query, "lower");
        required_int(query, "upper");
        optional_int(query, "tol");
        optional_int(query, "ends");
        break;
      case MolPickler::QUERY_SET:
        if (const auto *set = member(query, "set")) {
          for (const auto &value : array_of(*set, "a query's set")) {
            int_of(value, "a query's set member");
          }
        }
        break;
      default:
        rdkit_throws();
    }
  }

  // finishQuery's description, which it reads without checking that it is
  // there for a recursive or ring query, and its type and negation.
  static std::string check_finish(const bj::object &query) {
    auto descr =
        string_of(present(member(query, "descr"), "a query's description"),
                  "a query's description");
    if (const auto *type = member(query, "type")) {
      string_of(*type, "a query's type");
    }
    if (const auto *negated = member(query, "negated")) {
      bool_of(*negated, "a query's negation");
    }
    return descr;
  }

  // The values RDKit reads for the query's tag (readBaseQuery), each within
  // U7's bound.
  static void check_values(const bj::object &query, int tag) {
    auto bounded = [&](const char *key) {
      const auto *value = member(query, key);
      if (value && value->is_int64() &&
          (value->get_int64() > kMaxQueryValue ||
           value->get_int64() < -kMaxQueryValue)) {
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
  static void check_children(const bj::object &query, int tag, F &&check) {
    const auto *children = member(query, "children");
    const bj::array *list =
        children ? &array_of(*children, "a query's children") : nullptr;
    std::size_t n = list ? list->size() : 0;
    if (is_boolean_tag(tag) && n < 2) {
      refuse("gives an AND, OR, or XOR query fewer than two children");
    }
    if (!is_boolean_tag(tag) && n > 0) {
      refuse("gives children to a query that is not AND, OR, or XOR");
    }
    if (list) {
      for (const auto &child : *list) {
        // readQuery looks for the child's tag without checking that it is
        // an object.
        check(object_of(child, "a query's child"));
      }
    }
  }

  // As readQuery for atoms.
  void check_atom_query(const bj::object &query) {
    int tag = query_tag(query);
    if (tag == MolPickler::QUERY_RECURSIVE) {
      // RDKit reads the molecule before the description.
      const auto *mol = member(query, "subquery");
      if (!mol) {
        rdkit_throws();
      }
      check_mol(*mol);
      // RecursiveStructureQuery's matcher reads the molecule's first atom.
      if (mol->get_object().at("atoms").get_array().empty()) {
        refuse("gives a recursive query a molecule without atoms");
      }
    } else if (tag == MolPickler::QUERY_ATOMRING) {
      required_int(query, "val");
      optional_int(query, "tol");
    } else {
      check_base_query(query, tag);
    }
    auto descr = check_finish(query);
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
      known = is_atom_value_query(base) && tag == prefixed;
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
    // RDKit 2024.09's finalizeQueryFromDescription casts an AtomRingSize
    // query to EqualityQuery to read its value, which a RangeQuery is not.
    if (base == "AtomRingSize" && tag == MolPickler::QUERY_RANGE) {
      refuse(
          "gives an atom a range query of ring size, which RDKit 2024.09 "
          "reads as an equality query");
    }
    check_values(query, tag);
    check_children(query, tag,
                   [&](const bj::object &child) { check_atom_query(child); });
  }

  // As readQuery for bonds. A child of an AND or OR query is not a null
  // query: SmartsWrite.cpp casts each child that is not AND or OR to
  // EqualityQuery. XOR bond queries are not read: QueryBond.cpp's matcher
  // casts the other bond's query to EqualityQuery whatever its class.
  // finalizeQueryFromDescription casts a BondRingSize query to
  // EqualityQuery, which every class allowed it here is.
  void check_bond_query(const bj::object &query, bool in_and_or) {
    int tag = query_tag(query);
    check_base_query(query, tag);
    auto descr = check_finish(query);
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
    check_values(query, tag);
    check_children(query, tag, [&](const bj::object &child) {
      check_bond_query(child, true);
    });
  }

  const MI::JSONParseParameters &d_params;
  const bj::object *d_bond_defaults = nullptr;
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

// RDKit 2024.09's readSubstanceGroups leaves an attachment point's lvIdx
// uninitialized when the point has none, as RDKit's writer leaves out -1
// (still so in RDKit 2026.09.1). The checked document gets 4294967295
// there, which RDKit reads as -1 into AttachPoint's int.
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

// RDKit 2024.09 parses the document it is given with rapidjson 1.1.0's
// default flags (JSONDataToMols: no comments, NaN, or infinities), which it
// does not install. This parses it with boost::json's default options,
// which accept no more (they also refuse invalid UTF-8 and nesting past 32
// levels, which bounds RDKit's recursion through the document, and keep the
// last of duplicate keys), checks it, and hands RDKit the checked document,
// serialized again, so RDKit parses exactly what was checked (U2), whatever
// its rapidjson makes of the original text. rapidjson reads boost::json's
// serialization as the same values: strings escaped as JSON, integers as
// digits, doubles with an exponent (rj_int), and an infinity, which
// boost::json parses from a number past double's range, as 1e99999, past
// rapidjson's, which it refuses to parse, as RDKit then reports.
std::unique_ptr<MolList> json_to_mols(rust::Str json_block,
                                      const JSONParseParametersData &params) {
  auto p = parse_params(params);
  std::vector<boost::shared_ptr<RDKit::ROMol>> mols;
  try {
    bj::value doc =
        bj::parse(std::string_view(json_block.data(), json_block.size()));
    DocumentChecker check(doc, p);
    fill_leaving_atoms(doc);
    std::string checked = bj::serialize(doc);
    // rapidjson counts a string's characters and an array's or object's
    // members in 32 bits (SizeType), each fewer than the text's characters.
    if (checked.size() >= std::numeric_limits<std::uint32_t>::max()) {
      refuse("is too long for RDKit's parser, which counts in 32 bits");
    }
    mols = MI::JSONDataToMols(checked, p);
  } catch (const boost::system::system_error &e) {
    // boost::json's errors, for text that is not JSON, are parse errors
    // here (E4).
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
