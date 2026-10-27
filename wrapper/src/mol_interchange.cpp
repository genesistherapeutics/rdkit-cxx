// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_interchange.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cstdint>
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

// The wrapper reads the JSON with boost::json, header-only, its sources
// compiled once, here. (RDKit 2024.03 reads it with rapidjson, whose headers
// conda's RDKit does not ship.)
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

// Where RDKit throws for what it reads, the checker stops: RDKit reads no
// further, and its own error reaches the caller.
struct RDKitThrows {};

[[noreturn]] void rdkit_throws() { throw RDKitThrows{}; }

// RDKit reads the JSON's strings with rapidjson's GetString, as C strings,
// so they end at a NUL.
std::string text(const bj::string &s) { return std::string(s.c_str()); }

const bj::value *member(const bj::value &v, std::string_view key) {
  const auto *obj = v.if_object();
  if (!obj) {
    return nullptr;
  }
  auto it = obj->find(key);
  return it == obj->end() ? nullptr : &it->value();
}

bj::value *member(bj::value &v, std::string_view key) {
  return const_cast<bj::value *>(member(std::as_const(v), key));
}

// rapidjson's type tests, for the document as rapidjson reads boost::json's
// serialization of it: boost::json keeps an integer in int64 when it fits
// (otherwise in uint64, and past that in a double) and writes every double
// with an exponent, so rapidjson reads integers as integers and doubles as
// doubles. IsInt holds for integers in int's range, IsUint for those in
// unsigned's, IsDouble for doubles only, and IsNumber for all three.
bool rj_is_int(const bj::value &v) {
  return v.is_int64() && v.as_int64() >= INT_MIN && v.as_int64() <= INT_MAX;
}

bool rj_is_uint(const bj::value &v) {
  return v.is_int64() && v.as_int64() >= 0 && v.as_int64() <= UINT_MAX;
}

int rj_int(const bj::value &v) { return static_cast<int>(v.as_int64()); }

unsigned rj_uint(const bj::value &v) {
  return static_cast<unsigned>(v.as_int64());
}

// A number boost::json reads past double's range is infinity, which it
// writes as 1e99999; rapidjson's default parse refuses that
// (kParseErrorNumberTooBig) and leaves RDKit a null document, for which
// DocToMols throws.
bool has_infinity(const bj::value &v) {
  if (v.is_double()) {
    double d = v.as_double();
    return d > DBL_MAX || d < -DBL_MAX;
  }
  if (const auto *list = v.if_array()) {
    for (const auto &e : *list) {
      if (has_infinity(e)) {
        return true;
      }
    }
  } else if (const auto *obj = v.if_object()) {
    for (const auto &e : *obj) {
      if (has_infinity(e.value())) {
        return true;
      }
    }
  }
  return false;
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
// document, what it trusts, following it read by read. RDKit 2024.03 reads
// the document with rapidjson, and defines no RAPIDJSON_ASSERT, so
// rapidjson's asserts are compiled out of its release build: a getter
// (GetInt, GetUint, GetDouble, GetFloat, GetBool, GetString, GetArray,
// GetObject) on a value of another type, a member lookup (operator[],
// HasMember, FindMember) on a value that is not an object, operator[] for a
// member that is missing, and an index past an array's end are undefined.
// Whatever RDKit reads so without testing it first must be there, with the
// type its getter needs, and is refused with ValueErrorException otherwise;
// where RDKit tests what it reads and throws, the check stops (RDKitThrows)
// and leaves the error to RDKit, which reads no further. The checker stops
// only where RDKit is sure to throw, and goes on where it cannot tell.
class DocumentChecker {
 public:
  DocumentChecker(bj::value &doc, const MI::JSONParseParameters &params)
      : d_params(params) {
    try {
      check_document(doc);
    } catch (const RDKitThrows &) {
    }
  }

 private:
  // As DocToMols.
  void check_document(bj::value &doc) {
    if (has_infinity(doc) || !doc.is_object()) {
      rdkit_throws();
    }
    if (const auto *header = member(doc, "commonchem")) {
      const auto *version = member(*header, "version");
      if (!header->is_object() || !version) {
        rdkit_throws();
      }
      int v = header_version(*version);
      if (v != MI::currentMolJSONVersion) {
        rdkit_throws();
      }
    } else if (const auto *header = member(doc, "rdkitjson")) {
      const auto *version = member(*header, "version");
      if (!header->is_object() || !version) {
        rdkit_throws();
      }
      int v = header_version(*version);
      if (v > MI::currentRDKitJSONVersion || v < 10) {
        rdkit_throws();
      }
    } else {
      rdkit_throws();
    }

    // RDKit finds the defaults with GetValueByPointer, which finds nothing
    // under a "defaults" that is not an object, and throws for defaults
    // that are not objects.
    bj::value atom_defaults = bj::object();
    bj::value bond_defaults = bj::object();
    if (const auto *defaults = member(doc, "defaults")) {
      if (const auto *atom = member(*defaults, "atom")) {
        if (!atom->is_object()) {
          rdkit_throws();
        }
        atom_defaults = *atom;
      }
      if (const auto *bond = member(*defaults, "bond")) {
        if (!bond->is_object()) {
          rdkit_throws();
        }
        bond_defaults = *bond;
      }
    }
    // DefaultValueCache searches the defaults for every value an atom or
    // bond leaves out, and without them searches a null value, whose
    // FindMember asserts it is an object. RDKit reads empty defaults as it
    // reads none (0, false, ""), so the document RDKit reads has them.
    bj::object both;
    both["atom"] = std::move(atom_defaults);
    both["bond"] = std::move(bond_defaults);
    auto &defaults = doc.as_object()["defaults"];
    defaults = std::move(both);
    d_atom_defaults = &defaults.as_object()["atom"];
    d_bond_defaults = &defaults.as_object()["bond"];

    if (auto *mols = member(doc, "molecules")) {
      auto *list = mols->if_array();
      if (!list) {
        rdkit_throws();
      }
      for (auto &mol : *list) {
        check_mol(mol);
      }
    }
  }

  // The header's version, which RDKit reads with GetInt.
  static int header_version(const bj::value &version) {
    if (!rj_is_int(version)) {
      refuse("gives a version that is not an int");
    }
    return rj_int(version);
  }

  // What getIntDefaultValue, getStringDefaultValue, and getBoolDefaultValue
  // read: the member, or else the default (FindMember finds the first of
  // a key, and the document RDKit reads, written again by boost::json, which
  // keeps the last, has one of each).
  static const bj::value *with_default(const bj::value &from,
                                       const bj::value &defaults,
                                       std::string_view key) {
    const auto *value = member(from, key);
    return value ? value : member(defaults, key);
  }

  // getIntDefaultValue, which throws for a value that is not an int.
  static int default_int(const bj::value &from, const bj::value &defaults,
                         std::string_view key) {
    const auto *value = with_default(from, defaults, key);
    if (!value) {
      return 0;
    }
    if (!rj_is_int(*value)) {
      rdkit_throws();
    }
    return rj_int(*value);
  }

  // getStringDefaultValue, which throws for a value that is not a string.
  static std::string default_string(const bj::value &from,
                                    const bj::value &defaults,
                                    std::string_view key) {
    const auto *value = with_default(from, defaults, key);
    if (!value) {
      return std::string();
    }
    const auto *s = value->if_string();
    if (!s) {
      rdkit_throws();
    }
    return text(*s);
  }

  // An array RDKit reads with GetArray, unchecked.
  static const bj::array &array_of(const bj::value &v,
                                   const std::string &what) {
    const auto *list = v.if_array();
    if (!list) {
      refuse("gives " + what + " that is not an array");
    }
    return *list;
  }

  static bj::array &array_of(bj::value &v, const std::string &what) {
    return const_cast<bj::array &>(array_of(std::as_const(v), what));
  }

  // An object RDKit reads members of, unchecked.
  static void need_object(const bj::value &v, const std::string &what) {
    if (!v.is_object()) {
      refuse("gives " + what + " that is not an object");
    }
  }

  // A value RDKit reads with GetInt, unchecked.
  static int int_of(const bj::value *v, const std::string &what) {
    if (!v || !rj_is_int(*v)) {
      refuse("gives " + what + " that is missing or not an int");
    }
    return rj_int(*v);
  }

  // A value RDKit reads with GetUint, unchecked.
  static unsigned uint_of(const bj::value *v, const std::string &what) {
    if (!v || !rj_is_uint(*v)) {
      refuse("gives " + what + " that is missing or not an unsigned int");
    }
    return rj_uint(*v);
  }

  // A value RDKit reads with GetDouble, unchecked.
  static double number_of(const bj::value &v, const std::string &what) {
    if (!v.is_number()) {
      refuse("gives " + what + " that is not a number");
    }
    return v.to_number<double>();
  }

  // A value RDKit reads with GetString, unchecked.
  static std::string string_of(const bj::value *v, const std::string &what) {
    const auto *s = v ? v->if_string() : nullptr;
    if (!s) {
      refuse("gives " + what + " that is missing or not a string");
    }
    return text(*s);
  }

  // A point of three coordinates, which RDKit reads with GetDouble after
  // testing the size (readSubstanceGroups).
  static void check_point(const bj::value &v, const std::string &what) {
    const auto &point = array_of(v, what);
    if (point.size() != 3) {
      rdkit_throws();
    }
    for (const auto &x : point) {
      number_of(x, "a coordinate");
    }
  }

  // As processMol, for a molecule and for a recursive query's.
  void check_mol(bj::value &mol) {
    need_object(mol, "a molecule");
    if (const auto *name = member(mol, "name")) {
      string_of(name, "a molecule a name");
    }
    const auto *atoms = member(mol, "atoms");
    const auto *bonds = member(mol, "bonds");
    if (!atoms || !bonds) {
      rdkit_throws();
    }
    const auto &atom_list = array_of(*atoms, "a molecule atoms");
    for (const auto &atom : atom_list) {
      check_atom(atom);
    }
    const std::size_t num_atoms = atom_list.size();

    // readBond, then readBondStereo for every bond once one has stereo.
    const auto &bond_list = array_of(*bonds, "a molecule bonds");
    std::set<std::pair<unsigned, unsigned>> bonded;
    std::vector<std::string> stereo;
    bool need_stereo_loop = false;
    for (const auto &bond : bond_list) {
      stereo.push_back(check_bond(bond, num_atoms, bonded));
      need_stereo_loop |= stereo.back() != "unspecified";
    }
    const std::size_t num_bonds = bond_list.size();
    if (need_stereo_loop) {
      for (std::size_t i = 0; i < num_bonds; ++i) {
        check_bond_stereo(bond_list[i], stereo[i]);
      }
    }

    if (const auto *groups = member(mol, "stereoGroups")) {
      check_stereo_groups(*groups, num_atoms);
    }
    if (auto *sgroups = member(mol, "substanceGroups")) {
      check_substance_groups(*sgroups, num_atoms, num_bonds);
    }
    if (d_params.parseConformers) {
      if (const auto *confs = member(mol, "conformers")) {
        for (const auto &conf : array_of(*confs, "a molecule conformers")) {
          check_conformer(conf, num_atoms);
        }
      }
    }
    // parseProperties reads each value only after testing its type.
    if (d_params.parseProperties) {
      if (const auto *props = member(mol, "properties")) {
        need_object(*props, "a molecule properties");
      }
    }

    // The extensions, in order.
    if (auto *extensions = member(mol, "extensions")) {
      bool rings = false;
      for (auto &extension : array_of(*extensions, "a molecule extensions")) {
        need_object(extension, "an extension");
        const auto *name_value = member(extension, "name");
        if (!name_value) {
          rdkit_throws();
        }
        auto name = string_of(name_value, "an extension a name");
        if (name == "rdkitRepresentation") {
          if (format_version(extension) <=
              MI::currentRDKitRepresentationVersion) {
            check_representation(extension, num_atoms, num_bonds, bonded,
                                 rings);
          }
        } else if (name == "partialCharges") {
          if (format_version(extension) <=
              MI::currentChargeRepresentationVersion) {
            check_partial_charges(extension, num_atoms);
          }
        } else if (name == "rdkitQueries") {
          if (format_version(extension) <=
              MI::currentQueryRepresentationVersion) {
            check_queries(extension, num_atoms, num_bonds);
          }
        }
      }
    }
  }

  // As readAtom, whose reads of the atom's values all throw for the wrong
  // type.
  void check_atom(const bj::value &atom) const {
    need_object(atom, "an atom");
    default_int(atom, *d_atom_defaults, "z");
    if (d_params.useHCounts) {
      default_int(atom, *d_atom_defaults, "impHs");
    }
    default_int(atom, *d_atom_defaults, "chg");
    default_int(atom, *d_atom_defaults, "nRad");
    default_int(atom, *d_atom_defaults, "isotope");
    auto stereo = default_string(atom, *d_atom_defaults, "stereo");
    if (MI::chilookup.find(stereo) == MI::chilookup.end()) {
      rdkit_throws();
    }
  }

  // As readBond: the bond's stereo.
  std::string check_bond(
      const bj::value &bond, std::size_t num_atoms,
      std::set<std::pair<unsigned, unsigned>> &bonded) const {
    need_object(bond, "a bond");
    // RDKit reads bond["atoms"][0] and [1] with GetInt.
    const auto *atoms = member(bond, "atoms");
    if (!atoms) {
      refuse("gives a bond no atoms");
    }
    const auto &ends = array_of(*atoms, "a bond atoms");
    if (ends.size() < 2) {
      refuse("gives a bond fewer than two atoms");
    }
    auto begin = static_cast<unsigned>(int_of(&ends[0], "a bond an atom"));
    auto end = static_cast<unsigned>(int_of(&ends[1], "a bond an atom"));
    // addBond checks the atoms.
    if (begin >= num_atoms || end >= num_atoms || begin == end ||
        !bonded.insert(std::minmax(begin, end)).second) {
      rdkit_throws();
    }
    auto order =
        static_cast<unsigned>(default_int(bond, *d_bond_defaults, "bo"));
    if (MI::bolookup.find(order) == MI::bolookup.end()) {
      rdkit_throws();
    }
    return default_string(bond, *d_bond_defaults, "stereo");
  }

  // As readBondStereo.
  static void check_bond_stereo(const bj::value &bond,
                                const std::string &stereo) {
    if (stereo == "unspecified") {
      return;
    }
    if (MI::stereoBondlookup.find(stereo) == MI::stereoBondlookup.end()) {
      rdkit_throws();
    }
    if (const auto *atoms = member(bond, "stereoAtoms")) {
      // RDKit reads stereoAtoms[0] and [1] with GetInt; setStereoAtoms
      // checks the atoms.
      const auto &list = array_of(*atoms, "a bond stereo atoms");
      if (list.size() < 2) {
        refuse("gives a bond fewer than two stereo atoms");
      }
      int_of(&list[0], "a bond a stereo atom");
      int_of(&list[1], "a bond a stereo atom");
    } else if (stereo != "either") {
      rdkit_throws();
    }
  }

  // As readStereoGroups, which reads each group's atoms and, in 2024.03,
  // ignores its bonds.
  static void check_stereo_groups(const bj::value &groups,
                                  std::size_t num_atoms) {
    for (const auto &group : array_of(groups, "stereo groups")) {
      need_object(group, "a stereo group");
      const auto *type = member(group, "type");
      const auto *atoms = member(group, "atoms");
      if (!type || !atoms) {
        rdkit_throws();
      }
      auto name = string_of(type, "a stereo group a type");
      auto it = MI::stereoGrouplookup.find(name);
      if (it == MI::stereoGrouplookup.end()) {
        rdkit_throws();
      }
      if (it->second != RDKit::StereoGroupType::STEREO_ABSOLUTE) {
        if (const auto *id = member(group, "id")) {
          uint_of(id, "a stereo group an id");
        }
      }
      for (const auto &atom : array_of(*atoms, "a stereo group atoms")) {
        // getAtomWithIdx checks the atom.
        if (uint_of(&atom, "a stereo group an atom") >= num_atoms) {
          rdkit_throws();
        }
      }
    }
  }

  // As readSubstanceGroups. RDKit checks the atoms, bonds, and parent atoms
  // as it sets them, but not the brackets, the bonds of the cstates, or the
  // atoms of the attachment points, which the crate trusts as
  // SubstanceGroup's own setters keep them (substance_group.cpp).
  static void check_substance_groups(bj::value &sgroups, std::size_t num_atoms,
                                     std::size_t num_bonds) {
    for (auto &sgroup : array_of(sgroups, "substance groups")) {
      need_object(sgroup, "a substance group");
      const auto *props = member(sgroup, "properties");
      if (!props) {
        rdkit_throws();
      }
      // RDKit tests the properties for TYPE with HasMember, which needs an
      // object, and reads TYPE with GetString.
      need_object(*props, "a substance group properties");
      const auto *type = member(*props, "TYPE");
      if (!type) {
        rdkit_throws();
      }
      string_of(type, "a substance group a TYPE");
      for (const char *key : {"atoms", "bonds", "parentAtoms"}) {
        if (const auto *ids = member(sgroup, key)) {
          for (const auto &id :
               array_of(*ids, std::string("a substance group ") + key)) {
            uint_of(&id, std::string("a substance group ") + key);
          }
        }
      }
      if (const auto *brackets = member(sgroup, "brackets")) {
        for (const auto &bracket :
             array_of(*brackets, "a substance group brackets")) {
          const auto &points = array_of(bracket, "a bracket");
          for (std::size_t i = 0; i < points.size(); ++i) {
            check_point(points[i], "a bracket point");
            // RDKit writes each point into a Bracket, an array of three,
            // past its end for a fourth.
            if (i == 3) {
              refuse(
                  "gives a substance group a bracket of more than three "
                  "points");
            }
          }
        }
      }
      if (const auto *cstates = member(sgroup, "cstates")) {
        for (const auto &cstate :
             array_of(*cstates, "a substance group cstates")) {
          need_object(cstate, "a cstate");
          auto bond = uint_of(member(cstate, "bond"), "a cstate a bond");
          if (bond >= num_bonds) {
            refuse(
                "gives a substance group a cstate for a bond the molecule "
                "lacks");
          }
          if (const auto *vector = member(cstate, "vector")) {
            check_point(*vector, "a cstate vector");
          }
        }
      }
      if (auto *points = member(sgroup, "attachPoints")) {
        for (auto &point :
             array_of(*points, "a substance group attachment points")) {
          need_object(point, "an attachment point");
          auto atom =
              uint_of(member(point, "aIdx"), "an attachment point an atom");
          if (atom >= num_atoms) {
            refuse(
                "gives a substance group an attachment point at an atom the "
                "molecule lacks");
          }
          // RDKit leaves an AttachPoint's lvIdx uninitialized when the
          // point has none, and reads it with GetUint when it has one. Its
          // writer leaves out an lvIdx of -1, so a point without one, or
          // with -1, is handed to RDKit with UINT_MAX, which RDKit stores
          // in AttachPoint's int as -1.
          auto *lv_idx = member(point, "lvIdx");
          if (!lv_idx || (lv_idx->is_int64() && lv_idx->as_int64() == -1)) {
            point.as_object()["lvIdx"] = std::int64_t{UINT_MAX};
          } else {
            // As RDKit stores it, in AttachPoint's int.
            auto lv = static_cast<int>(
                uint_of(lv_idx, "an attachment point a leaving atom"));
            if (lv != -1 &&
                (lv < 0 || static_cast<std::size_t>(lv) >= num_atoms)) {
              refuse(
                  "gives a substance group an attachment point that leaves "
                  "an atom the molecule lacks");
            }
          }
          if (const auto *id = member(point, "id")) {
            string_of(id, "an attachment point an id");
          }
        }
      }
    }
  }

  // As readConformer, which reads the coordinates with GetFloat, a cast of
  // GetDouble to float, undefined for a double past float's range.
  static void check_conformer(const bj::value &conf, std::size_t num_atoms) {
    need_object(conf, "a conformer");
    const auto *dim_value = member(conf, "dim");
    if (!dim_value) {
      rdkit_throws();
    }
    // As RDKit stores it, in a size_t.
    auto dim =
        static_cast<std::size_t>(int_of(dim_value, "a conformer a dimension"));
    if (dim != 2 && dim != 3) {
      rdkit_throws();
    }
    const auto *coords = member(conf, "coords");
    if (!coords) {
      rdkit_throws();
    }
    const auto &points = array_of(*coords, "a conformer coordinates");
    for (const auto &point : points) {
      const auto &xyz = array_of(point, "a conformer a point");
      if (xyz.size() != dim) {
        rdkit_throws();
      }
      for (const auto &x : xyz) {
        double d = number_of(x, "a conformer a coordinate");
        if (d > FLT_MAX || d < -FLT_MAX) {
          refuse("gives a conformer a coordinate beyond float's range");
        }
      }
    }
    if (points.size() != num_atoms) {
      rdkit_throws();
    }
  }

  // An extension's formatVersion, which RDKit tests for and reads with
  // GetInt.
  static int format_version(const bj::value &extension) {
    const auto *version = member(extension, "formatVersion");
    if (!version) {
      rdkit_throws();
    }
    return int_of(version, "an extension a format version");
  }

  // As readRDKitRepresentation.
  static void check_representation(
      const bj::value &rep, std::size_t num_atoms, std::size_t num_bonds,
      const std::set<std::pair<unsigned, unsigned>> &bonded, bool &rings) {
    // RDKit tests each index with IsInt, and getAtomWithIdx and
    // getBondWithIdx check it.
    auto indices = [&](const char *key, std::size_t n) {
      if (const auto *ids = member(rep, key)) {
        for (const auto &id : array_of(*ids, std::string("an ") + key)) {
          if (!rj_is_int(id) || static_cast<unsigned>(rj_int(id)) >= n) {
            rdkit_throws();
          }
        }
      }
    };
    indices("aromaticAtoms", num_atoms);
    indices("aromaticBonds", num_bonds);
    if (const auto *ranks = member(rep, "cipRanks")) {
      std::size_t i = 0;
      for (const auto &rank : array_of(*ranks, "CIP ranks")) {
        if (!rj_is_int(rank) || i++ >= num_atoms) {
          rdkit_throws();
        }
      }
    }
    if (const auto *codes = member(rep, "cipCodes")) {
      for (const auto &code : array_of(*codes, "CIP codes")) {
        // RDKit tests that the code is an array, and reads its [0] with
        // GetInt and, once getAtomWithIdx checks it, its [1] with
        // GetString.
        const auto *pair = code.if_array();
        if (!pair) {
          rdkit_throws();
        }
        if (pair->size() < 2) {
          refuse("gives a CIP code fewer than two values");
        }
        if (static_cast<unsigned>(int_of(&(*pair)[0], "a CIP code an atom")) >=
            num_atoms) {
          rdkit_throws();
        }
        string_of(&(*pair)[1], "a CIP code a label");
      }
    }
    if (const auto *ring_list = member(rep, "atomRings")) {
      // CHECK_INVARIANT: RDKit reads the rings once.
      if (rings) {
        rdkit_throws();
      }
      rings = true;
      for (const auto &ring : array_of(*ring_list, "atom rings")) {
        const auto *atoms = ring.if_array();
        if (!atoms) {
          rdkit_throws();
        }
        // RDKit reads an empty ring's element sz - 1, past the end.
        if (atoms->empty()) {
          refuse("gives an empty atom ring");
        }
        // RDKit reads each pair of neighbors with GetInt, and
        // getBondBetweenAtoms and CHECK_INVARIANT check them.
        std::size_t sz = atoms->size();
        for (std::size_t i = 0; i < sz; ++i) {
          auto a = static_cast<unsigned>(
              int_of(&(*atoms)[i], "an atom ring an atom"));
          auto b = static_cast<unsigned>(
              int_of(&(*atoms)[(i + 1) % sz], "an atom ring an atom"));
          if (a >= num_atoms || b >= num_atoms ||
              bonded.count(std::minmax(a, b)) == 0) {
            rdkit_throws();
          }
        }
      }
    }
  }

  // As readPartialCharges, which tests the count and each value's type.
  static void check_partial_charges(const bj::value &rep,
                                    std::size_t num_atoms) {
    if (const auto *values = member(rep, "values")) {
      array_of(*values, "partial charges");
    }
  }

  // As readQueries.
  void check_queries(bj::value &rep, std::size_t num_atoms,
                     std::size_t num_bonds) {
    if (auto *queries = member(rep, "atomQueries")) {
      std::size_t idx = 0;
      for (auto &query : array_of(*queries, "atom queries")) {
        if (!query.is_object()) {
          rdkit_throws();
        }
        if (!member(query, "tag")) {
          continue;
        }
        if (idx++ >= num_atoms) {
          rdkit_throws();
        }
        check_atom_query(query);
      }
    }
    if (auto *queries = member(rep, "bondQueries")) {
      std::size_t idx = 0;
      for (auto &query : array_of(*queries, "bond queries")) {
        if (!query.is_object()) {
          rdkit_throws();
        }
        if (!member(query, "tag")) {
          continue;
        }
        if (idx++ >= num_bonds) {
          rdkit_throws();
        }
        check_bond_query(query, false);
      }
    }
  }

  // The query's tag, which RDKit reads with GetInt once it tests for it.
  static int query_tag(const bj::value &query) {
    return int_of(member(query, "tag"), "a query a tag");
  }

  // The values readBaseQuery reads with GetInt for the query's tag, after
  // testing for the description; it throws for a tag it does not know.
  static void check_base_query(const bj::value &query, int tag) {
    if (!member(query, "descr")) {
      rdkit_throws();
    }
    auto optional_int = [&](const char *key) {
      if (const auto *value = member(query, key)) {
        int_of(value, std::string("a query a ") + key);
      }
    };
    if (is_boolean_tag(tag) || tag == MolPickler::QUERY_NULL) {
      return;
    }
    if (is_equality_tag(tag)) {
      int_of(member(query, "val"), "a query a val");
      optional_int("tol");
    } else if (tag == MolPickler::QUERY_RANGE) {
      int_of(member(query, "lower"), "a query a lower");
      int_of(member(query, "upper"), "a query an upper");
      optional_int("tol");
      optional_int("ends");
    } else if (tag == MolPickler::QUERY_SET) {
      if (const auto *set = member(query, "set")) {
        for (const auto &value : array_of(*set, "a query a set")) {
          int_of(&value, "a query a set member");
        }
      }
    } else {
      rdkit_throws();
    }
  }

  // What finishQuery reads unchecked: the description, which the recursive
  // and ring queries never test for, with GetString, the type with
  // GetString, and negated with GetBool.
  static std::string check_finish(const bj::value &query) {
    auto descr = string_of(member(query, "descr"), "a query a description");
    if (const auto *type = member(query, "type")) {
      string_of(type, "a query a type");
    }
    if (const auto *negated = member(query, "negated")) {
      if (!negated->is_bool()) {
        refuse("gives a query a negation that is not a bool");
      }
    }
    return descr;
  }

  // The values RDKit reads for the query's tag (readBaseQuery), within U7's
  // bound; the type checks have shown each is an int.
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
  // for the queries RDKit makes without them. readQuery tests each child
  // for a tag with HasMember, which needs an object.
  template <class F>
  static void check_children(bj::value &query, int tag, F &&check) {
    auto *children = member(query, "children");
    bj::array *list =
        children ? &array_of(*children, "a query children") : nullptr;
    std::size_t n = list ? list->size() : 0;
    if (is_boolean_tag(tag) && n < 2) {
      refuse("gives an AND, OR, or XOR query fewer than two children");
    }
    if (!is_boolean_tag(tag) && n > 0) {
      refuse("gives children to a query that is not AND, OR, or XOR");
    }
    if (list) {
      for (auto &child : *list) {
        need_object(child, "a query a child");
        if (!member(child, "tag")) {
          rdkit_throws();
        }
        check(child);
      }
    }
  }

  // As readQuery for atoms.
  void check_atom_query(bj::value &query) {
    int tag = query_tag(query);
    if (tag == MolPickler::QUERY_RECURSIVE) {
      // RDKit reads the molecule before the description.
      auto *mol = member(query, "subquery");
      if (!mol) {
        rdkit_throws();
      }
      check_mol(*mol);
      // RecursiveStructureQuery's matcher reads the molecule's first atom.
      if (member(*mol, "atoms")->as_array().empty()) {
        refuse("gives a recursive query a molecule without atoms");
      }
    } else if (tag == MolPickler::QUERY_ATOMRING) {
      int_of(member(query, "val"), "a query a val");
      if (const auto *tol = member(query, "tol")) {
        int_of(tol, "a query a tol");
      }
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
      // 2024.03's finalizeQueryFromDescription casts an AtomRingSize query
      // to EqualityQuery to read its value, which a range query is not.
      known = is_atom_value_query(base) && tag == prefixed &&
              !(base == "AtomRingSize" && tag == MolPickler::QUERY_RANGE);
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
    check_values(query, tag);
    check_children(query, tag,
                   [&](bj::value &child) { check_atom_query(child); });
  }

  // As readQuery for bonds. A child of an AND or OR query is not a null
  // query: SmartsWrite.cpp casts each child that is not AND or OR to
  // EqualityQuery. XOR bond queries are not read: QueryBond.cpp's matcher
  // casts the other bond's query to EqualityQuery whatever its class.
  void check_bond_query(bj::value &query, bool in_and_or) {
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
    check_children(query, tag,
                   [&](bj::value &child) { check_bond_query(child, true); });
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

// The wrapper parses the text with boost::json's default options and checks
// that document; RDKit is handed the checked document, serialized again, so
// it parses exactly what was checked (U2), whatever rapidjson would make of
// the original text: one value per key (boost::json keeps the last of a
// key, rapidjson's FindMember finds the first), no NUL outside an escape
// (RDKit parses the text up to its first NUL), integers and doubles as
// boost::json read them, and no nesting deeper than boost::json's limit
// (rapidjson's default parse recurses, as do RDKit's reads of queries).
// rapidjson's default flags read neither comments nor NaN and Infinity,
// which boost::json's default options refuse too.
std::unique_ptr<MolList> json_to_mols(rust::Str json_block,
                                      const JSONParseParametersData &params) {
  auto p = parse_params(params);
  std::vector<boost::shared_ptr<RDKit::ROMol>> mols;
  try {
    bj::value doc =
        bj::parse(std::string_view(json_block.data(), json_block.size()));
    DocumentChecker check(doc, p);
    mols = MI::JSONDataToMols(bj::serialize(doc), p);
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
