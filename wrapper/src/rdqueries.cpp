// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/rdqueries.h"

#include <any>
#include <cstdint>
#include <string>
#include <type_traits>

#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/rdqueries.rs.h"

namespace rdkit_cxx {

namespace {

// Each value query matches by subtracting the atom's value from the query's
// in int (Queries::queryCmp). The values these queries read are below 2^30
// in magnitude (atom fields are 8 or 16 bits; degrees and ring counts are
// far smaller), except the mass, which RDKit scales by 1000 and which stays
// below 6.6e7 for 16-bit isotopes. Keeping the query value within 2^30, after
// that scaling, keeps the difference inside int.
constexpr std::int64_t max_query_value = std::int64_t{1} << 30;

void check_query_value(std::int64_t val) {
  if (val > max_query_value || val < -max_query_value) {
    throw ValueErrorException(
        "query values must lie within 2^30 (after mass's scaling by 1000), "
        "where RDKit's int comparison cannot overflow");
  }
}

// As QAFUNC1 in Code/GraphMol/Wrap/Queries.cpp: "Less" builds RDKit's
// greater-than query, and "Greater" its less-than query, because they read
// as the atom's value compared to val.
#define RDKIT_CXX_VALUE_QUERY(name, func)                                 \
  case AtomValueQuery::name:                                              \
    switch (comparison) {                                                 \
      case QueryComparison::Equals:                                       \
        return RDKit::func(val);                                          \
      case QueryComparison::Less:                                         \
        return RDKit::func<RDKit::ATOM_GREATER_QUERY>(val, #name "Less"); \
      case QueryComparison::Greater:                                      \
        return RDKit::func<RDKit::ATOM_LESS_QUERY>(val, #name "Greater"); \
    }                                                                     \
    break;

RDKit::Atom::QUERYATOM_QUERY *make_value_query(AtomValueQuery kind,
                                               QueryComparison comparison,
                                               int val) {
  switch (kind) {
    RDKIT_CXX_VALUE_QUERY(AtomNum, makeAtomNumQuery)
    RDKIT_CXX_VALUE_QUERY(ExplicitValence, makeAtomExplicitValenceQuery)
    RDKIT_CXX_VALUE_QUERY(TotalValence, makeAtomTotalValenceQuery)
    RDKIT_CXX_VALUE_QUERY(ExplicitDegree, makeAtomExplicitDegreeQuery)
    RDKIT_CXX_VALUE_QUERY(TotalDegree, makeAtomTotalDegreeQuery)
    RDKIT_CXX_VALUE_QUERY(NonHydrogenDegree, makeAtomNonHydrogenDegreeQuery)
    RDKIT_CXX_VALUE_QUERY(HCount, makeAtomHCountQuery)
    RDKIT_CXX_VALUE_QUERY(Mass, makeAtomMassQuery)
    RDKIT_CXX_VALUE_QUERY(Isotope, makeAtomIsotopeQuery)
    RDKIT_CXX_VALUE_QUERY(FormalCharge, makeAtomFormalChargeQuery)
    RDKIT_CXX_VALUE_QUERY(Hybridization, makeAtomHybridizationQuery)
    RDKIT_CXX_VALUE_QUERY(InNRings, makeAtomInNRingsQuery)
    RDKIT_CXX_VALUE_QUERY(MinRingSize, makeAtomMinRingSizeQuery)
    RDKIT_CXX_VALUE_QUERY(RingBondCount, makeAtomRingBondCountQuery)
    RDKIT_CXX_VALUE_QUERY(NumRadicalElectrons, makeAtomNumRadicalElectronsQuery)
    RDKIT_CXX_VALUE_QUERY(NumHeteroatomNeighbors,
                          makeAtomNumHeteroatomNbrsQuery)
    RDKIT_CXX_VALUE_QUERY(NumAliphaticHeteroatomNeighbors,
                          makeAtomNumAliphaticHeteroatomNbrsQuery)
  }
  throw ValueErrorException("unknown AtomValueQuery");
}

#undef RDKIT_CXX_VALUE_QUERY

RDKit::Atom::QUERYATOM_QUERY *make_flag_query(AtomFlagQuery kind) {
  switch (kind) {
    case AtomFlagQuery::IsUnsaturated:
      return RDKit::makeAtomUnsaturatedQuery();
    case AtomFlagQuery::IsAromatic:
      return RDKit::makeAtomAromaticQuery();
    case AtomFlagQuery::IsAliphatic:
      return RDKit::makeAtomAliphaticQuery();
    case AtomFlagQuery::IsInRing:
      return RDKit::makeAtomInRingQuery();
    case AtomFlagQuery::HasChiralTag:
      return RDKit::makeAtomHasChiralTagQuery();
    case AtomFlagQuery::MissingChiralTag:
      return RDKit::makeAtomMissingChiralTagQuery();
    case AtomFlagQuery::IsBridgehead:
      return RDKit::makeAtomIsBridgeheadQuery();
    case AtomFlagQuery::AAtom:
      return RDKit::makeAAtomQuery();
    case AtomFlagQuery::AHAtom:
      return RDKit::makeAHAtomQuery();
    case AtomFlagQuery::QAtom:
      return RDKit::makeQAtomQuery();
    case AtomFlagQuery::QHAtom:
      return RDKit::makeQHAtomQuery();
    case AtomFlagQuery::XAtom:
      return RDKit::makeXAtomQuery();
    case AtomFlagQuery::XHAtom:
      return RDKit::makeXHAtomQuery();
    case AtomFlagQuery::MAtom:
      return RDKit::makeMAtomQuery();
    case AtomFlagQuery::MHAtom:
      return RDKit::makeMHAtomQuery();
  }
  throw ValueErrorException("unknown AtomFlagQuery");
}

// The query atom or bond Python builds: an empty one given the query, which
// it negates on request.
template <class Result, class Query>
std::unique_ptr<Result> with_query(Query *query, bool negate) {
  std::unique_ptr<Query> owned(query);
  if (negate) {
    owned->setNegation(true);
  }
  if constexpr (std::is_same_v<Result, Atom>) {
    auto result = std::make_unique<RDKit::QueryAtom>();
    result->setQuery(owned.release());
    return result;
  } else {
    auto result = std::make_unique<RDKit::QueryBond>();
    result->setQuery(owned.release());
    return result;
  }
}

// RDKit's HasPropWithValueQuery<int> compares with queryCmp, whose
// subtraction overflows int for property values, which the target sets, far
// enough from the query's (and negates the tolerance, which overflows for
// INT_MIN). This one compares in 64 bits. It pickles and describes itself as
// the base class.
template <class Target>
class IntPropQuery : public RDKit::HasPropWithValueQuery<const Target *, int> {
 public:
  IntPropQuery(const std::string &propname, int val, int tolerance)
      : RDKit::HasPropWithValueQuery<const Target *, int>(propname, val,
                                                          tolerance),
        d_propname(propname),
        d_val(val),
        d_tolerance(tolerance) {}

  bool Match(const Target *what) const override {
    bool res = false;
    if (what->hasProp(d_propname)) {
      try {
        std::int64_t diff =
            std::int64_t{what->template getProp<int>(d_propname)} - d_val;
        res = diff <= d_tolerance && diff >= -std::int64_t{d_tolerance};
      } catch (KeyErrorException &) {
        res = false;
      } catch (std::bad_any_cast &) {
        res = false;
      }
    }
    return this->getNegation() ? !res : res;
  }

  Queries::Query<int, const Target *, true> *copy() const override {
    auto *res = new IntPropQuery(d_propname, d_val, d_tolerance);
    res->setNegation(this->getNegation());
    res->setDescription(this->getDescription());
    return res;
  }

 private:
  std::string d_propname;
  int d_val;
  int d_tolerance;
};

}  // namespace

std::unique_ptr<Atom> value_query_atom(AtomValueQuery kind,
                                       QueryComparison comparison,
                                       std::int32_t val, bool negate) {
  check_query_value(kind == AtomValueQuery::Mass ? std::int64_t{val} * 1000
                                                 : val);
  return with_query<Atom>(make_value_query(kind, comparison, val), negate);
}

std::unique_ptr<Atom> flag_query_atom(AtomFlagQuery kind, bool negate) {
  return with_query<Atom>(make_flag_query(kind), negate);
}

std::unique_ptr<Atom> has_prop_query_atom(rust::Str propname, bool negate) {
  return with_query<Atom>(
      RDKit::makeHasPropQuery<RDKit::Atom>(std::string(propname)), negate);
}

std::unique_ptr<Atom> has_int_prop_with_value_query_atom(
    rust::Str propname, std::int32_t val, bool negate, std::int32_t tolerance) {
  return with_query<Atom>(
      new IntPropQuery<RDKit::Atom>(std::string(propname), val, tolerance),
      negate);
}

std::unique_ptr<Atom> has_bool_prop_with_value_query_atom(rust::Str propname,
                                                          bool val,
                                                          bool negate) {
  return with_query<Atom>(
      RDKit::makePropQuery<RDKit::Atom, bool>(std::string(propname), val),
      negate);
}

std::unique_ptr<Atom> has_string_prop_with_value_query_atom(rust::Str propname,
                                                            rust::Str val,
                                                            bool negate) {
  return with_query<Atom>(RDKit::makePropQuery<RDKit::Atom, std::string>(
                              std::string(propname), std::string(val)),
                          negate);
}

std::unique_ptr<Atom> has_double_prop_with_value_query_atom(rust::Str propname,
                                                            double val,
                                                            bool negate,
                                                            double tolerance) {
  return with_query<Atom>(RDKit::makePropQuery<RDKit::Atom, double>(
                              std::string(propname), val, tolerance),
                          negate);
}

// The tolerance is the allowed Tanimoto distance.
std::unique_ptr<Atom> has_bit_vect_prop_with_value_query_atom(
    rust::Str propname, const BitVect &val, bool negate, double tolerance) {
  return with_query<Atom>(
      RDKit::makePropQuery<RDKit::Atom>(std::string(propname), as_explicit(val),
                                        static_cast<float>(tolerance)),
      negate);
}

std::unique_ptr<Bond> has_prop_query_bond(rust::Str propname, bool negate) {
  return with_query<Bond>(
      RDKit::makeHasPropQuery<RDKit::Bond>(std::string(propname)), negate);
}

std::unique_ptr<Bond> has_int_prop_with_value_query_bond(
    rust::Str propname, std::int32_t val, bool negate, std::int32_t tolerance) {
  return with_query<Bond>(
      new IntPropQuery<RDKit::Bond>(std::string(propname), val, tolerance),
      negate);
}

std::unique_ptr<Bond> has_bool_prop_with_value_query_bond(rust::Str propname,
                                                          bool val,
                                                          bool negate) {
  return with_query<Bond>(
      RDKit::makePropQuery<RDKit::Bond, bool>(std::string(propname), val),
      negate);
}

std::unique_ptr<Bond> has_string_prop_with_value_query_bond(rust::Str propname,
                                                            rust::Str val,
                                                            bool negate) {
  return with_query<Bond>(RDKit::makePropQuery<RDKit::Bond, std::string>(
                              std::string(propname), std::string(val)),
                          negate);
}

std::unique_ptr<Bond> has_double_prop_with_value_query_bond(rust::Str propname,
                                                            double val,
                                                            bool negate,
                                                            double tolerance) {
  return with_query<Bond>(RDKit::makePropQuery<RDKit::Bond, double>(
                              std::string(propname), val, tolerance),
                          negate);
}

// As replaceAtomWithQueryAtomHelper in Code/GraphMol/Wrap/Queries.cpp,
// with the atom given by index, so it is the molecule's own.
void replace_atom_with_query_atom(RWMol &mol, std::uint32_t idx) {
  check_index(idx, mol.getNumAtoms());
  RDKit::QueryOps::replaceAtomWithQueryAtom(&mol, mol.getAtomWithIdx(idx));
}

}  // namespace rdkit_cxx
