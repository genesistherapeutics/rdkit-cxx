/*
 * Copyright (c) 2026 Genesis Molecular AI, Inc.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/query.h"

#include <algorithm>
#include <mutex>
#include <vector>

#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"

namespace rdkit_cxx {

namespace {

Queries::CompositeQueryType composite(int how) {
  switch (how) {
    case Queries::COMPOSITE_AND:
    case Queries::COMPOSITE_OR:
    case Queries::COMPOSITE_XOR:
      return static_cast<Queries::CompositeQueryType>(how);
  }
  throw ValueErrorException("unknown CompositeQueryType");
}

template <class Query, class T>
Query &as_query(T &obj, const char *what) {
  auto *query = dynamic_cast<Query *>(&obj);
  if (!query) {
    throw ValueErrorException(std::string("the ") + what + " is not a query " +
                              what);
  }
  return *query;
}

}  // namespace

std::unique_ptr<std::string> atom_describe_query(const Atom &atom) {
  return std::make_unique<std::string>(RDKit::describeQuery(&atom));
}

std::unique_ptr<std::string> atom_get_query_type(const Atom &atom) {
  // QueryAtom::getQueryType reads its query without checking for one.
  if (dynamic_cast<const RDKit::QueryAtom *>(&atom) && !atom.hasQuery()) {
    return std::make_unique<std::string>();
  }
  return std::make_unique<std::string>(atom.getQueryType());
}

bool atom_match(const Atom &atom, const Atom &other) {
  return atom.Match(&other);
}

void atom_expand_query(Atom &atom, const Atom &other, int how,
                       bool maintain_order) {
  auto &self = as_query<RDKit::QueryAtom>(atom, "atom");
  const auto &query = as_query<const RDKit::QueryAtom>(other, "atom");
  auto how_ = composite(how);
  if (query.hasQuery()) {
    // QueryAtom::expandQuery reads this query's description; RDKit 2026.09
    // refuses a null one only with a PRECONDITION (an Invariant).
    if (!self.hasQuery()) {
      throw ValueErrorException("the query atom has no query to expand");
    }
    self.expandQuery(query.getQuery()->copy(), how_, maintain_order);
  }
}

void atom_set_query(Atom &atom, const Atom &other) {
  auto &self = as_query<RDKit::QueryAtom>(atom, "atom");
  const auto &query = as_query<const RDKit::QueryAtom>(other, "atom");
  if (query.hasQuery()) {
    self.setQuery(query.getQuery()->copy());
  }
}

std::unique_ptr<std::string> bond_describe_query(const Bond &bond) {
  return std::make_unique<std::string>(RDKit::describeQuery(&bond));
}

bool bond_match(const Bond &bond, const Bond &other) {
  return bond.Match(&other);
}

void bond_expand_query(Bond &bond, const Bond &other, int how,
                       bool maintain_order) {
  auto &self = as_query<RDKit::QueryBond>(bond, "bond");
  const auto &query = as_query<const RDKit::QueryBond>(other, "bond");
  auto how_ = composite(how);
  if (query.hasQuery()) {
    // QueryBond::expandQuery (QueryBond.cpp) reads this query's description
    // without checking that there is one, in every release: a null
    // dereference for a query bond without a query.
    if (!self.hasQuery()) {
      throw ValueErrorException("the query bond has no query to expand");
    }
    self.expandQuery(query.getQuery()->copy(), how_, maintain_order);
  }
}

void bond_set_query(Bond &bond, const Bond &other) {
  auto &self = as_query<RDKit::QueryBond>(bond, "bond");
  const auto &query = as_query<const RDKit::QueryBond>(other, "bond");
  if (query.hasQuery()) {
    self.setQuery(query.getQuery()->copy());
  }
}

std::unique_ptr<Atom> query_atom_without_query_for_tests() {
  return std::make_unique<RDKit::QueryAtom>();
}

std::unique_ptr<Bond> query_bond_without_query_for_tests() {
  return std::make_unique<RDKit::QueryBond>();
}

rust::Vec<std::uint32_t> mol_get_atoms_matching_query(const RWMol &mol,
                                                      const Atom &query) {
  const auto &qa = as_query<const RDKit::QueryAtom>(query, "atom");
  // Python matches the recursive queries first with detail::MatchSubqueries,
  // whose header RDKit does not install. SubstructMatch on a one-atom query
  // does the same, and with default parameters compares atoms with Match.
  RDKit::RWMol query_mol;
  query_mol.addAtom(qa.copy(), false, true);
  RDKit::SubstructMatchParameters params;
  params.maxMatches = std::max(1u, mol.getNumAtoms());
  std::vector<std::uint32_t> indices;
  // Every match constructs a ControlCHandler (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  for (const auto &match : RDKit::SubstructMatch(mol, query_mol, params)) {
    indices.push_back(match[0].second);
  }
  std::sort(indices.begin(), indices.end());
  rust::Vec<std::uint32_t> result;
  for (auto idx : indices) {
    result.push_back(idx);
  }
  return result;
}

}  // namespace rdkit_cxx
