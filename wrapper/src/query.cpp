/*
 * Copyright (c) 2026 Genesis Molecular AI
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
  // QueryAtom::expandQuery checks its own query with a PRECONDITION; this
  // refuses it as bond_expand_query does.
  if (!self.hasQuery()) {
    throw ValueErrorException("the atom has no query to expand");
  }
  if (query.hasQuery()) {
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
  // QueryBond::expandQuery reads its own query's description without
  // checking that it has one (QueryBond.cpp, in every release through
  // 2026.09.1), and a default-constructed QueryBond has none.
  if (!self.hasQuery()) {
    throw ValueErrorException("the bond has no query to expand");
  }
  if (query.hasQuery()) {
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
