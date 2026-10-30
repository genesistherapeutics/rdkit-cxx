// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_checks.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <GraphMol/MonomerInfo.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDValue.h>
#include <boost/shared_array.hpp>
#include <boost/tokenizer.hpp>

namespace rdkit_cxx {

// The matrices are the only properties RDKit stores as shared_arrays, under
// names a caller can prefix, so they are found by type rather than by name.
void clear_matrix_caches(const RDKit::ROMol &mol) {
  std::vector<std::string> keys;
  for (const auto &pair : mol.getDict()) {
    if (RDKit::rdvalue_is<boost::shared_array<double>>(pair.val) ||
        RDKit::rdvalue_is<boost::shared_array<int>>(pair.val)) {
      keys.push_back(pair.key);
    }
  }
  for (const auto &key : keys) {
    mol.clearProp(key);
  }
}

// The names RDKit caches descriptors under: Crippen.cpp, MolSurf.cpp (the
// TPSA names end in its includeSandP flag), ConnectivityDescriptors.cpp,
// PMI.cpp, and the principal axes and moments in MolTransforms.cpp, which
// PMI.cpp and PBF.cpp read.
void clear_descriptor_caches(const RDKit::ROMol &mol) {
  static const char *const names[] = {"_connectivityHKDeltas",
                                      "_connectivityNVals",
                                      "_crippenLogP",
                                      "_crippenMR",
                                      "_crippenLogPContribs",
                                      "_crippenMRContribs",
                                      "_labuteASA",
                                      "_labuteAtomContribs",
                                      "_labuteAtomHContrib",
                                      "_tpsa-0",
                                      "_tpsa-1",
                                      "_tpsaAtomContribs-0",
                                      "_tpsaAtomContribs-1",
                                      "_PMI1",
                                      "_PMI2",
                                      "_PMI3",
                                      "_PMI1_mass",
                                      "_PMI2_mass",
                                      "_PMI3_mass",
                                      "_PMI1_cov",
                                      "_PMI2_cov",
                                      "_PMI3_cov",
                                      "_PMI1_mass_cov",
                                      "_PMI2_mass_cov",
                                      "_PMI3_mass_cov",
                                      "_principalAxes",
                                      "_principalMoments",
                                      "_principalAxes_noH",
                                      "_principalMoments_noH",
                                      "_principalAxes_cov",
                                      "_principalMoments_cov",
                                      "_principalAxes_noH_cov",
                                      "_principalMoments_noH_cov"};
  for (const char *name : names) {
    if (mol.hasProp(name)) {
      mol.clearProp(name);
    }
  }
}

void check_matrix_size(const RDKit::ROMol &mol) {
  if (mol.getNumAtoms() > 46340) {
    throw ValueErrorException(
        "RDKit computes atom matrices in int, which overflows above 46340 "
        "atoms");
  }
}

void check_formal_charge_sum(const RDKit::ROMol &mol, std::int64_t headroom) {
  std::int64_t total = headroom;
  for (const auto *atom : mol.atoms()) {
    total += std::abs(atom->getFormalCharge());
  }
  if (total > std::numeric_limits<int>::max()) {
    throw ValueErrorException(
        "RDKit sums formal charges in int, which overflows for this molecule");
  }
}

namespace {

void check_monomer_info(const RDKit::Atom &atom, const std::string &which) {
  using RDKit::AtomMonomerInfo;
  const auto *info = atom.getMonomerInfo();
  if (!info) {
    return;
  }
  int type = info->getMonomerType();
  if (type < AtomMonomerInfo::UNKNOWN || type > AtomMonomerInfo::OTHER) {
    throw ValueErrorException(which + std::to_string(atom.getIdx()) +
                              " has unknown AtomMonomerType " +
                              std::to_string(type));
  }
  if (type == AtomMonomerInfo::PDBRESIDUE &&
      !dynamic_cast<const RDKit::AtomPDBResidueInfo *>(info)) {
    throw ValueErrorException(
        which + std::to_string(atom.getIdx()) +
        " has type PDBRESIDUE, which only AtomPDBResidueInfo can have");
  }
}

// MolPickler casts each byte it reads for an enum without checking it
// (static_cast in _addAtomFromPickle, _addBondFromPickle, _depickleStereo),
// and the Rust enums, which mirror RDKit's, have no variant for a value
// past the last one: their conversion panics. RDKit's enums count from 0.
void check_enum(int value, int last, const std::string &what) {
  if (value < 0 || value > last) {
    throw ValueErrorException(what + " has unknown value " +
                              std::to_string(value));
  }
}

void check_enum_values(const RDKit::ROMol &mol, const std::string &which) {
  for (const auto *atom : mol.atoms()) {
    const auto at = which + "atom " + std::to_string(atom->getIdx());
    check_enum(atom->getChiralTag(), RDKit::Atom::CHI_OCTAHEDRAL,
               at + "'s ChiralType");
    check_enum(atom->getHybridization(), RDKit::Atom::OTHER,
               at + "'s HybridizationType");
  }
  for (const auto *bond : mol.bonds()) {
    const auto at = which + "bond " + std::to_string(bond->getIdx());
    check_enum(bond->getBondType(), RDKit::Bond::ZERO, at + "'s BondType");
    check_enum(bond->getBondDir(), RDKit::Bond::UNKNOWN, at + "'s BondDir");
    check_enum(bond->getStereo(), RDKit::Bond::STEREOATROPCCW,
               at + "'s BondStereo");
  }
  for (const auto &group : mol.getStereoGroups()) {
    check_enum(static_cast<int>(group.getGroupType()),
               static_cast<int>(RDKit::StereoGroupType::STEREO_AND),
               which + "stereo group's StereoGroupType");
  }
}

// MolPickler pickles a recursive query's molecule as a whole molecule, and
// reads it back the same way. Calls visit(m, nested) for the molecule and
// each such molecule, and visit_query(q) for each query of their atoms and
// each query under one, each once. Queries can nest without bound, so the
// walk keeps its own stack.
template <class Visit, class VisitQuery>
void walk_pickled(const RDKit::ROMol &mol, Visit visit,
                  VisitQuery visit_query) {
  using AtomQuery = RDKit::Atom::QUERYATOM_QUERY;
  std::vector<const RDKit::ROMol *> mols{&mol};
  std::vector<const AtomQuery *> queries;
  std::unordered_set<const void *> seen{&mol};
  while (!mols.empty()) {
    const auto *current = mols.back();
    mols.pop_back();
    visit(*current, current != &mol);
    for (const auto *atom : current->atoms()) {
      if (atom->hasQuery() && seen.insert(atom->getQuery()).second) {
        queries.push_back(atom->getQuery());
      }
    }
    while (!queries.empty()) {
      const auto *query = queries.back();
      queries.pop_back();
      visit_query(*query);
      if (const auto *recursive =
              dynamic_cast<const RDKit::RecursiveStructureQuery *>(query)) {
        const auto *inner = recursive->getQueryMol();
        if (inner && seen.insert(inner).second) {
          mols.push_back(inner);
        }
      }
      for (auto child = query->beginChildren(); child != query->endChildren();
           ++child) {
        if (seen.insert(child->get()).second) {
          queries.push_back(child->get());
        }
      }
    }
  }
}

// The largest value MolPickler's one-byte fields hold.
constexpr std::size_t kPickleByteMax = 255;

void check_pickle_byte(std::size_t value, const std::string &what) {
  if (value > kPickleByteMax) {
    throw ValueErrorException(what + " is " + std::to_string(value) +
                              ", more than a molecule pickle holds (255)");
  }
}

template <class Query>
void check_query_children(const Query &query, const std::string &what) {
  check_pickle_byte(
      static_cast<std::size_t>(query.endChildren() - query.beginChildren()),
      what + "'s number of children");
}

// A bond query's tree, which no walk above visits, with its own stack.
void check_bond_query(const RDKit::Bond::QUERYBOND_QUERY *root,
                      const std::string &what) {
  std::vector<const RDKit::Bond::QUERYBOND_QUERY *> queries{root};
  std::unordered_set<const void *> seen{root};
  while (!queries.empty()) {
    const auto *query = queries.back();
    queries.pop_back();
    check_query_children(*query, what);
    for (auto child = query->beginChildren(); child != query->endChildren();
         ++child) {
      if (seen.insert(child->get()).second) {
        queries.push_back(child->get());
      }
    }
  }
}

template <class T>
std::size_t num_distinct(const std::vector<T> &values) {
  return std::unordered_set<T>(values.begin(), values.end()).size();
}

// Whether some part of at most max_part_bonds bonds can hold all of the
// substance group's distinct atoms, parent atoms, and bonds.
bool fits_part(const RDKit::SubstanceGroup &sgroup, unsigned max_part_bonds) {
  auto atoms = sgroup.getAtoms();
  const auto &parents = sgroup.getParentAtoms();
  atoms.insert(atoms.end(), parents.begin(), parents.end());
  return num_distinct(atoms) <= std::size_t{max_part_bonds} + 1 &&
         num_distinct(sgroup.getBonds()) <= max_part_bonds;
}

// The fields _pickle<unsigned char> writes in one byte, which MolPickler
// uses for a molecule of at most 255 atoms whatever its bond count (RDKit
// 2026.09 also counts bonds, rdkit/rdkit#9563): every atom index fits, but
// a bond index or a list's size may not. For parts (check_picklable_parts),
// whatever the molecule's size, with bonds named as a part of at most
// max_part_bonds bonds may number them.
void check_pickle_bytes(const RDKit::ROMol &mol, const std::string &where,
                        bool parts = false, unsigned max_part_bonds = 0) {
  if (!parts && mol.getNumAtoms() > kPickleByteMax) {
    return;
  }
  // A part's bonds are numbered from 0, so past 255 only in a part of more
  // than 256 bonds.
  const bool any_bond = parts && max_part_bonds > kPickleByteMax + 1;
  auto check_bond = [&](std::size_t bond, const std::string &what) {
    if (any_bond) {
      throw ValueErrorException(
          what + " " + std::to_string(bond) +
          " may be numbered past 255 in a part of more than 256 bonds, more "
          "than a molecule pickle holds");
    }
    if (!parts) {
      check_pickle_byte(bond, what);
    }
  };
  const auto &sgroups = RDKit::getSubstanceGroups(mol);
  for (std::size_t i = 0; i < sgroups.size(); ++i) {
    const auto &sgroup = sgroups[i];
    // copyMolSubset keeps a substance group only in a part that holds all
    // its atoms, parent atoms, and bonds, and a part of k bonds has at most
    // k + 1 atoms, so a group no part has room for is never pickled. (It
    // keeps a stereo group that shares one atom with the part, so every
    // stereo group counts.)
    if (parts && !fits_part(sgroup, max_part_bonds)) {
      continue;
    }
    const auto at = where + "substance group " + std::to_string(i) + "'s ";
    check_pickle_byte(sgroup.getAtoms().size(), at + "number of atoms");
    check_pickle_byte(sgroup.getParentAtoms().size(),
                      at + "number of parent atoms");
    check_pickle_byte(sgroup.getBonds().size(), at + "number of bonds");
    for (auto bond : sgroup.getBonds()) {
      check_bond(bond, at + "bond");
    }
    check_pickle_byte(sgroup.getBrackets().size(), at + "number of brackets");
    check_pickle_byte(sgroup.getCStates().size(), at + "number of CStates");
    for (const auto &cstate : sgroup.getCStates()) {
      check_bond(cstate.bondIdx, at + "CState's bond");
    }
    check_pickle_byte(sgroup.getAttachPoints().size(),
                      at + "number of attachment points");
  }
  const auto &groups = mol.getStereoGroups();
  check_pickle_byte(groups.size(), where + "number of stereo groups");
  for (std::size_t i = 0; i < groups.size(); ++i) {
    const auto at = where + "stereo group " + std::to_string(i) + "'s ";
    check_pickle_byte(groups[i].getAtoms().size(), at + "number of atoms");
    check_pickle_byte(groups[i].getBonds().size(), at + "number of bonds");
    for (const auto *bond : groups[i].getBonds()) {
      check_bond(bond->getIdx(), at + "bond");
    }
  }
}

// streamWriteProps counts the properties it writes in a std::uint16_t for
// the molecule, each conformer, atom and bond (with the property flags, or
// for a recursive query's molecule the default flags), and each substance
// group (always): past 65535 the count wraps, and the reader takes the rest
// of the properties for the pickle's next fields. Counts every property,
// whether or not the flags write it, which is at least as many.
constexpr std::size_t kPickleCountMax = 65535;

void check_prop_count(const RDKit::RDProps &props, const std::string &what) {
  auto count = props.getDict().size();
  if (count > kPickleCountMax) {
    throw ValueErrorException(
        what + " has " + std::to_string(count) +
        " properties, more than a molecule pickle counts (65535)");
  }
}

void check_prop_counts(const RDKit::ROMol &mol, const std::string &where) {
  check_prop_count(mol, where.substr(0, where.size() - 3));
  for (auto conf = mol.beginConformers(); conf != mol.endConformers(); ++conf) {
    check_prop_count(**conf,
                     where + "conformer " + std::to_string((*conf)->getId()));
  }
  for (const auto *atom : mol.atoms()) {
    check_prop_count(*atom, where + "atom " + std::to_string(atom->getIdx()));
  }
  for (const auto *bond : mol.bonds()) {
    check_prop_count(*bond, where + "bond " + std::to_string(bond->getIdx()));
  }
  const auto &sgroups = RDKit::getSubstanceGroups(mol);
  for (std::size_t i = 0; i < sgroups.size(); ++i) {
    check_prop_count(sgroups[i],
                     where + "substance group " + std::to_string(i));
  }
}

}  // namespace

void check_pickled_mol(const RDKit::ROMol &mol) {
  walk_pickled(
      mol,
      [](const RDKit::ROMol &current, bool nested) {
        check_enum_values(current,
                          nested ? "a recursive query's molecule's " : "");
        const std::string which =
            nested ? "atom of a recursive query's molecule, " : "atom ";
        for (const auto *atom : current.atoms()) {
          check_monomer_info(*atom, which);
        }
      },
      [](const auto &) {});
}

namespace {

// check_picklable, and for parts, check_picklable_parts: the molecule's own
// fields as parts' are, its recursive queries' molecules, which parts copy
// whole, as themselves.
void check_picklable(const RDKit::ROMol &mol, bool parts,
                     unsigned max_part_bonds) {
  walk_pickled(
      mol,
      [&](const RDKit::ROMol &current, bool nested) {
        const std::string where =
            nested ? "a recursive query's molecule's " : "the molecule's ";
        check_pickle_bytes(current, where, parts && !nested, max_part_bonds);
        check_prop_counts(current, where);
        for (const auto *bond : current.bonds()) {
          if (bond->hasQuery()) {
            check_bond_query(
                bond->getQuery(),
                where + "bond " + std::to_string(bond->getIdx()) + "'s query");
          }
        }
      },
      [](const auto &query) { check_query_children(query, "an atom query"); });
}

}  // namespace

void check_picklable(const RDKit::ROMol &mol) {
  check_picklable(mol, false, 0);
}

void check_picklable_parts(const RDKit::ROMol &mol, unsigned max_part_bonds) {
  check_picklable(mol, true, max_part_bonds);
}

namespace {

// Whether std::stod reads the token, and RDKit's conversion of the result to
// unsigned int (which truncates toward zero) is defined.
bool converts_to_unsigned(const std::string &token) {
  const char *begin = token.c_str();
  char *end = nullptr;
  errno = 0;
  double value = std::strtod(begin, &end);
  if (end == begin || errno == ERANGE) {
    return false;
  }
  return std::isfinite(value) && value > -1.0 && value < 4294967296.0;
}

}  // namespace

// Reads each value as RWMol.cpp does.
void check_bond_end_points(const RDKit::ROMol &mol) {
  std::string value;
  for (const auto bond : mol.bonds()) {
    if (!bond->getPropIfPresent(RDKit::common_properties::_MolFileBondEndPts,
                                value)) {
      continue;
    }
    if (value.empty()) {
      throw ValueErrorException("bond " + std::to_string(bond->getIdx()) +
                                " has an empty _MolFileBondEndPts");
    }
    if (value.front() != '(' || value.back() != ')') {
      continue;
    }
    std::string inner = value.substr(1, value.length() - 2);
    boost::char_separator<char> sep(" ");
    boost::tokenizer<boost::char_separator<char>> tokens(inner, sep);
    if (tokens.begin() == tokens.end()) {
      throw ValueErrorException("bond " + std::to_string(bond->getIdx()) +
                                " has a _MolFileBondEndPts with no count");
    }
    for (const auto &token : tokens) {
      if (!converts_to_unsigned(token)) {
        throw ValueErrorException(
            "bond " + std::to_string(bond->getIdx()) +
            " has a _MolFileBondEndPts value that is not an unsigned int: " +
            token);
      }
    }
  }
}

// Follows the chains as removeSubstanceGroupsReferencing (SubstanceGroup.cpp)
// does: a later group with the same index replaces an earlier one.
void check_substance_group_parents(const RDKit::ROMol &mol) {
  const auto &sgroups = RDKit::getSubstanceGroups(mol);
  std::map<unsigned, std::size_t> by_index;
  bool parents = false;
  for (std::size_t i = 0; i < sgroups.size(); ++i) {
    unsigned index;
    if (sgroups[i].getPropIfPresent("index", index)) {
      by_index[index] = i;
    }
    unsigned parent;
    parents = sgroups[i].getPropIfPresent("PARENT", parent) || parents;
  }
  if (!parents) {
    return;
  }
  for (std::size_t i = 0; i < sgroups.size(); ++i) {
    std::vector<char> seen(sgroups.size(), 0);
    std::size_t at = i;
    unsigned length = 0;
    unsigned parent;
    while (sgroups[at].getPropIfPresent("PARENT", parent)) {
      auto next = by_index.find(parent);
      if (next == by_index.end()) {
        break;
      }
      seen[at] = 1;
      at = next->second;
      if (seen[at]) {
        throw ValueErrorException(
            "substance group " + std::to_string(i) +
            "'s PARENT chain is a cycle, which RDKit follows forever when it "
            "removes an atom or bond");
      }
      if (++length > MAX_PARENT_CHAIN) {
        throw ValueErrorException("substance group " + std::to_string(i) +
                                  "'s PARENT chain is longer than " +
                                  std::to_string(MAX_PARENT_CHAIN) +
                                  ", which RDKit follows recursively");
      }
    }
  }
}

void check_removals(const RDKit::ROMol &mol) {
  check_bond_end_points(mol);
  check_substance_group_parents(mol);
}

// A depth-first search of the graph whose edges lead from each group to
// every group with its PARENT's index, with its own stack: a path in it is
// a chain RDKit may follow.
void check_removals(const std::vector<const RDKit::ROMol *> &mols) {
  std::vector<const RDKit::SubstanceGroup *> groups;
  for (const auto *mol : mols) {
    check_removals(*mol);
    for (const auto &sgroup : RDKit::getSubstanceGroups(*mol)) {
      groups.push_back(&sgroup);
    }
  }
  std::map<unsigned, std::vector<std::size_t>> by_index;
  std::vector<const std::vector<std::size_t> *> next(groups.size(), nullptr);
  std::vector<unsigned> parents(groups.size());
  std::vector<char> has_parent(groups.size(), 0);
  for (std::size_t i = 0; i < groups.size(); ++i) {
    unsigned index;
    if (groups[i]->getPropIfPresent("index", index)) {
      by_index[index].push_back(i);
    }
    has_parent[i] = groups[i]->getPropIfPresent("PARENT", parents[i]);
  }
  for (std::size_t i = 0; i < groups.size(); ++i) {
    if (has_parent[i]) {
      auto found = by_index.find(parents[i]);
      if (found != by_index.end()) {
        next[i] = &found->second;
      }
    }
  }
  // 0: not reached, 1: on the stack, 2: done, with its longest chain.
  std::vector<char> state(groups.size(), 0);
  std::vector<unsigned> longest(groups.size(), 0);
  std::vector<std::pair<std::size_t, std::size_t>> stack;
  for (std::size_t start = 0; start < groups.size(); ++start) {
    if (state[start]) {
      continue;
    }
    state[start] = 1;
    stack.emplace_back(start, 0);
    while (!stack.empty()) {
      const auto node = stack.back().first;
      const auto pos = stack.back().second;
      if (next[node] && pos < next[node]->size()) {
        stack.back().second = pos + 1;
        const auto child = (*next[node])[pos];
        if (state[child] == 1) {
          throw ValueErrorException(
              "the substance groups' PARENT chains can form a cycle once "
              "RDKit combines the molecules, which it would follow forever "
              "when it removes an atom or bond");
        }
        if (state[child] == 0) {
          state[child] = 1;
          stack.emplace_back(child, 0);
        }
        continue;
      }
      unsigned length = 0;
      if (next[node]) {
        for (auto child : *next[node]) {
          length = std::max(length, longest[child] + 1);
        }
      }
      if (length > MAX_PARENT_CHAIN) {
        throw ValueErrorException(
            "the substance groups' PARENT chains can grow longer than " +
            std::to_string(MAX_PARENT_CHAIN) +
            " once RDKit combines the molecules, which RDKit follows "
            "recursively");
      }
      longest[node] = length;
      state[node] = 2;
      stack.pop_back();
    }
  }
}

namespace {

// The group, holding `bonds` instead of its own, with its ids.
RDKit::StereoGroup with_bonds(const RDKit::StereoGroup &group,
                              std::vector<RDKit::Bond *> bonds) {
  RDKit::StereoGroup res(group.getGroupType(), group.getAtoms(), bonds,
                         group.getReadId());
  res.setWriteId(group.getWriteId());
  return res;
}

}  // namespace

// replaceBond deletes the old bond, so the groups' pointers to it are found
// first, by position, and never read once it is gone.
void replace_bond(RDKit::RWMol &mol, unsigned idx, RDKit::Bond *bond,
                  bool preserve_props, bool keep_s_groups) {
  const RDKit::Bond *old_bond =
      idx < mol.getNumBonds() ? mol.getBondWithIdx(idx) : nullptr;
  std::vector<std::vector<bool>> held;
  bool any = false;
  for (const auto &group : mol.getStereoGroups()) {
    std::vector<bool> positions;
    for (const auto *member : group.getBonds()) {
      positions.push_back(member == old_bond);
      any = any || member == old_bond;
    }
    held.push_back(std::move(positions));
  }
  mol.replaceBond(idx, bond, preserve_props, keep_s_groups);
  if (!any) {
    return;
  }
  auto *new_bond = mol.getBondWithIdx(idx);
  const auto &groups = mol.getStereoGroups();
  std::vector<RDKit::StereoGroup> updated;
  updated.reserve(groups.size());
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const auto &members = groups[g].getBonds();
    std::vector<RDKit::Bond *> bonds;
    bonds.reserve(members.size());
    for (std::size_t i = 0; i < members.size(); ++i) {
      bonds.push_back(held[g][i] ? new_bond : members[i]);
    }
    updated.push_back(with_bonds(groups[g], std::move(bonds)));
  }
  mol.setStereoGroups(std::move(updated));
}

bool has_stereo_group_bonds(const RDKit::ROMol &mol) {
  for (const auto &group : mol.getStereoGroups()) {
    if (!group.getBonds().empty()) {
      return true;
    }
  }
  return false;
}

void clear_stereo_group_bonds(RDKit::ROMol &mol) {
  if (!has_stereo_group_bonds(mol)) {
    return;
  }
  std::vector<RDKit::StereoGroup> stripped;
  for (const auto &group : mol.getStereoGroups()) {
    stripped.push_back(with_bonds(group, {}));
  }
  mol.setStereoGroups(std::move(stripped));
}

StereoGroupBondGuard::StereoGroupBondGuard(RDKit::ROMol &mol) : d_mol(mol) {
  if (!has_stereo_group_bonds(mol)) {
    return;
  }
  for (const auto &group : mol.getStereoGroups()) {
    std::vector<unsigned> indices;
    for (const auto *bond : group.getBonds()) {
      indices.push_back(bond->getIdx());
    }
    d_bonds.push_back(std::move(indices));
  }
  clear_stereo_group_bonds(mol);
}

// The guarded calls keep the groups and the bonds' indices; were they to
// change, the groups would be left without bonds, which is safe.
StereoGroupBondGuard::~StereoGroupBondGuard() {
  const auto &groups = d_mol.getStereoGroups();
  if (d_bonds.empty() || groups.size() != d_bonds.size()) {
    return;
  }
  std::vector<RDKit::StereoGroup> restored;
  restored.reserve(groups.size());
  for (std::size_t g = 0; g < groups.size(); ++g) {
    std::vector<RDKit::Bond *> bonds;
    for (auto idx : d_bonds[g]) {
      if (idx < d_mol.getNumBonds()) {
        bonds.push_back(d_mol.getBondWithIdx(idx));
      }
    }
    restored.push_back(with_bonds(groups[g], std::move(bonds)));
  }
  d_mol.setStereoGroups(std::move(restored));
}

}  // namespace rdkit_cxx
