// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_checks.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <GraphMol/MonomerInfo.h>
#include <GraphMol/PeriodicTable.h>
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
  for (const auto &pair : mol.getDict().getData()) {
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

// The fields _pickle<unsigned char> writes in one byte, which MolPickler
// uses for a molecule of at most 255 atoms whatever its bond count (RDKit
// 2026.09 also counts bonds, rdkit/rdkit#9563): every atom index fits, but
// a bond index or a list's size may not. For parts (check_picklable_parts),
// whatever the molecule's size, with bonds named as a part of at most
// max_part_bonds bonds may number them. A part keeps a substance group only
// if it holds all of the group's atoms, parent atoms and bonds, so a group
// with more distinct atoms or bonds than a connected part of max_part_bonds
// bonds has (max_part_bonds + 1 atoms) is in none; it keeps a stereo group
// that shares any atom with it (copySelectedStereoGroups in RDKit 2025.09's
// Subset.cpp), so every stereo group is counted.
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
  auto in_a_part = [&](const RDKit::SubstanceGroup &sgroup) {
    auto atoms = sgroup.getAtoms();
    const auto &parents = sgroup.getParentAtoms();
    atoms.insert(atoms.end(), parents.begin(), parents.end());
    const auto &bonds = sgroup.getBonds();
    return std::unordered_set<unsigned>(atoms.begin(), atoms.end()).size() <=
               std::size_t{max_part_bonds} + 1 &&
           std::unordered_set<unsigned>(bonds.begin(), bonds.end()).size() <=
               max_part_bonds;
  };
  const auto &sgroups = RDKit::getSubstanceGroups(mol);
  for (std::size_t i = 0; i < sgroups.size(); ++i) {
    const auto &sgroup = sgroups[i];
    if (parts && !in_a_part(sgroup)) {
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
// of the properties for the pickle's next fields (RDKit 2025.09.6's Python:
// "Bad pickle format: ENDPROPS tag not found" for a molecule, "Bond index out
// of range" for a substance group). Counts every property, whether or not
// the flags write it, which is at least as many.
constexpr std::size_t kPickleCountMax = 65535;

void check_prop_count(const RDKit::RDProps &props, const std::string &what) {
  auto count = props.getDict().getData().size();
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

void check_stereo_group_duplicates(const RDKit::ROMol &mol) {
  for (const auto &group : mol.getStereoGroups()) {
    std::set<const RDKit::Atom *> atoms(group.getAtoms().begin(),
                                        group.getAtoms().end());
    std::set<const RDKit::Bond *> bonds(group.getBonds().begin(),
                                        group.getBonds().end());
    if (atoms.size() != group.getAtoms().size() ||
        bonds.size() != group.getBonds().size()) {
      throw ValueErrorException(
          "a stereo group holds an atom or bond twice, which RDKit 2025.09 "
          "would keep a pointer to once it removed it");
    }
  }
}

// Bridges by Tarjan's low-link, with its own stack.
void check_kekulize_questions(const RDKit::ROMol &mol) {
  const auto n = mol.getNumAtoms();
  bool aromatic = false;
  unsigned dummies = 0;
  for (const auto *atom : mol.atoms()) {
    aromatic = aromatic || atom->getIsAromatic();
    dummies += atom->getAtomicNum() == 0;
  }
  for (const auto *bond : mol.bonds()) {
    aromatic = aromatic || bond->getIsAromatic() ||
               bond->getBondType() == RDKit::Bond::AROMATIC;
  }
  if (!aromatic || dummies < MAX_KEKULIZE_QUESTIONS) {
    return;
  }
  std::vector<unsigned> order(n, 0), low(n, 0);
  std::vector<char> bridge(mol.getNumBonds(), 0);
  unsigned counter = 0;
  // An atom, the bond it was reached by, and its bonds left to follow.
  struct Frame {
    unsigned atom;
    int via;
    std::vector<const RDKit::Bond *> bonds;
    std::size_t pos;
  };
  for (unsigned start = 0; start < n; ++start) {
    if (order[start]) {
      continue;
    }
    std::vector<Frame> stack;
    auto push = [&](unsigned atom, int via) {
      order[atom] = low[atom] = ++counter;
      Frame frame{atom, via, {}, 0};
      for (const auto *bond : mol.atomBonds(mol.getAtomWithIdx(atom))) {
        frame.bonds.push_back(bond);
      }
      stack.push_back(std::move(frame));
    };
    push(start, -1);
    while (!stack.empty()) {
      auto &top = stack.back();
      if (top.pos < top.bonds.size()) {
        const auto *bond = top.bonds[top.pos++];
        if (static_cast<int>(bond->getIdx()) == top.via) {
          continue;
        }
        const auto other = bond->getOtherAtomIdx(top.atom);
        if (order[other]) {
          low[top.atom] = std::min(low[top.atom], order[other]);
        } else {
          push(other, static_cast<int>(bond->getIdx()));
        }
        continue;
      }
      const auto atom = top.atom;
      const auto via = top.via;
      stack.pop_back();
      if (!stack.empty()) {
        auto &parent = stack.back();
        low[parent.atom] = std::min(low[parent.atom], low[atom]);
        if (low[atom] > order[parent.atom]) {
          bridge[via] = 1;
        }
      }
    }
  }
  // Union the atoms across ring bonds.
  std::vector<unsigned> root(n);
  for (unsigned i = 0; i < n; ++i) {
    root[i] = i;
  }
  auto find = [&root](unsigned i) {
    while (root[i] != i) {
      root[i] = root[root[i]];
      i = root[i];
    }
    return i;
  };
  std::vector<char> in_ring(n, 0);
  for (const auto *bond : mol.bonds()) {
    if (!bridge[bond->getIdx()]) {
      in_ring[bond->getBeginAtomIdx()] = in_ring[bond->getEndAtomIdx()] = 1;
      root[find(bond->getBeginAtomIdx())] = find(bond->getEndAtomIdx());
    }
  }
  std::map<unsigned, unsigned> dummy_count;
  std::map<unsigned, bool> has_non_dummy;
  for (const auto *atom : mol.atoms()) {
    const auto idx = atom->getIdx();
    if (!in_ring[idx]) {
      continue;
    }
    const auto r = find(idx);
    if (atom->getAtomicNum() == 0) {
      ++dummy_count[r];
    } else {
      has_non_dummy[r] = true;
    }
  }
  for (const auto &[r, count] : dummy_count) {
    if (count >= MAX_KEKULIZE_QUESTIONS && has_non_dummy[r]) {
      throw ValueErrorException(
          "a ring system holds " + std::to_string(count) +
          " dummy atoms, at least " + std::to_string(MAX_KEKULIZE_QUESTIONS) +
          ", for which RDKit 2025.09's kekulization shifts past the width of "
          "an unsigned int");
    }
  }
}

void check_removals(const RDKit::ROMol &mol) {
  check_bond_end_points(mol);
  check_substance_group_parents(mol);
  check_stereo_group_duplicates(mol);
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

namespace {

// As ProximityEntry and IsBonded in ProximityBonds.cpp.
struct ProximityAtom {
  float x, y, z, r;
  unsigned elem;
};

bool proximity_bonded(const ProximityAtom &p, const ProximityAtom &q,
                      unsigned flags) {
  if (flags & kCtdIgnoreHHContacts && p.elem == 1 && q.elem == 1) {
    return false;
  }
  constexpr double kExtDist = 0.45;
  constexpr double kMinDist2 = 0.16;
  constexpr double kMaxDist2 = 29.7025;
  double dx = (double)p.x - (double)q.x;
  double dy = (double)p.y - (double)q.y;
  double dz = (double)p.z - (double)q.z;
  double dist2 = dx * dx + dy * dy + dz * dz;
  if (dist2 > kMaxDist2 || dist2 < kMinDist2) {
    return false;
  }
  double radius = (double)p.r + (double)q.r + kExtDist;
  return dist2 <= radius * radius;
}

const RDKit::AtomPDBResidueInfo *pdb_info(const RDKit::Atom *atom) {
  auto *info =
      dynamic_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
  return info && info->getMonomerType() == RDKit::AtomMonomerInfo::PDBRESIDUE
             ? info
             : nullptr;
}

// As IsBlacklistedPair in ProximityBonds.cpp.
bool proximity_blacklisted(const RDKit::Atom *a, const RDKit::Atom *b) {
  auto *ai = pdb_info(a);
  auto *bi = pdb_info(b);
  if (!ai || !bi) {
    return false;
  }
  if (ai->getResidueNumber() == bi->getResidueNumber() &&
      ai->getResidueName() == bi->getResidueName() &&
      ai->getChainId() == bi->getChainId() &&
      ai->getInsertionCode() == bi->getInsertionCode()) {
    return false;
  }
  auto blacklisted = [](const RDKit::Atom *atom) {
    int elem = atom->getAtomicNum();
    return !((5 <= elem && elem <= 8) || (14 <= elem && elem <= 16) ||
             (32 <= elem && elem <= 34) || (51 <= elem && elem <= 52));
  };
  return blacklisted(a) || blacklisted(b) || ai->getResidueName() == "HOH" ||
         bi->getResidueName() == "HOH";
}

}  // namespace

void check_proximity_hydrogens(const RDKit::ROMol &mol, unsigned flags,
                               bool bonds_removed) {
  if (!mol.getNumConformers()) {
    return;
  }
  const auto &conf = mol.getConformer();
  auto *table = RDKit::PeriodicTable::getTable();
  const unsigned n = mol.getNumAtoms();
  std::vector<ProximityAtom> atoms(n);
  // Atoms RDKit bonds are within 5.45 of each other, so in neighboring cells
  // of that size.
  constexpr double kCell = 5.45;
  std::map<std::array<std::int64_t, 3>, std::vector<unsigned>> cells;
  auto cell = [&](const ProximityAtom &a) {
    return std::array<std::int64_t, 3>{
        static_cast<std::int64_t>(std::floor(a.x / kCell)),
        static_cast<std::int64_t>(std::floor(a.y / kCell)),
        static_cast<std::int64_t>(std::floor(a.z / kCell))};
  };
  bool any_h = false;
  for (unsigned i = 0; i < n; ++i) {
    const auto &p = conf.getAtomPos(i);
    unsigned elem = mol.getAtomWithIdx(i)->getAtomicNum();
    atoms[i] = {(float)p.x, (float)p.y, (float)p.z,
                (float)table->getRcovalent(elem), elem};
    if (!std::isfinite(atoms[i].x) || !std::isfinite(atoms[i].y) ||
        !std::isfinite(atoms[i].z)) {
      continue;  // never bonded (and refused elsewhere)
    }
    cells[cell(atoms[i])].push_back(i);
    any_h = any_h || elem == 1;
  }
  if (!any_h) {
    return;
  }
  for (unsigned i = 0; i < n; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    if (atoms[i].elem != 1 || !std::isfinite(atoms[i].x) ||
        !std::isfinite(atoms[i].y) || !std::isfinite(atoms[i].z)) {
      continue;
    }
    unsigned degree = 0;
    bool to_h = false;
    bool to_other = false;
    auto count = [&](unsigned j) {
      ++degree;
      (atoms[j].elem == 1 ? to_h : to_other) = true;
    };
    if (!bonds_removed) {
      for (const auto *nbr : mol.atomNeighbors(atom)) {
        count(nbr->getIdx());
      }
    }
    const auto c = cell(atoms[i]);
    for (std::int64_t dx = -1; dx <= 1; ++dx) {
      for (std::int64_t dy = -1; dy <= 1; ++dy) {
        for (std::int64_t dz = -1; dz <= 1; ++dz) {
          auto it = cells.find({c[0] + dx, c[1] + dy, c[2] + dz});
          if (it == cells.end()) {
            continue;
          }
          for (unsigned j : it->second) {
            if (j == i || (!bonds_removed && mol.getBondBetweenAtoms(i, j))) {
              continue;
            }
            if (proximity_bonded(atoms[i], atoms[j], flags) &&
                !proximity_blacklisted(atom, mol.getAtomWithIdx(j))) {
              count(j);
            }
          }
        }
      }
    }
    if (flags & kCtdQuickRemoveHHContacts && to_h && to_other) {
      --degree;
    }
    if (degree >= 3) {
      throw ValueErrorException(
          "hydrogen " + std::to_string(i) +
          " would be bonded to three or more atoms by proximity, which RDKit "
          "2025.09 reads past its neighbor list for when it keeps the "
          "closest");
    }
  }
}

}  // namespace rdkit_cxx
