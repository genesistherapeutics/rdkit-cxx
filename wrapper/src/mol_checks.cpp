// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_checks.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDValue.h>
#include <RingDecomposerLib.h>
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

namespace {

// As findRingFamilies (FindRings.cpp) chooses the bonds it gives
// RingDecomposerLib.
bool ring_family_bond(const RDKit::Bond &bond, bool include_dative_bonds,
                      bool include_hydrogen_bonds, bool include_zero_bonds) {
  auto type = bond.getBondType();
  return (include_zero_bonds || type != RDKit::Bond::ZERO) &&
         (include_dative_bonds || !RDKit::isDative(type)) &&
         (include_hydrogen_bonds || type != RDKit::Bond::HYDROGEN);
}

void check_relevant_cycle_count(double count) {
  // Also refuses RDL_INVALID_RC_COUNT and NaN.
  if (!(count <= std::numeric_limits<int>::max())) {
    throw ValueErrorException(
        "more relevant cycles than RDKit can enumerate: it converts their "
        "count to size_t unchecked and indexes rings in int");
  }
}

// Whether the graph RDKit gives RingDecomposerLib (every bond but zero-order
// ones, and dative and hydrogen bonds only if included) has at most
// MAX_INDEPENDENT_CYCLES independent cycles.
bool few_independent_cycles(const RDKit::ROMol &mol, bool include_dative_bonds,
                            bool include_hydrogen_bonds,
                            bool include_zero_bonds) {
  // There are no more than the bonds, which settles most molecules.
  if (mol.getNumBonds() <= MAX_INDEPENDENT_CYCLES) {
    return true;
  }
  std::vector<unsigned> parent(mol.getNumAtoms());
  std::iota(parent.begin(), parent.end(), 0u);
  auto root = [&parent](unsigned idx) {
    while (parent[idx] != idx) {
      idx = parent[idx] = parent[parent[idx]];
    }
    return idx;
  };
  // A bond that joins atoms already connected closes an independent cycle.
  unsigned cycles = 0;
  for (const auto *bond : mol.bonds()) {
    if (!ring_family_bond(*bond, include_dative_bonds, include_hydrogen_bonds,
                          include_zero_bonds)) {
      continue;
    }
    auto begin = root(bond->getBeginAtomIdx());
    auto end = root(bond->getEndAtomIdx());
    if (begin == end) {
      ++cycles;
    } else {
      parent[begin] = end;
    }
  }
  return cycles <= MAX_INDEPENDENT_CYCLES;
}

}  // namespace

void check_relevant_cycles(const RDKit::ROMol &mol, bool include_dative_bonds,
                           bool include_hydrogen_bonds, bool recalc_sssr) {
  const auto *info = mol.getRingInfo();
  // symmetrizeSSSR resets an initialized ring info's families only with
  // recalcSSSR, and finds new ones only when there are none.
  if (info->areRingFamiliesInitialized() &&
      (!recalc_sssr || !info->isInitialized())) {
    check_relevant_cycle_count(RDL_getNofRC(info->dp_urfData.get()));
    return;
  }
  if (few_independent_cycles(mol, include_dative_bonds, include_hydrogen_bonds,
                             false)) {
    return;
  }
  // As findRingFamilies, which gives an empty graph one atom.
  RDL_graph *graph = RDL_initNewGraph(std::max(mol.getNumAtoms(), 1u));
  for (const auto *bond : mol.bonds()) {
    if (ring_family_bond(*bond, include_dative_bonds, include_hydrogen_bonds,
                         false)) {
      RDL_addUEdge(graph, bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }
  }
  // The data owns the graph once calculated.
  RDL_data *data = RDL_calculate(graph);
  if (data == nullptr) {
    RDL_deleteGraph(graph);
    return;
  }
  double count = RDL_getNofRC(data);
  RDL_deleteData(data);
  check_relevant_cycle_count(count);
}

void check_sanitize_relevant_cycles(const RDKit::ROMol &mol,
                                    unsigned sanitize_ops) {
  using RDKit::MolOps::SanitizeFlags;
  // The steps before SANITIZE_SYMMRINGS, which find no relevant cycles.
  constexpr unsigned before = SanitizeFlags::SANITIZE_CLEANUP |
                              SanitizeFlags::SANITIZE_CLEANUP_ORGANOMETALLICS |
                              SanitizeFlags::SANITIZE_PROPERTIES |
                              SanitizeFlags::SANITIZE_KEKULIZE;
  if ((sanitize_ops & ~before) == 0 ||
      few_independent_cycles(mol, true, true, true)) {
    return;
  }
  // sanitizeMol clears the ring info first, so it finds the families of the
  // graph those steps leave.
  RDKit::RWMol copy(mol);
  unsigned failed = 0;
  try {
    RDKit::MolOps::sanitizeMol(copy, failed, sanitize_ops & before);
  } catch (...) {
    // sanitizeMol fails the same way before it finds rings.
    return;
  }
  check_relevant_cycles(copy);
}

void check_smiles_relevant_cycles(const std::string &text,
                                  const RDKit::SmilesParserParams &params) {
  if (!params.sanitize && !params.removeHs) {
    return;
  }
  // As preprocessSmiles (SmilesParse.cpp) finds the SMILES before a name or
  // CXSMILES extensions. Replacements can write ring closures.
  if (params.replacements == nullptr || params.replacements->empty()) {
    std::size_t end = text.size();
    if (params.parseName || params.allowCXSMILES) {
      auto sidx = text.find_first_of(" \t");
      if (sidx != std::string::npos && sidx != 0) {
        end = sidx;
      }
    }
    // A SMILES's bonds but its ring closures form a spanning forest, and
    // each closure is written as two labels of digits or a '%' outside
    // brackets, so half their characters bound its independent cycles.
    std::size_t label_chars = 0;
    bool in_bracket = false;
    for (std::size_t i = 0; i < end; ++i) {
      char ch = text[i];
      if (in_bracket) {
        in_bracket = ch != ']';
      } else if (ch == '[') {
        in_bracket = true;
      } else if (ch == '%' || (ch >= '0' && ch <= '9')) {
        ++label_chars;
      }
    }
    if (label_chars / 2 <= MAX_INDEPENDENT_CYCLES) {
      return;
    }
  }
  auto unsanitized = params;
  unsanitized.sanitize = false;
  unsanitized.removeHs = false;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(text, unsanitized));
  // Text that does not parse fails the same way when sanitized.
  if (mol) {
    check_sanitize_relevant_cycles(*mol,
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
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

using AtomQuery = RDKit::Atom::QUERYATOM_QUERY;
using BondQuery = RDKit::Bond::QUERYBOND_QUERY;

// MolPickler pickles a recursive query's molecule as a whole molecule, and
// reads it back the same way, so the checks of a pickled molecule cover those
// molecules too. Their queries can nest without bound, so the walk keeps its
// own stack. It calls visit(mol, nested) on the molecule (nested false) and
// on each recursive query's molecule, and visit_query(query) on each node of
// their atoms' queries, each once.
template <class Visit, class VisitQuery>
void walk_pickled(const RDKit::ROMol &mol, Visit visit,
                  VisitQuery visit_query) {
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
      [](const AtomQuery &) {});
}

namespace {

// The largest count MolPickler writes in one byte.
constexpr std::size_t MAX_ONE_BYTE = 255;

void check_one_byte(std::size_t count, const std::string &what) {
  if (count > MAX_ONE_BYTE) {
    throw ValueErrorException(
        what + " (" + std::to_string(count) +
        ") cannot be pickled: MolPickler writes the count in one byte, and "
        "cannot read back the pickle it writes past 255");
  }
}

// pickleQuery writes a node's number of children in one byte, in a molecule
// of any size.
template <class Query>
void check_children(const Query &query, const std::string &where) {
  check_one_byte(query.endChildren() - query.beginChildren(),
                 where + "query's number of children");
}

// Bond queries hold no recursive queries, so only their children are counted.
void check_bond_queries(const RDKit::ROMol &mol, const std::string &where) {
  std::vector<const BondQuery *> queries;
  std::unordered_set<const BondQuery *> seen;
  for (const auto *bond : mol.bonds()) {
    if (bond->hasQuery() && seen.insert(bond->getQuery()).second) {
      queries.push_back(bond->getQuery());
    }
  }
  while (!queries.empty()) {
    const auto *query = queries.back();
    queries.pop_back();
    check_children(*query, where + "bond ");
    for (auto child = query->beginChildren(); child != query->endChildren();
         ++child) {
      if (seen.insert(child->get()).second) {
        queries.push_back(child->get());
      }
    }
  }
}

// MolPickler pickles a molecule with one-byte indices and counts (its
// _pickle<unsigned char>) unless it has more than 255 atoms or bonds.
bool pickles_one_byte(const RDKit::ROMol &mol) {
  return mol.getNumAtoms() <= MAX_ONE_BYTE && mol.getNumBonds() <= MAX_ONE_BYTE;
}

template <class T>
std::size_t num_distinct(const std::vector<T> &values) {
  return std::unordered_set<T>(values.begin(), values.end()).size();
}

// The counts _pickleSubstanceGroup and _pickleStereo write with one-byte
// indices: for the substance groups that keep(atoms, bonds) says are
// pickled, given their numbers of distinct atoms and bonds, and for every
// stereo group, whose pickled lists hold at most atom_room atoms and
// bond_room bonds.
template <class Keep>
void check_one_byte_groups(const RDKit::ROMol &mol, const std::string &where,
                           Keep keep, std::size_t atom_room = SIZE_MAX,
                           std::size_t bond_room = SIZE_MAX) {
  const auto &sgroups = RDKit::getSubstanceGroups(mol);
  for (std::size_t i = 0; i < sgroups.size(); ++i) {
    const auto &sgroup = sgroups[i];
    auto atoms = sgroup.getAtoms();
    const auto &parents = sgroup.getParentAtoms();
    atoms.insert(atoms.end(), parents.begin(), parents.end());
    if (!keep(num_distinct(atoms), num_distinct(sgroup.getBonds()))) {
      continue;
    }
    const auto at = where + "substance group " + std::to_string(i) + "'s ";
    check_one_byte(sgroup.getAtoms().size(), at + "number of atoms");
    check_one_byte(parents.size(), at + "number of parent atoms");
    check_one_byte(sgroup.getBonds().size(), at + "number of bonds");
    check_one_byte(sgroup.getBrackets().size(), at + "number of brackets");
    check_one_byte(sgroup.getCStates().size(), at + "number of CStates");
    check_one_byte(sgroup.getAttachPoints().size(),
                   at + "number of attachment points");
  }
  const auto &groups = mol.getStereoGroups();
  check_one_byte(groups.size(), where + "number of stereo groups");
  for (const auto &group : groups) {
    // A stereo group holds each atom and bond once.
    check_one_byte(std::min(group.getAtoms().size(), atom_room),
                   where + "stereo group's number of atoms");
    check_one_byte(std::min(group.getBonds().size(), bond_room),
                   where + "stereo group's number of bonds");
  }
}

// streamWriteProps (RDGeneral/StreamOps.h), as MolPickler calls it, counts
// the properties it writes in a std::uint16_t: for the molecule, each
// conformer, atom, and bond (with the property flags, or the default flags
// for a recursive query's molecule), and each substance group (always).
// Past 65535 the count wraps, and the reader takes the remaining property
// bytes for the pickle's next fields: with 65536 properties every read
// fails, at whatever tag those bytes make ("Bond index out of range"), and
// 65535 user properties beside a computed one are silently dropped. This
// counts every property, whether or not the flags write it, which is at
// least as many. (As rdkit-2024.09's wrapper.)
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
  check_prop_count(mol, where + "molecule");
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

// What check_picklable checks of one molecule, apart from its atom queries.
void check_picklable_mol(const RDKit::ROMol &mol, bool nested) {
  const std::string where = nested ? "a recursive query's molecule's " : "";
  if (pickles_one_byte(mol)) {
    check_one_byte_groups(mol, where,
                          [](std::size_t, std::size_t) { return true; });
  }
  check_bond_queries(mol, where);
  check_prop_counts(mol, where);
}

void check_atom_query(const AtomQuery &query) {
  check_children(query, "an atom ");
}

}  // namespace

void check_picklable(const RDKit::ROMol &mol) {
  walk_pickled(mol, check_picklable_mol, check_atom_query);
}

// A part of k bonds is a connected subgraph, so it has at most k + 1 atoms.
void check_picklable_parts(const RDKit::ROMol &mol, unsigned max_part_bonds) {
  walk_pickled(
      mol,
      [max_part_bonds](const RDKit::ROMol &current, bool nested) {
        if (nested) {
          check_picklable_mol(current, nested);
          return;
        }
        check_one_byte_groups(
            current, "",
            [max_part_bonds](std::size_t atoms, std::size_t bonds) {
              return atoms <= std::size_t{max_part_bonds} + 1 &&
                     bonds <= max_part_bonds;
            },
            std::size_t{max_part_bonds} + 1, max_part_bonds);
        check_bond_queries(current, "");
        check_prop_counts(current, "");
      },
      check_atom_query);
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
