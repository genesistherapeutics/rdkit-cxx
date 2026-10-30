// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_enumerator.h"

#include <algorithm>
#include <cstddef>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/Exceptions.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/mol_enumerator.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

namespace ME = RDKit::MolEnumerator;
using tokenizer = boost::tokenizer<boost::char_separator<char>>;

// getVariations (MolEnumerator.cpp) recurses once per variation point of an
// operation, with about 320 bytes of stack a level, so an operation with
// more points is refused: the stack U12 allows std::regex.
constexpr std::size_t MAX_VARIATION_POINTS = 4000;

// LinkNodeOp labels atom i with the isotope 1000 * (i + 1), which Atom keeps
// in 16 bits, and finds each label again with a reaction's SMARTS.
constexpr unsigned MAX_LINK_NODE_ATOM = 64;

// As opFromName in Wrap/rdMolEnumerator.cpp, where an unknown type raises
// ValueError.
std::shared_ptr<ME::MolEnumeratorOp> op_from_type(EnumeratorType type) {
  switch (type) {
    case EnumeratorType::LinkNode:
      return std::make_shared<ME::LinkNodeOp>();
    case EnumeratorType::PositionVariation:
      return std::make_shared<ME::PositionVariationOp>();
    case EnumeratorType::RepeatUnit:
      return std::make_shared<ME::RepeatUnitOp>();
  }
  throw std::invalid_argument("unrecognized operator type");
}

void check_points(std::size_t count, const std::string &what) {
  if (count > MAX_VARIATION_POINTS) {
    throw ValueErrorException(
        "RDKit enumerates each of an operation's variation points one "
        "recursion deeper: more than " +
        std::to_string(MAX_VARIATION_POINTS) + " " + what);
  }
}

// PositionVariationOp::initFromMol's variation points: each bond with
// _MolFileBondEndPts and an _MolFileBondAttach of ANY, and the atom of it
// that operator() removes, chosen as it chooses.
std::vector<unsigned> position_variation_dummies(const RDKit::ROMol &mol) {
  std::vector<unsigned> dummies;
  for (const auto bond : mol.bonds()) {
    std::string endpts;
    std::string attach;
    if (!bond->getPropIfPresent(RDKit::common_properties::_MolFileBondEndPts,
                                endpts) ||
        !bond->getPropIfPresent(RDKit::common_properties::_MolFileBondAttach,
                                attach) ||
        attach != "ANY") {
      continue;
    }
    const RDKit::Atom *atom = bond->getBeginAtom();
    if (atom->getAtomicNum() == 0) {
      atom = bond->getEndAtom();
      if (atom->getAtomicNum() == 0 && atom->hasQuery() &&
          atom->getQuery()->getDescription() == "AtomNull" &&
          bond->getBeginAtom()->hasQuery() &&
          bond->getBeginAtom()->getQuery()->getDescription() != "AtomNull") {
        atom = bond->getBeginAtom();
      }
    }
    dummies.push_back(bond->getOtherAtomIdx(atom->getIdx()));
  }
  return dummies;
}

// RepeatUnitOp::initFromMol's variation points: the SRU groups it
// enumerates. It marks each one's atoms in bitsets sized for the molecule
// without a range check, so an atom out of range is refused.
std::size_t repeat_units(const RDKit::ROMol &mol) {
  std::size_t count = 0;
  for (const auto &sg : RDKit::getSubstanceGroups(mol)) {
    std::string type;
    if (!sg.getPropIfPresent("TYPE", type) || type != "SRU") {
      continue;
    }
    std::string connect;
    sg.getPropIfPresent("CONNECT", connect);
    if (!connect.empty() && connect != "HT" && connect != "HH") {
      continue;
    }
    ++count;
    for (auto idx : sg.getAtoms()) {
      if (idx >= mol.getNumAtoms()) {
        throw ValueErrorException(
            "an SRU substance group's atom index is out of range: " +
            std::to_string(idx));
      }
    }
  }
  return count;
}

// LinkNodeOp::operator() (LinkNode.cpp) runs a reaction per link node that
// finds the node's atom and its two neighbors by their isotope labels, across
// single bonds, and repeats the atom between them; it takes the reaction's
// first product without checking that there is one (ASSERT_INVARIANT is an
// assert, compiled out in release builds), so a reaction that finds nothing
// reads past the end of an empty vector. This refuses the nodes it cannot
// show will be found:
// - a node atom whose label does not fit 16 bits (index above 64);
// - a bond to the node atom that is not single;
// - a node whose neighbors are one atom;
// - two nodes on one atom, since the first's reaction may leave no copy of
//   the atom next to both of the second's neighbors.
// Otherwise each reaction finds its atoms: earlier reactions copy only their
// own node atoms, and keep every other atom, and a bond of the node atom to a
// neighbor that an earlier reaction copied, or an earlier node's atom, goes
// to a copy with the same label, singly. User isotopes equal to a label only
// add matches. Nodes RDKit rejects while it reads them (getMolLinkNodes)
// stop the check, since RDKit then runs no reaction. removed holds the atoms
// position variation removes first, which renumber the rest.
void check_link_nodes(const RDKit::ROMol &mol,
                      const std::vector<unsigned> &removed) {
  std::string text;
  if (!mol.getPropIfPresent(RDKit::common_properties::molFileLinkNodes, text)) {
    return;
  }
  const std::set<unsigned> gone(removed.begin(), removed.end());
  const unsigned num_atoms = mol.getNumAtoms();
  struct Node {
    unsigned atom, first, second;
  };
  std::vector<Node> nodes;
  boost::char_separator<char> pipe("|");
  boost::char_separator<char> space(" ");
  for (const auto &node_text : tokenizer(text, pipe)) {
    std::vector<unsigned> data;
    try {
      for (const auto &token : tokenizer(node_text, space)) {
        data.push_back(boost::lexical_cast<unsigned>(token));
      }
    } catch (const boost::bad_lexical_cast &) {
      return;
    }
    if (data.size() < 5 || data.size() < 3 + 2 * data[2] || data[0] == 0 ||
        data[1] < data[0] || data[2] != 2 || data[3] != data[5]) {
      return;
    }
    for (unsigned i = 3; i <= 6; ++i) {
      unsigned idx = data[i] - 1;
      if (idx >= num_atoms || gone.count(idx)) {
        return;  // not in the atom map
      }
    }
    Node node{data[3] - 1, data[4] - 1, data[6] - 1};
    if (node.first == node.atom || node.second == node.atom) {
      return;  // no bond to itself, which RDKit requires
    }
    nodes.push_back(node);
  }
  std::set<unsigned> atoms;
  for (const auto &node : nodes) {
    for (auto idx : {node.atom, node.first, node.second}) {
      auto shifted = idx - static_cast<unsigned>(std::distance(
                               gone.begin(), gone.lower_bound(idx)));
      if (shifted > MAX_LINK_NODE_ATOM) {
        throw ValueErrorException(
            "RDKit labels link node atoms with isotopes that only fit for "
            "the first " +
            std::to_string(MAX_LINK_NODE_ATOM + 1) + " atoms, not atom " +
            std::to_string(idx));
      }
    }
    for (auto nbr : {node.first, node.second}) {
      // A bond position variation adds is single.
      const auto *bond = mol.getBondBetweenAtoms(node.atom, nbr);
      if (bond && bond->getBondType() != RDKit::Bond::BondType::SINGLE) {
        throw ValueErrorException(
            "RDKit's link node reaction matches only single bonds, not the "
            "bond between atoms " +
            std::to_string(node.atom) + " and " + std::to_string(nbr));
      }
    }
    if (node.first == node.second) {
      throw ValueErrorException(
          "a link node's two bonds go to one atom, which RDKit's reaction "
          "cannot match");
    }
    if (!atoms.insert(node.atom).second) {
      throw ValueErrorException("two link nodes repeat atom " +
                                std::to_string(node.atom) +
                                ", which RDKit's reactions may then not find");
    }
  }
}

}  // namespace

// The checks of an enumeration with these operations: their variation
// points, and the atoms they remove (check_removals).
void check_enumerate(const RDKit::ROMol &mol, bool position_variation,
                     bool repeat_unit, bool link_node) {
  std::vector<unsigned> dummies;
  if (position_variation) {
    dummies = position_variation_dummies(mol);
    check_points(dummies.size(), "position variation bonds");
  }
  if (repeat_unit) {
    check_points(repeat_units(mol), "SRU substance groups");
  }
  if (position_variation || repeat_unit) {
    check_removals(mol);
  }
  if (link_node) {
    std::string text;
    if (mol.getPropIfPresent(RDKit::common_properties::molFileLinkNodes,
                             text)) {
      boost::char_separator<char> pipe("|");
      tokenizer nodes(text, pipe);
      check_points(std::distance(nodes.begin(), nodes.end()), "link nodes");
    }
    check_link_nodes(mol, dummies);
  }
}

namespace {

std::unique_ptr<MolList> to_list(const RDKit::MolBundle &bundle) {
  auto list = std::make_unique<MolList>();
  for (const auto &mol : bundle.getMols()) {
    auto copy = std::make_unique<RWMol>(*mol);
    clear_matrix_caches(*copy);
    list->push(std::move(copy));
  }
  return list;
}

}  // namespace

MolEnumeratorParamsData mol_enumerator_params_default() {
  ME::MolEnumeratorParams params;
  MolEnumeratorParamsData data{};
  data.sanitize = params.sanitize;
  data.max_to_enumerate = params.maxToEnumerate;
  data.do_random = params.doRandom;
  data.random_seed = params.randomSeed;
  data.has_operation = false;
  data.operation = EnumeratorType::LinkNode;
  return data;
}

// Position variation runs first, then repeat units, then link nodes, whose
// reactions match substructures, and every match constructs a
// ControlCHandler (control_c_lock.h). Their reactions also block RDKit's
// logs with RDLog::LogStateSetter, which changes process-wide state without
// a lock, as BlockLogs does.
std::unique_ptr<MolList> mol_enumerator_enumerate(
    const RWMol &mol, std::uint32_t max_per_operation) {
  check_enumerate(mol, true, true, true);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return to_list(ME::enumerate(mol, max_per_operation));
}

std::unique_ptr<MolList> mol_enumerator_enumerate_with_params(
    const RWMol &mol, const MolEnumeratorParamsData &data) {
  ME::MolEnumeratorParams params;
  params.sanitize = data.sanitize;
  params.maxToEnumerate = data.max_to_enumerate;
  params.doRandom = data.do_random;
  params.randomSeed = data.random_seed;
  std::unique_lock<ControlCMutex> control_c;
  if (data.has_operation) {
    params.dp_operation = op_from_type(data.operation);
    auto type = data.operation;
    check_enumerate(mol, type == EnumeratorType::PositionVariation,
                    type == EnumeratorType::RepeatUnit,
                    type == EnumeratorType::LinkNode);
    if (type == EnumeratorType::LinkNode) {
      control_c = std::unique_lock<ControlCMutex>(control_c_mutex());
    }
  }
  return to_list(ME::enumerate(mol, params));
}

}  // namespace rdkit_cxx
