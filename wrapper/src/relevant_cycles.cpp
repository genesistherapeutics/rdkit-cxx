// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/relevant_cycles.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/FileParsers/CDXMLParser.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/MarvinParse/MarvinParser.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/StreamOps.h>
#include <RingDecomposerLib.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/tokenizer.hpp>

namespace rdkit_cxx {

namespace FP = RDKit::v2::FileParsers;
namespace FPU = RDKit::FileParserUtils;

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

}  // namespace

namespace {

// Whether the sum over the biconnected blocks of the bonds' graph of
// 2^c - 1, for a block of c independent cycles, is at most INT_MAX. The
// blocks are found as Tarjan's algorithm finds them, keeping its own stack.
bool few_block_cycles(const RDKit::ROMol &mol, bool include_dative_bonds,
                      bool include_hydrogen_bonds, bool include_zero_bonds) {
  const unsigned n = mol.getNumAtoms();
  std::vector<std::vector<std::pair<unsigned, unsigned>>> adjacent(n);
  for (const auto *bond : mol.bonds()) {
    if (ring_family_bond(*bond, include_dative_bonds, include_hydrogen_bonds,
                         include_zero_bonds)) {
      auto a = bond->getBeginAtomIdx();
      auto b = bond->getEndAtomIdx();
      adjacent[a].emplace_back(b, bond->getIdx());
      adjacent[b].emplace_back(a, bond->getIdx());
    }
  }
  constexpr unsigned kUnseen = UINT_MAX;
  std::vector<unsigned> order(n, kUnseen), low(n, 0), stamp(n, kUnseen);
  std::vector<std::pair<unsigned, unsigned>> edges;  // the open blocks'
  struct Frame {
    unsigned atom, via_bond;
    std::size_t next = 0;
  };
  std::vector<Frame> stack;
  unsigned seen = 0, blocks = 0;
  std::uint64_t total = 0;
  for (unsigned root = 0; root < n; ++root) {
    if (order[root] != kUnseen) {
      continue;
    }
    order[root] = low[root] = seen++;
    stack.push_back({root, kUnseen});
    while (!stack.empty()) {
      auto &frame = stack.back();
      const auto v = frame.atom;
      if (frame.next < adjacent[v].size()) {
        auto [w, bond] = adjacent[v][frame.next++];
        if (bond == frame.via_bond) {
          continue;
        }
        if (order[w] == kUnseen) {
          edges.emplace_back(v, w);
          order[w] = low[w] = seen++;
          stack.push_back({w, bond});
        } else if (order[w] < order[v]) {
          edges.emplace_back(v, w);
          low[v] = std::min(low[v], order[w]);
        }
        continue;
      }
      stack.pop_back();
      if (stack.empty()) {
        break;
      }
      const auto u = stack.back().atom;
      low[u] = std::min(low[u], low[v]);
      if (low[v] < order[u]) {
        continue;
      }
      // The edges from (u, v) on are a block.
      std::uint64_t num_edges = 0, num_atoms = 0;
      for (;;) {
        auto [a, b] = edges.back();
        edges.pop_back();
        ++num_edges;
        for (auto atom : {a, b}) {
          if (stamp[atom] != blocks) {
            stamp[atom] = blocks;
            ++num_atoms;
          }
        }
        if (a == u && b == v) {
          break;
        }
      }
      ++blocks;
      const auto cycles = num_edges - num_atoms + 1;
      if (cycles >= 31) {
        return false;
      }
      total += (std::uint64_t{1} << cycles) - 1;
      if (total > static_cast<std::uint64_t>(INT_MAX)) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace

bool few_relevant_cycles(const RDKit::ROMol &mol, bool include_dative_bonds,
                         bool include_hydrogen_bonds, bool include_zero_bonds) {
  // There are no more independent cycles than bonds, which settles most
  // molecules.
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
  return cycles <= MAX_INDEPENDENT_CYCLES ||
         few_block_cycles(mol, include_dative_bonds, include_hydrogen_bonds,
                          include_zero_bonds);
}

std::uint64_t independent_cycles(const RDKit::ROMol &mol) {
  std::vector<unsigned> parent(mol.getNumAtoms());
  std::iota(parent.begin(), parent.end(), 0u);
  auto root = [&parent](unsigned idx) {
    while (parent[idx] != idx) {
      idx = parent[idx] = parent[parent[idx]];
    }
    return idx;
  };
  std::uint64_t cycles = 0;
  for (const auto *bond : mol.bonds()) {
    auto begin = root(bond->getBeginAtomIdx());
    auto end = root(bond->getEndAtomIdx());
    if (begin == end) {
      ++cycles;
    } else {
      parent[begin] = end;
    }
  }
  return cycles;
}

std::uint64_t CycleCounter::root(std::uint64_t node) {
  auto it = d_parent.try_emplace(node, node).first;
  while (it->second != node) {
    auto up = d_parent.find(it->second);
    // Path halving.
    it->second = up->second;
    node = it->second;
    it = d_parent.find(node);
  }
  return node;
}

void CycleCounter::add_bond(std::uint64_t begin, std::uint64_t end) {
  auto b = root(begin);
  auto e = root(end);
  if (b == e) {
    ++d_cycles;
  } else {
    d_parent[b] = e;
  }
}

bool reaction_keeps_cycles(const RDKit::ChemicalReaction &rxn, bool exactly) {
  auto pair = [](const RDKit::Bond &bond) {
    int a = bond.getBeginAtom()->getAtomMapNum();
    int b = bond.getEndAtom()->getAtomMapNum();
    return std::make_pair(std::min(a, b), std::max(a, b));
  };
  std::set<std::pair<int, int>> reactant_bonds;
  std::set<int> reactant_atoms;
  for (const auto &templ : rxn.getReactants()) {
    for (const auto *atom : templ->atoms()) {
      if (exactly && !atom->getAtomMapNum()) {
        return false;  // RDKit deletes it
      }
      reactant_atoms.insert(atom->getAtomMapNum());
    }
    for (const auto *bond : templ->bonds()) {
      auto p = pair(*bond);
      if (p.first) {
        reactant_bonds.insert(p);
      }
    }
  }
  std::set<std::pair<int, int>> product_bonds;
  std::set<int> product_atoms;
  for (const auto &templ : rxn.getProducts()) {
    CycleCounter cycles;
    const std::uint64_t mapped = templ->getNumAtoms();
    for (const auto *atom : templ->atoms()) {
      // A mapped atom written twice is copied with its other bonds.
      if (atom->getAtomMapNum() &&
          !product_atoms.insert(atom->getAtomMapNum()).second) {
        return false;
      }
    }
    for (const auto *bond : templ->bonds()) {
      if (exactly && !ring_family_bond(*bond, false, false, false)) {
        return false;
      }
      auto p = pair(*bond);
      if (p.first) {
        if (!reactant_bonds.count(p)) {
          return false;
        }
        product_bonds.insert(p);
        continue;
      }
      int a = bond->getBeginAtom()->getAtomMapNum();
      int b = bond->getEndAtom()->getAtomMapNum();
      cycles.add_bond(a ? mapped : bond->getBeginAtomIdx(),
                      b ? mapped : bond->getEndAtomIdx());
    }
    if (cycles.cycles() > 0) {
      return false;
    }
  }
  if (exactly) {
    for (int num : reactant_atoms) {
      if (!product_atoms.count(num)) {
        return false;
      }
    }
    return product_bonds == reactant_bonds;
  }
  return true;
}

TransformCycles transform_cycles(
    const std::vector<std::shared_ptr<RDKit::ChemicalReaction>> &rxns) {
  TransformCycles result;
  for (const auto &rxn : rxns) {
    result.keeps = result.keeps && reaction_keeps_cycles(*rxn, false);
    result.keeps_exactly =
        result.keeps_exactly && reaction_keeps_cycles(*rxn, true);
  }
  return result;
}

void check_standardize_relevant_cycles(const RDKit::ROMol &mol,
                                       TransformCycles normalizations,
                                       bool custom_metals) {
  if (RDKit::MolOps::getUseLegacyRingFinding()) {
    return;
  }
  if (!normalizations.keeps) {
    throw ValueErrorException(
        "normalizations that bond atoms anew cannot be checked for the "
        "relevant cycles of what they make before RDKit perceives them");
  }
  if (independent_cycles(mol) <= MAX_INDEPENDENT_CYCLES) {
    return;
  }
  // The atoms QueryOps::isMetal (and so cleanUpOrganometallics) does not
  // take for metals; no metal pattern of MetalDisconnector's names them.
  static const std::set<int> nonmetals{0,  1,  2,  5,  6,  7,  8,  9,
                                       10, 14, 15, 16, 17, 18, 33, 34,
                                       35, 36, 52, 53, 54, 85, 86};
  bool kept = normalizations.keeps_exactly && !custom_metals;
  for (const auto *bond : mol.bonds()) {
    if (!kept) {
      break;
    }
    kept = ring_family_bond(*bond, false, false, false) &&
           nonmetals.count(bond->getBeginAtom()->getAtomicNum()) &&
           nonmetals.count(bond->getEndAtom()->getAtomicNum());
  }
  if (!kept) {
    throw ValueErrorException(
        "the standardizer refuses a molecule with more than " +
        std::to_string(MAX_INDEPENDENT_CYCLES) +
        " independent cycles and bonds to metals, zero-order, dative or "
        "hydrogen bonds, or normalizations that change its rings: the "
        "relevant cycles of what it makes cannot be counted before RDKit "
        "perceives them");
  }
  check_graph_relevant_cycles(mol);
}

void check_graph_relevant_cycles(const RDKit::ROMol &mol,
                                 bool include_dative_bonds,
                                 bool include_hydrogen_bonds) {
  if (few_relevant_cycles(mol, include_dative_bonds, include_hydrogen_bonds,
                          false) ||
      // symmetrizeSSSR's default algorithm is then the legacy one, which
      // enumerates no relevant cycles.
      RDKit::MolOps::getUseLegacyRingFinding()) {
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

void check_relevant_cycles(const RDKit::ROMol &mol, bool include_dative_bonds,
                           bool include_hydrogen_bonds, bool recalc_sssr) {
  const auto *info = mol.getRingInfo();
  // symmetrizeSSSR resets an initialized ring info's families only with
  // recalcSSSR, and finds new ones only when there are none.
  if (info->areRingFamiliesInitialized() &&
      (!recalc_sssr || !info->isInitialized())) {
    if (!RDKit::MolOps::getUseLegacyRingFinding()) {
      check_relevant_cycle_count(RDL_getNofRC(info->dp_urfData.get()));
    }
    return;
  }
  check_graph_relevant_cycles(mol, include_dative_bonds,
                              include_hydrogen_bonds);
}

void check_perceive_relevant_cycles(const RDKit::ROMol &mol,
                                    bool unless_initialized) {
  const auto *info = mol.getRingInfo();
  if (unless_initialized ? info->isInitialized() : info->isSymmSssr()) {
    return;
  }
  check_relevant_cycles(mol);
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
      few_relevant_cycles(mol, true, true, true)) {
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

void check_remove_hs_relevant_cycles(
    const RDKit::ROMol &mol, const RDKit::MolOps::RemoveHsParameters &ps,
    bool sanitize) {
  if (!sanitize || few_relevant_cycles(mol, true, true, true)) {
    return;
  }
  RDKit::RWMol copy(mol);
  try {
    RDKit::MolOps::removeHs(copy, ps, false);
  } catch (...) {
    // removeHs fails the same way before it sanitizes.
    return;
  }
  check_sanitize_relevant_cycles(copy,
                                 RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
}

void check_stereo_relevant_cycles(const RDKit::ROMol &mol) {
  check_graph_relevant_cycles(mol);
}

namespace {

using AtomQuery = RDKit::Atom::QUERYATOM_QUERY;

// Checks the graph of each recursive query's molecule under `queries`, and
// of the recursive queries under those molecules' atoms. They nest without
// bound, so the walk keeps its own stack.
void check_recursive_queries(std::vector<const AtomQuery *> queries,
                             std::unordered_set<const void *> &seen) {
  std::vector<const RDKit::ROMol *> mols;
  while (!queries.empty() || !mols.empty()) {
    if (queries.empty()) {
      const auto *current = mols.back();
      mols.pop_back();
      check_graph_relevant_cycles(*current);
      for (const auto *atom : current->atoms()) {
        if (atom->hasQuery() && seen.insert(atom->getQuery()).second) {
          queries.push_back(atom->getQuery());
        }
      }
      continue;
    }
    const auto *query = queries.back();
    queries.pop_back();
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

}  // namespace

void check_smarts_relevant_cycles(const RDKit::ROMol &mol, bool include_mol) {
  if (include_mol) {
    check_graph_relevant_cycles(mol);
  }
  std::unordered_set<const void *> seen{&mol};
  std::vector<const AtomQuery *> queries;
  for (const auto *atom : mol.atoms()) {
    if (atom->hasQuery() && seen.insert(atom->getQuery()).second) {
      queries.push_back(atom->getQuery());
    }
  }
  check_recursive_queries(std::move(queries), seen);
}

void check_atom_smarts_relevant_cycles(const RDKit::Atom &atom) {
  if (!atom.hasQuery()) {
    return;
  }
  std::unordered_set<const void *> seen{atom.getQuery()};
  check_recursive_queries({atom.getQuery()}, seen);
}

void check_parsed_relevant_cycles(const RDKit::ROMol &raw, bool sanitize,
                                  bool remove_hs) {
  if (!sanitize) {
    return;
  }
  if (remove_hs) {
    check_remove_hs_relevant_cycles(raw, RDKit::MolOps::RemoveHsParameters(),
                                    true);
  } else {
    check_sanitize_relevant_cycles(raw,
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
  }
}

void check_remove_all_hs_relevant_cycles(const RDKit::ROMol &mol,
                                         bool sanitize) {
  if (!sanitize || few_relevant_cycles(mol, true, true, true)) {
    return;
  }
  RDKit::RWMol copy(mol);
  try {
    RDKit::MolOps::removeAllHs(copy, false);
  } catch (...) {
    return;  // removeAllHs fails the same way before it sanitizes
  }
  check_sanitize_relevant_cycles(copy,
                                 RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
}

void check_frags_sanitize_relevant_cycles(const RDKit::ROMol &mol) {
  if (few_relevant_cycles(mol, true, true, true)) {
    return;
  }
  for (const auto &frag : RDKit::MolOps::getMolFrags(mol, false)) {
    check_sanitize_relevant_cycles(*frag,
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
  }
}

void check_products_relevant_cycles(const RDKit::ChemicalReaction &rxn,
                                    const RDKit::MOL_SPTR_VECT &reactants,
                                    unsigned max_products, bool sanitize,
                                    unsigned reactant_idx) {
  auto directed = [](const RDKit::ROMOL_SPTR &mol) {
    for (const auto *bond : mol->bonds()) {
      if (bond->getBondType() == RDKit::Bond::SINGLE &&
          bond->getBondDir() > RDKit::Bond::NONE) {
        return true;
      }
    }
    return false;
  };
  const bool directions =
      std::any_of(reactants.begin(), reactants.end(), directed);
  if (!directions && !sanitize) {
    return;
  }
  std::uint64_t bound = 0;
  for (const auto &templ : rxn.getProducts()) {
    bound += templ->getNumBonds();
  }
  for (std::size_t i = 0; i < reactants.size(); ++i) {
    const auto idx = reactant_idx == UINT_MAX ? i : reactant_idx;
    bound += independent_cycles(*reactants[i]);
    if (idx < rxn.getNumReactantTemplates()) {
      bound += rxn.getReactants()[idx]->getNumAtoms();
    }
  }
  if (bound <= MAX_INDEPENDENT_CYCLES) {
    return;
  }
  RDKit::MOL_SPTR_VECT copies;
  for (const auto &mol : reactants) {
    auto *copy = new RDKit::RWMol(*mol);
    copies.emplace_back(copy);
    for (auto *bond : copy->bonds()) {
      if (bond->getBondType() == RDKit::Bond::SINGLE) {
        bond->setBondDir(RDKit::Bond::NONE);
      }
    }
  }
  auto products = reactant_idx == UINT_MAX
                      ? rxn.runReactants(copies, max_products)
                      : rxn.runReactant(copies[0], reactant_idx);
  for (const auto &set : products) {
    for (const auto &product : set) {
      if (directions) {
        check_graph_relevant_cycles(*product);
      }
      if (sanitize) {
        check_sanitize_relevant_cycles(
            *product, RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
      }
    }
  }
}

std::size_t smiles_closure_bound(std::string_view smiles) {
  // A SMILES's bonds but its ring closures form a spanning forest, and each
  // closure is written as two labels of digits or a '%' outside brackets, so
  // half their characters bound its independent cycles.
  std::size_t label_chars = 0;
  bool in_bracket = false;
  for (char ch : smiles) {
    if (in_bracket) {
      in_bracket = ch != ']';
    } else if (ch == '[') {
      in_bracket = true;
    } else if (ch == '%' || (ch >= '0' && ch <= '9')) {
      ++label_chars;
    }
  }
  return label_chars / 2;
}

std::size_t smiles_cycle_bound(std::string_view text) {
  return static_cast<std::size_t>(std::count_if(
             text.begin(), text.end(),
             [](char c) { return c == '%' || (c >= '0' && c <= '9'); })) /
         2;
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
    if (smiles_closure_bound(std::string_view(text).substr(0, end)) <=
        MAX_INDEPENDENT_CYCLES) {
      return;
    }
  }
  auto unsanitized = params;
  unsanitized.sanitize = false;
  unsanitized.removeHs = false;
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(text, unsanitized));
  // Text that does not parse fails the same way when sanitized.
  if (!mol) {
    return;
  }
  if (!params.removeHs) {
    check_sanitize_relevant_cycles(*mol,
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
    return;
  }
  // MolFromSmiles removes hydrogens (sanitizing if asked), then perceives
  // bond and atom stereochemistry, which symmetrizes the rings of what is
  // left when it was not sanitized.
  RDKit::MolOps::RemoveHsParameters ps;
  ps.updateExplicitCount = true;
  if (params.sanitize) {
    check_remove_hs_relevant_cycles(*mol, ps, true);
    return;
  }
  if (few_relevant_cycles(*mol, true, true, true)) {
    return;
  }
  try {
    RDKit::MolOps::removeHs(*mol, ps, false);
  } catch (...) {
    return;  // removeHs fails the same way
  }
  check_perceive_relevant_cycles(*mol);
}

void check_smiles_relevant_cycles(
    const std::string &text,
    const RDKit::v2::SmilesParse::SmilesParserParams &params) {
  RDKit::SmilesParserParams v1;
  v1.sanitize = params.sanitize;
  v1.allowCXSMILES = params.allowCXSMILES;
  v1.strictCXSMILES = params.strictCXSMILES;
  v1.parseName = params.parseName;
  v1.removeHs = params.removeHs;
  v1.skipCleanup = params.skipCleanup;
  auto replacements = params.replacements;
  if (!replacements.empty()) {
    v1.replacements = &replacements;
  }
  check_smiles_relevant_cycles(text, v1);
}

// Mol blocks

namespace {

// As MolFromMolDataStream (MolFileParser.cpp) reads a block up to
// finishMolProcessing: the header, the counts line, and the CTAB, with
// RDKit's own CTAB readers, which perceive nothing. Throws as RDKit's reader
// does, and returns null for a block it skips (an unhandled CTAB feature) or
// an input that ends at once.
std::unique_ptr<RDKit::RWMol> read_ctab(std::istream &in,
                                        const FP::MolFileParserParams &params,
                                        bool &chirality_possible) {
  unsigned int line = 1;
  std::string text = RDKit::getLine(in);
  if (in.eof()) {
    return nullptr;
  }
  auto res = std::make_unique<RDKit::RWMol>();
  res->setProp(RDKit::common_properties::_Name, text);
  text = RDKit::getLine(in);
  if (text.length() >= 22) {
    std::string dim = text.substr(20, 2);
    if (dim == "3d" || dim == "3D") {
      res->setProp(RDKit::common_properties::_3DConf, 1);
    }
  }
  RDKit::getLine(in);
  line = 4;
  text = RDKit::getLine(in);
  if (text.size() < 6) {
    throw RDKit::FileParseException("Counts line too short");
  }
  unsigned int n_atoms = 0;
  unsigned int n_bonds = 0;
  unsigned int chiral_flag = 0;
  try {
    n_atoms = FPU::toUnsigned(text.substr(0, 3), true);
    n_bonds = FPU::toUnsigned(text.substr(3, 3), true);
  } catch (const boost::bad_lexical_cast &) {
    throw RDKit::FileParseException("Cannot convert counts");
  }
  try {
    if (text.size() >= 15) {
      chiral_flag = FPU::toUnsigned(text.substr(12, 3), true);
    }
  } catch (const boost::bad_lexical_cast &) {
    // RDKit ignores the rest of the counts line once it fails.
  }
  unsigned version = 2000;
  if (text.size() > 35) {
    if (text.size() < 39 || text[34] != 'V') {
      if (params.strictParsing) {
        throw RDKit::FileParseException("CTAB version string invalid");
      }
    } else if (text.substr(34, 5) == "V3000") {
      version = 3000;
    } else if (text.substr(34, 5) != "V2000") {
      if (params.strictParsing) {
        throw RDKit::FileParseException("Unsupported CTAB version");
      }
    } else if (params.parsingSCSRMol) {
      throw RDKit::FileParseException("SCSR Mol files is not V3000");
    }
  }
  res->setProp(RDKit::common_properties::_MolFileChiralFlag, chiral_flag);
  RDKit::Conformer *conf = nullptr;
  bool complete = false;
  try {
    if (version == 2000) {
      complete =
          FPU::ParseV2000CTAB(&in, line, res.get(), conf, chirality_possible,
                              n_atoms, n_bonds, params.strictParsing);
    } else {
      if ((n_atoms != 0 || n_bonds != 0) && params.strictParsing) {
        throw RDKit::FileParseException("V3000 counts line not 0");
      }
      complete = FPU::ParseV3000CTAB(
          &in, line, res.get(), conf, chirality_possible, n_atoms, n_bonds,
          params.strictParsing, !params.parsingSCSRMol, params.parsingSCSRMol);
    }
  } catch (const RDKit::MolFileUnhandledFeatureException &) {
    delete conf;
    return nullptr;
  } catch (...) {
    delete conf;
    throw;
  }
  if (!complete) {
    throw RDKit::FileParseException("M  END missing");
  }
  return res;
}

// Whether cleanUpOrganometallics could make a bond of the molecule dative:
// it only makes single bonds to metals dative.
bool has_single_bond_to_metal(const RDKit::ROMol &mol) {
  for (const auto *bond : mol.bonds()) {
    if (bond->getBondType() == RDKit::Bond::SINGLE &&
        (RDKit::QueryOps::isMetal(*bond->getBeginAtom()) ||
         RDKit::QueryOps::isMetal(*bond->getEndAtom()))) {
      return true;
    }
  }
  return false;
}

// The sanitizing branches of finishMolProcessing, and of the Marvin reader's
// copy of it (MarvinCMLReader::parseMolecule), for a molecule read and
// finished without sanitizing, whose graph is checked: sanitization
// symmetrizes the graph cleanUpOrganometallics leaves, after removeHs, by
// default, takes only hydrogens with one bond. That has the same relevant
// cycles unless a bond to a metal becomes dative.
void check_finished_sanitize(const RDKit::ROMol &mol, bool remove_hs) {
  if (!has_single_bond_to_metal(mol)) {
    return;
  }
  RDKit::RWMol copy(mol);
  if (!remove_hs) {
    check_sanitize_relevant_cycles(copy,
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
    return;
  }
  try {
    unsigned failed = 0;
    RDKit::MolOps::sanitizeMol(copy, failed,
                               RDKit::MolOps::SanitizeFlags::SANITIZE_CLEANUP);
    RDKit::MolOps::removeHs(copy, RDKit::MolOps::RemoveHsParameters(), false);
  } catch (...) {
    return;  // RDKit's reader fails the same way before it sanitizes
  }
  check_sanitize_relevant_cycles(copy,
                                 RDKit::MolOps::SanitizeFlags::SANITIZE_ALL);
}

}  // namespace

namespace {

// Checks a molecule RDKit's CTAB reader read, which finishMolProcessing
// finishes with `params` (rxn templates are finished the same way).
void check_ctab_relevant_cycles(const RDKit::RWMol &raw,
                                bool chirality_possible,
                                const FP::MolFileParserParams &params) {
  // The graph detectBondStereochemistry symmetrizes, unsanitized or after
  // SANITIZE_CLEANUP, which changes only single and double bonds.
  check_graph_relevant_cycles(raw);
  if (!params.sanitize || !has_single_bond_to_metal(raw)) {
    return;
  }
  // The molecule finished as RDKit's reader finishes it without sanitizing,
  // which symmetrizes the graph just checked.
  auto unsanitized = params;
  unsanitized.sanitize = false;
  RDKit::RWMol mol(raw);
  try {
    FPU::finishMolProcessing(&mol, chirality_possible, unsanitized);
  } catch (...) {
    // The steps it shares with the sanitizing reader come first; past them,
    // the molecule as read stands in for it.
    mol = raw;
    for (auto *atom : mol.atoms()) {
      atom->calcExplicitValence(false);
    }
  }
  check_finished_sanitize(mol, params.removeHs);
}

// A V3000 CTAB as the V3000 rxn reader reads each template
// (ParseV3000RxnBlock: ParseV3000CTAB without "M  END", then
// finishMolProcessing), checked.
void check_rxn_ctab(std::istream &in, unsigned int &line, bool strict,
                    const FP::MolFileParserParams &params) {
  RDKit::RWMol templ;
  RDKit::Conformer *conf = nullptr;
  unsigned int atoms = 0;
  unsigned int bonds = 0;
  bool possible = false;
  try {
    FPU::ParseV3000CTAB(&in, line, &templ, conf, possible, atoms, bonds, strict,
                        false);
  } catch (...) {
    delete conf;
    throw;
  }
  check_ctab_relevant_cycles(templ, possible, params);
}

void expect_v3000(std::istream &in, unsigned int &line, std::string_view what) {
  auto text = FPU::getV3000Line(&in, line);
  boost::to_upper(text);
  if (text.substr(0, what.size()) != what) {
    throw RDKit::FileParseException("missing V3000 line");
  }
}

}  // namespace

void check_mol_block_relevant_cycles(std::string_view block,
                                     const FP::MolFileParserParams &params) {
  std::istringstream in{std::string(block)};
  bool chirality_possible = false;
  std::unique_ptr<RDKit::RWMol> raw;
  try {
    raw = read_ctab(in, params, chirality_possible);
  } catch (...) {
    // RDKit's reader fails the same way before it perceives anything.
    return;
  }
  if (raw) {
    check_ctab_relevant_cycles(*raw, chirality_possible, params);
  }
}

void check_rxn_block_relevant_cycles(const std::string &text,
                                     const FP::MolFileParserParams &params) {
  // As ReactionFromRxnDataStream (MDLParser.cpp) reads the templates.
  std::istringstream in(text);
  try {
    auto header = RDKit::getLine(in);
    if (in.eof() || header.substr(0, 4) != "$RXN") {
      return;
    }
    if (header.size() >= 10 && header.substr(5, 5) == "V3000") {
      unsigned int line = 4;
      for (int i = 0; i < 3; ++i) {
        RDKit::getLine(in);
      }
      auto counts = FPU::getV3000Line(&in, line);
      boost::to_upper(counts);
      counts = boost::trim_copy(counts);
      std::vector<std::string> tokens;
      boost::split(tokens, counts, boost::is_any_of(" \t"),
                   boost::token_compress_on);
      if (tokens.size() < 3 || tokens[0] != "COUNTS") {
        return;
      }
      auto reactants = FPU::stripSpacesAndCast<unsigned int>(tokens[1]);
      auto products = FPU::stripSpacesAndCast<unsigned int>(tokens[2]);
      unsigned int agents = 0;
      if (tokens.size() > 3) {
        agents = FPU::stripSpacesAndCast<unsigned int>(tokens[3]);
      }
      expect_v3000(in, line, "BEGIN REACTANT");
      for (unsigned int i = 0; i < reactants; ++i) {
        check_rxn_ctab(in, line, params.strictParsing, params);
      }
      expect_v3000(in, line, "END REACTANT");
      expect_v3000(in, line, "BEGIN PRODUCT");
      for (unsigned int i = 0; i < products; ++i) {
        check_rxn_ctab(in, line, params.strictParsing, params);
      }
      expect_v3000(in, line, "END PRODUCT");
      if (agents) {
        expect_v3000(in, line, "BEGIN AGENT");
      }
      for (unsigned int i = 0; i < agents; ++i) {
        check_rxn_ctab(in, line, true, params);
      }
      return;
    }
    std::string counts;
    for (int i = 0; i < 4; ++i) {
      counts = RDKit::getLine(in);
    }
    if (in.eof() || counts.size() < 6) {
      return;
    }
    auto reactants = FPU::stripSpacesAndCast<unsigned int>(counts.substr(0, 3));
    auto products = FPU::stripSpacesAndCast<unsigned int>(counts.substr(3, 3));
    unsigned int agents = 0;
    if (counts.size() > 6 && !boost::trim_copy(counts.substr(6, 3)).empty()) {
      agents = FPU::stripSpacesAndCast<unsigned int>(counts.substr(6, 3));
    }
    // Agents are read neither sanitized nor with the strictness asked.
    FP::MolFileParserParams agent_params;
    agent_params.sanitize = false;
    const auto *agents_ps = &agent_params;
    for (auto [n, ps] : {std::pair{reactants + products, &params},
                         std::pair{agents, agents_ps}}) {
      for (unsigned int i = 0; i < n; ++i) {
        auto mol_line = RDKit::getLine(in);
        if (in.eof() || mol_line.substr(0, 4) != "$MOL") {
          return;
        }
        bool possible = false;
        auto templ = read_ctab(in, *ps, possible);
        if (!templ) {
          return;
        }
        check_ctab_relevant_cycles(*templ, possible, *ps);
      }
    }
  } catch (const ValueErrorException &) {
    throw;
  } catch (...) {
    // RDKit's reader fails the same way before it perceives anything more.
  }
}

namespace {

constexpr std::string_view kV30 = "M  V30 ";

std::string_view strip_cr(std::string_view text) {
  if (!text.empty() && text.back() == '\r') {
    text.remove_suffix(1);
  }
  return text;
}

bool starts_with(std::string_view text, std::string_view prefix) {
  return text.substr(0, prefix.size()) == prefix;
}

// The first four tokens of a V3000 line, as tokenizeV3000Line
// (MolFileParser.cpp) splits it; fewer if it has fewer.
std::vector<std::string_view> v3000_tokens(std::string_view line) {
  std::vector<std::string_view> tokens;
  bool in_quotes = false;
  unsigned paren_depth = 0;
  std::size_t start = 0;
  std::size_t pos = 0;
  while (pos < line.size() && tokens.size() < 4) {
    char c = line[pos];
    if (c == ' ' || c == '\t') {
      if (start == pos) {
        ++start;
        ++pos;
      } else if (!in_quotes && paren_depth == 0) {
        tokens.push_back(line.substr(start, pos - start));
        ++pos;
        start = pos;
      } else {
        ++pos;
      }
    } else if (c == ')' && paren_depth > 0) {
      --paren_depth;
      ++pos;
    } else if (c == '(' && !in_quotes) {
      ++paren_depth;
      ++pos;
    } else if (c == '"' && paren_depth == 0) {
      if (pos + 1 < line.size() && line[pos + 1] == '"') {
        pos += 2;
      } else if (in_quotes) {
        tokens.push_back(line.substr(start + 1, pos - start - 1));
        ++pos;
        start = pos;
        in_quotes = false;
      } else {
        ++pos;
        in_quotes = true;
      }
    } else {
      ++pos;
    }
  }
  if (start != pos && tokens.size() < 4) {
    tokens.push_back(line.substr(start));
  }
  return tokens;
}

// As ParseV3000BondBlock reads an atom number: std::from_chars, which leaves
// 0 when the token does not start with a digit.
std::uint64_t v3000_atom(std::string_view token) {
  unsigned value = 0;
  std::from_chars(token.data(), token.data() + token.size(), value);
  return value;
}

}  // namespace

MolRecordScanner::MolRecordScanner(const FP::MolFileParserParams &params,
                                   bool sd)
    : d_params(params), d_sd(sd) {
  d_records.emplace_back();
}

std::uint64_t MolRecordScanner::held_from() const {
  std::uint64_t held = UINT64_MAX;
  for (const auto &record : d_records) {
    held = std::min(held, record.held_from);
  }
  return held;
}

bool MolRecordScanner::line(std::string_view text, std::uint64_t at) {
  if (!d_error.empty()) {
    return false;
  }
  for (auto it = d_records.begin(); it != d_records.end();) {
    bool more = feed(*it, text, at);
    if (!d_error.empty()) {
      return false;
    }
    it = more ? it + 1 : d_records.erase(it);
  }
  // ForwardSDMolSupplier and the others compare the line as read, "\r" and
  // all.
  if (d_sd && starts_with(text, "$$$$")) {
    d_records.emplace_back();
  }
  return true;
}

// The record's bonds are read: one within the bound is done with.
void MolRecordScanner::end_bonds(Record &record, std::uint64_t next_line) {
  record.bonds_done = true;
  if (record.cycles.cycles() > MAX_INDEPENDENT_CYCLES) {
    record.held_from = next_line;
  } else {
    record.keep = false;
  }
}

bool MolRecordScanner::end_ctab(Record &record) {
  try {
    check_mol_block_relevant_cycles(record.text, d_params);
  } catch (const std::exception &e) {
    d_error = e.what();
    d_refuse_at = record.held_from;
  }
  return false;
}

bool MolRecordScanner::feed(Record &r, std::string_view raw, std::uint64_t at) {
  r.text.append(raw);
  r.text.push_back('\n');
  const std::string_view text = strip_cr(raw);
  const std::uint64_t next_line = at + raw.size() + 1;
  switch (r.state) {
    case State::kHeader:
      if (++r.lines == 3) {
        r.state = State::kCounts;
      }
      return true;
    case State::kCounts: {
      // As MolFromMolDataStream reads the counts line; a line RDKit fails on
      // ends the record.
      if (text.size() < 6) {
        return false;
      }
      try {
        r.atoms = FPU::toUnsigned(text.substr(0, 3), true);
        r.bonds = FPU::toUnsigned(text.substr(3, 3), true);
      } catch (const boost::bad_lexical_cast &) {
        return false;
      }
      bool v3000 = false;
      if (text.size() > 35) {
        if (text.size() < 39 || text[34] != 'V') {
          if (d_params.strictParsing) {
            return false;
          }
        } else if (text.substr(34, 5) == "V3000") {
          v3000 = true;
        } else if (text.substr(34, 5) != "V2000") {
          if (d_params.strictParsing) {
            return false;
          }
        } else if (d_params.parsingSCSRMol) {
          return false;
        }
      }
      r.lines = 0;
      if (v3000) {
        if ((r.atoms != 0 || r.bonds != 0) && d_params.strictParsing) {
          return false;
        }
        r.state = State::kV3000;
        return true;
      }
      // A block's independent cycles are at most its bonds.
      if (r.bonds <= MAX_INDEPENDENT_CYCLES) {
        return false;
      }
      r.state = r.atoms > 0 ? State::kAtoms : State::kBonds;
      return true;
    }
    case State::kAtoms:
      if (++r.lines == r.atoms) {
        r.lines = 0;
        r.state = State::kBonds;
      }
      return true;
    case State::kBonds: {
      // As ParseMolFileBondLine reads a bond's atoms.
      if (text.size() < 9) {
        return false;
      }
      try {
        auto begin = FPU::toUnsigned(text.substr(0, 3));
        auto end = FPU::toUnsigned(text.substr(3, 3));
        FPU::toUnsigned(text.substr(6, 3));
        r.cycles.add_bond(begin, end);
      } catch (const boost::bad_lexical_cast &) {
        return false;
      }
      if (++r.lines == r.bonds) {
        end_bonds(r, next_line);
        r.state = State::kProperties;
        return r.keep;
      }
      return true;
    }
    case State::kProperties: {
      // As ParseMolBlockProperties reads up to "M  END".
      if (r.skip > 0) {
        --r.skip;
        return true;
      }
      if (r.first_property) {
        r.first_property = false;
        if (text.empty()) {
          // A blank line before the properties, which strict parsing
          // refuses.
          return !d_params.strictParsing;
        }
      }
      if (starts_with(text, "M  END")) {
        return end_ctab(r);
      }
      if (starts_with(text, "$$$$")) {
        return false;  // RDKit's reader fails: no "M  END"
      }
      if (!text.empty() && (text[0] == 'A' || text[0] == 'G')) {
        r.skip = 1;
      } else if (starts_with(text, "S  SKP")) {
        try {
          int n = FPU::toInt(std::string(text.substr(6, 3)));
          if (n < 0) {
            return false;
          }
          r.skip = static_cast<std::uint64_t>(n);
        } catch (const boost::bad_lexical_cast &) {
          return false;
        }
      }
      return true;
    }
    case State::kV3000: {
      // As getV3000Line joins lines. RDKit fails on one that does not start
      // "M  V30 ", but for the line after the CTAB, which ends it.
      if (!starts_with(text, kV30)) {
        return r.bonds_done && end_ctab(r);
      }
      if (text.back() == '-') {
        r.joined.append(text.substr(7, text.size() - 8));
        return true;
      }
      r.joined.append(text.substr(7));
      std::string logical = std::move(r.joined);
      r.joined.clear();
      if (r.bonds_done) {
        return true;  // held until the CTAB ends
      }
      ++r.logical;
      if (r.logical == 2) {
        // As ParseV3000CTAB reads its COUNTS line; the first is BEGIN CTAB.
        boost::to_upper(logical);
        if (logical.size() < 8 || logical.substr(0, 7) != "COUNTS ") {
          return false;
        }
        std::string trimmed = boost::trim_copy(logical.substr(7));
        std::vector<std::string> counts;
        boost::split(counts, trimmed, boost::is_any_of(" \t"),
                     boost::token_compress_on);
        if (counts.size() < 2) {
          return false;
        }
        try {
          r.atoms = FPU::toUnsigned(counts[0]);
          r.bonds = FPU::toUnsigned(counts[1]);
        } catch (const boost::bad_lexical_cast &) {
          return false;
        }
        return r.bonds > MAX_INDEPENDENT_CYCLES;
      }
      // BEGIN ATOM, the atoms, END ATOM; BEGIN BOND, then the bonds.
      const std::uint64_t first_bond = 4 + (r.atoms > 0 ? r.atoms + 2 : 0);
      if (r.logical >= first_bond) {
        auto tokens = v3000_tokens(FPU::strip(logical));
        if (tokens.size() < 4) {
          return false;
        }
        r.cycles.add_bond(v3000_atom(tokens[2]), v3000_atom(tokens[3]));
        if (r.logical == first_bond + r.bonds - 1) {
          end_bonds(r, next_line);
          return r.keep;
        }
      }
      return true;
    }
  }
  return true;
}

void check_mol_text_relevant_cycles(std::string_view text,
                                    const FP::MolFileParserParams &params) {
  MolRecordScanner scanner(params, true);
  std::uint64_t at = 0;
  while (at < text.size()) {
    auto newline = text.find('\n', at);
    auto end = newline == std::string_view::npos ? text.size() : newline;
    if (!scanner.line(text.substr(at, end - at), at)) {
      throw ValueErrorException(scanner.error());
    }
    at = end + 1;
  }
}

// SCSR

namespace {

// `text` with the type of every V3000 bond (in a BEGIN BOND block, as
// ParseV3000BondBlock reads them) made dative (9), but hydrogen bonds (10),
// which the SCSR expansion treats apart. Sets `dative` if a bond was dative
// already.
std::string scsr_with_dative_bonds(const std::string &text, bool &dative) {
  std::istringstream in(text);
  std::string out;
  bool in_bonds = false;
  std::string joined;
  for (std::string raw; std::getline(in, raw);) {
    const bool newline = !in.eof();
    std::string_view line = strip_cr(raw);
    if (!starts_with(line, kV30)) {
      out += raw;
      if (newline) {
        out += '\n';
      }
      continue;
    }
    if (!line.empty() && line.back() == '-') {
      joined.append(line.substr(7, line.size() - 8));
      continue;
    }
    joined.append(line.substr(7));
    std::string logical = std::move(joined);
    joined.clear();
    if (starts_with(logical, "END BOND")) {
      in_bonds = false;
    } else if (in_bonds) {
      auto stripped = FPU::strip(logical);
      auto tokens = v3000_tokens(stripped);
      if (tokens.size() >= 2) {
        unsigned type = 0;
        std::from_chars(tokens[1].data(), tokens[1].data() + tokens[1].size(),
                        type);
        dative |= type == 9;
        const auto at =
            static_cast<std::size_t>(tokens[1].data() - logical.data());
        logical.replace(at, tokens[1].size(), type == 10 ? "10" : "9");
      }
    } else if (starts_with(logical, "BEGIN BOND")) {
      in_bonds = true;
    }
    out += kV30;
    out += logical;
    if (newline) {
      out += '\n';
    }
  }
  return out;
}

}  // namespace

void check_scsr_relevant_cycles(const std::string &text,
                                const FP::MolFileParserParams &params,
                                const FP::MolFromSCSRParams &scsr) {
  auto main = params;
  main.parsingSCSRMol = true;
  check_mol_block_relevant_cycles(text, main);
  // Each template, as SCSRMolFromSCSRDataStream reads it after the main
  // molecule and finishes it without sanitizing.
  try {
    std::istringstream in(text);
    bool chirality_possible = false;
    if (!read_ctab(in, main, chirality_possible)) {
      return;
    }
    unsigned int line = 0;
    if (FPU::getV3000Line(&in, line) != "BEGIN TEMPLATE") {
      return;
    }
    for (auto t = FPU::getV3000Line(&in, line); t.substr(0, 8) == "TEMPLATE";
         t = FPU::getV3000Line(&in, line)) {
      RDKit::RWMol templ;
      RDKit::Conformer *conf = nullptr;
      unsigned int atoms = 0;
      unsigned int bonds = 0;
      bool possible = false;
      bool complete = false;
      try {
        complete = FPU::ParseV3000CTAB(&in, line, &templ, conf, possible, atoms,
                                       bonds, params.strictParsing, false);
      } catch (...) {
        delete conf;
        throw;
      }
      if (!complete) {
        return;
      }
      check_graph_relevant_cycles(templ);
    }
  } catch (const ValueErrorException &) {
    throw;
  } catch (...) {
    return;  // RDKit fails the same way before it expands anything
  }
  // The expansion, on the copy with dative bonds. RDKit treats all but
  // hydrogen bonds alike when it expands; without base hydrogen bonds it
  // adds none.
  bool dative = false;
  std::istringstream in(scsr_with_dative_bonds(text, dative));
  auto probe_params = params;
  probe_params.sanitize = false;
  probe_params.removeHs = false;
  auto probe_scsr = scsr;
  probe_scsr.scsrBaseHbondOptions = FP::SCSRBaseHbondOptions::Ignore;
  std::unique_ptr<RDKit::RWMol> probe;
  try {
    unsigned int line = 0;
    probe = FP::MolFromSCSRDataStream(in, line, probe_params, probe_scsr);
  } catch (...) {
    return;  // RDKit fails the same way
  }
  if (!probe) {
    return;
  }
  // The expansion's graph: the copy's dative bonds.
  RDKit::RWMol graph(*probe);
  for (auto *bond : graph.bonds()) {
    bond->setBondType(bond->getBondType() == RDKit::Bond::DATIVE
                          ? RDKit::Bond::SINGLE
                          : RDKit::Bond::ZERO);
  }
  if (few_relevant_cycles(graph, false, false, false)) {
    return;
  }
  if (dative || (params.sanitize && has_single_bond_to_metal(graph))) {
    throw ValueErrorException(
        "an SCSR block whose expansion has more than 31 independent cycles, "
        "with dative bonds or bonds to metals, cannot be checked for "
        "relevant cycles RDKit could not enumerate");
  }
  check_graph_relevant_cycles(graph);
}

// CDXML

namespace {

// The order a CDXML bond's Order attribute (`value`, or null when it has
// none) gives it in the copy check_cdxml_relevant_cycles reads: "hydrogen"
// for an order RDKit's ChemDraw reader turns into a bond in the rings'
// graph (chemdraw/bond.cpp: the single orders, and its three query
// combinations and Any), "dative" for a dative or hydrogen bond, which is
// not, and the value as it is for one the reader skips the fragment for.
std::string cdxml_order(const std::string *value) {
  if (value == nullptr) {
    return "hydrogen";  // single
  }
  std::istringstream in(*value);
  std::set<std::string> tokens;
  for (std::string token; in >> token;) {
    tokens.insert(token);
  }
  static const std::set<std::string> singles = {
      "1", "2", "3", "4", "5", "6", "1.5", "2.5", "3.5", "4.5", "5.5", "ionic"};
  static const std::set<std::set<std::string>> queries = {
      {"1", "2"},
      {"1", "1.5"},
      {"2", "1.5"},
      {"1", "2", "3", "4", "5", "6", "0.5", "1.5", "2.5", "3.5", "4.5", "5.5",
       "dative", "ionic", "hydrogen", "threecenter"}};
  if (tokens.size() == 1 && singles.count(*tokens.begin())) {
    return "hydrogen";
  }
  if (queries.count(tokens)) {
    return "hydrogen";
  }
  if (tokens == std::set<std::string>{"dative"} ||
      tokens == std::set<std::string>{"hydrogen"}) {
    return "dative";
  }
  return *value;
}

// `text` with every bond element's Order set as cdxml_order says. Comments,
// CDATA sections, and processing instructions are copied as they are; a
// quoted attribute value may hold '>'.
std::string cdxml_with_hydrogen_bonds(const std::string &text) {
  std::string out;
  out.reserve(text.size() + text.size() / 8);
  std::size_t pos = 0;
  auto copy_through = [&](std::string_view end) {
    auto stop = text.find(end, pos);
    stop = stop == std::string::npos ? text.size() : stop + end.size();
    out.append(text, pos, stop - pos);
    pos = stop;
  };
  while (pos < text.size()) {
    auto lt = text.find('<', pos);
    if (lt == std::string::npos) {
      out.append(text, pos, std::string::npos);
      break;
    }
    out.append(text, pos, lt - pos);
    pos = lt;
    if (text.compare(pos, 4, "<!--") == 0) {
      copy_through("-->");
      continue;
    }
    if (text.compare(pos, 9, "<![CDATA[") == 0) {
      copy_through("]]>");
      continue;
    }
    if (text.compare(pos, 2, "<?") == 0 || text.compare(pos, 2, "<!") == 0) {
      copy_through(">");
      continue;
    }
    // The tag, to its '>' outside quotes.
    std::size_t end = pos + 1;
    char quote = 0;
    while (end < text.size() && (quote || text[end] != '>')) {
      if (quote) {
        quote = text[end] == quote ? 0 : quote;
      } else if (text[end] == '"' || text[end] == '\'') {
        quote = text[end];
      }
      ++end;
    }
    end = std::min(end + 1, text.size());
    std::string tag = text.substr(pos, end - pos);
    pos = end;
    const bool bond = tag.size() > 2 && tag[1] == 'b' &&
                      (tag[2] == '>' || tag[2] == '/' ||
                       std::isspace(static_cast<unsigned char>(tag[2])));
    if (!bond) {
      out += tag;
      continue;
    }
    // Its attributes: name, '=', quoted value.
    std::string rewritten = "<b";
    bool had_order = false;
    std::size_t i = 2;
    while (i < tag.size()) {
      auto c = tag[i];
      if (std::isspace(static_cast<unsigned char>(c))) {
        rewritten += c;
        ++i;
        continue;
      }
      if (c == '/' || c == '>') {
        break;
      }
      auto eq = tag.find('=', i);
      if (eq == std::string::npos) {
        break;
      }
      std::string name = boost::trim_copy(tag.substr(i, eq - i));
      auto open = tag.find_first_of("\"'", eq + 1);
      if (open == std::string::npos) {
        break;
      }
      auto close = tag.find(tag[open], open + 1);
      if (close == std::string::npos) {
        break;
      }
      if (name == "Order") {
        std::string value = tag.substr(open + 1, close - open - 1);
        rewritten += "Order=\"" + cdxml_order(&value) + "\"";
        had_order = true;
      } else {
        rewritten += tag.substr(i, close + 1 - i);
      }
      i = close + 1;
    }
    if (!had_order) {
      rewritten += " Order=\"" + cdxml_order(nullptr) + "\"";
    }
    rewritten += tag.substr(i);
    out += rewritten;
  }
  return out;
}

}  // namespace

void check_cdxml_relevant_cycles(const std::string &text, bool sanitize,
                                 bool remove_hs) {
  (void)remove_hs;  // removeHs takes only hydrogens with one bond
  if (text.compare(0, 8, "VjCD0100") == 0) {
    throw ValueErrorException(
        "binary CDX input cannot be checked for relevant cycles RDKit could "
        "not enumerate; it is read only as CDXML");
  }
  std::istringstream in(cdxml_with_hydrogen_bonds(text));
  RDKit::v2::CDXMLParser::CDXMLParserParams params(
      false, false, RDKit::v2::CDXMLParser::CDXMLFormat::CDXML);
  std::vector<std::unique_ptr<RDKit::RWMol>> frags;
  try {
    frags = RDKit::v2::CDXMLParser::MolsFromCDXMLDataStream(in, params);
  } catch (...) {
    return;  // RDKit fails the same way
  }
  for (const auto &frag : frags) {
    // A fragment's graph: the copy's hydrogen bonds.
    RDKit::RWMol graph(*frag);
    for (auto *bond : graph.bonds()) {
      bond->setBondType(bond->getBondType() == RDKit::Bond::HYDROGEN
                            ? RDKit::Bond::SINGLE
                            : RDKit::Bond::ZERO);
    }
    if (few_relevant_cycles(graph, false, false, false)) {
      continue;
    }
    if (sanitize && has_single_bond_to_metal(graph)) {
      throw ValueErrorException(
          "a CDXML fragment with more than 31 independent cycles and bonds to "
          "metals cannot be checked for relevant cycles RDKit could not "
          "enumerate");
    }
    check_graph_relevant_cycles(graph);
  }
}

// Mol2
// Mol2

std::uint64_t mol2_cycle_bound(const std::string &text) {
  // As MolFromMol2DataStream finds the first molecule's blocks: the lines
  // after its "@<TRIPOS>" markers, up to the next MOLECULE, skipping a last
  // line without a newline.
  using tokenizer = boost::tokenizer<boost::char_separator<char>>;
  boost::char_separator<char> sep(" \t\n");
  std::vector<std::string> lines;
  std::size_t mol_start = 0;
  std::size_t atom_start = 0;
  std::size_t bond_start = 0;
  {
    std::istringstream in(text);
    while (!in.eof() && !in.fail()) {
      std::string line = RDKit::getLine(in);
      if (in.eof()) {
        break;
      }
      lines.push_back(line);
      if (!line.empty() && line[0] == '@') {
        tokenizer tokens(line, sep);
        const std::string first = *tokens.begin();
        if (first == "@<TRIPOS>MOLECULE") {
          if (mol_start != 0) {
            break;
          }
          mol_start = lines.size();
        } else if (first == "@<TRIPOS>ATOM") {
          atom_start = lines.size();
        } else if (first == "@<TRIPOS>BOND") {
          bond_start = lines.size();
        }
      }
    }
  }
  if (mol_start == 0 || atom_start == 0) {
    return 0;  // RDKit fails
  }
  // RDKit seeks back to the blocks and reads on, through the last line.
  std::vector<std::string> all;
  bool unterminated = false;
  {
    std::istringstream in(text);
    while (true) {
      std::string line = RDKit::getLine(in);
      if (in.eof() && line.empty()) {
        break;
      }
      all.push_back(line);
      if (in.eof()) {
        unterminated = true;
        break;
      }
    }
  }
  if (mol_start + 1 >= all.size()) {
    return 0;
  }
  tokenizer counts(all[mol_start + 1], sep);
  if (counts.begin() == counts.end()) {
    return 0;
  }
  std::uint64_t n_bonds = 0;
  try {
    auto it = counts.begin();
    if (boost::lexical_cast<unsigned int>(*it) == 0) {
      return 0;
    }
    if (++it != counts.end()) {
      n_bonds = boost::lexical_cast<unsigned int>(*it);
    }
  } catch (const boost::bad_lexical_cast &) {
    return 0;
  }
  if (n_bonds <= MAX_INDEPENDENT_CYCLES) {
    return n_bonds;
  }
  if (bond_start == 0) {
    // RDKit would read bonds from the start of the input.
    return UINT64_MAX;
  }
  CycleCounter cycles;
  for (std::uint64_t i = 0; i < n_bonds; ++i) {
    // RDKit fails at the end of the input, or on a line without the
    // bond's atoms (ParseMol2FileBondLine).
    const auto idx = bond_start + i;
    if (idx >= all.size() || (idx + 1 == all.size() && unterminated)) {
      return 0;
    }
    tokenizer tokens(all[bond_start + i], sep);
    auto it = tokens.begin();
    std::uint64_t ends[2];
    try {
      for (auto &end : ends) {
        if (it == tokens.end() || ++it == tokens.end()) {
          return 0;
        }
        end = boost::lexical_cast<unsigned int>(*it);
      }
    } catch (const boost::bad_lexical_cast &) {
      return 0;
    }
    cycles.add_bond(ends[0], ends[1]);
  }
  return cycles.cycles();
}

// Marvin

namespace {

namespace pt = boost::property_tree;

// The bound on the independent cycles of the molecules RDKit reads from a
// Marvin document: those of every bondArray bond, over atoms named by their
// ids (an id two molecules share joins them, which only adds cycles), or
// UINT64_MAX if a MultipleSgroup has RDKit copy bonds.
std::uint64_t mrv_cycle_bound(const pt::ptree &tree) {
  CycleCounter cycles;
  bool copies = false;
  std::vector<const pt::ptree *> stack{&tree};
  while (!stack.empty() && !copies) {
    const auto *node = stack.back();
    stack.pop_back();
    for (const auto &[key, child] : *node) {
      if (key == "<xmlattr>") {
        copies |= child.get<std::string>("role", "") == "MultipleSgroup";
        continue;
      }
      if (key == "bondArray") {
        for (const auto &[bkey, bond] : child) {
          if (bkey != "bond") {
            continue;
          }
          std::vector<std::string> refs;
          boost::algorithm::split(
              refs, bond.get<std::string>("<xmlattr>.atomRefs2", ""),
              boost::algorithm::is_space(), boost::token_compress_on);
          refs.erase(std::remove(refs.begin(), refs.end(), ""), refs.end());
          if (refs.size() >= 2) {
            cycles.add_bond(std::hash<std::string>()(refs[0]),
                            std::hash<std::string>()(refs[1]));
          }
        }
      }
      stack.push_back(&child);
    }
  }
  return copies ? UINT64_MAX : cycles.cycles();
}

// The document without its atoms' coordinates, which RDKit reads into
// molecules without conformers: detectBondStereochemistry then perceives
// nothing.
std::string without_coordinates(pt::ptree tree) {
  std::vector<pt::ptree *> stack{&tree};
  while (!stack.empty()) {
    auto *node = stack.back();
    stack.pop_back();
    for (auto &[key, child] : *node) {
      if (key == "<xmlattr>") {
        for (const char *coord : {"x2", "y2", "x3", "y3", "z3"}) {
          child.erase(coord);
        }
      } else {
        stack.push_back(&child);
      }
    }
  }
  std::ostringstream out;
  pt::write_xml(out, tree);
  return out.str();
}

// Checks a molecule RDKit's Marvin reader read without coordinates: its
// graph, which detectBondStereochemistry symmetrizes (unsanitized, or after
// SANITIZE_CLEANUP), and the graphs sanitization symmetrizes.
void check_mrv_molecule(const RDKit::ROMol &raw, bool sanitize,
                        bool remove_hs) {
  check_graph_relevant_cycles(raw);
  if (sanitize) {
    check_finished_sanitize(raw, remove_hs);
  }
}

pt::ptree read_mrv(const std::string &text) {
  std::istringstream in(text);
  pt::ptree tree;
  pt::read_xml(in, tree);
  return tree;
}

}  // namespace

void check_mrv_relevant_cycles(const std::string &text, bool sanitize,
                               bool remove_hs) {
  pt::ptree tree;
  try {
    tree = read_mrv(text);
  } catch (const pt::ptree_error &) {
    return;  // RDKit fails to read it the same way
  }
  if (mrv_cycle_bound(tree) <= MAX_INDEPENDENT_CYCLES) {
    return;
  }
  std::istringstream in(without_coordinates(std::move(tree)));
  RDKit::v2::MarvinParser::MrvParserParams params;
  params.sanitize = false;
  params.removeHs = false;
  std::unique_ptr<RDKit::RWMol> raw;
  try {
    raw = RDKit::v2::MarvinParser::MolFromMrvDataStream(in, params);
  } catch (...) {
    return;  // RDKit fails the same way before it perceives anything
  }
  if (raw) {
    check_mrv_molecule(*raw, sanitize, remove_hs);
  }
}

void check_mrv_reaction_relevant_cycles(const std::string &text, bool sanitize,
                                        bool remove_hs) {
  pt::ptree tree;
  try {
    tree = read_mrv(text);
  } catch (const pt::ptree_error &) {
    return;
  }
  if (mrv_cycle_bound(tree) <= MAX_INDEPENDENT_CYCLES) {
    return;
  }
  std::istringstream in(without_coordinates(std::move(tree)));
  RDKit::v2::MarvinParser::MrvParserParams params;
  params.sanitize = false;
  params.removeHs = false;
  std::unique_ptr<RDKit::ChemicalReaction> raw;
  try {
    raw = RDKit::v2::MarvinParser::ReactionFromMrvDataStream(in, params);
  } catch (...) {
    return;
  }
  if (!raw) {
    return;
  }
  for (const auto &templates :
       {raw->getReactants(), raw->getAgents(), raw->getProducts()}) {
    for (const auto &templ : templates) {
      check_mrv_molecule(*templ, sanitize, remove_hs);
    }
  }
}

// SVG

void check_svg_relevant_cycles(const std::string &svg, bool sanitize,
                               bool remove_hs) {
  if (!sanitize) {
    return;  // removeHs alone perceives nothing
  }
  pt::ptree tree;
  try {
    tree = read_mrv(svg);
  } catch (const pt::ptree_error &) {
    return;
  }
  pt::ptree empty;
  for (const auto &[key, mol] : tree.get_child("svg.metadata", empty)) {
    if (key != "rdkit:mol") {
      continue;
    }
    // As ptreeToMol (SVGParser.cpp) reads each bond's atoms.
    CycleCounter cycles;
    for (const auto &[bkey, bond] : mol) {
      if (bkey != "rdkit:bond") {
        continue;
      }
      auto begin = bond.get_optional<int>("<xmlattr>.begin-atom-idx");
      auto end = bond.get_optional<int>("<xmlattr>.end-atom-idx");
      if (begin && end) {
        cycles.add_bond(static_cast<std::uint32_t>(*begin),
                        static_cast<std::uint32_t>(*end));
      }
    }
    if (cycles.cycles() <= MAX_INDEPENDENT_CYCLES) {
      continue;
    }
    // The molecule alone, read without sanitizing.
    pt::ptree one;
    one.add_child("svg.metadata.rdkit:mol", mol);
    std::ostringstream out;
    pt::write_xml(out, one);
    std::unique_ptr<RDKit::RWMol> raw;
    try {
      raw.reset(RDKit::RDKitSVGToMol(out.str(), false, false));
    } catch (...) {
      return;  // RDKit fails the same way before it sanitizes
    }
    if (raw && raw->getNumAtoms()) {
      check_parsed_relevant_cycles(*raw, true, remove_hs);
    }
  }
}

// SMILES files
// SMILES files

void check_smiles_line_relevant_cycles(
    std::string_view line, const FP::SmilesMolSupplierParams &params) {
  // As SmilesMolSupplier::processLine (and
  // MultithreadedSmilesMolSupplier::processMoleculeRecord) find the SMILES.
  std::string text(strip_cr(line));
  boost::char_separator<char> sep(params.delimiter.c_str(), "",
                                  boost::keep_empty_tokens);
  boost::tokenizer<boost::char_separator<char>> tokens(text, sep);
  int column = 0;
  for (const auto &token : tokens) {
    if (column++ == params.smilesColumn) {
      check_smiles_relevant_cycles(std::string(FPU::strip(token)),
                                   params.parseParameters);
      return;
    }
  }
}

bool SmilesLineScanner::line(std::string_view text) {
  if (!d_error.empty()) {
    return false;
  }
  const bool held = d_held_from != UINT64_MAX;
  const auto at = d_held_from;
  d_closures = 0;
  d_held_from = UINT64_MAX;
  if (d_tdt) {
    d_column = 0;
    d_tdt_line = true;
    if (!held) {
      return true;
    }
    try {
      // As TDTMolSupplier::parseMol.
      auto open = text.find('<');
      auto close = text.find_last_of('>');
      check_smiles_relevant_cycles(
          std::string(text.substr(open + 1, close - open - 1)), d_tdt_params);
    } catch (const std::exception &e) {
      d_error = e.what();
      d_refuse_at = at;
      return false;
    }
    return true;
  }
  // Blank lines and comments are skipped, then the title line.
  const auto stripped = strip_cr(text);
  if (stripped.empty() || stripped[0] == '#' ||
      stripped.find_first_not_of(" \t\r\n") == std::string_view::npos) {
    return true;
  }
  if (d_title) {
    d_title = false;
    return true;
  }
  if (!held) {
    return true;
  }
  try {
    check_smiles_line_relevant_cycles(text, d_params);
  } catch (const std::exception &e) {
    d_error = e.what();
    d_refuse_at = at;
    return false;
  }
  return true;
}

namespace {

// The atoms but hydrogens of each component of an InChI formula layer, in
// the order the connection layer numbers them, with the components repeated
// as their multipliers say.
std::vector<std::vector<int>> inchi_formula_atoms(std::string_view formula) {
  std::vector<std::vector<int>> components;
  const auto *table = RDKit::PeriodicTable::getTable();
  std::size_t pos = 0;
  while (pos <= formula.size()) {
    auto end = std::min(formula.find('.', pos), formula.size());
    auto part = formula.substr(pos, end - pos);
    std::size_t i = 0;
    std::uint64_t repeat = 0;
    while (i < part.size() &&
           std::isdigit(static_cast<unsigned char>(part[i]))) {
      repeat =
          std::min<std::uint64_t>(repeat * 10 + (part[i++] - '0'), 1u << 20);
    }
    std::vector<int> atoms;
    while (i < part.size()) {
      if (!std::isupper(static_cast<unsigned char>(part[i]))) {
        ++i;
        continue;
      }
      auto start = i++;
      while (i < part.size() &&
             std::islower(static_cast<unsigned char>(part[i]))) {
        ++i;
      }
      std::string symbol(part.substr(start, i - start));
      std::uint64_t count = 0;
      bool digits = false;
      while (i < part.size() &&
             std::isdigit(static_cast<unsigned char>(part[i]))) {
        count =
            std::min<std::uint64_t>(count * 10 + (part[i++] - '0'), 1u << 20);
        digits = true;
      }
      if (symbol == "H") {
        continue;
      }
      int num = -1;
      try {
        num = table->getAtomicNumber(symbol);
      } catch (...) {
        // An unknown symbol counts as a metal.
      }
      atoms.insert(atoms.end(), digits ? count : 1, num);
    }
    for (std::uint64_t k = 0; k < std::max<std::uint64_t>(repeat, 1); ++k) {
      components.push_back(atoms);
    }
    pos = end + 1;
  }
  return components;
}

// The bonds of each component of an InChI connection layer, as pairs of its
// atom numbers, with the components repeated as their multipliers say; and
// the bound on each one's independent cycles.
struct InchiComponent {
  std::vector<std::pair<std::uint64_t, std::uint64_t>> bonds;
  std::uint64_t cycles = 0;
};

std::vector<InchiComponent> inchi_connections(std::string_view layer) {
  std::vector<InchiComponent> components;
  std::size_t pos = 0;
  while (pos <= layer.size()) {
    auto end = std::min(layer.find(';', pos), layer.size());
    auto part = layer.substr(pos, end - pos);
    std::uint64_t repeat = 1;
    if (auto star = part.find('*'); star != std::string_view::npos) {
      repeat = 0;
      for (char c : part.substr(0, star)) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
          repeat = std::min<std::uint64_t>(repeat * 10 + (c - '0'), 1u << 20);
        }
      }
      part = part.substr(star + 1);
    }
    InchiComponent component;
    std::vector<std::uint64_t> branches;
    std::unordered_set<std::uint64_t> atoms;
    std::uint64_t numbers = 0;
    std::uint64_t current = 0;  // atom numbers start at 1
    for (std::size_t i = 0; i < part.size();) {
      char c = part[i];
      if (std::isdigit(static_cast<unsigned char>(c))) {
        std::uint64_t n = 0;
        while (i < part.size() &&
               std::isdigit(static_cast<unsigned char>(part[i]))) {
          n = std::min<std::uint64_t>(n * 10 + (part[i++] - '0'), 1ull << 40);
        }
        ++numbers;
        atoms.insert(n);
        if (current) {
          component.bonds.emplace_back(current, n);
        }
        current = n;
        continue;
      }
      if (c == '(') {
        branches.push_back(current);
      } else if (c == ',' && !branches.empty()) {
        current = branches.back();
      } else if (c == ')' && !branches.empty()) {
        current = branches.back();
        branches.pop_back();
      }
      ++i;
    }
    component.cycles = numbers - atoms.size();
    for (std::uint64_t k = 0; k < repeat; ++k) {
      components.push_back(component);
    }
    pos = end + 1;
  }
  return components;
}

}  // namespace

void check_inchi_relevant_cycles(std::string_view inchi, bool sanitize) {
  if (RDKit::MolOps::getUseLegacyRingFinding()) {
    return;
  }
  std::vector<std::string_view> layers;
  boost::split(layers, inchi, boost::is_any_of("/"));
  std::string_view formula;
  for (std::size_t i = 1; i < layers.size(); ++i) {
    const auto layer = layers[i];
    if (i == 1) {
      formula = layer;
      continue;
    }
    if (!layer.empty() && layer[0] == 'r') {
      formula = layer.substr(1);
      continue;
    }
    if (layer.empty() || layer[0] != 'c') {
      continue;
    }
    // A quick bound before anything is built: every number past the first
    // of a component, which adds at most one cycle, follows one of these.
    if (static_cast<std::size_t>(std::count_if(layer.begin(), layer.end(),
                                               [](char c) {
                                                 return c == '-' || c == '(' ||
                                                        c == ')' || c == ',';
                                               })) <= MAX_INDEPENDENT_CYCLES &&
        layer.find('*') == std::string_view::npos) {
      continue;
    }
    auto components = inchi_connections(layer.substr(1));
    std::uint64_t cycles = 0;
    for (const auto &component : components) {
      cycles += component.cycles;
    }
    if (cycles <= MAX_INDEPENDENT_CYCLES) {
      continue;
    }
    auto elements = inchi_formula_atoms(formula);
    RDKit::RWMol graph;
    bool metal_bond = false;
    for (std::size_t k = 0; k < components.size(); ++k) {
      std::unordered_map<std::uint64_t, unsigned> index;
      auto atom = [&](std::uint64_t n) {
        auto [it, added] = index.emplace(n, graph.getNumAtoms());
        if (added) {
          int num = -1;
          if (k < elements.size() && n - 1 < elements[k].size()) {
            num = elements[k][n - 1];
          }
          // The element only says whether sanitization could make the
          // atom's bonds dative.
          RDKit::Atom a(std::max(num, 0));
          if (num < 0 || RDKit::QueryOps::isMetal(a)) {
            metal_bond = true;
          }
          graph.addAtom(new RDKit::Atom(0), false, true);
        }
        return it->second;
      };
      for (const auto &[a, b] : components[k].bonds) {
        auto i = atom(a);
        auto j = atom(b);
        if (i != j && !graph.getBondBetweenAtoms(i, j)) {
          graph.addBond(i, j, RDKit::Bond::SINGLE);
        }
      }
    }
    if (sanitize && metal_bond) {
      throw ValueErrorException(
          "an InChI structure with more than " +
          std::to_string(MAX_INDEPENDENT_CYCLES) +
          " independent cycles and bonds to metals cannot be checked for its "
          "relevant cycles before RDKit perceives them");
    }
    check_graph_relevant_cycles(graph);
  }
}

}  // namespace rdkit_cxx
