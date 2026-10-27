// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/abbreviations.h"

#include <algorithm>
#include <mutex>
#include <set>
#include <string>
#include <utility>

#include <GraphMol/MolOps.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/Exceptions.h>
#include <boost/tokenizer.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/abbreviations.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

namespace A = RDKit::Abbreviations;
using tokenizer = boost::tokenizer<boost::char_separator<char>>;

AbbreviationDefinitionData to_data(const A::AbbreviationDefinition &def) {
  AbbreviationDefinitionData data{};
  data.label = std::make_unique<std::string>(def.label);
  data.display_label = std::make_unique<std::string>(def.displayLabel);
  data.display_label_w = std::make_unique<std::string>(def.displayLabelW);
  if (def.mol) {
    auto copy = std::make_unique<RWMol>(*def.mol);
    clear_matrix_caches(*copy);
    data.mol = std::move(copy);
  }
  return data;
}

rust::Vec<AbbreviationDefinitionData> to_data(
    const std::vector<A::AbbreviationDefinition> &defs) {
  rust::Vec<AbbreviationDefinitionData> res;
  for (const auto &def : defs) {
    res.push_back(to_data(def));
  }
  return res;
}

A::AbbreviationDefinition definition(rust::Str label, rust::Str display_label,
                                     rust::Str display_label_w) {
  A::AbbreviationDefinition def;
  def.label = std::string(label);
  def.displayLabel = std::string(display_label);
  def.displayLabelW = std::string(display_label_w);
  return def;
}

// parseAbbreviations reads each line's label and SMARTS without checking that
// the line has two fields: on a shorter line it dereferences, and increments,
// boost's token_iterator past its end, whose precondition is that it is valid
// (in a release build it yields the last token again, so a lone label is read
// as its own SMARTS). Such lines are left out.
std::string lines_with_smarts(const std::string &text) {
  boost::char_separator<char> line_sep("\n");
  boost::char_separator<char> field_sep(" \t");
  std::string res;
  for (const auto &line : tokenizer(text, line_sep)) {
    tokenizer fields(line, field_sep);
    auto field = fields.begin();
    if (field == fields.end() || ++field == fields.end()) {
      continue;
    }
    res += line;
    res += '\n';
  }
  return res;
}

// applyMatches (Abbreviations.cpp) marks the atoms it removes, and the bonds
// of each, in bitsets sized for the molecule as it was, but it adds bonds as
// it goes: from each neighbor outside the match of a removed atom with more
// bonds than its query atom, and from each extra attachment atom, to the
// match's connecting atom. A later match that removes an atom such a bond
// reaches, or attaches to one, has it index its bond bitset with the new
// bond's index, past the end (boost::dynamic_bitset checks only with
// assertions). It also indexes its atom bitset with each match's atoms
// unchecked. This walks the matches as applyMatches will and throws
// ValueErrorException before either; it stops where RDKit throws first.
void check_apply_matches(const RDKit::ROMol &mol,
                         const std::vector<A::AbbreviationMatch> &matches) {
  const unsigned num_atoms = mol.getNumAtoms();
  auto in_range = [num_atoms](int idx) {
    return idx >= 0 && static_cast<unsigned>(idx) < num_atoms;
  };
  std::vector<char> removed(num_atoms, 0);
  // Atoms with a bond applyMatches has added.
  std::vector<char> touched(num_atoms, 0);
  std::set<std::pair<unsigned, unsigned>> added;
  auto has_bond = [&mol, &added](unsigned a, unsigned b) {
    return mol.getBondBetweenAtoms(a, b) != nullptr ||
           added.count(std::minmax(a, b)) != 0;
  };
  auto add_bond = [&added, &touched](unsigned a, unsigned b) {
    added.insert(std::minmax(a, b));
    touched[a] = 1;
    touched[b] = 1;
  };
  auto refuse = [](unsigned idx) {
    throw ValueErrorException(
        "RDKit would index past the end of its bond bitset: atom " +
        std::to_string(idx) +
        " is removed by, or attached to, one abbreviation after another "
        "abbreviation bonded it to its connecting atom");
  };
  for (const auto &amatch : matches) {
    const auto &match = amatch.match;
    // RDKit 2024.03 always connects through the match's second atom.
    const unsigned which = 1;
    if (match.size() <= which || !in_range(match[which].second)) {
      return;  // match.at() or getAtomWithIdx throws
    }
    auto connect = static_cast<unsigned>(match[which].second);
    for (std::size_t i = which + 1; i < match.size(); ++i) {
      if (!in_range(match[i].second)) {
        throw ValueErrorException(
            "an abbreviation's atom index is out of range: " +
            std::to_string(match[i].second));
      }
      auto idx = static_cast<unsigned>(match[i].second);
      if (removed[idx]) {
        return;  // the "overlapping matches" invariant throws
      }
      removed[idx] = 1;
      if (touched[idx]) {
        refuse(idx);
      }
      const auto &query = amatch.abbrev.mol;
      if (!query) {
        continue;
      }
      int query_idx = match[i].first;
      if (query_idx < 0 ||
          static_cast<unsigned>(query_idx) >= query->getNumAtoms()) {
        return;  // getAtomWithIdx throws
      }
      // The atom has no added bond, so its degree is the molecule's.
      const auto *atom = mol.getAtomWithIdx(idx);
      if (atom->getDegree() <= query->getAtomWithIdx(query_idx)->getDegree()) {
        continue;
      }
      for (const auto bond : mol.atomBonds(atom)) {
        auto nbr = bond->getOtherAtomIdx(idx);
        if (std::any_of(match.begin(), match.end(), [nbr](const auto &pr) {
              return pr.second == static_cast<int>(nbr);
            })) {
          continue;
        }
        if (has_bond(nbr, connect)) {
          return;  // addBond's precondition throws
        }
        add_bond(nbr, connect);
      }
    }
    for (auto other : amatch.abbrev.extraAttachAtoms) {
      if (other >= num_atoms) {
        return;  // getAtomWithIdx throws
      }
      if (touched[other]) {
        refuse(other);
      }
      unsigned marked = 0;
      for (const auto bond : mol.atomBonds(mol.getAtomWithIdx(other))) {
        if (removed[bond->getBeginAtomIdx()] ||
            removed[bond->getEndAtomIdx()]) {
          ++marked;
        }
      }
      if (marked != 1 || other == connect || has_bond(other, connect)) {
        return;  // an invariant or addBond's precondition throws
      }
      add_bond(other, connect);
    }
  }
}

// The matches condenseAbbreviationSubstanceGroups makes from the SUP groups,
// as it makes them; RDKit's accessors throw where it would.
std::vector<A::AbbreviationMatch> sup_matches(const RDKit::ROMol &mol) {
  std::vector<A::AbbreviationMatch> res;
  for (const auto &sg : RDKit::getSubstanceGroups(mol)) {
    if (sg.getProp<std::string>("TYPE") != "SUP") {
      continue;
    }
    A::AbbreviationMatch amatch;
    auto ats = sg.getAtoms();
    const auto &bnds = sg.getBonds();
    bool first_attach_found = false;
    for (auto bidx : bnds) {
      const auto *bnd = mol.getBondWithIdx(bidx);
      unsigned match_atom;
      unsigned other_atom;
      if (std::find(ats.begin(), ats.end(), bnd->getBeginAtomIdx()) !=
          ats.end()) {
        other_atom = bnd->getEndAtomIdx();
        match_atom = bnd->getBeginAtomIdx();
      } else if (std::find(ats.begin(), ats.end(), bnd->getEndAtomIdx()) !=
                 ats.end()) {
        other_atom = bnd->getBeginAtomIdx();
        match_atom = bnd->getEndAtomIdx();
      } else {
        continue;
      }
      if (!first_attach_found) {
        if (*ats.begin() != match_atom) {
          ats.erase(std::find(ats.begin(), ats.end(), match_atom));
          ats.insert(ats.begin(), match_atom);
        }
        ats.insert(ats.begin(), other_atom);
        first_attach_found = true;
      } else {
        amatch.abbrev.extraAttachAtoms.push_back(other_atom);
      }
    }
    for (unsigned i = 0; i < ats.size(); ++i) {
      amatch.match.push_back({i, ats[i]});
    }
    res.push_back(std::move(amatch));
  }
  return res;
}

}  // namespace

std::unique_ptr<AbbreviationList> abbreviation_list_new() {
  return std::make_unique<AbbreviationList>();
}

void abbreviation_list_push(AbbreviationList &list, rust::Str label,
                            rust::Str display_label,
                            rust::Str display_label_w) {
  list.defs.push_back(definition(label, display_label, display_label_w));
}

void abbreviation_list_push_with_mol(AbbreviationList &list, rust::Str label,
                                     rust::Str display_label,
                                     rust::Str display_label_w,
                                     const RWMol &mol) {
  auto def = definition(label, display_label, display_label_w);
  // RDKit only reads the definitions' molecules (and records recursive
  // queries' matches in them, under their own locks).
  def.mol.reset(const_cast<RWMol *>(&mol), [](RDKit::ROMol *) {});
  list.defs.push_back(std::move(def));
}

AbbreviationDefinitionData abbreviation_definition_default() {
  return to_data(A::AbbreviationDefinition());
}

// The definitions are made once, under a function-local static, and their
// molecules are only read here: nothing in the crate matches against them.
rust::Vec<AbbreviationDefinitionData> get_default_abbreviations() {
  return to_data(A::Utils::getDefaultAbbreviations());
}

rust::Vec<AbbreviationDefinitionData> get_default_linkers() {
  return to_data(A::Utils::getDefaultLinkers());
}

rust::Vec<AbbreviationDefinitionData> parse_abbreviations(
    rust::Str text, bool remove_extra_dummies,
    bool allow_connection_to_dummies) {
  return to_data(A::Utils::parseAbbreviations(
      lines_with_smarts(std::string(text)), remove_extra_dummies,
      allow_connection_to_dummies));
}

// As condenseMolAbbreviations, which finds the matches and applies them, with
// the matches checked in between. Finding them matches substructures, and
// every match constructs a ControlCHandler (control_c_lock.h); applying them
// removes atoms (check_removals).
std::unique_ptr<RWMol> condense_mol_abbreviations(
    const RWMol &mol, const AbbreviationList &abbrevs, double max_coverage,
    bool sanitize) {
  auto res = std::make_unique<RWMol>(mol);
  check_removals(*res);
  std::vector<A::AbbreviationMatch> matches;
  {
    std::lock_guard<ControlCMutex> control_c(control_c_mutex());
    matches =
        A::findApplicableAbbreviationMatches(*res, abbrevs.defs, max_coverage);
  }
  check_apply_matches(*res, matches);
  A::applyMatches(*res, matches);
  if (sanitize && !res->getRingInfo()->isSymmSssr()) {
    RDKit::MolOps::symmetrizeSSSR(*res);
  }
  return res;
}

// labelMatches only adds substance groups, reading each match's first two
// atoms, which findApplicableAbbreviationMatches checks are there.
std::unique_ptr<RWMol> label_mol_abbreviations(const RWMol &mol,
                                               const AbbreviationList &abbrevs,
                                               double max_coverage) {
  auto res = std::make_unique<RWMol>(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  A::labelMolAbbreviations(*res, abbrevs.defs, max_coverage);
  return res;
}

std::unique_ptr<RWMol> condense_abbreviation_substance_groups(
    const RWMol &mol) {
  auto res = std::make_unique<RWMol>(mol);
  check_removals(*res);
  check_apply_matches(*res, sup_matches(*res));
  A::condenseAbbreviationSubstanceGroups(*res);
  return res;
}

}  // namespace rdkit_cxx
