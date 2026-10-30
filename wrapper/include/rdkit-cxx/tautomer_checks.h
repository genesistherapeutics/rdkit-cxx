/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// RDKit's TautomerQuery::substructOf checks each match of the template
// against the tautomers in the template match's extraFinalCheck, which reads
// the match at the query's modified atoms and bonds. RDKit's
// RecursiveMatcher copies the parameters, check included, into the match of
// each recursive SMARTS, so the check reads that SMARTS's own, shorter,
// matches at the template's indices, past their end (Python's
// SubstructLibrary.GetMatches then reads atom 2532319792 of six).
// tautomer_query.cpp has its own matcher, which passes those matches;
// wrappers whose RDKit call matches a tautomer query with RDKit's own
// (SubstructLibrary's searches, RGroupDecomposition's doTautomers) refuse
// it instead.

#pragma once

#include <vector>

#include <GraphMol/ROMol.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

inline bool has_recursive_query(const RDKit::ROMol &mol) {
  for (const auto *atom : mol.atoms()) {
    if (!atom->hasQuery()) {
      continue;
    }
    std::vector<const RDKit::Atom::QUERYATOM_QUERY *> stack{atom->getQuery()};
    while (!stack.empty()) {
      const auto *query = stack.back();
      stack.pop_back();
      if (query->getDescription() == "RecursiveStructure") {
        return true;
      }
      for (auto child = query->beginChildren(); child != query->endChildren();
           ++child) {
        stack.push_back(child->get());
      }
    }
  }
  return false;
}

// For a template RDKit makes a tautomer query of: refuses any recursive
// query, since the modified atoms are not known yet.
inline void check_tautomer_template(const RDKit::ROMol &templ) {
  if (has_recursive_query(templ)) {
    throw ValueErrorException(
        "RDKit's tautomer matcher reads a recursive query's matches as the "
        "template's, past their end; this template has a recursive query");
  }
}

// A query whose check reads nothing (no modified atoms or bonds) is safe.
inline void check_rdkit_tautomer_match(const RDKit::TautomerQuery &query) {
  if (!query.getModifiedAtoms().empty() || !query.getModifiedBonds().empty()) {
    check_tautomer_template(query.getTemplateMolecule());
  }
}

}  // namespace rdkit_cxx
