// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/rgroup.h"

#include <chrono>
#include <climits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <RDGeneral/Exceptions.h>
#include <RDGeneral/types.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/rgroup.rs.h"
#include "rdkit-cxx/ffi/substruct.rs.h"
#include "rdkit-cxx/generalized_substruct.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/substruct.h"
#include "rdkit-cxx/tautomer_checks.h"

// RDKit exports relabelMappedDummies but does not install RGroupUtils.h,
// which declares it; this is RDKit 2026.09's declaration (CMake checks the
// release series).
namespace RDKit {
RDKIT_RGROUPDECOMPOSITION_EXPORT void relabelMappedDummies(
    ROMol &mol, unsigned int inputLabels, unsigned int outputLabels);
}  // namespace RDKit

namespace rdkit_cxx {

std::mutex &rgroup_mutex() {
  static std::mutex mutex;
  return mutex;
}

namespace {

// RDKit's own internal properties (RGroupDecomp.cpp), which a core's atoms
// may carry as input.
const std::string kRLabel = "tempRlabel";

// The locks for a call that can score, run the GA, or match (rgroup.h).
class Locks {
 public:
  Locks() : d_rgroup(rgroup_mutex()), d_control_c(control_c_mutex()) {}

 private:
  std::lock_guard<std::mutex> d_rgroup;
  std::lock_guard<std::mutex> d_control_c;
};

RDKit::RGroupDecompositionParameters params_from(
    const RGroupDecompositionParametersData &d) {
  RDKit::RGroupDecompositionParameters p;
  p.labels = d.labels;
  p.matchingStrategy = d.matching_strategy;
  p.scoreMethod = d.score_method;
  p.rgroupLabelling = d.rgroup_labelling;
  p.alignment = d.alignment;
  p.chunkSize = d.chunk_size;
  p.onlyMatchAtRGroups = d.only_match_at_r_groups;
  p.removeAllHydrogenRGroups = d.remove_all_hydrogen_r_groups;
  p.removeAllHydrogenRGroupsAndLabels =
      d.remove_all_hydrogen_r_groups_and_labels;
  p.removeHydrogensPostMatch = d.remove_hydrogens_post_match;
  p.allowNonTerminalRGroups = d.allow_non_terminal_r_groups;
  p.allowMultipleRGroupsOnUnlabelled = d.allow_multiple_r_groups_on_unlabelled;
  p.allowMultipleCoresInSameMol = d.allow_multiple_cores_in_same_mol;
  p.doTautomers = d.do_tautomers;
  p.doEnumeration = d.do_enumeration;
  p.includeTargetMolInResults = d.include_target_mol_in_results;
  p.timeout = d.timeout;
  p.gaPopulationSize = d.ga_population_size;
  p.gaMaximumOperations = d.ga_maximum_operations;
  p.gaNumberOperationsWithoutImprovement =
      d.ga_number_operations_without_improvement;
  p.gaRandomSeed = d.ga_random_seed;
  p.gaNumberRuns = d.ga_number_runs;
  // RDKit's parallel runs share the GA's random number generator and the
  // scoring static (rgroup.h); the runs go one after another here.
  p.gaParallelRuns = false;
  p.substructmatchParams = substruct_params_from_data(d.substruct_match_params);
  return p;
}

// RGroupDecomposition::add processes every chunkSize molecules with
// GreedyChunks (and not Greedy), taking the number of molecules modulo the
// chunk size, which is undefined (SIGFPE in Python) for 0.
void check_params(const RDKit::RGroupDecompositionParameters &p) {
  if (p.chunkSize == 0 && (p.matchingStrategy & RDKit::GreedyChunks) &&
      !(p.matchingStrategy & RDKit::Greedy)) {
    throw ValueErrorException("GreedyChunks needs a chunkSize above 0");
  }
}

// prepareCore turns a core's labels (an atom's tempRlabel property, and a
// dummy atom's MDL R label, isotope, or atom map number; isotopes have 16
// bits) into ints and
// keeps one past the largest in int (setLabel's maxLabel = label + 1, and
// prepareCore's rlabel + 1 for an aligned core), which overflows for
// INT_MAX; relabelling adds one per atom past that. Labels are refused at
// 2^30 and above: RDKit's R-group labels are below 1000 elsewhere (the MDL
// labels setAtomRLabel writes, atom map numbers).
constexpr std::int64_t kMaxLabel = std::int64_t{1} << 30;

void check_label(std::int64_t label) {
  if (label >= kMaxLabel) {
    throw ValueErrorException("R-group labels must be below 2^30");
  }
}

// MCS alignment runs findMCS on the cores with RDKit's default bond
// comparison, MCSBondCompareOrder, which indexes a table by both bonds'
// types up to ZERO; a molecule pickle can carry any type up to 127, which
// would be read past the table. doEnumeration enumerates the cores with
// MolEnumerator, whose link nodes check_link_nodes checks (Python segfaults
// on a link node over a double bond). doTautomers matches with RDKit's
// TautomerQuery::substructOf, which reads a recursive query's matches past
// their end (tautomer_checks.h; Python fails a precondition on a bond
// between atoms it read there); RDKit makes its tautomer query of the core
// less its terminal R groups, whose modified atoms are not known here, so
// any recursive query is refused.
void check_cores(const std::vector<boost::shared_ptr<RDKit::ROMol>> &cores,
                 const RDKit::RGroupDecompositionParameters &p) {
  bool mcs = !p.onlyMatchAtRGroups && (p.alignment & RDKit::MCS);
  for (const auto &core : cores) {
    for (const auto *atom : core->atoms()) {
      // A label of another type, which RDKit's own read refuses too, is
      // left to RDKit.
      int rlabel = 0;
      bool has_rlabel = false;
      try {
        has_rlabel = atom->getPropIfPresent(kRLabel, rlabel);
      } catch (const std::exception &) {
      }
      if (has_rlabel) {
        check_label(rlabel);
      }
      if (atom->getAtomicNum() == 0) {
        check_label(atom->getAtomMapNum());
        unsigned int mdl = 0;
        bool has_mdl = false;
        try {
          has_mdl = atom->getPropIfPresent(
              RDKit::common_properties::_MolFileRLabel, mdl);
        } catch (const std::exception &) {
        }
        if (has_mdl) {
          check_label(mdl);
        }
      }
    }
    if (p.doEnumeration) {
      check_link_nodes(*core);
    }
    if (p.doTautomers) {
      check_tautomer_template(*core);
    }
    if (mcs) {
      for (const auto *bond : core->bonds()) {
        if (bond->getBondType() > RDKit::Bond::ZERO) {
          throw ValueErrorException(
              "core bond " + std::to_string(bond->getIdx()) +
              " has a type the MCS bond comparison does not know");
        }
      }
    }
  }
}

std::unique_ptr<RGroupDecomposition> make_decomposition(
    BorrowedMols &cores, const RGroupDecompositionParametersData &data) {
  auto p = params_from(data);
  check_params(p);
  check_cores(cores.mols(), p);
  return std::make_unique<RGroupDecomposition>(cores.mols(), p);
}

std::unique_ptr<RWMol> copy_out(const RDKit::ROMOL_SPTR &mol) {
  if (!mol) {
    throw std::runtime_error("RDKit returned a null molecule");
  }
  auto copy = std::make_unique<RWMol>(*mol);
  clear_matrix_caches(*copy);
  return copy;
}

RGroupTable rows_table(const RDKit::RGroupRows &rows) {
  RGroupTable table{std::make_unique<MolList>(), {}};
  for (const auto &row : rows) {
    for (const auto &[label, mol] : row) {
      table.mols->push_named(label, copy_out(mol));
    }
    table.sizes.push_back(static_cast<std::uint32_t>(row.size()));
  }
  return table;
}

RGroupTable columns_table(const RDKit::RGroupColumns &columns) {
  RGroupTable table{std::make_unique<MolList>(), {}};
  for (const auto &[label, column] : columns) {
    for (const auto &mol : column) {
      table.mols->push_named(label, copy_out(mol));
    }
    table.sizes.push_back(static_cast<std::uint32_t>(column.size()));
  }
  return table;
}

}  // namespace

RGroupDecompositionParametersData r_group_decomposition_parameters_default() {
  RDKit::RGroupDecompositionParameters p;
  return RGroupDecompositionParametersData{
      p.labels,
      p.matchingStrategy,
      p.scoreMethod,
      p.rgroupLabelling,
      p.alignment,
      p.chunkSize,
      p.onlyMatchAtRGroups,
      p.removeAllHydrogenRGroups,
      p.removeAllHydrogenRGroupsAndLabels,
      p.removeHydrogensPostMatch,
      p.allowNonTerminalRGroups,
      p.allowMultipleRGroupsOnUnlabelled,
      p.allowMultipleCoresInSameMol,
      p.doTautomers,
      p.doEnumeration,
      p.includeTargetMolInResults,
      p.timeout,
      p.gaPopulationSize,
      p.gaMaximumOperations,
      p.gaNumberOperationsWithoutImprovement,
      p.gaRandomSeed,
      p.gaNumberRuns,
      p.gaParallelRuns,
      substruct_params_to_data(p.substructmatchParams)};
}

std::unique_ptr<RGroupDecomposition> r_group_decomposition_new(
    BorrowedMols &cores, const RGroupDecompositionParametersData &params) {
  Locks locks;
  return make_decomposition(cores, params);
}

std::int32_t r_group_decomposition_add(RGroupDecomposition &decomp,
                                       const RWMol &mol) {
  Locks locks;
  return decomp.add(mol);
}

MatchingCoreResult r_group_decomposition_get_matching_core_idx(
    RGroupDecomposition &decomp, const RWMol &mol) {
  std::vector<RDKit::MatchVectType> matches;
  int idx;
  {
    Locks locks;
    idx = decomp.getMatchingCoreIdx(mol, &matches);
  }
  MatchingCoreResult res{idx, {}};
  for (const auto &match : matches) {
    IndexList pairs;
    for (const auto &[core_idx, mol_idx] : match) {
      pairs.indices.push_back(static_cast<std::uint32_t>(core_idx));
      pairs.indices.push_back(static_cast<std::uint32_t>(mol_idx));
    }
    res.matches.push_back(std::move(pairs));
  }
  return res;
}

bool r_group_decomposition_process(RGroupDecomposition &decomp) {
  Locks locks;
  return decomp.process();
}

ProcessAndScoreResult r_group_decomposition_process_and_score(
    RGroupDecomposition &decomp) {
  Locks locks;
  auto result = decomp.processAndScore();
  return ProcessAndScoreResult{result.success, result.score};
}

std::unique_ptr<std::vector<std::string>>
r_group_decomposition_get_r_group_labels(RGroupDecomposition &decomp) {
  return std::make_unique<std::vector<std::string>>(decomp.getRGroupLabels());
}

RGroupTable r_group_decomposition_get_r_groups_as_rows(
    RGroupDecomposition &decomp) {
  return rows_table(decomp.getRGroupsAsRows());
}

RGroupTable r_group_decomposition_get_r_groups_as_columns(
    RGroupDecomposition &decomp) {
  return columns_table(decomp.getRGroupsAsColumns());
}

// As RGroupDecomp in rdRGroupComposition.cpp: the decomposition Python
// makes, its adds with the timeout checked after each, and its rows or
// columns.
RGroupDecomposeResult r_group_decompose(
    BorrowedMols &cores, BorrowedMols &mols,
    const RGroupDecompositionParametersData &params, bool as_rows) {
  auto t0 = std::chrono::steady_clock::now();
  RGroupDecomposeResult res{{std::make_unique<MolList>(), {}}, {}};
  Locks locks;
  auto decomp = make_decomposition(cores, params);
  std::uint32_t idx = 0;
  for (const auto &mol : mols.mols()) {
    if (decomp->add(*mol) == -1) {
      res.unmatched.push_back(idx);
    }
    ++idx;
    RDKit::checkForTimeout(t0, params.timeout);
  }
  decomp->process();
  res.table = as_rows ? rows_table(decomp->getRGroupsAsRows())
                      : columns_table(decomp->getRGroupsAsColumns());
  return res;
}

// relabelMappedDummies takes abs() of an atom map number, a property any
// int can be set in, which is undefined for INT_MIN.
void relabel_mapped_dummies(RWMol &mol, std::uint32_t input_labels,
                            std::uint32_t output_labels) {
  if (input_labels & RDKit::AtomMap) {
    for (const auto *atom : mol.atoms()) {
      if (!atom->getAtomicNum() && atom->getAtomMapNum() == INT_MIN) {
        throw ValueErrorException("atom " + std::to_string(atom->getIdx()) +
                                  " has atom map number INT_MIN");
      }
    }
  }
  RDKit::relabelMappedDummies(mol, input_labels, output_labels);
}

}  // namespace rdkit_cxx
