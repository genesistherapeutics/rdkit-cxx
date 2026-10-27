// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/standardize_tautomer.h"

#include <climits>
#include <cstdlib>
#include <iterator>
#include <mutex>
#include <utility>

#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/standardize.rs.h"
#include "rdkit-cxx/ffi/standardize_tautomer.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"
#include "rdkit-cxx/standardize.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace MS = RDKit::MolStandardize;
namespace TSF = RDKit::MolStandardize::TautomerScoringFunctions;

namespace {

// A tautomer enumerator callback whose operator() is a Rust closure. RDKit
// calls it from the thread that enumerates; enumerators copied from one
// another share it, which is why it is Send + Sync. It gets its own copy of
// the molecule, and a copy of the result so far that shares its tautomers
// with RDKit's (RDKit waits while the closure runs, and the copy is dropped
// before it resumes). A panic in it becomes an exception in Rust's half of
// the bridge, which leaves the enumeration. The closure runs without the
// enumeration's control_c_mutex, so it can call wrappers that take it; RDKit
// calls it between matches.
class RustTautomerCallback : public MS::TautomerEnumeratorCallback {
 public:
  explicit RustTautomerCallback(rust::Box<TautomerCallback> callback)
      : d_callback(std::move(callback)) {}

  bool operator()(const RDKit::ROMol &mol,
                  const TautomerEnumeratorResult &res) override {
    auto copy = std::make_unique<RWMol>(mol);
    clear_matrix_caches(*copy);
    clear_descriptor_caches(*copy);
    ReleaseControlC release;
    return call_tautomer_callback(
        *d_callback, std::move(copy),
        std::make_unique<TautomerEnumeratorResult>(res));
  }

 private:
  rust::Box<TautomerCallback> d_callback;
};

// A scoring function that calls a Rust closure with a copy of each
// tautomer. pickCanonical starts its search at INT_MIN and dereferences the
// best tautomer, null when none scores above that, so INT_MIN counts as
// INT_MIN + 1. The closure runs without control_c_mutex, as the callback's
// does.
boost::function<int(const RDKit::ROMol &)> scorer(const ScoreFunction &score) {
  return [&score](const RDKit::ROMol &mol) -> int {
    auto copy = std::make_unique<RWMol>(mol);
    clear_matrix_caches(*copy);
    clear_descriptor_caches(*copy);
    ReleaseControlC release;
    int value = call_score_function(score, std::move(copy));
    return value == INT_MIN ? INT_MIN + 1 : value;
  };
}

// The default score is checked for each tautomer it scores. scoreRings
// gives a molecule the symmetrized SSSR when it lacks it; that is done here
// first, as RDKit would, so its rings can be counted.
void check_default_score(const RDKit::ROMol &mol) {
  check_tautomer_score(mol, 0);
  if (!mol.getRingInfo()->isSymmSssr()) {
    check_relevant_cycles(mol);
    RDKit::MolOps::symmetrizeSSSR(const_cast<RDKit::ROMol &>(mol));
  }
  check_tautomer_score(mol, mol.getRingInfo()->numRings());
}

std::unique_ptr<RWMol> copy_mol(const RDKit::ROMol &mol) {
  auto result = std::make_unique<RWMol>(mol);
  clear_matrix_caches(*result);
  clear_descriptor_caches(*result);
  return result;
}

const MS::Tautomer &tautomer_at(const TautomerEnumeratorResult &result,
                                std::size_t i) {
  check_index(i, result.size());
  return std::next(result.smilesTautomerMap().begin(),
                   static_cast<std::ptrdiff_t>(i))
      ->second;
}

rust::Vec<std::uint32_t> indices(const boost::dynamic_bitset<> &bits) {
  rust::Vec<std::uint32_t> result;
  for (auto i = bits.find_first(); i != boost::dynamic_bitset<>::npos;
       i = bits.find_next(i)) {
    result.push_back(static_cast<std::uint32_t>(i));
  }
  return result;
}

// The enumerator's getters do not change it, but are not const.
TautomerEnumerator &mut(const TautomerEnumerator &enumerator) {
  return const_cast<TautomerEnumerator &>(enumerator);
}

}  // namespace

std::unique_ptr<TautomerEnumerator> tautomer_enumerator_default() {
  return std::make_unique<TautomerEnumerator>();
}

std::unique_ptr<TautomerEnumerator> tautomer_enumerator_new(
    const CleanupParametersData &params) {
  return checked_enumerator(to_cleanup_parameters(params));
}

std::unique_ptr<TautomerEnumerator> tautomer_enumerator_v1() {
  return std::unique_ptr<TautomerEnumerator>(MS::getV1TautomerEnumerator());
}

std::unique_ptr<TautomerEnumerator> tautomer_enumerator_clone(
    const TautomerEnumerator &enumerator) {
  return std::make_unique<TautomerEnumerator>(enumerator);
}

std::uint32_t tautomer_enumerator_get_max_tautomers(
    const TautomerEnumerator &enumerator) {
  return mut(enumerator).getMaxTautomers();
}

void tautomer_enumerator_set_max_tautomers(TautomerEnumerator &enumerator,
                                           std::uint32_t max_tautomers) {
  enumerator.setMaxTautomers(max_tautomers);
}

std::uint32_t tautomer_enumerator_get_max_transforms(
    const TautomerEnumerator &enumerator) {
  return mut(enumerator).getMaxTransforms();
}

void tautomer_enumerator_set_max_transforms(TautomerEnumerator &enumerator,
                                            std::uint32_t max_transforms) {
  enumerator.setMaxTransforms(max_transforms);
}

bool tautomer_enumerator_get_remove_sp3_stereo(
    const TautomerEnumerator &enumerator) {
  return mut(enumerator).getRemoveSp3Stereo();
}

void tautomer_enumerator_set_remove_sp3_stereo(TautomerEnumerator &enumerator,
                                               bool value) {
  enumerator.setRemoveSp3Stereo(value);
}

bool tautomer_enumerator_get_remove_bond_stereo(
    const TautomerEnumerator &enumerator) {
  return mut(enumerator).getRemoveBondStereo();
}

void tautomer_enumerator_set_remove_bond_stereo(TautomerEnumerator &enumerator,
                                                bool value) {
  enumerator.setRemoveBondStereo(value);
}

bool tautomer_enumerator_get_reassign_stereo(
    const TautomerEnumerator &enumerator) {
  return mut(enumerator).getReassignStereo();
}

void tautomer_enumerator_set_reassign_stereo(TautomerEnumerator &enumerator,
                                             bool value) {
  enumerator.setReassignStereo(value);
}

void tautomer_enumerator_set_callback(TautomerEnumerator &enumerator,
                                      const TautomerCallbackHandle &callback) {
  enumerator.setCallback(callback.ptr);
}

void tautomer_enumerator_clear_callback(TautomerEnumerator &enumerator) {
  enumerator.setCallback(std::shared_ptr<MS::TautomerEnumeratorCallback>());
}

// Enumerating, and the default scores, match substructures
// (control_c_lock.h); a Rust callback or score runs without the lock.
// Enumerating writes the molecule and its tautomers, copies of it, as
// SMILES, and pickCanonical writes the tautomers it ties, either of which
// removes atoms from copies of a molecule of several fragments
// (check_removals).
std::unique_ptr<TautomerEnumeratorResult> tautomer_enumerator_enumerate(
    const TautomerEnumerator &enumerator, const RWMol &mol) {
  check_removals(mol);
  HoldControlC control_c;
  return std::make_unique<TautomerEnumeratorResult>(enumerator.enumerate(mol));
}

// canonicalize scores tautomers that share the molecule's graph, and so its
// rings, which are counted on a copy.
std::unique_ptr<RWMol> tautomer_enumerator_canonicalize(
    const TautomerEnumerator &enumerator, const RWMol &mol) {
  check_tautomer_score(mol);
  check_removals(mol);
  HoldControlC control_c;
  return adopt_mol(enumerator.canonicalize(mol));
}

std::unique_ptr<RWMol> tautomer_enumerator_canonicalize_with(
    const TautomerEnumerator &enumerator, const RWMol &mol,
    const ScoreFunction &score) {
  check_removals(mol);
  HoldControlC control_c;
  return adopt_mol(enumerator.canonicalize(mol, scorer(score)));
}

namespace {

std::vector<RDKit::ROMOL_SPTR> &tautomers(BorrowedMols &mols) {
  if (mols.mols().empty()) {
    throw ValueErrorException("no tautomers to pick from");
  }
  for (const auto &mol : mols.mols()) {
    check_removals(*mol);
  }
  return mols.mols();
}

const TautomerEnumeratorResult &nonempty(
    const TautomerEnumeratorResult &result) {
  if (result.empty()) {
    throw ValueErrorException("no tautomers to pick from");
  }
  return result;
}

}  // namespace

// pickCanonical scores the tautomers only when there are several.
std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical(
    const TautomerEnumerator &enumerator, BorrowedMols &mols) {
  auto &list = tautomers(mols);
  if (list.size() > 1) {
    for (const auto &mol : list) {
      check_default_score(*mol);
    }
  }
  HoldControlC control_c;
  return adopt_mol(enumerator.pickCanonical(list));
}

std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical_with(
    const TautomerEnumerator &enumerator, BorrowedMols &mols,
    const ScoreFunction &score) {
  auto &list = tautomers(mols);
  HoldControlC control_c;
  return adopt_mol(enumerator.pickCanonical(list, scorer(score)));
}

std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical_result(
    const TautomerEnumerator &enumerator,
    const TautomerEnumeratorResult &result) {
  if (nonempty(result).size() > 1) {
    for (const auto &mol : result.tautomers()) {
      check_default_score(*mol);
    }
  }
  HoldControlC control_c;
  return adopt_mol(enumerator.pickCanonical(result));
}

std::unique_ptr<RWMol> tautomer_enumerator_pick_canonical_result_with(
    const TautomerEnumerator &enumerator,
    const TautomerEnumeratorResult &result, const ScoreFunction &score) {
  const auto &checked = nonempty(result);
  HoldControlC control_c;
  return adopt_mol(enumerator.pickCanonical(checked, scorer(score)));
}

std::unique_ptr<TautomerCallbackHandle> tautomer_callback_new(
    rust::Box<TautomerCallback> callback) {
  auto handle = std::make_unique<TautomerCallbackHandle>();
  handle->ptr = std::make_shared<RustTautomerCallback>(std::move(callback));
  return handle;
}

std::unique_ptr<TautomerCallbackHandle> tautomer_callback_clone(
    const TautomerCallbackHandle &handle) {
  return std::make_unique<TautomerCallbackHandle>(handle);
}

bool tautomer_callback_call(const TautomerCallbackHandle &handle,
                            const RWMol &mol,
                            const TautomerEnumeratorResult &result) {
  return (*handle.ptr)(mol, result);
}

std::size_t tautomer_result_size(const TautomerEnumeratorResult &result) {
  return result.size();
}

std::unique_ptr<RWMol> tautomer_result_tautomer(
    const TautomerEnumeratorResult &result, std::size_t i) {
  return copy_mol(*tautomer_at(result, i).tautomer);
}

std::unique_ptr<RWMol> tautomer_result_kekulized(
    const TautomerEnumeratorResult &result, std::size_t i) {
  const auto &tautomer = tautomer_at(result, i);
  if (tautomer.kekulized) {
    return copy_mol(*tautomer.kekulized);
  }
  // As Tautomer::getKekulized, whose cache would be written through const.
  auto kekulized = copy_mol(*tautomer.tautomer);
  RDKit::MolOps::Kekulize(*kekulized, false, true);
  return kekulized;
}

std::unique_ptr<std::vector<std::string>> tautomer_result_smiles(
    const TautomerEnumeratorResult &result) {
  return std::make_unique<std::vector<std::string>>(result.smiles());
}

TautomerEnumeratorStatus tautomer_result_status(
    const TautomerEnumeratorResult &result) {
  return result.status();
}

rust::Vec<std::uint32_t> tautomer_result_modified_atoms(
    const TautomerEnumeratorResult &result) {
  return indices(result.modifiedAtoms());
}

rust::Vec<std::uint32_t> tautomer_result_modified_bonds(
    const TautomerEnumeratorResult &result) {
  return indices(result.modifiedBonds());
}

// scoreRings adds 100 or 250 per ring in int, and the symmetrized SSSR's
// relevant cycles can be exponentially many.
std::int32_t score_rings(const RWMol &mol) {
  if (!mol.getRingInfo()->isSymmSssr()) {
    check_relevant_cycles(mol);
    RDKit::MolOps::symmetrizeSSSR(const_cast<RWMol &>(mol));
  }
  if (250 * static_cast<std::uint64_t>(mol.getRingInfo()->numRings()) >
      static_cast<std::uint64_t>(INT_MAX)) {
    throw ValueErrorException(
        "RDKit adds ring scores in int, which overflows for this many rings");
  }
  return TSF::scoreRings(mol);
}

// scoreHeteroHs subtracts unsigned counts, so its sum wraps but does not
// overflow.
std::int32_t score_hetero_hs(const RWMol &mol) {
  return TSF::scoreHeteroHs(mol);
}

// The default terms' scores times their at most 1000 matches stay in int.
std::int32_t score_substructs(const RWMol &mol) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return TSF::scoreSubstructs(mol);
}

// scoreSubstructs adds each term's score times its match count (at most
// SubstructMatchParameters' 1000) in int.
std::int32_t score_substructs_with_terms(
    const RWMol &mol, rust::Slice<const SubstructTermData> data) {
  std::vector<TSF::SubstructTerm> terms;
  std::int64_t bound = 0;
  for (const auto &term : data) {
    bound += 1000 * std::llabs(static_cast<long long>(term.score));
    if (bound > INT_MAX) {
      throw ValueErrorException(
          "RDKit adds substructure scores times their match counts (up to "
          "1000) in int, which overflows for these terms");
    }
    terms.emplace_back(to_std(rust::Slice<const std::uint8_t>(
                           term.name.data(), term.name.size())),
                       to_std(rust::Slice<const std::uint8_t>(
                           term.smarts.data(), term.smarts.size())),
                       term.score);
  }
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return TSF::scoreSubstructs(mol, terms);
}

std::int32_t score_tautomer(const RWMol &mol) {
  check_default_score(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return TSF::scoreTautomer(mol);
}

rust::Vec<SubstructTermData> default_tautomer_score_substructs() {
  rust::Vec<SubstructTermData> result;
  for (const auto &term : TSF::getDefaultTautomerScoreSubstructs()) {
    result.push_back(SubstructTermData{to_bytes(term.name),
                                       to_bytes(term.smarts), term.score});
  }
  return result;
}

std::unique_ptr<std::string> tautomer_score_version() {
  return std::make_unique<std::string>(TSF::tautomerScoringVersion);
}

}  // namespace rdkit_cxx
