// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/fmcs.h"

#include <chrono>
#include <cmath>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/fmcs.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

MCSAtomCompareParametersData atom_params_to_data(
    const RDKit::MCSAtomCompareParameters &p) {
  MCSAtomCompareParametersData data{};
  data.match_valences = p.MatchValences;
  data.match_chiral_tag = p.MatchChiralTag;
  data.match_formal_charge = p.MatchFormalCharge;
  data.ring_matches_ring_only = p.RingMatchesRingOnly;
  data.complete_rings_only = p.CompleteRingsOnly;
  data.match_isotope = p.MatchIsotope;
  data.max_distance = p.MaxDistance;
  return data;
}

RDKit::MCSAtomCompareParameters atom_params_from_data(
    const MCSAtomCompareParametersData &data) {
  RDKit::MCSAtomCompareParameters p;
  p.MatchValences = data.match_valences;
  p.MatchChiralTag = data.match_chiral_tag;
  p.MatchFormalCharge = data.match_formal_charge;
  p.RingMatchesRingOnly = data.ring_matches_ring_only;
  p.CompleteRingsOnly = data.complete_rings_only;
  p.MatchIsotope = data.match_isotope;
  p.MaxDistance = data.max_distance;
  return p;
}

MCSBondCompareParametersData bond_params_to_data(
    const RDKit::MCSBondCompareParameters &p) {
  MCSBondCompareParametersData data{};
  data.ring_matches_ring_only = p.RingMatchesRingOnly;
  data.complete_rings_only = p.CompleteRingsOnly;
  data.match_fused_rings = p.MatchFusedRings;
  data.match_fused_rings_strict = p.MatchFusedRingsStrict;
  data.match_stereo = p.MatchStereo;
  return data;
}

RDKit::MCSBondCompareParameters bond_params_from_data(
    const MCSBondCompareParametersData &data) {
  RDKit::MCSBondCompareParameters p;
  p.RingMatchesRingOnly = data.ring_matches_ring_only;
  p.CompleteRingsOnly = data.complete_rings_only;
  p.MatchFusedRings = data.match_fused_rings;
  p.MatchFusedRingsStrict = data.match_fused_rings_strict;
  p.MatchStereo = data.match_stereo;
  return p;
}

// The comparator for one of RDKit's typers, as Python's getMCSAtomTyper
// maps them; false for any other function (a trampoline).
bool predefined_atom_typer(RDKit::MCSAtomCompareFunction typer,
                           AtomComparator &comparator) {
  if (typer == RDKit::MCSAtomCompareAny) {
    comparator = RDKit::AtomCompareAny;
  } else if (typer == RDKit::MCSAtomCompareElements) {
    comparator = RDKit::AtomCompareElements;
  } else if (typer == RDKit::MCSAtomCompareIsotopes) {
    comparator = RDKit::AtomCompareIsotopes;
  } else if (typer == RDKit::MCSAtomCompareAnyHeavyAtom) {
    comparator = RDKit::AtomCompareAnyHeavyAtom;
  } else {
    return false;
  }
  return true;
}

bool predefined_bond_typer(RDKit::MCSBondCompareFunction typer,
                           BondComparator &comparator) {
  if (typer == RDKit::MCSBondCompareAny) {
    comparator = RDKit::BondCompareAny;
  } else if (typer == RDKit::MCSBondCompareOrder) {
    comparator = RDKit::BondCompareOrder;
  } else if (typer == RDKit::MCSBondCompareOrderExact) {
    comparator = RDKit::BondCompareOrderExact;
  } else {
    return false;
  }
  return true;
}

MCSParametersData params_to_data(const RDKit::MCSParameters &p) {
  MCSParametersData data{};
  data.maximize_bonds = p.MaximizeBonds;
  data.threshold = p.Threshold;
  data.timeout = p.Timeout;
  data.verbose = p.Verbose;
  data.atom_compare_parameters = atom_params_to_data(p.AtomCompareParameters);
  data.bond_compare_parameters = bond_params_to_data(p.BondCompareParameters);
  data.atom_typer = RDKit::AtomCompareElements;
  data.custom_atom_typer = !predefined_atom_typer(p.AtomTyper, data.atom_typer);
  data.bond_typer = RDKit::BondCompareOrder;
  data.custom_bond_typer = !predefined_bond_typer(p.BondTyper, data.bond_typer);
  data.progress_callback = p.ProgressCallback != nullptr;
  data.final_match_checker = p.FinalMatchChecker != nullptr;
  data.should_accept_mcs = p.ShouldAcceptMCS != nullptr;
  // The seed is RDKit's empty default or the UTF-8 Rust set, so this
  // conversion does not throw.
  data.initial_seed = rust::String(p.InitialSeed);
  data.store_all = p.StoreAll;
  return data;
}

// The state the trampolines share for one call, which RDKit's user-data
// pointers point to.
struct FmcsState {
  explicit FmcsState(const McsCallbacks &callbacks) : callbacks(callbacks) {}

  const McsCallbacks &callbacks;
  // Each molecule's index in the call, which is how Rust finds it.
  std::unordered_map<const RDKit::ROMol *, std::size_t> indices;
  // Whether Rust set a progress callback, or the trampoline stands in for
  // RDKit's timeout check.
  bool user_progress = false;
  std::chrono::steady_clock::time_point start;
};

std::size_t mol_index(const FmcsState &state, const RDKit::ROMol &mol) {
  auto it = state.indices.find(&mol);
  if (it == state.indices.end()) {
    throw std::runtime_error(
        "FindMCS called back with a molecule it was not given");
  }
  return it->second;
}

// Runs a trampoline's body. Nothing throws back into RDKit's search: an
// exception is recorded for Rust to return once FindMCS is done, and the
// callback answers false, as every callback does once one has failed (the
// progress trampoline's false cancels the search).
template <typename Body>
bool guarded(const FmcsState &state, Body &&body) {
  if (mcs_callbacks_failed(state.callbacks)) {
    return false;
  }
  try {
    return body();
  } catch (...) {
    mcs_callbacks_fail(state.callbacks, encode_current_exception());
    return false;
  }
}

const FmcsState &state_of(const void *user_data) {
  return *static_cast<const FmcsState *>(user_data);
}

// Calls a closure RDKit calls with no ControlCHandler alive (the typers,
// from MaximumCommonSubgraph::init; the progress callback, between seeds;
// the acceptance check, after its matches), with control_c_mutex released
// (ReleaseControlC) so the closure may call wrappers that take it.
template <typename Call>
bool unlocked(Call &&call) {
  ReleaseControlC release;
  return call();
}

// Marks a closure RDKit may call from inside vf2, with a ControlCHandler
// alive, so the lock stays held: wrappers the closure calls that would take
// it fail (control_c_lock.h).
class InControlCCallback {
 public:
  InControlCCallback() : d_outer(in_control_c_callback) {
    in_control_c_callback = true;
  }
  ~InControlCCallback() { in_control_c_callback = d_outer; }
  InControlCCallback(const InControlCCallback &) = delete;
  InControlCCallback &operator=(const InControlCCallback &) = delete;

 private:
  bool d_outer;
};

rust::Slice<const IndexPair> slice(const std::vector<IndexPair> &pairs) {
  return rust::Slice<const IndexPair>(pairs.data(), pairs.size());
}

bool atom_typer_trampoline(const RDKit::MCSAtomCompareParameters &p,
                           const RDKit::ROMol &mol1, unsigned int atom1,
                           const RDKit::ROMol &mol2, unsigned int atom2,
                           void *user_data) {
  const auto &state = state_of(user_data);
  return guarded(state, [&] {
    auto data = atom_params_to_data(p);
    auto index1 = mol_index(state, mol1);
    auto index2 = mol_index(state, mol2);
    return unlocked([&] {
      return mcs_call_atom_typer(state.callbacks, data, index1, atom1, index2,
                                 atom2);
    });
  });
}

bool bond_typer_trampoline(const RDKit::MCSBondCompareParameters &p,
                           const RDKit::ROMol &mol1, unsigned int bond1,
                           const RDKit::ROMol &mol2, unsigned int bond2,
                           void *user_data) {
  const auto &state = state_of(user_data);
  return guarded(state, [&] {
    auto data = bond_params_to_data(p);
    auto index1 = mol_index(state, mol1);
    auto index2 = mol_index(state, mol2);
    return unlocked([&] {
      return mcs_call_bond_typer(state.callbacks, data, index1, bond1, index2,
                                 bond2);
    });
  });
}

bool progress_trampoline(const RDKit::MCSProgressData &stat,
                         const RDKit::MCSParameters &params, void *user_data) {
  const auto &state = state_of(user_data);
  return guarded(state, [&] {
    if (!state.user_progress) {
      // RDKit's MCSProgressCallbackTimeout, which it installs only when no
      // callback is set: this one is there to cancel the search once a
      // closure failed, so it keeps the timeout too.
      auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
                         std::chrono::steady_clock::now() - state.start)
                         .count();
      return !params.Timeout || static_cast<unsigned long long>(elapsed) <=
                                    params.Timeout * 1000000ULL;
    }
    auto data = params_to_data(params);
    return unlocked([&] {
      return mcs_call_progress(state.callbacks, stat.NumAtoms, stat.NumBonds,
                               stat.SeedProcessed, data);
    });
  });
}

// As MCSFinalMatchCheckPyFunc in rdFMCS.cpp, with its
// buildAtomIdxMatchTuple and buildBondIdxMatchTuple.
bool final_match_check_trampoline(const std::uint32_t c1[],
                                  const std::uint32_t c2[],
                                  const RDKit::ROMol &mol1,
                                  const RDKit::FMCS::Graph &query,
                                  const RDKit::ROMol &mol2,
                                  const RDKit::FMCS::Graph &target,
                                  const RDKit::MCSParameters *params) {
  if (!params) {
    return false;
  }
  const auto &state = state_of(params->FinalMatchCheckerUserData);
  return guarded(state, [&] {
    std::vector<IndexPair> atoms;
    auto num_atoms = boost::num_vertices(query);
    for (unsigned int i = 0; i < num_atoms; ++i) {
      atoms.push_back(IndexPair{query[c1[i]], target[c2[i]]});
    }
    std::vector<IndexPair> bonds;
    auto [edge, end] = boost::edges(query);
    for (; edge != end; ++edge) {
      auto source = boost::source(*edge, query);
      auto dest = boost::target(*edge, query);
      const auto *query_bond =
          mol1.getBondBetweenAtoms(query[c1[source]], query[c1[dest]]);
      const auto *target_bond =
          mol2.getBondBetweenAtoms(target[c2[source]], target[c2[dest]]);
      if (!query_bond || !target_bond) {
        throw std::runtime_error("FindMCS matched a bond that is not there");
      }
      bonds.push_back(IndexPair{query_bond->getIdx(), target_bond->getIdx()});
    }
    auto index1 = mol_index(state, mol1);
    auto index2 = mol_index(state, mol2);
    auto data = params_to_data(*params);
    InControlCCallback in_callback;
    return mcs_call_final_match_checker(state.callbacks, index1, index2,
                                        slice(atoms), slice(bonds), data);
  });
}

std::vector<IndexPair> index_pairs(
    const std::vector<std::pair<int, int>> &matches) {
  std::vector<IndexPair> pairs;
  pairs.reserve(matches.size());
  for (const auto &[first, second] : matches) {
    pairs.push_back(IndexPair{static_cast<std::uint32_t>(first),
                              static_cast<std::uint32_t>(second)});
  }
  return pairs;
}

bool should_accept_trampoline(
    const RDKit::ROMol &query, const RDKit::ROMol &target,
    const std::vector<std::pair<int, int>> &atom_idx_match,
    const std::vector<std::pair<int, int>> &bond_idx_match,
    const RDKit::MCSParameters *params) {
  if (!params) {
    return false;
  }
  const auto &state = state_of(params->ShouldAcceptMCSUserData);
  return guarded(state, [&] {
    auto atoms = index_pairs(atom_idx_match);
    auto bonds = index_pairs(bond_idx_match);
    auto query_index = mol_index(state, query);
    auto target_index = mol_index(state, target);
    auto data = params_to_data(*params);
    return unlocked([&] {
      return mcs_call_should_accept_mcs(state.callbacks, query_index,
                                        target_index, slice(atoms),
                                        slice(bonds), data);
    });
  });
}

// findMCS gives molecules without ring info an empty one for the search
// and resets it before returning, but not when it throws (a callback's
// error is recorded, not thrown, but RDKit's own checks throw: a missing
// conformer for MaxDistance, valences it has not computed), which would
// leave such a molecule reporting that it has no rings. This resets them
// either way.
class RingInfoGuard {
 public:
  explicit RingInfoGuard(const std::vector<RDKit::ROMOL_SPTR> &mols) {
    for (const auto &mol : mols) {
      if (!mol->getRingInfo()->isInitialized()) {
        d_mols.push_back(mol.get());
      }
    }
  }
  RingInfoGuard(const RingInfoGuard &) = delete;
  RingInfoGuard &operator=(const RingInfoGuard &) = delete;
  ~RingInfoGuard() {
    for (auto *mol : d_mols) {
      mol->getRingInfo()->reset();
    }
  }

 private:
  std::vector<RDKit::ROMol *> d_mols;
};

std::unique_ptr<RWMol> copy_out(const RDKit::ROMol &mol) {
  auto copy = std::make_unique<RWMol>(mol);
  clear_matrix_caches(*copy);
  return copy;
}

void check_index_args(std::uint32_t index1, unsigned int len1,
                      std::uint32_t index2, unsigned int len2) {
  check_index(index1, len1);
  check_index(index2, len2);
}

}  // namespace

MCSAtomCompareParametersData mcs_atom_compare_parameters_default() {
  return atom_params_to_data(RDKit::MCSAtomCompareParameters());
}

MCSBondCompareParametersData mcs_bond_compare_parameters_default() {
  return bond_params_to_data(RDKit::MCSBondCompareParameters());
}

MCSParametersData mcs_parameters_default() {
  return params_to_data(RDKit::MCSParameters());
}

MCSResultData find_mcs(BorrowedMols &mols, const MCSParametersData &data,
                       const McsCallbacks &callbacks) {
  RDKit::MCSParameters params;
  params.MaximizeBonds = data.maximize_bonds;
  params.Threshold = data.threshold;
  params.Timeout = data.timeout;
  params.Verbose = data.verbose;
  params.AtomCompareParameters =
      atom_params_from_data(data.atom_compare_parameters);
  params.BondCompareParameters =
      bond_params_from_data(data.bond_compare_parameters);
  params.setMCSAtomTyperFromEnum(data.atom_typer);
  params.setMCSBondTyperFromEnum(data.bond_typer);
  params.InitialSeed = std::string(data.initial_seed);
  params.StoreAll = data.store_all;
  // RDKit 2024.09's MaximumCommonSubgraph::find matches every target against
  // res.QueryMol when Verbose is set, but with StoreAll it never sets
  // QueryMol, so that is a null dereference (rdkit/rdkit#9040). (Its verbose
  // report of each MCS found reads McsIdx.Bonds.front(), which no MCS of one
  // atom reaches: growSeeds stops before growing when the seeds' bonds match
  // the MCS's, both 0 for a seed of one atom.)
  if (params.Verbose && params.StoreAll) {
    throw ValueErrorException(
        "Verbose with StoreAll dereferences a null query molecule in RDKit "
        "2024.09; set only one");
  }

  // MaximumCommonSubgraph::find converts ceil(n * Threshold) to int, which
  // is undefined for NaN (Python passes it through) and for values past
  // INT_MIN. Every threshold at or below zero asks for one matching target,
  // as zero does, so clamping at -1 keeps RDKit's result in range.
  if (std::isnan(params.Threshold)) {
    throw ValueErrorException("Threshold must not be NaN");
  }
  if (params.Threshold < -1.0) {
    params.Threshold = -1.0;
  }
  // makeInitialSeeds dereferences SmartsToMol's result without a check, so
  // a seed that does not parse would be a null dereference (Python
  // crashes); it is parsed here first, the same way.
  if (!params.InitialSeed.empty()) {
    check_smiles_cx(params.InitialSeed, true, true);
    std::unique_ptr<RDKit::RWMol> seed(RDKit::SmartsToMol(params.InitialSeed));
    if (!seed) {
      throw RDKit::SmilesParse::SmilesParseException(
          "InitialSeed is not valid SMARTS: " + params.InitialSeed);
    }
  }

  // MCSBondCompareOrder and MCSBondCompareOrderExact index a table by both
  // bonds' types, through ZERO; a molecule pickle can carry any type up to
  // 127, which would be read past the table. Unpickling now refuses one
  // (check_pickled_mol), so this checks again only in case another path
  // sets one.
  if (params.BondTyper == RDKit::MCSBondCompareOrder ||
      params.BondTyper == RDKit::MCSBondCompareOrderExact) {
    for (const auto &mol : mols.mols()) {
      for (const auto *bond : mol->bonds()) {
        if (bond->getBondType() > RDKit::Bond::ZERO) {
          throw ValueErrorException(
              "bond " + std::to_string(bond->getIdx()) +
              " has a type the bond comparison does not know");
        }
      }
    }
  }

  FmcsState state(callbacks);
  for (std::size_t i = 0; i < mols.mols().size(); ++i) {
    state.indices.emplace(mols.mols()[i].get(), i);
  }
  // As Python's MCSParameters sets them: the comparison user data only for
  // a custom typer (RDKit's init reads its presence).
  if (data.custom_atom_typer) {
    params.AtomTyper = atom_typer_trampoline;
    params.CompareFunctionsUserData = &state;
  }
  if (data.custom_bond_typer) {
    params.BondTyper = bond_typer_trampoline;
    params.CompareFunctionsUserData = &state;
  }
  if (data.final_match_checker) {
    params.FinalMatchChecker = final_match_check_trampoline;
    params.FinalMatchCheckerUserData = &state;
  }
  if (data.should_accept_mcs) {
    params.ShouldAcceptMCS = should_accept_trampoline;
    params.ShouldAcceptMCSUserData = &state;
  }
  if (data.custom_atom_typer || data.custom_bond_typer ||
      data.progress_callback || data.final_match_checker ||
      data.should_accept_mcs) {
    // With a progress callback of its own, RDKit enforces no timeout, as
    // in Python; otherwise the trampoline checks it.
    params.ProgressCallback = progress_trampoline;
    params.ProgressCallbackUserData = &state;
    state.user_progress = data.progress_callback;
  }

  // The final checks match in vf2, and the seed with SubstructMatch, each
  // with a ControlCHandler (control_c_lock.h); the trampolines release the
  // lock where none is alive.
  RDKit::MCSResult result;
  {
    HoldControlC control_c;
    RingInfoGuard guard(mols.mols());
    state.start = std::chrono::steady_clock::now();
    result = RDKit::findMCS(mols.mols(), &params);
  }

  MCSResultData out{};
  out.num_atoms = result.NumAtoms;
  out.num_bonds = result.NumBonds;
  out.smarts_string = std::make_unique<std::string>(result.SmartsString);
  out.canceled = result.Canceled;
  if (result.QueryMol) {
    out.query_mol = copy_out(*result.QueryMol);
  }
  out.degenerate_smarts_query_mol_dict = std::make_unique<MolList>();
  for (const auto &[smarts, mol] : result.DegenerateSmartsQueryMolDict) {
    if (mol) {
      out.degenerate_smarts_query_mol_dict->push_named(smarts, copy_out(*mol));
    }
  }
  return out;
}

// The checks take indices RDKit reads unchecked through ring info (or
// checks with a RuntimeError in getAtomWithIdx and getBondWithIdx), so each
// is checked first (E3).
bool mcs_atom_compare_check_atom_ring_match(
    const MCSAtomCompareParametersData &params, const RWMol &mol1,
    std::uint32_t atom1, const RWMol &mol2, std::uint32_t atom2) {
  check_index_args(atom1, mol1.getNumAtoms(), atom2, mol2.getNumAtoms());
  return RDKit::checkAtomRingMatch(atom_params_from_data(params), mol1, atom1,
                                   mol2, atom2);
}

bool mcs_atom_compare_check_atom_charge(
    const MCSAtomCompareParametersData &params, const RWMol &mol1,
    std::uint32_t atom1, const RWMol &mol2, std::uint32_t atom2) {
  check_index_args(atom1, mol1.getNumAtoms(), atom2, mol2.getNumAtoms());
  return RDKit::checkAtomCharge(atom_params_from_data(params), mol1, atom1,
                                mol2, atom2);
}

bool mcs_atom_compare_check_atom_chirality(
    const MCSAtomCompareParametersData &params, const RWMol &mol1,
    std::uint32_t atom1, const RWMol &mol2, std::uint32_t atom2) {
  check_index_args(atom1, mol1.getNumAtoms(), atom2, mol2.getNumAtoms());
  return RDKit::checkAtomChirality(atom_params_from_data(params), mol1, atom1,
                                   mol2, atom2);
}

bool mcs_bond_compare_check_bond_stereo(
    const MCSBondCompareParametersData &params, const RWMol &mol1,
    std::uint32_t bond1, const RWMol &mol2, std::uint32_t bond2) {
  check_index_args(bond1, mol1.getNumBonds(), bond2, mol2.getNumBonds());
  return RDKit::checkBondStereo(bond_params_from_data(params), mol1, bond1,
                                mol2, bond2);
}

bool mcs_bond_compare_check_bond_ring_match(
    const MCSBondCompareParametersData &params, const RWMol &mol1,
    std::uint32_t bond1, const RWMol &mol2, std::uint32_t bond2) {
  check_index_args(bond1, mol1.getNumBonds(), bond2, mol2.getNumBonds());
  return RDKit::checkBondRingMatch(bond_params_from_data(params), mol1, bond1,
                                   mol2, bond2);
}

}  // namespace rdkit_cxx
