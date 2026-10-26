// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_align.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <Geometry/Transform3D.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolAlign/O3AAlignMolecules.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Vector.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/mol_align.rs.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/symmetrize_lock.h"
#include "rdkit-cxx/threads.h"

namespace rdkit_cxx {

namespace {

namespace MolAlign = RDKit::MolAlign;
using RDKit::MatchVectType;
using RDKit::ROMol;

rust::Vec<double> to_values(const RDGeom::Transform3D &trans) {
  rust::Vec<double> values;
  const double *data = trans.getData();
  for (unsigned i = 0; i < 16; ++i) {
    values.push_back(data[i]);
  }
  return values;
}

// Each pair is (probe atom, reference atom); RDKit reads their positions
// with a range check that, on RDKit's threads, calls std::terminate (E3).
MatchVectType to_match(const IndexPair *pairs, std::size_t n, const ROMol &prb,
                       const ROMol &ref) {
  MatchVectType match;
  for (std::size_t i = 0; i < n; ++i) {
    check_index(pairs[i].first, prb.getNumAtoms());
    check_index(pairs[i].second, ref.getNumAtoms());
    match.emplace_back(static_cast<int>(pairs[i].first),
                       static_cast<int>(pairs[i].second));
  }
  return match;
}

// As translateDoubleSeq in Code/RDBoost/Wrap.cpp, which
// Code/GraphMol/MolAlign/Wrap/rdMolAlign.cpp uses: no values are no weights.
std::unique_ptr<RDNumeric::DoubleVector> to_vector(
    rust::Slice<const double> values) {
  if (values.empty()) {
    return nullptr;
  }
  auto vec = std::make_unique<RDNumeric::DoubleVector>(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    vec->setVal(i, values[i]);
  }
  return vec;
}

// RDKit's own threads (getBestRMSInternal, getAllConformerBestRMS,
// getAllConformerBestRMSToRef) call std::terminate on any exception,
// getNumThreadsToUse negates its argument in int, which is undefined for
// INT_MIN, and nothing bounds the count. The count is computed here, at
// most the hardware's.
int thread_count(std::int32_t target) {
  unsigned n = num_threads_to_use(target);
  unsigned hardware = std::max(1u, std::thread::hardware_concurrency());
  return static_cast<int>(std::min(n, hardware));
}

// As getAllMatchesPrbRef in AlignMolecules.cpp, which RDKit calls when no
// map is given. symmetrizeTerminalAtoms writes a static (symmetrize_mutex),
// so the matches are found here, under the lock only for that, and passed
// to RDKit as the map. Both steps match substructures, so both hold
// control_c_mutex (control_c_lock.h), taken first.
std::vector<MatchVectType> all_matches(const ROMol &prb, const ROMol &ref,
                                       int max_matches, bool symmetrize,
                                       bool ignore_hs) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::unique_ptr<RDKit::RWMol> symmetrized;
  if (symmetrize) {
    symmetrized = std::make_unique<RDKit::RWMol>(prb);
    // symmetrizeTerminalAtoms replaces bonds with RWMol::replaceBond, which
    // leaves stereo groups pointing at the bonds it deletes (replace_bond);
    // the copy is only matched, without its groups' bonds.
    clear_stereo_group_bonds(*symmetrized);
    std::lock_guard<std::mutex> lock(symmetrize_mutex());
    MolAlign::details::symmetrizeTerminalAtoms(*symmetrized);
  }
  const ROMol &for_match = symmetrized ? *symmetrized : prb;
  std::vector<MatchVectType> matches;
  RDKit::SubstructMatch(ref, for_match, matches, false, true, false, false,
                        max_matches);
  if (matches.empty()) {
    throw MolAlign::MolAlignException(
        "No sub-structure match found between the reference and probe mol");
  }
  if (matches.size() > 1e6) {
    std::string name;
    prb.getPropIfPresent(RDKit::common_properties::_Name, name);
    BOOST_LOG(rdWarningLog)
        << "Warning in getAllMatchesPrbRef: " << matches.size()
        << " matches detected for molecule " << name << ", this may "
        << "lead to a performance slowdown." << std::endl;
  }
  if (ignore_hs) {
    for (auto &match : matches) {
      std::erase_if(match, [&for_match](const auto &mi) {
        return for_match.getAtomWithIdx(mi.first)->getAtomicNum() == 1;
      });
    }
  }
  return matches;
}

// AlignPoints and calcMSDInternal check the number of weights against each
// match with a precondition, which on RDKit's threads calls std::terminate.
void check_weights(const std::vector<MatchVectType> &matches,
                   const RDNumeric::DoubleVector *weights) {
  if (!weights) {
    return;
  }
  for (const auto &match : matches) {
    if (match.size() != weights->size()) {
      throw ValueErrorException("Mismatch between number of weights (" +
                                std::to_string(weights->size()) +
                                ") and number of atoms in the alignment map (" +
                                std::to_string(match.size()) + ")");
    }
  }
}

// MolAlign::BestAlignmentParams, with its weights, and its map checked or
// found here (all_matches).
struct NativeBestParams {
  std::unique_ptr<RDNumeric::DoubleVector> weights;
  MolAlign::BestAlignmentParams params;
};

NativeBestParams to_native(const BestAlignmentParamsData &data,
                           const ROMol &prb, const ROMol &ref, bool ignore_hs) {
  std::size_t total = 0;
  for (auto count : data.map_counts) {
    total += count;
  }
  if (total != data.map_pairs.size()) {
    throw ValueErrorException("the atom maps' lengths do not add up");
  }
  NativeBestParams native;
  std::vector<MatchVectType> map;
  std::size_t at = 0;
  for (auto count : data.map_counts) {
    map.push_back(to_match(data.map_pairs.data() + at, count, prb, ref));
    at += count;
  }
  native.weights = to_vector(
      rust::Slice<const double>(data.weights.data(), data.weights.size()));
  if (map.empty()) {
    map = all_matches(prb, ref, data.max_matches,
                      data.symmetrize_conjugated_terminal_groups, ignore_hs);
  }
  check_weights(map, native.weights.get());
  native.params.maxMatches = data.max_matches;
  native.params.symmetrizeConjugatedTerminalGroups =
      data.symmetrize_conjugated_terminal_groups;
  native.params.ignoreHs = data.ignore_hs;
  native.params.numThreads = thread_count(data.num_threads);
  native.params.map = std::move(map);
  native.params.weights = native.weights.get();
  return native;
}

// Open3DAlign trusts more of its input than the alignments do.
// MolHistogram and LAP::computeCostMatrix (O3AAlignMolecules.cpp) count
// atoms whose atomic number is not 1 as heavy, but size their tables with
// getNumHeavyAtoms, which counts those above 1, so a dummy atom makes them
// write past the tables' ends; Python's heap breaks.
void check_o3a_atoms(const ROMol &mol) {
  check_matrix_size(mol);
  for (const auto atom : mol.atoms()) {
    if (atom->getAtomicNum() < 1) {
      throw ValueErrorException(
          "O3A cannot align molecules with dummy atoms: it writes past its "
          "tables for them");
    }
  }
}

// MolHistogram casts each distance to unsigned int, which is undefined
// unless it is finite and below 2^32; coordinates within 1e9 keep every
// distance below 3.5e9. Its sorts compare distances, which NaN breaks.
void check_o3a_conformer(const ROMol &mol, int cid) {
  for (const auto &pos : mol.getConformer(cid).getPositions()) {
    for (double v : {pos.x, pos.y, pos.z}) {
      if (!(std::abs(v) <= 1e9)) {
        throw ValueErrorException("O3A needs finite coordinates within 1e9");
      }
    }
  }
}

void check_o3a_conformers(const ROMol &prb, const ROMol &ref,
                          bool for_probe_confs, const O3AOptionsData &opts) {
  check_o3a_atoms(prb);
  check_o3a_atoms(ref);
  if (for_probe_confs) {
    for (auto it = prb.beginConformers(); it != prb.endConformers(); ++it) {
      check_o3a_conformer(prb, (*it)->getId());
    }
  } else {
    check_o3a_conformer(prb, opts.prb_cid);
  }
  check_o3a_conformer(ref, opts.ref_cid);
}

// As parseConstraints in rdMolAlign.cpp. The constraints are sorted by
// weight (O3AConstraintVect::append, compareSDMDist), which NaN breaks.
struct Constraints {
  std::unique_ptr<MatchVectType> map;
  std::unique_ptr<RDNumeric::DoubleVector> weights;
};

Constraints to_constraints(const O3AOptionsData &opts, const ROMol &prb,
                           const ROMol &ref) {
  Constraints result;
  if (opts.constraint_map.empty()) {
    return result;
  }
  result.map = std::make_unique<MatchVectType>(to_match(
      opts.constraint_map.data(), opts.constraint_map.size(), prb, ref));
  if (opts.has_constraint_weights) {
    result.weights = to_vector(rust::Slice<const double>(
        opts.constraint_weights.data(), opts.constraint_weights.size()));
  }
  if (result.weights) {
    if (result.weights->size() != result.map->size()) {
      throw ValueErrorException(
          "The number of weights should match the number of constraints");
    }
    for (std::size_t i = 0; i < result.weights->size(); ++i) {
      if (std::isnan((*result.weights)[i])) {
        throw ValueErrorException("constraint weights must not be NaN");
      }
    }
  }
  for (const auto &[p, r] : *result.map) {
    if (prb.getAtomWithIdx(p)->getAtomicNum() == 1 ||
        ref.getAtomWithIdx(r)->getAtomicNum() == 1) {
      throw ValueErrorException("Constrained atoms must be heavy atoms");
    }
  }
  return result;
}

std::unique_ptr<O3AList> run_o3a(RWMol &prb, const RWMol &ref, void *prb_prop,
                                 void *ref_prop,
                                 MolAlign::O3A::AtomTypeScheme scheme,
                                 bool for_probe_confs,
                                 const O3AOptionsData &opts,
                                 const Constraints &constraints) {
  auto list = std::make_unique<O3AList>();
  if (for_probe_confs) {
    // RDKit's threads here run through std::async, which hands their
    // exceptions back to this thread.
    MolAlign::getO3AForProbeConfs(prb, ref, prb_prop, ref_prop, list->items,
                                  thread_count(opts.num_threads), scheme,
                                  opts.ref_cid, opts.reflect, opts.max_iters,
                                  opts.options, constraints.map.get(),
                                  constraints.weights.get());
  } else {
    list->items.push_back(boost::make_shared<MolAlign::O3A>(
        prb, ref, prb_prop, ref_prop, scheme, opts.prb_cid, opts.ref_cid,
        opts.reflect, opts.max_iters, opts.options, constraints.map.get(),
        constraints.weights.get()));
  }
  return list;
}

// As getMMFFO3A in rdMolAlign.cpp, which checks only the properties it
// computes. The cost functions index MMFF's similarity matrix by atom type
// minus 1, so an untyped atom (type 0) reads before it, and read the
// properties by atom index. Missing properties are computed on a copy, so
// MMFF's sanitization leaves the molecule as it is.
const RDKit::MMFF::MMFFMolProperties *mmff_props(
    const RWMol &mol, const MMFFMolProperties *given,
    std::unique_ptr<RDKit::MMFF::MMFFMolProperties> &owned, const char *which) {
  const RDKit::MMFF::MMFFMolProperties *props;
  if (given) {
    if (given->num_atoms() != mol.getNumAtoms()) {
      throw ValueErrorException(std::string("the MMFF properties for the ") +
                                which +
                                " molecule have another number of "
                                "atoms");
    }
    props = &given->props();
  } else {
    RDKit::RWMol copy(mol);
    owned = std::make_unique<RDKit::MMFF::MMFFMolProperties>(copy);
    props = owned.get();
  }
  if (!props->isValid()) {
    throw ValueErrorException(std::string("missing MMFF94 parameters for ") +
                              which + " molecule");
  }
  return props;
}

// As getCrippenO3A in rdMolAlign.cpp, which uses the given contributions
// only when there is one per atom, and computes them otherwise. The cost
// function rounds them, scaled, to int, and the assignment solver sums those
// costs in int, so given ones stay within the range RDKit's span.
std::vector<double> crippen_contribs(const RWMol &mol, bool has_given,
                                     rust::Slice<const double> given) {
  unsigned n = mol.getNumAtoms();
  if (has_given && given.size() == n) {
    for (double v : given) {
      if (!(std::abs(v) <= 100.0)) {
        throw ValueErrorException(
            "Crippen contributions must be finite and within 100");
      }
    }
    return std::vector<double>(given.begin(), given.end());
  }
  std::vector<double> logp(n), mr(n);
  std::vector<unsigned int> types(n);
  std::vector<std::string> labels(n);
  // getCrippenAtomContribs caches them in the molecule's properties (U8),
  // and types the atoms by substructure matches (control_c_lock.h).
  DescriptorCacheGuard guard(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  RDKit::Descriptors::getCrippenAtomContribs(mol, logp, mr, true, &types,
                                             &labels);
  return logp;
}

}  // namespace

BestAlignmentParamsData best_alignment_params_default() {
  MolAlign::BestAlignmentParams params;
  BestAlignmentParamsData data{};
  data.max_matches = params.maxMatches;
  data.symmetrize_conjugated_terminal_groups =
      params.symmetrizeConjugatedTerminalGroups;
  data.ignore_hs = params.ignoreHs;
  data.num_threads = params.numThreads;
  data.prb_cid = -1;
  data.ref_cid = -1;
  data.reflect = false;
  data.max_iters = 50;
  return data;
}

// As getMolAlignTransform in rdMolAlign.cpp, which checks the number of
// weights.
double get_alignment_transform(const RWMol &prb_mol, const RWMol &ref_mol,
                               std::int32_t prb_cid, std::int32_t ref_cid,
                               rust::Slice<const IndexPair> atom_map,
                               rust::Slice<const double> weights, bool reflect,
                               std::uint32_t max_iters,
                               rust::Vec<double> &trans) {
  auto map = to_match(atom_map.data(), atom_map.size(), prb_mol, ref_mol);
  auto weights_vec = to_vector(weights);
  if (weights_vec) {
    unsigned n = map.empty() ? prb_mol.getNumAtoms() : map.size();
    if (weights_vec->size() != n) {
      throw ValueErrorException("Incorrect number of weights specified");
    }
  }
  RDGeom::Transform3D transform;
  // RDKit matches the probe to the reference itself (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  double rmsd = MolAlign::getAlignmentTransform(
      prb_mol, ref_mol, transform, prb_cid, ref_cid,
      map.empty() ? nullptr : &map, weights_vec.get(), reflect, max_iters);
  trans = to_values(transform);
  return rmsd;
}

// As alignMolecule in rdMolAlign.cpp.
double align_mol(RWMol &prb_mol, const RWMol &ref_mol, std::int32_t prb_cid,
                 std::int32_t ref_cid, rust::Slice<const IndexPair> atom_map,
                 rust::Slice<const double> weights, bool reflect,
                 std::uint32_t max_iters) {
  auto map = to_match(atom_map.data(), atom_map.size(), prb_mol, ref_mol);
  auto weights_vec = to_vector(weights);
  if (weights_vec) {
    unsigned n = map.empty() ? prb_mol.getNumAtoms() : map.size();
    if (weights_vec->size() != n) {
      throw ValueErrorException("Incorrect number of weights specified");
    }
  }
  // RDKit matches the probe to the reference itself (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return MolAlign::alignMol(prb_mol, ref_mol, prb_cid, ref_cid,
                            map.empty() ? nullptr : &map, weights_vec.get(),
                            reflect, max_iters);
}

double get_best_alignment_transform(const RWMol &prb_mol, const RWMol &ref_mol,
                                    const BestAlignmentParamsData &params,
                                    rust::Vec<double> &trans,
                                    rust::Vec<IndexPair> &best_match) {
  auto native = to_native(params, prb_mol, ref_mol, params.ignore_hs);
  RDGeom::Transform3D transform;
  MatchVectType match;
  double rmsd = MolAlign::getBestAlignmentTransform(
      prb_mol, ref_mol, transform, match, native.params, params.prb_cid,
      params.ref_cid, params.reflect, params.max_iters);
  trans = to_values(transform);
  best_match.clear();
  for (const auto &[p, r] : match) {
    best_match.push_back(IndexPair{static_cast<std::uint32_t>(p),
                                   static_cast<std::uint32_t>(r)});
  }
  return rmsd;
}

double get_best_rms(RWMol &prb_mol, const RWMol &ref_mol,
                    const BestAlignmentParamsData &params) {
  auto native = to_native(params, prb_mol, ref_mol, params.ignore_hs);
  return MolAlign::getBestRMS(prb_mol, ref_mol, native.params, params.prb_cid,
                              params.ref_cid);
}

rust::Vec<double> get_all_conformer_best_rms(
    const RWMol &mol, const BestAlignmentParamsData &params) {
  auto native = to_native(params, mol, mol, params.ignore_hs);
  rust::Vec<double> result;
  for (double v : MolAlign::getAllConformerBestRMS(mol, native.params)) {
    result.push_back(v);
  }
  return result;
}

rust::Vec<double> get_all_conformer_best_rms_to_ref(
    const RWMol &prb_mol, const RWMol &ref_mol,
    const BestAlignmentParamsData &params) {
  auto native = to_native(params, prb_mol, ref_mol, params.ignore_hs);
  rust::Vec<double> result;
  for (double v :
       MolAlign::getAllConformerBestRMSToRef(prb_mol, ref_mol, native.params)) {
    result.push_back(v);
  }
  return result;
}

// CalcRMS finds its matches with the hydrogens. It takes the probe as
// non-const, but only reads it.
double calc_rms(const RWMol &prb_mol, const RWMol &ref_mol,
                const BestAlignmentParamsData &params) {
  auto native = to_native(params, prb_mol, ref_mol, false);
  return MolAlign::CalcRMS(const_cast<RWMol &>(prb_mol), ref_mol,
                           params.prb_cid, params.ref_cid, native.params.map,
                           params.max_matches, false, native.weights.get());
}

// As alignMolConfs in rdMolAlign.cpp, returning what Python's RMSlist
// collects.
rust::Vec<double> align_mol_conformers(
    RWMol &mol, rust::Slice<const std::uint32_t> atom_ids,
    rust::Slice<const std::uint32_t> conf_ids,
    rust::Slice<const double> weights, bool reflect, std::uint32_t max_iters) {
  std::vector<unsigned int> atoms;
  for (auto idx : atom_ids) {
    check_index(idx, mol.getNumAtoms());
    atoms.push_back(idx);
  }
  std::vector<unsigned int> confs(conf_ids.begin(), conf_ids.end());
  auto weights_vec = to_vector(weights);
  std::vector<double> rms;
  MolAlign::alignMolConformers(mol, atoms.empty() ? nullptr : &atoms,
                               confs.empty() ? nullptr : &confs,
                               weights_vec.get(), reflect, max_iters, &rms);
  rust::Vec<double> result;
  for (double v : rms) {
    result.push_back(v);
  }
  return result;
}

// randomTransform seeds and draws from RDKit's global generator
// (random_mutex), and reads the bounds of a conformer without atoms
// uninitialized; that conformer has nothing to move.
void random_transform(RWMol &mol, std::int32_t cid, std::int32_t seed) {
  if (!mol.getConformer(cid).getNumAtoms()) {
    return;
  }
  std::lock_guard<std::mutex> random(random_mutex());
  MolAlign::randomTransform(mol, cid, seed);
}

std::unique_ptr<O3AList> get_mmff_o3a(RWMol &prb_mol, const RWMol &ref_mol,
                                      const MMFFMolProperties *prb_props,
                                      const MMFFMolProperties *ref_props,
                                      bool for_probe_confs,
                                      const O3AOptionsData &options) {
  auto constraints = to_constraints(options, prb_mol, ref_mol);
  check_o3a_conformers(prb_mol, ref_mol, for_probe_confs, options);
  std::unique_ptr<RDKit::MMFF::MMFFMolProperties> prb_owned, ref_owned;
  const auto *prb = mmff_props(prb_mol, prb_props, prb_owned, "probe");
  const auto *ref = mmff_props(ref_mol, ref_props, ref_owned, "reference");
  return run_o3a(prb_mol, ref_mol,
                 const_cast<RDKit::MMFF::MMFFMolProperties *>(prb),
                 const_cast<RDKit::MMFF::MMFFMolProperties *>(ref),
                 MolAlign::O3A::MMFF94, for_probe_confs, options, constraints);
}

std::unique_ptr<O3AList> get_crippen_o3a(RWMol &prb_mol, const RWMol &ref_mol,
                                         bool has_prb_contribs,
                                         rust::Slice<const double> prb_contribs,
                                         bool has_ref_contribs,
                                         rust::Slice<const double> ref_contribs,
                                         bool for_probe_confs,
                                         const O3AOptionsData &options) {
  auto constraints = to_constraints(options, prb_mol, ref_mol);
  check_o3a_conformers(prb_mol, ref_mol, for_probe_confs, options);
  auto prb = crippen_contribs(prb_mol, has_prb_contribs, prb_contribs);
  auto ref = crippen_contribs(ref_mol, has_ref_contribs, ref_contribs);
  return run_o3a(prb_mol, ref_mol, &prb, &ref, MolAlign::O3A::CRIPPEN,
                 for_probe_confs, options, constraints);
}

std::size_t o3a_list_size(const O3AList &list) { return list.items.size(); }

std::unique_ptr<O3AHolder> o3a_list_take(O3AList &list, std::size_t i) {
  return std::make_unique<O3AHolder>(std::move(list.items.at(i)));
}

double o3a_align(O3AHolder &o3a) { return o3a.get().align(); }

double o3a_trans(const O3AHolder &o3a, rust::Vec<double> &trans) {
  RDGeom::Transform3D transform;
  double rmsd = o3a.get().trans(transform);
  trans = to_values(transform);
  return rmsd;
}

double o3a_score(const O3AHolder &o3a) { return o3a.get().score(); }

rust::Vec<IndexPair> o3a_matches(const O3AHolder &o3a) {
  rust::Vec<IndexPair> result;
  for (const auto &[p, r] : *o3a.get().matches()) {
    result.push_back(IndexPair{static_cast<std::uint32_t>(p),
                               static_cast<std::uint32_t>(r)});
  }
  return result;
}

rust::Vec<double> o3a_weights(const O3AHolder &o3a) {
  rust::Vec<double> result;
  const auto *weights = o3a.get().weights();
  for (unsigned i = 0; i < weights->size(); ++i) {
    result.push_back((*weights)[i]);
  }
  return result;
}

}  // namespace rdkit_cxx
