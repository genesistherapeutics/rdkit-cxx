// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/fingerprint_generator.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <exception>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/fingerprint_generator.rs.h"
#include "rdkit-cxx/ffi/fingerprints.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/threads.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

using RDKit::AtomPair::AtomPairArguments;
using RDKit::MorganFingerprint::MorganArguments;
using RDKit::RDKitFP::RDKitFPArguments;
using RDKit::TopologicalTorsion::TopologicalTorsionArguments;

// The count bounds a factory passes RDKit, whose constructor asserts that
// count simulation has some; the check here comes first, as a ValueError.
std::vector<std::uint32_t> count_bounds(
    bool count_simulation, bool has_count_bounds,
    const rust::Vec<std::uint32_t> &bounds) {
  // As the generator wrappers in Code/GraphMol/Fingerprints/Wrap, except that
  // they treat an empty countBounds list as None; here it is used as given.
  if (!has_count_bounds) {
    return {1, 2, 4, 8};
  }
  if (count_simulation && bounds.empty()) {
    throw ValueErrorException("count simulation needs count bounds");
  }
  return std::vector<std::uint32_t>(bounds.begin(), bounds.end());
}

// The torsion code shifts each atom's code by codeSize (plus the chiral bits)
// per atom into 64 bits, and the result size by the whole length; a shift of
// 64 or more is undefined.
void check_torsion_atom_count(std::uint32_t count, bool include_chirality) {
  std::uint64_t bits =
      RDKit::AtomPairs::codeSize +
      (include_chirality ? RDKit::AtomPairs::numChiralBits : 0);
  if (std::uint64_t{count} * bits > 63) {
    throw ValueErrorException(
        "torsionAtomCount is too large: its codes must fit in 63 bits (7 "
        "atoms, or 5 with chirality)");
  }
}

// What every fingerprint call trusts about the generator's options. They
// are checked when a generator is made or given new options.
void check_option_values(std::uint32_t fp_size, bool count_simulation,
                         bool no_count_bounds,
                         std::uint32_t num_bits_per_feature) {
  // U3: the bit vectors hold fpSize bits, and their int arithmetic needs no
  // more than INT_MAX.
  if (fp_size > static_cast<std::uint32_t>(INT_MAX)) {
    throw ValueErrorException("fpSize must be at most 2^31 - 1");
  }
  // getSparseFingerprint divides by the number of count bounds.
  if (count_simulation && no_count_bounds) {
    throw ValueErrorException("count simulation needs count bounds");
  }
  // As FingerprintArguments' constructor asserts.
  if (num_bits_per_feature == 0) {
    throw ValueErrorException("numBitsPerFeature must be at least 1");
  }
}

void check_options(const RDKit::FingerprintArguments &args) {
  check_option_values(args.d_fpSize, args.df_countSimulation,
                      args.d_countBounds.empty(), args.d_numBitsPerFeature);
  if (auto *tt = dynamic_cast<const TopologicalTorsionArguments *>(&args)) {
    check_torsion_atom_count(tt->d_torsionAtomCount, tt->df_includeChirality);
  }
}

std::unique_ptr<FpGenerator> checked(FpGenerator *gen) {
  std::unique_ptr<FpGenerator> owned(gen);
  check_options(*owned->getOptions());
  return owned;
}

// Python's empty lists mean "none" (convertPyArguments in
// FingerprintGeneratorWrapper.cpp).
std::unique_ptr<std::vector<std::uint32_t>> optional_list(
    bool has, const rust::Vec<std::uint32_t> &values) {
  if (!has || values.empty()) {
    return nullptr;
  }
  return std::make_unique<std::vector<std::uint32_t>>(values.begin(),
                                                      values.end());
}

bool uses_3d_distances(const FpGenerator &gen) {
  auto *ap = dynamic_cast<const AtomPairArguments *>(gen.getOptions());
  return ap && !ap->df_use2D;
}

// The atom-pair generator casts each distance to unsigned int; for 3D
// distances that is undefined unless the distance is finite and below 2^32.
// Coordinates within 1e9 keep every distance below 3.5e9.
void check_3d_distances(const RDKit::ROMol &mol, int conf_id) {
  const auto &conf = mol.getConformer(conf_id);
  for (const auto &pos : conf.getPositions()) {
    for (double v : {pos.x, pos.y, pos.z}) {
      if (!(std::abs(v) <= 1e9)) {
        throw ValueErrorException(
            "atom-pair fingerprints from 3D distances need finite coordinates "
            "within 1e9");
      }
    }
  }
}

// What one call trusts about the molecule and the arguments. RDKit sets bits
// in boost::dynamic_bitset<> for fromAtoms and ignoreAtoms without checking
// them (Morgan, RDKit, torsion), and indexes the custom invariants by atom
// and bond (the atom-pair check is off by one).
void check_call(const FpGenerator &gen, const RDKit::ROMol &mol,
                const std::vector<std::uint32_t> *from_atoms,
                const std::vector<std::uint32_t> *ignore_atoms,
                const std::vector<std::uint32_t> *atom_invariants,
                const std::vector<std::uint32_t> *bond_invariants,
                int conf_id) {
  for (const auto *atoms : {from_atoms, ignore_atoms}) {
    if (atoms) {
      for (auto idx : *atoms) {
        check_index(idx, mol.getNumAtoms());
      }
    }
  }
  if (atom_invariants && atom_invariants->size() < mol.getNumAtoms()) {
    throw ValueErrorException("customAtomInvariants needs one per atom");
  }
  if (bond_invariants && bond_invariants->size() < mol.getNumBonds()) {
    throw ValueErrorException("customBondInvariants needs one per bond");
  }
  // Atom pairs read the n*n distance matrix.
  if (dynamic_cast<const AtomPairArguments *>(gen.getOptions())) {
    check_matrix_size(mol);
  }
  if (uses_3d_distances(gen)) {
    check_3d_distances(mol, conf_id);
  }
}

// Python's AdditionalOutput, rebuilt for each call.
class Outputs {
 public:
  explicit Outputs(const AdditionalOutputData &data) {
    if (data.collect_atom_to_bits) {
      d_ao.allocateAtomToBits();
    }
    if (data.collect_bit_info_map) {
      d_ao.allocateBitInfoMap();
    }
    if (data.collect_bit_paths) {
      d_ao.allocateBitPaths();
    }
    if (data.collect_atom_counts) {
      d_ao.allocateAtomCounts();
    }
    d_any = data.collect_atom_to_bits || data.collect_bit_info_map ||
            data.collect_bit_paths || data.collect_atom_counts;
  }

  RDKit::AdditionalOutput *get() { return d_any ? &d_ao : nullptr; }

  void copy_to(AdditionalOutputData &data) const {
    if (d_ao.atomToBits) {
      for (const auto &bits : *d_ao.atomToBits) {
        data.atom_to_bits.push_back(U64List{to_vec(bits)});
      }
    }
    if (d_ao.bitInfoMap) {
      for (const auto &[bit, infos] : *d_ao.bitInfoMap) {
        data.bit_info_map_keys.push_back(bit);
        data.bit_info_map_counts.push_back(
            static_cast<std::uint32_t>(infos.size()));
        for (const auto &[first, second] : infos) {
          data.bit_info_map_values.push_back(IndexPair{first, second});
        }
      }
    }
    copy_lists(d_ao.bitPaths, data.bit_paths_keys, data.bit_paths_counts,
               data.bit_paths_values);
    if (d_ao.atomCounts) {
      data.atom_counts = to_vec(*d_ao.atomCounts);
    }
  }

 private:
  template <class Map>
  static void copy_lists(const Map *map, rust::Vec<std::uint64_t> &keys,
                         rust::Vec<std::uint32_t> &counts,
                         rust::Vec<IndexList> &values) {
    if (!map) {
      return;
    }
    for (const auto &[bit, lists] : *map) {
      keys.push_back(bit);
      counts.push_back(static_cast<std::uint32_t>(lists.size()));
      for (const auto &list : lists) {
        rust::Vec<std::uint32_t> indices;
        for (auto idx : list) {
          indices.push_back(static_cast<std::uint32_t>(idx));
        }
        values.push_back(IndexList{std::move(indices)});
      }
    }
  }

  RDKit::AdditionalOutput d_ao;
  bool d_any = false;
};

// Runs `call` with the arguments RDKit takes, after the checks, and with the
// molecule's matrix caches cleared around it (U6).
template <class Call>
auto with_args(const FpGenerator &gen, const RWMol &mol,
               const FingerprintArgsData &data, AdditionalOutputData &outputs,
               Call call) {
  auto from_atoms = optional_list(data.has_from_atoms, data.from_atoms);
  auto ignore_atoms = optional_list(data.has_ignore_atoms, data.ignore_atoms);
  auto atom_invariants = optional_list(data.has_custom_atom_invariants,
                                       data.custom_atom_invariants);
  auto bond_invariants = optional_list(data.has_custom_bond_invariants,
                                       data.custom_bond_invariants);
  check_call(gen, mol, from_atoms.get(), ignore_atoms.get(),
             atom_invariants.get(), bond_invariants.get(), data.conf_id);
  Outputs ao(outputs);
  RDKit::FingerprintFuncArguments args(
      from_atoms.get(), ignore_atoms.get(), data.conf_id, ao.get(),
      atom_invariants.get(), bond_invariants.get());
  MatrixCacheGuard guard(mol);
  auto control_c = lock_control_c_for(gen);
  auto result = call(args);
  ao.copy_to(outputs);
  return result;
}

std::unique_ptr<BitVect> bits(const FpGenerator &gen, const RDKit::ROMol &mol,
                              FingerprintOutput output,
                              RDKit::FingerprintFuncArguments &args) {
  switch (output) {
    case FingerprintOutput::Bits:
      return wrap_bit_vect(gen.getFingerprint(mol, args));
    case FingerprintOutput::SparseBits:
      return wrap_bit_vect(gen.getSparseFingerprint(mol, args));
    default:
      throw ValueErrorException("not a bit-vector fingerprint");
  }
}

std::unique_ptr<SparseIntVect> counts(const FpGenerator &gen,
                                      const RDKit::ROMol &mol,
                                      FingerprintOutput output,
                                      RDKit::FingerprintFuncArguments &args) {
  switch (output) {
    case FingerprintOutput::Counts:
      return wrap_sparse_int_vect(gen.getCountFingerprint(mol, args));
    case FingerprintOutput::SparseCounts:
      return wrap_sparse_int_vect(gen.getSparseCountFingerprint(mol, args));
    default:
      throw ValueErrorException("not a count fingerprint");
  }
}

bool is_bits(FingerprintOutput output) {
  return output == FingerprintOutput::Bits ||
         output == FingerprintOutput::SparseBits;
}

// As mtgetFingerprints in FingerprintGenerator.cpp, which runs RDKit's
// threads without catching: an exception in one (an invalid conformer, a
// count-bounds check) would call std::terminate. Here each thread keeps its
// first exception, and the call rethrows it. Each molecule is read by one
// thread only, as RDKit's own split does. RDKit's threads would each
// construct a ControlCHandler when they match feature patterns or run the
// CIP labeler, so a call that may (lock_control_c_for) runs on this thread
// alone.
std::unique_ptr<FingerprintList> run_bulk(const FpGenerator &gen,
                                          BorrowedMols &mols,
                                          FingerprintOutput output,
                                          int num_threads) {
  const auto &list = mols.mols();
  for (const auto &mol : list) {
    check_call(gen, *mol, nullptr, nullptr, nullptr, nullptr, -1);
  }
  auto result = std::make_unique<FingerprintList>();
  if (is_bits(output)) {
    result->bit_vects.resize(list.size());
  } else {
    result->int_vects.resize(list.size());
  }
  auto work = [&](std::size_t i) {
    RDKit::FingerprintFuncArguments args;
    MatrixCacheGuard guard(*list[i]);
    if (is_bits(output)) {
      result->bit_vects[i] = bits(gen, *list[i], output, args);
    } else {
      result->int_vects[i] = counts(gen, *list[i], output, args);
    }
  };
  auto control_c = lock_control_c_for(gen);
  // RDKit's getNumThreadsToUse negates a negative count in int, undefined
  // for INT_MIN (num_threads_to_use). Threads past the molecules would have
  // nothing to do, and run_on_threads joins the threads it started when
  // another cannot start, where destroying a joinable std::thread would
  // call std::terminate.
  std::size_t n = control_c.owns_lock() ? 1 : num_threads_to_use(num_threads);
  n = std::min(n, list.size());
  if (n <= 1) {
    for (std::size_t i = 0; i < list.size(); ++i) {
      work(i);
    }
    return result;
  }
  run_on_threads(static_cast<unsigned>(n), [&](unsigned t) {
    for (std::size_t i = t; i < list.size(); i += n) {
      work(i);
    }
  });
  return result;
}

}  // namespace

// Holds control_c_mutex while the generator may construct a ControlCHandler:
// always for feature invariants, which MorganFeatureAtomInvGenerator finds
// by substructure matches (a generator that has one names it in its info
// string), and otherwise when it may run the CIP labeler
// (lock_control_c_for_cip).
std::unique_lock<ControlCMutex> lock_control_c_for(const FpGenerator &gen) {
  if (gen.infoString().find("MorganFeatureInvariantGenerator") !=
      std::string::npos) {
    return std::unique_lock<ControlCMutex>(control_c_mutex());
  }
  return lock_control_c_for_cip();
}

bool fingerprint_reads_distance_matrix(const FpGenerator &gen) {
  return dynamic_cast<const AtomPairArguments *>(gen.getOptions()) != nullptr;
}

std::unique_ptr<FpGenerator> get_morgan_generator(
    const MorganGeneratorData &data, std::unique_ptr<AtomInvGen> atom_inv_gen,
    std::unique_ptr<BondInvGen> bond_inv_gen) {
  // As getMorganGenerator in Code/GraphMol/Fingerprints/Wrap/
  // MorganWrapper.cpp, which ignores includeRingMembership.
  return checked(RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(
      data.radius, data.count_simulation, data.include_chirality,
      data.use_bond_types, data.only_nonzero_invariants,
      data.include_redundant_environments, atom_inv_gen.release(),
      bond_inv_gen.release(), data.fp_size,
      count_bounds(data.count_simulation, data.has_count_bounds,
                   data.count_bounds),
      true, true));
}

std::unique_ptr<FpGenerator> get_atom_pair_generator(
    const AtomPairGeneratorData &data,
    std::unique_ptr<AtomInvGen> atom_inv_gen) {
  return checked(RDKit::AtomPair::getAtomPairGenerator<std::uint64_t>(
      data.min_distance, data.max_distance, data.include_chirality, data.use_2d,
      atom_inv_gen.release(), data.count_simulation, data.fp_size,
      count_bounds(data.count_simulation, data.has_count_bounds,
                   data.count_bounds),
      true));
}

std::unique_ptr<FpGenerator> get_rdkit_fp_generator(
    const RDKitFPGeneratorData &data,
    std::unique_ptr<AtomInvGen> atom_inv_gen) {
  if (data.num_bits_per_feature == 0) {
    throw ValueErrorException("numBitsPerFeature must be at least 1");
  }
  return checked(RDKit::RDKitFP::getRDKitFPGenerator<std::uint64_t>(
      data.min_path, data.max_path, data.use_hs, data.branched_paths,
      data.use_bond_order, atom_inv_gen.release(), data.count_simulation,
      count_bounds(data.count_simulation, data.has_count_bounds,
                   data.count_bounds),
      data.fp_size, data.num_bits_per_feature, true));
}

// Python's wrapper passes ownsAtomInvGen = false for the clone it makes,
// which leaks; this one owns it.
std::unique_ptr<FpGenerator> get_topological_torsion_generator(
    const TopologicalTorsionGeneratorData &data,
    std::unique_ptr<AtomInvGen> atom_inv_gen) {
  // Check the length before RDKit computes anything from it.
  check_torsion_atom_count(data.torsion_atom_count, data.include_chirality);
  return checked(
      RDKit::TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
          data.include_chirality, data.torsion_atom_count,
          atom_inv_gen.release(), data.count_simulation, data.fp_size,
          count_bounds(data.count_simulation, data.has_count_bounds,
                       data.count_bounds),
          true));
}

std::unique_ptr<AtomInvGen> get_morgan_atom_inv_gen(
    bool include_ring_membership) {
  return std::make_unique<RDKit::MorganFingerprint::MorganAtomInvGenerator>(
      include_ring_membership);
}

// The generator copies the patterns.
std::unique_ptr<AtomInvGen> get_morgan_feature_atom_inv_gen(
    bool has_patterns, BorrowedMols &patterns) {
  if (!has_patterns) {
    return std::make_unique<
        RDKit::MorganFingerprint::MorganFeatureAtomInvGenerator>(nullptr);
  }
  std::vector<const RDKit::ROMol *> mols;
  for (const auto &mol : patterns.mols()) {
    mols.push_back(mol.get());
  }
  return std::make_unique<
      RDKit::MorganFingerprint::MorganFeatureAtomInvGenerator>(&mols);
}

std::unique_ptr<BondInvGen> get_morgan_bond_inv_gen(bool use_bond_types,
                                                    bool use_chirality) {
  return std::make_unique<RDKit::MorganFingerprint::MorganBondInvGenerator>(
      use_bond_types, use_chirality);
}

std::unique_ptr<AtomInvGen> get_atom_pair_atom_inv_gen(bool include_chirality) {
  return std::make_unique<RDKit::AtomPair::AtomPairAtomInvGenerator>(
      include_chirality);
}

std::unique_ptr<AtomInvGen> get_rdkit_atom_inv_gen() {
  return std::make_unique<RDKit::RDKitFP::RDKitFPAtomInvGenerator>();
}

std::unique_ptr<AtomInvGen> atom_inv_gen_clone(const AtomInvGen &gen) {
  return std::unique_ptr<AtomInvGen>(gen.clone());
}

std::unique_ptr<BondInvGen> bond_inv_gen_clone(const BondInvGen &gen) {
  return std::unique_ptr<BondInvGen>(gen.clone());
}

std::unique_ptr<std::string> atom_inv_gen_info_string(const AtomInvGen &gen) {
  return std::make_unique<std::string>(gen.infoString());
}

std::unique_ptr<std::string> bond_inv_gen_info_string(const BondInvGen &gen) {
  return std::make_unique<std::string>(gen.infoString());
}

FingerprintOptionsData fp_generator_get_options(const FpGenerator &gen) {
  const auto *args = gen.getOptions();
  FingerprintOptionsData data{};
  data.count_simulation = args->df_countSimulation;
  data.include_chirality = args->df_includeChirality;
  data.fp_size = args->d_fpSize;
  data.num_bits_per_feature = args->d_numBitsPerFeature;
  data.count_bounds = to_vec(args->d_countBounds);
  if (auto *m = dynamic_cast<const MorganArguments *>(args)) {
    data.kind = FPType::MorganFP;
    data.radius = m->d_radius;
    data.only_nonzero_invariants = m->df_onlyNonzeroInvariants;
    data.include_redundant_environments = m->df_includeRedundantEnvironments;
  } else if (auto *ap = dynamic_cast<const AtomPairArguments *>(args)) {
    data.kind = FPType::AtomPairFP;
    data.min_distance = ap->d_minDistance;
    data.max_distance = ap->d_maxDistance;
    data.use_2d = ap->df_use2D;
  } else if (auto *r = dynamic_cast<const RDKitFPArguments *>(args)) {
    data.kind = FPType::RDKitFP;
    data.min_path = r->d_minPath;
    data.max_path = r->d_maxPath;
    data.use_hs = r->df_useHs;
    data.branched_paths = r->df_branchedPaths;
    data.use_bond_order = r->df_useBondOrder;
  } else if (auto *tt =
                 dynamic_cast<const TopologicalTorsionArguments *>(args)) {
    data.kind = FPType::TopologicalTorsionFP;
    data.torsion_atom_count = tt->d_torsionAtomCount;
    data.only_shortest_paths = tt->df_onlyShortestPaths;
  } else {
    throw ValueErrorException("unknown fingerprint options type");
  }
  return data;
}

// Python's GetOptions returns the generator's own options to change in
// place; here the new options are checked as a whole first. RDKit 2024.03's
// FingerprintArguments cannot be copied, so the values are checked before
// any is set.
void fp_generator_set_options(FpGenerator &gen,
                              const FingerprintOptionsData &data) {
  auto *args = gen.getOptions();
  auto *morgan = dynamic_cast<MorganArguments *>(args);
  auto *atom_pair = dynamic_cast<AtomPairArguments *>(args);
  auto *rdkit_fp = dynamic_cast<RDKitFPArguments *>(args);
  auto *torsion = dynamic_cast<TopologicalTorsionArguments *>(args);
  bool matches = (data.kind == FPType::MorganFP && morgan) ||
                 (data.kind == FPType::AtomPairFP && atom_pair) ||
                 (data.kind == FPType::RDKitFP && rdkit_fp) ||
                 (data.kind == FPType::TopologicalTorsionFP && torsion);
  if (!matches) {
    throw ValueErrorException("the options are for another fingerprint");
  }
  check_option_values(data.fp_size, data.count_simulation,
                      data.count_bounds.empty(), data.num_bits_per_feature);
  if (torsion) {
    check_torsion_atom_count(data.torsion_atom_count, data.include_chirality);
  }
  args->df_countSimulation = data.count_simulation;
  args->df_includeChirality = data.include_chirality;
  args->d_fpSize = data.fp_size;
  args->d_numBitsPerFeature = data.num_bits_per_feature;
  args->d_countBounds.assign(data.count_bounds.begin(),
                             data.count_bounds.end());
  if (morgan) {
    morgan->d_radius = data.radius;
    morgan->df_onlyNonzeroInvariants = data.only_nonzero_invariants;
    morgan->df_includeRedundantEnvironments =
        data.include_redundant_environments;
  } else if (atom_pair) {
    atom_pair->d_minDistance = data.min_distance;
    atom_pair->d_maxDistance = data.max_distance;
    atom_pair->df_use2D = data.use_2d;
  } else if (rdkit_fp) {
    rdkit_fp->d_minPath = data.min_path;
    rdkit_fp->d_maxPath = data.max_path;
    rdkit_fp->df_useHs = data.use_hs;
    rdkit_fp->df_branchedPaths = data.branched_paths;
    rdkit_fp->df_useBondOrder = data.use_bond_order;
  } else {
    torsion->d_torsionAtomCount = data.torsion_atom_count;
    torsion->df_onlyShortestPaths = data.only_shortest_paths;
  }
}

std::unique_ptr<std::string> fp_generator_info_string(const FpGenerator &gen) {
  return std::make_unique<std::string>(gen.infoString());
}

std::unique_ptr<BitVect> fp_generator_bits(const FpGenerator &gen,
                                           const RWMol &mol,
                                           FingerprintOutput output,
                                           const FingerprintArgsData &args,
                                           AdditionalOutputData &outputs) {
  return with_args(gen, mol, args, outputs,
                   [&](RDKit::FingerprintFuncArguments &a) {
                     return bits(gen, mol, output, a);
                   });
}

std::unique_ptr<SparseIntVect> fp_generator_counts(
    const FpGenerator &gen, const RWMol &mol, FingerprintOutput output,
    const FingerprintArgsData &args, AdditionalOutputData &outputs) {
  return with_args(gen, mol, args, outputs,
                   [&](RDKit::FingerprintFuncArguments &a) {
                     return counts(gen, mol, output, a);
                   });
}

// As getNumPyFingerprint and getNumPyCountFingerprint in
// FingerprintGeneratorWrapper.cpp: one entry per bit.
rust::Vec<std::uint32_t> fp_generator_values(const FpGenerator &gen,
                                             const RWMol &mol,
                                             FingerprintOutput output,
                                             const FingerprintArgsData &args,
                                             AdditionalOutputData &outputs) {
  return with_args(
      gen, mol, args, outputs, [&](RDKit::FingerprintFuncArguments &a) {
        rust::Vec<std::uint32_t> values;
        if (output == FingerprintOutput::Bits) {
          auto fp = gen.getFingerprint(mol, a);
          values.reserve(fp->size());
          for (unsigned i = 0; i < fp->size(); ++i) {
            values.push_back((*fp)[i] ? 1 : 0);
          }
        } else if (output == FingerprintOutput::Counts) {
          auto fp = gen.getCountFingerprint(mol, a);
          auto size = static_cast<unsigned>(fp->getLength());
          values.reserve(size);
          for (unsigned i = 0; i < size; ++i) {
            values.push_back(static_cast<std::uint32_t>(fp->getVal(i)));
          }
        } else {
          throw ValueErrorException("NumPy fingerprints are bits or counts");
        }
        return values;
      });
}

std::unique_ptr<FingerprintList> fp_generator_bulk(const FpGenerator &gen,
                                                   BorrowedMols &mols,
                                                   FingerprintOutput output,
                                                   std::int32_t num_threads) {
  return run_bulk(gen, mols, output, num_threads);
}

std::unique_ptr<FingerprintList> fp_type_bulk(BorrowedMols &mols,
                                              FPType fp_type,
                                              FingerprintOutput output) {
  std::unique_ptr<FpGenerator> gen;
  switch (fp_type) {
    case FPType::AtomPairFP:
      gen.reset(RDKit::AtomPair::getAtomPairGenerator<std::uint64_t>());
      break;
    case FPType::MorganFP:
      gen.reset(RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
      break;
    case FPType::RDKitFP:
      gen.reset(RDKit::RDKitFP::getRDKitFPGenerator<std::uint64_t>());
      break;
    case FPType::TopologicalTorsionFP:
      gen.reset(RDKit::TopologicalTorsion::getTopologicalTorsionGenerator<
                std::uint64_t>());
      break;
    default:
      throw ValueErrorException("unknown fingerprint type");
  }
  return run_bulk(*gen, mols, output, 1);
}

void check_fingerprint_mol(const FpGenerator &gen, const RDKit::ROMol &mol) {
  check_call(gen, mol, nullptr, nullptr, nullptr, nullptr, -1, false);
}

std::size_t fingerprint_list_size(const FingerprintList &list) {
  return list.bit_vects.empty() ? list.int_vects.size() : list.bit_vects.size();
}

std::unique_ptr<BitVect> fingerprint_list_take_bits(FingerprintList &list,
                                                    std::size_t i) {
  return std::move(list.bit_vects.at(i));
}

std::unique_ptr<SparseIntVect> fingerprint_list_take_counts(
    FingerprintList &list, std::size_t i) {
  return std::move(list.int_vects.at(i));
}

}  // namespace rdkit_cxx
