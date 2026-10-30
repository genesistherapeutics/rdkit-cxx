// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/embedder.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <mutex>
#include <tuple>
#include <utility>
#include <vector>

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"
#include "rdkit-cxx/ffi/embedder.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/force_field.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/symmetrize_lock.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

namespace DG = RDKit::DGeomHelpers;

// Each n * n matrix is indexed in unsigned int (U3).
constexpr std::uint32_t MAX_BOUNDS_MAT_SIZE = 46340;

// initBoundsMat's default upper bound, which no bound from a molecule's
// topology exceeds; triangle smoothing only lowers upper bounds.
constexpr double MAX_TOPOLOGICAL_BOUND = 1000.0;

// Copies the fields Python exposes as attributes.
void set_fields(EmbedParametersData &d, const DG::EmbedParameters &p) {
  d.max_iterations = p.maxIterations;
  d.num_threads = p.numThreads;
  d.timeout = p.timeout;
  d.random_seed = p.randomSeed;
  d.clear_confs = p.clearConfs;
  d.use_random_coords = p.useRandomCoords;
  d.box_size_mult = p.boxSizeMult;
  d.rand_neg_eig = p.randNegEig;
  d.num_zero_fail = p.numZeroFail;
  d.optimizer_force_tol = p.optimizerForceTol;
  d.basin_thresh = p.basinThresh;
  d.ignore_smoothing_failures = p.ignoreSmoothingFailures;
  d.enforce_chirality = p.enforceChirality;
  d.use_exp_torsion_angle_prefs = p.useExpTorsionAnglePrefs;
  d.use_basic_knowledge = p.useBasicKnowledge;
  d.et_version = p.ETversion;
  d.verbose = p.verbose;
  d.prune_rms_thresh = p.pruneRmsThresh;
  d.only_heavy_atoms_for_rms = p.onlyHeavyAtomsForRMS;
  d.embed_fragments_separately = p.embedFragmentsSeparately;
  d.use_small_ring_torsions = p.useSmallRingTorsions;
  d.use_macrocycle_torsions = p.useMacrocycleTorsions;
  d.use_macrocycle14config = p.useMacrocycle14config;
  d.bounds_mat_force_scaling = p.boundsMatForceScaling;
  d.use_symmetry_for_pruning = p.useSymmetryForPruning;
  d.force_trans_amides = p.forceTransAmides;
  d.track_failures = p.trackFailures;
  d.enable_sequential_random_seeds = p.enableSequentialRandomSeeds;
  d.symmetrize_conjugated_terminal_groups_for_pruning =
      p.symmetrizeConjugatedTerminalGroupsForPruning;
  d.failures = to_vec(p.failures);
}

// EmbedParameters points at its coordinate map, so the two live together.
struct NativeParams {
  DG::EmbedParameters params;
  std::map<int, RDGeom::Point3D> coord_map;
};

// As PyEmbedParameters' setters in Code/GraphMol/DistGeomHelpers/Wrap/
// rdDistGeom.cpp, except that an empty map is no map here, where Python sets
// an empty one.
std::unique_ptr<NativeParams> to_native(const EmbedParametersData &d) {
  auto native = std::make_unique<NativeParams>();
  auto &p = native->params;
  p.maxIterations = d.max_iterations;
  p.numThreads = d.num_threads;
  p.timeout = d.timeout;
  p.randomSeed = d.random_seed;
  p.clearConfs = d.clear_confs;
  p.useRandomCoords = d.use_random_coords;
  p.boxSizeMult = d.box_size_mult;
  p.randNegEig = d.rand_neg_eig;
  p.numZeroFail = d.num_zero_fail;
  p.optimizerForceTol = d.optimizer_force_tol;
  p.basinThresh = d.basin_thresh;
  p.ignoreSmoothingFailures = d.ignore_smoothing_failures;
  p.enforceChirality = d.enforce_chirality;
  p.useExpTorsionAnglePrefs = d.use_exp_torsion_angle_prefs;
  p.useBasicKnowledge = d.use_basic_knowledge;
  p.ETversion = d.et_version;
  p.verbose = d.verbose;
  p.pruneRmsThresh = d.prune_rms_thresh;
  p.onlyHeavyAtomsForRMS = d.only_heavy_atoms_for_rms;
  p.embedFragmentsSeparately = d.embed_fragments_separately;
  p.useSmallRingTorsions = d.use_small_ring_torsions;
  p.useMacrocycleTorsions = d.use_macrocycle_torsions;
  p.useMacrocycle14config = d.use_macrocycle14config;
  p.boundsMatForceScaling = d.bounds_mat_force_scaling;
  p.useSymmetryForPruning = d.use_symmetry_for_pruning;
  p.forceTransAmides = d.force_trans_amides;
  p.trackFailures = d.track_failures;
  p.enableSequentialRandomSeeds = d.enable_sequential_random_seeds;
  p.symmetrizeConjugatedTerminalGroupsForPruning =
      d.symmetrize_conjugated_terminal_groups_for_pruning;
  p.failures.assign(d.failures.begin(), d.failures.end());

  for (std::size_t i = 0; i < d.coord_map_atoms.size(); ++i) {
    const auto &point = d.coord_map_points.at(i);
    native->coord_map[static_cast<int>(d.coord_map_atoms[i])] =
        RDGeom::Point3D(point.x, point.y, point.z);
  }
  if (!native->coord_map.empty()) {
    p.coordMap = &native->coord_map;
  }

  if (!d.cpci_pairs.empty()) {
    p.CPCI = std::make_shared<
        std::map<std::pair<unsigned int, unsigned int>, double>>();
    for (std::size_t i = 0; i < d.cpci_pairs.size(); ++i) {
      const auto &pair = d.cpci_pairs[i];
      (*p.CPCI)[{pair.first, pair.second}] = d.cpci_values.at(i);
    }
  }

  // As setBoundsMatrix, whose nrows * nrows wraps in unsigned int. An
  // infinite bound becomes an infinite distance, which RDKit casts to int
  // (check_initial_coords_cast).
  if (d.has_bounds_mat) {
    std::uint32_t n = d.bounds_mat_size;
    if (n == 0) {
      throw ValueErrorException("The array has to have a nonzero size");
    }
    if (n > MAX_BOUNDS_MAT_SIZE ||
        d.bounds_mat.size() != static_cast<std::size_t>(n) * n) {
      throw ValueErrorException("The array has to be square");
    }
    DistGeom::BoundsMatrix::DATA_SPTR data(new double[d.bounds_mat.size()]);
    for (std::size_t i = 0; i < d.bounds_mat.size(); ++i) {
      if (!std::isfinite(d.bounds_mat[i])) {
        throw ValueErrorException("the bounds matrix has non-finite values");
      }
      data[i] = d.bounds_mat[i];
    }
    p.boundsMat = boost::shared_ptr<const DistGeom::BoundsMatrix>(
        new DistGeom::BoundsMatrix(n, data));
  }
  return native;
}

// The embedder computes each piece's int distance matrix and minimizes it in
// four dimensions (check_matrix_size, check_minimize_size). Returns the
// number of pieces and the size of the largest. It embeds fragments as
// copies of the molecule with the other atoms removed (getMolFrags), and
// prunes by RMS on a copy without its hydrogens, both of which trust what
// check_removals checks.
std::pair<unsigned, unsigned> check_embed_pieces(
    const RDKit::ROMol &mol, const DG::EmbedParameters &params) {
  unsigned num_atoms = mol.getNumAtoms();
  if (num_atoms == 0) {
    throw ValueErrorException("molecule has no atoms");
  }
  check_matrix_size(mol);
  check_removals(mol, Removes::FragmentsAndHydrogens);
  unsigned num_frags = 1;
  unsigned largest = num_atoms;
  if (params.embedFragmentsSeparately) {
    std::vector<int> mapping;
    num_frags = RDKit::MolOps::getMolFrags(mol, mapping);
    std::vector<unsigned> sizes(num_frags);
    for (auto frag : mapping) {
      ++sizes[frag];
    }
    largest = *std::max_element(sizes.begin(), sizes.end());
  }
  check_minimize_size(largest, 4);
  return {num_frags, largest};
}

// With distance-geometry starting coordinates (without useRandomCoords) the
// embedder casts a sum of squared distances to int
// (check_initial_coords_cast).
void check_initial_coords(const DG::EmbedParameters &params, unsigned largest,
                          double max_bound) {
  if (!params.useRandomCoords) {
    check_initial_coords_cast(largest, max_bound);
  }
}

// The embedder trusts several of its inputs:
// - the sizes check_embed_pieces checks;
// - CPCI's indices are checked by an invariant inside the minimizers, and
//   the coordinate map's by RDKit with ValueError (E3 puts the wrapper's
//   check first); with more than one fragment RDKit drops the map, except
//   that random coordinates still index each fragment's points with it,
//   through vector::at;
// - the bounds behind the initial coordinates (check_initial_coords).
void check_embed(const RWMol &mol, const EmbedParametersData &d,
                 const DG::EmbedParameters &params) {
  unsigned num_atoms = mol.getNumAtoms();
  auto [num_frags, largest] = check_embed_pieces(mol, params);

  double max_bound = MAX_TOPOLOGICAL_BOUND;
  if (!d.coord_map_atoms.empty()) {
    if (num_frags == 1) {
      // RDKit's own checks, with the index check first (E3).
      for (std::size_t i = 0; i < d.coord_map_atoms.size(); ++i) {
        check_index(d.coord_map_atoms[i], num_atoms);
        const auto &p = d.coord_map_points[i];
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) {
          throw ValueErrorException("coordMap contains non-finite coordinates");
        }
      }
      // The map's distances become bounds.
      for (const auto &[i, pi] : *params.coordMap) {
        for (const auto &[j, pj] : *params.coordMap) {
          if (i < j) {
            max_bound = std::max(max_bound, (pi - pj).length());
          }
        }
      }
    }
  }
  for (const auto &pair : d.cpci_pairs) {
    check_index(pair.first, num_atoms);
    check_index(pair.second, num_atoms);
  }
  if (params.boundsMat) {
    const auto *data = params.boundsMat->getData();
    for (std::size_t i = 0; i < params.boundsMat->getDataSize(); ++i) {
      max_bound = std::max(max_bound, std::abs(data[i]));
    }
  }
  check_initial_coords(params, largest, max_bound);
}

}  // namespace

void check_embed_mol(const RDKit::ROMol &mol,
                     const DG::EmbedParameters &params) {
  auto largest = check_embed_pieces(mol, params).second;
  check_initial_coords(params, largest, MAX_TOPOLOGICAL_BOUND);
}

void check_initial_coords_cast(unsigned n, double max_bound) {
  double worst = (static_cast<double>(n) - 1.0) * max_bound * max_bound / 2.0;
  if (!(worst < static_cast<double>(INT_MAX))) {
    throw ValueErrorException(
        "RDKit casts the sum of the squared embedding distances to int, which "
        "overflows for this many atoms and bounds this long");
  }
}

EmbedParametersData embed_parameters_preset(EmbedPreset preset) {
  const DG::EmbedParameters *params = nullptr;
  DG::EmbedParameters defaults;
  switch (preset) {
    case EmbedPreset::Default:
      params = &defaults;
      break;
    case EmbedPreset::ETKDG:
      params = &DG::ETKDG;
      break;
    case EmbedPreset::ETKDGv2:
      params = &DG::ETKDGv2;
      break;
    case EmbedPreset::ETKDGv3:
      params = &DG::ETKDGv3;
      break;
    case EmbedPreset::SrETKDGv3:
      params = &DG::srETKDGv3;
      break;
    case EmbedPreset::ETDG:
      params = &DG::ETDG;
      break;
    case EmbedPreset::KDG:
      params = &DG::KDG;
      break;
  }
  EmbedParametersData data{};
  set_fields(data, params ? *params : defaults);
  return data;
}

// RDKit's own threads race on the success flag of each conformer, which
// share the words of a boost::dynamic_bitset; on the parameters, whose
// maxIterations and basinThresh each thread sets; and, without a seed, on
// RDKit's global random generator. The embedding therefore runs on the
// calling thread, which gives the conformers RDKit's threads would, since
// each conformer's seed depends only on its index. Each call also swaps the
// process's SIGINT handler (control_c_mutex), and pruning symmetrizes
// terminal groups (symmetrize_mutex).
rust::Vec<std::int32_t> embed_multiple_confs(RWMol &mol,
                                             std::uint32_t num_confs,
                                             EmbedParametersData &data) {
  auto native = to_native(data);
  auto &params = native->params;
  check_embed(mol, data, params);
  params.numThreads = 1;
  MatrixCacheGuard guard(mol);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  // Pruning by RMS symmetrizes terminal groups (getMolSelfMatches).
  std::unique_lock<std::mutex> symmetrize(symmetrize_mutex(), std::defer_lock);
  if (params.pruneRmsThresh > 0.0 && params.useSymmetryForPruning &&
      params.symmetrizeConjugatedTerminalGroupsForPruning) {
    symmetrize.lock();
  }
  std::unique_lock<std::mutex> random(random_mutex(), std::defer_lock);
  if (params.randomSeed < 0) {
    random.lock();
  }
  RDKit::INT_VECT res;
  DG::EmbedMultipleConfs(mol, res, num_confs, params);
  // Python's parameters show what the embedder changed: maxIterations,
  // basinThresh, initialEmbeddingMode, and the failure counts.
  auto num_threads = data.num_threads;
  set_fields(data, params);
  data.num_threads = num_threads;
  rust::Vec<std::int32_t> ids;
  for (auto id : res) {
    ids.push_back(id);
  }
  return ids;
}

// As getMolBoundsMatrix in rdDistGeom.cpp, which sizes its copy in unsigned
// int. setTopolBounds reads the distance matrix (U6).
rust::Vec<double> get_molecule_bounds_matrix(const RWMol &mol, bool set15bounds,
                                             bool scale_vdw,
                                             bool do_triangle_smoothing,
                                             bool use_macrocycle14config) {
  check_matrix_size(mol);
  MatrixCacheGuard guard(mol);
  unsigned n = mol.getNumAtoms();
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(n));
  DG::initBoundsMat(mat);
  DG::setTopolBounds(mol, mat, set15bounds, scale_vdw, use_macrocycle14config);
  if (do_triangle_smoothing) {
    DistGeom::triangleSmoothBounds(mat);
  }
  rust::Vec<double> values;
  const double *matrix = mat->getData();
  std::size_t size = static_cast<std::size_t>(n) * n;
  values.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    values.push_back(matrix[i]);
  }
  return values;
}

// As getExpTorsHelper in rdDistGeom.cpp. RDKit raises ValueError for
// versions other than 1 and 2 when it reads the torsion preferences.
rust::Vec<ExperimentalTorsionData> get_experimental_torsions(
    const RWMol &mol, bool use_exp_torsion_angle_prefs,
    bool use_small_ring_torsions, bool use_macrocycle_torsions,
    bool use_basic_knowledge, std::uint32_t et_version, bool verbose) {
  namespace CFF = ForceFields::CrystalFF;
  CFF::CrystalFFDetails details;
  std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                         const CFF::ExpTorsionAngle *>>
      torsions;
  // It matches the torsion patterns (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  CFF::getExperimentalTorsions(mol, details, torsions,
                               use_exp_torsion_angle_prefs,
                               use_small_ring_torsions, use_macrocycle_torsions,
                               use_basic_knowledge, et_version, verbose);
  rust::Vec<ExperimentalTorsionData> result;
  for (const auto &[bond, atoms, angle] : torsions) {
    ExperimentalTorsionData item{};
    item.bond_index = bond;
    item.atom_indices = to_vec(atoms);
    item.torsion_index = angle->torsionIdx;
    item.smarts = std::make_unique<std::string>(angle->smarts);
    item.v = to_vec(angle->V);
    item.signs = to_vec(angle->signs);
    result.push_back(std::move(item));
  }
  return result;
}

}  // namespace rdkit_cxx
