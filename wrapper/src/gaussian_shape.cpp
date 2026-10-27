// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/gaussian_shape.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>
#include <boost/dynamic_bitset.hpp>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/gaussian_shape.rs.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/relevant_cycles.h"

namespace rdkit_cxx {

namespace {

namespace GS = RDKit::GaussianShape;

GS::ShapeInputOptions to_options(const GaussianShapeInputOptionsData &data) {
  std::size_t features = 0;
  for (auto count : data.custom_feature_counts) {
    features += count;
  }
  std::size_t atoms = 0;
  for (auto count : data.custom_feature_atom_counts) {
    atoms += count;
  }
  if (data.atom_radii_atoms.size() != data.atom_radii_values.size() ||
      data.custom_feature_types.size() != features ||
      data.custom_feature_positions.size() != features ||
      data.custom_feature_radii.size() != features ||
      data.custom_feature_atom_counts.size() != features ||
      data.custom_feature_atoms.size() != atoms) {
    throw ValueErrorException("mismatched ShapeInputOptions vectors");
  }
  GS::ShapeInputOptions opts;
  opts.useColors = data.use_colors;
  std::size_t f = 0;
  std::size_t a = 0;
  for (auto count : data.custom_feature_counts) {
    std::vector<GS::CustomFeature> conf_features;
    for (std::uint32_t i = 0; i < count; ++i, ++f) {
      const auto &p = data.custom_feature_positions[f];
      std::vector<unsigned int> feature_atoms(
          data.custom_feature_atoms.begin() + a,
          data.custom_feature_atoms.begin() + a +
              data.custom_feature_atom_counts[f]);
      a += data.custom_feature_atom_counts[f];
      conf_features.emplace_back(data.custom_feature_types[f],
                                 RDGeom::Point3D(p.x, p.y, p.z),
                                 data.custom_feature_radii[f], feature_atoms);
    }
    opts.customFeatures.push_back(std::move(conf_features));
  }
  opts.atomSubset.assign(data.atom_subset.begin(), data.atom_subset.end());
  for (std::size_t i = 0; i < data.atom_radii_atoms.size(); ++i) {
    opts.atomRadii.emplace_back(data.atom_radii_atoms[i],
                                data.atom_radii_values[i]);
  }
  opts.allCarbonRadii = data.all_carbon_radii;
  opts.shapePruneThreshold = data.shape_prune_threshold;
  opts.sortShapes = data.sort_shapes;
  opts.includeDummies = data.include_dummies;
  return opts;
}

GS::ShapeOverlayOptions to_overlay(const ShapeOverlayOptionsData &data) {
  GS::ShapeOverlayOptions opts;
  opts.startMode = data.start_mode;
  opts.optimMode = data.optim_mode;
  opts.simAlpha = data.sim_alpha;
  opts.simBeta = data.sim_beta;
  opts.optParam = data.opt_param;
  opts.nSteps = data.n_steps;
  opts.normalize = data.normalize;
  opts.useDistCutoff = data.use_dist_cutoff;
  opts.distCutoff = data.dist_cutoff;
  opts.shapeConvergenceCriterion = data.shape_convergence_criterion;
  return opts;
}

// The overlaps' Gaussians have widths KAPPA / r^2 for radius r, which they
// multiply, add, and divide: radii in [1e-3, 1e3] in magnitude and
// coordinates within 1e8 keep every volume finite. A NaN volume would reach
// sortShapesByVolumes' std::sort, whose comparison is then not a strict
// weak order (undefined; libstdc++ can read past the vector).
void check_radius(double radius) {
  if (!(std::abs(radius) >= 1e-3 && std::abs(radius) <= 1e3)) {
    throw ValueErrorException(
        "ShapeInput needs finite radii from 0.001 to 1000 in magnitude");
  }
}

void check_point(const RDGeom::Point3D &p) {
  for (double v : {p.x, p.y, p.z}) {
    if (!(std::abs(v) <= 1e8)) {
      throw ValueErrorException(
          "ShapeInput needs finite coordinates within 1e8");
    }
  }
}

// What ShapeInput's constructor (ShapeInput.cpp) trusts:
// - atomRadii's atoms, which it looks up, and atomSubset's, which it sets in
//   a dynamic_bitset unchecked (a write past its end);
// - custom features, which it appends to each conformer's coordinates but
//   counts once, from the first conformer, so other conformers with fewer
//   are read past their ends (a segfault in Python);
// - an atom with volume: without one, every overlay score is 0 / 0, which
//   its std::sort of starting poses then compares, and its best score is
//   returned uninitialized;
// - the finite radii and coordinates above.
// It also stores atomRadii as properties on the molecule, through a const
// reference, where later shapes of the molecule find them; the wrapper hands
// it a copy instead, which this returns.
std::unique_ptr<RWMol> checked_mol(const RWMol &mol, int conf_id,
                                   const GS::ShapeInputOptions &opts) {
  // ShapeInput removes atoms from copies (extractSubset, removeHs) and
  // writes SMILES, which trusts what check_removals checks.
  check_removals(mol);
  if (!mol.getNumConformers()) {
    throw ValueErrorException(
        "ShapeInput needs the molecule to have conformers");
  }
  auto num_atoms = mol.getNumAtoms();
  for (const auto &[idx, radius] : opts.atomRadii) {
    check_index(idx, num_atoms);
    check_radius(radius);
  }
  for (auto idx : opts.atomSubset) {
    check_index(idx, num_atoms);
  }
  std::vector<const RDKit::Conformer *> confs;
  if (conf_id >= 0) {
    confs.push_back(&mol.getConformer(conf_id));
  } else {
    for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
      confs.push_back(it->get());
    }
  }
  for (const auto *conf : confs) {
    for (const auto &pos : conf->getPositions()) {
      check_point(pos);
    }
  }
  bool has_volume = false;
  auto counts = [&](const RDKit::Atom *atom) {
    auto num = atom->getAtomicNum();
    return num != 1 && (opts.includeDummies || num != 0);
  };
  if (opts.atomSubset.empty()) {
    for (const auto *atom : mol.atoms()) {
      has_volume = has_volume || counts(atom);
    }
  } else {
    for (auto idx : opts.atomSubset) {
      has_volume = has_volume || counts(mol.getAtomWithIdx(idx));
    }
  }
  if (!has_volume) {
    throw ValueErrorException(
        "ShapeInput needs an atom other than hydrogen (or a dummy, with "
        "includeDummies)");
  }
  if (opts.useColors && !opts.customFeatures.empty()) {
    auto used = std::min(confs.size(), opts.customFeatures.size());
    for (std::size_t k = 0; k < used; ++k) {
      if (opts.customFeatures[k].size() != opts.customFeatures[0].size()) {
        throw ValueErrorException(
            "ShapeInput needs the same number of custom features for every "
            "conformer");
      }
    }
    for (const auto &conf_features : opts.customFeatures) {
      for (const auto &feature : conf_features) {
        check_point(feature.pos);
        check_radius(feature.rad);
      }
    }
  }
  auto copy = std::make_unique<RWMol>(mol);
  if (!opts.atomRadii.empty()) {
    // As the constructor will, and then the radii it will read.
    for (const auto &[idx, radius] : opts.atomRadii) {
      copy->getAtomWithIdx(idx)->setProp<double>("BespokeRadius", radius);
    }
    for (const auto *atom : copy->atoms()) {
      if (atom->hasProp("BespokeRadius")) {
        check_radius(atom->getProp<double>("BespokeRadius"));
      }
    }
  }
  return copy;
}

// The constructor writes the molecule, or what is left of it after
// removing the atoms outside atomSubset (extractSubset in ShapeInput.cpp),
// with MolToSmiles, which perceives its stereochemistry.
void check_shape_relevant_cycles(const RWMol &mol,
                                 const GS::ShapeInputOptions &opts) {
  if (opts.atomSubset.empty() ||
      independent_cycles(mol) <= MAX_INDEPENDENT_CYCLES) {
    check_stereo_relevant_cycles(mol);
    return;
  }
  RWMol subset(mol);
  boost::dynamic_bitset<> keep(mol.getNumAtoms());
  for (auto idx : opts.atomSubset) {
    keep.set(idx);
  }
  subset.beginBatchEdit();
  for (const auto *atom : subset.atoms()) {
    if (!keep[atom->getIdx()]) {
      subset.removeAtom(atom->getIdx());
    }
  }
  subset.commitBatchEdit();
  check_stereo_relevant_cycles(subset);
}

// shapeToMol, which the A_LA_PUBCHEM start mode reaches through
// calcMomentsOfInertia, parses the shape's SMILES and dereferences the
// molecule it gets, null if the SMILES does not parse.
void check_smiles(const GaussianShapeInput &shape) {
  RDKit::v2::SmilesParse::SmilesParserParams params;
  params.sanitize = false;
  if (!RDKit::v2::SmilesParse::MolFromSmiles(shape.getSmiles(), params)) {
    throw ValueErrorException("the shape's SMILES does not parse");
  }
}

GaussianShapeInput make_shape(const RWMol &mol, int conf_id,
                              const GS::ShapeInputOptions &opts,
                              const GS::ShapeOverlayOptions &overlay_opts) {
  auto copy = checked_mol(mol, conf_id, opts);
  check_shape_relevant_cycles(*copy, opts);
  std::lock_guard<std::mutex> lock(gaussian_shape_mutex());
  // Feature patterns are matched (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  if (opts.shapePruneThreshold > 0.0 && conf_id < 0 &&
      copy->getNumConformers() > 1) {
    // Pruning overlays the conformers' shapes in the A_LA_PUBCHEM mode, so
    // the SMILES must parse; it is checked on a shape made without pruning.
    auto unpruned = opts;
    unpruned.shapePruneThreshold = -1.0;
    unpruned.sortShapes = false;
    check_smiles(GaussianShapeInput(*copy, conf_id, unpruned, overlay_opts));
  }
  return GaussianShapeInput(*copy, conf_id, opts, overlay_opts);
}

// What alignShape (GaussianShape.cpp) trusts, besides the shapes'
// construction: it sorts its starting poses by score with std::sort and
// returns its best score uninitialized unless some pose scores above -1, so
// every score must be a number of at least 0. With volumes above 0, finite
// overlaps of at most half their sum, Tversky weights in [0, 1] (not both 0)
// and a color weight in [0, 1], each score lies in [0, 1]. It runs
// nSteps - 10 further steps, in unsigned arithmetic, which wraps for
// nSteps below 10.
void check_overlay(const GS::ShapeOverlayOptions &opts) {
  if (!(opts.simAlpha >= 0.0 && opts.simAlpha <= 1.0 && opts.simBeta >= 0.0 &&
        opts.simBeta <= 1.0 && (opts.simAlpha > 0.0 || opts.simBeta > 0.0))) {
    throw ValueErrorException(
        "shape overlays need simAlpha and simBeta in [0, 1], not both 0");
  }
  if (!(opts.optParam >= 0.0 && opts.optParam <= 1.0)) {
    throw ValueErrorException("shape overlays need optParam in [0, 1]");
  }
  if (opts.nSteps < 10) {
    throw ValueErrorException("shape overlays need nSteps of at least 10");
  }
}

void check_volume(double volume) {
  if (!(volume > 0.0 && std::isfinite(volume))) {
    throw ValueErrorException("shape overlays need shapes with volume");
  }
}

void check_alignable(const GaussianShapeInput &shape,
                     const GS::ShapeOverlayOptions &opts, bool all_shapes) {
  if (all_shapes) {
    for (unsigned i = 0; i < shape.getNumShapes(); ++i) {
      check_volume(shape.getShapeVolume(i));
    }
  } else {
    check_volume(shape.getShapeVolume());
  }
  if (opts.startMode == GS::StartMode::A_LA_PUBCHEM) {
    check_smiles(shape);
  }
}

std::array<double, 3> align(const GaussianShapeInput &ref,
                            GaussianShapeInput &fit, RDGeom::Transform3D &xform,
                            const GS::ShapeOverlayOptions &opts) {
  check_overlay(opts);
  check_alignable(ref, opts, false);
  check_alignable(fit, opts, false);
  return GS::AlignShape(ref, fit, &xform, opts);
}

void push_transform(rust::Vec<double> &values,
                    const RDGeom::Transform3D &xform) {
  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      values.push_back(xform.getValUnchecked(i, j));
    }
  }
}

template <std::size_t N>
rust::Vec<double> to_vec(const std::array<double, N> &values) {
  rust::Vec<double> result;
  for (double v : values) {
    result.push_back(v);
  }
  return result;
}

// As AlignMolecule(refShape, fit, ...) in GaussianShape.cpp.
std::array<double, 3> align_molecule(const GaussianShapeInput &ref_shape,
                                     RWMol &fit,
                                     const GS::ShapeInputOptions &fit_opts,
                                     const GS::ShapeOverlayOptions &opts,
                                     int fit_conf_id) {
  auto fit_shape = make_shape(fit, fit_conf_id, fit_opts, opts);
  RDGeom::Transform3D xform;
  auto scores = align(ref_shape, fit_shape, xform, opts);
  MolTransforms::transformConformer(fit.getConformer(fit_conf_id), xform);
  return scores;
}

rust::Vec<double> score_with_volumes(const GaussianShapeInput &ref,
                                     const GaussianShapeInput &fit,
                                     const GS::ShapeOverlayOptions &opts) {
  std::pair<double, double> volumes;
  auto result = to_vec(GS::ScoreShape(ref, fit, opts, &volumes));
  result.push_back(volumes.first);
  result.push_back(volumes.second);
  return result;
}

}  // namespace

std::mutex &gaussian_shape_mutex() {
  check_not_in_control_c_callback();
  static std::mutex mutex;
  return mutex;
}

void check_shape_coordinates(const RDKit::ROMol &mol) {
  for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
    for (const auto &pos : (*it)->getPositions()) {
      check_point(pos);
    }
  }
}

GaussianShapeInputOptionsData gaussian_shape_input_options_default() {
  GS::ShapeInputOptions opts;
  GaussianShapeInputOptionsData data;
  data.use_colors = opts.useColors;
  data.all_carbon_radii = opts.allCarbonRadii;
  data.shape_prune_threshold = opts.shapePruneThreshold;
  data.sort_shapes = opts.sortShapes;
  data.include_dummies = opts.includeDummies;
  return data;
}

ShapeOverlayOptionsData shape_overlay_options_default() {
  GS::ShapeOverlayOptions opts;
  return ShapeOverlayOptionsData{
      opts.startMode,  opts.optimMode,
      opts.simAlpha,   opts.simBeta,
      opts.optParam,   opts.nSteps,
      opts.normalize,  opts.useDistCutoff,
      opts.distCutoff, opts.shapeConvergenceCriterion};
}

std::unique_ptr<GaussianShapeInput> gaussian_shape_input_new(
    const RWMol &mol, std::int32_t conf_id,
    const GaussianShapeInputOptionsData &shape_opt,
    const ShapeOverlayOptionsData &overlay_opts) {
  return std::make_unique<GaussianShapeInput>(make_shape(
      mol, conf_id, to_options(shape_opt), to_overlay(overlay_opts)));
}

std::unique_ptr<GaussianShapeInput> gaussian_shape_input_clone(
    const GaussianShapeInput &shape) {
  return std::make_unique<GaussianShapeInput>(shape);
}

std::unique_ptr<std::string> gaussian_shape_input_get_smiles(
    const GaussianShapeInput &shape) {
  return std::make_unique<std::string>(shape.getSmiles());
}

void gaussian_shape_input_set_active_shape(GaussianShapeInput &shape,
                                           std::uint32_t shape_num) {
  check_index(shape_num, shape.getNumShapes());
  shape.setActiveShape(shape_num);
}

std::uint32_t gaussian_shape_input_get_active_shape(
    const GaussianShapeInput &shape) {
  return shape.getActiveShape();
}

std::uint32_t gaussian_shape_input_num_atoms(const GaussianShapeInput &shape) {
  return shape.getNumAtoms();
}

std::uint32_t gaussian_shape_input_num_features(
    const GaussianShapeInput &shape) {
  return shape.getNumFeatures();
}

std::uint32_t gaussian_shape_input_num_shapes(const GaussianShapeInput &shape) {
  return shape.getNumShapes();
}

double gaussian_shape_input_shape_volume(const GaussianShapeInput &shape) {
  return shape.getShapeVolume();
}

double gaussian_shape_input_color_volume(const GaussianShapeInput &shape) {
  return shape.getColorVolume();
}

void gaussian_shape_input_normalize_coords(GaussianShapeInput &shape) {
  shape.normalizeCoords();
}

std::unique_ptr<RWMol> gaussian_shape_input_shape_to_mol(
    const GaussianShapeInput &shape, bool include_colors, bool with_bonds) {
  if (with_bonds) {
    check_smiles(shape);
  }
  return shape.shapeToMol(include_colors, with_bonds);
}

rust::Vec<double> gaussian_shape_input_best_similarity(
    GaussianShapeInput &shape, const GaussianShapeInput &fit_shape,
    double threshold, const ShapeOverlayOptionsData &overlay_opts) {
  auto opts = to_overlay(overlay_opts);
  check_overlay(opts);
  check_alignable(shape, opts, true);
  check_alignable(fit_shape, opts, true);
  unsigned best_this = 0;
  unsigned best_fit = 0;
  RDGeom::Transform3D xform;
  auto scores = shape.bestSimilarity(fit_shape, best_this, best_fit, xform,
                                     threshold, opts);
  auto result = to_vec(scores);
  result.push_back(best_this);
  result.push_back(best_fit);
  push_transform(result, xform);
  return result;
}

double gaussian_shape_input_max_possible_similarity(
    const GaussianShapeInput &shape, const GaussianShapeInput &fit_shape,
    const ShapeOverlayOptionsData &overlay_opts) {
  return shape.maxPossibleSimilarity(fit_shape, to_overlay(overlay_opts));
}

// As AlignMolecule(ref, fit, ...).
rust::Vec<double> gaussian_align_mol(
    const RWMol &ref, RWMol &fit, const GaussianShapeInputOptionsData &ref_opts,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t ref_conf_id,
    std::int32_t fit_conf_id) {
  auto opts = to_overlay(overlay_opts);
  check_overlay(opts);
  auto ref_shape = make_shape(ref, ref_conf_id, to_options(ref_opts), opts);
  return to_vec(
      align_molecule(ref_shape, fit, to_options(fit_opts), opts, fit_conf_id));
}

rust::Vec<double> gaussian_align_mol_shape_input(
    const GaussianShapeInput &ref_shape, RWMol &fit,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t fit_conf_id) {
  auto opts = to_overlay(overlay_opts);
  check_overlay(opts);
  return to_vec(
      align_molecule(ref_shape, fit, to_options(fit_opts), opts, fit_conf_id));
}

rust::Vec<double> gaussian_align_shapes(
    const GaussianShapeInput &ref_shape, GaussianShapeInput &fit_shape,
    const ShapeOverlayOptionsData &overlay_opts) {
  RDGeom::Transform3D xform;
  auto result =
      to_vec(align(ref_shape, fit_shape, xform, to_overlay(overlay_opts)));
  push_transform(result, xform);
  return result;
}

// As ScoreMolecule(ref, fit, ...): the shapes are scored where they are.
rust::Vec<double> gaussian_score_mol(
    const RWMol &ref, const RWMol &fit,
    const GaussianShapeInputOptionsData &ref_opts,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t ref_conf_id,
    std::int32_t fit_conf_id) {
  auto opts = to_overlay(overlay_opts);
  opts.normalize = false;
  opts.startMode = GS::StartMode::ROTATE_0;
  auto ref_shape = make_shape(ref, ref_conf_id, to_options(ref_opts), opts);
  auto fit_shape = make_shape(fit, fit_conf_id, to_options(fit_opts), opts);
  return score_with_volumes(ref_shape, fit_shape, opts);
}

// As ScoreMolecule(refShape, fit, ...).
rust::Vec<double> gaussian_score_mol_shape_input(
    const GaussianShapeInput &ref_shape, const RWMol &fit,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts, std::int32_t fit_conf_id) {
  auto opts = to_overlay(overlay_opts);
  auto fit_shape = make_shape(fit, fit_conf_id, to_options(fit_opts), opts);
  return score_with_volumes(ref_shape, fit_shape, opts);
}

rust::Vec<double> gaussian_score_shape(
    const GaussianShapeInput &ref_shape, const GaussianShapeInput &fit_shape,
    const ShapeOverlayOptionsData &overlay_opts) {
  return score_with_volumes(ref_shape, fit_shape, to_overlay(overlay_opts));
}

// As ScoreMoleculeAllConformers, which leaves its best conformers
// uninitialized unless some score beats -1.
rust::Vec<double> gaussian_score_molecule_all_conformers(
    const RWMol &ref, const RWMol &fit,
    const GaussianShapeInputOptionsData &ref_opts,
    const GaussianShapeInputOptionsData &fit_opts,
    const ShapeOverlayOptionsData &overlay_opts) {
  auto opts = to_overlay(overlay_opts);
  check_overlay(opts);
  auto ref_options = to_options(ref_opts);
  ref_options.shapePruneThreshold = -1;
  ref_options.sortShapes = false;
  auto fit_options = to_options(fit_opts);
  fit_options.shapePruneThreshold = -1;
  fit_options.sortShapes = false;
  auto ref_shape = make_shape(ref, -1, ref_options, opts);
  auto fit_shape = make_shape(fit, -1, fit_options, opts);
  check_alignable(ref_shape, opts, true);
  check_alignable(fit_shape, opts, true);
  unsigned num_ref = ref_shape.getNumShapes();
  unsigned num_fit = fit_shape.getNumShapes();
  rust::Vec<double> result;
  result.push_back(num_ref);
  result.push_back(num_fit);
  double best = -1.0;
  int best_ref = -1;
  int best_fit = -1;
  RDGeom::Transform3D best_xform;
  for (unsigned i = 0; i < num_ref; ++i) {
    ref_shape.setActiveShape(i);
    for (unsigned j = 0; j < num_fit; ++j) {
      fit_shape.setActiveShape(j);
      RDGeom::Transform3D xform;
      auto scores = GS::AlignShape(ref_shape, fit_shape, &xform, opts);
      result.push_back(scores[0]);
      if (scores[0] > best) {
        best = scores[0];
        best_ref = static_cast<int>(i);
        best_fit = static_cast<int>(j);
        best_xform = xform;
      }
    }
  }
  result.push_back(best_ref);
  result.push_back(best_fit);
  push_transform(result, best_xform);
  return result;
}

}  // namespace rdkit_cxx
