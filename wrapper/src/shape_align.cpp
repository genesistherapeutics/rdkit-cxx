// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/shape_align.h"

#include <cfloat>
#include <cmath>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Geometry/point.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/shape_align.rs.h"

namespace rdkit_cxx {

namespace {

// PubChemShape.cpp and pubchem-align3d keep state in globals they write
// without a lock: PrepareConformer and AlignMolecule call
// Align3D::setUseCutOff(true), which assigns the dCutOff and d2CutOff that
// every overlap reads (shape_functions2.cpp), and PrepareConformer fills
// getPh4Patterns' static list of patterns on first use. Every call into them
// holds this mutex (CONVENTIONS.md U4), before control_c_mutex(), so a
// closure that holds that one must not wait for it.
std::mutex &pubchem_shape_mutex() {
  check_not_in_control_c_callback();
  static std::mutex mutex;
  return mutex;
}

constexpr const char *kFeaturesProp = "PUBCHEM_PHARMACOPHORE_FEATURES";

ShapeInputOptions to_options(const ShapeInputOptionsData &data) {
  if (data.atom_radii_atoms.size() != data.atom_radii_values.size() ||
      data.custom_feature_types.size() !=
          data.custom_feature_positions.size() ||
      data.custom_feature_types.size() != data.custom_feature_radii.size()) {
    throw ValueErrorException("mismatched ShapeInputOptions vectors");
  }
  ShapeInputOptions opts;
  opts.useColors = data.use_colors;
  opts.includeDummies = data.include_dummies;
  opts.dummyRadius = data.dummy_radius;
  opts.atomSubset.assign(data.atom_subset.begin(), data.atom_subset.end());
  opts.notColorAtoms.assign(data.not_color_atoms.begin(),
                            data.not_color_atoms.end());
  for (std::size_t i = 0; i < data.atom_radii_atoms.size(); ++i) {
    opts.atomRadii.emplace_back(data.atom_radii_atoms[i],
                                data.atom_radii_values[i]);
  }
  for (std::size_t i = 0; i < data.custom_feature_types.size(); ++i) {
    const auto &p = data.custom_feature_positions[i];
    opts.customFeatures.emplace_back(data.custom_feature_types[i],
                                     RDGeom::Point3D(p.x, p.y, p.z),
                                     data.custom_feature_radii[i]);
  }
  opts.normalize = data.normalize;
  return opts;
}

void check_coordinate(double v) {
  // PrepareConformer stores the coordinates as float; a double beyond
  // float's range does not convert.
  if (!(std::abs(v) <= 1e8)) {
    throw ValueErrorException(
        "PrepareConformer needs finite coordinates within 1e8");
  }
}

// extractFeatures (PubChemShape.cpp) reads the counts in a molecule's
// PUBCHEM_PHARMACOPHORE_FEATURES property with stoul and resizes vectors to
// them before it checks them against the lines and tokens that follow, so a
// count of billions allocates gigabytes (bad_alloc in Python, or the OOM
// killer). This counts the lines and tokens first, and throws the
// ValueErrorException RDKit would throw after allocating.
void check_features(const RDKit::ROMol &mol) {
  std::string features;
  if (!mol.getPropIfPresent(kFeaturesProp, features)) {
    return;
  }
  std::istringstream iss(features);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(iss, line)) {
    lines.push_back(line);
  }
  if (lines.empty()) {
    return;
  }
  unsigned long count;
  try {
    count = std::stoul(lines[0]);
  } catch (...) {
    return;  // RDKit throws the same
  }
  if (count > lines.size() - 1) {
    throw ValueErrorException("Wrong number of features");
  }
  for (std::size_t f = 1; f <= count; ++f) {
    std::istringstream tokens_in(lines[f]);
    std::vector<std::string> tokens;
    std::string token;
    while (std::getline(tokens_in, token, ' ')) {
      tokens.push_back(token);
    }
    if (tokens.empty()) {
      return;  // RDKit throws before it allocates
    }
    unsigned long atoms;
    try {
      atoms = std::stoul(tokens[0]);
    } catch (...) {
      return;
    }
    if (tokens.size() < 2 || atoms > tokens.size() - 2) {
      throw ValueErrorException("Wrong number of tokens in feature");
    }
  }
}

PubChemShapeInput prepare(const RDKit::ROMol &mol, int conf_id,
                          const ShapeInputOptions &opts) {
  if (opts.useColors) {
    check_features(mol);
  }
  for (const auto &pos : mol.getConformer(conf_id).getPositions()) {
    check_coordinate(pos.x);
    check_coordinate(pos.y);
    check_coordinate(pos.z);
  }
  for (const auto &feature : opts.customFeatures) {
    const auto &pos = std::get<1>(feature);
    check_coordinate(pos.x);
    check_coordinate(pos.y);
    check_coordinate(pos.z);
  }
  std::lock_guard<std::mutex> lock(pubchem_shape_mutex());
  // The color features come from substructure matches (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return PrepareConformer(mol, conf_id, opts);
}

void check_shape(const PubChemShapeInput &shape) {
  if (shape.alpha_vector.empty()) {
    throw ValueErrorException(
        "shape alignment needs shapes with atoms or features");
  }
  for (float c : shape.coord) {
    if (!std::isfinite(c)) {
      throw ValueErrorException(
          "shape alignment needs shapes with finite coordinates");
    }
  }
  for (double alpha : shape.alpha_vector) {
    if (!(alpha > 0.0 && alpha <= 100.0)) {
      throw ValueErrorException(
          "shape alignment needs radii of at least 0.16 (Gaussian widths in "
          "(0, 100])");
    }
  }
}

// Align3D::Neighbor_Conformers (shape_neighbor.cpp) copies three
// coordinates per alpha into a work vector and reads its first point, past
// the end for a shape without atoms or features (a segfault in Python). It
// returns the pose in best_qom, which it sets only when some pose scores
// above 0, so for shapes that never overlap it reads uninitialized memory
// (valgrind traces Python's AlignShapes matrix there). The check: its
// first four starting poses, in both of its branches, are the identity and
// three axis flips; its line searches only increase the shape overlap, and
// a pose without overlap does not move; with radii giving Gaussian widths in
// (0, 100] (so each overlap term lies in (0, 1]), the color weight opt_param
// in (0, 1], and the overlap of two shapes below the sum of their self
// overlaps, as Gaussians keep it, one of those four poses scoring above 0
// makes some pose score above 0 at the end. ScoreShape computes those
// scores with the same overlaps, exactly, since the flips only negate
// coordinates.
void check_alignable(const PubChemShapeInput &ref, const PubChemShapeInput &fit,
                     double opt_param) {
  check_shape(ref);
  check_shape(fit);
  std::set<unsigned> ref_types(ref.atom_type_vector.begin(),
                               ref.atom_type_vector.end());
  bool colors = false;
  for (unsigned type : fit.atom_type_vector) {
    if (type > 0 && ref_types.count(type)) {
      colors = true;
    }
  }
  if (colors && !(opt_param > 0.0 && opt_param <= 1.0)) {
    throw ValueErrorException(
        "shape alignment with colors needs opt_param in (0, 1]");
  }
  static const float flips[4][3] = {
      {1, 1, 1}, {1, -1, -1}, {-1, -1, 1}, {-1, 1, -1}};
  auto positive = [](double v) {
    return std::isfinite(v) && v > 0.0 ? v : 0.0;
  };
  for (const auto &flip : flips) {
    auto posed = fit;
    for (std::size_t i = 0; i < posed.coord.size(); ++i) {
      posed.coord[i] *= flip[i % 3];
    }
    auto [shape, color] = ScoreShape(ref, posed, colors);
    double fst = positive(shape);
    double fct = positive(color);
    double score = colors ? opt_param * (fst - fct) + fct : fst;
    if (score > 0.0) {
      return;
    }
  }
  throw ValueErrorException(
      "shape alignment needs shapes that overlap in a starting pose");
}

// AlignShape, on a copy of the probe that replaces it only on success.
std::pair<double, double> align(const PubChemShapeInput &ref,
                                PubChemShapeInput &fit,
                                std::vector<float> &matrix, double opt_param,
                                unsigned max_preiters, unsigned max_postiters) {
  std::lock_guard<std::mutex> lock(pubchem_shape_mutex());
  check_alignable(ref, fit, opt_param);
  auto posed = fit;
  auto scores =
      AlignShape(ref, posed, matrix, opt_param, max_preiters, max_postiters);
  // Unreachable after check_alignable, short of an overlap above the sum of
  // the self overlaps: no pose scored, and the matrix is not one.
  if (scores.first == 0.0 && scores.second == 0.0) {
    throw ValueErrorException("shape alignment found no overlapping pose");
  }
  fit = std::move(posed);
  return scores;
}

void check_size(std::size_t size, std::size_t expected, const char *what) {
  if (size != expected) {
    throw ValueErrorException(std::string("The ") + what + " must have " +
                              std::to_string(expected) + " values.  It had " +
                              std::to_string(size) + ".");
  }
}

// As AlignMolecule(refShape, fit, ...) in PubChemShape.cpp, step by step,
// with the checks above between the steps.
std::pair<double, double> align_molecule(
    const PubChemShapeInput &ref_shape, RWMol &fit,
    const ShapeInputOptions &opts, int fit_conf_id, double opt_param,
    unsigned max_preiters, unsigned max_postiters, bool apply_ref_shift) {
  auto fit_shape = prepare(fit, fit_conf_id, opts);
  std::vector<float> matrix(12, 0.0f);
  auto scores = align(ref_shape, fit_shape, matrix, opt_param, max_preiters,
                      max_postiters);
  auto &conf = fit.getConformer(fit_conf_id);
  std::vector<double> final_trans{0.0, 0.0, 0.0};
  std::vector<double> final_rot{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  if (apply_ref_shift) {
    const auto &rot = ref_shape.inertialRot;
    final_trans = {-ref_shape.shift[0], -ref_shape.shift[1],
                   -ref_shape.shift[2]};
    final_rot = {rot[0], rot[3], rot[6], rot[1], rot[4],
                 rot[7], rot[2], rot[5], rot[8]};
  }
  TransformConformer(final_trans, final_rot, matrix, fit_shape, conf);
  fit.setProp("shape_align_shape_tanimoto", scores.first);
  fit.setProp("shape_align_color_tanimoto", scores.second);
  return scores;
}

rust::Vec<double> pair_vec(const std::pair<double, double> &scores) {
  rust::Vec<double> result;
  result.push_back(scores.first);
  result.push_back(scores.second);
  return result;
}

template <class T, class U>
rust::Vec<T> copy_vec(const std::vector<U> &v) {
  rust::Vec<T> result;
  result.reserve(v.size());
  for (const auto &x : v) {
    result.push_back(x);
  }
  return result;
}

}  // namespace

ShapeInputOptionsData pubchem_shape_input_options_default() {
  ShapeInputOptions opts;
  ShapeInputOptionsData data;
  data.use_colors = opts.useColors;
  data.include_dummies = opts.includeDummies;
  data.dummy_radius = opts.dummyRadius;
  data.normalize = opts.normalize;
  return data;
}

std::unique_ptr<PubChemShapeInput> pubchem_prepare_conformer(
    const RWMol &mol, std::int32_t conf_id, const ShapeInputOptionsData &opts) {
  return std::make_unique<PubChemShapeInput>(
      prepare(mol, conf_id, to_options(opts)));
}

std::unique_ptr<PubChemShapeInput> pubchem_shape_input_clone(
    const PubChemShapeInput &shape) {
  return std::make_unique<PubChemShapeInput>(shape);
}

rust::Vec<float> pubchem_shape_input_coord(const PubChemShapeInput &shape) {
  return copy_vec<float>(shape.coord);
}

rust::Vec<double> pubchem_shape_input_alpha_vector(
    const PubChemShapeInput &shape) {
  return copy_vec<double>(shape.alpha_vector);
}

rust::Vec<std::uint32_t> pubchem_shape_input_atom_type_vector(
    const PubChemShapeInput &shape) {
  return copy_vec<std::uint32_t>(shape.atom_type_vector);
}

rust::Vec<std::uint32_t> pubchem_shape_input_volume_atom_index_vector(
    const PubChemShapeInput &shape) {
  return copy_vec<std::uint32_t>(shape.volumeAtomIndexVector);
}

rust::Vec<double> pubchem_shape_input_shift(const PubChemShapeInput &shape) {
  return copy_vec<double>(shape.shift);
}

rust::Vec<double> pubchem_shape_input_inertial_rot(
    const PubChemShapeInput &shape) {
  return copy_vec<double>(shape.inertialRot);
}

double pubchem_shape_input_sov(const PubChemShapeInput &shape) {
  return shape.sov;
}

double pubchem_shape_input_sof(const PubChemShapeInput &shape) {
  return shape.sof;
}

// As AlignMolecule(ref, fit, matrix, refShapeOpts, fitShapeOpts, ...).
rust::Vec<double> pubchem_align_mol(
    const RWMol &ref, RWMol &probe, const ShapeInputOptionsData &ref_opts,
    const ShapeInputOptionsData &probe_opts, std::int32_t ref_conf_id,
    std::int32_t probe_conf_id, double opt_param, std::uint32_t max_preiters,
    std::uint32_t max_postiters) {
  auto ref_options = to_options(ref_opts);
  auto probe_options = to_options(probe_opts);
  if (ref_options.useColors != probe_options.useColors) {
    BOOST_LOG(rdWarningLog)
        << "useColor values inconsistent between the reference and fit "
           "molecules. Colors will not be used in the alignment."
        << std::endl;
  }
  auto ref_shape = prepare(ref, ref_conf_id, ref_options);
  return pair_vec(align_molecule(ref_shape, probe, probe_options, probe_conf_id,
                                 opt_param, max_preiters, max_postiters, true));
}

rust::Vec<double> pubchem_align_mol_shape_input(
    const PubChemShapeInput &ref_shape, RWMol &probe,
    std::int32_t probe_conf_id, bool use_colors, double opt_param,
    std::uint32_t max_preiters, std::uint32_t max_postiters,
    bool apply_ref_shift) {
  ShapeInputOptions opts;
  opts.useColors = use_colors;
  return pair_vec(align_molecule(ref_shape, probe, opts, probe_conf_id,
                                 opt_param, max_preiters, max_postiters,
                                 apply_ref_shift));
}

rust::Vec<double> pubchem_align_shapes(const PubChemShapeInput &ref_shape,
                                       PubChemShapeInput &probe_shape,
                                       double opt_param,
                                       std::uint32_t max_preiters,
                                       std::uint32_t max_postiters) {
  std::vector<float> matrix(12, 0.0f);
  auto result = pair_vec(align(ref_shape, probe_shape, matrix, opt_param,
                               max_preiters, max_postiters));
  for (float v : matrix) {
    result.push_back(v);
  }
  return result;
}

void pubchem_transform_conformer_to_shape(rust::Slice<const double> final_trans,
                                          rust::Slice<const double> final_rot,
                                          rust::Slice<const double> matrix,
                                          const PubChemShapeInput &probe_shape,
                                          Conformer &probe_conformer) {
  // As Python's TransformConformer, which checks the sizes RDKit indexes.
  check_size(matrix.size(), 12, "transformation matrix");
  check_size(final_trans.size(), 3, "final translation vector");
  check_size(final_rot.size(), 9, "final rotation vector");
  std::vector<float> values;
  for (double v : matrix) {
    // A finite double beyond float's range does not convert.
    if (std::isfinite(v) && std::abs(v) > FLT_MAX) {
      throw ValueErrorException(
          "The transformation matrix's values must fit in float.");
    }
    values.push_back(static_cast<float>(v));
  }
  TransformConformer(
      std::vector<double>(final_trans.begin(), final_trans.end()),
      std::vector<double>(final_rot.begin(), final_rot.end()), values,
      probe_shape, probe_conformer);
}

// As ScoreMolecule(mol1, mol2, ...): shapes left where they are.
rust::Vec<double> pubchem_score_mol(const RWMol &mol1, const RWMol &mol2,
                                    const ShapeInputOptionsData &mol1_opts,
                                    const ShapeInputOptionsData &mol2_opts,
                                    std::int32_t mol1_conf_id,
                                    std::int32_t mol2_conf_id) {
  auto opts1 = to_options(mol1_opts);
  auto opts2 = to_options(mol2_opts);
  bool use_colors = opts1.useColors && opts2.useColors;
  opts1.normalize = false;
  opts2.normalize = false;
  auto shape1 = prepare(mol1, mol1_conf_id, opts1);
  auto shape2 = prepare(mol2, mol2_conf_id, opts2);
  std::lock_guard<std::mutex> lock(pubchem_shape_mutex());
  return pair_vec(ScoreShape(shape1, shape2, use_colors));
}

// As ScoreMolecule(shape, mol, ...).
rust::Vec<double> pubchem_score_mol_shape_input(
    const PubChemShapeInput &shape, const RWMol &mol,
    const ShapeInputOptionsData &mol_opts, std::int32_t mol_conf_id) {
  auto opts = to_options(mol_opts);
  opts.normalize = false;
  auto shape2 = prepare(mol, mol_conf_id, opts);
  std::lock_guard<std::mutex> lock(pubchem_shape_mutex());
  return pair_vec(ScoreShape(shape, shape2, opts.useColors));
}

rust::Vec<double> pubchem_score_shape(const PubChemShapeInput &shape1,
                                      const PubChemShapeInput &shape2,
                                      bool use_colors) {
  auto copy = shape2;  // ScoreShape takes it by non-const reference
  std::lock_guard<std::mutex> lock(pubchem_shape_mutex());
  return pair_vec(ScoreShape(shape1, copy, use_colors));
}

}  // namespace rdkit_cxx
