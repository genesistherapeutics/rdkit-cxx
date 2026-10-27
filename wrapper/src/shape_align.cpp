// Copyright (c) 2026 Steven Kearnes
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

// pubchem-align3d does not install its header (shape_functions.hpp); these
// are its declarations at the commit RDKit 2024.09 builds (daefab3), whose
// library exports them.
namespace Align3D {
void setUseCutOff(bool useCutOff);
double ComputeShapeOverlap(const float *ref_coord,
                           const std::vector<double> &alpha_ref_vector,
                           const std::vector<unsigned> &ref_volume_atom_index,
                           const float *fit_coord,
                           const std::vector<double> &alpha_fit_vector,
                           const std::vector<unsigned> &fit_volume_atom_index);
}  // namespace Align3D

namespace rdkit_cxx {

namespace {

// PubChemShape.cpp and pubchem-align3d keep state in globals they write
// without a lock: AlignMolecule calls Align3D::setUseCutOff(true), which
// assigns the dCutOff and d2CutOff that every overlap reads
// (shape_functions2.cpp), and PrepareConformer fills getPh4Patterns' static
// list of patterns on first use. Every call into them
// holds this mutex (CONVENTIONS.md U4), before control_c_mutex(), so a
// closure that holds that one must not wait for it.
std::mutex &pubchem_shape_mutex() {
  check_not_in_control_c_callback();
  static std::mutex mutex;
  return mutex;
}

constexpr const char *kFeaturesProp = "PUBCHEM_PHARMACOPHORE_FEATURES";

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

// PrepareConformer writes each heavy atom's coordinates at the atom's own
// index in an array sized for the heavy atoms and the features, so a heavy
// atom after a hydrogen or a dummy atom writes past the heavy atoms (onto a
// feature, or past the array's end).
void check_heavy_atoms_first(const RDKit::ROMol &mol) {
  const unsigned heavy = mol.getNumHeavyAtoms();
  for (const auto *atom : mol.atoms()) {
    if (atom->getAtomicNum() > 1 && atom->getIdx() >= heavy) {
      throw ValueErrorException(
          "PrepareConformer needs the heavy atoms before any hydrogen or "
          "dummy atom; atom " +
          std::to_string(atom->getIdx()) + " follows one");
    }
  }
}

// What PrepareConformer trusts, checked. RDKit 2024.09's PrepareConformer
// computes the self overlaps under whatever cutoff the last alignment left
// (none before the first); the wrapper sets the alignment's first, so a
// shape is the same before and after other alignments, as AlignMol makes it.
PubChemShapeInput prepare(const RDKit::ROMol &mol, int conf_id,
                          bool use_colors) {
  if (use_colors) {
    check_features(mol);
  }
  for (const auto &pos : mol.getConformer(conf_id).getPositions()) {
    check_coordinate(pos.x);
    check_coordinate(pos.y);
    check_coordinate(pos.z);
  }
  check_heavy_atoms_first(mol);
  std::lock_guard<std::mutex> lock(pubchem_shape_mutex());
  Align3D::setUseCutOff(true);
  // The color features come from substructure matches (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return PrepareConformer(mol, conf_id, use_colors);
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
// makes some pose score above 0 at the end. A pose's score is opt_param
// times its shape Tanimoto plus the rest times its color Tanimoto, which is
// not negative, so a positive shape overlap in one of the four poses is
// enough; ComputeShapeOverlap computes it with the same terms, exactly,
// since the flips only negate coordinates, under the cutoff the alignment
// sets.
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
    double overlap = Align3D::ComputeShapeOverlap(
        ref.coord.data(), ref.alpha_vector, ref.volumeAtomIndexVector,
        posed.coord.data(), posed.alpha_vector, posed.volumeAtomIndexVector);
    if (positive(overlap) > 0.0) {
      return;
    }
  }
  throw ValueErrorException(
      "shape alignment needs shapes that overlap in a starting pose");
}

// As AlignMolecule(refShape, fit, ...) in PubChemShape.cpp, after the
// checks: the fit's shape is prepared first, as AlignMolecule prepares it,
// and checked against the reference's.
std::pair<double, double> align_molecule(const PubChemShapeInput &ref_shape,
                                         RWMol &fit, int fit_conf_id,
                                         bool use_colors, double opt_param,
                                         unsigned max_preiters,
                                         unsigned max_postiters) {
  auto fit_shape = prepare(fit, fit_conf_id, use_colors);
  std::lock_guard<std::mutex> lock(pubchem_shape_mutex());
  check_alignable(ref_shape, fit_shape, opt_param);
  std::vector<float> matrix(12, 0.0f);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return AlignMolecule(ref_shape, fit, matrix, fit_conf_id, use_colors,
                       opt_param, max_preiters, max_postiters);
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

std::unique_ptr<PubChemShapeInput> pubchem_prepare_conformer(
    const RWMol &mol, std::int32_t conf_id, bool use_colors) {
  return std::make_unique<PubChemShapeInput>(
      prepare(mol, conf_id, use_colors));
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

double pubchem_shape_input_sov(const PubChemShapeInput &shape) {
  return shape.sov;
}

double pubchem_shape_input_sof(const PubChemShapeInput &shape) {
  return shape.sof;
}

// As AlignMolecule(ref, fit, ...), which prepares the reference's shape and
// aligns the fit to it.
rust::Vec<double> pubchem_align_mol(const RWMol &ref, RWMol &probe,
                                    std::int32_t ref_conf_id,
                                    std::int32_t probe_conf_id,
                                    bool use_colors, double opt_param,
                                    std::uint32_t max_preiters,
                                    std::uint32_t max_postiters) {
  auto ref_shape = prepare(ref, ref_conf_id, use_colors);
  return pair_vec(align_molecule(ref_shape, probe, probe_conf_id, use_colors,
                                 opt_param, max_preiters, max_postiters));
}

rust::Vec<double> pubchem_align_mol_shape_input(
    const PubChemShapeInput &ref_shape, RWMol &probe,
    std::int32_t probe_conf_id, bool use_colors, double opt_param,
    std::uint32_t max_preiters, std::uint32_t max_postiters) {
  return pair_vec(align_molecule(ref_shape, probe, probe_conf_id, use_colors,
                                 opt_param, max_preiters, max_postiters));
}

}  // namespace rdkit_cxx
