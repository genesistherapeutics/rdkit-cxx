// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/free_sasa.h"

#include <climits>
#include <cmath>
#include <string>
#include <vector>

#include <Geometry/point.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/RDFreeSASA.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/free_sasa.rs.h"

namespace rdkit_cxx {

namespace {

FreeSASA::SASAOpts to_opts(const SASAOptsData &data) {
  FreeSASA::SASAOpts opts;
  switch (data.algorithm) {
    case FreeSASA::SASAOpts::LeeRichards:
    case FreeSASA::SASAOpts::ShrakeRupley:
      opts.algorithm =
          static_cast<FreeSASA::SASAOpts::Algorithm>(data.algorithm);
      break;
    default:
      throw ValueErrorException("unknown SASAAlgorithm");
  }
  switch (data.classifier) {
    case FreeSASA::SASAOpts::Protor:
    case FreeSASA::SASAOpts::NACCESS:
    case FreeSASA::SASAOpts::OONS:
      opts.classifier =
          static_cast<FreeSASA::SASAOpts::Classifier>(data.classifier);
      break;
    default:
      throw ValueErrorException("unknown SASAClassifier");
  }
  opts.probeRadius = data.probe_radius;
  return opts;
}

// freesasa's classifiers look a residue's and an atom's names up
// (find_string in classifier.c) after trimming them with sscanf's "%s",
// which writes nothing for a name of only whitespace, so strcmp then reads
// an uninitialized buffer, past its end unless it holds a NUL. RDKit passes
// the names as C strings, which end at an embedded NUL. A name without a
// printable ASCII character there matches no classifier entry, so while
// classifyAtoms runs, each such name is replaced by one that matches none
// either, and restored afterwards, even when it throws.
bool blank(const std::string &name) {
  for (char c : name) {
    if (c == '\0') {
      break;
    }
    if (c > ' ' && c < 0x7f) {
      return false;
    }
  }
  return true;
}

class BlankNameGuard {
 public:
  explicit BlankNameGuard(RWMol &mol) {
    for (auto *atom : mol.atoms()) {
      // classifyAtoms reads the names only of PDB residue infos, which
      // U10 guarantees have that class.
      auto *info =
          dynamic_cast<RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
      if (!info) {
        continue;
      }
      bool name = blank(info->getName());
      bool residue = blank(info->getResidueName());
      if (!name && !residue) {
        continue;
      }
      d_saved.push_back(
          Saved{info, info->getName(), info->getResidueName(), name, residue});
      if (name) {
        info->setName(kStandIn);
      }
      if (residue) {
        info->setResidueName(kStandIn);
      }
    }
  }
  ~BlankNameGuard() {
    for (auto &saved : d_saved) {
      if (saved.name_blank) {
        saved.info->setName(saved.name);
      }
      if (saved.residue_blank) {
        saved.info->setResidueName(saved.residue);
      }
    }
  }
  BlankNameGuard(const BlankNameGuard &) = delete;
  BlankNameGuard &operator=(const BlankNameGuard &) = delete;

 private:
  static constexpr const char *kStandIn = "?";
  struct Saved {
    RDKit::AtomPDBResidueInfo *info;
    std::string name;
    std::string residue;
    bool name_blank;
    bool residue_blank;
  };
  std::vector<Saved> d_saved;
};

// freesasa's neighbor lists (freesasa_nb_new in nb.c) divide the box around
// the coordinates, padded by half a cell, into cubic cells twice the largest
// radius plus probe radius wide. They cast each axis's number of cells, and
// each atom's cell along it, to int, multiply the numbers in int, and index
// the cells with the atom's cell, all unchecked. A NaN or infinite
// coordinate, a largest radius plus probe that is not positive and finite
// (as when every radius plus the probe is at most 0), a box of more than
// INT_MAX cells, or a cell too small to change coordinates that large makes
// that undefined, or writes past the cells. This repeats the computation in
// the same double arithmetic (fmin and fmax skip NaN, as in nb.c) and
// refuses the input unless the cells fit in int and every atom lies at
// least a quarter cell inside the box, where nb.c puts it half a cell in.
void check_cell_list(const RDGeom::POINT3D_VECT &positions,
                     const std::vector<double> &radii, double probe_radius) {
  double max_radius = 0.0;
  for (double radius : radii) {
    max_radius = std::fmax(radius + probe_radius, max_radius);
  }
  double d = 2 * max_radius;
  if (!(d > 0.0 && std::isfinite(d))) {
    throw ValueErrorException(
        "CalcSASA needs the largest radius plus the probe radius to be "
        "positive and finite");
  }
  double cells = 1.0;
  for (int axis = 0; axis < 3; ++axis) {
    double lo = positions[0][axis];
    double hi = lo;
    for (const auto &pos : positions) {
      lo = std::fmin(pos[axis], lo);
      hi = std::fmax(pos[axis], hi);
    }
    double min = lo - d / 2.;
    double max = hi + d / 2.;
    double num = std::ceil((max - min) / d);
    if (!(num >= 1.0 && num <= INT_MAX)) {
      throw ValueErrorException(
          "CalcSASA needs finite coordinates spanning at most INT_MAX of "
          "freesasa's cells");
    }
    cells *= num;
    for (const auto &pos : positions) {
      double t = (pos[axis] - min) / d;
      if (!(t >= 0.25 && t <= num - 0.25)) {
        throw ValueErrorException(
            "CalcSASA needs finite coordinates that freesasa's cells, "
            "twice the largest radius plus probe radius wide, can tell "
            "apart");
      }
    }
  }
  if (cells > INT_MAX) {
    throw ValueErrorException(
        "CalcSASA needs finite coordinates spanning at most INT_MAX of "
        "freesasa's cells");
  }
}

std::unique_ptr<Atom> own(const RDKit::QueryAtom *atom) {
  return std::unique_ptr<Atom>(const_cast<RDKit::QueryAtom *>(atom));
}

}  // namespace

SASAOptsData sasa_opts_default() {
  FreeSASA::SASAOpts opts;
  return SASAOptsData{static_cast<std::int32_t>(opts.algorithm),
                      static_cast<std::int32_t>(opts.classifier),
                      opts.probeRadius};
}

rust::Vec<double> classify_atoms(RWMol &mol, const SASAOptsData &opts) {
  auto options = to_opts(opts);
  std::vector<double> radii;
  bool classified;
  {
    BlankNameGuard guard(mol);
    classified = FreeSASA::classifyAtoms(mol, radii, options);
  }
  rust::Vec<double> result;
  if (classified) {
    result.reserve(radii.size());
    for (double radius : radii) {
      result.push_back(radius);
    }
  }
  return result;
}

double calc_sasa(const RWMol &mol, rust::Slice<const double> radii,
                 std::int32_t conf_idx, const Atom *query,
                 const SASAOptsData &opts) {
  auto options = to_opts(opts);
  const RDKit::QueryAtom *query_atom = nullptr;
  if (query) {
    query_atom = dynamic_cast<const RDKit::QueryAtom *>(query);
    if (!query_atom) {
      throw ValueErrorException("Query is not a query atom!");
    }
  }
  std::vector<double> vradii(radii.begin(), radii.end());
  // Without conformers or atoms, RDKit's preconditions throw first.
  if (mol.getNumConformers() && mol.getNumAtoms()) {
    const auto &conf = mol.getConformer(conf_idx);
    // RDKit hands freesasa &radii[0], which reads a radius for every atom.
    if (vradii.size() != mol.getNumAtoms()) {
      throw ValueErrorException("CalcSASA needs a radius for each atom");
    }
    check_cell_list(conf.getPositions(), vradii, options.probeRadius);
  }
  return FreeSASA::calcSASA(mol, vradii, conf_idx, query_atom, options);
}

std::unique_ptr<Atom> make_free_sasa_a_polar_atom_query() {
  return own(FreeSASA::makeFreeSasaAPolarAtomQuery());
}

std::unique_ptr<Atom> make_free_sasa_polar_atom_query() {
  return own(FreeSASA::makeFreeSasaPolarAtomQuery());
}

}  // namespace rdkit_cxx
