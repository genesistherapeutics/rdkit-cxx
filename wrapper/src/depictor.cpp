// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/depictor.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <Geometry/point.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDConfig.h>
#include <unistd.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/depictor.rs.h"
#include "rdkit-cxx/file_checks.h"
#include "rdkit-cxx/global_state.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

RDDepict::ConstrainedDepictionParams to_params(
    const ConstrainedDepictionParamsData &data) {
  RDDepict::ConstrainedDepictionParams params;
  params.acceptFailure = data.accept_failure;
  params.forceRDKit = data.force_rdkit;
  params.allowRGroups = data.allow_r_groups;
  params.alignOnly = data.align_only;
  params.adjustMolBlockWedging = data.adjust_mol_block_wedging;
  params.existingConfId = data.existing_conf_id;
  params.useRingTemplates = data.use_ring_templates;
  return params;
}

// What a constrained depiction reads besides the molecule's graph: the
// reference's coordinates (a coordinate map, or the target of an alignment)
// and, to align or to decide about wedging, the molecule's own. A missing
// conformer is left for RDKit, which raises only if it reads it.
void check_constrained_inputs(const RDKit::ROMol &mol,
                              const RDKit::ROMol &reference, int conf_id,
                              const RDDepict::ConstrainedDepictionParams &p) {
  check_compute_2d_coords_mol(mol, p.forceRDKit);
  if (const auto *conf = find_conformer(reference, conf_id)) {
    check_depiction_conformer(*conf, "the reference's coordinates");
  }
  if (const auto *conf = find_conformer(mol, p.existingConfId)) {
    check_depiction_conformer(*conf, "the molecule's coordinates");
  }
}

const RDKit::Conformer *checked_bonded_conformer(const RDKit::ROMol &mol,
                                                 int conf_id) {
  // Both return before reading the conformer when there is no bond.
  if (mol.getNumBonds() == 0) {
    return nullptr;
  }
  const auto *conf = find_conformer(mol, conf_id);
  if (conf) {
    check_depiction_conformer(*conf, "the coordinates");
  }
  return conf;
}

// A file of checked text that is removed when this goes out of scope.
class TemporaryFile {
 public:
  explicit TemporaryFile(const std::string &text) {
    auto pattern =
        (std::filesystem::temp_directory_path() / "rdkit-cxx-XXXXXX").string();
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    int fd = mkstemp(name.data());
    if (fd < 0) {
      throw RDDepict::DepictException("Could not create a temporary file");
    }
    d_path = name.data();
    std::size_t written = 0;
    while (written < text.size()) {
      auto n = write(fd, text.data() + written, text.size() - written);
      if (n < 0 && errno == EINTR) {
        continue;
      }
      if (n <= 0) {
        close(fd);
        std::remove(d_path.c_str());
        throw RDDepict::DepictException("Could not write a temporary file");
      }
      written += static_cast<std::size_t>(n);
    }
    close(fd);
  }
  ~TemporaryFile() { std::remove(d_path.c_str()); }
  TemporaryFile(const TemporaryFile &) = delete;
  TemporaryFile &operator=(const TemporaryFile &) = delete;
  const std::string &path() const { return d_path; }

 private:
  std::string d_path;
};

std::string replace_all(std::string text, const std::string &from,
                        const std::string &to) {
  for (auto pos = text.find(from); pos != std::string::npos;
       pos = text.find(from, pos + to.size())) {
    text.replace(pos, from.size(), to);
  }
  return text;
}

// RDKit copies a ring system template's coordinates into depictions
// unchecked (EmbeddedFrag::matchToTemplate), and CXSMILES coordinates are
// parsed with atof, so "nan" is one. This reads the file, checks each
// template RDKit will load, and hands RDKit the checked text in a temporary
// file, so RDKit never reads a file that changed after the check.
template <typename Load>
void load_checked_templates(rust::Slice<const std::uint8_t> path, Load load) {
  const auto name = to_path(path);
  std::ifstream in(name, std::ios::binary);
  if (!in) {
    throw RDDepict::DepictException("Could not open file " + name);
  }
  std::stringstream buffer;
  buffer << in.rdbuf();
  const auto text = buffer.str();
  std::istringstream lines(text);
  std::string line;
  // As CoordinateTemplates::loadTemplatesFromPath, which stops at the first
  // line it cannot parse.
  while (std::getline(lines, line)) {
    check_smiles_cx(line);
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(line));
    if (!mol) {
      break;
    }
    for (auto it = mol->beginConformers(); it != mol->endConformers(); ++it) {
      check_depiction_conformer(**it, "ring system template coordinates");
    }
  }
  TemporaryFile checked(text);
  try {
    load(checked.path());
  } catch (const RDDepict::DepictException &e) {
    throw RDDepict::DepictException(
        replace_all(e.what(), checked.path(), name));
  }
}

}  // namespace

bool is_coord_gen_support_available() {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  return true;
#else
  return false;
#endif
}

// As setPreferCoordGen and getPreferCoordGen in Wrap/rdDepictor.cpp, which
// ignore the setting without CoordGen support.
void set_prefer_coord_gen(bool val) { replace_prefer_coord_gen(val); }

bool get_prefer_coord_gen() {
  auto lock = lock_depiction();
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  return RDDepict::preferCoordGen;
#else
  return false;
#endif
}

bool replace_prefer_coord_gen(bool val) {
  auto lock = lock_depiction();
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool old = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = val;
  return old;
#else
  (void)val;
  return false;
#endif
}

void set_ring_system_templates(rust::Slice<const std::uint8_t> path) {
  load_checked_templates(path, [](const std::string &checked) {
    auto lock = lock_depiction();
    RDDepict::setRingSystemTemplates(checked);
  });
}

void add_ring_system_templates(rust::Slice<const std::uint8_t> path) {
  load_checked_templates(path, [](const std::string &checked) {
    auto lock = lock_depiction();
    RDDepict::addRingSystemTemplates(checked);
  });
}

void load_default_ring_system_templates() {
  auto lock = lock_depiction();
  RDDepict::loadDefaultRingSystemTemplates();
}

ConstrainedDepictionParamsData constrained_depiction_params_default() {
  RDDepict::ConstrainedDepictionParams p;
  return ConstrainedDepictionParamsData{
      p.acceptFailure,         p.forceRDKit,
      p.allowRGroups,          p.alignOnly,
      p.adjustMolBlockWedging, p.existingConfId,
      p.useRingTemplates,      -1,
  };
}

// As Compute2DCoords in Wrap/rdDepictor.cpp.
unsigned compute_2d_coords(RWMol &mol, const Compute2DCoordsData &data) {
  RDGeom::INT_POINT2D_MAP coord_map;
  const auto n =
      std::min(data.coord_map_keys.size(), data.coord_map_values.size());
  for (std::size_t i = 0; i < n; ++i) {
    // Python checks the atom index; RDKit's depictor would check it only
    // for a map of two or more atoms (CoordGen not at all).
    check_index(data.coord_map_keys[i], mol.getNumAtoms());
    const auto &p = data.coord_map_values[i];
    check_depiction_coordinate(p.x, "coordMap coordinates");
    check_depiction_coordinate(p.y, "coordMap coordinates");
    coord_map[static_cast<int>(data.coord_map_keys[i])] =
        RDGeom::Point2D(p.x, p.y);
  }
  check_bond_length(data.bond_length);
  auto lock = lock_depiction();
  // Ring-system templates are matched to the molecule (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  // Random bond flips (nSample and nFlipsPerSample) draw from RDKit's
  // generator.
  std::lock_guard<std::mutex> random_lock(random_mutex());
  check_compute_2d_coords_mol(mol, data.force_rdkit);
  MatrixCacheGuard guard(mol);
  BondLengthOverride bond_length(data.bond_length);
  RDDepict::Compute2DCoordParameters params;
  params.coordMap = coord_map.empty() ? nullptr : &coord_map;
  params.canonOrient = data.canon_orient;
  params.clearConfs = data.clear_confs;
  params.nFlipsPerSample = data.n_flips_per_sample;
  params.nSamples = data.n_sample;
  params.sampleSeed = data.sample_seed;
  params.permuteDeg4Nodes = data.permute_deg4_nodes != 0;
  params.forceRDKit = data.force_rdkit;
  params.useRingTemplates = data.use_ring_templates;
  return RDDepict::compute2DCoords(mol, params);
}

// As Compute2DCoordsMimicDistmat in Wrap/rdDepictor.cpp, whose
// n * (n - 1) / 2 is unsigned and wraps for large molecules; RDKit indexes
// the matrix that way too, so the check comes first.
unsigned compute_2d_coords_mimic_distmat(
    RWMol &mol, rust::Slice<const double> dist_mat,
    const Compute2DCoordsMimicDistmatData &data) {
  check_matrix_size(mol);
  const std::size_t n = mol.getNumAtoms();
  const std::size_t items = n * (n - 1) / 2;
  if (dist_mat.size() != items) {
    throw ValueErrorException(
        "The array size does not match the number of atoms in the molecule");
  }
  check_bond_length(data.bond_length);
  auto lock = lock_depiction();
  std::lock_guard<std::mutex> random_lock(random_mutex());
  // compute2DCoordsMimicDistMat always uses RDKit's depictor.
  check_rdkit_depictor_mol(mol);
  MatrixCacheGuard guard(mol);
  BondLengthOverride bond_length(data.bond_length);
  RDDepict::DOUBLE_SMART_PTR dmat(new double[items]);
  std::copy(dist_mat.begin(), dist_mat.end(), dmat.get());
  return RDDepict::compute2DCoordsMimicDistMat(
      mol, &dmat, data.canon_orient, data.clear_confs, data.weight_dist_mat,
      data.n_flips_per_sample, data.n_sample, data.sample_seed,
      data.permute_deg4_nodes, data.force_rdkit);
}

rust::Vec<IndexPair> generate_depiction_matching_2d_structure(
    RWMol &mol, const RWMol &reference, const RWMol &ref_patt,
    bool has_ref_patt, const ConstrainedDepictionParamsData &data) {
  auto params = to_params(data);
  auto lock = lock_depiction();
  // It matches the reference, or its pattern, to the molecule
  // (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  check_constrained_inputs(mol, reference, data.conf_id, params);
  // With allowRGroups, RDKit adjusts a copy of the pattern, or of the
  // reference, with adjustQueryProperties, whose RWMol::replaceBond leaves
  // the copy's stereo groups pointing at the bonds it deletes, and then
  // copies the copy, which reads them (prepareTemplateForRGroups in
  // DepictUtils.cpp). The depictor reads nothing else of those groups, so it
  // gets a copy whose groups hold no bonds.
  const RDKit::ROMol &query = has_ref_patt ? ref_patt : reference;
  if (params.allowRGroups) {
    // It then removes the copy's terminal R groups (check_removals).
    check_removals(query);
  }
  std::unique_ptr<RDKit::RWMol> stripped;
  if (params.allowRGroups && has_stereo_group_bonds(query)) {
    stripped = std::make_unique<RDKit::RWMol>(query);
    clear_stereo_group_bonds(*stripped);
  }
  const RDKit::ROMol &ref = stripped && !has_ref_patt ? *stripped : reference;
  const RDKit::ROMol *patt = nullptr;
  if (has_ref_patt) {
    patt = stripped ? stripped.get() : &ref_patt;
  }
  MatrixCacheGuard mol_guard(mol);
  MatrixCacheGuard reference_guard(ref);
  auto match = RDDepict::generateDepictionMatching2DStructure(
      mol, ref, data.conf_id, patt, params);
  rust::Vec<IndexPair> pairs;
  for (const auto &[ref_idx, mol_idx] : match) {
    pairs.push_back(IndexPair{static_cast<std::uint32_t>(ref_idx),
                              static_cast<std::uint32_t>(mol_idx)});
  }
  return pairs;
}

void generate_depiction_matching_2d_structure_atom_map(
    RWMol &mol, const RWMol &reference, rust::Slice<const IndexPair> atom_map,
    const ConstrainedDepictionParamsData &data) {
  // RDKit checks the pairs off by one (an index equal to the atom count
  // passes) and then sets bits by the molecule index unchecked, past the end
  // of a boost::dynamic_bitset.
  RDKit::MatchVectType match;
  for (const auto &pair : atom_map) {
    check_index(pair.first, reference.getNumAtoms());
    check_index(pair.second, mol.getNumAtoms());
    match.emplace_back(static_cast<int>(pair.first),
                       static_cast<int>(pair.second));
  }
  auto params = to_params(data);
  auto lock = lock_depiction();
  check_constrained_inputs(mol, reference, data.conf_id, params);
  MatrixCacheGuard mol_guard(mol);
  MatrixCacheGuard reference_guard(reference);
  RDDepict::generateDepictionMatching2DStructure(mol, reference, match,
                                                 data.conf_id, params);
}

void generate_depiction_matching_3d_structure(
    RWMol &mol, const RWMol &reference, int conf_id, const RWMol &ref_patt,
    bool has_ref_patt, bool accept_failure, bool force_rdkit) {
  auto lock = lock_depiction();
  // It matches the reference to the molecule (control_c_lock.h).
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  std::lock_guard<std::mutex> random_lock(random_mutex());
  // It mimics the reference's distances with RDKit's depictor, or with
  // acceptFailure falls back to compute2DCoords without forceRDKit.
  check_rdkit_depictor_mol(mol);
  if (accept_failure) {
    check_compute_2d_coords_mol(mol, false);
  }
  if (const auto *conf = find_conformer(reference, conf_id)) {
    check_depiction_conformer(*conf, "the reference's coordinates");
  }
  MatrixCacheGuard mol_guard(mol);
  MatrixCacheGuard reference_guard(reference);
  // RDKit takes the pattern as non-const but only matches it.
  auto *pattern = has_ref_patt ? const_cast<RWMol *>(&ref_patt) : nullptr;
  RDDepict::generateDepictionMatching3DStructure(
      mol, reference, conf_id, pattern, accept_failure, force_rdkit);
}

// straightenDepiction and normalizeDepiction read no depiction setting.
void straighten_depiction(RWMol &mol, int conf_id, bool minimize_rotation) {
  // It bins each bond's angle by a cast to int, undefined for NaN.
  checked_bonded_conformer(mol, conf_id);
  RDDepict::straightenDepiction(mol, conf_id, minimize_rotation);
}

double normalize_depiction(RWMol &mol, int conf_id, int canonicalize,
                           double scale_factor) {
  // It bins bond lengths times 10 by a cast to int, undefined beyond 2^31.
  checked_bonded_conformer(mol, conf_id);
  return RDDepict::normalizeDepiction(mol, conf_id, canonicalize, scale_factor);
}

}  // namespace rdkit_cxx
