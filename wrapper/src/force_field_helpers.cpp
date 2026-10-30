// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/force_field_helpers.h"

#include <algorithm>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <ForceField/UFF/Params.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/force_field.rs.h"
#include "rdkit-cxx/ffi/force_field_helpers.rs.h"
#include "rdkit-cxx/mol_checks.h"
#include "rdkit-cxx/threads.h"

namespace rdkit_cxx {

namespace {

// The builders size their neighbor matrices in unsigned int, and the MMFF
// one reads the int distance matrix (U3, U6).
void check_builder_size(const RDKit::ROMol &mol) { check_matrix_size(mol); }

// The builders match torsion patterns (addTorsions in UFF/Builder.cpp and
// MMFF/Builder.cpp), so they run under control_c_mutex (control_c_lock.h);
// minimizing, which may use several threads, matches nothing.
// Ignoring interfragment interactions, UFF's builder splits the molecule
// into its fragments (getMolFrags), removing atoms from copies, which trusts
// what check_removals checks.
template <class... Args>
ForceFields::ForceField *build_uff(RDKit::ROMol &mol, Args &&...args) {
  check_removals(mol, Removes::Fragments);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return RDKit::UFF::constructForceField(mol, std::forward<Args>(args)...);
}

template <class... Args>
ForceFields::ForceField *build_mmff(Args &&...args) {
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return RDKit::MMFF::constructForceField(std::forward<Args>(args)...);
}

// Minimizes one conformer, as OptimizeMoleculeConfsST.
ConfOptimization optimize_conf(ForceFields::ForceField &field,
                               RDKit::Conformer &conf, unsigned num_atoms,
                               unsigned max_iters) {
  for (unsigned aidx = 0; aidx < num_atoms; ++aidx) {
    field.positions()[aidx] = &conf.getAtomPos(aidx);
  }
  field.initialize();
  int needs_more = field.minimize(max_iters);
  double energy = field.calcEnergy();
  return ConfOptimization{needs_more, energy};
}

std::unique_ptr<ForceField> wrap(ForceFields::ForceField *field) {
  auto ff = std::make_unique<ForceField>(
      std::unique_ptr<ForceFields::ForceField>(field));
  ff->field().initialize();
  return ff;
}

TermParams not_found() { return TermParams{false, 0, {}}; }

TermParams found(std::initializer_list<double> values) {
  rust::Vec<double> result;
  for (auto value : values) {
    result.push_back(value);
  }
  return TermParams{true, 0, std::move(result)};
}

void check_indices(const RWMol &mol, std::initializer_list<unsigned> indices) {
  for (auto idx : indices) {
    check_index(idx, mol.getNumAtoms());
  }
}

}  // namespace

// RDKit's single-threaded loop writes the first numAtoms positions, past the
// end of a force field with fewer points, and leaves the force field pointing
// at the last conformer. Its threads copy the force field, which keeps its
// contributions but not its points or fixed points (ForceField's copy
// constructor), resize the copy to numAtoms points, and call std::terminate
// when minimizing throws (a NaN coordinate fails an invariant in
// linearSearch). Here one thread runs RDKit's loop on the force field and
// then restores its points; more threads run in run_on_threads, each on a
// copy made as RDKit makes it.
rust::Vec<ConfOptimization> optimize_confs(RDKit::ROMol &mol,
                                           ForceFields::ForceField &field,
                                           int num_threads,
                                           unsigned max_iters) {
  unsigned num_atoms = mol.getNumAtoms();
  std::vector<RDKit::Conformer *> confs;
  for (auto it = mol.beginConformers(); it != mol.endConformers(); ++it) {
    confs.push_back(it->get());
  }
  std::vector<ConfOptimization> results(confs.size());
  unsigned threads = num_threads_to_use(num_threads);
  if (threads == 1) {
    if (field.positions().size() < num_atoms) {
      throw ValueErrorException(
          "the force field has fewer points than the molecule has atoms");
    }
    check_minimize_size(field.positions().size(), field.dimension());
    auto saved = field.positions();
    auto restore = [&]() {
      field.positions() = saved;
      field.initialize();
    };
    try {
      for (std::size_t i = 0; i < confs.size(); ++i) {
        results[i] = optimize_conf(field, *confs[i], num_atoms, max_iters);
      }
    } catch (...) {
      restore();
      throw;
    }
    restore();
  } else {
    // A contribution may use a point past numAtoms, which a copy lacks.
    if (field.positions().size() != num_atoms) {
      throw ValueErrorException(
          "with more than one thread, RDKit's copies of the force field keep "
          "only one point per atom, and this one has others");
    }
    check_minimize_size(num_atoms, field.dimension());
    // Conformers are independent, so threads past their number have no
    // work.
    unsigned n = static_cast<unsigned>(
        std::min<std::size_t>(threads, std::max<std::size_t>(confs.size(), 1)));
    run_on_threads(n, [&](unsigned t) {
      ForceFields::ForceField copy(field);
      copy.positions().resize(num_atoms);
      for (std::size_t i = t; i < confs.size(); i += n) {
        results[i] = optimize_conf(copy, *confs[i], num_atoms, max_iters);
      }
    });
  }
  rust::Vec<ConfOptimization> res;
  for (const auto &result : results) {
    res.push_back(result);
  }
  return res;
}

// As UFFHelper in Code/GraphMol/ForceFieldHelpers/Wrap/rdForceFields.cpp.
int uff_optimize_molecule(RWMol &mol, unsigned max_iters, double vdw_thresh,
                          int conf_id, bool ignore_interfrag_interactions) {
  check_builder_size(mol);
  check_minimize_size(mol.getNumAtoms(), 3);
  std::unique_ptr<ForceFields::ForceField> ff(
      build_uff(mol, vdw_thresh, conf_id, ignore_interfrag_interactions));
  return RDKit::ForceFieldsHelper::OptimizeMolecule(*ff, max_iters).first;
}

rust::Vec<ConfOptimization> uff_optimize_molecule_confs(
    RWMol &mol, int num_threads, unsigned max_iters, double vdw_thresh,
    bool ignore_interfrag_interactions) {
  check_builder_size(mol);
  std::unique_ptr<ForceFields::ForceField> ff(
      build_uff(mol, vdw_thresh, -1, ignore_interfrag_interactions));
  return optimize_confs(mol, *ff, num_threads, max_iters);
}

std::unique_ptr<ForceField> uff_get_molecule_force_field(
    RWMol &mol, double vdw_thresh, int conf_id,
    bool ignore_interfrag_interactions) {
  check_builder_size(mol);
  return wrap(
      build_uff(mol, vdw_thresh, conf_id, ignore_interfrag_interactions));
}

bool uff_has_all_molecule_params(const RWMol &mol) {
  return RDKit::UFF::getAtomTypes(mol).second;
}

// As MMFFOptimizeMolecule in rdForceFields.cpp. The builder reads the
// distance matrix (U6).
int mmff_optimize_molecule(RWMol &mol, rust::Str mmff_variant,
                           unsigned max_iters, double non_bonded_thresh,
                           int conf_id, bool ignore_interfrag_interactions) {
  check_builder_size(mol);
  check_minimize_size(mol.getNumAtoms(), 3);
  MatrixCacheGuard guard(mol);
  RDKit::MMFF::MMFFMolProperties props(mol, std::string(mmff_variant));
  if (!props.isValid()) {
    return -1;
  }
  std::unique_ptr<ForceFields::ForceField> ff(build_mmff(
      mol, &props, non_bonded_thresh, conf_id, ignore_interfrag_interactions));
  ff->initialize();
  return ff->minimize(max_iters);
}

unsigned mmff_sanitize_molecule(RWMol &mol) {
  return RDKit::MMFF::sanitizeMMFFMol(mol);
}

std::unique_ptr<MMFFMolProperties> mmff_get_molecule_properties(
    RWMol &mol, rust::Str mmff_variant, MMFFVerbosity verbosity) {
  auto props = std::make_unique<MMFFMolProperties>(
      mol, std::string(mmff_variant), static_cast<std::uint8_t>(verbosity));
  if (!props->props().isValid()) {
    return nullptr;
  }
  return props;
}

// The builder looks up each atom's type by its index; properties for a
// smaller molecule fail RDKit's range check, and a larger one's types belong
// to other atoms, so the sizes must match. It takes the properties by
// non-const pointer only to reach their output stream.
std::unique_ptr<ForceField> mmff_get_molecule_force_field(
    RWMol &mol, const MMFFMolProperties &props, double non_bonded_thresh,
    int conf_id, bool ignore_interfrag_interactions) {
  if (props.num_atoms() != mol.getNumAtoms()) {
    throw ValueErrorException(
        "the MMFF properties are for a molecule with another number of atoms");
  }
  check_builder_size(mol);
  MatrixCacheGuard guard(mol);
  auto &mmff_props =
      const_cast<RDKit::MMFF::MMFFMolProperties &>(props.props());
  return wrap(build_mmff(mol, &mmff_props, non_bonded_thresh, conf_id,
                         ignore_interfrag_interactions));
}

std::unique_ptr<ForceField> create_empty_force_field_for_mol(RWMol &mol,
                                                             int conf_id) {
  auto ff = std::make_unique<ForceField>(
      RDKit::ForceFieldsHelper::createEmptyForceFieldForMol(mol, conf_id));
  ff->field().initialize();
  return ff;
}

// As MMFFHasAllMoleculeParams in rdForceFields.cpp, which types a copy.
bool mmff_has_all_molecule_params(const RWMol &mol) {
  RDKit::RWMol copy(mol);
  RDKit::MMFF::MMFFMolProperties props(copy);
  return props.isValid();
}

// As MMFF::MMFFOptimizeMoleculeConfs (ForceFieldHelpers/MMFF/MMFF.h), whose
// result for untypeable atoms is (-1, -1) for each conformer.
rust::Vec<ConfOptimization> mmff_optimize_molecule_confs(
    RWMol &mol, int num_threads, unsigned max_iters, rust::Str mmff_variant,
    double non_bonded_thresh, bool ignore_interfrag_interactions) {
  check_builder_size(mol);
  MatrixCacheGuard guard(mol);
  RDKit::MMFF::MMFFMolProperties props(mol, std::string(mmff_variant));
  if (!props.isValid()) {
    rust::Vec<ConfOptimization> res;
    for (unsigned i = 0; i < mol.getNumConformers(); ++i) {
      res.push_back(ConfOptimization{-1, -1.0});
    }
    return res;
  }
  std::unique_ptr<ForceFields::ForceField> ff(build_mmff(
      mol, &props, non_bonded_thresh, -1, ignore_interfrag_interactions));
  return optimize_confs(mol, *ff, num_threads, max_iters);
}

// The getters read atoms through range-checked accessors; the wrapper's
// check comes first (E3).
TermParams get_uff_bond_stretch_params(const RWMol &mol, unsigned idx1,
                                       unsigned idx2) {
  check_indices(mol, {idx1, idx2});
  ForceFields::UFF::UFFBond params;
  if (!RDKit::UFF::getUFFBondStretchParams(mol, idx1, idx2, params)) {
    return not_found();
  }
  return found({params.kb, params.r0});
}

TermParams get_uff_angle_bend_params(const RWMol &mol, unsigned idx1,
                                     unsigned idx2, unsigned idx3) {
  check_indices(mol, {idx1, idx2, idx3});
  ForceFields::UFF::UFFAngle params;
  if (!RDKit::UFF::getUFFAngleBendParams(mol, idx1, idx2, idx3, params)) {
    return not_found();
  }
  return found({params.ka, params.theta0});
}

TermParams get_uff_torsion_params(const RWMol &mol, unsigned idx1,
                                  unsigned idx2, unsigned idx3, unsigned idx4) {
  check_indices(mol, {idx1, idx2, idx3, idx4});
  ForceFields::UFF::UFFTor params;
  if (!RDKit::UFF::getUFFTorsionParams(mol, idx1, idx2, idx3, idx4, params)) {
    return not_found();
  }
  return found({params.V});
}

TermParams get_uff_inversion_params(const RWMol &mol, unsigned idx1,
                                    unsigned idx2, unsigned idx3,
                                    unsigned idx4) {
  check_indices(mol, {idx1, idx2, idx3, idx4});
  ForceFields::UFF::UFFInv params;
  if (!RDKit::UFF::getUFFInversionParams(mol, idx1, idx2, idx3, idx4, params)) {
    return not_found();
  }
  return found({params.K});
}

TermParams get_uff_vdw_params(const RWMol &mol, unsigned idx1, unsigned idx2) {
  check_indices(mol, {idx1, idx2});
  ForceFields::UFF::UFFVdW params;
  if (!RDKit::UFF::getUFFVdWParams(mol, idx1, idx2, params)) {
    return not_found();
  }
  return found({params.x_ij, params.D_ij});
}

// As FFHelper in rdForceFields.cpp.
int optimize_molecule(ForceField &ff, unsigned max_iters) {
  check_minimize_size(ff.field().positions().size(), ff.field().dimension());
  return RDKit::ForceFieldsHelper::OptimizeMolecule(ff.field(), max_iters)
      .first;
}

rust::Vec<ConfOptimization> optimize_molecule_confs(RWMol &mol, ForceField &ff,
                                                    int num_threads,
                                                    unsigned max_iters) {
  return optimize_confs(mol, ff.field(), num_threads, max_iters);
}

}  // namespace rdkit_cxx
