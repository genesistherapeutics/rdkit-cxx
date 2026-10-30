/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdForceFieldHelpers: UFF and MMFF force fields for molecules.

#pragma once

#include <cstdint>
#include <memory>

#include <ForceField/ForceField.h>

#include "rdkit-cxx/force_field.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridges (src/ffi/force_field_helpers.rs,
// force_field.rs).
struct ConfOptimization;
struct TermParams;

// Minimizes each of the molecule's conformers with the force field, as
// ForceFieldsHelper::OptimizeMoleculeConfs (ForceFieldHelpers/
// FFConvenience.h), and leaves the force field's positions as they were.
rust::Vec<ConfOptimization> optimize_confs(RDKit::ROMol &mol,
                                           ForceFields::ForceField &field,
                                           int num_threads, unsigned max_iters);

int uff_optimize_molecule(RWMol &mol, unsigned max_iters, double vdw_thresh,
                          int conf_id, bool ignore_interfrag_interactions);
rust::Vec<ConfOptimization> uff_optimize_molecule_confs(
    RWMol &mol, int num_threads, unsigned max_iters, double vdw_thresh,
    bool ignore_interfrag_interactions);
std::unique_ptr<ForceField> uff_get_molecule_force_field(
    RWMol &mol, double vdw_thresh, int conf_id,
    bool ignore_interfrag_interactions);
bool uff_has_all_molecule_params(const RWMol &mol);
int mmff_optimize_molecule(RWMol &mol, rust::Str mmff_variant,
                           unsigned max_iters, double non_bonded_thresh,
                           int conf_id, bool ignore_interfrag_interactions);
unsigned mmff_sanitize_molecule(RWMol &mol);
// Null when MMFF cannot type every atom.
std::unique_ptr<MMFFMolProperties> mmff_get_molecule_properties(
    RWMol &mol, rust::Str mmff_variant, MMFFVerbosity verbosity);
std::unique_ptr<ForceField> mmff_get_molecule_force_field(
    RWMol &mol, const MMFFMolProperties &props, double non_bonded_thresh,
    int conf_id, bool ignore_interfrag_interactions);
std::unique_ptr<ForceField> create_empty_force_field_for_mol(RWMol &mol,
                                                             int conf_id);
bool mmff_has_all_molecule_params(const RWMol &mol);
rust::Vec<ConfOptimization> mmff_optimize_molecule_confs(
    RWMol &mol, int num_threads, unsigned max_iters, rust::Str mmff_variant,
    double non_bonded_thresh, bool ignore_interfrag_interactions);
TermParams get_uff_bond_stretch_params(const RWMol &mol, unsigned idx1,
                                       unsigned idx2);
TermParams get_uff_angle_bend_params(const RWMol &mol, unsigned idx1,
                                     unsigned idx2, unsigned idx3);
TermParams get_uff_torsion_params(const RWMol &mol, unsigned idx1,
                                  unsigned idx2, unsigned idx3, unsigned idx4);
TermParams get_uff_inversion_params(const RWMol &mol, unsigned idx1,
                                    unsigned idx2, unsigned idx3,
                                    unsigned idx4);
TermParams get_uff_vdw_params(const RWMol &mol, unsigned idx1, unsigned idx2);
int optimize_molecule(ForceField &ff, unsigned max_iters);
rust::Vec<ConfOptimization> optimize_molecule_confs(RWMol &mol, ForceField &ff,
                                                    int num_threads,
                                                    unsigned max_iters);

}  // namespace rdkit_cxx
