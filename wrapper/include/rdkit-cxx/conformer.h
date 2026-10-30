/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#pragma once

#include <memory>

#include <GraphMol/Conformer.h>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using Conformer = RDKit::Conformer;

std::unique_ptr<Conformer> conformer_new();
std::unique_ptr<Conformer> conformer_with_num_atoms(unsigned num_atoms);
std::unique_ptr<Conformer> conformer_clone(const Conformer &conf);

unsigned conformer_get_id(const Conformer &conf);
unsigned conformer_get_num_atoms(const Conformer &conf);
bool conformer_is_3d(const Conformer &conf);
bool conformer_has_owning_mol(const Conformer &conf);
// The position functions check their atom index (E3).
Point3D conformer_get_atom_position(const Conformer &conf, unsigned aid);
rust::Vec<Point3D> conformer_get_positions(const Conformer &conf);
void conformer_set_atom_position(Conformer &conf, unsigned aid,
                                 const Point3D &loc);
void conformer_set_positions(Conformer &conf,
                             rust::Slice<const Point3D> positions);
void conformer_set_3d(Conformer &conf, bool v);
void conformer_set_id(Conformer &conf, unsigned id);

unsigned mol_get_num_conformers(const RWMol &mol);
const Conformer &mol_get_conformer(const RWMol &mol, int id);
const Conformer &mol_find_conformer(const RWMol &mol, unsigned id);
Conformer &mol_find_conformer_mut(RWMol &mol, unsigned id);
void mol_set_conformer_id(RWMol &mol, unsigned current, unsigned id);
// The conformer at a position in the molecule's list, for iteration.
const Conformer &mol_get_conformer_at(const RWMol &mol, unsigned position);
unsigned mol_add_conformer(RWMol &mol, const Conformer &conf, bool assign_id);
void mol_remove_conformer(RWMol &mol, unsigned id);
void mol_remove_all_conformers(RWMol &mol);

}  // namespace rdkit_cxx
