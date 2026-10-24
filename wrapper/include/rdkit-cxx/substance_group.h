/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Substance groups, which live in their molecule; the Rust views reach them by
// molecule and index (src/chem/rdchem/substance_group.rs).

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/SubstanceGroup.h>

#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using SubstanceGroup = RDKit::SubstanceGroup;

// Shared structs, defined by the bridges (src/ffi/geometry.rs,
// src/ffi/substance_group.rs).
struct Point3D;
struct SubstanceGroupCStateData;
struct SubstanceGroupAttachData;

std::uint32_t mol_get_num_substance_groups(const RWMol &mol);
const SubstanceGroup &mol_get_substance_group(const RWMol &mol,
                                              std::uint32_t idx);
SubstanceGroup &mol_get_substance_group_mut(RWMol &mol, std::uint32_t idx);
void mol_clear_substance_groups(RWMol &mol);
// As createMolSubstanceGroup, createMolDataSubstanceGroup, and
// addMolSubstanceGroup in Code/GraphMol/nbWrap/SubstanceGroup.cpp; each
// returns the new group's index. The copy's indices are checked against the
// molecule it joins (CONVENTIONS.md E3).
std::uint32_t mol_create_substance_group(RWMol &mol, rust::Str type);
std::uint32_t mol_create_data_substance_group(RWMol &mol, rust::Str field_name,
                                              rust::Str value);
std::uint32_t mol_add_substance_group(RWMol &mol, const SubstanceGroup &sgroup);

std::uint32_t substance_group_get_index_in_mol(const SubstanceGroup &sgroup);
rust::Vec<std::uint32_t> substance_group_get_atoms(
    const SubstanceGroup &sgroup);
rust::Vec<std::uint32_t> substance_group_get_parent_atoms(
    const SubstanceGroup &sgroup);
rust::Vec<std::uint32_t> substance_group_get_bonds(
    const SubstanceGroup &sgroup);
rust::Vec<SubstanceGroupCStateData> substance_group_get_cstates(
    const SubstanceGroup &sgroup);
rust::Vec<SubstanceGroupAttachData> substance_group_get_attach_points(
    const SubstanceGroup &sgroup);
// Three points per bracket; two-point brackets end with the origin.
rust::Vec<Point3D> substance_group_get_brackets(const SubstanceGroup &sgroup);

void substance_group_set_atoms(SubstanceGroup &sgroup,
                               rust::Slice<const std::uint32_t> atoms);
void substance_group_set_parent_atoms(SubstanceGroup &sgroup,
                                      rust::Slice<const std::uint32_t> patoms);
void substance_group_set_bonds(SubstanceGroup &sgroup,
                               rust::Slice<const std::uint32_t> bonds);
void substance_group_add_atom_with_idx(SubstanceGroup &sgroup,
                                       std::uint32_t idx);
void substance_group_add_parent_atom_with_idx(SubstanceGroup &sgroup,
                                              std::uint32_t idx);
void substance_group_add_bond_with_idx(SubstanceGroup &sgroup,
                                       std::uint32_t idx);
void substance_group_add_atom_with_bookmark(SubstanceGroup &sgroup, int mark);
void substance_group_add_parent_atom_with_bookmark(SubstanceGroup &sgroup,
                                                   int mark);
void substance_group_add_bond_with_bookmark(SubstanceGroup &sgroup, int mark);
void substance_group_add_cstate(SubstanceGroup &sgroup, std::uint32_t bond_idx,
                                const Point3D &vector);
void substance_group_add_attach_point(SubstanceGroup &sgroup,
                                      std::uint32_t a_idx, int lv_idx,
                                      rust::Str id_str);
// As addBracketHelper in Code/GraphMol/nbWrap/SubstanceGroup.cpp.
void substance_group_add_bracket(SubstanceGroup &sgroup,
                                 rust::Slice<const Point3D> pts);
void substance_group_clear_brackets(SubstanceGroup &sgroup);
void substance_group_clear_cstates(SubstanceGroup &sgroup);
void substance_group_clear_attach_points(SubstanceGroup &sgroup);

SubstanceGroupCStateData substance_group_cstate_default();
SubstanceGroupAttachData substance_group_attach_default();

}  // namespace rdkit_cxx
