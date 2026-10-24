/*
 * Copyright (c) 2026 Steven Kearnes. All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/substance_group.h"

#include <stdexcept>
#include <vector>

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"
#include "rdkit-cxx/ffi/substance_group.rs.h"

namespace rdkit_cxx {

namespace {

Point3D to_point(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

RDGeom::Point3D from_point(const Point3D &p) {
  return RDGeom::Point3D(p.x, p.y, p.z);
}

rust::Vec<std::uint32_t> to_vec(const std::vector<unsigned> &v) {
  rust::Vec<std::uint32_t> result;
  for (auto idx : v) {
    result.push_back(idx);
  }
  return result;
}

std::vector<unsigned> checked(rust::Slice<const std::uint32_t> indices,
                              unsigned len) {
  std::vector<unsigned> result;
  for (auto idx : indices) {
    check_index(idx, len);
    result.push_back(idx);
  }
  return result;
}

unsigned num_atoms(const SubstanceGroup &sgroup) {
  return sgroup.getOwningMol().getNumAtoms();
}

unsigned num_bonds(const SubstanceGroup &sgroup) {
  return sgroup.getOwningMol().getNumBonds();
}

void check_lv_idx(int lv_idx, unsigned len) {
  // -1 marks an attachment with no leaving atom.
  if (lv_idx != -1) {
    check_index(static_cast<unsigned>(lv_idx), len);
  }
}

// ROMol::getUniqueAtomWithBookmark reads the first entry without checking
// that there is one.
void check_atom_bookmark(RDKit::ROMol &mol, int mark) {
  if (!mol.hasAtomBookmark(mark) || mol.getAllAtomsWithBookmark(mark).empty()) {
    throw std::runtime_error("bookmark not found");
  }
}

void check_bond_bookmark(RDKit::ROMol &mol, int mark) {
  if (!mol.hasBondBookmark(mark) || mol.getAllBondsWithBookmark(mark).empty()) {
    throw std::runtime_error("bookmark not found");
  }
}

}  // namespace

std::uint32_t mol_get_num_substance_groups(const RWMol &mol) {
  return RDKit::getSubstanceGroups(mol).size();
}

const SubstanceGroup &mol_get_substance_group(const RWMol &mol,
                                              std::uint32_t idx) {
  const auto &sgroups = RDKit::getSubstanceGroups(mol);
  check_index(idx, sgroups.size());
  return sgroups[idx];
}

SubstanceGroup &mol_get_substance_group_mut(RWMol &mol, std::uint32_t idx) {
  auto &sgroups = RDKit::getSubstanceGroups(mol);
  check_index(idx, sgroups.size());
  return sgroups[idx];
}

void mol_clear_substance_groups(RWMol &mol) {
  RDKit::getSubstanceGroups(mol).clear();
}

std::uint32_t mol_create_substance_group(RWMol &mol, rust::Str type) {
  return RDKit::addSubstanceGroup(mol, SubstanceGroup(&mol, std::string(type)));
}

std::uint32_t mol_create_data_substance_group(RWMol &mol, rust::Str field_name,
                                              rust::Str value) {
  SubstanceGroup sgroup(&mol, "DAT");
  sgroup.setProp("FIELDNAME", std::string(field_name));
  RDKit::STR_VECT data_fields{std::string(value)};
  sgroup.setProp("DATAFIELDS", data_fields);
  return RDKit::addSubstanceGroup(mol, sgroup);
}

std::uint32_t mol_add_substance_group(RWMol &mol,
                                      const SubstanceGroup &sgroup) {
  auto atoms = mol.getNumAtoms();
  auto bonds = mol.getNumBonds();
  for (const auto *list : {&sgroup.getAtoms(), &sgroup.getParentAtoms()}) {
    for (auto idx : *list) {
      check_index(idx, atoms);
    }
  }
  for (auto idx : sgroup.getBonds()) {
    check_index(idx, bonds);
  }
  for (const auto &cstate : sgroup.getCStates()) {
    check_index(cstate.bondIdx, bonds);
  }
  for (const auto &ap : sgroup.getAttachPoints()) {
    check_index(ap.aIdx, atoms);
    check_lv_idx(ap.lvIdx, atoms);
  }
  return RDKit::addSubstanceGroup(mol, sgroup);
}

std::uint32_t substance_group_get_index_in_mol(const SubstanceGroup &sgroup) {
  return sgroup.getIndexInMol();
}

rust::Vec<std::uint32_t> substance_group_get_atoms(
    const SubstanceGroup &sgroup) {
  return to_vec(sgroup.getAtoms());
}

rust::Vec<std::uint32_t> substance_group_get_parent_atoms(
    const SubstanceGroup &sgroup) {
  return to_vec(sgroup.getParentAtoms());
}

rust::Vec<std::uint32_t> substance_group_get_bonds(
    const SubstanceGroup &sgroup) {
  return to_vec(sgroup.getBonds());
}

rust::Vec<SubstanceGroupCStateData> substance_group_get_cstates(
    const SubstanceGroup &sgroup) {
  rust::Vec<SubstanceGroupCStateData> result;
  for (const auto &cstate : sgroup.getCStates()) {
    result.push_back(
        SubstanceGroupCStateData{cstate.bondIdx, to_point(cstate.vector)});
  }
  return result;
}

rust::Vec<SubstanceGroupAttachData> substance_group_get_attach_points(
    const SubstanceGroup &sgroup) {
  rust::Vec<SubstanceGroupAttachData> result;
  for (const auto &ap : sgroup.getAttachPoints()) {
    result.push_back(SubstanceGroupAttachData{
        ap.aIdx, ap.lvIdx, std::make_unique<std::string>(ap.id)});
  }
  return result;
}

rust::Vec<Point3D> substance_group_get_brackets(const SubstanceGroup &sgroup) {
  rust::Vec<Point3D> result;
  for (const auto &bracket : sgroup.getBrackets()) {
    for (const auto &p : bracket) {
      result.push_back(to_point(p));
    }
  }
  return result;
}

void substance_group_set_atoms(SubstanceGroup &sgroup,
                               rust::Slice<const std::uint32_t> atoms) {
  sgroup.setAtoms(checked(atoms, num_atoms(sgroup)));
}

void substance_group_set_parent_atoms(SubstanceGroup &sgroup,
                                      rust::Slice<const std::uint32_t> patoms) {
  sgroup.setParentAtoms(checked(patoms, num_atoms(sgroup)));
}

void substance_group_set_bonds(SubstanceGroup &sgroup,
                               rust::Slice<const std::uint32_t> bonds) {
  sgroup.setBonds(checked(bonds, num_bonds(sgroup)));
}

void substance_group_add_atom_with_idx(SubstanceGroup &sgroup,
                                       std::uint32_t idx) {
  check_index(idx, num_atoms(sgroup));
  sgroup.addAtomWithIdx(idx);
}

void substance_group_add_parent_atom_with_idx(SubstanceGroup &sgroup,
                                              std::uint32_t idx) {
  check_index(idx, num_atoms(sgroup));
  sgroup.addParentAtomWithIdx(idx);
}

void substance_group_add_bond_with_idx(SubstanceGroup &sgroup,
                                       std::uint32_t idx) {
  check_index(idx, num_bonds(sgroup));
  sgroup.addBondWithIdx(idx);
}

void substance_group_add_atom_with_bookmark(SubstanceGroup &sgroup, int mark) {
  check_atom_bookmark(sgroup.getOwningMol(), mark);
  sgroup.addAtomWithBookmark(mark);
}

void substance_group_add_parent_atom_with_bookmark(SubstanceGroup &sgroup,
                                                   int mark) {
  check_atom_bookmark(sgroup.getOwningMol(), mark);
  sgroup.addParentAtomWithBookmark(mark);
}

void substance_group_add_bond_with_bookmark(SubstanceGroup &sgroup, int mark) {
  check_bond_bookmark(sgroup.getOwningMol(), mark);
  sgroup.addBondWithBookmark(mark);
}

void substance_group_add_cstate(SubstanceGroup &sgroup, std::uint32_t bond_idx,
                                const Point3D &vector) {
  check_index(bond_idx, num_bonds(sgroup));
  sgroup.addCState(bond_idx, from_point(vector));
}

void substance_group_add_attach_point(SubstanceGroup &sgroup,
                                      std::uint32_t a_idx, int lv_idx,
                                      rust::Str id_str) {
  check_index(a_idx, num_atoms(sgroup));
  check_lv_idx(lv_idx, num_atoms(sgroup));
  sgroup.addAttachPoint(a_idx, lv_idx, std::string(id_str));
}

void substance_group_add_bracket(SubstanceGroup &sgroup,
                                 rust::Slice<const Point3D> pts) {
  if (pts.size() != 2 && pts.size() != 3) {
    throw ValueErrorException("pts object have a length of 2 or 3");
  }
  SubstanceGroup::Bracket bracket;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    bracket[i] = from_point(pts[i]);
  }
  sgroup.addBracket(bracket);
}

void substance_group_clear_brackets(SubstanceGroup &sgroup) {
  sgroup.clearBrackets();
}

void substance_group_clear_cstates(SubstanceGroup &sgroup) {
  sgroup.clearCStates();
}

void substance_group_clear_attach_points(SubstanceGroup &sgroup) {
  sgroup.clearAttachPoints();
}

SubstanceGroupCStateData substance_group_cstate_default() {
  SubstanceGroup::CState cstate{};
  return SubstanceGroupCStateData{cstate.bondIdx, to_point(cstate.vector)};
}

SubstanceGroupAttachData substance_group_attach_default() {
  SubstanceGroup::AttachPoint ap{};
  return SubstanceGroupAttachData{ap.aIdx, ap.lvIdx,
                                  std::make_unique<std::string>(ap.id)};
}

}  // namespace rdkit_cxx
