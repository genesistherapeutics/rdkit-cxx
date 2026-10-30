/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/stereo.h"

#include <vector>

#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/stereo.rs.h"
#include "rdkit-cxx/mol_checks.h"

namespace rdkit_cxx {

namespace {

StereoGroupData to_data(const RDKit::StereoGroup &group) {
  StereoGroupData data;
  data.group_type = group.getGroupType();
  for (const auto *atom : group.getAtoms()) {
    data.atoms.push_back(atom->getIdx());
  }
  for (const auto *bond : group.getBonds()) {
    data.bonds.push_back(bond->getIdx());
  }
  data.read_id = group.getReadId();
  data.write_id = group.getWriteId();
  return data;
}

RDKit::StereoGroup from_data(RWMol &mol, StereoGroupType group_type,
                             rust::Slice<const std::uint32_t> atom_ids,
                             rust::Slice<const std::uint32_t> bond_ids,
                             unsigned read_id) {
  std::vector<RDKit::Atom *> atoms;
  for (auto idx : atom_ids) {
    check_index(idx, mol.getNumAtoms());
    atoms.push_back(mol.getAtomWithIdx(idx));
  }
  std::vector<RDKit::Bond *> bonds;
  for (auto idx : bond_ids) {
    check_index(idx, mol.getNumBonds());
    bonds.push_back(mol.getBondWithIdx(idx));
  }
  if (atoms.empty() && bonds.empty()) {
    throw ValueErrorException(
        "New StereoGroup must contain at least one atom or bond.");
  }
  // RDKit 2025.03 removes an atom or bond from a group once (RWMol's
  // removeAtom and removeBond), leaving a duplicate to dangle; 2026.03's
  // constructor refuses them, with these messages.
  check_stereo_group_members(atoms, bonds);
  return RDKit::StereoGroup(group_type, atoms, bonds, read_id);
}

rust::Slice<const std::uint32_t> slice(const rust::Vec<std::uint32_t> &v) {
  return rust::Slice<const std::uint32_t>(v.data(), v.size());
}

}  // namespace

rust::Vec<StereoGroupData> mol_get_stereo_groups(const RWMol &mol) {
  rust::Vec<StereoGroupData> groups;
  for (const auto &group : mol.getStereoGroups()) {
    groups.push_back(to_data(group));
  }
  return groups;
}

void rwmol_set_stereo_groups(RWMol &mol,
                             rust::Slice<const StereoGroupData> groups) {
  std::vector<RDKit::StereoGroup> result;
  for (const auto &data : groups) {
    auto group = from_data(mol, data.group_type, slice(data.atoms),
                           slice(data.bonds), data.read_id);
    group.setWriteId(data.write_id);
    result.push_back(std::move(group));
  }
  mol.setStereoGroups(std::move(result));
}

StereoGroupData create_stereo_group(StereoGroupType group_type,
                                    const RWMol &mol,
                                    rust::Slice<const std::uint32_t> atom_ids,
                                    rust::Slice<const std::uint32_t> bond_ids,
                                    unsigned read_id) {
  // The group is built only to read back RDKit's checks and ids.
  return to_data(from_data(const_cast<RWMol &>(mol), group_type, atom_ids,
                           bond_ids, read_id));
}

void forward_stereo_group_ids(RWMol &mol) { RDKit::forwardStereoGroupIds(mol); }

rust::Vec<StereoInfoData> find_potential_stereo(RWMol &mol, bool clean_it,
                                                bool flag_possible) {
  rust::Vec<StereoInfoData> result;
  for (const auto &info :
       RDKit::Chirality::findPotentialStereo(mol, clean_it, flag_possible)) {
    StereoInfoData data;
    data.type_ = info.type;
    data.specified = info.specified;
    data.centered_on = info.centeredOn;
    data.descriptor = info.descriptor;
    data.permutation = info.permutation;
    for (auto idx : info.controllingAtoms) {
      data.controlling_atoms.push_back(idx);
    }
    result.push_back(std::move(data));
  }
  return result;
}

}  // namespace rdkit_cxx
