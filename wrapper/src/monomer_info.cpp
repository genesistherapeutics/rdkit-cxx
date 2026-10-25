/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/monomer_info.h"

#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

namespace {

std::unique_ptr<std::string> str(const std::string &s) {
  return std::make_unique<std::string>(s);
}

AtomMonomerInfo::AtomMonomerType monomer_type(int type) {
  if (type < AtomMonomerInfo::UNKNOWN || type > AtomMonomerInfo::OTHER) {
    throw ValueErrorException("unknown AtomMonomerType " +
                              std::to_string(type));
  }
  return static_cast<AtomMonomerInfo::AtomMonomerType>(type);
}

// RDKit takes an info's type for its class: code that finds PDBRESIDUE
// (MolPickler.cpp, PDBWriter.cpp, AddHs.cpp, ProximityBonds.cpp) casts the
// info to AtomPDBResidueInfo and reads past the end of a plain
// AtomMonomerInfo. Only an AtomPDBResidueInfo may have the type
// (CONVENTIONS.md U10).
void check_not_pdb(AtomMonomerInfo::AtomMonomerType type) {
  if (type == AtomMonomerInfo::PDBRESIDUE) {
    throw ValueErrorException(
        "only AtomPDBResidueInfo can have AtomMonomerType PDBRESIDUE");
  }
}

}  // namespace

std::unique_ptr<AtomMonomerInfo> monomer_info_new() {
  return std::make_unique<AtomMonomerInfo>();
}

std::unique_ptr<AtomMonomerInfo> monomer_info_with(int type, rust::Str name,
                                                   rust::Str residue_name,
                                                   int res_num,
                                                   rust::Str chain_id,
                                                   rust::Str monomer_class) {
  auto typ = monomer_type(type);
  check_not_pdb(typ);
  return std::make_unique<AtomMonomerInfo>(
      typ, std::string(name), std::string(residue_name), res_num,
      std::string(chain_id), std::string(monomer_class));
}

std::unique_ptr<AtomMonomerInfo> monomer_info_clone(
    const AtomMonomerInfo &info) {
  return std::unique_ptr<AtomMonomerInfo>(info.copy());
}

std::unique_ptr<std::string> monomer_info_get_name(
    const AtomMonomerInfo &info) {
  return str(info.getName());
}

int monomer_info_get_monomer_type(const AtomMonomerInfo &info) {
  return info.getMonomerType();
}

std::unique_ptr<std::string> monomer_info_get_residue_name(
    const AtomMonomerInfo &info) {
  return str(info.getResidueName());
}

int monomer_info_get_residue_number(const AtomMonomerInfo &info) {
  return info.getResidueNumber();
}

std::unique_ptr<std::string> monomer_info_get_chain_id(
    const AtomMonomerInfo &info) {
  return str(info.getChainId());
}

std::unique_ptr<std::string> monomer_info_get_monomer_class(
    const AtomMonomerInfo &info) {
  return str(info.getMonomerClass());
}

void monomer_info_set_name(AtomMonomerInfo &info, rust::Str nm) {
  info.setName(std::string(nm));
}

void monomer_info_set_monomer_type(AtomMonomerInfo &info, int typ) {
  auto type = monomer_type(typ);
  // A PDB residue's info, or a plain copy of one, may have its type back.
  if (!dynamic_cast<AtomPDBResidueInfo *>(&info)) {
    check_not_pdb(type);
  }
  info.setMonomerType(type);
}

void monomer_info_set_residue_name(AtomMonomerInfo &info, rust::Str val) {
  info.setResidueName(std::string(val));
}

void monomer_info_set_residue_number(AtomMonomerInfo &info, int val) {
  info.setResidueNumber(val);
}

void monomer_info_set_chain_id(AtomMonomerInfo &info, rust::Str val) {
  info.setChainId(std::string(val));
}

void monomer_info_set_monomer_class(AtomMonomerInfo &info, rust::Str val) {
  info.setMonomerClass(std::string(val));
}

std::unique_ptr<AtomPDBResidueInfo> pdb_residue_info_new() {
  return std::make_unique<AtomPDBResidueInfo>();
}

std::unique_ptr<AtomPDBResidueInfo> pdb_residue_info_with(
    rust::Str atom_name, int serial_number, rust::Str alt_loc,
    rust::Str residue_name, int residue_number, rust::Str chain_id,
    rust::Str insertion_code, double occupancy, double temp_factor,
    bool is_hetero_atom, unsigned secondary_structure, unsigned segment_number,
    rust::Str monomer_class) {
  return std::make_unique<AtomPDBResidueInfo>(
      std::string(atom_name), serial_number, std::string(alt_loc),
      std::string(residue_name), residue_number, std::string(chain_id),
      std::string(insertion_code), occupancy, temp_factor, is_hetero_atom,
      secondary_structure, segment_number, std::string(monomer_class));
}

std::unique_ptr<AtomPDBResidueInfo> pdb_residue_info_clone(
    const AtomPDBResidueInfo &info) {
  return std::make_unique<AtomPDBResidueInfo>(info);
}

const AtomMonomerInfo &pdb_residue_info_as_monomer_info(
    const AtomPDBResidueInfo &info) {
  return info;
}

AtomMonomerInfo &pdb_residue_info_as_monomer_info_mut(
    AtomPDBResidueInfo &info) {
  return info;
}

int pdb_residue_info_get_serial_number(const AtomPDBResidueInfo &info) {
  return info.getSerialNumber();
}

std::unique_ptr<std::string> pdb_residue_info_get_alt_loc(
    const AtomPDBResidueInfo &info) {
  return str(info.getAltLoc());
}

std::unique_ptr<std::string> pdb_residue_info_get_insertion_code(
    const AtomPDBResidueInfo &info) {
  return str(info.getInsertionCode());
}

double pdb_residue_info_get_occupancy(const AtomPDBResidueInfo &info) {
  return info.getOccupancy();
}

double pdb_residue_info_get_temp_factor(const AtomPDBResidueInfo &info) {
  return info.getTempFactor();
}

bool pdb_residue_info_get_is_hetero_atom(const AtomPDBResidueInfo &info) {
  return info.getIsHeteroAtom();
}

unsigned pdb_residue_info_get_secondary_structure(
    const AtomPDBResidueInfo &info) {
  return info.getSecondaryStructure();
}

unsigned pdb_residue_info_get_segment_number(const AtomPDBResidueInfo &info) {
  return info.getSegmentNumber();
}

void pdb_residue_info_set_serial_number(AtomPDBResidueInfo &info, int val) {
  info.setSerialNumber(val);
}

void pdb_residue_info_set_alt_loc(AtomPDBResidueInfo &info, rust::Str val) {
  info.setAltLoc(std::string(val));
}

void pdb_residue_info_set_insertion_code(AtomPDBResidueInfo &info,
                                         rust::Str val) {
  info.setInsertionCode(std::string(val));
}

void pdb_residue_info_set_occupancy(AtomPDBResidueInfo &info, double val) {
  info.setOccupancy(val);
}

void pdb_residue_info_set_temp_factor(AtomPDBResidueInfo &info, double val) {
  info.setTempFactor(val);
}

void pdb_residue_info_set_is_hetero_atom(AtomPDBResidueInfo &info, bool val) {
  info.setIsHeteroAtom(val);
}

void pdb_residue_info_set_secondary_structure(AtomPDBResidueInfo &info,
                                              unsigned val) {
  info.setSecondaryStructure(val);
}

void pdb_residue_info_set_segment_number(AtomPDBResidueInfo &info,
                                         unsigned val) {
  info.setSegmentNumber(val);
}

std::unique_ptr<AtomMonomerInfo> atom_get_monomer_info(const Atom &atom) {
  const auto *info = atom.getMonomerInfo();
  return info ? std::unique_ptr<AtomMonomerInfo>(info->copy()) : nullptr;
}

std::unique_ptr<AtomPDBResidueInfo> atom_get_pdb_residue_info(
    const Atom &atom) {
  const auto *info = atom.getMonomerInfo();
  if (!info) {
    return nullptr;
  }
  // Python casts; only an AtomPDBResidueInfo can have the type (U10), but the
  // class is checked too.
  const auto *pdb = dynamic_cast<const AtomPDBResidueInfo *>(info);
  if (info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE || !pdb) {
    throw ValueErrorException("MonomerInfo is not a PDB Residue");
  }
  return std::make_unique<AtomPDBResidueInfo>(*pdb);
}

void atom_set_monomer_info(Atom &atom, const AtomMonomerInfo &info) {
  atom.setMonomerInfo(info.copy());
}

}  // namespace rdkit_cxx
