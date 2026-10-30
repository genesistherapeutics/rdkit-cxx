/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Atom monomer and PDB residue information. Rust holds copies; the atom
// accessors copy in and out.

#pragma once

#include <memory>
#include <string>

#include <GraphMol/MonomerInfo.h>

#include "rdkit-cxx/atom.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using AtomMonomerInfo = RDKit::AtomMonomerInfo;
using AtomPDBResidueInfo = RDKit::AtomPDBResidueInfo;

// AtomMonomerType has no fixed underlying type, so it crosses as int. Only
// an AtomPDBResidueInfo may have type PDBRESIDUE (CONVENTIONS.md U10), so
// monomer_info_with refuses it, and monomer_info_set_monomer_type refuses it
// for an info of another class, with ValueError.
std::unique_ptr<AtomMonomerInfo> monomer_info_new();
std::unique_ptr<AtomMonomerInfo> monomer_info_with(int type, rust::Str name,
                                                   rust::Str residue_name,
                                                   int res_num,
                                                   rust::Str chain_id,
                                                   rust::Str monomer_class);
std::unique_ptr<AtomMonomerInfo> monomer_info_clone(
    const AtomMonomerInfo &info);
std::unique_ptr<std::string> monomer_info_get_name(const AtomMonomerInfo &info);
int monomer_info_get_monomer_type(const AtomMonomerInfo &info);
std::unique_ptr<std::string> monomer_info_get_residue_name(
    const AtomMonomerInfo &info);
int monomer_info_get_residue_number(const AtomMonomerInfo &info);
std::unique_ptr<std::string> monomer_info_get_chain_id(
    const AtomMonomerInfo &info);
std::unique_ptr<std::string> monomer_info_get_monomer_class(
    const AtomMonomerInfo &info);
void monomer_info_set_name(AtomMonomerInfo &info, rust::Str nm);
void monomer_info_set_monomer_type(AtomMonomerInfo &info, int typ);
void monomer_info_set_residue_name(AtomMonomerInfo &info, rust::Str val);
void monomer_info_set_residue_number(AtomMonomerInfo &info, int val);
void monomer_info_set_chain_id(AtomMonomerInfo &info, rust::Str val);
void monomer_info_set_monomer_class(AtomMonomerInfo &info, rust::Str val);

std::unique_ptr<AtomPDBResidueInfo> pdb_residue_info_new();
std::unique_ptr<AtomPDBResidueInfo> pdb_residue_info_with(
    rust::Str atom_name, int serial_number, rust::Str alt_loc,
    rust::Str residue_name, int residue_number, rust::Str chain_id,
    rust::Str insertion_code, double occupancy, double temp_factor,
    bool is_hetero_atom, unsigned secondary_structure, unsigned segment_number,
    rust::Str monomer_class);
std::unique_ptr<AtomPDBResidueInfo> pdb_residue_info_clone(
    const AtomPDBResidueInfo &info);
const AtomMonomerInfo &pdb_residue_info_as_monomer_info(
    const AtomPDBResidueInfo &info);
AtomMonomerInfo &pdb_residue_info_as_monomer_info_mut(AtomPDBResidueInfo &info);
int pdb_residue_info_get_serial_number(const AtomPDBResidueInfo &info);
std::unique_ptr<std::string> pdb_residue_info_get_alt_loc(
    const AtomPDBResidueInfo &info);
std::unique_ptr<std::string> pdb_residue_info_get_insertion_code(
    const AtomPDBResidueInfo &info);
double pdb_residue_info_get_occupancy(const AtomPDBResidueInfo &info);
double pdb_residue_info_get_temp_factor(const AtomPDBResidueInfo &info);
bool pdb_residue_info_get_is_hetero_atom(const AtomPDBResidueInfo &info);
unsigned pdb_residue_info_get_secondary_structure(
    const AtomPDBResidueInfo &info);
unsigned pdb_residue_info_get_segment_number(const AtomPDBResidueInfo &info);
void pdb_residue_info_set_serial_number(AtomPDBResidueInfo &info, int val);
void pdb_residue_info_set_alt_loc(AtomPDBResidueInfo &info, rust::Str val);
void pdb_residue_info_set_insertion_code(AtomPDBResidueInfo &info,
                                         rust::Str val);
void pdb_residue_info_set_occupancy(AtomPDBResidueInfo &info, double val);
void pdb_residue_info_set_temp_factor(AtomPDBResidueInfo &info, double val);
void pdb_residue_info_set_is_hetero_atom(AtomPDBResidueInfo &info, bool val);
void pdb_residue_info_set_secondary_structure(AtomPDBResidueInfo &info,
                                              unsigned val);
void pdb_residue_info_set_segment_number(AtomPDBResidueInfo &info,
                                         unsigned val);

// Copies of the atom's info, or null; the PDB form throws ValueError for
// other monomer types, as AtomGetPDBResidueInfo in Wrap/Atom.cpp does, and
// checks the info's class with dynamic_cast where Python casts.
std::unique_ptr<AtomMonomerInfo> atom_get_monomer_info(const Atom &atom);
std::unique_ptr<AtomPDBResidueInfo> atom_get_pdb_residue_info(const Atom &atom);
void atom_set_monomer_info(Atom &atom, const AtomMonomerInfo &info);

}  // namespace rdkit_cxx
