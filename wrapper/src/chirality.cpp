/*
 * Copyright (c) 2026 Steven Kearnes. All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

#include "rdkit-cxx/chirality.h"

#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SanitException.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/chirality.rs.h"

namespace rdkit_cxx {

namespace {

// Throws the problem as its most derived type, so encode_current_exception
// sees the atom indices.
[[noreturn]] void rethrow(const RDKit::MolSanitizeException &e) {
  if (auto *x = dynamic_cast<const RDKit::AtomValenceException *>(&e)) {
    throw *x;
  }
  if (auto *x = dynamic_cast<const RDKit::AtomKekulizeException *>(&e)) {
    throw *x;
  }
  if (auto *x = dynamic_cast<const RDKit::KekulizeException *>(&e)) {
    throw *x;
  }
  if (auto *x = dynamic_cast<const RDKit::AtomSanitizeException *>(&e)) {
    throw *x;
  }
  throw e;
}

}  // namespace

void assign_stereochemistry(RWMol &mol, bool clean_it, bool force,
                            bool flag_possible_stereo_centers) {
  RDKit::MolOps::assignStereochemistry(mol, clean_it, force,
                                       flag_possible_stereo_centers);
}

void assign_stereochemistry_from_3d(RWMol &mol, std::int32_t conf_id,
                                    bool replace_existing_tags) {
  RDKit::MolOps::assignStereochemistryFrom3D(mol, conf_id,
                                             replace_existing_tags);
}

void assign_chiral_types_from_bond_dirs(RWMol &mol, std::int32_t conf_id,
                                        bool replace_existing_tags) {
  RDKit::MolOps::assignChiralTypesFromBondDirs(mol, conf_id,
                                               replace_existing_tags);
}

void assign_atom_chiral_tags_from_structure(RWMol &mol, std::int32_t conf_id,
                                            bool replace_existing_tags) {
  RDKit::MolOps::assignChiralTypesFrom3D(mol, conf_id, replace_existing_tags);
}

void assign_atom_chiral_tags_from_mol_parity(RWMol &mol,
                                             bool replace_existing_tags) {
  RDKit::MolOps::assignChiralTypesFromMolParity(mol, replace_existing_tags);
}

void find_potential_stereo_bonds(RWMol &mol, bool clean_it) {
  RDKit::MolOps::findPotentialStereoBonds(mol, clean_it);
}

void remove_stereochemistry(RWMol &mol) {
  RDKit::MolOps::removeStereochemistry(mol);
}

void remove_non_explicit_3d_chirality(RWMol &mol) {
  RDKit::Chirality::removeNonExplicit3DChirality(mol);
}

void detect_bond_stereochemistry(RWMol &mol, std::int32_t conf_id) {
  RDKit::MolOps::detectBondStereochemistry(mol, conf_id);
}

void set_double_bond_neighbor_directions(RWMol &mol) {
  RDKit::MolOps::setDoubleBondNeighborDirections(mol, nullptr);
}

// RDKit reads the conformer's position for each atom, through its own range
// check.
void set_double_bond_neighbor_directions_with_conf(RWMol &mol,
                                                   const Conformer &conf) {
  RDKit::MolOps::setDoubleBondNeighborDirections(mol, &conf);
}

void set_bond_stereo_from_directions(RWMol &mol) {
  RDKit::MolOps::setBondStereoFromDirections(mol);
}

rust::Vec<std::uint32_t> compute_atom_cip_ranks(RWMol &mol) {
  RDKit::UINT_VECT ranks;
  RDKit::Chirality::assignAtomCIPRanks(mol, ranks);
  rust::Vec<std::uint32_t> result;
  for (auto rank : ranks) {
    result.push_back(rank);
  }
  return result;
}

rust::Vec<IndexPair> find_meso_centers(const RWMol &mol, bool include_isotopes,
                                       bool include_atom_maps) {
  rust::Vec<IndexPair> result;
  for (const auto &[first, second] : RDKit::Chirality::findMesoCenters(
           mol, include_isotopes, include_atom_maps)) {
    result.push_back(IndexPair{first, second});
  }
  return result;
}

// As canonicalizeStereoGroupsHelper in Code/GraphMol/nbWrap/MolOps.cpp: a
// copy gets the canonical groups.
std::unique_ptr<RWMol> canonicalize_stereo_groups(
    const RWMol &mol, StereoGroupAbsOptions output_absolute_groups,
    std::uint32_t max_stereo_groups) {
  auto copy = std::make_unique<RDKit::ROMol>(mol);
  RDKit::canonicalizeStereoGroups(copy, output_absolute_groups,
                                  max_stereo_groups);
  return std::make_unique<RWMol>(*copy);
}

void cleanup_stereo_groups(RWMol &mol) {
  RDKit::Chirality::cleanupStereoGroups(mol);
}

void simplify_enhanced_stereo(RWMol &mol, bool remove_affected_stereo_groups) {
  RDKit::Chirality::simplifyEnhancedStereo(mol, remove_affected_stereo_groups);
}

void translate_chiral_flag_to_stereo_groups(
    RWMol &mol, StereoGroupType zero_flag_group_type) {
  RDKit::translateChiralFlagToStereoGroups(mol, zero_flag_group_type);
}

void add_stereo_annotations(RWMol &mol, rust::Str abs_label, rust::Str or_label,
                            rust::Str and_label, rust::Str cip_label,
                            rust::Str bond_label) {
  RDKit::Chirality::addStereoAnnotations(
      mol, std::string(abs_label), std::string(or_label),
      std::string(and_label), std::string(cip_label), std::string(bond_label));
}

void add_wavy_bonds_for_stereo_any(RWMol &mol, bool clear_double_bond_flags,
                                   std::uint32_t add_when_impossible) {
  RDKit::addWavyBondsForStereoAny(mol, clear_double_bond_flags,
                                  add_when_impossible);
}

void wedge_bond(Bond &bond, std::uint32_t from_atom_idx, std::int32_t conf_id) {
  const auto &conf = bond.getOwningMol().getConformer(conf_id);
  RDKit::Chirality::wedgeBond(&bond, from_atom_idx, &conf);
}

BondWedgingParametersData bond_wedging_parameters_default() {
  RDKit::Chirality::BondWedgingParameters params;
  return BondWedgingParametersData{params.wedgeTwoBondsIfPossible};
}

void wedge_mol_bonds(RWMol &mol, const Conformer &conformer,
                     const BondWedgingParametersData &data) {
  RDKit::Chirality::BondWedgingParameters params;
  params.wedgeTwoBondsIfPossible = data.wedge_two_bonds_if_possible;
  RDKit::Chirality::wedgeMolBonds(mol, &conformer, &params);
}

void reapply_mol_block_wedging(RWMol &mol, bool all_bond_types, bool verify) {
  RDKit::Chirality::reapplyMolBlockWedging(mol, all_bond_types, verify);
}

void set_allow_nontetrahedral_chirality(bool val) {
  RDKit::Chirality::setAllowNontetrahedralChirality(val);
}

bool get_allow_nontetrahedral_chirality() {
  return RDKit::Chirality::getAllowNontetrahedralChirality();
}

void set_use_legacy_stereo_perception(bool val) {
  RDKit::Chirality::setUseLegacyStereoPerception(val);
}

bool get_use_legacy_stereo_perception() {
  return RDKit::Chirality::getUseLegacyStereoPerception();
}

void expand_attachment_points(RWMol &mol, bool add_as_queries,
                              bool add_coords) {
  RDKit::MolOps::expandAttachmentPoints(mol, add_as_queries, add_coords);
}

void collapse_attachment_points(RWMol &mol, bool marked_only) {
  RDKit::MolOps::collapseAttachmentPoints(mol, marked_only);
}

std::unique_ptr<std::vector<std::string>> detect_chemistry_problems(
    const RWMol &mol, std::uint32_t sanitize_ops) {
  auto result = std::make_unique<std::vector<std::string>>();
  for (const auto &problem :
       RDKit::MolOps::detectChemistryProblems(mol, sanitize_ops)) {
    try {
      rethrow(*problem);
    } catch (...) {
      result->push_back(encode_current_exception());
    }
  }
  return result;
}

}  // namespace rdkit_cxx
