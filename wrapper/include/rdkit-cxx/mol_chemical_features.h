/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolChemicalFeatures: feature factories and the features they find in
// molecules.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

#include "rdkit-cxx/geometry.h"
#include "rdkit-cxx/mol.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using MolChemicalFeatureFactory = RDKit::MolChemicalFeatureFactory;

// A feature the factory found. It points at the molecule it was found in
// (and that molecule's atoms), at the factory, and at the factory's feature
// definition; its Rust handle borrows the molecule and the factory for as
// long as it lives.
class MolChemicalFeature {
 public:
  RDKit::FeatSPtr feat;
};

// As buildFeatureFactory; FeatureFileParseException becomes
// FileParseException, with its line number.
std::unique_ptr<MolChemicalFeatureFactory> build_feature_factory_from_string(
    rust::Str fdef);

std::int32_t factory_get_num_feature_defs(
    const MolChemicalFeatureFactory &factory);
// The definitions' families, types, and SMARTS, in order (E3).
std::unique_ptr<std::string> factory_feature_def_family(
    const MolChemicalFeatureFactory &factory, std::uint32_t idx);
std::unique_ptr<std::string> factory_feature_def_type(
    const MolChemicalFeatureFactory &factory, std::uint32_t idx);
std::unique_ptr<std::string> factory_feature_def_smarts(
    const MolChemicalFeatureFactory &factory, std::uint32_t idx);
// As getFeaturesForMol, which matches each definition's pattern.
std::int32_t factory_get_num_mol_features(
    const MolChemicalFeatureFactory &factory, const RWMol &mol,
    rust::Str include_only);
std::unique_ptr<MolChemicalFeature> factory_get_mol_feature(
    const MolChemicalFeatureFactory &factory, const RWMol &mol,
    std::int32_t idx, rust::Str include_only, std::int32_t conf_id);

std::int32_t mol_chemical_feature_get_id(const MolChemicalFeature &feat);
std::unique_ptr<std::string> mol_chemical_feature_get_family(
    const MolChemicalFeature &feat);
std::unique_ptr<std::string> mol_chemical_feature_get_type(
    const MolChemicalFeature &feat);
Point3D mol_chemical_feature_get_pos(const MolChemicalFeature &feat);
Point3D mol_chemical_feature_get_pos_with_conf_id(
    const MolChemicalFeature &feat, std::int32_t conf_id);
rust::Vec<std::uint32_t> mol_chemical_feature_get_atom_ids(
    const MolChemicalFeature &feat);
std::int32_t mol_chemical_feature_get_active_conformer(
    const MolChemicalFeature &feat);
void mol_chemical_feature_set_active_conformer(MolChemicalFeature &feat,
                                               std::int32_t conf_id);
void mol_chemical_feature_clear_cache(MolChemicalFeature &feat);

}  // namespace rdkit_cxx
