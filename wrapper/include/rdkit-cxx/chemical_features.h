/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdChemicalFeatures: free chemical features, not tied to a molecule.

#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include <ChemicalFeatures/FreeChemicalFeature.h>

#include "rdkit-cxx/geometry.h"
#include "rust/cxx.h"

namespace rdkit_cxx {

using FreeChemicalFeature = ChemicalFeatures::FreeChemicalFeature;

std::unique_ptr<FreeChemicalFeature> free_chemical_feature_new();
std::unique_ptr<FreeChemicalFeature> free_chemical_feature_from_parts(
    rust::Str family, rust::Str type, const Point3D &loc, std::int32_t id);
// RDKit reads a pickle's lengths without checking them; the pickle is checked
// first (see chemical_features.cpp).
std::unique_ptr<FreeChemicalFeature> free_chemical_feature_from_pickle(
    rust::Slice<const std::uint8_t> pickle);
std::unique_ptr<FreeChemicalFeature> free_chemical_feature_clone(
    const FreeChemicalFeature &feat);
rust::Vec<std::uint8_t> free_chemical_feature_to_pickle(
    const FreeChemicalFeature &feat);

std::int32_t free_chemical_feature_get_id(const FreeChemicalFeature &feat);
std::unique_ptr<std::string> free_chemical_feature_get_family(
    const FreeChemicalFeature &feat);
std::unique_ptr<std::string> free_chemical_feature_get_type(
    const FreeChemicalFeature &feat);
Point3D free_chemical_feature_get_pos(const FreeChemicalFeature &feat);
void free_chemical_feature_set_id(FreeChemicalFeature &feat, std::int32_t id);
void free_chemical_feature_set_family(FreeChemicalFeature &feat,
                                      rust::Str family);
void free_chemical_feature_set_type(FreeChemicalFeature &feat, rust::Str type);
void free_chemical_feature_set_pos(FreeChemicalFeature &feat,
                                   const Point3D &loc);

}  // namespace rdkit_cxx
