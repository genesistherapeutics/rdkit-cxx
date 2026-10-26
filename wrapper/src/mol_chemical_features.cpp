// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/mol_chemical_features.h"

#include <iterator>
#include <mutex>

#include <GraphMol/MolChemicalFeatures/FeatureParser.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureDef.h>
#include <RDGeneral/FileParseException.h>

#include "rdkit-cxx/control_c_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/geometry.rs.h"

namespace rdkit_cxx {

namespace {

Point3D to_point(const RDGeom::Point3D &p) { return Point3D{p.x, p.y, p.z}; }

const RDKit::MolChemicalFeatureDef &feature_def(
    const MolChemicalFeatureFactory &factory, std::uint32_t idx) {
  auto n = static_cast<unsigned>(factory.getNumFeatureDefs());
  check_index(idx, n);
  auto it = factory.beginFeatureDefs();
  std::advance(it, idx);
  return **it;
}

// getFeaturesForMol matches every definition's pattern, and each match
// constructs a ControlCHandler (control_c_lock.h, CONVENTIONS.md U4).
RDKit::FeatSPtrList features(const MolChemicalFeatureFactory &factory,
                             const RWMol &mol, rust::Str include_only,
                             std::int32_t conf_id) {
  // Python passes includeOnly as a C string, which ends at a NUL.
  std::string only(include_only);
  std::lock_guard<ControlCMutex> control_c(control_c_mutex());
  return factory.getFeaturesForMol(mol, only.c_str(), conf_id);
}

}  // namespace

// The parser checks every pattern as SMARTS before it makes a definition
// from it, so no definition has a null pattern, which getFeaturesForMol
// would dereference.
std::unique_ptr<MolChemicalFeatureFactory> build_feature_factory_from_string(
    rust::Str fdef) {
  std::unique_ptr<MolChemicalFeatureFactory> factory;
  try {
    factory.reset(RDKit::buildFeatureFactory(std::string(fdef)));
  } catch (const RDKit::FeatureFileParseException &e) {
    // Python's FeatureFileParseException is a ValueError subclass; it is a
    // parse error here (E4).
    throw RDKit::FileParseException(std::string(e.what()) + " (line " +
                                    std::to_string(e.lineNo()) + ")");
  }
  if (!factory) {
    throw RDKit::FileParseException("could not build a feature factory");
  }
  return factory;
}

std::int32_t factory_get_num_feature_defs(
    const MolChemicalFeatureFactory &factory) {
  return factory.getNumFeatureDefs();
}

std::unique_ptr<std::string> factory_feature_def_family(
    const MolChemicalFeatureFactory &factory, std::uint32_t idx) {
  return std::make_unique<std::string>(feature_def(factory, idx).getFamily());
}

std::unique_ptr<std::string> factory_feature_def_type(
    const MolChemicalFeatureFactory &factory, std::uint32_t idx) {
  return std::make_unique<std::string>(feature_def(factory, idx).getType());
}

std::unique_ptr<std::string> factory_feature_def_smarts(
    const MolChemicalFeatureFactory &factory, std::uint32_t idx) {
  return std::make_unique<std::string>(feature_def(factory, idx).getSmarts());
}

std::int32_t factory_get_num_mol_features(
    const MolChemicalFeatureFactory &factory, const RWMol &mol,
    rust::Str include_only) {
  return static_cast<std::int32_t>(
      features(factory, mol, include_only, -1).size());
}

// As getMolFeature in nbWrap/MolChemicalFeatureFactory.cpp, which keeps the
// features in a static list that recompute=False reads back: a list shared by
// every thread, holding features of molecules Python may have freed. The
// features are computed each time here.
std::unique_ptr<MolChemicalFeature> factory_get_mol_feature(
    const MolChemicalFeatureFactory &factory, const RWMol &mol,
    std::int32_t idx, rust::Str include_only, std::int32_t conf_id) {
  auto feats = features(factory, mol, include_only, conf_id);
  if (idx < 0 || static_cast<std::size_t>(idx) >= feats.size()) {
    throw IndexError(idx, feats.size());
  }
  auto it = feats.begin();
  std::advance(it, idx);
  auto result = std::make_unique<MolChemicalFeature>();
  result->feat = *it;
  return result;
}

std::int32_t mol_chemical_feature_get_id(const MolChemicalFeature &feat) {
  return feat.feat->getId();
}

std::unique_ptr<std::string> mol_chemical_feature_get_family(
    const MolChemicalFeature &feat) {
  return std::make_unique<std::string>(feat.feat->getFamily());
}

std::unique_ptr<std::string> mol_chemical_feature_get_type(
    const MolChemicalFeature &feat) {
  return std::make_unique<std::string>(feat.feat->getType());
}

Point3D mol_chemical_feature_get_pos(const MolChemicalFeature &feat) {
  return to_point(feat.feat->getPos());
}

Point3D mol_chemical_feature_get_pos_with_conf_id(
    const MolChemicalFeature &feat, std::int32_t conf_id) {
  return to_point(feat.feat->getPos(conf_id));
}

rust::Vec<std::uint32_t> mol_chemical_feature_get_atom_ids(
    const MolChemicalFeature &feat) {
  rust::Vec<std::uint32_t> ids;
  for (const auto *atom : feat.feat->getAtoms()) {
    ids.push_back(atom->getIdx());
  }
  return ids;
}

std::int32_t mol_chemical_feature_get_active_conformer(
    const MolChemicalFeature &feat) {
  return feat.feat->getActiveConformer();
}

void mol_chemical_feature_set_active_conformer(MolChemicalFeature &feat,
                                               std::int32_t conf_id) {
  feat.feat->setActiveConformer(conf_id);
}

void mol_chemical_feature_clear_cache(MolChemicalFeature &feat) {
  feat.feat->clearCache();
}

}  // namespace rdkit_cxx
