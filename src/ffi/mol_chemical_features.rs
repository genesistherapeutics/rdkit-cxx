// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_chemical_features.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type MolChemicalFeatureFactory;
        type MolChemicalFeature;

        fn build_feature_factory_from_string(
            fdef: &str,
        ) -> Result<UniquePtr<MolChemicalFeatureFactory>>;
        fn factory_get_num_feature_defs(factory: &MolChemicalFeatureFactory) -> Result<i32>;
        fn factory_feature_def_family(
            factory: &MolChemicalFeatureFactory,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn factory_feature_def_type(
            factory: &MolChemicalFeatureFactory,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn factory_feature_def_smarts(
            factory: &MolChemicalFeatureFactory,
            idx: u32,
        ) -> Result<UniquePtr<CxxString>>;
        fn factory_get_num_mol_features(
            factory: &MolChemicalFeatureFactory,
            mol: &RWMol,
            include_only: &str,
        ) -> Result<i32>;
        fn factory_get_mol_feature(
            factory: &MolChemicalFeatureFactory,
            mol: &RWMol,
            idx: i32,
            include_only: &str,
            conf_id: i32,
        ) -> Result<UniquePtr<MolChemicalFeature>>;

        fn mol_chemical_feature_get_id(feat: &MolChemicalFeature) -> Result<i32>;
        fn mol_chemical_feature_get_family(
            feat: &MolChemicalFeature,
        ) -> Result<UniquePtr<CxxString>>;
        fn mol_chemical_feature_get_type(feat: &MolChemicalFeature)
        -> Result<UniquePtr<CxxString>>;
        fn mol_chemical_feature_get_pos(feat: &MolChemicalFeature) -> Result<Point3D>;
        fn mol_chemical_feature_get_pos_with_conf_id(
            feat: &MolChemicalFeature,
            conf_id: i32,
        ) -> Result<Point3D>;
        fn mol_chemical_feature_get_atom_ids(feat: &MolChemicalFeature) -> Result<Vec<u32>>;
        fn mol_chemical_feature_get_active_conformer(feat: &MolChemicalFeature) -> Result<i32>;
        fn mol_chemical_feature_set_active_conformer(
            feat: Pin<&mut MolChemicalFeature>,
            conf_id: i32,
        ) -> Result<()>;
        fn mol_chemical_feature_clear_cache(feat: Pin<&mut MolChemicalFeature>) -> Result<()>;
    }
}

// SAFETY: a factory owns its feature definitions, whose strings and patterns
// nothing changes after parsing. Its only call that touches the patterns,
// getFeaturesForMol, matches them holding control_c_mutex, so no two threads
// match them at once, and RDKit locks each recursive query while it records
// matches in it (RDK_BUILD_THREADSAFE_SSS, CONVENTIONS.md S5).
unsafe impl Send for ffi::MolChemicalFeatureFactory {}
unsafe impl Sync for ffi::MolChemicalFeatureFactory {}
