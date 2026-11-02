// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/chemical_features.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type Point3D = crate::ffi::geometry::ffi::Point3D;
        type FreeChemicalFeature;

        fn free_chemical_feature_new() -> Result<UniquePtr<FreeChemicalFeature>>;
        fn free_chemical_feature_from_parts(
            family: &str,
            type_: &str,
            loc: &Point3D,
            id: i32,
        ) -> Result<UniquePtr<FreeChemicalFeature>>;
        fn free_chemical_feature_from_pickle(
            pickle: &[u8],
        ) -> Result<UniquePtr<FreeChemicalFeature>>;
        fn free_chemical_feature_clone(
            feat: &FreeChemicalFeature,
        ) -> Result<UniquePtr<FreeChemicalFeature>>;
        fn free_chemical_feature_to_pickle(feat: &FreeChemicalFeature) -> Result<Vec<u8>>;

        fn free_chemical_feature_get_id(feat: &FreeChemicalFeature) -> Result<i32>;
        fn free_chemical_feature_get_family(
            feat: &FreeChemicalFeature,
        ) -> Result<UniquePtr<CxxString>>;
        fn free_chemical_feature_get_type(
            feat: &FreeChemicalFeature,
        ) -> Result<UniquePtr<CxxString>>;
        fn free_chemical_feature_get_pos(feat: &FreeChemicalFeature) -> Result<Point3D>;
        fn free_chemical_feature_set_id(feat: Pin<&mut FreeChemicalFeature>, id: i32)
        -> Result<()>;
        fn free_chemical_feature_set_family(
            feat: Pin<&mut FreeChemicalFeature>,
            family: &str,
        ) -> Result<()>;
        fn free_chemical_feature_set_type(
            feat: Pin<&mut FreeChemicalFeature>,
            type_: &str,
        ) -> Result<()>;
        fn free_chemical_feature_set_pos(
            feat: Pin<&mut FreeChemicalFeature>,
            loc: &Point3D,
        ) -> Result<()>;
    }
}

// SAFETY: a FreeChemicalFeature owns its strings and position and refers to
// nothing else; its const methods only read them.
unsafe impl Send for ffi::FreeChemicalFeature {}
unsafe impl Sync for ffi::FreeChemicalFeature {}
