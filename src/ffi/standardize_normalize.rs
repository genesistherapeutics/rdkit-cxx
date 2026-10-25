// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolStandardize::MetalDisconnectorOptions` (A7).
    struct MetalDisconnectorOptionsData {
        split_grignards: bool,
        split_aromatic_c: bool,
        adjust_charges: bool,
        remove_haptic_dummies: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/standardize_normalize.h");
        include!("rdkit-cxx/ffi/standardize.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type CleanupParametersData = crate::ffi::standardize::ffi::CleanupParametersData;
        type Normalizer;
        type MetalDisconnector;

        fn normalizer_default() -> Result<UniquePtr<Normalizer>>;
        fn normalizer_new(
            normalize_filename: &[u8],
            max_restarts: u32,
        ) -> Result<UniquePtr<Normalizer>>;
        fn normalizer_from_data(
            param_data: &str,
            params: &CleanupParametersData,
        ) -> Result<UniquePtr<Normalizer>>;
        fn normalizer_from_params(params: &CleanupParametersData) -> Result<UniquePtr<Normalizer>>;
        fn normalizer_normalize(normalizer: &Normalizer, mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn normalizer_normalize_in_place(
            normalizer: &Normalizer,
            mol: Pin<&mut RWMol>,
        ) -> Result<()>;

        fn metal_disconnector_options_default() -> Result<MetalDisconnectorOptionsData>;
        fn metal_disconnector_new(
            options: &MetalDisconnectorOptionsData,
        ) -> Result<UniquePtr<MetalDisconnector>>;
        fn metal_disconnector_get_metal_nof(
            disconnector: &MetalDisconnector,
        ) -> Result<UniquePtr<CxxString>>;
        fn metal_disconnector_get_metal_non(
            disconnector: &MetalDisconnector,
        ) -> Result<UniquePtr<CxxString>>;
        fn metal_disconnector_set_metal_nof(
            disconnector: Pin<&mut MetalDisconnector>,
            mol: &RWMol,
        ) -> Result<()>;
        fn metal_disconnector_set_metal_non(
            disconnector: Pin<&mut MetalDisconnector>,
            mol: &RWMol,
        ) -> Result<()>;
        fn metal_disconnector_disconnect(
            disconnector: &MetalDisconnector,
            mol: &RWMol,
        ) -> Result<UniquePtr<RWMol>>;
        fn metal_disconnector_disconnect_in_place(
            disconnector: &MetalDisconnector,
            mol: Pin<&mut RWMol>,
        ) -> Result<()>;

        fn disconnect_organometallics(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn disconnect_organometallics_with_options(
            mol: &RWMol,
            options: &MetalDisconnectorOptionsData,
        ) -> Result<UniquePtr<RWMol>>;
        fn disconnect_organometallics_in_place(mol: Pin<&mut RWMol>) -> Result<()>;
        fn disconnect_organometallics_in_place_with_options(
            mol: Pin<&mut RWMol>,
            options: &MetalDisconnectorOptionsData,
        ) -> Result<()>;
    }
}

// SAFETY: a Normalizer owns its catalog, a deep copy of the transformations
// (TransformCatalogParams' copy constructor copies each reaction), and a
// MetalDisconnector owns its queries.
unsafe impl Send for ffi::Normalizer {}
unsafe impl Send for ffi::MetalDisconnector {}
