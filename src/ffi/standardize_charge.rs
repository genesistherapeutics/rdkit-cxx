// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolStandardize::ChargeCorrection` (A7), its strings as bytes.
    struct ChargeCorrectionData {
        name: Vec<u8>,
        smarts: Vec<u8>,
        charge: i32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/standardize_charge.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Reionizer;
        type Uncharger;

        fn charge_corrections() -> Result<Vec<ChargeCorrectionData>>;

        fn reionizer_default() -> Result<UniquePtr<Reionizer>>;
        fn reionizer_new(
            acidbase_file: &[u8],
            ccs: &[ChargeCorrectionData],
            default_ccs: bool,
        ) -> Result<UniquePtr<Reionizer>>;
        fn reionizer_from_data(
            param_data: &str,
            ccs: &[ChargeCorrectionData],
        ) -> Result<UniquePtr<Reionizer>>;
        fn reionizer_reionize(reionizer: &Reionizer, mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn reionizer_reionize_in_place(reionizer: &Reionizer, mol: Pin<&mut RWMol>) -> Result<()>;

        fn uncharger_new(
            canonical_order: bool,
            force: bool,
            protonation_only: bool,
        ) -> Result<UniquePtr<Uncharger>>;
        fn uncharger_uncharge(uncharger: &Uncharger, mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn uncharger_uncharge_in_place(uncharger: &Uncharger, mol: Pin<&mut RWMol>) -> Result<()>;
    }
}

// SAFETY: a Reionizer owns its catalog, whose acid-base pairs it shares
// (through shared_ptrs, with atomic counts) with RDKit's cache of parsed
// files, which only reads them; an Uncharger owns its patterns.
unsafe impl Send for ffi::Reionizer {}
unsafe impl Send for ffi::Uncharger {}
