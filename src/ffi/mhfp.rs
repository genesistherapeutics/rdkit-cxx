// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// The shingling arguments every `MHFPEncoder` method takes.
    struct ShinglingOptions {
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: bool,
        min_radius: u8,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mhfp.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type MHFPEncoder;

        fn mhfp_encoder_new(n_permutations: u32, seed: u32) -> Result<UniquePtr<MHFPEncoder>>;
        fn mhfp_from_string_array(encoder: &MHFPEncoder, vec: &[&str]) -> Result<Vec<u32>>;
        fn mhfp_from_array(encoder: &MHFPEncoder, vec: &[u32]) -> Result<Vec<u32>>;
        fn mhfp_shingling_from_smiles(
            encoder: &MHFPEncoder,
            smiles: &str,
            options: &ShinglingOptions,
        ) -> Result<Vec<String>>;
        fn mhfp_shingling_from_mol(
            encoder: &MHFPEncoder,
            mol: &RWMol,
            options: &ShinglingOptions,
        ) -> Result<Vec<String>>;
        fn mhfp_encode_smiles(
            encoder: &MHFPEncoder,
            smiles: &str,
            options: &ShinglingOptions,
        ) -> Result<Vec<u32>>;
        fn mhfp_encode_mol(
            encoder: &MHFPEncoder,
            mol: &RWMol,
            options: &ShinglingOptions,
        ) -> Result<Vec<u32>>;
        fn mhfp_encode_secfp_smiles(
            encoder: &MHFPEncoder,
            smiles: &str,
            options: &ShinglingOptions,
            length: usize,
        ) -> Result<UniquePtr<BitVect>>;
        fn mhfp_encode_secfp_mol(
            encoder: &MHFPEncoder,
            mol: &RWMol,
            options: &ShinglingOptions,
            length: usize,
        ) -> Result<UniquePtr<BitVect>>;
        fn mhfp_distance(a: &[u32], b: &[u32]) -> Result<f64>;
    }
}

// SAFETY: an encoder owns its permutations and has no thread affinity.
unsafe impl Send for ffi::MHFPEncoder {}
// SAFETY: every method the bridge binds only reads the encoder (the C++
// methods are non-const but write nothing), so it can be shared.
unsafe impl Sync for ffi::MHFPEncoder {}
