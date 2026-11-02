// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::EnumerateStereoisomers::StereoEnumerationOptions`.
    struct StereoEnumerationOptionsData {
        try_embedding: bool,
        only_unassigned: bool,
        only_stereo_groups: bool,
        unique: bool,
        max_isomers: u64,
        random_seed: i32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/enumerate_stereoisomers.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type StereoisomerEnumerator;

        fn stereo_enumeration_options_default() -> Result<StereoEnumerationOptionsData>;
        fn stereoisomer_enumerator_new(
            mol: &RWMol,
            options: &StereoEnumerationOptionsData,
            verbose: bool,
        ) -> Result<UniquePtr<StereoisomerEnumerator>>;
        fn stereoisomer_enumerator_get_stereoisomer_count(
            enumerator: &StereoisomerEnumerator,
        ) -> Result<u64>;
        fn stereoisomer_enumerator_next(
            enumerator: Pin<&mut StereoisomerEnumerator>,
        ) -> Result<UniquePtr<RWMol>>;
    }
}

// SAFETY: an enumerator owns its copy of the molecule, which its flippers
// point into, its random generator, and the sets of isomers it has made;
// nothing outside it shares them, so it may move between threads. It is not
// Sync: next() changes all of them.
unsafe impl Send for ffi::StereoisomerEnumerator {}
