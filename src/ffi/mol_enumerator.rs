// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// The operations Wrap/rdMolEnumerator.cpp names (`EnumeratorTypes`,
    /// which RDKit's headers lack).
    #[repr(u8)]
    enum EnumeratorType {
        LinkNode,
        PositionVariation,
        RepeatUnit,
    }

    /// `RDKit::MolEnumerator::MolEnumeratorParams`, with its operation, if
    /// any, by type.
    struct MolEnumeratorParamsData {
        sanitize: bool,
        max_to_enumerate: usize,
        do_random: bool,
        random_seed: i32,
        has_operation: bool,
        operation: EnumeratorType,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_enumerator.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;

        fn mol_enumerator_params_default() -> Result<MolEnumeratorParamsData>;
        fn mol_enumerator_enumerate(
            mol: &RWMol,
            max_per_operation: u32,
        ) -> Result<UniquePtr<MolList>>;
        fn mol_enumerator_enumerate_with_params(
            mol: &RWMol,
            params: &MolEnumeratorParamsData,
        ) -> Result<UniquePtr<MolList>>;
    }
}
