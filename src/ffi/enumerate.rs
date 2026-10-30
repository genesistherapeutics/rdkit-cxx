// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `EnumerationParams`.
    struct EnumerationParamsData {
        reagent_max_match_count: i32,
        sane_partial_products: bool,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/enumerate.h");

        type ChemicalReaction = crate::ffi::chem_reactions::ffi::ChemicalReaction;
        type ProductSets = crate::ffi::chem_reactions::ffi::ProductSets;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type EnumerationStrategyBase;
        type EnumerateLibrary;

        fn enumeration_strategy_new(kind: u8) -> Result<UniquePtr<EnumerationStrategyBase>>;
        fn enumeration_strategy_clone(
            strategy: &EnumerationStrategyBase,
        ) -> Result<UniquePtr<EnumerationStrategyBase>>;
        fn enumeration_strategy_initialize(
            strategy: Pin<&mut EnumerationStrategyBase>,
            rxn: &ChemicalReaction,
            mols: Pin<&mut BorrowedMols>,
            group_sizes: &[u32],
        ) -> Result<()>;
        fn enumeration_strategy_get_position(
            strategy: &EnumerationStrategyBase,
        ) -> Result<Vec<u64>>;
        fn enumeration_strategy_get_num_permutations(
            strategy: &EnumerationStrategyBase,
        ) -> Result<u64>;
        fn enumeration_strategy_type(
            strategy: &EnumerationStrategyBase,
        ) -> Result<UniquePtr<CxxString>>;
        fn enumeration_strategy_has_next(strategy: &EnumerationStrategyBase) -> Result<bool>;
        fn enumeration_strategy_next(
            strategy: Pin<&mut EnumerationStrategyBase>,
        ) -> Result<Vec<u64>>;
        fn enumeration_strategy_skip(
            strategy: Pin<&mut EnumerationStrategyBase>,
            skip_count: u64,
        ) -> Result<bool>;
        fn even_sample_pairs_strategy_stats(
            strategy: &EnumerationStrategyBase,
        ) -> Result<UniquePtr<CxxString>>;

        fn enumeration_params_default() -> Result<EnumerationParamsData>;
        fn enumerate_library_can_serialize() -> Result<bool>;

        fn enumerate_library_new() -> Result<UniquePtr<EnumerateLibrary>>;
        fn enumerate_library_new_with_params(
            rxn: &ChemicalReaction,
            mols: Pin<&mut BorrowedMols>,
            group_sizes: &[u32],
            params: &EnumerationParamsData,
        ) -> Result<UniquePtr<EnumerateLibrary>>;
        fn enumerate_library_new_with_enumerator(
            rxn: &ChemicalReaction,
            mols: Pin<&mut BorrowedMols>,
            group_sizes: &[u32],
            strategy: &EnumerationStrategyBase,
            params: &EnumerationParamsData,
        ) -> Result<UniquePtr<EnumerateLibrary>>;
        fn enumerate_library_get_reagents(lib: &EnumerateLibrary)
        -> Result<UniquePtr<ProductSets>>;
        fn enumerate_library_get_reaction(
            lib: &EnumerateLibrary,
        ) -> Result<UniquePtr<ChemicalReaction>>;
        fn enumerate_library_get_enumerator(
            lib: &EnumerateLibrary,
        ) -> Result<UniquePtr<EnumerationStrategyBase>>;
        fn enumerate_library_get_position(lib: &EnumerateLibrary) -> Result<Vec<u64>>;
        fn enumerate_library_has_next(lib: &EnumerateLibrary) -> Result<bool>;
        fn enumerate_library_next(
            lib: Pin<&mut EnumerateLibrary>,
        ) -> Result<UniquePtr<ProductSets>>;
        fn enumerate_library_get_state(lib: &EnumerateLibrary) -> Result<UniquePtr<CxxString>>;
        fn enumerate_library_set_state(lib: Pin<&mut EnumerateLibrary>, state: &str) -> Result<()>;
        fn enumerate_library_reset_state(lib: Pin<&mut EnumerateLibrary>) -> Result<()>;
        fn enumerate_library_serialize(lib: &EnumerateLibrary) -> Result<Vec<u8>>;
        fn enumerate_library_is_empty(lib: &EnumerateLibrary) -> Result<bool>;
        fn enumerate_library_from_string(data: &[u8]) -> Result<UniquePtr<EnumerateLibrary>>;
    }
}

// SAFETY: a strategy owns its state (positions, counts, a random number
// generator), and reads nothing else once initialized.
unsafe impl Send for ffi::EnumerationStrategyBase {}

// SAFETY: a library owns its reaction (a copy), its strategies, and its
// building blocks (copies the wrappers make, or molecules read from its
// pickle), so moving it moves everything it reads. It is not Sync: running
// it writes its building blocks' bookmarks and its reaction's templates, as
// running a reaction does (CONVENTIONS.md S4).
unsafe impl Send for ffi::EnumerateLibrary {}
