// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::RGroupLabels`; cxx checks the values.
    #[repr(u32)]
    enum RGroupLabels {
        IsotopeLabels = 0x01,
        AtomMapLabels = 0x02,
        AtomIndexLabels = 0x04,
        RelabelDuplicateLabels = 0x08,
        MDLRGroupLabels = 0x10,
        DummyAtomLabels = 0x20,
        AutoDetect = 0xFF,
    }

    /// `RDKit::RGroupMatching`; cxx checks the values.
    #[repr(u32)]
    enum RGroupMatching {
        Greedy = 0x01,
        GreedyChunks = 0x02,
        Exhaustive = 0x04,
        NoSymmetrization = 0x08,
        GA = 0x10,
    }

    /// `RDKit::RGroupLabelling`; cxx checks the values.
    #[repr(u32)]
    enum RGroupLabelling {
        AtomMap = 0x01,
        Isotope = 0x02,
        MDLRGroup = 0x04,
    }

    /// `RDKit::RGroupCoreAlignment`; cxx checks the values.
    #[repr(u32)]
    enum RGroupCoreAlignment {
        NoAlignment = 0x0,
        MCS = 0x01,
    }

    /// `RDKit::RGroupScore`; cxx checks the values.
    #[repr(u32)]
    enum RGroupScore {
        Match = 0x1,
        FingerprintVariance = 0x4,
    }

    /// `RGroupDecompositionParameters`.
    struct RGroupDecompositionParametersData {
        labels: u32,
        matching_strategy: u32,
        score_method: u32,
        rgroup_labelling: u32,
        alignment: u32,
        chunk_size: u32,
        only_match_at_r_groups: bool,
        remove_all_hydrogen_r_groups: bool,
        remove_all_hydrogen_r_groups_and_labels: bool,
        remove_hydrogens_post_match: bool,
        allow_non_terminal_r_groups: bool,
        allow_multiple_r_groups_on_unlabelled: bool,
        allow_multiple_cores_in_same_mol: bool,
        do_tautomers: bool,
        do_enumeration: bool,
        include_target_mol_in_results: bool,
        timeout: f64,
        ga_population_size: i32,
        ga_maximum_operations: i32,
        ga_number_operations_without_improvement: i32,
        ga_random_seed: i32,
        ga_number_runs: i32,
        ga_parallel_runs: bool,
        substruct_match_params: SubstructMatchParametersData,
    }

    /// Rows or columns of molecules: each named by its label, with the
    /// number in each row or column.
    struct RGroupTable {
        mols: UniquePtr<MolList>,
        sizes: Vec<u32>,
    }

    struct RGroupDecomposeResult {
        table: RGroupTable,
        unmatched: Vec<u32>,
    }

    struct ProcessAndScoreResult {
        success: bool,
        score: f64,
    }

    /// A core's index and its matches, each as core and molecule atom
    /// index pairs, flattened.
    struct MatchingCoreResult {
        core_idx: i32,
        matches: Vec<IndexList>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/rgroup.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");
        include!("rdkit-cxx/ffi/substruct.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type SubstructMatchParametersData =
            crate::ffi::substruct::ffi::SubstructMatchParametersData;
        type RGroupLabels;
        type RGroupMatching;
        type RGroupLabelling;
        type RGroupCoreAlignment;
        type RGroupScore;
        type RGroupDecomposition;

        fn r_group_decomposition_parameters_default() -> Result<RGroupDecompositionParametersData>;
        fn r_group_decomposition_new(
            cores: Pin<&mut BorrowedMols>,
            params: &RGroupDecompositionParametersData,
        ) -> Result<UniquePtr<RGroupDecomposition>>;
        fn r_group_decomposition_add(
            decomp: Pin<&mut RGroupDecomposition>,
            mol: &RWMol,
        ) -> Result<i32>;
        fn r_group_decomposition_get_matching_core_idx(
            decomp: Pin<&mut RGroupDecomposition>,
            mol: &RWMol,
        ) -> Result<MatchingCoreResult>;
        fn r_group_decomposition_process(decomp: Pin<&mut RGroupDecomposition>) -> Result<bool>;
        fn r_group_decomposition_process_and_score(
            decomp: Pin<&mut RGroupDecomposition>,
        ) -> Result<ProcessAndScoreResult>;
        fn r_group_decomposition_get_r_group_labels(
            decomp: Pin<&mut RGroupDecomposition>,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn r_group_decomposition_get_r_groups_as_rows(
            decomp: Pin<&mut RGroupDecomposition>,
        ) -> Result<RGroupTable>;
        fn r_group_decomposition_get_r_groups_as_columns(
            decomp: Pin<&mut RGroupDecomposition>,
        ) -> Result<RGroupTable>;
        fn r_group_decompose(
            cores: Pin<&mut BorrowedMols>,
            mols: Pin<&mut BorrowedMols>,
            params: &RGroupDecompositionParametersData,
            as_rows: bool,
        ) -> Result<RGroupDecomposeResult>;
        fn relabel_mapped_dummies(
            mol: Pin<&mut RWMol>,
            input_labels: u32,
            output_labels: u32,
        ) -> Result<()>;
    }
}

// SAFETY: a decomposition owns copies of its cores and molecules, its
// parameters, and its matches, so moving it moves everything it reads;
// the process-wide state it uses (the GA's random number generator, the
// scoring static) is used under rgroup_mutex. It is not Sync: its getters
// rewrite its matched cores in place (outputCoreMolecule), so they take
// &mut self (CONVENTIONS.md S4).
unsafe impl Send for ffi::RGroupDecomposition {}
