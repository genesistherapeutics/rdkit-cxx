// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolStandardize::PipelineOptions` (A7), its strings as bytes.
    struct PipelineOptionsData {
        strict_parsing: bool,
        report_all_failures: bool,
        allow_empty_molecules: bool,
        allow_enhanced_stereo: bool,
        allow_aromatic_bond_type: bool,
        allow_dative_bond_type: bool,
        is_2d_zero_threshold: f64,
        atom_clash_limit: f64,
        min_median_bond_length: f64,
        bond_length_limit: f64,
        allow_long_bonds_in_rings: bool,
        allow_atom_bond_clash_exemption: bool,
        metal_nof: Vec<u8>,
        metal_non: Vec<u8>,
        normalizer_data: Vec<u8>,
        normalizer_max_restarts: u32,
        scaled_median_bond_length: f64,
        output_v2000: bool,
    }

    /// `RDKit::MolStandardize::PipelineLogEntry`, its detail as bytes.
    struct PipelineLogEntryData {
        status: u32,
        detail: Vec<u8>,
    }

    /// `RDKit::MolStandardize::PipelineResult`, its strings as bytes.
    struct PipelineResultData {
        status: u32,
        stage: u32,
        log: Vec<PipelineLogEntryData>,
        input_mol_data: Vec<u8>,
        output_mol_data: Vec<u8>,
        parent_mol_data: Vec<u8>,
    }

    /// `RDKit::MolStandardize::PipelineStatus`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum PipelineStatus {
        NO_EVENT = 0,
        INPUT_ERROR = 0x1,
        PREPARE_FOR_VALIDATION_ERROR = 0x2,
        FEATURES_VALIDATION_ERROR = 0x4,
        BASIC_VALIDATION_ERROR = 0x8,
        IS2D_VALIDATION_ERROR = 0x10,
        LAYOUT2D_VALIDATION_ERROR = 0x20,
        STEREO_VALIDATION_ERROR = 0x40,
        VALIDATION_ERROR = 0x7c,
        PREPARE_FOR_STANDARDIZATION_ERROR = 0x80,
        METAL_STANDARDIZATION_ERROR = 0x100,
        NORMALIZER_STANDARDIZATION_ERROR = 0x200,
        FRAGMENT_STANDARDIZATION_ERROR = 0x400,
        CHARGE_STANDARDIZATION_ERROR = 0x800,
        STANDARDIZATION_ERROR = 0xf00,
        OUTPUT_ERROR = 0x1000,
        PIPELINE_ERROR = 0x1fff,
        METALS_DISCONNECTED = 0x80_0000,
        NORMALIZATION_APPLIED = 0x100_0000,
        FRAGMENTS_REMOVED = 0x200_0000,
        PROTONATION_CHANGED = 0x400_0000,
        STRUCTURE_MODIFICATION = 0x780_0000,
    }

    /// `RDKit::MolStandardize::PipelineStage`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum PipelineStage {
        NOT_STARTED = 0,
        PARSING_INPUT,
        PREPARE_FOR_VALIDATION,
        VALIDATION,
        PREPARE_FOR_STANDARDIZATION,
        STANDARDIZATION,
        REAPPLY_WEDGING,
        CLEANUP_2D,
        MAKE_PARENT,
        SERIALIZING_OUTPUT,
        COMPLETED,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/standardize_pipeline.h");

        type Pipeline;
        type PipelineStatus;
        type PipelineStage;

        fn pipeline_options_default() -> Result<PipelineOptionsData>;
        fn pipeline_new(options: &PipelineOptionsData) -> Result<UniquePtr<Pipeline>>;
        fn pipeline_run(pipeline: &Pipeline, mol_data: &str) -> Result<PipelineResultData>;
    }
}

// SAFETY: a Pipeline holds its options and its steps (function pointers),
// which run only reads; each run makes its own molecules and objects, and
// the wrapper's serialization step locks 2D depiction's globals.
unsafe impl Send for ffi::Pipeline {}
unsafe impl Sync for ffi::Pipeline {}
