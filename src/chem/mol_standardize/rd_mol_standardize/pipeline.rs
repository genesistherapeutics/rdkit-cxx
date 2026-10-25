// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The validation and standardization pipeline.

use std::fmt;

use cxx::UniquePtr;

use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::expect_ok;
use crate::ffi::standardize_pipeline::ffi;
use crate::{Error, Result};

/// What happened in a pipeline run, as bits: errors, by step, and changes
/// to the structure. Mirrors `rdMolStandardize.PipelineStatus`
/// (CONVENTIONS.md F1).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct PipelineStatus(u32);

macro_rules! pipeline_status {
    ($($name:ident),* $(,)?) => {
        impl PipelineStatus {
            $(
                #[doc = concat!("Mirrors `PipelineStatus.", stringify!($name), "`.")]
                pub const $name: PipelineStatus = PipelineStatus(ffi::PipelineStatus::$name.repr);
            )*

            const NAMES: &[(&str, PipelineStatus)] =
                &[$((stringify!($name), PipelineStatus::$name)),*];
        }
    };
}

pipeline_status!(
    NO_EVENT,
    INPUT_ERROR,
    PREPARE_FOR_VALIDATION_ERROR,
    FEATURES_VALIDATION_ERROR,
    BASIC_VALIDATION_ERROR,
    IS2D_VALIDATION_ERROR,
    LAYOUT2D_VALIDATION_ERROR,
    STEREO_VALIDATION_ERROR,
    VALIDATION_ERROR,
    PREPARE_FOR_STANDARDIZATION_ERROR,
    METAL_STANDARDIZATION_ERROR,
    NORMALIZER_STANDARDIZATION_ERROR,
    FRAGMENT_STANDARDIZATION_ERROR,
    CHARGE_STANDARDIZATION_ERROR,
    STANDARDIZATION_ERROR,
    OUTPUT_ERROR,
    PIPELINE_ERROR,
    METALS_DISCONNECTED,
    NORMALIZATION_APPLIED,
    FRAGMENTS_REMOVED,
    PROTONATION_CHANGED,
    STRUCTURE_MODIFICATION,
);

crate::flags::flag_ops!(PipelineStatus);

rdkit_enum! {
    /// The step a pipeline run reached. Mirrors
    /// `rdMolStandardize.PipelineStage`, which lacks the C++ stages
    /// `NotStarted`, `ReapplyWedging`, `Cleanup2d`, and `MakeParent`.
    PipelineStage: ffi::PipelineStage {
        NotStarted = NOT_STARTED,
        ParsingInput = PARSING_INPUT,
        PrepareForValidation = PREPARE_FOR_VALIDATION,
        Validation = VALIDATION,
        PrepareForStandardization = PREPARE_FOR_STANDARDIZATION,
        Standardization = STANDARDIZATION,
        ReapplyWedging = REAPPLY_WEDGING,
        Cleanup2d = CLEANUP_2D,
        MakeParent = MAKE_PARENT,
        SerializingOutput = SERIALIZING_OUTPUT,
        Completed = COMPLETED,
    }
}

/// Options for a pipeline. Mirrors `rdMolStandardize.PipelineOptions`.
#[derive(Clone, Debug, PartialEq)]
pub struct PipelineOptions {
    /// Parses the mol block strictly.
    pub strict_parsing: bool,
    /// Goes on after a failure, to report them all.
    pub report_all_failures: bool,
    /// Allows molecules with no atoms.
    pub allow_empty_molecules: bool,
    /// Allows enhanced stereo.
    pub allow_enhanced_stereo: bool,
    /// Allows bonds of aromatic type.
    pub allow_aromatic_bond_type: bool,
    /// Allows dative bonds.
    pub allow_dative_bond_type: bool,
    /// `Is2DValidation`'s threshold.
    pub is_2d_zero_threshold: f64,
    /// `Layout2DValidation`'s clash limit.
    pub atom_clash_limit: f64,
    /// `Layout2DValidation`'s shortest median bond length.
    pub min_median_bond_length: f64,
    /// `Layout2DValidation`'s bond length limit.
    pub bond_length_limit: f64,
    /// Allows long bonds in rings.
    pub allow_long_bonds_in_rings: bool,
    /// Exempts long ring bonds from the atom-bond clash check.
    pub allow_atom_bond_clash_exemption: bool,
    /// The SMARTS of the metals to disconnect from nitrogen, oxygen, and
    /// fluorine, and the atom bonded to each: two atoms joined by a bond,
    /// or empty.
    pub metal_nof: String,
    /// The SMARTS of the metals to disconnect from other non-metals, as
    /// `metal_nof`.
    pub metal_non: String,
    /// The normalizations, as `normalizer_from_data` reads them; empty for
    /// RDKit's built-in ones.
    pub normalizer_data: String,
    /// The most times to restart the normalizations.
    pub normalizer_max_restarts: u32,
    /// The median bond length the output is scaled to; 0 or less leaves it.
    pub scaled_median_bond_length: f64,
    /// Writes V2000 mol blocks rather than V3000.
    pub output_v2000: bool,
}

/// Mirrors `PipelineOptions()`.
impl Default for PipelineOptions {
    fn default() -> PipelineOptions {
        let data = expect_ok(ffi::pipeline_options_default());
        let text = |bytes: Vec<u8>| String::from_utf8(bytes).expect("RDKit's options are UTF-8");
        PipelineOptions {
            strict_parsing: data.strict_parsing,
            report_all_failures: data.report_all_failures,
            allow_empty_molecules: data.allow_empty_molecules,
            allow_enhanced_stereo: data.allow_enhanced_stereo,
            allow_aromatic_bond_type: data.allow_aromatic_bond_type,
            allow_dative_bond_type: data.allow_dative_bond_type,
            is_2d_zero_threshold: data.is_2d_zero_threshold,
            atom_clash_limit: data.atom_clash_limit,
            min_median_bond_length: data.min_median_bond_length,
            bond_length_limit: data.bond_length_limit,
            allow_long_bonds_in_rings: data.allow_long_bonds_in_rings,
            allow_atom_bond_clash_exemption: data.allow_atom_bond_clash_exemption,
            metal_nof: text(data.metal_nof),
            metal_non: text(data.metal_non),
            normalizer_data: text(data.normalizer_data),
            normalizer_max_restarts: data.normalizer_max_restarts,
            scaled_median_bond_length: data.scaled_median_bond_length,
            output_v2000: data.output_v2000,
        }
    }
}

impl PipelineOptions {
    fn to_ffi(&self) -> ffi::PipelineOptionsData {
        ffi::PipelineOptionsData {
            strict_parsing: self.strict_parsing,
            report_all_failures: self.report_all_failures,
            allow_empty_molecules: self.allow_empty_molecules,
            allow_enhanced_stereo: self.allow_enhanced_stereo,
            allow_aromatic_bond_type: self.allow_aromatic_bond_type,
            allow_dative_bond_type: self.allow_dative_bond_type,
            is_2d_zero_threshold: self.is_2d_zero_threshold,
            atom_clash_limit: self.atom_clash_limit,
            min_median_bond_length: self.min_median_bond_length,
            bond_length_limit: self.bond_length_limit,
            allow_long_bonds_in_rings: self.allow_long_bonds_in_rings,
            allow_atom_bond_clash_exemption: self.allow_atom_bond_clash_exemption,
            metal_nof: self.metal_nof.as_bytes().to_vec(),
            metal_non: self.metal_non.as_bytes().to_vec(),
            normalizer_data: self.normalizer_data.as_bytes().to_vec(),
            normalizer_max_restarts: self.normalizer_max_restarts,
            scaled_median_bond_length: self.scaled_median_bond_length,
            output_v2000: self.output_v2000,
        }
    }
}

/// A message a pipeline run logged. Mirrors
/// `rdMolStandardize.PipelineLogEntry`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PipelineLogEntry {
    /// What happened. Mirrors `PipelineLogEntry.status`.
    pub status: PipelineStatus,
    /// The message. Mirrors `PipelineLogEntry.detail`.
    pub detail: String,
}

/// What a pipeline run did. Mirrors `rdMolStandardize.PipelineResult`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PipelineResult {
    /// Every status the run logged, combined. Mirrors
    /// `PipelineResult.status`.
    pub status: PipelineStatus,
    /// The step the run reached. Mirrors `PipelineResult.stage`.
    pub stage: PipelineStage,
    /// The messages. Mirrors `PipelineResult.log`, a `PipelineLog`.
    pub log: Vec<PipelineLogEntry>,
    /// The input. Mirrors `PipelineResult.inputMolData`.
    pub input_mol_data: String,
    /// The standardized molecule, as a mol block. Mirrors
    /// `PipelineResult.outputMolData`.
    pub output_mol_data: String,
    /// Its parent (uncharged, without dative bonds), as a mol block.
    /// Mirrors `PipelineResult.parentMolData`.
    pub parent_mol_data: String,
}

fn utf8(bytes: Vec<u8>) -> Result<String> {
    String::from_utf8(bytes).map_err(|err| Error::Value {
        message: format!("invalid UTF-8 in a pipeline result: {err}"),
    })
}

impl PipelineResult {
    fn from_ffi(data: ffi::PipelineResultData) -> Result<PipelineResult> {
        Ok(PipelineResult {
            status: PipelineStatus(data.status),
            stage: PipelineStage::from_ffi(ffi::PipelineStage { repr: data.stage }),
            log: data
                .log
                .into_iter()
                .map(|entry| {
                    Ok(PipelineLogEntry {
                        status: PipelineStatus(entry.status),
                        detail: utf8(entry.detail)?,
                    })
                })
                .collect::<Result<_>>()?,
            input_mol_data: utf8(data.input_mol_data)?,
            output_mol_data: utf8(data.output_mol_data)?,
            parent_mol_data: utf8(data.parent_mol_data)?,
        })
    }
}

/// Validates and standardizes molecules read from mol blocks, as
/// Novartis's pipeline does: parses, validates a copy, disconnects metals,
/// normalizes, keeps the largest fragment, reapplies the input's wedging,
/// scales the 2D layout, makes the uncharged parent, and writes both.
/// Mirrors `rdMolStandardize.Pipeline`.
pub struct Pipeline {
    /// Never null.
    ptr: UniquePtr<ffi::Pipeline>,
}

impl fmt::Debug for Pipeline {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Pipeline").finish_non_exhaustive()
    }
}

/// Mirrors `Pipeline()`.
impl Default for Pipeline {
    fn default() -> Pipeline {
        Pipeline::new(&PipelineOptions::default()).expect("the default options are valid")
    }
}

impl Pipeline {
    /// A pipeline with these options. `Error::Value` for metal SMARTS that
    /// do not parse, which RDKit would dereference as null, or that are not
    /// two atoms joined by a bond (or empty), whose matches RDKit would read
    /// past (`MetalDisconnector::set_metal_nof`). Mirrors
    /// `Pipeline(options)`.
    pub fn new(options: &PipelineOptions) -> Result<Pipeline> {
        let ptr = ffi::pipeline_new(&options.to_ffi())?;
        assert!(!ptr.is_null(), "RDKit returned a null Pipeline");
        Ok(Pipeline { ptr })
    }

    /// Runs the pipeline on a mol block. Its problems are in the result's
    /// status and log; `Error::Value` for a mol block that breaks the
    /// reader's limits (CONVENTIONS.md U12), for a molecule whose charges
    /// could overflow RDKit's `int` sums (U9), and for text in the result
    /// that is not UTF-8, which Python raises on reading. Mirrors
    /// `Pipeline.run`.
    pub fn run(&self, mol_data: &str) -> Result<PipelineResult> {
        PipelineResult::from_ffi(ffi::pipeline_run(self.ptr.as_ref().unwrap(), mol_data)?)
    }
}
