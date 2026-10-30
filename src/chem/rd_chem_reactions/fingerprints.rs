// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Reaction fingerprints.

use crate::Result;
use crate::chem::rd_chem_reactions::ChemicalReaction;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::data_structs::{ExplicitBitVect, UIntSparseIntVect};
use crate::error::expect_ok;
use crate::ffi::chem_reactions::ffi;

rdkit_enum! {
    /// The molecular fingerprint a reaction fingerprint is built from.
    /// Mirrors `rdChemReactions.FingerprintType`.
    FingerprintType: ffi::FingerprintType {
        AtomPairFp = AtomPairFP,
        TopologicalTorsion = TopologicalTorsion,
        MorganFp = MorganFP,
        RdkitFp = RDKitFP,
        PatternFp = PatternFP,
    }
}

/// How a reaction fingerprint is computed. Mirrors
/// `rdChemReactions.ReactionFingerprintParams`.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionFingerprintParams {
    /// Includes the agents.
    pub include_agents: bool,
    /// The agents' share of a structural fingerprint's bits.
    pub bit_ratio_agents: f64,
    /// A difference fingerprint's weight for reactants and products, with
    /// agents.
    pub non_agent_weight: u32,
    /// A difference fingerprint's weight for agents.
    pub agent_weight: i32,
    /// The number of bits.
    pub fp_size: u32,
    /// The molecular fingerprint. Difference fingerprints take only atom
    /// pairs, topological torsions, and Morgan.
    pub fp_type: FingerprintType,
}

impl ReactionFingerprintParams {
    /// Mirrors `ReactionFingerprintParams(includeAgents, bitRatioAgents,
    /// nonAgentWeight, agentWeight, fpSize, fpType)`.
    pub fn new(
        include_agents: bool,
        bit_ratio_agents: f64,
        non_agent_weight: u32,
        agent_weight: i32,
        fp_size: u32,
        fp_type: FingerprintType,
    ) -> ReactionFingerprintParams {
        ReactionFingerprintParams {
            include_agents,
            bit_ratio_agents,
            non_agent_weight,
            agent_weight,
            fp_size,
            fp_type,
        }
    }

    fn from_ffi(data: ffi::ReactionFingerprintParamsData) -> ReactionFingerprintParams {
        ReactionFingerprintParams {
            include_agents: data.include_agents,
            bit_ratio_agents: data.bit_ratio_agents,
            non_agent_weight: data.non_agent_weight,
            agent_weight: data.agent_weight,
            fp_size: data.fp_size,
            fp_type: FingerprintType::from_ffi(data.fp_type),
        }
    }

    fn to_ffi(&self) -> ffi::ReactionFingerprintParamsData {
        ffi::ReactionFingerprintParamsData {
            include_agents: self.include_agents,
            bit_ratio_agents: self.bit_ratio_agents,
            non_agent_weight: self.non_agent_weight,
            agent_weight: self.agent_weight,
            fp_size: self.fp_size,
            fp_type: self.fp_type.to_ffi(),
        }
    }
}

/// RDKit's defaults. Mirrors `ReactionFingerprintParams()`.
impl Default for ReactionFingerprintParams {
    fn default() -> ReactionFingerprintParams {
        ReactionFingerprintParams::from_ffi(expect_ok(ffi::reaction_fingerprint_params_default()))
    }
}

/// The products' count fingerprint less the reactants', 2048 atom-pair
/// bits with agents weighted in. Mirrors
/// `rdChemReactions.CreateDifferenceFingerprintForReaction`.
pub fn create_difference_fingerprint_for_reaction(
    reaction: &ChemicalReaction,
) -> Result<UIntSparseIntVect> {
    let params =
        ReactionFingerprintParams::from_ffi(expect_ok(ffi::default_difference_fp_params()));
    create_difference_fingerprint_for_reaction_with_reaction_finger_print_params(reaction, &params)
}

/// RDKit counts in `i32`, so counts that would overflow one after
/// weighting, or a `non_agent_weight` above `i32::MAX`, are `Error::Value`;
/// a zero `fp_size` or a fingerprint type other than atom pairs, torsions,
/// or Morgan is `Error::Runtime`. Mirrors
/// `rdChemReactions.CreateDifferenceFingerprintForReaction(reaction,
/// ReactionFingerPrintParams)`.
pub fn create_difference_fingerprint_for_reaction_with_reaction_finger_print_params(
    reaction: &ChemicalReaction,
    reaction_finger_print_params: &ReactionFingerprintParams,
) -> Result<UIntSparseIntVect> {
    Ok(UIntSparseIntVect::from_ffi(
        ffi::create_difference_fingerprint_for_reaction(
            reaction.as_ffi(),
            &reaction_finger_print_params.to_ffi(),
        )?,
    ))
}

/// The reactants' fingerprint, then the products', then the agents', 4096
/// pattern-fingerprint bits in all. Mirrors
/// `rdChemReactions.CreateStructuralFingerprintForReaction`.
pub fn create_structural_fingerprint_for_reaction(
    reaction: &ChemicalReaction,
) -> Result<ExplicitBitVect> {
    let params =
        ReactionFingerprintParams::from_ffi(expect_ok(ffi::default_structural_fp_params()));
    create_structural_fingerprint_for_reaction_with_reaction_finger_print_params(reaction, &params)
}

/// Parameters whose parts would hold more than `i32::MAX` bits, or whose
/// `fp_size * bit_ratio_agents` does not fit an `i32`, are `Error::Value`; a
/// zero `fp_size` is `Error::Runtime`. Mirrors
/// `rdChemReactions.CreateStructuralFingerprintForReaction(reaction,
/// ReactionFingerPrintParams)`.
pub fn create_structural_fingerprint_for_reaction_with_reaction_finger_print_params(
    reaction: &ChemicalReaction,
    reaction_finger_print_params: &ReactionFingerprintParams,
) -> Result<ExplicitBitVect> {
    Ok(ExplicitBitVect::from_ffi(
        ffi::create_structural_fingerprint_for_reaction(
            reaction.as_ffi(),
            &reaction_finger_print_params.to_ffi(),
        )?,
    ))
}
