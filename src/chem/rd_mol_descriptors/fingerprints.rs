// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's MACCS keys, atom-pair codes, and Morgan invariants.
//! RDKit deprecates the module's atom-pair, torsion, and Morgan fingerprint
//! functions in favor of `rd_fingerprint_generator`, so they are not bound
//! (CONVENTIONS.md I2).

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::AsAtom;
use crate::data_structs::ExplicitBitVect;
use crate::error::{expect_ok, to_string};
use crate::ffi::descriptors::ffi;

/// The 166 MACCS keys, as bits 1 to 166 of 167. Mirrors
/// `rdMolDescriptors.GetMACCSKeysFingerprint`.
pub fn get_maccs_keys_fingerprint(mol: &Mol) -> Result<ExplicitBitVect> {
    Ok(ExplicitBitVect::from_ffi(ffi::get_maccs_keys_fingerprint(
        mol.as_ffi(),
    )?))
}

/// Parameters for `get_atom_pair_atom_code_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetAtomPairAtomCodeParams {
    /// Subtracted from the atom's degree, as the torsion codes do.
    pub branch_subtract: u32,
    /// Adds two bits for the atom's CIP label, which RDKit assigns to the
    /// molecule first when it has none.
    pub include_chirality: bool,
}

/// An atom's code in the atom-pair scheme, from its degree, pi electrons,
/// and element. The atom must be in a molecule: an owned `Atom` is
/// `Error::Runtime`. Mirrors `rdMolDescriptors.GetAtomPairAtomCode`.
pub fn get_atom_pair_atom_code(atom: &impl AsAtom) -> Result<u32> {
    get_atom_pair_atom_code_with_params(atom, &GetAtomPairAtomCodeParams::default())
}

/// Mirrors `rdMolDescriptors.GetAtomPairAtomCode(atom, branchSubtract,
/// includeChirality)`.
pub fn get_atom_pair_atom_code_with_params(
    atom: &impl AsAtom,
    params: &GetAtomPairAtomCodeParams,
) -> Result<u32> {
    Ok(ffi::get_atom_pair_atom_code(
        atom.atom(),
        params.branch_subtract,
        params.include_chirality,
    )?)
}

/// The code for two atoms' codes `distance` bonds apart; `Error::Runtime`
/// for a distance above 30. Mirrors `rdMolDescriptors.GetAtomPairCode`.
pub fn get_atom_pair_code(atom1_code: u32, atom2_code: u32, distance: u32) -> Result<u32> {
    get_atom_pair_code_with_include_chirality(atom1_code, atom2_code, distance, false)
}

/// With room for the chirality bits of `include_chirality` codes. Mirrors
/// `rdMolDescriptors.GetAtomPairCode(atom1Code, atom2Code, distance,
/// includeChirality)`.
pub fn get_atom_pair_code_with_include_chirality(
    atom1_code: u32,
    atom2_code: u32,
    distance: u32,
    include_chirality: bool,
) -> Result<u32> {
    Ok(ffi::get_atom_pair_code(
        atom1_code,
        atom2_code,
        distance,
        include_chirality,
    )?)
}

/// Each atom's ECFP-like invariant: a hash of its element, degree,
/// hydrogens, charge, mass difference, and ring membership. Mirrors
/// `rdMolDescriptors.GetConnectivityInvariants`.
pub fn get_connectivity_invariants(mol: &Mol) -> Result<Vec<u32>> {
    get_connectivity_invariants_with_include_ring_membership(mol, true)
}

/// Without ring membership unless `include_ring_membership`. Mirrors
/// `rdMolDescriptors.GetConnectivityInvariants(mol, includeRingMembership)`.
pub fn get_connectivity_invariants_with_include_ring_membership(
    mol: &Mol,
    include_ring_membership: bool,
) -> Result<Vec<u32>> {
    Ok(ffi::get_connectivity_invariants(
        mol.as_ffi(),
        include_ring_membership,
    )?)
}

/// Each atom's FCFP-like invariant: a bit for each pharmacophoric feature it
/// is in. Mirrors `rdMolDescriptors.GetFeatureInvariants`.
pub fn get_feature_invariants(mol: &Mol) -> Result<Vec<u32>> {
    Ok(ffi::get_feature_invariants(mol.as_ffi())?)
}

/// The constants of RDKit's atom-pair codes, from `AtomPairs.h`. Mirrors
/// `rdMolDescriptors.AtomPairsParameters`, whose properties belong to the
/// class; like it, this has no values.
#[derive(Debug)]
pub enum AtomPairsParameters {}

fn constants() -> ffi::AtomPairsConstants {
    expect_ok(ffi::atom_pairs_constants())
}

impl AtomPairsParameters {
    /// The atom-pair code's version. Mirrors
    /// `AtomPairsParameters.version`.
    pub fn version() -> String {
        to_string(expect_ok(ffi::atom_pairs_version())).expect("RDKit's version is ASCII")
    }

    /// The bits for an atom's element. Mirrors
    /// `AtomPairsParameters.numTypeBits`.
    pub fn num_type_bits() -> u32 {
        constants().num_type_bits
    }

    /// The bits for an atom's pi electrons. Mirrors
    /// `AtomPairsParameters.numPiBits`.
    pub fn num_pi_bits() -> u32 {
        constants().num_pi_bits
    }

    /// The bits for an atom's degree. Mirrors
    /// `AtomPairsParameters.numBranchBits`.
    pub fn num_branch_bits() -> u32 {
        constants().num_branch_bits
    }

    /// The bits for an atom's CIP label. Mirrors
    /// `AtomPairsParameters.numChiralBits`.
    pub fn num_chiral_bits() -> u32 {
        constants().num_chiral_bits
    }

    /// The bits in an atom's code without chirality. Mirrors
    /// `AtomPairsParameters.codeSize`.
    pub fn code_size() -> u32 {
        constants().code_size
    }

    /// The elements the type bits encode, by type; other elements take the
    /// last type. Mirrors `AtomPairsParameters.atomTypes`.
    pub fn atom_types() -> Vec<u32> {
        constants().atom_types
    }

    /// The bits for a pair's distance. Mirrors
    /// `AtomPairsParameters.numPathBits`.
    pub fn num_path_bits() -> u32 {
        constants().num_path_bits
    }

    /// The bits in a pair's code without chirality. Mirrors
    /// `AtomPairsParameters.numAtomPairFingerprintBits`.
    pub fn num_atom_pair_fingerprint_bits() -> u32 {
        constants().num_atom_pair_fingerprint_bits
    }
}
