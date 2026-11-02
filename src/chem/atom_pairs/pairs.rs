// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.AtomPairs.Pairs`, a Python module: atom-pair
//! fingerprints (Carhart, Smith, and Venkataraghavan, JCICS 25, 64-73
//! (1985)).
//!
//! Python's `GetAtomPairFingerprint`, `GetHashedAtomPairFingerprint`, and
//! `GetAtomPairFingerprintAsIntVect` are rdMolDescriptors' deprecated
//! functions, which `rdFingerprintGenerator`'s atom-pair generator
//! replaces ([`get_atom_pair_generator`]); its unfolded count fingerprint
//! has the same codes.

use crate::Result;
use crate::chem::Mol;
use crate::chem::atom_pairs::utils::{
    AtomCodeExplanation, ChiralAtomCodeExplanation, explain_atom_code,
    explain_atom_code_include_chirality,
};
use crate::chem::rd_fingerprint_generator::get_atom_pair_generator;
use crate::chem::rd_mol_descriptors::{
    AtomPairsParameters, GetAtomPairAtomCodeParams, get_atom_pair_atom_code_with_params,
};
use crate::chem::rdchem::AsAtom;
use crate::data_structs::SparseBitVect;

/// The bits for a pair's distance. Mirrors `Pairs.numPathBits`.
pub fn num_path_bits() -> u32 {
    AtomPairsParameters::num_path_bits()
}

/// The bits in a pair's code. Mirrors `Pairs.numFpBits`.
pub fn num_fp_bits() -> u32 {
    num_path_bits() + 2 * AtomPairsParameters::code_size()
}

/// The number of possible pair codes. Mirrors `Pairs.fpLen`.
pub fn fp_len() -> u32 {
    1 << num_fp_bits()
}

/// Parameters for [`py_score_pair_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PyScorePairParams {
    /// The two atoms' codes, used instead of computing them.
    pub atom_codes: Option<(u32, u32)>,
    /// Includes the atoms' CIP labels in their codes.
    pub include_chirality: bool,
}

/// The code of two atoms `dist` bonds apart, from their atom codes, the
/// smaller first; the distance is taken modulo 31. The atoms must be in a
/// molecule. Mirrors `Pairs.pyScorePair`.
pub fn py_score_pair(at1: &impl AsAtom, at2: &impl AsAtom, dist: u32) -> Result<u64> {
    py_score_pair_with_params(at1, at2, dist, &PyScorePairParams::default())
}

/// Mirrors `Pairs.pyScorePair(at1, at2, dist, atomCodes, includeChirality)`.
pub fn py_score_pair_with_params(
    at1: &impl AsAtom,
    at2: &impl AsAtom,
    dist: u32,
    params: &PyScorePairParams,
) -> Result<u64> {
    let (code1, code2) = match params.atom_codes {
        Some(codes) => codes,
        None => {
            let code_params = GetAtomPairAtomCodeParams {
                branch_subtract: 0,
                include_chirality: params.include_chirality,
            };
            (
                get_atom_pair_atom_code_with_params(at1, &code_params)?,
                get_atom_pair_atom_code_with_params(at2, &code_params)?,
            )
        }
    };
    let mut code_size = AtomPairsParameters::code_size();
    if params.include_chirality {
        code_size += AtomPairsParameters::num_chiral_bits();
    }
    let max_path_len = (1u64 << num_path_bits()) - 1;
    let mut accum = u64::from(dist) % max_path_len;
    accum |= u64::from(code1.min(code2)) << num_path_bits();
    accum |= u64::from(code1.max(code2)) << (code_size + num_path_bits());
    Ok(accum)
}

fn split_pair_score(score: u64, include_chirality: bool) -> (u64, u32, u64) {
    let mut code_size = AtomPairsParameters::code_size();
    if include_chirality {
        code_size += AtomPairsParameters::num_chiral_bits();
    }
    let code_mask = (1u64 << code_size) - 1;
    let path_mask = (1u64 << num_path_bits()) - 1;
    let dist = (score & path_mask) as u32;
    let score = score >> num_path_bits();
    (score & code_mask, dist, (score >> code_size) & code_mask)
}

/// What a pair code means: its two atoms' codes explained, and their
/// distance. Mirrors `Pairs.ExplainPairScore`.
pub fn explain_pair_score(score: u64) -> (AtomCodeExplanation, u32, AtomCodeExplanation) {
    let (code1, dist, code2) = split_pair_score(score, false);
    (explain_atom_code(code1), dist, explain_atom_code(code2))
}

/// For a code made with `includeChirality`, with the atoms' CIP labels;
/// `Error::Key` for a label code of 3. Mirrors
/// `Pairs.ExplainPairScore(score, includeChirality=True)` (A4).
pub fn explain_pair_score_include_chirality(
    score: u64,
) -> Result<(ChiralAtomCodeExplanation, u32, ChiralAtomCodeExplanation)> {
    let (code1, dist, code2) = split_pair_score(score, true);
    Ok((
        explain_atom_code_include_chirality(code1)?,
        dist,
        explain_atom_code_include_chirality(code2)?,
    ))
}

/// The atom-pair fingerprint as a [`SparseBitVect`] of [`fp_len`] bits, one
/// for each pair code the molecule has (not their counts, as the standard
/// definition has). Mirrors `Pairs.GetAtomPairFingerprintAsBitVect`.
pub fn get_atom_pair_fingerprint_as_bit_vect(mol: &Mol) -> Result<SparseBitVect> {
    let fp = get_atom_pair_generator().get_sparse_count_fingerprint(mol)?;
    let mut res = SparseBitVect::new(fp_len());
    for &code in fp.get_nonzero_elements().keys() {
        let code = u32::try_from(code).expect("atom-pair codes are below fp_len");
        res.set_bit(code)?;
    }
    Ok(res)
}
