// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.AtomPairs.Torsions`, a Python module:
//! topological-torsion fingerprints (Nilakantan, Bauman, Dixon, and
//! Venkataraghavan, JCICS 27, 82-85 (1987)).
//!
//! Python's `GetTopologicalTorsionFingerprint`,
//! `GetHashedTopologicalTorsionFingerprint`, and
//! `GetTopologicalTorsionFingerprintAsIntVect` are rdMolDescriptors'
//! deprecated functions, which `rdFingerprintGenerator`'s torsion generator
//! replaces ([`get_topological_torsion_generator`](crate::chem::rd_fingerprint_generator::get_topological_torsion_generator)); its unfolded count
//! fingerprint has the same codes.

use crate::chem::Mol;
use crate::chem::atom_pairs::utils::{AtomCodeExplanation, explain_atom_code};
use crate::chem::rd_fingerprint_generator::{
    GetTopologicalTorsionGeneratorParams, get_topological_torsion_generator_with_params,
};
use crate::chem::rd_mol_descriptors::{
    AtomPairsParameters, GetAtomPairAtomCodeParams, get_atom_pair_atom_code_with_params,
};
use crate::{Error, Result};

/// The number subtracted from the degree of a path's `i`th atom of `size`.
fn branch_subtract(i: usize, size: usize) -> u32 {
    if i == 0 || i + 1 == size { 1 } else { 2 }
}

/// The code of a path of `size` atoms (`path`'s first `size`), from their
/// atom codes, each less its bonds in the path, in the direction that puts
/// the smaller code first; the same either way along the path. A path of
/// more than 7 atoms is `Error::Value`, since its code would not fit in 64
/// bits (Python's integers grow), as is a `path` shorter than `size`.
/// Mirrors `Torsions.pyScorePath`.
pub fn py_score_path(mol: &Mol, path: &[u32], size: usize) -> Result<u64> {
    py_score_path_with_atom_codes(mol, path, size, None)
}

/// With each atom's code given by index in `atom_codes` (without the
/// subtraction, which this does), instead of computed. Mirrors
/// `Torsions.pyScorePath(mol, path, size, atomCodes)`.
pub fn py_score_path_with_atom_codes(
    mol: &Mol,
    path: &[u32],
    size: usize,
    atom_codes: Option<&[u32]>,
) -> Result<u64> {
    let code_size = AtomPairsParameters::code_size() as usize;
    if size * code_size > 64 {
        return Err(Error::Value {
            message: format!("a path of {size} atoms has a code wider than 64 bits"),
        });
    }
    if path.len() < size {
        return Err(Error::Value {
            message: format!("the path has {} atoms, fewer than {size}", path.len()),
        });
    }
    // Python's `if not atomCodes` also computes the codes for an empty list.
    let atom_codes = atom_codes.filter(|codes| !codes.is_empty());
    let mut codes = Vec::with_capacity(size);
    for (i, &idx) in path[..size].iter().enumerate() {
        let sub = branch_subtract(i, size);
        let code = match atom_codes {
            Some(atom_codes) => {
                let code = *atom_codes.get(idx as usize).ok_or(Error::Index {
                    index: idx,
                    len: atom_codes.len() as u32,
                })?;
                // Python's integers go negative here; RDKit's codes do not.
                code.checked_sub(sub).ok_or_else(|| Error::Value {
                    message: format!("atom code {code} is smaller than {sub}"),
                })?
            }
            None => get_atom_pair_atom_code_with_params(
                &mol.get_atom_with_idx(idx)?,
                &GetAtomPairAtomCodeParams {
                    branch_subtract: sub,
                    include_chirality: false,
                },
            )?,
        };
        codes.push(code);
    }
    // "Canonicalize" the codes, as Python does.
    let (mut beg, mut end) = (0, codes.len().saturating_sub(1));
    while beg < end {
        if codes[beg] == codes[end] {
            beg += 1;
            end -= 1;
        } else {
            if codes[beg] > codes[end] {
                codes.reverse();
            }
            break;
        }
    }
    let mut accum = 0u64;
    for (i, &code) in codes.iter().enumerate() {
        accum |= u64::from(code)
            .checked_shl((code_size * i) as u32)
            .unwrap_or_default();
    }
    Ok(accum)
}

/// What a path code means: each atom's code explained, with its branches
/// counting its bonds in the path. Mirrors `Torsions.ExplainPathScore`.
pub fn explain_path_score(score: u64) -> Vec<AtomCodeExplanation> {
    explain_path_score_with_size(score, 4)
}

/// For a path of `size` atoms. Mirrors `Torsions.ExplainPathScore(score,
/// size)`.
pub fn explain_path_score_with_size(score: u64, size: usize) -> Vec<AtomCodeExplanation> {
    let code_size = AtomPairsParameters::code_size();
    let code_mask = (1u64 << code_size) - 1;
    let mut score = score;
    (0..size)
        .map(|i| {
            let code = score & code_mask;
            score = score.checked_shr(code_size).unwrap_or_default();
            let (symbol, n_branch, n_pi) = explain_atom_code(code);
            (symbol, n_branch + branch_subtract(i, size), n_pi)
        })
        .collect()
}

/// The torsion codes of the molecule's unfolded topological-torsion
/// fingerprint, sorted, each repeated by its count. Mirrors
/// `Torsions.GetTopologicalTorsionFingerprintAsIds`.
pub fn get_topological_torsion_fingerprint_as_ids(mol: &Mol) -> Result<Vec<u64>> {
    get_topological_torsion_fingerprint_as_ids_with_target_size(mol, 4)
}

/// For torsions of `target_size` atoms. Mirrors
/// `Torsions.GetTopologicalTorsionFingerprintAsIds(mol, targetSize)`.
pub fn get_topological_torsion_fingerprint_as_ids_with_target_size(
    mol: &Mol,
    target_size: u32,
) -> Result<Vec<u64>> {
    let generator =
        get_topological_torsion_generator_with_params(&GetTopologicalTorsionGeneratorParams {
            torsion_atom_count: target_size,
            ..Default::default()
        })?;
    let fp = generator.get_sparse_count_fingerprint(mol)?;
    let mut res = Vec::new();
    for (code, count) in fp.get_nonzero_elements() {
        res.extend(std::iter::repeat_n(code, count.max(0) as usize));
    }
    Ok(res)
}
