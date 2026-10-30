// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm2D.Generate`, a Python module: 2D
//! pharmacophore fingerprints.
//!
//! An N-point pharmacophore is N features and the topological distances
//! between them; an N-point proto-pharmacophore is N feature families
//! without distances; and an N-point scaffold is the binned distances
//! without the families.

use std::collections::BTreeMap;

use crate::chem::pharm_2d::sig_factory::{SigFactory, Signature};
use crate::chem::pharm_2d::utils::{get_index_combinations, get_unique_combinations, point_dists};
use crate::chem::{GetDistanceMatrixParams, Mol, get_distance_matrix_with_params};
use crate::{Error, Result};

/// The atoms of each feature of each pharmacophore that set a bit, by bit:
/// Python's `bitInfo`.
pub type BitInfo = BTreeMap<u64, Vec<Vec<Vec<u32>>>>;

/// Sets the bit of a match's pharmacophore from its features' shortest
/// distances, as Python's `_ShortestPathsMatch`; `None` when a distance
/// falls outside the bins.
fn shortest_paths_match(
    m: &[Vec<u32>],
    feature_set: &[usize],
    sig: &mut Signature,
    d_mat: &[Vec<f64>],
    factory: &SigFactory,
) -> Result<Option<u64>> {
    let dists_to_check = point_dists(m.len())?;
    let bins = factory.get_bins();
    let (Some(first), Some(last)) = (bins.first(), bins.last()) else {
        return Err(Error::index(0, 0));
    };
    let (min_d, max_d) = (f64::from(first.0), f64::from(last.1));
    let entry = |i: u32, j: u32| -> Result<f64> {
        let row = d_mat.get(i as usize).ok_or(Error::index(i, d_mat.len()))?;
        row.get(j as usize)
            .copied()
            .ok_or(Error::index(j, row.len()))
    };
    let mut dist = Vec::with_capacity(dists_to_check.len());
    for &(pt0, pt1) in dists_to_check {
        let mut min_seen = max_d;
        for &idx1 in &m[pt0] {
            for &idx2 in &m[pt1] {
                min_seen = min_seen.min(entry(idx1, idx2)?);
                if min_seen == 0.0 || min_seen < min_d {
                    return Ok(None);
                }
            }
        }
        // Python's int() truncates toward zero.
        let d = min_seen.trunc();
        if d == 0.0 || d < min_d || d >= max_d {
            return Ok(None);
        }
        dist.push(d as u32);
    }
    let idx = factory.get_bit_idx_with_sort_indices(feature_set, &dist, false)?;
    match sig {
        Signature::Bits(bits) => {
            bits.set_bit(idx as u32)?;
        }
        Signature::Counts(counts) => {
            let i = idx as i32;
            counts.set_val(i, counts.get_val(i)? + 1)?;
        }
        Signature::LongCounts(counts) => {
            let i = idx as i64;
            counts.set_val(i, counts.get_val(i)? + 1)?;
        }
    }
    Ok(Some(idx))
}

/// Parameters for [`gen_2d_fingerprint_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Gen2DFingerprintParams {
    /// The proto-pharmacophores to look for, as sorted feature indices;
    /// `None` for all of the factory's.
    pub perms: Option<Vec<Vec<usize>>>,
    /// The molecule's distance matrix; `None` to compute it.
    pub d_mat: Option<Vec<Vec<f64>>>,
}

/// The molecule's 2D pharmacophore fingerprint from an initialized
/// factory: a bit (or count) for each pharmacophore of its features, from
/// their shortest topological distances. Pharmacophores of more than 3
/// points are not supported, so a factory's larger `max_point_count` is
/// taken as 3 (Python logs a warning). A distance between the bins'
/// bounds that falls in no bin is `Error::Index`, as Python raises. Mirrors
/// `Generate.Gen2DFingerprint`.
pub fn gen_2d_fingerprint(mol: &Mol, sig_factory: &SigFactory) -> Result<Signature> {
    gen_2d_fingerprint_with_params(mol, sig_factory, &Gen2DFingerprintParams::default())
}

/// Mirrors `Generate.Gen2DFingerprint(mol, sigFactory, perms, dMat)`.
pub fn gen_2d_fingerprint_with_params(
    mol: &Mol,
    sig_factory: &SigFactory,
    params: &Gen2DFingerprintParams,
) -> Result<Signature> {
    Ok(gen_2d_fingerprint_with_outputs(mol, sig_factory, params)?.0)
}

/// The fingerprint, and the features' atoms in each match that set each
/// bit (A9). Mirrors `Generate.Gen2DFingerprint(..., bitInfo)`.
pub fn gen_2d_fingerprint_with_outputs(
    mol: &Mol,
    sig_factory: &SigFactory,
    params: &Gen2DFingerprintParams,
) -> Result<(Signature, BitInfo)> {
    let n_feats = sig_factory.get_feat_families().len();
    let min_count = sig_factory.min_point_count;
    let max_count = sig_factory.max_point_count.min(3);
    let computed;
    let d_mat = match &params.d_mat {
        Some(d_mat) => d_mat,
        None => {
            let dm_params = GetDistanceMatrixParams {
                use_bo: sig_factory.include_bond_order,
                use_atom_wts: false,
            };
            computed = get_distance_matrix_with_params(mol, &dm_params)?;
            &computed
        }
    };
    let perms = match &params.perms {
        Some(perms) => perms.clone(),
        None => (min_count..=max_count)
            .flat_map(|count| get_index_combinations(n_feats, count))
            .collect(),
    };
    let feat_matches = sig_factory.get_mol_feats(mol)?;
    let mut sig = sig_factory.get_signature()?;
    let mut bit_info = BitInfo::new();
    for perm in &perms {
        let mut feat_classes = vec![0; perm.len()];
        for i in 1..perm.len() {
            feat_classes[i] = if perm[i] == perm[i - 1] {
                feat_classes[i - 1]
            } else {
                feat_classes[i - 1] + 1
            };
        }
        let match_perms = perm
            .iter()
            .map(|&x| {
                feat_matches
                    .get(x)
                    .cloned()
                    .ok_or(Error::index(x, feat_matches.len()))
            })
            .collect::<Result<Vec<_>>>()?;
        let matches_to_map = get_unique_combinations(&match_perms, &feat_classes)?;
        for entry in matches_to_map {
            let m: Vec<Vec<u32>> = entry.into_iter().map(|(_, atoms)| atoms).collect();
            if sig_factory.shortest_paths_only
                && let Some(idx) = shortest_paths_match(&m, perm, &mut sig, d_mat, sig_factory)?
            {
                bit_info.entry(idx).or_default().push(m);
            }
        }
    }
    Ok((sig, bit_info))
}
