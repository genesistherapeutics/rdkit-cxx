// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm2D.Matcher`, a pure-Python module: the atoms
//! that match a 2D pharmacophore's bit.
//!
//! Python's `MatchError` is an exception class nothing raises, so it is
//! not ported.

use crate::chem::pharm_2d::sig_factory::SigFactory;
use crate::chem::pharm_2d::utils::{get_all_combinations, point_dists};
use crate::chem::{GetDistanceMatrixParams, Mol, get_distance_matrix_with_params};
use crate::{Error, Result};

/// Parameters for [`get_atoms_matching_bit_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetAtomsMatchingBitParams {
    /// The molecule's distance matrix; `None` to compute it.
    pub d_mat: Option<Vec<Vec<f64>>>,
    /// Returns only the first match.
    pub just_one: bool,
    /// The atoms of each of the molecule's features, by family, as
    /// [`SigFactory::get_mol_feats`] gives them; `None` to find them.
    pub matching_atoms: Option<Vec<Vec<Vec<u32>>>>,
}

/// The features' atoms in each pharmacophore of the molecule that matches
/// the bit, each sorted, with a pharmacophore once. The distances are
/// those of each feature's first atom. A factory not limited to shortest
/// paths is `Error::Value` (Python's `AssertionError`). Mirrors
/// `Matcher.GetAtomsMatchingBit`.
pub fn get_atoms_matching_bit(
    sig_factory: &SigFactory,
    bit_idx: u64,
    mol: &Mol,
) -> Result<Vec<Vec<Vec<u32>>>> {
    get_atoms_matching_bit_with_params(
        sig_factory,
        bit_idx,
        mol,
        &GetAtomsMatchingBitParams::default(),
    )
}

/// Mirrors `Matcher.GetAtomsMatchingBit(sigFactory, bitIdx, mol, dMat,
/// justOne, matchingAtoms)`.
pub fn get_atoms_matching_bit_with_params(
    sig_factory: &SigFactory,
    bit_idx: u64,
    mol: &Mol,
    params: &GetAtomsMatchingBitParams,
) -> Result<Vec<Vec<Vec<u32>>>> {
    if !sig_factory.shortest_paths_only {
        return Err(Error::Value {
            message: "not implemented for non-shortest path signatures".to_owned(),
        });
    }
    let (n_pts, feat_combo, scaffold) = sig_factory.get_bit_info(bit_idx)?;
    let found;
    let matching_atoms = match &params.matching_atoms {
        Some(atoms) => atoms,
        None => {
            found = sig_factory.get_mol_feats(mol)?;
            &found
        }
    };
    let mut choices = Vec::new();
    for &feat_idx in &feat_combo {
        let tmp = matching_atoms.get(feat_idx).ok_or(Error::Index {
            index: feat_idx as u32,
            len: matching_atoms.len() as u32,
        })?;
        if tmp.is_empty() {
            return Ok(Vec::new());
        }
        choices.push(tmp.clone());
    }
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
    let entry = |i: u32, j: u32| -> Result<f64> {
        let row = d_mat.get(i as usize).ok_or(Error::Index {
            index: i,
            len: d_mat.len() as u32,
        })?;
        row.get(j as usize).copied().ok_or(Error::Index {
            index: j,
            len: row.len() as u32,
        })
    };
    let dists_to_check = point_dists(n_pts)?;
    let bins = sig_factory.get_bins();
    let mut res: Vec<Vec<Vec<u32>>> = Vec::new();
    for mut proto_pharm in get_all_combinations(&choices) {
        let mut passes = true;
        for (i, &(a1, a2)) in dists_to_check.iter().enumerate() {
            let (d_low, d_high) = *bins.get(scaffold[i]).ok_or(Error::Index {
                index: scaffold[i] as u32,
                len: bins.len() as u32,
            })?;
            let first = |atoms: &Vec<u32>| -> Result<u32> {
                atoms
                    .first()
                    .copied()
                    .ok_or(Error::Index { index: 0, len: 0 })
            };
            let dist = entry(first(&proto_pharm[a1])?, first(&proto_pharm[a2])?)?;
            if dist < f64::from(d_low) || dist >= f64::from(d_high) {
                passes = false;
                break;
            }
        }
        if passes {
            proto_pharm.sort();
            if !res.contains(&proto_pharm) {
                res.push(proto_pharm);
                if params.just_one {
                    break;
                }
            }
        }
    }
    Ok(res)
}
