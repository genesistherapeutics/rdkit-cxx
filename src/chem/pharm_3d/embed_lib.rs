// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm3D.EmbedLib`, a Python module: matching 3D
//! pharmacophores to molecules and embedding molecules to satisfy them.
//!
//! Bounds matrices are rows, as `get_molecule_bounds_matrix` returns them:
//! upper bounds above the diagonal and lower bounds below. Python's verbose
//! output and timing are not ported, nor are its prints to standard error
//! of excluded volumes' positions.

use std::collections::BTreeSet;

use crate::chem::pharm_3d::excluded_volume::{ExcludedVolume, ExcludedVolumeParams};
use crate::chem::pharm_3d::pharmacophore::Pharmacophore;
use crate::chem::rd_chemical_features::FreeChemicalFeature;
use crate::chem::rd_dist_geom::get_molecule_bounds_matrix;
use crate::chem::rd_force_field_helpers::uff_get_molecule_force_field;
use crate::chem::rd_mol_chemical_features::{
    GetMolFeatureParams, MolChemicalFeature, MolChemicalFeatureFactory,
};
use crate::chem::rdchem::{AtomRef, Conformer, Mol};
use crate::chem::{
    GetDistanceMatrixParams, assign_stereochemistry, find_mol_chiral_centers,
    get_distance_matrix_with_params,
};
use crate::distance_geometry::{
    EmbedBoundsMatrixParams, do_triangle_smoothing, embed_bounds_matrix_with_params,
};
use crate::{Error, Result};

/// The length of feature directions in bounds matrices. Mirrors
/// `EmbedLib.defaultFeatLength`.
pub fn default_feat_length() -> f64 {
    2.0
}

fn value_error(message: &str) -> Error {
    Error::Value {
        message: message.to_owned(),
    }
}

fn index_error(index: usize, len: usize) -> Error {
    Error::Index {
        index: index as u32,
        len: len as u32,
    }
}

/// `bm[i][j]`, or `Error::Index` outside the matrix.
fn at(bm: &[Vec<f64>], i: usize, j: usize) -> Result<f64> {
    let row = bm.get(i).ok_or_else(|| index_error(i, bm.len()))?;
    row.get(j).copied().ok_or_else(|| index_error(j, row.len()))
}

fn set(bm: &mut [Vec<f64>], i: usize, j: usize, val: f64) -> Result<()> {
    let len = bm.len();
    let row = bm.get_mut(i).ok_or_else(|| index_error(i, len))?;
    let len = row.len();
    *row.get_mut(j).ok_or_else(|| index_error(j, len))? = val;
    Ok(())
}

/// The atom's neighbors other than hydrogens. Mirrors
/// `EmbedLib.GetAtomHeavyNeighbors`.
pub fn get_atom_heavy_neighbors<'m>(atom: &AtomRef<'m>) -> Vec<AtomRef<'m>> {
    atom.get_neighbors()
        .into_iter()
        .filter(|n| n.get_atomic_num() != 1)
        .collect()
}

/// Parameters for [`replace_group_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReplaceGroupParams {
    /// The fraction by which the bounds widen.
    pub slop: f64,
    /// Adds a point for the feature's direction too.
    pub use_dirs: bool,
    /// The direction's length.
    pub dir_length: f64,
}

impl Default for ReplaceGroupParams {
    fn default() -> ReplaceGroupParams {
        ReplaceGroupParams {
            slop: 0.01,
            use_dirs: false,
            dir_length: default_feat_length(),
        }
    }
}

/// A copy of the bounds matrix with a row and column added at the end for
/// the center of a multi-atom feature, whose atoms form a regular polygon
/// in order, and the new row's index. Mirrors `EmbedLib.ReplaceGroup`.
pub fn replace_group(match_: &[u32], bounds: &[Vec<f64>]) -> Result<(Vec<Vec<f64>>, usize)> {
    replace_group_with_params(match_, bounds, &ReplaceGroupParams::default())
}

/// Mirrors `EmbedLib.ReplaceGroup(match, bounds, slop, useDirs,
/// dirLength)`.
pub fn replace_group_with_params(
    match_: &[u32],
    bounds: &[Vec<f64>],
    params: &ReplaceGroupParams,
) -> Result<(Vec<Vec<f64>>, usize)> {
    let mut max_val = -1000.0f64;
    let mut min_val = 1e8f64;
    let n_pts = match_.len();
    for i in 0..n_pts {
        let mut idx0 = match_[i] as usize;
        let mut idx1 = match_[(i + 1) % n_pts] as usize;
        if idx1 < idx0 {
            std::mem::swap(&mut idx0, &mut idx1);
        }
        min_val = min_val.min(at(bounds, idx1, idx0)?);
        max_val = max_val.max(at(bounds, idx0, idx1)?);
    }
    let slop = params.slop;
    max_val *= 1.0 + slop;
    min_val *= 1.0 - slop;
    let scale_fact = 1.0 / (2.0 * (std::f64::consts::PI / n_pts as f64).sin());
    min_val *= scale_fact;
    max_val *= scale_fact;

    let replace_idx = bounds.len();
    let size = bounds.len() + 1 + usize::from(params.use_dirs);
    let mut bm = vec![vec![0.0; size]; size];
    for (i, row) in bounds.iter().enumerate() {
        bm[i][..row.len()].copy_from_slice(row);
    }
    for row in bm.iter_mut().take(replace_idx) {
        row[replace_idx] = 1000.0;
    }
    if params.use_dirs {
        for row in bm.iter_mut().take(replace_idx + 1) {
            row[replace_idx + 1] = 1000.0;
        }
        bm[replace_idx][replace_idx + 1] = params.dir_length + slop;
        bm[replace_idx + 1][replace_idx] = params.dir_length - slop;
    }
    for &idx1 in match_ {
        let idx1 = idx1 as usize;
        bm[idx1][replace_idx] = max_val;
        bm[replace_idx][idx1] = min_val;
        if params.use_dirs {
            let up = bm[replace_idx][replace_idx + 1];
            let low = bm[replace_idx + 1][replace_idx];
            bm[idx1][replace_idx + 1] = (up * up + max_val * max_val).sqrt();
            bm[replace_idx + 1][idx1] = (low * low + min_val * min_val).sqrt();
        }
    }
    Ok((bm, replace_idx))
}

/// Parameters for [`embed_mol_with_params`], with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct EmbedMolParams {
    /// The pharmacophore's atoms; its length sets which bounds get
    /// `weight` (Python weights the pairs of the first that many rows).
    pub atom_match: Option<Vec<Vec<u32>>>,
    /// The weight of the pharmacophore's and the volumes' bounds.
    pub weight: f64,
    /// The random seed; −1 for none.
    pub random_seed: i32,
    /// Excluded volumes in the bounds matrix, whose positions the
    /// embedding sets (an in-out argument, A9).
    pub excluded_volumes: Option<Vec<ExcludedVolume>>,
}

impl Default for EmbedMolParams {
    fn default() -> EmbedMolParams {
        EmbedMolParams {
            atom_match: None,
            weight: 2.0,
            random_seed: -1,
            excluded_volumes: None,
        }
    }
}

/// Embeds the molecule from a bounds matrix (which may have rows past its
/// atoms) and adds the coordinates as conformer 0. A matrix that does not
/// embed is `Error::Value`. Mirrors `EmbedLib.EmbedMol`.
pub fn embed_mol(mol: &mut Mol, bm: &[Vec<f64>]) -> Result<()> {
    embed_mol_with_params(mol, bm, &EmbedMolParams::default())
}

/// Mirrors `EmbedLib.EmbedMol(mol, bm, atomMatch, weight, randomSeed,
/// excludedVolumes)`.
pub fn embed_mol_with_params(
    mol: &mut Mol,
    bm: &[Vec<f64>],
    params: &EmbedMolParams,
) -> Result<()> {
    embed_mol_with_outputs(mol, bm, params).map(|_| ())
}

/// The excluded volumes with their positions in the embedding (A9).
/// Mirrors `EmbedLib.EmbedMol(..., excludedVolumes)`, which sets them.
pub fn embed_mol_with_outputs(
    mol: &mut Mol,
    bm: &[Vec<f64>],
    params: &EmbedMolParams,
) -> Result<Vec<ExcludedVolume>> {
    let n_ats = mol.get_num_atoms();
    let mut weights = Vec::new();
    if let Some(atom_match) = params.atom_match.as_ref().filter(|m| !m.is_empty()) {
        let size = atom_match.len() as u32;
        for i in 0..size {
            for j in i + 1..size {
                weights.push((i, j, params.weight));
            }
        }
    }
    let mut volumes = params.excluded_volumes.clone().unwrap_or_default();
    for vol in &volumes {
        let idx = vol
            .index
            .ok_or_else(|| value_error("the excluded volume has no index"))?;
        for i in 0..n_ats {
            weights.push((i, idx as u32, params.weight));
        }
    }
    let coords = embed_bounds_matrix_with_params(
        bm,
        &EmbedBoundsMatrixParams {
            num_zero_fail: 1,
            weights,
            random_seed: params.random_seed,
            ..Default::default()
        },
    )?;
    let mut conf = Conformer::new(n_ats);
    conf.set_id(0)?;
    for i in 0..n_ats {
        let p = *coords
            .get(i as usize)
            .ok_or_else(|| index_error(i as usize, coords.len()))?;
        conf.set_atom_position(i, p)?;
    }
    for vol in &mut volumes {
        let idx = vol.index.expect("checked above");
        vol.pos = Some(
            *coords
                .get(idx)
                .ok_or_else(|| index_error(idx, coords.len()))?,
        );
    }
    mol.add_conformer(&conf)?;
    Ok(volumes)
}

/// A copy of the bounds matrix with a row and column for each excluded
/// volume, which gets that index; with `smooth_it`, triangle-smoothed. An
/// atom past the matrix is `Error::Index`. Mirrors
/// `EmbedLib.AddExcludedVolumes`.
pub fn add_excluded_volumes(
    bm: &[Vec<f64>],
    excluded_volumes: &mut [ExcludedVolume],
) -> Result<Vec<Vec<f64>>> {
    add_excluded_volumes_with_smooth_it(bm, excluded_volumes, true)
}

/// Without smoothing unless `smooth_it`. Mirrors
/// `EmbedLib.AddExcludedVolumes(bm, excludedVolumes, smoothIt)`.
#[allow(clippy::needless_range_loop)] // writes a row and a column together
pub fn add_excluded_volumes_with_smooth_it(
    bm: &[Vec<f64>],
    excluded_volumes: &mut [ExcludedVolume],
    smooth_it: bool,
) -> Result<Vec<Vec<f64>>> {
    let o_dim = bm.len();
    let dim = o_dim + excluded_volumes.len();
    let mut res = vec![vec![0.0; dim]; dim];
    for (i, row) in bm.iter().enumerate() {
        res[i][..row.len()].copy_from_slice(row);
    }
    for (i, vol) in excluded_volumes.iter_mut().enumerate() {
        let bm_idx = o_dim + i;
        vol.index = Some(bm_idx);
        for j in 0..bm_idx {
            res[bm_idx][j] = vol.exclusion_dist;
            res[j][bm_idx] = 1000.0;
        }
        for (indices, min_v, max_v) in &vol.feat_info {
            for &index in indices {
                set(&mut res, bm_idx, index as usize, *min_v)?;
                set(&mut res, index as usize, bm_idx, *max_v)?;
            }
        }
        for j in bm_idx + 1..dim {
            res[bm_idx][j] = 0.0;
            res[j][bm_idx] = 1000.0;
        }
    }
    if smooth_it {
        do_triangle_smoothing(&mut res)?;
    }
    Ok(res)
}

/// Parameters for [`update_pharmacophore_bounds_with_params`], with
/// Python's defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UpdatePharmacophoreBoundsParams {
    /// Adds points for multi-atom features' directions.
    pub use_dirs: bool,
    /// Unused, as in Python.
    pub dir_length: f64,
}

impl Default for UpdatePharmacophoreBoundsParams {
    fn default() -> UpdatePharmacophoreBoundsParams {
        UpdatePharmacophoreBoundsParams {
            use_dirs: false,
            dir_length: default_feat_length(),
        }
    }
}

/// The bounds matrix with the pharmacophore's bounds between its features'
/// atoms, and a center row for each multi-atom feature
/// ([`replace_group`]). Mirrors `EmbedLib.UpdatePharmacophoreBounds`.
pub fn update_pharmacophore_bounds(
    bm: Vec<Vec<f64>>,
    atom_match: &[Vec<u32>],
    pcophore: &Pharmacophore,
) -> Result<Vec<Vec<f64>>> {
    update_pharmacophore_bounds_with_params(bm, atom_match, pcophore, &Default::default())
}

/// Mirrors `EmbedLib.UpdatePharmacophoreBounds(bm, atomMatch, pcophore,
/// useDirs, dirLength)`.
pub fn update_pharmacophore_bounds_with_params(
    mut bm: Vec<Vec<f64>>,
    atom_match: &[Vec<u32>],
    pcophore: &Pharmacophore,
    params: &UpdatePharmacophoreBoundsParams,
) -> Result<Vec<Vec<f64>>> {
    let mut replace_map = vec![None; atom_match.len()];
    for (i, match_i) in atom_match.iter().enumerate() {
        if match_i.len() > 1 {
            // Python passes useDirs, but not dirLength, on.
            let group = ReplaceGroupParams {
                use_dirs: params.use_dirs,
                ..Default::default()
            };
            let (new_bm, idx) = replace_group_with_params(match_i, &bm, &group)?;
            bm = new_bm;
            replace_map[i] = Some(idx);
        }
    }
    let first = |i: usize| -> Result<usize> {
        match replace_map[i] {
            Some(idx) => Ok(idx),
            None => atom_match[i]
                .first()
                .map(|&a| a as usize)
                .ok_or_else(|| index_error(0, 0)),
        }
    };
    for i in 0..atom_match.len() {
        let mi = first(i)?;
        for j in i + 1..atom_match.len() {
            let mj = first(j)?;
            let (idx0, idx1) = if mi < mj { (mi, mj) } else { (mj, mi) };
            set(&mut bm, idx0, idx1, pcophore.get_upper_bound(i, j)?)?;
            set(&mut bm, idx1, idx0, pcophore.get_lower_bound(i, j)?)?;
        }
    }
    Ok(bm)
}

/// The volume of the tetrahedron of a center's neighbors (and the center
/// itself for three), ordered by CIP rank: negative for R, positive for S,
/// and 0 for an atom without a CIP code.
fn chiral_volume(mol: &mut Mol, center_idx: u32, conf_id: i32) -> Result<f64> {
    assign_stereochemistry(mol)?;
    let conf = mol.get_conformer_with_id(conf_id)?;
    let center = mol.get_atom_with_idx(center_idx)?;
    if !center.has_prop("_CIPCode") {
        return Ok(0.0);
    }
    let mut nbr_ranks = Vec::new();
    for nbr in center.get_neighbors() {
        let rank = i64::from(nbr.get_unsigned_prop("_CIPRank")?);
        nbr_ranks.push((rank, conf.get_atom_position(nbr.get_idx())?));
    }
    if nbr_ranks.len() == 3 {
        nbr_ranks.push((-1, conf.get_atom_position(center_idx)?));
    }
    if nbr_ranks.len() != 4 {
        return Err(index_error(nbr_ranks.len(), 4));
    }
    nbr_ranks.sort_by_key(|&(rank, _)| rank);
    let ps: Vec<_> = nbr_ranks.into_iter().map(|(_, p)| p).collect();
    let v1 = ps[0] - ps[3];
    let v2 = ps[1] - ps[3];
    let v3 = ps[2] - ps[3];
    Ok(v1.dot_product(&v2.cross_product(&v3)))
}

/// The chiral volume of an atom (Blaney and Dixon, Rev. Comp. Chem. V,
/// 299-335 (1994), Figure 7), from the CIP ranks of its neighbors, which
/// `assign_stereochemistry` sets with legacy stereo perception: negative
/// for R, positive for S, and 0 for an atom without a CIP code. A
/// neighbor without a rank is `Error::Key`. Mirrors
/// `EmbedLib.ComputeChiralVolume`.
pub fn compute_chiral_volume(mol: &mut Mol, center_idx: u32) -> Result<f64> {
    chiral_volume(mol, center_idx, -1)
}

/// In the conformer with `conf_id`. Mirrors
/// `EmbedLib.ComputeChiralVolume(mol, centerIdx, confId)`.
pub fn compute_chiral_volume_with_conf_id(
    mol: &mut Mol,
    center_idx: u32,
    conf_id: i32,
) -> Result<f64> {
    chiral_volume(mol, center_idx, conf_id)
}

/// Parameters for [`embed_pharmacophore_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct EmbedPharmacophoreParams {
    /// The seed of the first attempt, the others' 10 apart; 0 or less for
    /// 1.
    pub random_seed: i32,
    /// The most attempts.
    pub count: u32,
    /// Triangle-smooths the molecule's bounds first.
    pub smooth_first: bool,
    /// Python logs failures unless `silent`; the port does not log.
    pub silent: bool,
    /// The molecule's bounds matrix; `None` to compute it. Python smooths a
    /// given matrix in place; this smooths a copy.
    pub bounds: Option<Vec<Vec<f64>>>,
    /// Excluded volumes to add to the bounds (an in-out argument, A9).
    pub excluded_volumes: Option<Vec<ExcludedVolume>>,
    /// The most embeddings to make; 0 or less for `count`.
    pub target_number: i32,
    /// Adds points for multi-atom features' directions.
    pub use_dirs: bool,
}

impl Default for EmbedPharmacophoreParams {
    fn default() -> EmbedPharmacophoreParams {
        EmbedPharmacophoreParams {
            random_seed: -1,
            count: 10,
            smooth_first: true,
            silent: false,
            bounds: None,
            excluded_volumes: None,
            target_number: -1,
            use_dirs: false,
        }
    }
}

/// What [`embed_pharmacophore`] returns: the bounds matrix, the
/// embeddings, and the number of failed attempts.
pub type EmbeddedPharmacophore = (Vec<Vec<f64>>, Vec<Mol>, u32);

/// The pharmacophore's embeddings of the molecule: up to `count` attempts,
/// keeping those whose specified stereocenters keep their chirality. Returns
/// the bounds matrix with the pharmacophore's bounds, the embeddings (each
/// a copy of the molecule with one conformer), and the number of attempts
/// that failed. Bounds that cannot be smoothed are `Error::Value`. Mirrors
/// `EmbedLib.EmbedPharmacophore`.
pub fn embed_pharmacophore(
    mol: &Mol,
    atom_match: &[Vec<u32>],
    pcophore: &Pharmacophore,
) -> Result<EmbeddedPharmacophore> {
    embed_pharmacophore_with_params(mol, atom_match, pcophore, &Default::default())
}

/// Mirrors `EmbedLib.EmbedPharmacophore(mol, atomMatch, pcophore,
/// randomSeed, count, smoothFirst, silent, bounds, excludedVolumes,
/// targetNumber, useDirs)`.
pub fn embed_pharmacophore_with_params(
    mol: &Mol,
    atom_match: &[Vec<u32>],
    pcophore: &Pharmacophore,
    params: &EmbedPharmacophoreParams,
) -> Result<EmbeddedPharmacophore> {
    embed_pharmacophore_with_outputs(mol, atom_match, pcophore, params).map(|(res, _)| res)
}

/// With the excluded volumes as the call leaves them: indexed in the
/// bounds, at their positions in the last embedding (A9). Mirrors
/// `EmbedLib.EmbedPharmacophore(..., excludedVolumes)`, which changes them.
pub fn embed_pharmacophore_with_outputs(
    mol: &Mol,
    atom_match: &[Vec<u32>],
    pcophore: &Pharmacophore,
    params: &EmbedPharmacophoreParams,
) -> Result<(EmbeddedPharmacophore, Vec<ExcludedVolume>)> {
    let chiral_centers = find_mol_chiral_centers(&mut mol.clone())?;
    let mut bounds = match &params.bounds {
        Some(bounds) => bounds.clone(),
        None => get_molecule_bounds_matrix(mol)?,
    };
    if params.smooth_first {
        do_triangle_smoothing(&mut bounds)?;
    }
    let update = UpdatePharmacophoreBoundsParams {
        use_dirs: params.use_dirs,
        ..Default::default()
    };
    let mut bm = update_pharmacophore_bounds_with_params(bounds, atom_match, pcophore, &update)?;
    let mut volumes = params.excluded_volumes.clone().unwrap_or_default();
    if !volumes.is_empty() {
        bm = add_excluded_volumes_with_smooth_it(&bm, &mut volumes, false)?;
    }
    if !do_triangle_smoothing(&mut bm)? {
        return Err(value_error("could not smooth bounds matrix"));
    }
    let target_number = if params.target_number <= 0 {
        params.count as usize
    } else {
        params.target_number as usize
    };
    let mut n_failed = 0;
    let mut res = Vec::new();
    for i in 0..params.count {
        let mut m2 = mol.clone();
        let seed = if params.random_seed <= 0 {
            i as i32 * 10 + 1
        } else {
            i as i32 * 10 + params.random_seed
        };
        let embed = EmbedMolParams {
            atom_match: Some(atom_match.to_vec()),
            random_seed: seed,
            excluded_volumes: (!volumes.is_empty()).then(|| volumes.clone()),
            ..Default::default()
        };
        match embed_mol_with_outputs(&mut m2, &bm, &embed) {
            Err(Error::Value { .. }) => {
                n_failed += 1;
                continue;
            }
            Err(err) => return Err(err),
            Ok(placed) => {
                if !volumes.is_empty() {
                    volumes = placed;
                }
            }
        }
        let mut keep_it = true;
        for (idx, stereo) in &chiral_centers {
            if stereo == "R" || stereo == "S" {
                let vol = compute_chiral_volume(&mut m2, *idx)?;
                if (stereo == "R" && vol >= 0.0) || (stereo == "S" && vol <= 0.0) {
                    keep_it = false;
                    break;
                }
            }
        }
        if keep_it {
            res.push(m2);
        }
        if res.len() == target_number {
            break;
        }
    }
    Ok(((bm, res, n_failed), volumes))
}

/// Python's `isNaN`, which is false for every value except NaN on Windows.
/// Mirrors `EmbedLib.isNaN`.
pub fn is_nan(v: f64) -> bool {
    cfg!(windows) && v.is_nan()
}

/// Parameters for [`optimize_mol_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct OptimizeMolParams {
    /// The pharmacophore's atoms, whose pairs across features keep their
    /// bounds.
    pub atom_matches: Option<Vec<Vec<u32>>>,
    /// Excluded volumes, placed by an embedding, that atoms keep away from.
    pub excluded_volumes: Option<Vec<ExcludedVolume>>,
    /// The force constant of the constraints.
    pub force_constant: f64,
    /// The most minimizations after the first.
    pub max_passes: u32,
}

impl Default for OptimizeMolParams {
    fn default() -> OptimizeMolParams {
        OptimizeMolParams {
            atom_matches: None,
            excluded_volumes: None,
            force_constant: 1200.0,
            max_passes: 5,
        }
    }
}

/// Minimizes the molecule's conformer with UFF, subject to the bounds of
/// the pharmacophore's atoms, and returns its energy (with the constraints)
/// before and after; `(-1, -1)` when UFF cannot be set up for it. Mirrors
/// `EmbedLib.OptimizeMol`.
pub fn optimize_mol(mol: &mut Mol, bm: &[Vec<f64>]) -> Result<(f64, f64)> {
    optimize_mol_with_params(mol, bm, &OptimizeMolParams::default())
}

/// Mirrors `EmbedLib.OptimizeMol(mol, bm, atomMatches, excludedVolumes,
/// forceConstant, maxPasses)`. A volume without a position or index is
/// `Error::Value`.
pub fn optimize_mol_with_params(
    mol: &mut Mol,
    bm: &[Vec<f64>],
    params: &OptimizeMolParams,
) -> Result<(f64, f64)> {
    let n_ats = mol.get_num_atoms();
    let positions: Vec<_> = {
        let conf = mol.get_conformer();
        match conf {
            Ok(conf) => (0..n_ats)
                .map(|i| conf.get_atom_position(i))
                .collect::<Result<_>>()?,
            Err(_) => Vec::new(),
        }
    };
    let Ok(mut ff) = uff_get_molecule_force_field(mol) else {
        return Ok((-1.0, -1.0));
    };
    let force_constant = params.force_constant;
    if let Some(matches) = params.atom_matches.as_ref() {
        for k in 0..matches.len() {
            for &i in &matches[k] {
                for other in &matches[k + 1..] {
                    for &j in other {
                        let (i, j) = if j < i { (j, i) } else { (i, j) };
                        let (iu, ju) = (i as usize, j as usize);
                        ff.add_distance_constraint(
                            i,
                            j,
                            at(bm, ju, iu)?,
                            at(bm, iu, ju)?,
                            force_constant,
                        )?;
                    }
                }
            }
        }
    }
    for (idx, ex_vol) in (n_ats..).zip(params.excluded_volumes.iter().flatten()) {
        let pos = ex_vol
            .pos
            .ok_or_else(|| value_error("the excluded volume has no position"))?;
        let vol_idx = ex_vol
            .index
            .ok_or_else(|| value_error("the excluded volume has no index"))?;
        ff.add_extra_point_with_fixed(pos.x, pos.y, pos.z, true);
        let indices: BTreeSet<u32> = ex_vol
            .feat_info
            .iter()
            .flat_map(|(ids, _, _)| ids.iter().copied())
            .collect();
        for i in 0..n_ats {
            if !indices.contains(&i) {
                let v = positions[i as usize] - pos;
                let d = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
                if d < 5.0 {
                    ff.add_distance_constraint(
                        i,
                        idx,
                        ex_vol.exclusion_dist,
                        1000.0,
                        force_constant,
                    )?;
                }
            } else {
                let iu = i as usize;
                ff.add_distance_constraint(
                    i,
                    idx,
                    at(bm, vol_idx, iu)?,
                    at(bm, iu, vol_idx)?,
                    force_constant,
                )?;
            }
        }
    }
    ff.initialize()?;
    let e1 = ff.calc_energy()?;
    if is_nan(e1) {
        return Err(value_error("bogus energy"));
    }
    let mut needs_more = ff.minimize()?;
    let mut n_passes = 0;
    while needs_more != 0 && n_passes < params.max_passes {
        needs_more = ff.minimize()?;
        n_passes += 1;
    }
    let e2 = ff.calc_energy()?;
    if is_nan(e2) {
        return Err(value_error("bogus energy"));
    }
    Ok((e1, e2))
}

/// Python's `rdkit.ML.Data.Stats.MeanAndDev`: the mean and sample standard
/// deviation, summed left to right as Python's `sum` does over NumPy values.
fn mean_and_dev(vect: &[f64]) -> (f64, f64) {
    let n = vect.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mut total = 0.0;
    for &x in vect {
        total += x;
    }
    let mean = total / n as f64;
    if n == 1 {
        return (mean, 0.0);
    }
    let mut squares = 0.0;
    for &x in vect {
        let v = x - mean;
        squares += v * v;
    }
    (mean, (squares / (n - 1) as f64).sqrt())
}

/// The statistics of a molecule's embeddings in [`embed_one`]: the mean
/// and sample standard deviation of the energies E1 (of each embedding,
/// with constraints), E2 (optimized, with constraints), E3 (of that
/// geometry, without them), and E4 (optimized without them), and the number
/// of embeddings that failed.
pub type EmbedOneStats = (f64, f64, f64, f64, f64, f64, f64, f64, u32);

/// Parameters for [`embed_one_with_params`], with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct EmbedOneParams {
    /// The most attempts at embedding.
    pub count: u32,
    /// Prints the statistics, as Python does, unless `silent`.
    pub silent: bool,
    /// The rest of [`EmbedPharmacophoreParams`] (Python's keyword
    /// arguments); its `count` and `silent` are these.
    pub embed: EmbedPharmacophoreParams,
}

impl Default for EmbedOneParams {
    fn default() -> EmbedOneParams {
        EmbedOneParams {
            count: 1,
            silent: false,
            embed: EmbedPharmacophoreParams::default(),
        }
    }
}

/// Embeds the molecule to match the features, optimizes each embedding
/// with and then without the pharmacophore's constraints, and returns the
/// energies' statistics; prints them as `name(failures): E1(dev) -> ...`
/// unless silent. Mirrors `EmbedLib.EmbedOne`.
pub fn embed_one(
    mol: &Mol,
    name: &str,
    match_: &[MolChemicalFeature<'_>],
    pcophore: &Pharmacophore,
) -> Result<EmbedOneStats> {
    embed_one_with_params(mol, name, match_, pcophore, &EmbedOneParams::default())
}

/// Mirrors `EmbedLib.EmbedOne(mol, name, match, pcophore, count, silent,
/// **kwargs)`.
pub fn embed_one_with_params(
    mol: &Mol,
    name: &str,
    match_: &[MolChemicalFeature<'_>],
    pcophore: &Pharmacophore,
    params: &EmbedOneParams,
) -> Result<EmbedOneStats> {
    let atom_match: Vec<Vec<u32>> = match_.iter().map(|x| x.get_atom_ids()).collect();
    let embed = EmbedPharmacophoreParams {
        count: params.count,
        silent: params.silent,
        ..params.embed.clone()
    };
    let (bm, ms, n_failed) = embed_pharmacophore_with_params(mol, &atom_match, pcophore, &embed)?;
    let (mut e1s, mut e2s, mut e3s, mut e4s) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for mut m in ms {
        let constrained = OptimizeMolParams {
            atom_matches: Some(atom_match.clone()),
            ..Default::default()
        };
        let (e1, e2) = match optimize_mol_with_params(&mut m, &bm, &constrained) {
            Err(Error::Value { .. }) => continue,
            other => other?,
        };
        e1s.push(e1);
        e2s.push(e2);
        let (e3, e4) = match optimize_mol(&mut m, &bm) {
            Err(Error::Value { .. }) => continue,
            other => other?,
        };
        e3s.push(e3);
        e4s.push(e4);
    }
    let (e1, e1d) = mean_and_dev(&e1s);
    let (e2, e2d) = mean_and_dev(&e2s);
    let (e3, e3d) = mean_and_dev(&e3s);
    let (e4, e4d) = mean_and_dev(&e4s);
    if !params.silent {
        println!(
            "{name}({n_failed}): {e1:.2}({e1d:.2}) -> {e2:.2}({e2d:.2}) : {e3:.2}({e3d:.2}) -> {e4:.2}({e4d:.2})"
        );
    }
    Ok((e1, e1d, e2, e2d, e3, e3d, e4, e4d, n_failed))
}

/// The molecule's features of each family, as Python's `_getFeatDict` and
/// `GetFeaturesForMol(includeOnly=family)` find them.
fn features_of<'a>(
    mol: &'a Mol,
    feat_factory: &'a MolChemicalFeatureFactory,
    family: &str,
) -> Result<Vec<MolChemicalFeature<'a>>> {
    let count = feat_factory.get_num_mol_features_with_include_only(mol, family)?;
    let params = GetMolFeatureParams {
        include_only: family.to_owned(),
        conf_id: -1,
    };
    (0..count)
        .map(|i| feat_factory.get_mol_feature_with_params(mol, i, &params))
        .collect()
}

/// Each feature's matches in the molecule, by the feature's family, or
/// `None` if a feature has none. Mirrors `EmbedLib.MatchFeatsToMol`, whose
/// first value is whether the second is not `None`.
pub fn match_feats_to_mol<'a>(
    mol: &'a Mol,
    feat_factory: &'a MolChemicalFeatureFactory,
    features: &[FreeChemicalFeature],
) -> Result<(bool, Option<Vec<Vec<MolChemicalFeature<'a>>>>)> {
    let mut res = Vec::with_capacity(features.len());
    for feat in features {
        let matches = features_of(mol, feat_factory, &feat.get_family())?;
        if matches.is_empty() {
            return Ok((false, None));
        }
        res.push(matches);
    }
    Ok((true, Some(res)))
}

/// Each of the pharmacophore's features' matches in the molecule. Mirrors
/// `EmbedLib.MatchPharmacophoreToMol`.
pub fn match_pharmacophore_to_mol<'a>(
    mol: &'a Mol,
    feat_factory: &'a MolChemicalFeatureFactory,
    pcophore: &Pharmacophore,
) -> Result<(bool, Option<Vec<Vec<MolChemicalFeature<'a>>>>)> {
    match_feats_to_mol(mol, feat_factory, pcophore.get_features())
}

/// Every combination of one element of each sequence, in order; one empty
/// combination for no sequences. Python's is a generator. Mirrors
/// `EmbedLib.CombiEnum`.
pub fn combi_enum<T: Clone>(sequence: &[Vec<T>]) -> Vec<Vec<T>> {
    let mut res = vec![Vec::new()];
    for choices in sequence {
        let mut next = Vec::with_capacity(res.len() * choices.len());
        for prefix in &res {
            for entry in choices {
                let mut combo = prefix.clone();
                combo.push(entry.clone());
                next.push(combo);
            }
        }
        res = next;
    }
    res
}

/// The bounds matrix without the rows (and columns) of points whose upper
/// bounds to every one of `indices` (looking only at later points) are at
/// least 4. Mirrors `EmbedLib.DownsampleBoundsMatrix`.
pub fn downsample_bounds_matrix(bm: &[Vec<f64>], indices: &[usize]) -> Result<Vec<Vec<f64>>> {
    downsample_bounds_matrix_with_max_thresh(bm, indices, 4.0)
}

/// With threshold `max_thresh`. Mirrors
/// `EmbedLib.DownsampleBoundsMatrix(bm, indices, maxThresh)`.
pub fn downsample_bounds_matrix_with_max_thresh(
    bm: &[Vec<f64>],
    indices: &[usize],
    max_thresh: f64,
) -> Result<Vec<Vec<f64>>> {
    let n_pts = bm.len();
    if indices.is_empty() {
        return Ok(Vec::new());
    }
    let mut mask = vec![false; n_pts];
    for &idx in indices.iter().collect::<BTreeSet<_>>() {
        *mask.get_mut(idx).ok_or_else(|| index_error(idx, n_pts))? = true;
        for (j, keep) in mask.iter_mut().enumerate().skip(idx + 1) {
            if at(bm, idx, j)? < max_thresh {
                *keep = true;
            }
        }
    }
    let keep: Vec<usize> = (0..n_pts).filter(|&i| mask[i]).collect();
    Ok(keep
        .iter()
        .map(|&i| keep.iter().map(|&j| bm[i][j]).collect())
        .collect())
}

/// Whether every pair of single-atom features can be as far apart as the
/// pharmacophore needs: the molecule's upper bound above the
/// pharmacophore's lower bound, and its lower bound below the
/// pharmacophore's upper bound. Mirrors `EmbedLib.CoarseScreenPharmacophore`.
pub fn coarse_screen_pharmacophore(
    atom_match: &[Vec<u32>],
    bounds: &[Vec<f64>],
    pcophore: &Pharmacophore,
) -> Result<bool> {
    for k in 0..atom_match.len() {
        if atom_match[k].len() != 1 {
            continue;
        }
        for l in k + 1..atom_match.len() {
            if atom_match[l].len() != 1 {
                continue;
            }
            let (a, b) = (atom_match[k][0] as usize, atom_match[l][0] as usize);
            let (idx0, idx1) = if b < a { (b, a) } else { (a, b) };
            if at(bounds, idx1, idx0)? >= pcophore.get_upper_bound(k, l)?
                || at(bounds, idx0, idx1)? <= pcophore.get_lower_bound(k, l)?
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Whether the features' atoms' shortest topological distances are within
/// the pharmacophore's 2D bounds. Mirrors `EmbedLib.Check2DBounds`.
pub fn check_2d_bounds(
    atom_match: &[Vec<u32>],
    mol: &Mol,
    pcophore: &Pharmacophore,
) -> Result<bool> {
    let dm = get_distance_matrix_with_params(
        mol,
        &GetDistanceMatrixParams {
            use_bo: false,
            use_atom_wts: false,
        },
    )?;
    let n_feats = atom_match.len();
    for i in 0..n_feats {
        for j in i + 1..n_feats {
            let lower_b = pcophore.get_lower_bound_2d(i, j)?;
            let upper_b = pcophore.get_upper_bound_2d(i, j)?;
            let mut dij = 10000.0f64;
            for &atom_i in &atom_match[i] {
                for &atom_j in &atom_match[j] {
                    dij = dij.min(at(&dm, atom_i as usize, atom_j as usize)?);
                }
            }
            if dij < lower_b as f64 || dij > upper_b as f64 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Python's `ChemicalFeatures.GetAtomMatch` of borrowed features: each
/// feature's atoms, or nothing if two share an atom.
fn atom_match_of(feats: &[&MolChemicalFeature<'_>]) -> Result<Vec<Vec<u32>>> {
    let mut seen = BTreeSet::new();
    let mut result = Vec::with_capacity(feats.len());
    for feat in feats {
        let ids = feat.get_atom_ids();
        for &idx in &ids {
            // RDKit's flags cover atoms below its maxAts of 1024.
            if idx >= 1024 {
                return Err(Error::Index {
                    index: idx,
                    len: 1024,
                });
            }
            if !seen.insert(idx) {
                return Ok(Vec::new());
            }
        }
        result.push(ids);
    }
    Ok(result)
}

/// The atom match of a partial assignment of features, if it passes the
/// screens, as Python's `_checkMatch`.
fn check_match(
    m: &[&MolChemicalFeature<'_>],
    mol: Option<&Mol>,
    bounds: &[Vec<f64>],
    pcophore: &Pharmacophore,
    use_2d_limits: bool,
) -> Result<Option<Vec<Vec<u32>>>> {
    let atom_match = atom_match_of(m)?;
    if atom_match.is_empty() {
        return Ok(None);
    }
    if use_2d_limits {
        let mol = mol.ok_or_else(|| value_error("2D limits need the molecule"))?;
        if !check_2d_bounds(&atom_match, mol, pcophore)? {
            return Ok(None);
        }
    }
    if !coarse_screen_pharmacophore(&atom_match, bounds, pcophore)? {
        return Ok(None);
    }
    Ok(Some(atom_match))
}

/// A feature assignment and its atoms.
pub type FeatureAssignment<'f, 'a> = (Vec<&'f MolChemicalFeature<'a>>, Vec<Vec<u32>>);

/// Parameters for [`constrained_enum_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, Default)]
pub struct ConstrainedEnumParams<'m> {
    /// Also screens assignments by the pharmacophore's topological bounds,
    /// which needs `mol`.
    pub use_2d_limits: bool,
    /// The molecule, for `use_2d_limits`.
    pub mol: Option<&'m Mol>,
}

/// Every assignment of one match to each feature, in order, whose atoms
/// are distinct and pass the bounds screens as it grows. Python's is a
/// generator whose recursion arguments are not bound. Mirrors
/// `EmbedLib.ConstrainedEnum`.
pub fn constrained_enum<'f, 'a>(
    matches: &'f [Vec<MolChemicalFeature<'a>>],
    pcophore: &Pharmacophore,
    bounds: &[Vec<f64>],
) -> Result<Vec<FeatureAssignment<'f, 'a>>> {
    constrained_enum_with_params(matches, pcophore, bounds, &Default::default())
}

/// Mirrors `EmbedLib.ConstrainedEnum(matches, mol, pcophore, bounds,
/// use2DLimits)`.
pub fn constrained_enum_with_params<'f, 'a>(
    matches: &'f [Vec<MolChemicalFeature<'a>>],
    pcophore: &Pharmacophore,
    bounds: &[Vec<f64>],
    params: &ConstrainedEnumParams<'_>,
) -> Result<Vec<FeatureAssignment<'f, 'a>>> {
    let mut res = Vec::new();
    if matches.is_empty() {
        res.push((Vec::new(), Vec::new()));
        return Ok(res);
    }
    let mut so_far = Vec::new();
    enumerate_assignments(matches, pcophore, bounds, params, 0, &mut so_far, &mut res)?;
    Ok(res)
}

fn enumerate_assignments<'f, 'a>(
    matches: &'f [Vec<MolChemicalFeature<'a>>],
    pcophore: &Pharmacophore,
    bounds: &[Vec<f64>],
    params: &ConstrainedEnumParams<'_>,
    index: usize,
    so_far: &mut Vec<&'f MolChemicalFeature<'a>>,
    res: &mut Vec<FeatureAssignment<'f, 'a>>,
) -> Result<()> {
    let last = index == matches.len() - 1;
    for entry in &matches[index] {
        so_far.push(entry);
        // Python screens each partial assignment from the second feature
        // on; a lone first feature only needs distinct atoms.
        let atom_match = if index != 0 {
            check_match(so_far, params.mol, bounds, pcophore, params.use_2d_limits)?
        } else if last {
            Some(atom_match_of(so_far)?).filter(|m| !m.is_empty())
        } else {
            None
        };
        if last {
            if let Some(atom_match) = atom_match {
                res.push((so_far.clone(), atom_match));
            }
        } else if index == 0 || atom_match.is_some() {
            enumerate_assignments(matches, pcophore, bounds, params, index + 1, so_far, res)?;
        }
        so_far.pop();
    }
    Ok(())
}

/// A match of features to a pharmacophore: the bounds, the features
/// assigned to the pharmacophore's, and the bounds' sizes before and after
/// downsampling.
pub type PharmacophoreMatch<'f, 'a> = (
    Vec<Vec<f64>>,
    Vec<&'f MolChemicalFeature<'a>>,
    (usize, usize),
);

/// Parameters for [`match_pharmacophore_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Debug, Default)]
pub struct MatchPharmacophoreParams<'m> {
    /// Drops the bounds' rows far from the pharmacophore before smoothing.
    pub use_downsampling: bool,
    /// Also screens by the pharmacophore's topological bounds, which needs
    /// `mol`.
    pub use_2d_limits: bool,
    /// The molecule, for `use_2d_limits`.
    pub mol: Option<&'m Mol>,
    /// Excluded volumes, with feature information for each of the
    /// pharmacophore's features.
    pub excluded_volumes: Option<Vec<ExcludedVolume>>,
    /// Adds points for multi-atom features' directions.
    pub use_dirs: bool,
}

/// The first assignment of matches to the pharmacophore's features whose
/// bounds (with the pharmacophore's) can be smoothed: the bounds, the
/// assignment, and the bounds' sizes before and after downsampling; `None`
/// if none can, where Python returns `(1, None, None, None)` (and
/// `(0, ...)` otherwise). Mirrors `EmbedLib.MatchPharmacophore`.
pub fn match_pharmacophore<'f, 'a>(
    matches: &'f [Vec<MolChemicalFeature<'a>>],
    bounds: &[Vec<f64>],
    pcophore: &Pharmacophore,
) -> Result<Option<PharmacophoreMatch<'f, 'a>>> {
    match_pharmacophore_with_params(matches, bounds, pcophore, &Default::default())
}

/// Mirrors `EmbedLib.MatchPharmacophore(matches, bounds, pcophore,
/// useDownsampling, use2DLimits, mol, excludedVolumes, useDirs)`.
pub fn match_pharmacophore_with_params<'f, 'a>(
    matches: &'f [Vec<MolChemicalFeature<'a>>],
    bounds: &[Vec<f64>],
    pcophore: &Pharmacophore,
    params: &MatchPharmacophoreParams<'_>,
) -> Result<Option<PharmacophoreMatch<'f, 'a>>> {
    let enum_params = ConstrainedEnumParams {
        use_2d_limits: params.use_2d_limits,
        mol: params.mol,
    };
    for (m, atom_match) in constrained_enum_with_params(matches, pcophore, bounds, &enum_params)? {
        let update = UpdatePharmacophoreBoundsParams {
            use_dirs: params.use_dirs,
            ..Default::default()
        };
        let mut bm = update_pharmacophore_bounds_with_params(
            bounds.to_vec(),
            &atom_match,
            pcophore,
            &update,
        )?;
        let mut local_evs = Vec::new();
        if let Some(volumes) = params.excluded_volumes.as_ref().filter(|v| !v.is_empty()) {
            for ev in volumes {
                let mut feat_info = Vec::new();
                for (i, entry) in atom_match.iter().enumerate() {
                    let (_, min_v, max_v) = ev
                        .feat_info
                        .get(i)
                        .ok_or_else(|| index_error(i, ev.feat_info.len()))?;
                    feat_info.push((entry.clone(), *min_v, *max_v));
                }
                local_evs.push(ExcludedVolume::new_with_params(
                    feat_info,
                    &ExcludedVolumeParams {
                        index: ev.index,
                        exclusion_dist: ev.exclusion_dist,
                    },
                )?);
            }
            bm = add_excluded_volumes_with_smooth_it(&bm, &mut local_evs, false)?;
        }
        let sz = bm.len();
        if params.use_downsampling {
            let mut indices: Vec<usize> =
                atom_match.iter().flatten().map(|&a| a as usize).collect();
            indices.extend(local_evs.iter().filter_map(|vol| vol.index));
            bm = downsample_bounds_matrix(&bm, &indices)?;
        }
        if do_triangle_smoothing(&mut bm)? {
            let n = bm.len();
            return Ok(Some((bm, m, (sz, n))));
        }
    }
    Ok(None)
}

/// Parameters for [`get_all_pharmacophore_matches_with_params`], with
/// Python's defaults (A2).
#[derive(Default)]
pub struct GetAllPharmacophoreMatchesParams<'m> {
    /// Drops the bounds' rows far from the pharmacophore before smoothing.
    pub use_downsampling: bool,
    /// Called with the number of assignments tried after each.
    pub progress_callback: Option<&'m dyn Fn(usize)>,
    /// Also screens by the pharmacophore's topological bounds, when `mol`
    /// is given.
    pub use_2d_limits: bool,
    /// The molecule, for `use_2d_limits`.
    pub mol: Option<&'m Mol>,
}

/// Every assignment of matches to the pharmacophore's features with
/// distinct atoms that passes the screens and whose bounds can be
/// smoothed. Mirrors `EmbedLib.GetAllPharmacophoreMatches`.
pub fn get_all_pharmacophore_matches<'f, 'a>(
    matches: &'f [Vec<MolChemicalFeature<'a>>],
    bounds: &[Vec<f64>],
    pcophore: &Pharmacophore,
) -> Result<Vec<Vec<&'f MolChemicalFeature<'a>>>> {
    get_all_pharmacophore_matches_with_params(matches, bounds, pcophore, &Default::default())
}

/// Mirrors `EmbedLib.GetAllPharmacophoreMatches(matches, bounds, pcophore,
/// useDownsampling, progressCallback, use2DLimits, mol)`.
pub fn get_all_pharmacophore_matches_with_params<'f, 'a>(
    matches: &'f [Vec<MolChemicalFeature<'a>>],
    bounds: &[Vec<f64>],
    pcophore: &Pharmacophore,
    params: &GetAllPharmacophoreMatchesParams<'_>,
) -> Result<Vec<Vec<&'f MolChemicalFeature<'a>>>> {
    let refs: Vec<Vec<&'f MolChemicalFeature<'a>>> =
        matches.iter().map(|m| m.iter().collect()).collect();
    let mut res = Vec::new();
    let mut n_done = 0;
    for m in combi_enum(&refs) {
        let atom_match = atom_match_of(&m)?;
        let pass_2d = match (atom_match.is_empty(), params.use_2d_limits, params.mol) {
            (false, true, Some(mol)) => check_2d_bounds(&atom_match, mol, pcophore)?,
            _ => true,
        };
        if !atom_match.is_empty()
            && pass_2d
            && coarse_screen_pharmacophore(&atom_match, bounds, pcophore)?
        {
            let mut bm = update_pharmacophore_bounds(bounds.to_vec(), &atom_match, pcophore)?;
            if params.use_downsampling {
                let indices: Vec<usize> =
                    atom_match.iter().flatten().map(|&a| a as usize).collect();
                bm = downsample_bounds_matrix(&bm, &indices)?;
            }
            if do_triangle_smoothing(&mut bm)? {
                res.push(m);
            }
            n_done += 1;
            if let Some(callback) = params.progress_callback {
                callback(n_done);
            }
        }
    }
    Ok(res)
}
