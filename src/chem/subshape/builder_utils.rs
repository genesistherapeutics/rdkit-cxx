// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Subshape.BuilderUtils`, a Python module: finding
//! the terminal and skeleton points of a shape.
//!
//! Python solves each point's moments with `numpy.linalg.eigh`; this uses
//! Jacobi rotations, which agree with it to within rounding, with
//! directions possibly of the opposite sign (which the aligner ignores).
//! Python's `AssignMolFeatsToPoints` also prints each point's features,
//! which this does not.

use crate::chem::Mol;
use crate::chem::rd_mol_chemical_features::MolChemicalFeatureFactory;
use crate::chem::rdchem::ConformerRef;
use crate::chem::subshape::subshape_objects::{ShapeWithSkeleton, SkeletonPoint};
use crate::geometry::rd_geometry::{
    Point3D, UniformGrid3D, compute_grid_centroid, find_grid_terminal_points,
};
use crate::{Error, Result};

fn value_error(message: &str) -> Error {
    Error::Value {
        message: message.to_owned(),
    }
}

/// The grid's values, as Python's `GetOccupancyVect` lists them.
pub(crate) fn occupancy(grid: &UniformGrid3D) -> Result<Vec<u32>> {
    let vect = grid.get_occupancy_vect();
    (0..vect.len()).map(|i| vect.get_val(i)).collect()
}

fn shape_grid(shape: &ShapeWithSkeleton) -> Result<&UniformGrid3D> {
    shape
        .grid
        .as_ref()
        .ok_or_else(|| value_error("the shape has no grid"))
}

/// The offsets of the grid indices within `win_rad` of a point, in grid
/// steps, as flat indices. Python caches them on the grid, so its later
/// calls with another radius get the first radius's offsets; this computes
/// them for the radius given. Mirrors `BuilderUtils.ComputeGridIndices`.
pub fn compute_grid_indices(shape_grid: &UniformGrid3D, win_rad: f64) -> Vec<i64> {
    let dx = i64::from(shape_grid.get_num_x());
    let dy = i64::from(shape_grid.get_num_y());
    let rad_in_grid = (win_rad / shape_grid.get_spacing()) as i64;
    let mut indices = Vec::new();
    for i in -rad_in_grid..=rad_in_grid {
        for j in -rad_in_grid..=rad_in_grid {
            for k in -rad_in_grid..=rad_in_grid {
                let d = ((i * i + j * j + k * k) as f64).sqrt() as i64;
                if d <= rad_in_grid {
                    indices.push((i * dy + j) * dx + k);
                }
            }
        }
    }
    indices
}

/// The weight and the weighted centroid of the grid's values near `pt`,
/// over [`compute_grid_indices`]' offsets; no weight is `Error::Value`.
/// Mirrors `BuilderUtils.ComputeShapeGridCentroid`.
pub fn compute_shape_grid_centroid(
    pt: &Point3D,
    shape_grid: &UniformGrid3D,
    win_rad: f64,
) -> Result<(u64, Point3D)> {
    let mut count = 0u64;
    let mut centroid = Point3D::new(0.0, 0.0, 0.0);
    let idx_i = i64::from(shape_grid.get_grid_point_index(pt));
    let values = occupancy(shape_grid)?;
    for idx_j in compute_grid_indices(shape_grid, win_rad) {
        let idx = idx_i + idx_j;
        if idx >= 0 && (idx as usize) < values.len() {
            let wt = values[idx as usize];
            let t_pt = shape_grid.get_grid_point_loc(idx as u32)?;
            centroid += t_pt * f64::from(wt);
            count += u64::from(wt);
        }
    }
    if count == 0 {
        return Err(value_error("found no weight in sphere"));
    }
    centroid /= count as f64;
    Ok((count, centroid))
}

/// The shape's terminal points, from its grid. Mirrors
/// `BuilderUtils.FindTerminalPtsFromShape`, whose `maxGridVal` is unused.
pub fn find_terminal_pts_from_shape(
    shape: &ShapeWithSkeleton,
    win_rad: f64,
    fraction: f64,
) -> Result<Vec<SkeletonPoint>> {
    let pts = find_grid_terminal_points(shape_grid(shape)?, win_rad, fraction)?;
    Ok(pts.into_iter().map(SkeletonPoint::from_location).collect())
}

/// The terminal points of a conformer's heavy atoms: for each with fewer
/// than `nbr_count` heavy atoms within `win_rad` (counting itself), their
/// weighted mean position, raising the count until there is one. A
/// conformer without heavy atoms has none (Python loops forever). Mirrors
/// `BuilderUtils.FindTerminalPtsFromConformer`.
pub fn find_terminal_pts_from_conformer(
    conf: &ConformerRef<'_>,
    win_rad: f64,
    nbr_count: usize,
) -> Result<Vec<SkeletonPoint>> {
    let mol = conf.get_owning_mol();
    let n_ats = conf.get_num_atoms();
    let mut nbr_lists: Vec<Vec<(usize, Point3D)>> = vec![Vec::new(); n_ats as usize];
    for i in 0..n_ats {
        if mol.get_atom_with_idx(i)?.get_atomic_num() <= 1 {
            continue;
        }
        let pi = conf.get_atom_position(i)?;
        nbr_lists[i as usize].push((i as usize, pi));
        for j in i + 1..n_ats {
            if mol.get_atom_with_idx(j)?.get_atomic_num() <= 1 {
                continue;
            }
            let pj = conf.get_atom_position(j)?;
            if pi.distance(&pj) < win_rad {
                nbr_lists[i as usize].push((j as usize, pj));
                nbr_lists[j as usize].push((i as usize, pi));
            }
        }
    }
    if nbr_lists.iter().all(Vec::is_empty) {
        return Ok(Vec::new());
    }
    let mut nbr_count = nbr_count;
    loop {
        let mut term_pts = Vec::new();
        for nbr_list in &nbr_lists {
            if nbr_list.is_empty() || nbr_list.len() >= nbr_count {
                continue;
            }
            let mut pos = Point3D::new(0.0, 0.0, 0.0);
            let mut tot_wt = 0.0;
            for &(nbr_j, pos_j) in nbr_list {
                let weight = 1.0 * nbr_list.len() as f64 / nbr_lists[nbr_j].len() as f64;
                pos += pos_j * weight;
                tot_wt += weight;
            }
            pos /= tot_wt;
            term_pts.push(SkeletonPoint::from_location(pos));
        }
        if !term_pts.is_empty() {
            return Ok(term_pts);
        }
        nbr_count += 1;
    }
}

/// The grid's centroid near the midpoint of two points, refined until it
/// moves less than a grid spacing. Mirrors
/// `BuilderUtils.FindGridPointBetweenPoints`.
pub fn find_grid_point_between_points(
    pt1: &Point3D,
    pt2: &Point3D,
    shape_grid: &UniformGrid3D,
    win_rad: f64,
) -> Result<Point3D> {
    let mut center = *pt1 + *pt2;
    center /= 2.0;
    let mut d = 1e8;
    while d > shape_grid.get_spacing() {
        let (_, centroid) = compute_grid_centroid(shape_grid, &center, win_rad)?;
        d = center.distance(&centroid);
        center = centroid;
    }
    Ok(center)
}

/// The points merged into clusters of points within `scale * win_rad` of
/// each cluster's first, each at its mean. Mirrors
/// `BuilderUtils.ClusterTerminalPts`.
pub fn cluster_terminal_pts(pts: &[SkeletonPoint], win_rad: f64, scale: f64) -> Vec<SkeletonPoint> {
    let mut res = Vec::new();
    let mut tagged: Vec<&SkeletonPoint> = pts.iter().collect();
    while !tagged.is_empty() {
        let head = tagged.remove(0);
        let mut curr_set = vec![head];
        let mut i = 0;
        while i < tagged.len() {
            if head.location.distance(&tagged[i].location) < scale * win_rad {
                curr_set.push(tagged.remove(i));
            } else {
                i += 1;
            }
        }
        let mut pt = Point3D::new(0.0, 0.0, 0.0);
        for o in &curr_set {
            pt += o.location;
        }
        pt /= curr_set.len() as f64;
        res.push(SkeletonPoint::from_location(pt));
    }
    res
}

/// Adds terminal points until there are 5: each the grid's centroid near
/// the full grid point farthest from the points so far. Mirrors
/// `BuilderUtils.GetMoreTerminalPoints`.
pub fn get_more_terminal_points(
    shape: &ShapeWithSkeleton,
    pts: &mut Vec<SkeletonPoint>,
    win_rad: f64,
    max_grid_val: f64,
) -> Result<()> {
    get_more_terminal_points_with_target_number(shape, pts, win_rad, max_grid_val, 5)
}

/// Until there are `target_number`. Mirrors
/// `BuilderUtils.GetMoreTerminalPoints(shape, pts, winRad, maxGridVal,
/// targetNumber)`. A grid with no point of `max_grid_val` is `Error::Value`.
pub fn get_more_terminal_points_with_target_number(
    shape: &ShapeWithSkeleton,
    pts: &mut Vec<SkeletonPoint>,
    win_rad: f64,
    max_grid_val: f64,
    target_number: usize,
) -> Result<()> {
    let grid = shape_grid(shape)?;
    let values = occupancy(grid)?;
    while pts.len() < target_number {
        let mut max_min = -1.0;
        let mut best_pt = None;
        for (i, &v) in values.iter().enumerate() {
            if f64::from(v) < max_grid_val {
                continue;
            }
            let pos_i = grid.get_grid_point_loc(i as u32)?;
            let mut min_val = 1e8f64;
            for curr_pt in pts.iter() {
                let dst = pos_i.distance(&curr_pt.location);
                if dst < min_val {
                    min_val = dst;
                }
            }
            if min_val > max_min {
                max_min = min_val;
                best_pt = Some(pos_i);
            }
        }
        let best_pt = best_pt.ok_or_else(|| value_error("the grid has no full point"))?;
        let (_, centroid) = compute_grid_centroid(grid, &best_pt, win_rad)?;
        pts.push(SkeletonPoint::from_location(centroid));
    }
    Ok(())
}

/// The grid's centroid near the full grid point farthest from `loc`; a
/// grid with no point of `max_grid_val` is `Error::Value`. Mirrors
/// `BuilderUtils.FindFarthestGridPoint`.
pub fn find_farthest_grid_point(
    shape: &ShapeWithSkeleton,
    loc: &Point3D,
    win_rad: f64,
    max_grid_val: f64,
) -> Result<Point3D> {
    let grid = shape_grid(shape)?;
    let mut d_max = -1.0;
    let mut res = None;
    for (i, &v) in occupancy(grid)?.iter().enumerate() {
        if f64::from(v) < max_grid_val {
            continue;
        }
        let pos_i = grid.get_grid_point_loc(i as u32)?;
        let dst = pos_i.distance(loc);
        if dst > d_max {
            d_max = dst;
            res = Some(pos_i);
        }
    }
    let res = res.ok_or_else(|| value_error("the grid has no full point"))?;
    Ok(compute_grid_centroid(grid, &res, win_rad)?.1)
}

/// Parameters for [`expand_terminal_pts_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExpandTerminalPtsParams {
    /// The value of full grid points.
    pub max_grid_val: f64,
    /// The number of terminal points to reach.
    pub target_num_pts: usize,
}

impl Default for ExpandTerminalPtsParams {
    fn default() -> ExpandTerminalPtsParams {
        ExpandTerminalPtsParams {
            max_grid_val: 3.0,
            target_num_pts: 5,
        }
    }
}

/// Adds terminal points until there are 5: for one point, the farthest
/// full one; for two, one between them; then by
/// [`get_more_terminal_points`]. Mirrors `BuilderUtils.ExpandTerminalPts`.
pub fn expand_terminal_pts(
    shape: &ShapeWithSkeleton,
    pts: &mut Vec<SkeletonPoint>,
    win_rad: f64,
) -> Result<()> {
    expand_terminal_pts_with_params(shape, pts, win_rad, &ExpandTerminalPtsParams::default())
}

/// Mirrors `BuilderUtils.ExpandTerminalPts(shape, pts, winRad, maxGridVal,
/// targetNumPts)`.
pub fn expand_terminal_pts_with_params(
    shape: &ShapeWithSkeleton,
    pts: &mut Vec<SkeletonPoint>,
    win_rad: f64,
    params: &ExpandTerminalPtsParams,
) -> Result<()> {
    if pts.len() == 1 {
        let pt2 = find_farthest_grid_point(shape, &pts[0].location, win_rad, params.max_grid_val)?;
        pts.push(SkeletonPoint::from_location(pt2));
    }
    if pts.len() == 2 {
        let center = find_grid_point_between_points(
            &pts[0].location,
            &pts[1].location,
            shape_grid(shape)?,
            win_rad,
        )?;
        pts.push(SkeletonPoint::from_location(center));
    }
    if pts.len() < params.target_num_pts {
        get_more_terminal_points_with_target_number(
            shape,
            pts,
            win_rad,
            params.max_grid_val,
            params.target_num_pts,
        )?;
    }
    Ok(())
}

/// Parameters for [`append_skeleton_points_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AppendSkeletonPointsParams {
    /// The value of full grid points.
    pub max_grid_val: f64,
    /// The farthest a point may move to its centroid.
    pub max_dist_c: f64,
    /// Unused, as in Python.
    pub dist_tol: f64,
    /// Points closer than this times the step distance are merged.
    pub sym_factor: f64,
}

impl Default for AppendSkeletonPointsParams {
    fn default() -> AppendSkeletonPointsParams {
        AppendSkeletonPointsParams {
            max_grid_val: 3.0,
            max_dist_c: 15.0,
            dist_tol: 1.5,
            sym_factor: 1.5,
        }
    }
}

/// The terminal points followed by skeleton points: the centroids of the
/// full grid points at least `step_dist` from them, thinned to one point
/// (the fullest) per neighborhood. Mirrors
/// `BuilderUtils.AppendSkeletonPoints`.
pub fn append_skeleton_points(
    shape_grid: &UniformGrid3D,
    term_pts: &[SkeletonPoint],
    win_rad: f64,
    step_dist: f64,
) -> Result<Vec<SkeletonPoint>> {
    append_skeleton_points_with_params(
        shape_grid,
        term_pts,
        win_rad,
        step_dist,
        &Default::default(),
    )
}

/// Mirrors `BuilderUtils.AppendSkeletonPoints(shapeGrid, termPts, winRad,
/// stepDist, maxGridVal, maxDistC, distTol, symFactor)`; Python's
/// `verbose` output is not ported.
pub fn append_skeleton_points_with_params(
    shape_grid: &UniformGrid3D,
    term_pts: &[SkeletonPoint],
    win_rad: f64,
    step_dist: f64,
    params: &AppendSkeletonPointsParams,
) -> Result<Vec<SkeletonPoint>> {
    let n_term_pts = term_pts.len();
    let mut skel_pts = Vec::new();
    for (i, &v) in occupancy(shape_grid)?.iter().enumerate() {
        if f64::from(v) < params.max_grid_val {
            continue;
        }
        let pos_i = shape_grid.get_grid_point_loc(i as u32)?;
        if term_pts
            .iter()
            .all(|pt| pos_i.distance(&pt.location) >= step_dist)
        {
            skel_pts.push(SkeletonPoint::from_location(pos_i));
        }
    }
    let grid_box_volume = shape_grid.get_spacing().powf(3.0);
    let max_vol = 4.0 * std::f64::consts::PI / 3.0 * win_rad.powf(3.0) * params.max_grid_val
        / grid_box_volume;
    let mut i = 0;
    while i < skel_pts.len() {
        let (count, centroid) = compute_grid_centroid(shape_grid, &skel_pts[i].location, win_rad)?;
        if centroid.distance(&skel_pts[i].location) > params.max_dist_c {
            skel_pts.remove(i);
        } else {
            skel_pts[i].frac_vol = count / max_vol;
            skel_pts[i].location = centroid;
            i += 1;
        }
    }
    let mut res: Vec<SkeletonPoint> = term_pts.to_vec();
    res.extend(skel_pts);
    let limit = params.sym_factor * step_dist;
    let mut i = 0;
    while i < res.len() {
        let mut p = None;
        let mut m_frac = 0.0;
        let pt_i = res[i].location;
        let start_j = (i + 1).max(n_term_pts);
        for (j, pt_j) in res.iter().enumerate().skip(start_j) {
            if pt_i.distance(&pt_j.location) < limit && pt_j.frac_vol > m_frac {
                p = Some(j);
                m_frac = pt_j.frac_vol;
            }
        }
        if let Some(p) = p {
            let pt_p = res.remove(p);
            let mut j = start_j;
            while j < res.len() {
                if pt_i.distance(&res[j].location) < limit {
                    res.remove(j);
                } else {
                    j += 1;
                }
            }
            res.push(pt_p);
        }
        i += 1;
    }
    Ok(res)
}

/// The eigenvalues and eigenvectors (as columns) of a symmetric 3x3 matrix,
/// by cyclic Jacobi rotations.
#[allow(clippy::needless_range_loop)] // rotates a row and a column together
fn symmetric_eigen(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..100 {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        if off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for k in 0..3 {
                let (akp, akq) = (a[k][p], a[k][q]);
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let (apk, aqk) = (a[p][k], a[q][k]);
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for row in &mut v {
                let (vkp, vkq) = (row[p], row[q]);
                row[p] = c * vkp - s * vkq;
                row[q] = s * vkp + c * vkq;
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

/// Sets the point's shape moments and directions: the eigenvalues and
/// eigenvectors, largest first, of the covariance of the grid's values
/// within `win_rad`. Mirrors `BuilderUtils.CalculateDirectionsAtPoint`.
pub fn calculate_directions_at_point(
    pt: &mut SkeletonPoint,
    shape_grid: &UniformGrid3D,
    win_rad: f64,
) -> Result<()> {
    let values = occupancy(shape_grid)?;
    let n_grid_pts = values.len() as i64;
    let tmp = win_rad / shape_grid.get_spacing();
    let rad_in_grid = tmp as i64;
    let rad_in_grid2 = (tmp * tmp) as i64;
    let mut cov = [[0.0f64; 3]; 3];
    let dx = i64::from(shape_grid.get_num_x());
    let dy = i64::from(shape_grid.get_num_y());
    let idx = i64::from(shape_grid.get_grid_point_index(&pt.location));
    // Python's // and % floor, as these do for the negative index of a point
    // outside the grid.
    let idx_z = idx.div_euclid(dx * dy);
    let rem = idx.rem_euclid(dx * dy);
    let (idx_y, idx_x) = (rem / dx, rem % dx);
    let mut tot_wt = 0.0;
    for i in -rad_in_grid..rad_in_grid {
        let xi = idx_x + i;
        for j in -rad_in_grid..rad_in_grid {
            let xj = idx_y + j;
            for k in -rad_in_grid..rad_in_grid {
                let xk = idx_z + k;
                let d2 = i * i + j * j + k * k;
                if d2 > rad_in_grid2 && ((d2 as f64).sqrt() as i64) > rad_in_grid {
                    continue;
                }
                let grid_idx = (xk * dy + xj) * dx + xi;
                if grid_idx >= 0 && grid_idx < n_grid_pts {
                    let wt = f64::from(values[grid_idx as usize]);
                    tot_wt += wt;
                    let p = shape_grid.get_grid_point_loc(grid_idx as u32)? - pt.location;
                    cov[0][0] += wt * (p.x * p.x);
                    cov[0][1] += wt * (p.x * p.y);
                    cov[0][2] += wt * (p.x * p.z);
                    cov[1][1] += wt * (p.y * p.y);
                    cov[1][2] += wt * (p.y * p.z);
                    cov[2][2] += wt * (p.z * p.z);
                }
            }
        }
    }
    for row in &mut cov {
        for x in row.iter_mut() {
            *x /= tot_wt;
        }
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];
    let (e_vals, e_vects) = symmetric_eigen(cov);
    let mut order = [0, 1, 2];
    order.sort_by(|&a, &b| e_vals[b].total_cmp(&e_vals[a]));
    let dir = |c: usize| Point3D::new(e_vects[0][c], e_vects[1][c], e_vects[2][c]);
    pt.shape_moments = order.map(|c| e_vals[c]);
    pt.shape_dirs = Some(order.map(dir));
    Ok(())
}

/// Adds to each point the families of the molecule's features within
/// `win_rad` of it. Mirrors `BuilderUtils.AssignMolFeatsToPoints`.
pub fn assign_mol_feats_to_points(
    pts: &mut [SkeletonPoint],
    mol: &Mol,
    feat_factory: &MolChemicalFeatureFactory,
    win_rad: f64,
) -> Result<()> {
    let mut feats = Vec::new();
    for i in 0..feat_factory.get_num_mol_features(mol)? {
        let feat = feat_factory.get_mol_feature(mol, i)?;
        feats.push((feat.get_pos()?, feat.get_family()));
    }
    for pt in pts.iter_mut() {
        for (pos, family) in &feats {
            if pos.distance(&pt.location) < win_rad && !pt.mol_features.contains(family) {
                pt.mol_features.push(family.clone());
            }
        }
    }
    Ok(())
}
