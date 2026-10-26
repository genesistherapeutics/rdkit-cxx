// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Subshape.SubshapeAligner`, a Python module:
//! aligning shapes by triangles of their skeleton points, pruned by
//! features, directions, and shape overlap.
//!
//! Python's generators return lists here, and its pruning statistics
//! (counts and the seconds each step took) come back from the
//! `_with_outputs` functions (A9). Its log messages are not ported.

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use crate::chem::rdchem::{Conformer, Mol};
use crate::chem::subshape::subshape_builder::SubshapeBuilder;
use crate::chem::subshape::subshape_objects::{ShapeWithSkeleton, SkeletonPoint};
use crate::geometry::rd_geometry::{Point3D, UniformGrid3D, protrude_distance, tanimoto_distance};
use crate::numerics::rd_alignment::get_alignment_transform;
use crate::{Error, Result};

/// Pruning statistics: how many alignments each step removed
/// (`"features"`, `"direction"`, `"coarseGrid"`, `"medGrid"`, `"fineGrid"`)
/// and the seconds steps took (`"gtm_time"`, ...), as Python's
/// `pruneStats` dictionary.
pub type PruneStats = BTreeMap<String, f64>;

fn count(stats: &mut PruneStats, key: &str) {
    *stats.entry(key.to_owned()).or_insert(0.0) += 1.0;
}

/// An alignment of a query shape to a target shape by a triangle of
/// skeleton points. Mirrors `SubshapeAligner.SubshapeAlignment`.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct SubshapeAlignment {
    /// The 4x4 transform of the query onto the target. Mirrors
    /// `SubshapeAlignment.transform`.
    pub transform: Vec<Vec<f64>>,
    /// The triangles' sum of squared deviations. Mirrors
    /// `SubshapeAlignment.triangleSSD`.
    pub triangle_ssd: f64,
    /// The target's skeleton points. Mirrors `SubshapeAlignment.targetTri`.
    pub target_tri: (usize, usize, usize),
    /// The query's skeleton points. Mirrors `SubshapeAlignment.queryTri`.
    pub query_tri: (usize, usize, usize),
    /// Mirrors `SubshapeAlignment.alignedConfId`.
    pub aligned_conf_id: i32,
    /// The sum of the points' directions' absolute dot products. Mirrors
    /// `SubshapeAlignment.dirMatch`.
    pub dir_match: f64,
    /// The shapes' distance at the finest grid. Mirrors
    /// `SubshapeAlignment.shapeDist`.
    pub shape_dist: f64,
}

/// The shape distance an aligner uses. Mirrors
/// `SubshapeAligner.SubshapeDistanceMetric`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum SubshapeDistanceMetric {
    /// The Tanimoto distance.
    Tanimoto,
    /// The protrude distance of the smaller shape from the larger.
    #[default]
    Protrude,
}

fn grid_of(shape: &ShapeWithSkeleton) -> Result<&UniformGrid3D> {
    shape.grid.as_ref().ok_or_else(|| Error::Value {
        message: "the shape has no grid".to_owned(),
    })
}

/// The distance between two shapes' grids. Mirrors
/// `SubshapeAligner.GetShapeShapeDistance`.
pub fn get_shape_shape_distance(
    s1: &ShapeWithSkeleton,
    s2: &ShapeWithSkeleton,
    dist_metric: SubshapeDistanceMetric,
) -> Result<f64> {
    let (g1, g2) = (grid_of(s1)?, grid_of(s2)?);
    match dist_metric {
        SubshapeDistanceMetric::Protrude => {
            if g1.get_occupancy_vect().get_total_val() < g2.get_occupancy_vect().get_total_val() {
                protrude_distance(g1, g2)
            } else {
                protrude_distance(g2, g1)
            }
        }
        SubshapeDistanceMetric::Tanimoto => tanimoto_distance(g1, g2),
    }
}

/// Python's `rdkit.ML.Cluster.Butina.ClusterData` of a lower-triangle
/// distance list: each cluster's centroid first.
#[allow(clippy::needless_range_loop)] // fills both triangles of the matrix
fn butina_cluster(dists: &[f64], n_pts: usize, dist_thresh: f64) -> Vec<Vec<usize>> {
    let mut matrix = vec![vec![0.0; n_pts]; n_pts];
    let mut k = 0;
    for i in 0..n_pts {
        for j in 0..i {
            matrix[i][j] = dists[k];
            matrix[j][i] = dists[k];
            k += 1;
        }
    }
    let neighbors: Vec<Vec<usize>> = matrix
        .iter()
        .map(|row| (0..n_pts).filter(|&j| row[j] <= dist_thresh).collect())
        .collect();
    let mut sorted: Vec<(usize, usize)> = neighbors
        .iter()
        .enumerate()
        .map(|(i, n)| (n.len(), i))
        .collect();
    sorted.sort_by(|a, b| b.cmp(a));
    let mut clusters = Vec::new();
    let mut seen = vec![false; n_pts];
    let mut pos = 0;
    while pos < sorted.len() && sorted[pos].0 > 1 {
        let idx = sorted[pos].1;
        pos += 1;
        if seen[idx] {
            continue;
        }
        let mut cluster = vec![idx];
        seen[idx] = true;
        for &nbr in &neighbors[idx] {
            if !seen[nbr] {
                cluster.push(nbr);
                seen[nbr] = true;
            }
        }
        clusters.push(cluster);
    }
    for &(_, idx) in &sorted[pos..] {
        if !seen[idx] {
            clusters.push(vec![idx]);
        }
    }
    clusters
}

/// Parameters for [`cluster_alignments_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClusterAlignmentsParams {
    /// Alignments whose shapes are this close are neighbors.
    pub neighbor_tol: f64,
    /// The shape distance.
    pub dist_metric: SubshapeDistanceMetric,
    /// The conformer ids the transformed molecule takes while it is
    /// compared (this and the next).
    pub temp_conf_id: i32,
}

impl Default for ClusterAlignmentsParams {
    fn default() -> ClusterAlignmentsParams {
        ClusterAlignmentsParams {
            neighbor_tol: 0.1,
            dist_metric: SubshapeDistanceMetric::Protrude,
            temp_conf_id: 1001,
        }
    }
}

/// The alignments at the centers of Butina clusters of the shapes they
/// give the molecule. Mirrors `SubshapeAligner.ClusterAlignments`.
pub fn cluster_alignments(
    mol: &mut Mol,
    alignments: &[SubshapeAlignment],
    builder: &SubshapeBuilder,
) -> Result<Vec<SubshapeAlignment>> {
    cluster_alignments_with_params(mol, alignments, builder, &Default::default())
}

/// Mirrors `SubshapeAligner.ClusterAlignments(mol, alignments, builder,
/// neighborTol, distMetric, tempConfId)`.
pub fn cluster_alignments_with_params(
    mol: &mut Mol,
    alignments: &[SubshapeAlignment],
    builder: &SubshapeBuilder,
    params: &ClusterAlignmentsParams,
) -> Result<Vec<SubshapeAlignment>> {
    let temp = params.temp_conf_id;
    let mut dists = Vec::new();
    for i in 0..alignments.len() {
        transform_mol_with_params(
            mol,
            &alignments[i].transform,
            &TransformMolParams {
                conf_id: -1,
                new_conf_id: temp,
            },
        )?;
        let shape_i = builder.shape_with_spacing(mol, temp, builder.grid_spacing)?;
        for alignment in &alignments[..i] {
            let params_j = TransformMolParams {
                conf_id: -1,
                new_conf_id: temp + 1,
            };
            transform_mol_with_params(mol, &alignment.transform, &params_j)?;
            let shape_j = builder.shape_with_spacing(mol, temp + 1, builder.grid_spacing)?;
            dists.push(get_shape_shape_distance(
                &shape_i,
                &shape_j,
                params.dist_metric,
            )?);
            mol.remove_conformer((temp + 1) as u32);
        }
        mol.remove_conformer(temp as u32);
    }
    let clusters = butina_cluster(&dists, alignments.len(), params.neighbor_tol);
    Ok(clusters.iter().map(|c| alignments[c[0]].clone()).collect())
}

/// Parameters for [`transform_mol_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransformMolParams {
    /// The conformer to transform; −1 for the default.
    pub conf_id: i32,
    /// The id of the transformed conformer, which replaces any with it.
    pub new_conf_id: i32,
}

impl Default for TransformMolParams {
    fn default() -> TransformMolParams {
        TransformMolParams {
            conf_id: -1,
            new_conf_id: 100,
        }
    }
}

/// Adds the conformer transformed by the 4x4 `tform` as conformer 100,
/// replacing any. NumPy's matrix product may round differently in the last
/// place. Mirrors `SubshapeAligner.TransformMol`.
pub fn transform_mol(mol: &mut Mol, tform: &[Vec<f64>]) -> Result<()> {
    transform_mol_with_params(mol, tform, &Default::default())
}

/// Mirrors `SubshapeAligner.TransformMol(mol, tform, confId, newConfId)`.
/// A transform that is not 4x4 is `Error::Value`, and a negative
/// `new_conf_id` too.
pub fn transform_mol_with_params(
    mol: &mut Mol,
    tform: &[Vec<f64>],
    params: &TransformMolParams,
) -> Result<()> {
    if tform.len() < 3 || tform.iter().take(3).any(|row| row.len() != 4) {
        return Err(Error::Value {
            message: "the transform must be 4x4".to_owned(),
        });
    }
    let new_conf_id = u32::try_from(params.new_conf_id).map_err(|_| Error::Value {
        message: format!("negative conformer id {}", params.new_conf_id),
    })?;
    let ref_conf = mol.get_conformer_with_id(params.conf_id)?;
    let n = ref_conf.get_num_atoms();
    let mut new_conf = Conformer::new(n);
    for i in 0..n {
        let p = ref_conf.get_atom_position(i)?;
        let pos = [p.x, p.y, p.z, 1.0];
        let row = |r: &Vec<f64>| {
            r.iter()
                .zip(pos)
                .map(|(a, b)| a * b)
                .fold(0.0, |acc, x| acc + x)
        };
        new_conf.set_atom_position(
            i,
            Point3D::new(row(&tform[0]), row(&tform[1]), row(&tform[2])),
        )?;
    }
    new_conf.set_id(new_conf_id);
    mol.remove_conformer(new_conf_id);
    mol.add_conformer(&new_conf)?;
    Ok(())
}

/// The triangles of point indices: (i, j, k) with i < j < k when ordered;
/// otherwise j may be below i, and k is above j.
fn all_triangles(n: usize, ordered: bool) -> Vec<(usize, usize, usize)> {
    let mut res = Vec::new();
    for i in 0..n {
        let j_start = if ordered { i + 1 } else { 0 };
        for j in j_start..n {
            if j == i {
                continue;
            }
            for k in j + 1..n {
                if k == i {
                    continue;
                }
                res.push((i, j, k));
            }
        }
    }
    res
}

/// Parameters for [`SubshapeAligner`]'s calls that take conformers, with
/// Python's defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlignerConfParams {
    /// The target's conformer; −1 for the default.
    pub tgt_conf: i32,
    /// The query's conformer; −1 for the default.
    pub query_conf: i32,
}

impl Default for AlignerConfParams {
    fn default() -> AlignerConfParams {
        AlignerConfParams {
            tgt_conf: -1,
            query_conf: -1,
        }
    }
}

/// Finds alignments of shapes. Python's class attributes are its fields.
/// Mirrors `SubshapeAligner.SubshapeAligner`.
#[derive(Clone, Debug, PartialEq)]
pub struct SubshapeAligner {
    /// The largest RMS deviation of matched triangles. Mirrors
    /// `SubshapeAligner.triangleRMSTol`.
    pub triangle_rms_tol: f64,
    /// The shape distance. Mirrors `SubshapeAligner.distMetric`.
    pub dist_metric: SubshapeDistanceMetric,
    /// The largest shape distance. Mirrors `SubshapeAligner.shapeDistTol`.
    pub shape_dist_tol: f64,
    /// The fewest points whose features match. Mirrors
    /// `SubshapeAligner.numFeatThresh`.
    pub num_feat_thresh: usize,
    /// The smallest sum of the directions' dot products. Mirrors
    /// `SubshapeAligner.dirThresh`.
    pub dir_thresh: f64,
    /// The largest difference of matched edges' squared lengths, as a
    /// length. Mirrors `SubshapeAligner.edgeTol`.
    pub edge_tol: f64,
    /// Scales the shape tolerance on the coarse grid. Mirrors
    /// `SubshapeAligner.coarseGridToleranceMult`.
    pub coarse_grid_tolerance_mult: f64,
    /// Scales the shape tolerance on the medium grid. Mirrors
    /// `SubshapeAligner.medGridToleranceMult`.
    pub med_grid_tolerance_mult: f64,
}

impl Default for SubshapeAligner {
    fn default() -> SubshapeAligner {
        SubshapeAligner {
            triangle_rms_tol: 1.0,
            dist_metric: SubshapeDistanceMetric::Protrude,
            shape_dist_tol: 0.2,
            num_feat_thresh: 3,
            dir_thresh: 2.6,
            edge_tol: 6.0,
            coarse_grid_tolerance_mult: 1.0,
            med_grid_tolerance_mult: 1.0,
        }
    }
}

fn tri(t: (usize, usize, usize)) -> [usize; 3] {
    [t.0, t.1, t.2]
}

fn point(pts: &[SkeletonPoint], i: usize) -> Result<&SkeletonPoint> {
    pts.get(i).ok_or(Error::Index {
        index: i as u32,
        len: pts.len() as u32,
    })
}

impl SubshapeAligner {
    /// The alignments of triangles of the query's skeleton points onto the
    /// target's whose edges have compatible lengths and whose points
    /// superimpose within the RMS tolerance. Mirrors
    /// `SubshapeAligner.GetTriangleMatches`.
    pub fn get_triangle_matches(
        &self,
        target: &ShapeWithSkeleton,
        query: &ShapeWithSkeleton,
    ) -> Result<Vec<SubshapeAlignment>> {
        let ssd_tol = self.triangle_rms_tol.powf(2.0) * 9.0;
        let (tgt_pts, query_pts) = (&target.skel_pts, &query.skel_pts);
        let lengths = |pts: &[SkeletonPoint]| {
            let mut ls = Vec::new();
            for i in 0..pts.len() {
                for j in i + 1..pts.len() {
                    ls.push(((i, j), (pts[i].location - pts[j].location).length_sq()));
                }
            }
            ls
        };
        let tol2 = self.edge_tol * self.edge_tol;
        let mut compat_edges = BTreeSet::new();
        for (tk, tv) in lengths(tgt_pts) {
            for (qk, qv) in lengths(query_pts) {
                if (tv - qv).abs() < tol2 {
                    compat_edges.insert((tk, qk));
                }
            }
        }
        let mut res = Vec::new();
        for tgt_tri in all_triangles(tgt_pts.len(), true) {
            let tgt_locs: Vec<Point3D> =
                tri(tgt_tri).iter().map(|&x| tgt_pts[x].location).collect();
            for query_tri in all_triangles(query_pts.len(), false) {
                let (t, q) = (tgt_tri, query_tri);
                if compat_edges.contains(&((t.0, t.1), (q.0, q.1)))
                    && compat_edges.contains(&((t.0, t.2), (q.0, q.2)))
                    && compat_edges.contains(&((t.1, t.2), (q.1, q.2)))
                {
                    let query_locs: Vec<Point3D> =
                        tri(q).iter().map(|&x| query_pts[x].location).collect();
                    let (ssd, tf) = get_alignment_transform(&tgt_locs, &query_locs)?;
                    if ssd <= ssd_tol {
                        res.push(SubshapeAlignment {
                            transform: tf,
                            triangle_ssd: ssd,
                            target_tri: t,
                            query_tri: q,
                            aligned_conf_id: -1,
                            dir_match: 0.0,
                            shape_dist: 0.0,
                        });
                    }
                }
            }
        }
        Ok(res)
    }

    fn check_match_features(
        &self,
        target_pts: &[SkeletonPoint],
        query_pts: &[SkeletonPoint],
        alignment: &SubshapeAlignment,
    ) -> Result<bool> {
        let mut n_matched = 0;
        for (t, q) in tri(alignment.target_tri)
            .into_iter()
            .zip(tri(alignment.query_tri))
        {
            let tgt_feats = &point(target_pts, t)?.mol_features;
            let q_feats = &point(query_pts, q)?.mol_features;
            if (tgt_feats.is_empty() && q_feats.is_empty())
                || tgt_feats.iter().any(|f| q_feats.contains(f))
            {
                n_matched += 1;
            }
            if n_matched >= self.num_feat_thresh {
                break;
            }
        }
        Ok(n_matched >= self.num_feat_thresh)
    }

    /// Removes the alignments whose points' features do not match: points
    /// match when neither has features or they share one. Mirrors
    /// `SubshapeAligner.PruneMatchesUsingFeatures`.
    pub fn prune_matches_using_features(
        &self,
        target: &ShapeWithSkeleton,
        query: &ShapeWithSkeleton,
        alignments: &mut Vec<SubshapeAlignment>,
    ) -> Result<()> {
        self.prune_matches_using_features_with_outputs(target, query, alignments)
            .map(|_| ())
    }

    /// With the number removed, under `"features"` (A9). Mirrors
    /// `SubshapeAligner.PruneMatchesUsingFeatures(..., pruneStats)`.
    pub fn prune_matches_using_features_with_outputs(
        &self,
        target: &ShapeWithSkeleton,
        query: &ShapeWithSkeleton,
        alignments: &mut Vec<SubshapeAlignment>,
    ) -> Result<PruneStats> {
        let mut stats = PruneStats::new();
        let mut i = 0;
        while i < alignments.len() {
            if self.check_match_features(&target.skel_pts, &query.skel_pts, &alignments[i])? {
                i += 1;
            } else {
                count(&mut stats, "features");
                alignments.remove(i);
            }
        }
        Ok(stats)
    }

    fn check_match_directions(
        &self,
        target_pts: &[SkeletonPoint],
        query_pts: &[SkeletonPoint],
        alignment: &mut SubshapeAlignment,
    ) -> Result<bool> {
        let no_dirs = || Error::Value {
            message: "the skeleton point has no directions".to_owned(),
        };
        let t = &alignment.transform;
        let mut dot = 0.0;
        for (ti, qi) in tri(alignment.target_tri)
            .into_iter()
            .zip(tri(alignment.query_tri))
        {
            let qv = point(query_pts, qi)?.shape_dirs.ok_or_else(no_dirs)?[0];
            let tv = point(target_pts, ti)?.shape_dirs.ok_or_else(no_dirs)?[0];
            let rot = |r: usize| t[r][0] * qv.x + t[r][1] * qv.y + t[r][2] * qv.z;
            let rot_v = [rot(0), rot(1), rot(2)];
            dot += (rot_v[0] * tv.x + rot_v[1] * tv.y + rot_v[2] * tv.z).abs();
            if dot >= self.dir_thresh {
                break;
            }
        }
        alignment.dir_match = dot;
        Ok(dot >= self.dir_thresh)
    }

    /// Removes the alignments whose points' principal directions do not
    /// align, and sets the others' `dir_match`. Mirrors
    /// `SubshapeAligner.PruneMatchesUsingDirection`.
    pub fn prune_matches_using_direction(
        &self,
        target: &ShapeWithSkeleton,
        query: &ShapeWithSkeleton,
        alignments: &mut Vec<SubshapeAlignment>,
    ) -> Result<()> {
        self.prune_matches_using_direction_with_outputs(target, query, alignments)
            .map(|_| ())
    }

    /// With the number removed, under `"direction"` (A9). Mirrors
    /// `SubshapeAligner.PruneMatchesUsingDirection(..., pruneStats)`.
    pub fn prune_matches_using_direction_with_outputs(
        &self,
        target: &ShapeWithSkeleton,
        query: &ShapeWithSkeleton,
        alignments: &mut Vec<SubshapeAlignment>,
    ) -> Result<PruneStats> {
        let mut stats = PruneStats::new();
        let mut i = 0;
        while i < alignments.len() {
            if self.check_match_directions(&target.skel_pts, &query.skel_pts, &mut alignments[i])? {
                i += 1;
            } else {
                count(&mut stats, "direction");
                alignments.remove(i);
            }
        }
        Ok(stats)
    }

    /// Gives the target its medium (1.5 times the spacing) and coarse (2
    /// times) grids, from the molecule or, without one, resampled.
    fn add_coarse_and_medium_grids(
        &self,
        mol: Option<&Mol>,
        tgt: &mut ShapeWithSkeleton,
        conf_id: i32,
        builder: &SubshapeBuilder,
    ) -> Result<()> {
        let o_space = builder.grid_spacing;
        let (med, coarse) = match mol {
            Some(mol) => (
                builder.shape_with_spacing(mol, conf_id, o_space * 1.5)?,
                builder.shape_with_spacing(mol, conf_id, o_space * 2.0)?,
            ),
            None => (
                builder.sample_subshape(tgt, o_space * 1.5)?,
                builder.sample_subshape(tgt, o_space * 2.0)?,
            ),
        };
        tgt.med_grid = Some(Box::new(med));
        tgt.coarse_grid = Some(Box::new(coarse));
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn check_match_shape(
        &self,
        target: &ShapeWithSkeleton,
        query_mol: &mut Mol,
        alignment: &mut SubshapeAlignment,
        builder: &SubshapeBuilder,
        query_conf: i32,
        stats: &mut PruneStats,
    ) -> Result<bool> {
        let t_conf_id = 1001;
        let params = TransformMolParams {
            conf_id: query_conf,
            new_conf_id: t_conf_id,
        };
        transform_mol_with_params(query_mol, &alignment.transform, &params)?;
        let o_space = builder.grid_spacing;
        let missing = || Error::Value {
            message: "the target has no coarse grids".to_owned(),
        };
        let result = (|| -> Result<bool> {
            let coarse = builder.shape_with_spacing(query_mol, t_conf_id, o_space * 2.0)?;
            let target_coarse = target.coarse_grid.as_deref().ok_or_else(missing)?;
            let d = get_shape_shape_distance(&coarse, target_coarse, self.dist_metric)?;
            if d > self.shape_dist_tol * self.coarse_grid_tolerance_mult {
                count(stats, "coarseGrid");
                return Ok(false);
            }
            let med = builder.shape_with_spacing(query_mol, t_conf_id, o_space * 1.5)?;
            let target_med = target.med_grid.as_deref().ok_or_else(missing)?;
            let d = get_shape_shape_distance(&med, target_med, self.dist_metric)?;
            if d > self.shape_dist_tol * self.med_grid_tolerance_mult {
                count(stats, "medGrid");
                return Ok(false);
            }
            let fine = builder.shape_with_spacing(query_mol, t_conf_id, o_space)?;
            let d = get_shape_shape_distance(&fine, target, self.dist_metric)?;
            let ok = d <= self.shape_dist_tol;
            if !ok {
                count(stats, "fineGrid");
            }
            alignment.shape_dist = d;
            Ok(ok)
        })();
        query_mol.remove_conformer(t_conf_id as u32);
        result
    }

    /// Removes the alignments whose transformed query shape is too far from
    /// the target's at the coarse, medium, and fine grids, and sets the
    /// others' `shape_dist`. Gives the target its coarse grids first if it
    /// has none. Mirrors `SubshapeAligner.PruneMatchesUsingShape`.
    pub fn prune_matches_using_shape(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        builder: &SubshapeBuilder,
        alignments: &mut Vec<SubshapeAlignment>,
    ) -> Result<()> {
        self.prune_matches_using_shape_with_params(
            target_mol,
            target,
            query_mol,
            builder,
            alignments,
            &Default::default(),
        )
    }

    /// Mirrors `SubshapeAligner.PruneMatchesUsingShape(targetMol, target,
    /// queryMol, query, builder, alignments, tgtConf, queryConf)`.
    pub fn prune_matches_using_shape_with_params(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        builder: &SubshapeBuilder,
        alignments: &mut Vec<SubshapeAlignment>,
        params: &AlignerConfParams,
    ) -> Result<()> {
        self.prune_matches_using_shape_with_outputs(
            target_mol, target, query_mol, builder, alignments, params,
        )
        .map(|_| ())
    }

    /// With the numbers removed at each grid (A9). Mirrors
    /// `SubshapeAligner.PruneMatchesUsingShape(..., pruneStats)`.
    pub fn prune_matches_using_shape_with_outputs(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        builder: &SubshapeBuilder,
        alignments: &mut Vec<SubshapeAlignment>,
        params: &AlignerConfParams,
    ) -> Result<PruneStats> {
        if target.med_grid.is_none() {
            self.add_coarse_and_medium_grids(Some(target_mol), target, params.tgt_conf, builder)?;
        }
        let mut stats = PruneStats::new();
        let mut i = 0;
        while i < alignments.len() {
            if self.check_match_shape(
                target,
                query_mol,
                &mut alignments[i],
                builder,
                params.query_conf,
                &mut stats,
            )? {
                i += 1;
            } else {
                alignments.remove(i);
            }
        }
        Ok(stats)
    }

    /// The query's alignments onto the target that survive pruning by
    /// features (with a feature factory), directions, and shape. Mirrors
    /// `SubshapeAligner.GetSubshapeAlignments`.
    pub fn get_subshape_alignments(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        query: &ShapeWithSkeleton,
        builder: &SubshapeBuilder,
    ) -> Result<Vec<SubshapeAlignment>> {
        self.get_subshape_alignments_with_params(
            target_mol,
            target,
            query_mol,
            query,
            builder,
            &Default::default(),
        )
    }

    /// Mirrors `SubshapeAligner.GetSubshapeAlignments(targetMol, target,
    /// queryMol, query, builder, tgtConf, queryConf)`.
    pub fn get_subshape_alignments_with_params(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        query: &ShapeWithSkeleton,
        builder: &SubshapeBuilder,
        params: &AlignerConfParams,
    ) -> Result<Vec<SubshapeAlignment>> {
        Ok(self
            .get_subshape_alignments_with_outputs(
                target_mol, target, query_mol, query, builder, params,
            )?
            .0)
    }

    /// With the pruning statistics, including each step's seconds (A9).
    /// Mirrors `SubshapeAligner.GetSubshapeAlignments(..., pruneStats)`.
    pub fn get_subshape_alignments_with_outputs(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        query: &ShapeWithSkeleton,
        builder: &SubshapeBuilder,
        params: &AlignerConfParams,
    ) -> Result<(Vec<SubshapeAlignment>, PruneStats)> {
        let mut stats = PruneStats::new();
        let merge = |stats: &mut PruneStats, more: PruneStats| {
            for (k, v) in more {
                *stats.entry(k).or_insert(0.0) += v;
            }
        };
        let t1 = Instant::now();
        let mut res = self.get_triangle_matches(target, query)?;
        stats.insert("gtm_time".to_owned(), t1.elapsed().as_secs_f64());
        if builder.feat_factory.is_some() {
            let t1 = Instant::now();
            let more = self.prune_matches_using_features_with_outputs(target, query, &mut res)?;
            merge(&mut stats, more);
            stats.insert("feats_time".to_owned(), t1.elapsed().as_secs_f64());
        }
        let t1 = Instant::now();
        let more = self.prune_matches_using_direction_with_outputs(target, query, &mut res)?;
        merge(&mut stats, more);
        stats.insert("direction_time".to_owned(), t1.elapsed().as_secs_f64());
        let t1 = Instant::now();
        let more = self.prune_matches_using_shape_with_outputs(
            target_mol, target, query_mol, builder, &mut res, params,
        )?;
        merge(&mut stats, more);
        stats.insert("shape_time".to_owned(), t1.elapsed().as_secs_f64());
        Ok((res, stats))
    }

    /// The alignments that pass each check in turn, as
    /// [`get_subshape_alignments`](SubshapeAligner::get_subshape_alignments)
    /// finds them. Mirrors `SubshapeAligner.__call__`.
    pub fn call(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        query: &ShapeWithSkeleton,
        builder: &SubshapeBuilder,
    ) -> Result<Vec<SubshapeAlignment>> {
        self.call_with_params(
            target_mol,
            target,
            query_mol,
            query,
            builder,
            &Default::default(),
        )
    }

    /// Mirrors `SubshapeAligner.__call__(targetMol, target, queryMol, query,
    /// builder, tgtConf, queryConf)`.
    pub fn call_with_params(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        query: &ShapeWithSkeleton,
        builder: &SubshapeBuilder,
        params: &AlignerConfParams,
    ) -> Result<Vec<SubshapeAlignment>> {
        Ok(self
            .call_with_outputs(target_mol, target, query_mol, query, builder, params)?
            .0)
    }

    /// With the pruning statistics (A9). Mirrors
    /// `SubshapeAligner.__call__(..., pruneStats)`.
    pub fn call_with_outputs(
        &self,
        target_mol: &Mol,
        target: &mut ShapeWithSkeleton,
        query_mol: &mut Mol,
        query: &ShapeWithSkeleton,
        builder: &SubshapeBuilder,
        params: &AlignerConfParams,
    ) -> Result<(Vec<SubshapeAlignment>, PruneStats)> {
        let mut stats = PruneStats::new();
        let mut res = Vec::new();
        for mut alignment in self.get_triangle_matches(target, query)? {
            if !self.check_match_features(&target.skel_pts, &query.skel_pts, &alignment)?
                && builder.feat_factory.is_some()
            {
                count(&mut stats, "features");
                continue;
            }
            if !self.check_match_directions(&target.skel_pts, &query.skel_pts, &mut alignment)? {
                count(&mut stats, "direction");
                continue;
            }
            if target.med_grid.is_none() {
                self.add_coarse_and_medium_grids(
                    Some(target_mol),
                    target,
                    params.tgt_conf,
                    builder,
                )?;
            }
            if self.check_match_shape(
                target,
                query_mol,
                &mut alignment,
                builder,
                params.query_conf,
                &mut stats,
            )? {
                res.push(alignment);
            }
        }
        Ok((res, stats))
    }
}
