// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FeatMaps.FeatMapUtils`, a Python module: merging
//! and combining feature maps.
//!
//! Python's mode classes hold integer constants, and their `valid` class
//! methods check that a value is one of them; here they are enums, whose
//! values all are.

use crate::Result;
use crate::chem::feat_maps::feat_map_point::{FeatMapPoint, FeatureLike};
use crate::chem::feat_maps::feat_maps::{FeatMap, FeatMapInitParams};

/// How two merged points' positions and weights combine. Mirrors
/// `FeatMapUtils.MergeMethod`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum MergeMethod {
    /// The weighted mean of the positions, and the mean weight.
    #[default]
    WeightedAverage,
    /// The mean of the positions and of the weights.
    Average,
    /// The position and weight of the point with the larger weight.
    UseLarger,
}

/// What decides which points merge. Mirrors `FeatMapUtils.MergeMetric`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum MergeMetric {
    /// No points merge.
    #[default]
    NoMerge,
    /// Points closer than the tolerance merge.
    Distance,
    /// Points whose overlap score is above the tolerance merge.
    Overlap,
}

/// How merged points' directions combine. Python accepts it and does not
/// use it. Mirrors `FeatMapUtils.DirMergeMode`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum DirMergeMode {
    /// Directions are not merged.
    #[default]
    NoMerge,
    /// Directions are added.
    Sum,
}

/// The matrix of the map's point pairs' merge metrics: the squared
/// distance for [`MergeMetric::Distance`] below the tolerance squared, the
/// negated overlap score (times the second point's weight) for
/// [`MergeMetric::Overlap`] below the tolerance, and 1e8 otherwise, or
/// for pairs `compat_func` rejects. Mirrors
/// `FeatMapUtils.GetFeatFeatDistMatrix`.
pub fn get_feat_feat_dist_matrix(
    fm: &FeatMap,
    merge_metric: MergeMetric,
    merge_tol: f64,
    _dir_merge_mode: DirMergeMode,
    compat_func: impl Fn(&FeatMapPoint, &FeatMapPoint) -> bool,
) -> Result<Vec<Vec<f64>>> {
    let n = fm.get_num_features();
    let mut dists = vec![vec![1e8; n]; n];
    let benchmark = match merge_metric {
        MergeMetric::NoMerge => return Ok(dists),
        MergeMetric::Distance => merge_tol * merge_tol,
        MergeMetric::Overlap => merge_tol,
    };
    let feats = fm.get_features();
    for i in 0..n {
        for j in i + 1..n {
            let (pt_i, pt_j) = (&feats[i], &feats[j]);
            if !compat_func(pt_i, pt_j) {
                continue;
            }
            let value = if merge_metric == MergeMetric::Distance {
                pt_i.get_dist2(pt_j)?
            } else {
                fm.get_feat_feat_score_with_type_match(pt_i, pt_j, false)? * (-pt_j.weight)
            };
            if value < benchmark {
                dists[i][j] = value;
                dists[j][i] = value;
            }
        }
    }
    Ok(dists)
}

/// Whether two features have the same family. Mirrors
/// `FeatMapUtils.familiesMatch`.
pub fn families_match(f1: &dyn FeatureLike, f2: &dyn FeatureLike) -> bool {
    f1.feature_family() == f2.feature_family()
}

/// Whether two values are within 1e-4. Mirrors `FeatMapUtils.feq`.
pub fn feq(v1: f64, v2: f64) -> bool {
    feq_with_tol(v1, v2, 1e-4)
}

/// Within `tol`. Mirrors `FeatMapUtils.feq(v1, v2, tol)`.
pub fn feq_with_tol(v1: f64, v2: f64, tol: f64) -> bool {
    (v1 - v2).abs() < tol
}

fn points_match(f1: &FeatMapPoint, f2: &FeatMapPoint) -> bool {
    families_match(f1, f2)
}

/// Parameters for [`merge_feat_points_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Copy, Debug)]
pub struct MergeFeatPointsParams {
    /// What decides which points merge.
    pub merge_metric: MergeMetric,
    /// The largest distance for distance-based merging, or the smallest
    /// score for score-based merging.
    pub merge_tol: f64,
    /// How directions merge (unused, as in Python).
    pub dir_merge_mode: DirMergeMode,
    /// How positions and weights merge.
    pub merge_method: MergeMethod,
    /// Which points may merge; by default, those of the same family.
    pub compat_func: fn(&FeatMapPoint, &FeatMapPoint) -> bool,
}

impl Default for MergeFeatPointsParams {
    fn default() -> MergeFeatPointsParams {
        MergeFeatPointsParams {
            merge_metric: MergeMetric::NoMerge,
            merge_tol: 1.5,
            dir_merge_mode: DirMergeMode::NoMerge,
            merge_method: MergeMethod::WeightedAverage,
            compat_func: points_match,
        }
    }
}

/// Merges the map's points that are each other's nearest under the merge
/// metric, pair by pair, and returns whether any merged. With the default
/// [`MergeMetric::NoMerge`], none do. Mirrors `FeatMapUtils.MergeFeatPoints`.
pub fn merge_feat_points(fm: &mut FeatMap) -> Result<bool> {
    merge_feat_points_with_params(fm, &MergeFeatPointsParams::default())
}

/// Mirrors `FeatMapUtils.MergeFeatPoints(fm, mergeMetric, mergeTol,
/// dirMergeMode, mergeMethod, compatFunc)`.
pub fn merge_feat_points_with_params(
    fm: &mut FeatMap,
    params: &MergeFeatPointsParams,
) -> Result<bool> {
    let mut res = false;
    if params.merge_metric == MergeMetric::NoMerge {
        return Ok(res);
    }
    let dists = get_feat_feat_dist_matrix(
        fm,
        params.merge_metric,
        params.merge_tol,
        params.dir_merge_mode,
        params.compat_func,
    )?;
    // Python sorts (distance, index) tuples; the distances are never NaN.
    let dist_orders: Vec<Vec<(f64, usize)>> = dists
        .iter()
        .map(|dist_v| {
            let mut order: Vec<(f64, usize)> = dist_v
                .iter()
                .enumerate()
                .filter(|&(_, &d)| d < params.merge_tol)
                .map(|(j, &d)| (d, j))
                .collect();
            order.sort_by(|a, b| a.partial_cmp(b).expect("the distances are not NaN"));
            order
        })
        .collect();
    // Python then tries to remove merged indices from these lists of
    // tuples, which never holds them, so the lists stay as they are.
    let mut feats_in_play: Vec<usize> = (0..fm.get_num_features()).collect();
    let mut feats_to_remove = Vec::new();
    while !feats_in_play.is_empty() {
        let mut merge = None;
        for fi in feats_in_play.clone() {
            let Some(&(dist, nbr)) = dist_orders[fi].first() else {
                feats_in_play.retain(|&f| f != fi);
                continue;
            };
            if !feats_in_play.contains(&nbr) {
                continue;
            }
            let nbr_order = &dist_orders[nbr];
            let mut merge_them = nbr_order[0].1 == fi;
            if !merge_them && feq(nbr_order[0].0, dist) {
                for &(dist_j, nbr_j) in &nbr_order[1..] {
                    if !feq(dist, dist_j) {
                        break;
                    }
                    if nbr_j == fi {
                        merge_them = true;
                        break;
                    }
                }
            }
            if merge_them {
                merge = Some((fi, nbr));
                break;
            }
        }
        let Some((fi, nbr)) = merge else {
            break;
        };
        res = true;
        let nbr_feat = fm.get_feature(nbr)?.clone();
        let feat_i = fm.get_feature_mut(fi)?;
        let (pos_i, pos_n) = (feat_i.get_pos(), nbr_feat.get_pos());
        let (w_i, w_n) = (feat_i.weight, nbr_feat.weight);
        let (new_pos, new_weight) = match params.merge_method {
            MergeMethod::WeightedAverage => {
                ((pos_i * w_i + pos_n * w_n) / (w_i + w_n), (w_i + w_n) / 2.0)
            }
            MergeMethod::Average => ((pos_i + pos_n) / 2.0, (w_i + w_n) / 2.0),
            MergeMethod::UseLarger => {
                if w_i > w_n {
                    (pos_i, w_i)
                } else {
                    (pos_n, w_n)
                }
            }
        };
        feat_i.set_pos(new_pos);
        feat_i.weight = new_weight;
        feats_to_remove.push(nbr);
        feats_in_play.retain(|&f| f != fi && f != nbr);
    }
    feats_to_remove.sort();
    for (i, f_idx) in feats_to_remove.into_iter().enumerate() {
        fm.drop_feature(f_idx - i)?;
    }
    Ok(res)
}

/// Parameters for [`combine_feat_maps_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CombineFeatMapsParams {
    /// What decides which points merge.
    pub merge_metric: MergeMetric,
    /// As in [`MergeFeatPointsParams`].
    pub merge_tol: f64,
    /// Unused, as in Python.
    pub dir_merge_mode: DirMergeMode,
}

impl Default for CombineFeatMapsParams {
    fn default() -> CombineFeatMapsParams {
        CombineFeatMapsParams {
            merge_metric: MergeMetric::NoMerge,
            merge_tol: 1.5,
            dir_merge_mode: DirMergeMode::NoMerge,
        }
    }
}

/// A map with both maps' points, and `fm1`'s parameters. Mirrors
/// `FeatMapUtils.CombineFeatMaps`.
pub fn combine_feat_maps(fm1: &FeatMap, fm2: &FeatMap) -> Result<FeatMap> {
    combine_feat_maps_with_params(fm1, fm2, &CombineFeatMapsParams::default())
}

/// Merging the points as [`merge_feat_points_with_params`] does, by the
/// default method. Mirrors `FeatMapUtils.CombineFeatMaps(fm1, fm2,
/// mergeMetric, mergeTol, dirMergeMode)`.
pub fn combine_feat_maps_with_params(
    fm1: &FeatMap,
    fm2: &FeatMap,
    params: &CombineFeatMapsParams,
) -> Result<FeatMap> {
    let mut res = FeatMap::new_with_params(&FeatMapInitParams {
        params: fm1.params.clone(),
        ..Default::default()
    })?;
    for feat in fm1.get_features().iter().chain(fm2.get_features()) {
        res.add_feat_point(feat.clone())?;
    }
    if params.merge_metric != MergeMetric::NoMerge {
        let merge = MergeFeatPointsParams {
            merge_metric: params.merge_metric,
            merge_tol: params.merge_tol,
            ..Default::default()
        };
        merge_feat_points_with_params(&mut res, &merge)?;
    }
    Ok(res)
}
