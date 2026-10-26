// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FeatMaps.FeatMaps`, a Python module: feature maps
//! and their scoring.
//!
//! Python's mode classes hold integer constants; here they are enums, and
//! Python's class-level defaults are those of [`FeatMap::default`].

use std::collections::BTreeMap;
use std::fmt;

use crate::chem::feat_maps::feat_map_point::{FeatMapPoint, FeatureLike};
use crate::chem::py_math::py_sum;
use crate::{Error, Result};

/// How [`FeatMap::score_feats`] scores each feature. Mirrors
/// `FeatMaps.FeatMapScoreMode`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum FeatMapScoreMode {
    /// Against every matching feature of the map.
    #[default]
    All,
    /// Against the closest matching feature of the map.
    Closest,
    /// Against the matching feature of the map that scores highest.
    Best,
}

/// How feature directions weigh in a score. Mirrors
/// `FeatMaps.FeatDirScoreMode`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum FeatDirScoreMode {
    /// Directions are ignored.
    #[default]
    Ignore,
    /// The dot product of the directions, negative when they are
    /// antiparallel.
    DotFullRange,
    /// The dot product of the directions scaled to 0 to 1.
    DotPosRange,
}

/// A feature's scoring profile. Mirrors `FeatMaps.FeatMapParams.FeatProfile`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum FeatProfile {
    /// A Gaussian of width `width`.
    #[default]
    Gaussian,
    /// A linear fall-off to 0 at `width`.
    Triangle,
    /// 1 within the radius.
    Box,
}

/// The scoring parameters of one feature family in a map. Mirrors
/// `FeatMaps.FeatMapParams`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FeatMapParams {
    /// The cutoff radius. Mirrors `FeatMapParams.radius`.
    pub radius: f64,
    /// The width of the profile (a Gaussian's sigma). Mirrors
    /// `FeatMapParams.width`.
    pub width: f64,
    /// The scoring profile. Mirrors `FeatMapParams.featProfile`.
    pub feat_profile: FeatProfile,
}

impl Default for FeatMapParams {
    fn default() -> FeatMapParams {
        FeatMapParams {
            radius: 2.5,
            width: 1.0,
            feat_profile: FeatProfile::Gaussian,
        }
    }
}

/// Parameters for [`FeatMap::new_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default)]
pub struct FeatMapInitParams {
    /// The scoring parameters, by family.
    pub params: BTreeMap<String, FeatMapParams>,
    /// The points to add, with the weights of `weights`.
    pub feats: Vec<FeatMapPoint>,
    /// One weight per point.
    pub weights: Vec<f64>,
}

/// The scores [`FeatMap::score_feats_with_outputs`] fills in, as Python's
/// output arguments (A9).
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct ScoreFeatsOutputs {
    /// Each map feature's score. Mirrors `mapScoreVect`.
    pub map_score_vect: Vec<f64>,
    /// Each scored feature's score. Mirrors `featsScoreVect`.
    pub feats_score_vect: Vec<f64>,
    /// The map features each scored feature was scored against. Mirrors
    /// `featsToFeatMapIdx`.
    pub feats_to_feat_map_idx: Vec<Vec<usize>>,
}

/// A feature map: weighted points and how they score other features.
/// Mirrors `FeatMaps.FeatMap`.
#[derive(Clone, Debug, Default)]
pub struct FeatMap {
    /// How directions weigh in scores. Mirrors `FeatMap.dirScoreMode`.
    pub dir_score_mode: FeatDirScoreMode,
    /// How features are scored. Mirrors `FeatMap.scoreMode`.
    pub score_mode: FeatMapScoreMode,
    /// The scoring parameters, by family; when there are any, every point
    /// must have one. Mirrors `FeatMap.params`.
    pub params: BTreeMap<String, FeatMapParams>,
    feats: Vec<FeatMapPoint>,
}

impl FeatMap {
    /// A map with these parameters and points; points and weights of
    /// different lengths are `Error::Value`, as is a point whose family
    /// the parameters lack. Mirrors `FeatMap(params, feats, weights)`.
    pub fn new_with_params(params: &FeatMapInitParams) -> Result<FeatMap> {
        let mut map = FeatMap {
            params: params.params.clone(),
            ..Default::default()
        };
        if !params.feats.is_empty() {
            if params.feats.len() != params.weights.len() {
                return Err(Error::Value {
                    message: "feats and weights lists must be the same length".to_owned(),
                });
            }
            for (feat, &weight) in params.feats.iter().zip(&params.weights) {
                map.add_feature_with_weight(feat, weight)?;
            }
        }
        Ok(map)
    }

    fn check_family(&self, family: &str) -> Result<()> {
        if !self.params.is_empty() && !self.params.contains_key(family) {
            return Err(Error::Value {
                message: format!("feature family {family} not found in params"),
            });
        }
        Ok(())
    }

    /// Adds a point copied from a feature, with weight 0 (Python's is
    /// `None`, which fails when the point is scored). Mirrors
    /// `FeatMap.AddFeature`.
    pub fn add_feature(&mut self, feat: &dyn FeatureLike) -> Result<()> {
        self.add_feature_with_weight(feat, 0.0)
    }

    /// With the point's weight. Mirrors `FeatMap.AddFeature(feat, weight)`.
    pub fn add_feature_with_weight(&mut self, feat: &dyn FeatureLike, weight: f64) -> Result<()> {
        self.check_family(&feat.feature_family())?;
        let mut new_feat = FeatMapPoint::default();
        new_feat.init_from_feat(feat)?;
        new_feat.weight = weight;
        self.add_feat_point(new_feat)
    }

    /// Adds a point; `Error::Value` if the parameters lack its family.
    /// Mirrors `FeatMap.AddFeatPoint`.
    pub fn add_feat_point(&mut self, feat_pt: FeatMapPoint) -> Result<()> {
        self.check_family(&feat_pt.get_family())?;
        self.feats.push(feat_pt);
        Ok(())
    }

    /// Mirrors `FeatMap.GetFeatures`.
    pub fn get_features(&self) -> &[FeatMapPoint] {
        &self.feats
    }

    /// Mirrors `FeatMap.GetNumFeatures`.
    pub fn get_num_features(&self) -> usize {
        self.feats.len()
    }

    fn index_error(&self, i: usize) -> Error {
        Error::Index {
            index: i as u32,
            len: self.feats.len() as u32,
        }
    }

    /// The point at `i`; `Error::Index` past the end. Mirrors
    /// `FeatMap.GetFeature`.
    pub fn get_feature(&self, i: usize) -> Result<&FeatMapPoint> {
        self.feats.get(i).ok_or_else(|| self.index_error(i))
    }

    /// The point at `i`, to change it in place, as Python's is shared.
    /// Mirrors `FeatMap.GetFeature`.
    pub fn get_feature_mut(&mut self, i: usize) -> Result<&mut FeatMapPoint> {
        let error = self.index_error(i);
        self.feats.get_mut(i).ok_or(error)
    }

    /// Removes the point at `i`; `Error::Index` past the end. Mirrors
    /// `FeatMap.DropFeature`.
    pub fn drop_feature(&mut self, i: usize) -> Result<()> {
        if i >= self.feats.len() {
            return Err(self.index_error(i));
        }
        self.feats.remove(i);
        Ok(())
    }

    /// The score of another feature against one of the map's points: 0 for
    /// another family or beyond the family's radius, else the family's
    /// profile at the distance, times the point's weight and, unless
    /// directions are ignored, their match. A family without parameters is
    /// `Error::Key`. Mirrors `FeatMap.GetFeatFeatScore`.
    pub fn get_feat_feat_score(
        &self,
        feat1: &FeatMapPoint,
        feat2: &dyn FeatureLike,
    ) -> Result<f64> {
        self.get_feat_feat_score_with_type_match(feat1, feat2, true)
    }

    /// Scores features of other families too unless `type_match`. Mirrors
    /// `FeatMap.GetFeatFeatScore(feat1, feat2, typeMatch)`.
    pub fn get_feat_feat_score_with_type_match(
        &self,
        feat1: &FeatMapPoint,
        feat2: &dyn FeatureLike,
        type_match: bool,
    ) -> Result<f64> {
        let family = feat1.get_family();
        if type_match && family != feat2.feature_family() {
            return Ok(0.0);
        }
        let d2 = feat1.get_dist2(feat2)?;
        let params = self.params.get(&family).ok_or(Error::Key { key: family })?;
        if d2 > params.radius * params.radius {
            return Ok(0.0);
        }
        let mut score = match params.feat_profile {
            FeatProfile::Gaussian => (-d2 / params.width).exp(),
            FeatProfile::Triangle => {
                let d = d2.sqrt();
                if d < params.width {
                    1.0 - d / params.width
                } else {
                    0.0
                }
            }
            FeatProfile::Box => 1.0,
        };
        score *= feat1.weight;
        if self.dir_score_mode != FeatDirScoreMode::Ignore {
            let mut dir_score = feat1.get_dir_match(feat2);
            if self.dir_score_mode == FeatDirScoreMode::DotPosRange {
                dir_score = (dir_score + 1.0) / 2.0;
            }
            score *= dir_score;
        }
        Ok(score)
    }

    /// The map's score of the features, as the score mode sets. Mirrors
    /// `FeatMap.ScoreFeats`.
    pub fn score_feats(&self, feats_to_score: &[&dyn FeatureLike]) -> Result<f64> {
        Ok(self.score_feats_with_outputs(feats_to_score)?.0)
    }

    /// The score, with each map point's and each feature's scores and the
    /// points each feature was scored against (A9). Mirrors
    /// `FeatMap.ScoreFeats(featsToScore, mapScoreVect, featsScoreVect,
    /// featsToFeatMapIdx)`.
    pub fn score_feats_with_outputs(
        &self,
        feats_to_score: &[&dyn FeatureLike],
    ) -> Result<(f64, ScoreFeatsOutputs)> {
        let n_to_score = feats_to_score.len();
        let mut map_score_vect = vec![0.0; self.feats.len()];
        let def_score = if self.score_mode == FeatMapScoreMode::Closest {
            1000.0
        } else {
            0.0
        };
        let mut feats_score_vect = vec![def_score; n_to_score];
        // Python's `[-1]` is `Some(None)`: a best or closest point not found yet.
        let mut best: Vec<Option<usize>> = vec![None; n_to_score];
        let mut all: Vec<Vec<usize>> = vec![Vec::new(); n_to_score];
        // Python's Closest pass reuses the loop variable, the last point
        // matched for any feature.
        let mut last_s_feat = None;
        for (o_idx, o_feat) in feats_to_score.iter().enumerate() {
            let family = o_feat.feature_family();
            for (s_idx, s_feat) in self.feats.iter().enumerate() {
                if s_feat.get_family() != family {
                    continue;
                }
                last_s_feat = Some(s_idx);
                match self.score_mode {
                    FeatMapScoreMode::Closest => {
                        let d = s_feat.get_dist2(*o_feat)?;
                        if d < feats_score_vect[o_idx] {
                            feats_score_vect[o_idx] = d;
                            best[o_idx] = Some(s_idx);
                        }
                    }
                    FeatMapScoreMode::Best => {
                        let l_score =
                            self.get_feat_feat_score_with_type_match(s_feat, *o_feat, false)?;
                        if l_score > feats_score_vect[o_idx] {
                            feats_score_vect[o_idx] = l_score;
                            best[o_idx] = Some(s_idx);
                        }
                    }
                    FeatMapScoreMode::All => {
                        let l_score =
                            self.get_feat_feat_score_with_type_match(s_feat, *o_feat, false)?;
                        feats_score_vect[o_idx] += l_score;
                        map_score_vect[s_idx] += l_score;
                        all[o_idx].push(s_idx);
                    }
                }
            }
        }
        let mut tot_score = 0.0;
        if self.score_mode == FeatMapScoreMode::Closest {
            for (o_idx, o_feat) in feats_to_score.iter().enumerate() {
                match (best[o_idx], last_s_feat) {
                    (Some(s_idx), Some(last)) => {
                        let l_score = self.get_feat_feat_score_with_type_match(
                            &self.feats[last],
                            *o_feat,
                            false,
                        )?;
                        feats_score_vect[o_idx] = l_score;
                        map_score_vect[s_idx] = l_score;
                        tot_score += l_score;
                    }
                    _ => feats_score_vect[o_idx] = 0.0,
                }
            }
        } else {
            tot_score = py_sum(feats_score_vect.iter().copied());
            if self.score_mode == FeatMapScoreMode::Best {
                for (o_idx, &l_score) in feats_score_vect.iter().enumerate() {
                    if let Some(s_idx) = best[o_idx] {
                        map_score_vect[s_idx] = l_score;
                    }
                }
            }
        }
        let feats_to_feat_map_idx = if self.score_mode == FeatMapScoreMode::All {
            all
        } else {
            best.into_iter().map(|b| b.into_iter().collect()).collect()
        };
        Ok((
            tot_score,
            ScoreFeatsOutputs {
                map_score_vect,
                feats_score_vect,
                feats_to_feat_map_idx,
            },
        ))
    }
}

/// Python's `'% 6.4f'`: a space for a sign that is not `-`.
fn space_signed(x: f64) -> String {
    let text = if x.is_nan() {
        "nan".to_owned()
    } else if x.is_infinite() {
        "inf".to_owned()
    } else {
        format!("{:.4}", x.abs())
    };
    let sign = if x.is_sign_negative() && !x.is_nan() {
        '-'
    } else {
        ' '
    };
    format!("{:>6}", format!("{sign}{text}"))
}

/// One line per point: its number, family, position, and weight. Mirrors
/// `FeatMap.__str__`.
impl fmt::Display for FeatMap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, feat) in self.feats.iter().enumerate() {
            let pos = feat.get_pos();
            writeln!(
                f,
                "{:>3} {:>12} {} {} {} {}",
                format!(" {}", i + 1),
                feat.get_family(),
                space_signed(pos.x),
                space_signed(pos.y),
                space_signed(pos.z),
                space_signed(feat.weight),
            )?;
        }
        Ok(())
    }
}
