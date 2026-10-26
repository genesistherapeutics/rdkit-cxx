// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FeatMaps.FeatMapPoint`, a pure-Python module: a
//! feature map's points.

use std::ops::{Deref, DerefMut};

use crate::Result;
use crate::chem::rd_chemical_features::FreeChemicalFeature;
use crate::chem::rd_mol_chemical_features::MolChemicalFeature;
use crate::geometry::rd_geometry::Point3D;

/// A chemical feature as the feature-map code reads one: a family, a
/// type, a position, and perhaps directions. Python reads these from any
/// feature object; its `featDirs` attribute, which only feature-map points
/// have, is [`feat_dirs`](FeatureLike::feat_dirs).
pub trait FeatureLike {
    /// The feature's family (`GetFamily`).
    fn feature_family(&self) -> String;
    /// The feature's type (`GetType`).
    fn feature_type(&self) -> String;
    /// The feature's position (`GetPos`); a molecule's feature needs a
    /// conformer.
    fn feature_pos(&self) -> Result<Point3D>;
    /// The feature's directions; `None` for features that have none.
    fn feat_dirs(&self) -> Option<&[Point3D]> {
        None
    }
}

impl FeatureLike for FreeChemicalFeature {
    fn feature_family(&self) -> String {
        FreeChemicalFeature::get_family(self)
    }

    fn feature_type(&self) -> String {
        FreeChemicalFeature::get_type(self)
    }

    fn feature_pos(&self) -> Result<Point3D> {
        Ok(FreeChemicalFeature::get_pos(self))
    }
}

impl FeatureLike for MolChemicalFeature<'_> {
    fn feature_family(&self) -> String {
        MolChemicalFeature::get_family(self)
    }

    fn feature_type(&self) -> String {
        MolChemicalFeature::get_type(self)
    }

    fn feature_pos(&self) -> Result<Point3D> {
        MolChemicalFeature::get_pos(self)
    }
}

/// A point of a feature map: a free chemical feature (which it
/// dereferences to, as Python's subclasses one) with a weight and
/// directions. Build the feature with `FreeChemicalFeature`'s constructors
/// and convert it with `From`. Mirrors `FeatMapPoint.FeatMapPoint`.
#[derive(Clone, Debug, Default)]
pub struct FeatMapPoint {
    feature: FreeChemicalFeature,
    /// The point's weight. Mirrors `FeatMapPoint.weight`.
    pub weight: f64,
    /// The point's directions, as unit vectors. Mirrors
    /// `FeatMapPoint.featDirs`.
    pub feat_dirs: Vec<Point3D>,
}

impl From<FreeChemicalFeature> for FeatMapPoint {
    fn from(feature: FreeChemicalFeature) -> FeatMapPoint {
        FeatMapPoint {
            feature,
            weight: 0.0,
            feat_dirs: Vec::new(),
        }
    }
}

impl Deref for FeatMapPoint {
    type Target = FreeChemicalFeature;

    fn deref(&self) -> &FreeChemicalFeature {
        &self.feature
    }
}

impl DerefMut for FeatMapPoint {
    fn deref_mut(&mut self) -> &mut FreeChemicalFeature {
        &mut self.feature
    }
}

impl FeatureLike for FeatMapPoint {
    fn feature_family(&self) -> String {
        self.feature.get_family()
    }

    fn feature_type(&self) -> String {
        self.feature.get_type()
    }

    fn feature_pos(&self) -> Result<Point3D> {
        Ok(self.feature.get_pos())
    }

    fn feat_dirs(&self) -> Option<&[Point3D]> {
        Some(&self.feat_dirs)
    }
}

impl FeatMapPoint {
    /// Takes the feature's family, type, position, and directions, if it
    /// has them. Mirrors `FeatMapPoint.initFromFeat`.
    pub fn init_from_feat(&mut self, feat: &dyn FeatureLike) -> Result<()> {
        self.feature.set_family(&feat.feature_family());
        self.feature.set_type(&feat.feature_type());
        self.feature.set_pos(feat.feature_pos()?);
        if let Some(dirs) = feat.feat_dirs() {
            self.feat_dirs = dirs.to_vec();
        }
        Ok(())
    }

    /// The squared distance to another feature. Mirrors
    /// `FeatMapPoint.GetDist2`.
    pub fn get_dist2(&self, other: &dyn FeatureLike) -> Result<f64> {
        Ok((self.feature.get_pos() - other.feature_pos()?).length_sq())
    }

    /// The best dot product of the points' directions with the other
    /// feature's; 1 when either has none. Mirrors `FeatMapPoint.GetDirMatch`.
    pub fn get_dir_match(&self, other: &dyn FeatureLike) -> f64 {
        self.get_dir_match_with_use_best(other, true)
    }

    /// The mean of the dot products instead unless `use_best`. Mirrors
    /// `FeatMapPoint.GetDirMatch(other, useBest)`.
    pub fn get_dir_match_with_use_best(&self, other: &dyn FeatureLike, use_best: bool) -> f64 {
        let other_dirs = other.feat_dirs().unwrap_or_default();
        if self.feat_dirs.is_empty() || other_dirs.is_empty() {
            return 1.0;
        }
        let mut accum = if use_best { -100000.0 } else { 0.0 };
        for s_dir in &self.feat_dirs {
            for o_dir in other_dirs {
                let d = s_dir.dot_product(o_dir);
                if use_best {
                    if d > accum {
                        accum = d;
                    }
                } else {
                    accum += d;
                }
            }
        }
        if !use_best {
            accum /= (self.feat_dirs.len() * other_dirs.len()) as f64;
        }
        accum
    }
}
