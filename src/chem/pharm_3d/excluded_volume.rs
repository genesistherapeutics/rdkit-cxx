// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm3D.ExcludedVolume`, a pure-Python module: points
//! a molecule's atoms must stay away from when it is embedded.

use crate::geometry::rd_geometry::Point3D;
use crate::{Error, Result};

/// A point that atoms must stay `exclusion_dist` from, placed by its
/// distances to features. Python's attributes are its fields. Mirrors
/// `ExcludedVolume.ExcludedVolume`.
#[derive(Clone, Debug, PartialEq)]
pub struct ExcludedVolume {
    /// Each feature's atoms, and the volume's least and greatest distance to
    /// them.
    pub feat_info: Vec<(Vec<u32>, f64, f64)>,
    /// The volume's row in a bounds matrix; `None` (Python's −1) until
    /// [`add_excluded_volumes`](crate::chem::pharm_3d::embed_lib::add_excluded_volumes)
    /// adds it.
    pub index: Option<usize>,
    /// The least distance from atoms not in the features.
    pub exclusion_dist: f64,
    /// The volume's position in the last embedding.
    pub pos: Option<Point3D>,
}

/// Parameters for [`ExcludedVolume::new_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExcludedVolumeParams {
    /// The volume's row in a bounds matrix; `None` for Python's −1.
    pub index: Option<usize>,
    /// The least distance from atoms not in the features.
    pub exclusion_dist: f64,
}

impl Default for ExcludedVolumeParams {
    fn default() -> ExcludedVolumeParams {
        ExcludedVolumeParams {
            index: None,
            exclusion_dist: 3.0,
        }
    }
}

impl ExcludedVolume {
    /// A volume placed by its distances to features; no features is
    /// `Error::Value`. Mirrors `ExcludedVolume(featInfo)`.
    pub fn new(feat_info: Vec<(Vec<u32>, f64, f64)>) -> Result<ExcludedVolume> {
        ExcludedVolume::new_with_params(feat_info, &ExcludedVolumeParams::default())
    }

    /// Mirrors `ExcludedVolume(featInfo, index, exclusionDist)`.
    pub fn new_with_params(
        feat_info: Vec<(Vec<u32>, f64, f64)>,
        params: &ExcludedVolumeParams,
    ) -> Result<ExcludedVolume> {
        if feat_info.is_empty() {
            return Err(Error::Value {
                message: "featInfo argument must non-empty".to_owned(),
            });
        }
        Ok(ExcludedVolume {
            feat_info,
            index: params.index,
            exclusion_dist: params.exclusion_dist,
            pos: None,
        })
    }
}
