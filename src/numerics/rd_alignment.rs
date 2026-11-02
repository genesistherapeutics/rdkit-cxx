// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Numerics.rdAlignment`: the optimal superposition of two
//! sets of points.

use crate::Result;
use crate::ffi::alignment::ffi;
use crate::geometry::Point3D;

/// Parameters for `get_alignment_transform_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetAlignmentTransformParams {
    /// A positive weight for each pair of points; `None` weighs them
    /// equally. Python treats an empty list as `None`.
    pub weights: Option<Vec<f64>>,
    /// Reflects the probe points first.
    pub reflect: bool,
    /// The most iterations of the eigensolver.
    pub max_iterations: u32,
}

impl Default for GetAlignmentTransformParams {
    fn default() -> GetAlignmentTransformParams {
        GetAlignmentTransformParams {
            weights: None,
            reflect: false,
            max_iterations: 50,
        }
    }
}

/// The sum of squared deviations after the probe points' optimal (least
/// RMSD) superposition on the reference points, and the 4 x 4 transform
/// that does it, as rows. Mirrors `rdAlignment.GetAlignmentTransform`.
/// Empty or mismatched lists are `Error::Value`.
pub fn get_alignment_transform(
    ref_points: &[Point3D],
    probe_points: &[Point3D],
) -> Result<(f64, Vec<Vec<f64>>)> {
    get_alignment_transform_with_params(
        ref_points,
        probe_points,
        &GetAlignmentTransformParams::default(),
    )
}

/// Mirrors `rdAlignment.GetAlignmentTransform(refPoints, probePoints,
/// weights, reflect, maxIterations)`. A weight count other than the number
/// of points is `Error::Value`, and a weight that is not positive
/// `Error::Runtime`.
pub fn get_alignment_transform_with_params(
    ref_points: &[Point3D],
    probe_points: &[Point3D],
    params: &GetAlignmentTransformParams,
) -> Result<(f64, Vec<Vec<f64>>)> {
    let mut transform = Vec::new();
    let ssd = ffi::align_point_pairs(
        ref_points,
        probe_points,
        params.weights.is_some(),
        params.weights.as_deref().unwrap_or_default(),
        params.reflect,
        params.max_iterations,
        &mut transform,
    )?;
    Ok((ssd, transform.chunks(4).map(<[f64]>::to_vec).collect()))
}
