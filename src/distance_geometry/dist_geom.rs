// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.DistanceGeometry.DistGeom`: triangle smoothing and
//! embedding of distance bounds matrices.
//!
//! A bounds matrix is square, given as rows, with lower bounds below the
//! diagonal and upper bounds above it, as `rd_dist_geom` returns them.

use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::distance_geometry::ffi;
use crate::geometry::Point3D;
use crate::{Error, Result};

/// The size of a square, non-empty matrix given as rows.
fn square_size(rows: &[Vec<f64>]) -> Result<u32> {
    if rows.is_empty() {
        return Err(Error::Value {
            message: "The array has to have a nonzero size".to_owned(),
        });
    }
    if rows.iter().any(|row| row.len() != rows.len()) {
        return Err(Error::Value {
            message: "The array has to be square".to_owned(),
        });
    }
    u32::try_from(rows.len()).map_err(|_| Error::Value {
        message: "The array has to be square".to_owned(),
    })
}

/// Triangle-smooths a bounds matrix in place, returning whether the bounds
/// are consistent. Mirrors `DistGeom.DoTriangleSmoothing`. A matrix that is
/// not square, is empty, or has more than 46340 rows is `Error::Value`.
pub fn do_triangle_smoothing(bounds_matrix: &mut [Vec<f64>]) -> Result<bool> {
    do_triangle_smoothing_with_tol(bounds_matrix, 0.0)
}

/// Mirrors `DistGeom.DoTriangleSmoothing(boundsMatrix, tol)`; a lower bound
/// above its upper bound by less than `tol`, relatively, raises the upper
/// bound.
pub fn do_triangle_smoothing_with_tol(bounds_matrix: &mut [Vec<f64>], tol: f64) -> Result<bool> {
    let n = square_size(bounds_matrix)?;
    let mut values = bounds_matrix.concat();
    let res = ffi::do_triangle_smoothing(&mut values, n, tol)?;
    for (row, smoothed) in bounds_matrix.iter_mut().zip(values.chunks(n as usize)) {
        row.copy_from_slice(smoothed);
    }
    Ok(res)
}

/// Parameters for `embed_bounds_matrix_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct EmbedBoundsMatrixParams {
    /// The most random distance matrices to try.
    pub max_iters: i32,
    /// Uses random coordinates when a distance matrix fails to embed.
    pub randomize_on_failure: bool,
    /// Fails with this many zero eigenvalues.
    pub num_zero_fail: i32,
    /// `(i, j, weight)` for bounds whose weight should change.
    pub weights: Vec<(u32, u32, f64)>,
    /// The random seed; −1 leaves RDKit's generator unseeded.
    pub random_seed: i32,
}

impl Default for EmbedBoundsMatrixParams {
    fn default() -> EmbedBoundsMatrixParams {
        EmbedBoundsMatrixParams {
            max_iters: 10,
            randomize_on_failure: false,
            num_zero_fail: 2,
            weights: Vec::new(),
            random_seed: -1,
        }
    }
}

/// Coordinates for each row of a smoothed bounds matrix. Mirrors
/// `DistGeom.EmbedBoundsMatrix`. `Error::Value` for a matrix that cannot be
/// embedded, is not square or finite, has more than 21845 rows (where the
/// minimizer overflows), or has bounds whose squares RDKit would sum past
/// `i32::MAX`; a lower bound above its upper bound is `Error::Runtime`.
pub fn embed_bounds_matrix(bounds_matrix: &[Vec<f64>]) -> Result<Vec<Point3D>> {
    embed_bounds_matrix_with_params(bounds_matrix, &EmbedBoundsMatrixParams::default())
}

/// Mirrors `DistGeom.EmbedBoundsMatrix(boundsMatrix, maxIters,
/// randomizeOnFailure, numZeroFail, weights, randomSeed)`.
pub fn embed_bounds_matrix_with_params(
    bounds_matrix: &[Vec<f64>],
    params: &EmbedBoundsMatrixParams,
) -> Result<Vec<Point3D>> {
    let n = square_size(bounds_matrix)?;
    let pairs: Vec<IndexPair> = params
        .weights
        .iter()
        .map(|&(first, second, _)| IndexPair { first, second })
        .collect();
    let weights: Vec<f64> = params.weights.iter().map(|&(_, _, w)| w).collect();
    Ok(ffi::embed_bounds_matrix(
        &bounds_matrix.concat(),
        n,
        params.max_iters,
        params.randomize_on_failure,
        params.num_zero_fail,
        &pairs,
        &weights,
        params.random_seed,
    )?)
}
