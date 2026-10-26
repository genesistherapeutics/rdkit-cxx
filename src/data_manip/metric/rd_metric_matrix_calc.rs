// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.DataManip.Metric.rdMetricMatrixCalc`: distance and
//! similarity matrices between items, as the lower triangles of the
//! symmetric matrices. The value for items `i > j` is at
//! `i * (i - 1) / 2 + j`.
//!
//! ```no_run
//! use rdkit_cxx::data_manip::metric::get_euclidean_dist_mat;
//!
//! let dists = get_euclidean_dist_mat(&[vec![0.0, 0.0], vec![3.0, 4.0]]).unwrap();
//! assert_eq!(dists, [5.0]);
//! ```

use crate::data_structs::BitVect;
use crate::error::expect_ok;
use crate::ffi::metric_matrix_calc::ffi;
use crate::{Error, Result};

/// Calls `f` with a C++ list of the vectors, which borrows them and is
/// dropped before this returns.
pub(crate) fn with_bit_vect_list<'a, B: BitVect + 'a, T>(
    bit_vects: impl IntoIterator<Item = &'a B>,
    f: impl FnOnce(&ffi::BitVectList) -> T,
) -> T {
    let mut list = expect_ok(ffi::bit_vect_list_new());
    for bv in bit_vects {
        expect_ok(ffi::bit_vect_list_push(list.pin_mut(), bv.bit_vect()));
    }
    f(&list)
}

fn too_many() -> Error {
    Error::Value {
        message: "at most 65537 items: RDKit indexes the matrix in unsigned int".to_owned(),
    }
}

/// The Euclidean distances between the rows of a descriptor matrix, one
/// row per item. Mirrors `rdMetricMatrixCalc.GetEuclideanDistMat`, which
/// takes a NumPy array or a sequence of sequences and computes in `f64`
/// either way.
///
/// `Error::Runtime` for no rows, as Python raises. `Error::Value` for rows
/// of different lengths, as Python raises, or for more than 65537 rows,
/// where RDKit's `unsigned int` index into the triangle wraps.
pub fn get_euclidean_dist_mat(descrip_mat: &[Vec<f64>]) -> Result<Vec<f64>> {
    let ncols = descrip_mat.first().map_or(0, Vec::len);
    if descrip_mat.iter().any(|row| row.len() != ncols) {
        return Err(Error::Value {
            message: "All subsequences must be the same length".to_owned(),
        });
    }
    let nrows = u32::try_from(descrip_mat.len()).map_err(|_| too_many())?;
    let ncols = u32::try_from(ncols).map_err(|_| Error::Value {
        message: "rows hold at most u32::MAX values: RDKit counts them in unsigned int".to_owned(),
    })?;
    Ok(ffi::get_euclidean_dist_mat(
        &descrip_mat.concat(),
        nrows,
        ncols,
    )?)
}

/// The Tanimoto distances (1 - similarity) between bit vectors, all
/// `ExplicitBitVect`s or all `SparseBitVect`s. Mirrors
/// `rdMetricMatrixCalc.GetTanimotoDistMat`.
///
/// RDKit folds the larger of two vectors of different sizes to the
/// smaller's size, so the sizes must divide evenly. `Error::Runtime` for
/// fewer than two vectors, as Python raises. `Error::Value` for sizes that
/// do not fold, as Python raises; for sizes that differ where one vector
/// has no bits (RDKit divides by zero) or more than `i32::MAX` (RDKit's
/// folding computes in `int`); or for more than 65537 vectors.
pub fn get_tanimoto_dist_mat<'a, B: BitVect + 'a>(
    bit_vect_list: impl IntoIterator<Item = &'a B>,
) -> Result<Vec<f64>> {
    with_bit_vect_list(bit_vect_list, |list| Ok(ffi::get_tanimoto_dist_mat(list)?))
}

/// The Tanimoto similarities between bit vectors, as
/// [`get_tanimoto_dist_mat`]. Mirrors
/// `rdMetricMatrixCalc.GetTanimotoSimMat`.
pub fn get_tanimoto_sim_mat<'a, B: BitVect + 'a>(
    bit_vect_list: impl IntoIterator<Item = &'a B>,
) -> Result<Vec<f64>> {
    with_bit_vect_list(bit_vect_list, |list| Ok(ffi::get_tanimoto_sim_mat(list)?))
}
