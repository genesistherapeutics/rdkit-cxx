// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Data.SplitData`, a Python module: splitting data
//! sets in two at random.
//!
//! The random choices draw from the generator `rdkit_cxx::rd_random`
//! shares, as Python's draw from `RDRandom`, so with the same seed they
//! are Python's. Python's `silent` arguments, which print the sizes of the
//! sets when false, are not taken.

use crate::rd_random;
use crate::{Error, Result};

/// Parameters for [`split_indices_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SplitIndicesParams {
    /// Puts each index in the first set when a random float is below the
    /// fraction, so the sets' sizes vary; the indices stay in order.
    pub legacy: bool,
    /// Draws the first set's indices with replacement, so it may repeat
    /// some; the second set holds those never drawn.
    pub replacement: bool,
}

fn check_frac(frac: f64) -> Result<()> {
    if !(0.0..=1.0).contains(&frac) {
        return Err(Error::Value {
            message: format!("frac must be between 0.0 and 1.0 (frac={frac:.6})"),
        });
    }
    Ok(())
}

/// Splits the indices `0..n_pts` in two, the first set holding `frac` of
/// them (rounded down), by a random permutation. Mirrors
/// `SplitData.SplitIndices(nPts, frac)`.
///
/// `Error::Value` for a fraction outside 0 to 1, as Python raises.
pub fn split_indices(n_pts: usize, frac: f64) -> Result<(Vec<usize>, Vec<usize>)> {
    split_indices_with_params(n_pts, frac, &SplitIndicesParams::default())
}

/// Mirrors `SplitData.SplitIndices(nPts, frac, silent, legacy,
/// replacement)`.
pub fn split_indices_with_params(
    n_pts: usize,
    frac: f64,
    params: &SplitIndicesParams,
) -> Result<(Vec<usize>, Vec<usize>)> {
    check_frac(frac)?;
    let (res_data, res_test);
    if params.replacement {
        let n_train = (n_pts as f64 * frac) as usize;
        let drawn: Vec<usize> = (0..n_train)
            .map(|_| {
                let val = (rd_random::random() * n_pts as f64) as usize;
                if val == n_pts { n_pts - 1 } else { val }
            })
            .collect();
        res_test = (0..n_pts).filter(|i| !drawn.contains(i)).collect();
        res_data = drawn;
    } else if params.legacy {
        let (mut data, mut test) = (Vec::new(), Vec::new());
        for i in 0..n_pts {
            if rd_random::random() < frac {
                data.push(i);
            } else {
                test.push(i);
            }
        }
        res_data = data;
        res_test = test;
    } else {
        let mut perm: Vec<usize> = (0..n_pts).collect();
        rd_random::shuffle_with_random(&mut perm, rd_random::random)?;
        let n_train = (n_pts as f64 * frac) as usize;
        res_test = perm.split_off(n_train);
        res_data = perm;
    }
    Ok((res_data, res_test))
}

/// Splits the examples in two, the first set holding `frac` of them, as
/// [`split_indices`] splits their indices. Mirrors
/// `SplitData.SplitDataSet`.
///
/// `Error::Value` for a fraction outside 0 to 1, as Python raises.
pub fn split_data_set<T: Clone>(data: &[T], frac: f64) -> Result<(Vec<T>, Vec<T>)> {
    check_frac(frac)?;
    let (train, test) = split_indices(data.len(), frac)?;
    Ok((
        train.into_iter().map(|x| data[x].clone()).collect(),
        test.into_iter().map(|x| data[x].clone()).collect(),
    ))
}
