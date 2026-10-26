// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Python built-ins the ports of RDKit's pure-Python descriptor modules
//! reproduce exactly, so their results match Python's, NaN included.

use crate::{Error, Result};

/// Python's `bisect.bisect_right`: the index after the last element of
/// `bins` that is not greater than `x`, found by Python's binary search, so
/// unsorted bins and NaN give Python's answers.
pub(crate) fn bisect_right(bins: &[f64], x: f64) -> usize {
    let (mut lo, mut hi) = (0, bins.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if x < bins[mid] {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Python's `max(values)`: the first value no later one exceeds, as Python
/// compares them (so a leading NaN wins). `Error::Value` for no values, as
/// Python raises.
pub(crate) fn py_max(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut values = values.into_iter();
    let first = values.next().ok_or_else(empty)?;
    Ok(values.fold(first, |best, x| if x > best { x } else { best }))
}

/// Python's `min(values)`, as `py_max`.
pub(crate) fn py_min(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut values = values.into_iter();
    let first = values.next().ok_or_else(empty)?;
    Ok(values.fold(first, |best, x| if x < best { x } else { best }))
}

fn empty() -> Error {
    Error::Value {
        message: "max() or min() of an empty sequence".to_owned(),
    }
}
