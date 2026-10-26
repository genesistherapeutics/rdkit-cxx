// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Python's formatting of floats, for ports that print what Python prints.

/// Python's `%f`: six decimals, with Python's `nan` and `inf`.
pub(crate) fn format_f(x: f64) -> String {
    if x.is_nan() {
        "nan".to_owned()
    } else {
        format!("{x:.6}")
    }
}
