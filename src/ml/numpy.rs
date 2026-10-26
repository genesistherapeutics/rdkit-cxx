// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! NumPy's floating-point reductions, in NumPy's order, so ports of Python
//! code that sums with NumPy get Python's results to the last bit.

/// NumPy's pairwise summation (`pairwise_sum` in NumPy's
/// `loops_utils.h.src`).
fn pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        a.iter().fold(-0.0, |res, &x| res + x)
    } else if n <= 128 {
        let mut r = [0.0; 8];
        r.copy_from_slice(&a[..8]);
        let mut i = 8;
        while i < n - n % 8 {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for &x in &a[i..] {
            res += x;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise_sum(&a[..n2]) + pairwise_sum(&a[n2..])
    }
}

/// `numpy.sum` of a contiguous float array: its identity, 0.0, plus the
/// pairwise sum.
pub(crate) fn sum(a: &[f64]) -> f64 {
    0.0 + pairwise_sum(a)
}

/// Python's `%f`: six decimals, with Python's `nan` and `inf`.
pub(crate) fn format_f(x: f64) -> String {
    if x.is_nan() {
        "nan".to_owned()
    } else {
        format!("{x:.6}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sum_follows_numpy_blocks() {
        // Values from numpy.sum in NumPy 2.5.
        let a: Vec<f64> = (0..300).map(|i| 1.0 / (i as f64 + 1.0)).collect();
        assert_eq!(sum(&a[..5]), 2.283333333333333);
        assert_eq!(sum(&a), 6.282663880299504);
        assert_eq!(sum(&[-0.0]), 0.0);
        assert!(sum(&[-0.0]).is_sign_positive());
    }
}
