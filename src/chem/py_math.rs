// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Python's and NumPy's floating-point reductions, for ports that must
//! compute what Python computes to the last bit (CONVENTIONS.md R2).

/// Python's `sum` of floats (3.12 and later): Neumaier's compensated sum,
/// left to right from 0.
pub(crate) fn py_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut c) = (0.0f64, 0.0f64);
    for x in values {
        let t = sum + x;
        if sum.abs() >= x.abs() {
            c += (sum - t) + x;
        } else {
            c += (x - t) + sum;
        }
        sum = t;
    }
    if c != 0.0 && c.is_finite() {
        sum += c;
    }
    sum
}

/// NumPy's pairwise sum of a contiguous `float64` array (`np.add.reduce`).
pub(crate) fn np_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - n % 8 {
            for (j, rj) in r.iter_mut().enumerate() {
                *rj += a[i + j];
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
        np_sum(&a[..n2]) + np_sum(&a[n2..])
    }
}

/// NumPy's `np.std` of a list of floats: the population standard deviation.
pub(crate) fn np_std(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    let mean = (0.0 + np_sum(a)) / n;
    let squares: Vec<f64> = a.iter().map(|&x| (x - mean) * (x - mean)).collect();
    ((0.0 + np_sum(&squares)) / n).sqrt()
}
