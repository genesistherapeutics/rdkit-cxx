// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Data.Quantize`, a Python module: finding
//! quantization bounds for a variable by the information gain they give
//! about a result.
//!
//! Python's module runs RDKit's C++ `cQuantize` extension when it is built,
//! as conda-forge's is; the search here follows that code, so the bounds
//! and gains match Python's to the last bit. Bounds are "less than": with
//! bounds `[1, 2]`, the values `[0.9, 1, 1.1, 2, 2.2]` quantize to
//! `[0, 1, 1, 2, 2]`.

use crate::{Error, Result};

/// Whether two floats differ by less than `1e-8`. Mirrors `Quantize.feq(v1,
/// v2)`.
pub fn feq(v1: f64, v2: f64) -> bool {
    feq_with_tol(v1, v2, 1e-8)
}

/// Whether two floats differ by less than `tol`. Mirrors `Quantize.feq(v1,
/// v2, tol)`.
pub fn feq_with_tol(v1: f64, v2: f64, tol: f64) -> bool {
    (v1 - v2).abs() < tol
}

/// RDKit's `InfoEntropy` for integer counts (`InfoGainFuncs.h`).
fn info_entropy(counts: &[i64]) -> f64 {
    let n_instances: i64 = counts.iter().sum();
    let mut accum = 0.0;
    if n_instances != 0 {
        for &c in counts {
            let d = c as f64 / n_instances as f64;
            if d != 0.0 {
                accum += -d * d.ln();
            }
        }
    }
    accum / 2f64.ln()
}

/// RDKit's `InfoEntropyGain` for a table of integer counts, `dim1` rows of
/// `dim2` (`InfoGainFuncs.h`).
fn info_entropy_gain(table: &[i64], dim1: usize, dim2: usize) -> f64 {
    let variable_res: Vec<i64> = (0..dim1)
        .map(|i| table[i * dim2..(i + 1) * dim2].iter().sum())
        .collect();
    let overall_res: Vec<i64> = (0..dim2)
        .map(|i| (0..dim1).map(|j| table[j * dim2 + i]).sum())
        .collect();
    let mut term2 = 0.0;
    for i in 0..dim1 {
        term2 += variable_res[i] as f64 * info_entropy(&table[i * dim2..(i + 1) * dim2]);
    }
    // RDKit sums the column totals in `int`.
    let t_sum: i64 = overall_res.iter().map(|&x| i64::from(x as i32)).sum();
    if t_sum != 0 {
        term2 /= t_sum as f64;
        info_entropy(&overall_res) - term2
    } else {
        0.0
    }
}

/// The counts of each result for each bin the cuts make (`GenVarTable` in
/// `cQuantize.cpp`).
fn gen_var_table(
    n_vals: usize,
    cuts: &[usize],
    starts: &[usize],
    results: &[usize],
    n_possible_res: usize,
) -> Vec<i64> {
    let mut var_table = vec![0i64; (cuts.len() + 1) * n_possible_res];
    let mut idx = 0;
    for (i, &cut) in cuts.iter().enumerate() {
        let i_tab = i * n_possible_res;
        while idx < starts[cut] {
            var_table[i_tab + results[idx]] += 1;
            idx += 1;
        }
    }
    let i_tab = cuts.len() * n_possible_res;
    while idx < n_vals {
        var_table[i_tab + results[idx]] += 1;
        idx += 1;
    }
    var_table
}

/// The best gain over the positions of cuts `which` onward, leaving the
/// best cuts in `cuts` (`RecurseHelper` in `cQuantize.cpp`).
fn recurse_helper(
    n_vals: usize,
    cuts: &mut [usize],
    which: usize,
    starts: &[usize],
    results: &[usize],
    n_possible_res: usize,
) -> f64 {
    let n_cuts = cuts.len();
    let mut max_gain = -1e6;
    let mut best_cuts = vec![0usize; n_cuts];
    let highest_cut_here = starts.len() - n_cuts + which;
    let mut var_table = gen_var_table(n_vals, cuts, starts, results, n_possible_res);
    while cuts[which] <= highest_cut_here {
        let gain_here = info_entropy_gain(&var_table, n_cuts + 1, n_possible_res);
        if gain_here > max_gain {
            max_gain = gain_here;
            best_cuts.copy_from_slice(cuts);
        }
        if which < n_cuts - 1 {
            let mut t_cuts = cuts.to_vec();
            let gain_here = recurse_helper(
                n_vals,
                &mut t_cuts,
                which + 1,
                starts,
                results,
                n_possible_res,
            );
            if gain_here > max_gain {
                max_gain = gain_here;
                best_cuts.copy_from_slice(&t_cuts);
            }
        }
        let old_cut = cuts[which];
        cuts[which] += 1;
        let bot = starts[old_cut];
        let top = if old_cut + 1 < starts.len() {
            starts[old_cut + 1]
        } else {
            starts[starts.len() - 1]
        };
        for &v in &results[bot..top] {
            var_table[which * n_possible_res + v] += 1;
            var_table[(which + 1) * n_possible_res + v] -= 1;
        }
        for i in which + 1..n_cuts {
            if cuts[i] == cuts[i - 1] {
                cuts[i] += 1;
            }
        }
    }
    cuts.copy_from_slice(&best_cuts);
    max_gain
}

/// The indices where the sorted values may be cut: where a block of equal
/// values (within 1e-8) starts whose results differ from the last block's
/// (`cQuantize_FindStartPoints`).
fn find_start_points(vals: &[f64], res: &[usize]) -> Vec<usize> {
    let n_data = vals.len();
    let mut start_pts = Vec::new();
    if n_data < 2 {
        return start_pts;
    }
    // RDKit marks a block of mixed results with -1, and the last block's
    // result starts at -2.
    let res = |i: usize| res[i] as i64;
    let mut first_block = true;
    let mut last_block_act = -2;
    let mut block_act = res(0);
    let mut last_div = 0;
    let tol = 1e-8;
    let mut i = 1;
    while i < n_data {
        while i < n_data && vals[i] - vals[i - 1] <= tol {
            if res(i) != block_act {
                block_act = -1;
            }
            i += 1;
        }
        if first_block {
            first_block = false;
            last_block_act = block_act;
            last_div = i;
        } else if block_act == -1 || last_block_act == -1 || block_act != last_block_act {
            start_pts.push(last_div);
            last_div = i;
            last_block_act = block_act;
        } else {
            last_div = i;
        }
        if i < n_data {
            block_act = res(i);
        }
        i += 1;
    }
    if block_act != last_block_act {
        start_pts.push(last_div);
    }
    start_pts
}

/// One quantization bound for a variable and its information gain. Mirrors
/// `Quantize.FindVarQuantBound`.
///
/// Errors as [`find_var_mult_quant_bounds`].
pub fn find_var_quant_bound(
    vals: &[f64],
    results: &[usize],
    n_possible_res: usize,
) -> Result<(f64, f64)> {
    let (bounds, gain) = find_var_mult_quant_bounds(vals, 1, results, n_possible_res)?;
    let bound = *bounds.first().ok_or(Error::index(0, 0))?;
    Ok((bound, gain))
}

/// Up to `n_bounds` quantization bounds for a variable, chosen for the
/// largest information gain about the results (integers below
/// `n_possible_res`), and that gain. With no values, no bounds and a gain
/// of -1e8; with no place to cut, the bound 0 and a gain of 0. Mirrors
/// `Quantize.FindVarMultQuantBounds`.
///
/// `Error::Value` for different numbers of values and results (Python's
/// assertion), and for a NaN value, which Python sorts among the others as
/// its sort's comparisons happen to fall, giving bounds no ordering
/// explains; `Error::Index` for a result not below `n_possible_res`.
pub fn find_var_mult_quant_bounds(
    vals: &[f64],
    n_bounds: usize,
    results: &[usize],
    n_possible_res: usize,
) -> Result<(Vec<f64>, f64)> {
    if vals.len() != results.len() {
        return Err(Error::Value {
            message: "vals/results length mismatch".to_owned(),
        });
    }
    let n_data = vals.len();
    if n_data == 0 {
        return Ok((Vec::new(), -1e8));
    }
    if let Some(&bad) = results.iter().find(|&&r| r >= n_possible_res) {
        return Err(Error::index(bad, n_possible_res));
    }
    if vals.iter().any(|v| v.is_nan()) {
        return Err(Error::Value {
            message: "the values must not be NaN".to_owned(),
        });
    }
    // Python sorts (value, result) pairs.
    let mut svs: Vec<(f64, usize)> = vals.iter().copied().zip(results.iter().copied()).collect();
    svs.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.cmp(&b.1))
    });
    let sort_vals: Vec<f64> = svs.iter().map(|p| p.0).collect();
    let sort_results: Vec<usize> = svs.iter().map(|p| p.1).collect();
    let start_next = find_start_points(&sort_vals, &sort_results);
    if start_next.is_empty() {
        return Ok((vec![0.0], 0.0));
    }
    let mut n_bounds = n_bounds;
    if start_next.len() < n_bounds {
        n_bounds = start_next.len() - 1;
    }
    if n_bounds == 0 {
        n_bounds = 1;
    }
    let mut cuts: Vec<usize> = (0..n_bounds).collect();
    let max_gain = recurse_helper(
        n_data,
        &mut cuts,
        0,
        &start_next,
        &sort_results,
        n_possible_res,
    );
    let n_vs = sort_vals.len();
    let quant_bounds = cuts
        .iter()
        .map(|&cut| {
            let idx = start_next[cut];
            if idx == n_vs {
                sort_vals[n_vs - 1]
            } else if idx == 0 {
                sort_vals[idx]
            } else {
                (sort_vals[idx] + sort_vals[idx - 1]) / 2.0
            }
        })
        .collect();
    Ok((quant_bounds, max_gain))
}

#[cfg(test)]
mod tests {
    use super::find_start_points;

    type Case<'a> = (&'a [(f64, usize)], &'a [usize]);

    /// Python's `test9NewSplits`, against `cQuantize._FindStartPoints`.
    #[test]
    fn start_points() {
        let cases: [Case; 11] = [
            (&[(0., 0), (1., 1), (2., 0)], &[1, 2]),
            (&[(0., 1), (1., 0), (2., 1)], &[1, 2]),
            (
                &[(0., 0), (0., 0), (1., 1), (1., 1), (2., 0), (2., 1)],
                &[2, 4],
            ),
            (
                &[(0., 0), (0., 1), (1., 1), (1., 1), (2., 0), (2., 1)],
                &[2, 4],
            ),
            (
                &[(0., 0), (0., 0), (1., 0), (1., 1), (2., 0), (2., 1)],
                &[2, 4],
            ),
            (
                &[(0., 0), (0., 0), (1., 0), (1., 0), (2., 1), (2., 1)],
                &[4],
            ),
            (
                &[(0., 0), (0., 0), (1., 1), (1., 1), (2., 1), (2., 1)],
                &[2],
            ),
            (&[(0., 0), (0., 0), (1., 0), (1., 0), (2., 0), (2., 0)], &[]),
            (
                &[(0., 0), (0., 1), (1., 0), (1., 1), (2., 0), (2., 0)],
                &[2, 4],
            ),
            (
                &[
                    (1., 0),
                    (2., 1),
                    (2., 1),
                    (3., 1),
                    (3., 1),
                    (3., 1),
                    (4., 0),
                    (4., 1),
                    (4., 1),
                ],
                &[1, 6],
            ),
            (
                &[
                    (1.65175902843, 0),
                    (1.89935600758, 0),
                    (1.89935600758, 1),
                    (1.89935600758, 1),
                    (2.7561609745, 1),
                    (2.7561609745, 1),
                    (2.7561609745, 1),
                    (2.7561609745, 1),
                    (3.53454303741, 1),
                    (3.53454303741, 1),
                    (3.53454303741, 1),
                    (3.53454303741, 1),
                    (3.53454303741, 1),
                ],
                &[1, 4],
            ),
        ];
        for (pairs, expected) in cases {
            let (vals, res): (Vec<f64>, Vec<usize>) = pairs.iter().copied().unzip();
            assert_eq!(find_start_points(&vals, &res), expected, "{pairs:?}");
        }
    }
}
