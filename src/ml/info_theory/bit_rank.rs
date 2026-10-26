// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.InfoTheory.BitRank`, a Python module: ranking bits
//! by their information gain about activities.
//!
//! In RDKit 2026.09.1 with NumPy 2, Python's functions raise `TypeError`
//! (they make arrays of the abstract `numpy.integer` type); these compute
//! what they mean to, with integer counts, and the gains with
//! `rdInfoTheory.InfoGain` ([`info_gain`]) as Python's do.
//!
//! The bit vectors are given as lists of values, each in `0..n_possible_bit_vals`
//! ([`to_list`](crate::data_structs::ExplicitBitVect::to_list) gives a bit
//! vector's), and the activities as integers below `n_possible_acts`.

use crate::data_structs::SparseBitVect;
use crate::ml::info_theory::rd_info_theory::info_gain;
use crate::{Error, Result};

fn index_error(index: usize, len: usize) -> Error {
    Error::Index {
        index: u32::try_from(index).unwrap_or(u32::MAX),
        len: u32::try_from(len).unwrap_or(u32::MAX),
    }
}

fn check_lengths(n_vects: usize, n_acts: usize) -> Result<()> {
    if n_vects != n_acts {
        return Err(Error::Value {
            message: "var and activity lists should be the same length".to_owned(),
        });
    }
    Ok(())
}

/// The counts of each activity (columns) for each value of bit `which_bit`
/// (rows), for bits of 0 or 1. Mirrors `BitRank.FormCounts(bitVects,
/// actVals, whichBit, nPossibleActs)`.
///
/// `Error::Value` for different numbers of vectors and activities, as
/// Python raises; `Error::Index` for a vector without the bit, or a value
/// or activity past the table.
pub fn form_counts(
    bit_vects: &[Vec<u32>],
    act_vals: &[u32],
    which_bit: usize,
    n_possible_acts: usize,
) -> Result<Vec<Vec<u64>>> {
    form_counts_with_n_possible_bit_vals(bit_vects, act_vals, which_bit, n_possible_acts, 2)
}

/// As [`form_counts`], for bits taking `n_possible_bit_vals` values.
/// Mirrors `BitRank.FormCounts(bitVects, actVals, whichBit, nPossibleActs,
/// nPossibleBitVals)`.
pub fn form_counts_with_n_possible_bit_vals(
    bit_vects: &[Vec<u32>],
    act_vals: &[u32],
    which_bit: usize,
    n_possible_acts: usize,
    n_possible_bit_vals: usize,
) -> Result<Vec<Vec<u64>>> {
    check_lengths(bit_vects.len(), act_vals.len())?;
    let mut res = vec![vec![0u64; n_possible_acts]; n_possible_bit_vals];
    for (v, &act) in bit_vects.iter().zip(act_vals) {
        let bit = *v.get(which_bit).ok_or(index_error(which_bit, v.len()))? as usize;
        let row = res
            .get_mut(bit)
            .ok_or(index_error(bit, n_possible_bit_vals))?;
        let count = row
            .get_mut(act as usize)
            .ok_or(index_error(act as usize, n_possible_acts))?;
        *count += 1;
    }
    Ok(res)
}

/// The information gain of each bit (as many as the first vector has).
/// Mirrors `BitRank.CalcInfoGains(bitVects, actVals, nPossibleActs)`.
///
/// Errors as [`form_counts`], and `Error::Index` for no vectors.
pub fn calc_info_gains(
    bit_vects: &[Vec<u32>],
    act_vals: &[u32],
    n_possible_acts: usize,
) -> Result<Vec<f64>> {
    calc_info_gains_with_n_possible_bit_vals(bit_vects, act_vals, n_possible_acts, 2)
}

/// As [`calc_info_gains`], for bits taking `n_possible_bit_vals` values.
/// Mirrors `BitRank.CalcInfoGains(bitVects, actVals, nPossibleActs,
/// nPossibleBitVals)`.
pub fn calc_info_gains_with_n_possible_bit_vals(
    bit_vects: &[Vec<u32>],
    act_vals: &[u32],
    n_possible_acts: usize,
    n_possible_bit_vals: usize,
) -> Result<Vec<f64>> {
    check_lengths(bit_vects.len(), act_vals.len())?;
    let n_bits = bit_vects.first().ok_or(index_error(0, 0))?.len();
    (0..n_bits)
        .map(|bit| {
            let counts = form_counts_with_n_possible_bit_vals(
                bit_vects,
                act_vals,
                bit,
                n_possible_acts,
                n_possible_bit_vals,
            )?;
            let counts: Vec<Vec<f64>> = counts
                .into_iter()
                .map(|row| row.into_iter().map(|c| c as f64).collect())
                .collect();
            info_gain(&counts)
        })
        .collect()
}

/// A metric for [`rank_bits_with_params`]: given the vectors, the
/// activities, the number of possible activities, and the number of
/// possible bit values, a value for each bit.
pub type MetricFunc = fn(&[Vec<u32>], &[u32], usize, usize) -> Result<Vec<f64>>;

/// Parameters for [`rank_bits_with_params`], with Python's defaults (A2).
#[derive(Clone, Copy, Debug)]
pub struct RankBitsParams {
    /// The number of values a bit takes.
    pub n_possible_bit_vals: usize,
    /// The metric bits are ranked by.
    pub metric_func: MetricFunc,
}

impl Default for RankBitsParams {
    fn default() -> RankBitsParams {
        RankBitsParams {
            n_possible_bit_vals: 2,
            metric_func: calc_info_gains_with_n_possible_bit_vals,
        }
    }
}

/// The indices ranked by their values, largest first. NumPy's `argsort`,
/// reversed: ties come in the reverse of their order, as a stable sort
/// orders them (NumPy's default sort leaves ties' order unspecified past
/// 16 values).
fn rank(metrics: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..metrics.len()).collect();
    order.sort_by(|&a, &b| metrics[a].total_cmp(&metrics[b]));
    order.reverse();
    order
}

/// The bits ranked by information gain, largest first, and each bit's
/// gain. Mirrors `BitRank.RankBits(bitVects, actVals)`.
///
/// Errors as [`calc_info_gains`].
pub fn rank_bits(bit_vects: &[Vec<u32>], act_vals: &[u32]) -> Result<(Vec<usize>, Vec<f64>)> {
    rank_bits_with_params(bit_vects, act_vals, &RankBitsParams::default())
}

/// As [`rank_bits`], by the metric and bit values of `params`. Mirrors
/// `BitRank.RankBits(bitVects, actVals, nPossibleBitVals, metricFunc)`.
pub fn rank_bits_with_params(
    bit_vects: &[Vec<u32>],
    act_vals: &[u32],
    params: &RankBitsParams,
) -> Result<(Vec<usize>, Vec<f64>)> {
    let n_possible_acts = act_vals
        .iter()
        .max()
        .map(|&m| m as usize + 1)
        .ok_or_else(|| Error::Value {
            message: "max() arg is an empty sequence".to_owned(),
        })?;
    let metrics = (params.metric_func)(
        bit_vects,
        act_vals,
        n_possible_acts,
        params.n_possible_bit_vals,
    )?;
    Ok((rank(&metrics), metrics))
}

/// A bit's analysis: the bit, its information gain, and the numbers of
/// actives and inactives that set it.
pub type BitInfo = (u32, f64, u64, u64);

/// For each bit set in any vector, its information gain about binary
/// activities (nonzero is active), with its counts; and the gains alone.
/// Mirrors `BitRank.AnalyzeSparseVects`.
///
/// `Error::Value` for different numbers of vectors and activities, as
/// Python raises; `Error::Index` for no vectors, or a vector longer than
/// the first.
pub fn analyze_sparse_vects(
    bit_vects: &[SparseBitVect],
    act_vals: &[u32],
) -> Result<(Vec<BitInfo>, Vec<f64>)> {
    let n_pts = bit_vects.len();
    check_lengths(n_pts, act_vals.len())?;
    let n_bits = bit_vects.first().ok_or(index_error(0, 0))?.get_num_bits() as usize;
    let mut actives = vec![0u64; n_bits];
    let mut inactives = vec![0u64; n_bits];
    for (sig, &act) in bit_vects.iter().zip(act_vals) {
        let counts = if act != 0 {
            &mut actives
        } else {
            &mut inactives
        };
        for bit in sig.get_on_bits() {
            *counts
                .get_mut(bit as usize)
                .ok_or(index_error(bit as usize, n_bits))? += 1;
        }
    }
    let mut res = Vec::new();
    let mut gains = Vec::new();
    for bit in 0..n_bits {
        let (n_act, n_inact) = (actives[bit], inactives[bit]);
        if n_act != 0 || n_inact != 0 {
            let n = n_pts as u64;
            // Python's table is integers, as these are.
            let table = vec![
                vec![n_act as f64, n_inact as f64],
                vec![(n - n_act) as f64, (n - n_inact) as f64],
            ];
            let gain = info_gain(&table)?;
            gains.push(gain);
            res.push((bit as u32, gain, n_act, n_inact));
        }
    }
    Ok((res, gains))
}

/// A metric for [`sparse_rank_bits_with_metric_func`]: the analysis of
/// each bit, and a value to rank each by.
pub type SparseMetricFunc = fn(&[SparseBitVect], &[u32]) -> Result<(Vec<BitInfo>, Vec<f64>)>;

/// The analyzed bits ranked by information gain, largest first (as
/// indices into the analysis), and the analysis of each, as
/// [`analyze_sparse_vects`] gives it. Mirrors
/// `BitRank.SparseRankBits(bitVects, actVals)`.
///
/// Errors as [`analyze_sparse_vects`].
pub fn sparse_rank_bits(
    bit_vects: &[SparseBitVect],
    act_vals: &[u32],
) -> Result<(Vec<usize>, Vec<BitInfo>)> {
    sparse_rank_bits_with_metric_func(bit_vects, act_vals, analyze_sparse_vects)
}

/// As [`sparse_rank_bits`], by another metric. Mirrors
/// `BitRank.SparseRankBits(bitVects, actVals, metricFunc)`.
pub fn sparse_rank_bits_with_metric_func(
    bit_vects: &[SparseBitVect],
    act_vals: &[u32],
    metric_func: SparseMetricFunc,
) -> Result<(Vec<usize>, Vec<BitInfo>)> {
    let (info, metrics) = metric_func(bit_vects, act_vals)?;
    Ok((rank(&metrics), info))
}
