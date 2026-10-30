// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports the Python code in `rdkit/DataStructs/__init__.py`, re-exported
//! from `rdkit_cxx::data_structs` as Python's is from `rdkit.DataStructs`.

use std::any::TypeId;

use crate::data_structs::{
    BitVect, SparseBitVect, all_bit_similarity, asymmetric_similarity, braun_blanquet_similarity,
    cosine_similarity, dice_similarity, fold_fingerprint_with_fold_factor, kulczynski_similarity,
    mc_connaughey_similarity, rogot_goldberg_similarity, russel_similarity, sokal_similarity,
    tanimoto_similarity,
};
use crate::error::expect_ok;
use crate::ffi::data_structs::ffi;
use crate::{Error, Result};

fn zero_division() -> Error {
    Error::Value {
        message: "division by zero: the fingerprint is empty".to_owned(),
    }
}

/// A similarity metric on bit vectors.
/// A Rust addition, naming the functions of `similarityFunctions`.
pub type SimilarityFunction<B> = fn(&B, &B) -> Result<f64>;

/// The similarity metrics on bit vectors, as `(name, function,
/// description)`; the descriptions are empty. Mirrors
/// `DataStructs.similarityFunctions`.
pub fn similarity_functions<B: BitVect>() -> Vec<(&'static str, SimilarityFunction<B>, &'static str)>
{
    vec![
        ("Tanimoto", tanimoto_similarity::<B>, ""),
        ("Dice", dice_similarity::<B>, ""),
        ("Cosine", cosine_similarity::<B>, ""),
        ("Sokal", sokal_similarity::<B>, ""),
        ("Russel", russel_similarity::<B>, ""),
        ("RogotGoldberg", rogot_goldberg_similarity::<B>, ""),
        ("AllBit", all_bit_similarity::<B>, ""),
        ("Kulczynski", kulczynski_similarity::<B>, ""),
        ("McConnaughey", mc_connaughey_similarity::<B>, ""),
        ("Asymmetric", asymmetric_similarity::<B>, ""),
        ("BraunBlanquet", braun_blanquet_similarity::<B>, ""),
    ]
}

pub(crate) fn num_bits<B: BitVect>(bv: &B) -> u32 {
    expect_ok(ffi::bit_vect_get_num_bits(bv.bit_vect()))
}

pub(crate) fn num_on_bits<B: BitVect>(bv: &B) -> u32 {
    expect_ok(ffi::bit_vect_get_num_on_bits(bv.bit_vect()))
}

pub(crate) fn on_bits<B: BitVect>(bv: &B) -> Vec<u32> {
    expect_ok(ffi::bit_vect_get_on_bits(bv.bit_vect()))
}

pub(crate) fn get_bit<B: BitVect>(bv: &B, which: u32) -> Result<bool> {
    Ok(ffi::bit_vect_get_bit(bv.bit_vect(), which)?)
}

pub(crate) fn set_bit<B: BitVect>(bv: &mut B, which: u32) -> Result<bool> {
    Ok(ffi::bit_vect_set_bit(bv.bit_vect_mut(), which)?)
}

/// A new, empty vector of the same type as `B`, as Python's
/// `bv.__class__(size)`.
pub(crate) fn new_like<B: BitVect + 'static>(size: u32) -> Result<B> {
    if size > i32::MAX as u32 {
        return Err(Error::Value {
            message: "bit vectors hold at most i32::MAX bits".to_owned(),
        });
    }
    let ptr = if TypeId::of::<B>() == TypeId::of::<SparseBitVect>() {
        ffi::bit_vect_new_sparse(size)?
    } else {
        ffi::bit_vect_new_explicit(size, false)?
    };
    Ok(B::from_ffi(ptr))
}

/// The Tanimoto similarity of two fingerprints, folding the longer to the
/// shorter's size first. Mirrors `DataStructs.FingerprintSimilarity(fp1,
/// fp2)`.
///
/// Errors as [`fingerprint_similarity_with_metric`].
pub fn fingerprint_similarity<B: BitVect>(fp1: &B, fp2: &B) -> Result<f64> {
    fingerprint_similarity_with_metric(fp1, fp2, tanimoto_similarity::<B>)
}

/// The similarity of two fingerprints by `metric`, folding the longer by
/// the integer ratio of their sizes first. Mirrors
/// `DataStructs.FingerprintSimilarity(fp1, fp2, metric)`.
///
/// `Error::Value` when one fingerprint is empty and the other is not, as
/// Python raises (`ZeroDivisionError`); errors from the folding and the
/// metric pass through.
pub fn fingerprint_similarity_with_metric<B: BitVect>(
    fp1: &B,
    fp2: &B,
    metric: impl Fn(&B, &B) -> Result<f64>,
) -> Result<f64> {
    let (sz1, sz2) = (num_bits(fp1), num_bits(fp2));
    if sz1 != sz2 && sz1.min(sz2) == 0 {
        return Err(zero_division());
    }
    if sz1 < sz2 {
        let fp2 = fold_fingerprint_with_fold_factor(fp2, sz2 / sz1)?;
        metric(fp1, &fp2)
    } else if sz2 < sz1 {
        let fp1 = fold_fingerprint_with_fold_factor(fp1, sz1 / sz2)?;
        metric(&fp1, fp2)
    } else {
        metric(fp1, fp2)
    }
}

/// Parameters for [`fold_to_target_density_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FoldToTargetDensityParams {
    /// The largest fraction of bits on that stops the folding.
    pub density: f64,
    /// The folding stops before a fold would make the fingerprint this
    /// short or shorter.
    pub min_length: u32,
}

impl Default for FoldToTargetDensityParams {
    fn default() -> FoldToTargetDensityParams {
        FoldToTargetDensityParams {
            density: 0.3,
            min_length: 64,
        }
    }
}

/// The fingerprint folded in half until at most 30% of its bits are on, or
/// another fold would leave it 64 bits or shorter. Mirrors
/// `DataStructs.FoldToTargetDensity(fp)`.
///
/// `Error::Value` for an empty fingerprint, as Python raises
/// (`ZeroDivisionError`); errors from the folding pass through.
pub fn fold_to_target_density<B: BitVect + Clone>(fp: &B) -> Result<B> {
    fold_to_target_density_with_params(fp, &FoldToTargetDensityParams::default())
}

/// Mirrors `DataStructs.FoldToTargetDensity(fp, density, minLength)`.
pub fn fold_to_target_density_with_params<B: BitVect + Clone>(
    fp: &B,
    params: &FoldToTargetDensityParams,
) -> Result<B> {
    if num_bits(fp) == 0 {
        return Err(zero_division());
    }
    let mut fp = fp.clone();
    while num_on_bits(&fp) as f64 / num_bits(&fp) as f64 > params.density
        && num_bits(&fp) / 2 > params.min_length
    {
        fp = fold_fingerprint_with_fold_factor(&fp, 2)?;
    }
    Ok(fp)
}

/// The number of items `n` whose strict lower (or upper) triangular matrix
/// has `matrix.len()` entries. Mirrors `DataStructs.getNForFlatMatrix`.
pub fn get_n_for_flat_matrix<T>(matrix: &[T]) -> usize {
    ((1 + 8 * matrix.len() as u64).isqrt() as usize).div_ceil(2)
}

/// Entry `(i, j)` of a symmetric matrix stored as its strict lower
/// triangle, row by row, with zeros on the diagonal. Mirrors
/// `DataStructs.getElementFromFlatMatrix`.
///
/// `Error::Index` for an entry past the matrix.
pub fn get_element_from_flat_matrix(matrix: &[f64], i: usize, j: usize) -> Result<f64> {
    if i == j {
        return Ok(0.0);
    }
    let (i, j) = if i > j { (j, i) } else { (i, j) };
    let idx = j
        .checked_mul(j - 1)
        .and_then(|x| (x / 2).checked_add(i))
        .unwrap_or(usize::MAX);
    matrix
        .get(idx)
        .copied()
        .ok_or(Error::index(idx, matrix.len()))
}
