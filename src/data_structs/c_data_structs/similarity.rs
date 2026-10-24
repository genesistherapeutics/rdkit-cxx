// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Similarity metrics, as `rdkit.DataStructs.cDataStructs` defines them.
//!
//! Each bit-vector metric `xxx_similarity` takes two `ExplicitBitVect`s or two
//! `SparseBitVect`s, with `_with_return_distance` for `1 - similarity`,
//! `_bytes` to compare against a pickle, `bulk_` for one vector against many,
//! and `_neighbors`/`_neighbors_sparse` for each query's nearest vector. The
//! Bulk and Neighbors loops, which Python's wrappers write in C++, are written
//! here in Rust.

use crate::Result;
use crate::data_structs::c_data_structs::{
    BitVect, ExplicitBitVect, SparseBitVect, SparseIntVectType,
};
use crate::ffi::data_structs::ffi;

// The metric numbers of bit_vect_similarity in
// wrapper/src/data_structs.cpp.
const TANIMOTO: u8 = 0;
const COSINE: u8 = 1;
const KULCZYNSKI: u8 = 2;
const DICE: u8 = 3;
const SOKAL: u8 = 4;
const MC_CONNAUGHEY: u8 = 5;
const ASYMMETRIC: u8 = 6;
const BRAUN_BLANQUET: u8 = 7;
const RUSSEL: u8 = 8;
const ROGOT_GOLDBERG: u8 = 9;
const ON_BIT: u8 = 10;
const ALL_BIT: u8 = 11;
const TVERSKY: u8 = 12;

fn similarity<B: BitVect>(bv1: &B, bv2: &B, metric: u8, a: f64, b: f64, rd: bool) -> Result<f64> {
    Ok(ffi::bit_vect_similarity(
        bv1.bit_vect(),
        bv2.bit_vect(),
        metric,
        a,
        b,
        rd,
    )?)
}

fn similarity_bytes<B: BitVect>(
    bv1: &B,
    pkl: &[u8],
    metric: u8,
    a: f64,
    b: f64,
    rd: bool,
) -> Result<f64> {
    Ok(ffi::bit_vect_similarity_pkl(
        bv1.bit_vect(),
        pkl,
        metric,
        a,
        b,
        rd,
    )?)
}

fn bulk<'a, B: BitVect + 'a>(
    bv1: &B,
    bv_list: impl IntoIterator<Item = &'a B>,
    metric: u8,
    a: f64,
    b: f64,
    rd: bool,
) -> Result<Vec<f64>> {
    bv_list
        .into_iter()
        .map(|bv2| similarity(bv1, bv2, metric, a, b, rd))
        .collect()
}

/// As NeighborWrapper in Code/DataStructs/nbWrap/wrap_BitOps.cpp: for each
/// query, the index and similarity of the first most similar vector.
fn neighbors<'a, 'b, B: BitVect + 'a + 'b>(
    bvqueries: impl IntoIterator<Item = &'a B>,
    bv_list: impl IntoIterator<Item = &'b B>,
    metric: u8,
) -> Result<Vec<(u32, f64)>> {
    let bv_list: Vec<&B> = bv_list.into_iter().collect();
    bvqueries
        .into_iter()
        .map(|query| {
            let mut closest = -1.0;
            let mut nbr = 0;
            for (j, bv) in bv_list.iter().enumerate() {
                let sim = similarity(query, *bv, metric, 0.0, 0.0, false)?;
                if sim > closest {
                    closest = sim;
                    nbr = u32::try_from(j).expect("fewer than 2^32 vectors");
                }
            }
            Ok((nbr, closest))
        })
        .collect()
}

/// Defines a metric's functions (BIG_DEF in wrap_BitOps.cpp).
macro_rules! big_metric {
    (
        $metric:expr, $label:literal,
        $sim:ident, $sim_rd:ident, $sim_bytes:ident, $sim_bytes_rd:ident,
        $bulk:ident, $bulk_rd:ident, $neighbors:ident, $neighbors_sparse:ident
    ) => {
        #[doc = concat!("The ", $label, " similarity of two bit vectors of one size. Mirrors `cDataStructs.", $label, "Similarity`.")]
        pub fn $sim<B: BitVect>(bv1: &B, bv2: &B) -> Result<f64> {
            similarity(bv1, bv2, $metric, 0.0, 0.0, false)
        }

        #[doc = concat!("As [`", stringify!($sim), "`], or 1 minus it when `return_distance`.")]
        pub fn $sim_rd<B: BitVect>(bv1: &B, bv2: &B, return_distance: bool) -> Result<f64> {
            similarity(bv1, bv2, $metric, 0.0, 0.0, return_distance)
        }

        #[doc = concat!("As [`", stringify!($sim), "`], against the pickle (`to_binary`) of a vector of the same type (A5).")]
        pub fn $sim_bytes<B: BitVect>(bv1: &B, pkl: &[u8]) -> Result<f64> {
            similarity_bytes(bv1, pkl, $metric, 0.0, 0.0, false)
        }

        #[doc = concat!("As [`", stringify!($sim_bytes), "`], or 1 minus it when `return_distance`.")]
        pub fn $sim_bytes_rd<B: BitVect>(bv1: &B, pkl: &[u8], return_distance: bool) -> Result<f64> {
            similarity_bytes(bv1, pkl, $metric, 0.0, 0.0, return_distance)
        }

        #[doc = concat!("[`", stringify!($sim), "`] against each vector. Mirrors `cDataStructs.Bulk", $label, "Similarity`.")]
        pub fn $bulk<'a, B: BitVect + 'a>(
            bv1: &B,
            bv_list: impl IntoIterator<Item = &'a B>,
        ) -> Result<Vec<f64>> {
            bulk(bv1, bv_list, $metric, 0.0, 0.0, false)
        }

        #[doc = concat!("As [`", stringify!($bulk), "`], or distances when `return_distance`.")]
        pub fn $bulk_rd<'a, B: BitVect + 'a>(
            bv1: &B,
            bv_list: impl IntoIterator<Item = &'a B>,
            return_distance: bool,
        ) -> Result<Vec<f64>> {
            bulk(bv1, bv_list, $metric, 0.0, 0.0, return_distance)
        }

        #[doc = concat!("For each query, the index of its most similar vector and their similarity. Mirrors `cDataStructs.", $label, "SimilarityNeighbors`.")]
        pub fn $neighbors<'a, 'b>(
            bvqueries: impl IntoIterator<Item = &'a ExplicitBitVect>,
            bv_list: impl IntoIterator<Item = &'b ExplicitBitVect>,
        ) -> Result<Vec<(u32, f64)>> {
            neighbors(bvqueries, bv_list, $metric)
        }

        #[doc = concat!("As [`", stringify!($neighbors), "`], for `SparseBitVect`s. Mirrors `cDataStructs.", $label, "SimilarityNeighbors_sparse`.")]
        pub fn $neighbors_sparse<'a, 'b>(
            bvqueries: impl IntoIterator<Item = &'a SparseBitVect>,
            bv_list: impl IntoIterator<Item = &'b SparseBitVect>,
        ) -> Result<Vec<(u32, f64)>> {
            neighbors(bvqueries, bv_list, $metric)
        }
    };
}

big_metric!(
    TANIMOTO,
    "Tanimoto",
    tanimoto_similarity,
    tanimoto_similarity_with_return_distance,
    tanimoto_similarity_bytes,
    tanimoto_similarity_bytes_with_return_distance,
    bulk_tanimoto_similarity,
    bulk_tanimoto_similarity_with_return_distance,
    tanimoto_similarity_neighbors,
    tanimoto_similarity_neighbors_sparse
);
big_metric!(
    COSINE,
    "Cosine",
    cosine_similarity,
    cosine_similarity_with_return_distance,
    cosine_similarity_bytes,
    cosine_similarity_bytes_with_return_distance,
    bulk_cosine_similarity,
    bulk_cosine_similarity_with_return_distance,
    cosine_similarity_neighbors,
    cosine_similarity_neighbors_sparse
);
big_metric!(
    KULCZYNSKI,
    "Kulczynski",
    kulczynski_similarity,
    kulczynski_similarity_with_return_distance,
    kulczynski_similarity_bytes,
    kulczynski_similarity_bytes_with_return_distance,
    bulk_kulczynski_similarity,
    bulk_kulczynski_similarity_with_return_distance,
    kulczynski_similarity_neighbors,
    kulczynski_similarity_neighbors_sparse
);
big_metric!(
    DICE,
    "Dice",
    dice_similarity,
    dice_similarity_with_return_distance,
    dice_similarity_bytes,
    dice_similarity_bytes_with_return_distance,
    bulk_dice_similarity,
    bulk_dice_similarity_with_return_distance,
    dice_similarity_neighbors,
    dice_similarity_neighbors_sparse
);
big_metric!(
    SOKAL,
    "Sokal",
    sokal_similarity,
    sokal_similarity_with_return_distance,
    sokal_similarity_bytes,
    sokal_similarity_bytes_with_return_distance,
    bulk_sokal_similarity,
    bulk_sokal_similarity_with_return_distance,
    sokal_similarity_neighbors,
    sokal_similarity_neighbors_sparse
);
big_metric!(
    MC_CONNAUGHEY,
    "McConnaughey",
    mc_connaughey_similarity,
    mc_connaughey_similarity_with_return_distance,
    mc_connaughey_similarity_bytes,
    mc_connaughey_similarity_bytes_with_return_distance,
    bulk_mc_connaughey_similarity,
    bulk_mc_connaughey_similarity_with_return_distance,
    mc_connaughey_similarity_neighbors,
    mc_connaughey_similarity_neighbors_sparse
);
big_metric!(
    ASYMMETRIC,
    "Asymmetric",
    asymmetric_similarity,
    asymmetric_similarity_with_return_distance,
    asymmetric_similarity_bytes,
    asymmetric_similarity_bytes_with_return_distance,
    bulk_asymmetric_similarity,
    bulk_asymmetric_similarity_with_return_distance,
    asymmetric_similarity_neighbors,
    asymmetric_similarity_neighbors_sparse
);
big_metric!(
    BRAUN_BLANQUET,
    "BraunBlanquet",
    braun_blanquet_similarity,
    braun_blanquet_similarity_with_return_distance,
    braun_blanquet_similarity_bytes,
    braun_blanquet_similarity_bytes_with_return_distance,
    bulk_braun_blanquet_similarity,
    bulk_braun_blanquet_similarity_with_return_distance,
    braun_blanquet_similarity_neighbors,
    braun_blanquet_similarity_neighbors_sparse
);
big_metric!(
    RUSSEL,
    "Russel",
    russel_similarity,
    russel_similarity_with_return_distance,
    russel_similarity_bytes,
    russel_similarity_bytes_with_return_distance,
    bulk_russel_similarity,
    bulk_russel_similarity_with_return_distance,
    russel_similarity_neighbors,
    russel_similarity_neighbors_sparse
);
big_metric!(
    ROGOT_GOLDBERG,
    "RogotGoldberg",
    rogot_goldberg_similarity,
    rogot_goldberg_similarity_with_return_distance,
    rogot_goldberg_similarity_bytes,
    rogot_goldberg_similarity_bytes_with_return_distance,
    bulk_rogot_goldberg_similarity,
    bulk_rogot_goldberg_similarity_with_return_distance,
    rogot_goldberg_similarity_neighbors,
    rogot_goldberg_similarity_neighbors_sparse
);

/// Defines a metric with no return distance or pickle forms (DBL_DEF).
macro_rules! dbl_metric {
    ($metric:expr, $label:literal, $sim:ident, $bulk:ident, $bulk_rd:ident) => {
        #[doc = concat!("The ", $label, " similarity of two bit vectors of one size. Mirrors `cDataStructs.", $label, "Similarity`.")]
        pub fn $sim<B: BitVect>(v1: &B, v2: &B) -> Result<f64> {
            similarity(v1, v2, $metric, 0.0, 0.0, false)
        }

        #[doc = concat!("[`", stringify!($sim), "`] against each vector. Mirrors `cDataStructs.Bulk", $label, "Similarity`.")]
        pub fn $bulk<'a, B: BitVect + 'a>(
            v1: &B,
            v2: impl IntoIterator<Item = &'a B>,
        ) -> Result<Vec<f64>> {
            bulk(v1, v2, $metric, 0.0, 0.0, false)
        }

        #[doc = concat!("As [`", stringify!($bulk), "`], or distances when `return_distance`.")]
        pub fn $bulk_rd<'a, B: BitVect + 'a>(
            v1: &B,
            v2: impl IntoIterator<Item = &'a B>,
            return_distance: bool,
        ) -> Result<Vec<f64>> {
            bulk(v1, v2, $metric, 0.0, 0.0, return_distance)
        }
    };
}

dbl_metric!(
    ON_BIT,
    "OnBit",
    on_bit_similarity,
    bulk_on_bit_similarity,
    bulk_on_bit_similarity_with_return_distance
);
dbl_metric!(
    ALL_BIT,
    "AllBit",
    all_bit_similarity,
    bulk_all_bit_similarity,
    bulk_all_bit_similarity_with_return_distance
);

/// The Tversky similarity of two bit vectors of one size, with weights `a`
/// and `b`. Mirrors `cDataStructs.TverskySimilarity`.
pub fn tversky_similarity<B: BitVect>(bv1: &B, bv2: &B, a: f64, b: f64) -> Result<f64> {
    similarity(bv1, bv2, TVERSKY, a, b, false)
}

/// As [`tversky_similarity`], or 1 minus it when `return_distance`.
pub fn tversky_similarity_with_return_distance<B: BitVect>(
    bv1: &B,
    bv2: &B,
    a: f64,
    b: f64,
    return_distance: bool,
) -> Result<f64> {
    similarity(bv1, bv2, TVERSKY, a, b, return_distance)
}

/// As [`tversky_similarity`], against a pickle (A5).
pub fn tversky_similarity_bytes<B: BitVect>(bv1: &B, pkl: &[u8], a: f64, b: f64) -> Result<f64> {
    similarity_bytes(bv1, pkl, TVERSKY, a, b, false)
}

/// As [`tversky_similarity_bytes`], or 1 minus it when `return_distance`.
pub fn tversky_similarity_bytes_with_return_distance<B: BitVect>(
    bv1: &B,
    pkl: &[u8],
    a: f64,
    b: f64,
    return_distance: bool,
) -> Result<f64> {
    similarity_bytes(bv1, pkl, TVERSKY, a, b, return_distance)
}

/// [`tversky_similarity`] against each vector. Mirrors
/// `cDataStructs.BulkTverskySimilarity`.
pub fn bulk_tversky_similarity<'a, B: BitVect + 'a>(
    bv1: &B,
    bv_list: impl IntoIterator<Item = &'a B>,
    a: f64,
    b: f64,
) -> Result<Vec<f64>> {
    bulk(bv1, bv_list, TVERSKY, a, b, false)
}

/// As [`bulk_tversky_similarity`], or distances when `return_distance`.
pub fn bulk_tversky_similarity_with_return_distance<'a, B: BitVect + 'a>(
    bv1: &B,
    bv_list: impl IntoIterator<Item = &'a B>,
    a: f64,
    b: f64,
    return_distance: bool,
) -> Result<Vec<f64>> {
    bulk(bv1, bv_list, TVERSKY, a, b, return_distance)
}

/// Parameters for the sparse integer vector similarities'
/// `_with_params` forms, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SparseIntVectSimilarityParams {
    /// Returns 1 minus the similarity.
    pub return_distance: bool,
    /// Returns 0 early when the similarity cannot reach this.
    pub bounds: f64,
}

fn siv_similarity<S: SparseIntVectType>(
    siv1: &S,
    siv2: &S,
    metric: u8,
    a: f64,
    b: f64,
    params: &SparseIntVectSimilarityParams,
) -> Result<f64> {
    Ok(ffi::sparse_int_vect_similarity(
        siv1.sparse_int_vect(),
        siv2.sparse_int_vect(),
        metric,
        a,
        b,
        params.return_distance,
        params.bounds,
    )?)
}

fn siv_bulk<'a, S: SparseIntVectType + 'a>(
    v1: &S,
    v2: impl IntoIterator<Item = &'a S>,
    metric: u8,
    a: f64,
    b: f64,
    return_distance: bool,
) -> Result<Vec<f64>> {
    let params = SparseIntVectSimilarityParams {
        return_distance,
        ..Default::default()
    };
    v2.into_iter()
        .map(|siv2| siv_similarity(v1, siv2, metric, a, b, &params))
        .collect()
}

/// The Dice similarity of two sparse integer vectors of one length. Mirrors
/// `cDataStructs.DiceSimilarity(siv1, siv2)` (A5).
///
/// RDKit sums the values' absolute values in `i32`, so vectors whose
/// absolute values sum past `i32::MAX` are `Error::Value`.
pub fn dice_similarity_sparse_int_vect<S: SparseIntVectType>(siv1: &S, siv2: &S) -> Result<f64> {
    dice_similarity_sparse_int_vect_with_params(siv1, siv2, &Default::default())
}

/// Mirrors `cDataStructs.DiceSimilarity(siv1, siv2, returnDistance, bounds)`.
pub fn dice_similarity_sparse_int_vect_with_params<S: SparseIntVectType>(
    siv1: &S,
    siv2: &S,
    params: &SparseIntVectSimilarityParams,
) -> Result<f64> {
    siv_similarity(siv1, siv2, 0, 0.0, 0.0, params)
}

/// The Tanimoto similarity of two sparse integer vectors of one length.
/// Mirrors `cDataStructs.TanimotoSimilarity(siv1, siv2)` (A5).
///
/// RDKit sums the values' absolute values in `i32`, so vectors whose
/// absolute values sum past `i32::MAX` are `Error::Value`.
pub fn tanimoto_similarity_sparse_int_vect<S: SparseIntVectType>(
    siv1: &S,
    siv2: &S,
) -> Result<f64> {
    tanimoto_similarity_sparse_int_vect_with_params(siv1, siv2, &Default::default())
}

/// Mirrors `cDataStructs.TanimotoSimilarity(siv1, siv2, returnDistance,
/// bounds)`.
pub fn tanimoto_similarity_sparse_int_vect_with_params<S: SparseIntVectType>(
    siv1: &S,
    siv2: &S,
    params: &SparseIntVectSimilarityParams,
) -> Result<f64> {
    siv_similarity(siv1, siv2, 1, 0.0, 0.0, params)
}

/// The Tversky similarity of two sparse integer vectors of one length.
/// Mirrors `cDataStructs.TverskySimilarity(siv1, siv2, a, b)` (A5).
///
/// The values' absolute values must sum to at most `i32::MAX`, as for
/// `dice_similarity_sparse_int_vect`.
pub fn tversky_similarity_sparse_int_vect<S: SparseIntVectType>(
    siv1: &S,
    siv2: &S,
    a: f64,
    b: f64,
) -> Result<f64> {
    tversky_similarity_sparse_int_vect_with_params(siv1, siv2, a, b, &Default::default())
}

/// Mirrors `cDataStructs.TverskySimilarity(siv1, siv2, a, b, returnDistance,
/// bounds)`.
pub fn tversky_similarity_sparse_int_vect_with_params<S: SparseIntVectType>(
    siv1: &S,
    siv2: &S,
    a: f64,
    b: f64,
    params: &SparseIntVectSimilarityParams,
) -> Result<f64> {
    siv_similarity(siv1, siv2, 2, a, b, params)
}

/// [`dice_similarity_sparse_int_vect`] against each vector. Mirrors
/// `cDataStructs.BulkDiceSimilarity(v1, v2)` (A5).
pub fn bulk_dice_similarity_sparse_int_vect<'a, S: SparseIntVectType + 'a>(
    v1: &S,
    v2: impl IntoIterator<Item = &'a S>,
) -> Result<Vec<f64>> {
    siv_bulk(v1, v2, 0, 0.0, 0.0, false)
}

/// As [`bulk_dice_similarity_sparse_int_vect`], or distances when
/// `return_distance`.
pub fn bulk_dice_similarity_sparse_int_vect_with_return_distance<'a, S: SparseIntVectType + 'a>(
    v1: &S,
    v2: impl IntoIterator<Item = &'a S>,
    return_distance: bool,
) -> Result<Vec<f64>> {
    siv_bulk(v1, v2, 0, 0.0, 0.0, return_distance)
}

/// [`tanimoto_similarity_sparse_int_vect`] against each vector. Mirrors
/// `cDataStructs.BulkTanimotoSimilarity(v1, v2)` (A5).
pub fn bulk_tanimoto_similarity_sparse_int_vect<'a, S: SparseIntVectType + 'a>(
    v1: &S,
    v2: impl IntoIterator<Item = &'a S>,
) -> Result<Vec<f64>> {
    siv_bulk(v1, v2, 1, 0.0, 0.0, false)
}

/// As [`bulk_tanimoto_similarity_sparse_int_vect`], or distances when
/// `return_distance`.
pub fn bulk_tanimoto_similarity_sparse_int_vect_with_return_distance<
    'a,
    S: SparseIntVectType + 'a,
>(
    v1: &S,
    v2: impl IntoIterator<Item = &'a S>,
    return_distance: bool,
) -> Result<Vec<f64>> {
    siv_bulk(v1, v2, 1, 0.0, 0.0, return_distance)
}

/// [`tversky_similarity_sparse_int_vect`] against each vector. Mirrors
/// `cDataStructs.BulkTverskySimilarity(v1, v2, a, b)` (A5).
pub fn bulk_tversky_similarity_sparse_int_vect<'a, S: SparseIntVectType + 'a>(
    v1: &S,
    v2: impl IntoIterator<Item = &'a S>,
    a: f64,
    b: f64,
) -> Result<Vec<f64>> {
    siv_bulk(v1, v2, 2, a, b, false)
}

/// As [`bulk_tversky_similarity_sparse_int_vect`], or distances when
/// `return_distance`.
pub fn bulk_tversky_similarity_sparse_int_vect_with_return_distance<
    'a,
    S: SparseIntVectType + 'a,
>(
    v1: &S,
    v2: impl IntoIterator<Item = &'a S>,
    a: f64,
    b: f64,
    return_distance: bool,
) -> Result<Vec<f64>> {
    siv_bulk(v1, v2, 2, a, b, return_distance)
}

/// The fraction of each vector's on bits the other shares. Mirrors
/// `cDataStructs.OnBitProjSimilarity`.
pub fn on_bit_proj_similarity<B: BitVect>(bv1: &B, bv2: &B) -> Result<Vec<f64>> {
    Ok(ffi::bit_vect_proj_similarity(
        bv1.bit_vect(),
        bv2.bit_vect(),
        true,
    )?)
}

/// The fraction of each vector's off bits the other shares. Mirrors
/// `cDataStructs.OffBitProjSimilarity`.
pub fn off_bit_proj_similarity<B: BitVect>(bv1: &B, bv2: &B) -> Result<Vec<f64>> {
    Ok(ffi::bit_vect_proj_similarity(
        bv1.bit_vect(),
        bv2.bit_vect(),
        false,
    )?)
}

/// The number of bits with the same value in both. Mirrors
/// `cDataStructs.NumBitsInCommon`.
pub fn num_bits_in_common<B: BitVect>(bv1: &B, bv2: &B) -> Result<i32> {
    Ok(ffi::bit_vect_num_bits_in_common(
        bv1.bit_vect(),
        bv2.bit_vect(),
    )?)
}

/// The bits on in both. Mirrors `cDataStructs.OnBitsInCommon`.
pub fn on_bits_in_common<B: BitVect>(bv1: &B, bv2: &B) -> Result<Vec<u32>> {
    Ok(ffi::bit_vect_bits_in_common(
        bv1.bit_vect(),
        bv2.bit_vect(),
        true,
    )?)
}

/// The bits off in both. Mirrors `cDataStructs.OffBitsInCommon`.
pub fn off_bits_in_common<B: BitVect>(bv1: &B, bv2: &B) -> Result<Vec<u32>> {
    Ok(ffi::bit_vect_bits_in_common(
        bv1.bit_vect(),
        bv2.bit_vect(),
        false,
    )?)
}

/// The vector folded in half, by OR. Mirrors `cDataStructs.FoldFingerprint`.
pub fn fold_fingerprint<B: BitVect>(bv: &B) -> Result<B> {
    fold_fingerprint_with_fold_factor(bv, 2)
}

/// The vector folded `fold_factor` times, which must be more than 0 and
/// less than its size. Mirrors `cDataStructs.FoldFingerprint(bv,
/// foldFactor)`.
pub fn fold_fingerprint_with_fold_factor<B: BitVect>(bv: &B, fold_factor: u32) -> Result<B> {
    Ok(B::from_ffi(ffi::bit_vect_fold(bv.bit_vect(), fold_factor)?))
}

/// Whether every bit on in `probe` is on in `ref_` (N8). Mirrors
/// `cDataStructs.AllProbeBitsMatch`.
pub fn all_probe_bits_match<B: BitVect>(probe: &B, ref_: &B) -> Result<bool> {
    Ok(ffi::bit_vect_all_probe_bits_match(
        probe.bit_vect(),
        ref_.bit_vect(),
    )?)
}

/// As [`all_probe_bits_match`], against the pickle of a vector of the same
/// type (A5).
pub fn all_probe_bits_match_bytes<B: BitVect>(probe: &B, ref_: &[u8]) -> Result<bool> {
    Ok(ffi::bit_vect_all_probe_bits_match_pkl(
        probe.bit_vect(),
        ref_,
    )?)
}
