// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.ML.InfoTheory.rdInfoTheory`: the entropy, information
//! gain, and chi-squared statistics of count matrices, a ranker of
//! fingerprint bits by them, and a matrix of how often bits are set
//! together.
//!
//! A count matrix has a row per value of a variable and a column per
//! class, each entry the count of the examples of that value and class.
//! RDKit computes in `f64`, as Python's does for every NumPy type but
//! `float32`.
//!
//! ```no_run
//! use rdkit_cxx::data_structs::ExplicitBitVect;
//! use rdkit_cxx::ml::info_theory::InfoBitRanker;
//!
//! let mut ranker = InfoBitRanker::new(4, 2);
//! let mut bv = ExplicitBitVect::new(4);
//! bv.set_bit(3).unwrap();
//! ranker.accumulate_votes(&bv, 0).unwrap();
//! let top = ranker.get_top_n(1).unwrap();
//! assert_eq!(top[0][0], 3.0);
//! ```

use std::fmt;
use std::path::Path;

use cxx::UniquePtr;

use crate::chem::rdchem::enums::rdkit_enum;
use crate::data_structs::BitVect;
use crate::error::expect_ok;
use crate::ffi::info_theory::ffi;
use crate::util::path_bytes;
use crate::{Error, Result};

/// A rectangular matrix's values by rows, with its shape.
fn flatten(res_arr: &[Vec<f64>]) -> Result<(Vec<f64>, u64, u64)> {
    let cols = res_arr.first().map_or(0, Vec::len);
    if res_arr.iter().any(|row| row.len() != cols) {
        return Err(Error::Value {
            message: "the matrix's rows must have the same length".to_owned(),
        });
    }
    Ok((res_arr.concat(), res_arr.len() as u64, cols as u64))
}

/// The informational entropy, in bits, of counts. Mirrors
/// `rdInfoTheory.InfoEntropy`.
///
/// `Error::Runtime` for no counts, as Python raises; `Error::Value` for
/// more than `i32::MAX`, which RDKit counts in `int`.
pub fn info_entropy(res_arr: &[f64]) -> Result<f64> {
    Ok(ffi::info_entropy(res_arr)?)
}

/// The information gain of a variable, given its count matrix. Mirrors
/// `rdInfoTheory.InfoGain`.
///
/// `Error::Value` for rows of different lengths, more than `i32::MAX`
/// columns, or column sums that do not truncate to `i32`, or whose
/// truncations' sum does not fit one: RDKit sums them in `int` (undefined
/// behavior, and garbage in Python).
pub fn info_gain(res_arr: &[Vec<f64>]) -> Result<f64> {
    let (values, rows, cols) = flatten(res_arr)?;
    Ok(ffi::info_gain(&values, rows, cols)?)
}

/// The chi-squared statistic of a variable, given its count matrix.
/// Mirrors `rdInfoTheory.ChiSquare`.
///
/// `Error::Value` for rows of different lengths, more than `i32::MAX`
/// values, which RDKit indexes in `int`, or row sums that do not truncate
/// to `i32`, or whose truncations' sum does not fit one, as for
/// [`info_gain`].
pub fn chi_square(res_arr: &[Vec<f64>]) -> Result<f64> {
    let (values, rows, cols) = flatten(res_arr)?;
    Ok(ffi::chi_square(&values, rows, cols)?)
}

rdkit_enum! {
    /// What a ranker ranks bits by: entropy (information gain) or the
    /// chi-squared statistic, biased or not. Mirrors
    /// `rdInfoTheory.InfoType`.
    InfoType: ffi::InfoType {
        Entropy = ENTROPY,
        Biasentropy = BIASENTROPY,
        Chisquare = CHISQUARE,
        Biaschisquare = BIASCHISQUARE,
    }
}

/// Ranks the bits of labeled fingerprints by how much they tell about the
/// labels. Mirrors `rdInfoTheory.InfoBitRanker`.
///
/// Each bit's counts are, per class, the fingerprints of that class that
/// set it; RDKit keeps them in 16 bits, so they wrap past 65535 votes, as
/// in Python.
pub struct InfoBitRanker {
    /// Never null.
    ptr: UniquePtr<ffi::InfoBitRanker>,
}

impl fmt::Debug for InfoBitRanker {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("InfoBitRanker").finish_non_exhaustive()
    }
}

impl InfoBitRanker {
    /// A ranker of `n_bits` bits by entropy, for `n_classes` classes.
    /// Mirrors `InfoBitRanker(nBits, nClasses)`.
    ///
    /// # Panics
    ///
    /// For more than `i32::MAX` bits, which RDKit numbers in `int`, or more
    /// than 32768 classes, whose counts RDKit adds in `int` (U3).
    pub fn new(n_bits: u32, n_classes: u32) -> InfoBitRanker {
        InfoBitRanker::new_with_info_type(n_bits, n_classes, InfoType::Entropy)
    }

    /// As [`new`](Self::new), ranking by `info_type`. Mirrors
    /// `InfoBitRanker(nBits, nClasses, infoType)`.
    ///
    /// # Panics
    ///
    /// As [`new`](Self::new).
    pub fn new_with_info_type(n_bits: u32, n_classes: u32, info_type: InfoType) -> InfoBitRanker {
        assert!(
            n_bits <= i32::MAX as u32,
            "an InfoBitRanker ranks at most i32::MAX bits"
        );
        assert!(
            n_classes <= 32768,
            "an InfoBitRanker takes at most 32768 classes"
        );
        let ptr = expect_ok(ffi::info_bit_ranker_new(
            n_bits,
            n_classes,
            info_type.to_ffi(),
        ));
        InfoBitRanker { ptr }
    }

    /// Counts the bits a fingerprint of class `label` sets, except those
    /// the mask leaves out. Mirrors `InfoBitRanker.AccumulateVotes`.
    ///
    /// `Error::Runtime` for a `label` past the classes or a fingerprint of
    /// another size than the ranker's, as Python raises.
    pub fn accumulate_votes<B: BitVect>(&mut self, bit_vect: &B, label: u32) -> Result<()> {
        Ok(ffi::info_bit_ranker_accumulate_votes(
            self.ptr.pin_mut(),
            bit_vect.bit_vect(),
            label,
        )?)
    }

    /// The top `num` bits, ranked, each a row of the bit, its score, and
    /// its count in each class. Mirrors `InfoBitRanker.GetTopN`.
    ///
    /// `Error::Value` for a `num` above the ranker's bits, as Python
    /// raises, for a `num` of 0 (RDKit reads its empty ranking, a segfault
    /// in Python), or for `num * (2 + n_classes)` above `i32::MAX`, which
    /// RDKit indexes in `int`; `Error::Runtime` for a `num` above the
    /// mask's bits, as Python raises.
    pub fn get_top_n(&mut self, num: u32) -> Result<Vec<Vec<f64>>> {
        let ncols = expect_ok(ffi::info_bit_ranker_get_num_classes(&self.ptr)) as usize + 2;
        let values = ffi::info_bit_ranker_get_top_n(self.ptr.pin_mut(), num)?;
        Ok(values.chunks(ncols).map(<[f64]>::to_vec).collect())
    }

    /// The classes the biased types favor: a bit counts only when its
    /// fraction of hits in one of these classes is at least its fraction
    /// in every other class. Mirrors `InfoBitRanker.SetBiasList`.
    ///
    /// `Error::Runtime` for as many classes as the ranker has or more, a
    /// class past them, or a class listed twice, as Python raises (whose
    /// ranker keeps such a list and reads past its classes' fractions).
    pub fn set_bias_list(&mut self, class_list: &[u32]) -> Result<()> {
        Ok(ffi::info_bit_ranker_set_bias_list(
            self.ptr.pin_mut(),
            class_list,
        )?)
    }

    /// Ranks only these bits, from votes counted after this call. Mirrors
    /// `InfoBitRanker.SetMaskBits`.
    ///
    /// `Error::Index` for a bit past the ranker's, which leaves the mask as
    /// it was.
    pub fn set_mask_bits(&mut self, mask_bits: &[u32]) -> Result<()> {
        Ok(ffi::info_bit_ranker_set_mask_bits(
            self.ptr.pin_mut(),
            mask_bits,
        )?)
    }

    /// Writes the bits [`get_top_n`](Self::get_top_n) last ranked, as a
    /// table, to a file. Mirrors `InfoBitRanker.WriteTopBitsToFile`.
    ///
    /// `Error::Os` when the file cannot be opened or written.
    pub fn write_top_bits_to_file<P: AsRef<Path>>(&self, file_name: P) -> Result<()> {
        Ok(ffi::info_bit_ranker_write_top_bits_to_file(
            &self.ptr,
            path_bytes(file_name.as_ref()),
        )?)
    }
}

/// Counts, for each pair of listed bits, the fingerprints that set both.
/// Mirrors `rdInfoTheory.BitCorrMatGenerator`.
pub struct BitCorrMatGenerator {
    /// Never null.
    ptr: UniquePtr<ffi::BitCorrMatGenerator>,
}

impl fmt::Debug for BitCorrMatGenerator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BitCorrMatGenerator")
            .finish_non_exhaustive()
    }
}

impl Default for BitCorrMatGenerator {
    /// A generator with no bits listed. Mirrors `BitCorrMatGenerator()`.
    fn default() -> BitCorrMatGenerator {
        BitCorrMatGenerator {
            ptr: expect_ok(ffi::bit_corr_mat_generator_new()),
        }
    }
}

impl BitCorrMatGenerator {
    /// Lists the bits to correlate, clearing the counts. Mirrors
    /// `BitCorrMatGenerator.SetBitList`.
    ///
    /// `Error::Value` for more than 46341 bits, whose matrix RDKit sizes
    /// in `int`.
    pub fn set_bit_list(&mut self, bit_list: &[u32]) -> Result<()> {
        Ok(ffi::bit_corr_mat_generator_set_bit_list(
            self.ptr.pin_mut(),
            bit_list,
        )?)
    }

    /// Counts a fingerprint: each pair of listed bits it sets. Mirrors
    /// `BitCorrMatGenerator.CollectVotes`.
    ///
    /// `Error::Index` for a listed bit past the fingerprint's end, as
    /// Python raises.
    pub fn collect_votes<B: BitVect>(&mut self, bit_vect: &B) -> Result<()> {
        Ok(ffi::bit_corr_mat_generator_collect_votes(
            self.ptr.pin_mut(),
            bit_vect.bit_vect(),
        )?)
    }

    /// The counts, as the lower triangle of the symmetric matrix by rows:
    /// the count for listed bits `i > j` is at `i * (i - 1) / 2 + j`.
    /// Empty before a bit list is set. Mirrors
    /// `BitCorrMatGenerator.GetCorrMatrix`.
    pub fn get_corr_matrix(&self) -> Vec<f64> {
        expect_ok(ffi::bit_corr_mat_generator_get_corr_matrix(&self.ptr))
    }
}
