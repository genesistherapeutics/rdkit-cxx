// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDInfoTheory::InfoBitRanker::InfoType`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum InfoType {
        ENTROPY = 1,
        BIASENTROPY,
        CHISQUARE,
        BIASCHISQUARE,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/info_theory.h");

        type InfoType;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type InfoBitRanker;
        type BitCorrMatGenerator;

        fn info_entropy(res_arr: &[f64]) -> Result<f64>;
        fn info_gain(values: &[f64], rows: u64, cols: u64) -> Result<f64>;
        fn chi_square(values: &[f64], rows: u64, cols: u64) -> Result<f64>;

        fn info_bit_ranker_new(
            n_bits: u32,
            n_classes: u32,
            info_type: InfoType,
        ) -> Result<UniquePtr<InfoBitRanker>>;
        fn info_bit_ranker_get_num_classes(ranker: &InfoBitRanker) -> Result<u32>;
        fn info_bit_ranker_accumulate_votes(
            ranker: Pin<&mut InfoBitRanker>,
            bit_vect: &BitVect,
            label: u32,
        ) -> Result<()>;
        fn info_bit_ranker_set_bias_list(
            ranker: Pin<&mut InfoBitRanker>,
            classes: &[u32],
        ) -> Result<()>;
        fn info_bit_ranker_set_mask_bits(
            ranker: Pin<&mut InfoBitRanker>,
            bits: &[u32],
        ) -> Result<()>;
        fn info_bit_ranker_get_top_n(ranker: Pin<&mut InfoBitRanker>, num: u32)
        -> Result<Vec<f64>>;
        fn info_bit_ranker_write_top_bits_to_file(
            ranker: &InfoBitRanker,
            path: &[u8],
        ) -> Result<()>;

        fn bit_corr_mat_generator_new() -> Result<UniquePtr<BitCorrMatGenerator>>;
        fn bit_corr_mat_generator_set_bit_list(
            generator: Pin<&mut BitCorrMatGenerator>,
            bits: &[u32],
        ) -> Result<()>;
        fn bit_corr_mat_generator_collect_votes(
            generator: Pin<&mut BitCorrMatGenerator>,
            bit_vect: &BitVect,
        ) -> Result<()>;
        fn bit_corr_mat_generator_get_corr_matrix(
            generator: &BitCorrMatGenerator,
        ) -> Result<Vec<f64>>;
    }
}

// SAFETY: the ranker and the generator are plain data, with no globals or
// thread affinity, and their const methods write nothing.
unsafe impl Send for ffi::InfoBitRanker {}
unsafe impl Sync for ffi::InfoBitRanker {}
unsafe impl Send for ffi::BitCorrMatGenerator {}
unsafe impl Sync for ffi::BitCorrMatGenerator {}
