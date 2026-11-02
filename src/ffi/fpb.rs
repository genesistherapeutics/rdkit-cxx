// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A fingerprint's similarity and index.
    struct FPBNeighbor {
        similarity: f64,
        idx: u32,
    }

    /// A fingerprint's similarity, index, and reader.
    struct MultiFPBNeighbor {
        similarity: f64,
        idx: u32,
        reader: u32,
    }

    /// A fingerprint's index and reader.
    struct MultiFPBContaining {
        idx: u32,
        reader: u32,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/fpb.h");

        type FPBReader;
        type MultiFPBReader;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;

        fn fpb_reader_new(path: &[u8], lazy: bool) -> Result<UniquePtr<FPBReader>>;
        fn fpb_reader_init(reader: Pin<&mut FPBReader>) -> Result<()>;
        fn fpb_reader_length(reader: &FPBReader) -> Result<u32>;
        fn fpb_reader_get_num_bits(reader: &FPBReader) -> Result<u32>;
        fn fpb_reader_get_fp(reader: &FPBReader, idx: u32) -> Result<UniquePtr<BitVect>>;
        fn fpb_reader_get_bytes(reader: &FPBReader, which: u32) -> Result<Vec<u8>>;
        fn fpb_reader_get_id(reader: &FPBReader, idx: u32) -> Result<UniquePtr<CxxString>>;
        fn fpb_reader_get_tanimoto(reader: &FPBReader, which: u32, bytes: &[u8]) -> Result<f64>;
        fn fpb_reader_get_tversky(
            reader: &FPBReader,
            which: u32,
            bytes: &[u8],
            ca: f64,
            cb: f64,
        ) -> Result<f64>;
        fn fpb_reader_get_tanimoto_neighbors(
            reader: &FPBReader,
            bv: &[u8],
            threshold: f64,
        ) -> Result<Vec<FPBNeighbor>>;
        fn fpb_reader_get_tversky_neighbors(
            reader: &FPBReader,
            bv: &[u8],
            ca: f64,
            cb: f64,
            threshold: f64,
        ) -> Result<Vec<FPBNeighbor>>;
        fn fpb_reader_get_containing_neighbors(reader: &FPBReader, bv: &[u8]) -> Result<Vec<u32>>;

        fn multi_fpb_reader_new(init_on_search: bool) -> Result<UniquePtr<MultiFPBReader>>;
        fn multi_fpb_reader_add_reader(
            multi: Pin<&mut MultiFPBReader>,
            reader: Pin<&mut FPBReader>,
        ) -> Result<u32>;
        fn multi_fpb_reader_init(multi: Pin<&mut MultiFPBReader>) -> Result<()>;
        fn multi_fpb_reader_get_num_bits(multi: &MultiFPBReader) -> Result<u32>;
        fn multi_fpb_reader_get_tanimoto_neighbors(
            multi: &MultiFPBReader,
            bv: &[u8],
            threshold: f64,
            num_threads: u32,
        ) -> Result<Vec<MultiFPBNeighbor>>;
        fn multi_fpb_reader_get_tversky_neighbors(
            multi: &MultiFPBReader,
            bv: &[u8],
            ca: f64,
            cb: f64,
            threshold: f64,
            num_threads: u32,
        ) -> Result<Vec<MultiFPBNeighbor>>;
        fn multi_fpb_reader_get_containing_neighbors(
            multi: &MultiFPBReader,
            bv: &[u8],
            num_threads: u32,
        ) -> Result<Vec<MultiFPBContaining>>;
    }
}

// SAFETY: a reader owns its file stream and buffers, with no thread affinity.
// It is not Sync: lazy readers seek their stream in const methods.
unsafe impl Send for ffi::FPBReader {}
unsafe impl Send for ffi::MultiFPBReader {}
