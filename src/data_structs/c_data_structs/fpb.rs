// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::path::Path;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::data_structs::c_data_structs::ExplicitBitVect;
use crate::error::to_string;
use crate::ffi::fpb::ffi;
use crate::util::path_bytes;
use crate::{Error, Result};

/// Reads fingerprints from an FPB file. Mirrors `cDataStructs.FPBReader`.
///
/// Call `init` before anything else; until then, the queries fail.
/// Fingerprints cross as their bytes (`get_bytes`, or `bit_vect_to_binary_text`
/// of an `ExplicitBitVect`), which must have the file's length.
pub struct FPBReader {
    /// Never null.
    ptr: UniquePtr<ffi::FPBReader>,
}

impl FPBReader {
    fn raw(&self) -> &ffi::FPBReader {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::FPBReader> {
        self.ptr.pin_mut()
    }

    /// A reader for `filename`, which it reads and checks now and parses in
    /// `init`. RDKit trusts every size and offset in the file, so the reader
    /// checks them first: a malformed file is `Error::Value`, and so is one
    /// RDKit would misread (8-byte id offsets, or 4 GiB of fingerprints or
    /// more). Mirrors `FPBReader(filename)`.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<FPBReader> {
        FPBReader::new_with_lazy(filename, false)
    }

    /// When `lazy`, `init` leaves the fingerprints undecoded. The reader
    /// keeps the bytes it checked rather than reading the file again, so a
    /// lazy reader holds the whole file in memory. Mirrors
    /// `FPBReader(filename, lazy)`.
    pub fn new_with_lazy<P: AsRef<Path>>(filename: P, lazy: bool) -> Result<FPBReader> {
        let ptr = ffi::fpb_reader_new(path_bytes(filename.as_ref()), lazy)?;
        Ok(FPBReader { ptr })
    }

    /// Reads the file's header, and its fingerprints unless lazy. Mirrors
    /// `FPBReader.Init`.
    pub fn init(&mut self) -> Result<()> {
        Ok(ffi::fpb_reader_init(self.raw_mut())?)
    }

    /// The number of fingerprints. Mirrors `len(FPBReader)`.
    pub fn len(&self) -> Result<u32> {
        Ok(ffi::fpb_reader_length(self.raw())?)
    }

    /// Whether the file has no fingerprints.
    pub fn is_empty(&self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// Mirrors `FPBReader.GetNumBits`.
    pub fn get_num_bits(&self) -> Result<u32> {
        Ok(ffi::fpb_reader_get_num_bits(self.raw())?)
    }

    /// Mirrors `FPBReader.GetFP`.
    pub fn get_fp(&self, idx: u32) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(ffi::fpb_reader_get_fp(
            self.raw(),
            idx,
        )?))
    }

    /// A fingerprint and its id. Mirrors `FPBReader[which]`.
    pub fn get(&self, which: u32) -> Result<(ExplicitBitVect, String)> {
        Ok((self.get_fp(which)?, self.get_id(which)?))
    }

    /// A fingerprint's bytes. Mirrors `FPBReader.GetBytes`.
    pub fn get_bytes(&self, which: u32) -> Result<Vec<u8>> {
        Ok(ffi::fpb_reader_get_bytes(self.raw(), which)?)
    }

    /// Mirrors `FPBReader.GetId`.
    pub fn get_id(&self, idx: u32) -> Result<String> {
        to_string(ffi::fpb_reader_get_id(self.raw(), idx)?)
    }

    /// The Tanimoto similarity of a fingerprint to `bytes`. Mirrors
    /// `FPBReader.GetTanimoto`.
    pub fn get_tanimoto(&self, which: u32, bytes: &[u8]) -> Result<f64> {
        Ok(ffi::fpb_reader_get_tanimoto(self.raw(), which, bytes)?)
    }

    /// The Tversky similarity of a fingerprint to `bytes`. Mirrors
    /// `FPBReader.GetTversky`.
    pub fn get_tversky(&self, which: u32, bytes: &[u8], ca: f64, cb: f64) -> Result<f64> {
        Ok(ffi::fpb_reader_get_tversky(
            self.raw(),
            which,
            bytes,
            ca,
            cb,
        )?)
    }

    /// The (similarity, index) of each fingerprint at least 0.7 similar to
    /// `bv`, most similar first. Mirrors `FPBReader.GetTanimotoNeighbors`.
    pub fn get_tanimoto_neighbors(&self, bv: &[u8]) -> Result<Vec<(f64, u32)>> {
        self.get_tanimoto_neighbors_with_threshold(bv, 0.7)
    }

    /// Mirrors `FPBReader.GetTanimotoNeighbors(bv, threshold)`.
    pub fn get_tanimoto_neighbors_with_threshold(
        &self,
        bv: &[u8],
        threshold: f64,
    ) -> Result<Vec<(f64, u32)>> {
        Ok(
            ffi::fpb_reader_get_tanimoto_neighbors(self.raw(), bv, threshold)?
                .into_iter()
                .map(|nbr| (nbr.similarity, nbr.idx))
                .collect(),
        )
    }

    /// As `get_tanimoto_neighbors`, by Tversky similarity. Mirrors
    /// `FPBReader.GetTverskyNeighbors`.
    pub fn get_tversky_neighbors(&self, bv: &[u8], ca: f64, cb: f64) -> Result<Vec<(f64, u32)>> {
        self.get_tversky_neighbors_with_threshold(bv, ca, cb, 0.7)
    }

    /// Mirrors `FPBReader.GetTverskyNeighbors(bv, ca, cb, threshold)`.
    pub fn get_tversky_neighbors_with_threshold(
        &self,
        bv: &[u8],
        ca: f64,
        cb: f64,
        threshold: f64,
    ) -> Result<Vec<(f64, u32)>> {
        Ok(
            ffi::fpb_reader_get_tversky_neighbors(self.raw(), bv, ca, cb, threshold)?
                .into_iter()
                .map(|nbr| (nbr.similarity, nbr.idx))
                .collect(),
        )
    }

    /// The fingerprints with every bit of `bv` on. Mirrors
    /// `FPBReader.GetContainingNeighbors`.
    pub fn get_containing_neighbors(&self, bv: &[u8]) -> Result<Vec<u32>> {
        Ok(ffi::fpb_reader_get_containing_neighbors(self.raw(), bv)?)
    }
}

impl fmt::Debug for FPBReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FPBReader")
            .field("len", &self.len().ok())
            .finish()
    }
}

/// Parameters for `MultiFPBReader`'s neighbor searches, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MultiFPBReaderNeighborsParams {
    pub threshold: f64,
    pub num_threads: u32,
}

impl Default for MultiFPBReaderNeighborsParams {
    fn default() -> MultiFPBReaderNeighborsParams {
        MultiFPBReaderNeighborsParams {
            threshold: 0.7,
            num_threads: 1,
        }
    }
}

/// Searches several FPB files with fingerprints of one size. Mirrors
/// `cDataStructs.MultiFPBReader`.
///
/// It owns its readers, where Python keeps each alive while the
/// MultiFPBReader is.
pub struct MultiFPBReader {
    // Declared first, so it drops before the readers it points into.
    ptr: UniquePtr<ffi::MultiFPBReader>,
    readers: Vec<FPBReader>,
}

impl MultiFPBReader {
    fn raw(&self) -> &ffi::MultiFPBReader {
        self.ptr.as_ref().unwrap()
    }

    /// Mirrors `MultiFPBReader()`.
    pub fn new() -> MultiFPBReader {
        MultiFPBReader::new_with_init_on_search(false)
    }

    /// When `init_on_search`, the readers are initialized in each search's
    /// threads. Mirrors `MultiFPBReader(initOnSearch)`.
    pub fn new_with_init_on_search(init_on_search: bool) -> MultiFPBReader {
        let ptr = crate::error::expect_ok(ffi::multi_fpb_reader_new(init_on_search));
        MultiFPBReader {
            ptr,
            readers: Vec::new(),
        }
    }

    /// Adds a reader, returning the number of readers. Mirrors
    /// `MultiFPBReader.AddReader`.
    pub fn add_reader(&mut self, mut rdr: FPBReader) -> u32 {
        // The reader's C++ object stays where it is when the reader moves.
        let count = crate::error::expect_ok(ffi::multi_fpb_reader_add_reader(
            self.ptr.pin_mut(),
            rdr.raw_mut(),
        ));
        self.readers.push(rdr);
        count
    }

    /// Mirrors `MultiFPBReader.GetReader`.
    pub fn get_reader(&self, which: u32) -> Result<&FPBReader> {
        self.readers.get(which as usize).ok_or(Error::Index {
            index: which,
            len: self.len(),
        })
    }

    /// Initializes the readers. Mirrors `MultiFPBReader.Init`.
    pub fn init(&mut self) -> Result<()> {
        Ok(ffi::multi_fpb_reader_init(self.ptr.pin_mut())?)
    }

    /// The number of readers. Mirrors `len(MultiFPBReader)`.
    pub fn len(&self) -> u32 {
        u32::try_from(self.readers.len()).expect("fewer than 2^32 readers")
    }

    /// Whether there are no readers.
    pub fn is_empty(&self) -> bool {
        self.readers.is_empty()
    }

    /// Mirrors `MultiFPBReader.GetNumBits`.
    pub fn get_num_bits(&self) -> Result<u32> {
        Ok(ffi::multi_fpb_reader_get_num_bits(self.raw())?)
    }

    /// The (similarity, index, reader) of each fingerprint at least 0.7
    /// similar to `bv`. Mirrors `MultiFPBReader.GetTanimotoNeighbors`.
    pub fn get_tanimoto_neighbors(&self, bv: &[u8]) -> Result<Vec<(f64, u32, u32)>> {
        self.get_tanimoto_neighbors_with_params(bv, &Default::default())
    }

    /// Mirrors `MultiFPBReader.GetTanimotoNeighbors(bv, threshold,
    /// numThreads)`.
    pub fn get_tanimoto_neighbors_with_params(
        &self,
        bv: &[u8],
        params: &MultiFPBReaderNeighborsParams,
    ) -> Result<Vec<(f64, u32, u32)>> {
        Ok(ffi::multi_fpb_reader_get_tanimoto_neighbors(
            self.raw(),
            bv,
            params.threshold,
            params.num_threads,
        )?
        .into_iter()
        .map(|nbr| (nbr.similarity, nbr.idx, nbr.reader))
        .collect())
    }

    /// As `get_tanimoto_neighbors`, by Tversky similarity. Mirrors
    /// `MultiFPBReader.GetTverskyNeighbors`.
    pub fn get_tversky_neighbors(
        &self,
        bv: &[u8],
        ca: f64,
        cb: f64,
    ) -> Result<Vec<(f64, u32, u32)>> {
        self.get_tversky_neighbors_with_params(bv, ca, cb, &Default::default())
    }

    /// Mirrors `MultiFPBReader.GetTverskyNeighbors(bv, ca, cb, threshold,
    /// numThreads)`.
    pub fn get_tversky_neighbors_with_params(
        &self,
        bv: &[u8],
        ca: f64,
        cb: f64,
        params: &MultiFPBReaderNeighborsParams,
    ) -> Result<Vec<(f64, u32, u32)>> {
        Ok(ffi::multi_fpb_reader_get_tversky_neighbors(
            self.raw(),
            bv,
            ca,
            cb,
            params.threshold,
            params.num_threads,
        )?
        .into_iter()
        .map(|nbr| (nbr.similarity, nbr.idx, nbr.reader))
        .collect())
    }

    /// The (index, reader) of each fingerprint with every bit of `bv` on.
    /// Mirrors `MultiFPBReader.GetContainingNeighbors`.
    pub fn get_containing_neighbors(&self, bv: &[u8]) -> Result<Vec<(u32, u32)>> {
        self.get_containing_neighbors_with_num_threads(bv, 1)
    }

    /// Mirrors `MultiFPBReader.GetContainingNeighbors(bv, numThreads)`.
    pub fn get_containing_neighbors_with_num_threads(
        &self,
        bv: &[u8],
        num_threads: u32,
    ) -> Result<Vec<(u32, u32)>> {
        Ok(
            ffi::multi_fpb_reader_get_containing_neighbors(self.raw(), bv, num_threads)?
                .into_iter()
                .map(|nbr| (nbr.idx, nbr.reader))
                .collect(),
        )
    }
}

impl Default for MultiFPBReader {
    fn default() -> MultiFPBReader {
        MultiFPBReader::new()
    }
}

impl fmt::Debug for MultiFPBReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MultiFPBReader")
            .field("readers", &self.readers)
            .finish()
    }
}
