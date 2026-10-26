// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.DataStructs.LazySignature`, a Python module: a
//! signature whose bits are computed when first read.

use std::collections::HashMap;
use std::fmt;

use crate::{Error, Result};

/// A signature of `sig_size` values, each computed by a function of its
/// index the first time it is read, and kept. Mirrors
/// `LazySignature.LazySig`.
pub struct LazySig<T = bool, F = fn(usize) -> T> {
    compute_func: F,
    size: usize,
    cache: HashMap<usize, T>,
}

impl<T, F> fmt::Debug for LazySig<T, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LazySig")
            .field("size", &self.size)
            .finish_non_exhaustive()
    }
}

impl<T: Clone, F: FnMut(usize) -> T> LazySig<T, F> {
    /// A signature of `sig_size` values computed by `compute_func`. Mirrors
    /// `LazySig(computeFunc, sigSize)`.
    ///
    /// `Error::Value` for a size of zero, as Python raises.
    pub fn new(compute_func: F, sig_size: usize) -> Result<LazySig<T, F>> {
        if sig_size == 0 {
            return Err(Error::Value {
                message: "zero size".to_owned(),
            });
        }
        Ok(LazySig {
            compute_func,
            size: sig_size,
            cache: HashMap::new(),
        })
    }

    /// The signature's size. Mirrors `len(sig)`.
    pub fn len(&self) -> usize {
        self.size
    }

    /// Always `false`: a signature has at least one value.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// The value at `which`, counting from the end when negative, computed
    /// on first use. Mirrors `sig[which]`.
    ///
    /// `Error::Index` outside the signature, and, as in Python, at index 0
    /// (and so at `-len`), which Python refuses too.
    pub fn get(&mut self, which: i64) -> Result<T> {
        let size = i64::try_from(self.size).unwrap_or(i64::MAX);
        let index = if which < 0 { size + which } else { which };
        if index <= 0 || index >= size {
            return Err(Error::index(which, self.size));
        }
        let which = index as usize;
        if let Some(v) = self.cache.get(&which) {
            return Ok(v.clone());
        }
        let v = (self.compute_func)(which);
        self.cache.insert(which, v.clone());
        Ok(v)
    }
}
