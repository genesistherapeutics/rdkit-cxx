// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.ChemUtils.DescriptorUtilities`, a Python module:
//! a wrapper that makes a descriptor function returning a vector seem one
//! function per entry.
//!
//! Python's `setDescriptorVersion` decorator and `VectorDescriptorNamespace`
//! class set attributes and fill namespaces of Python functions, which Rust
//! has no counterpart for.

use std::fmt;

use crate::chem::Mol;
use crate::{Error, Result};

/// A descriptor function that returns a vector.
type VectorFn = dyn Fn(&Mol) -> Result<Vec<f64>> + Send + Sync;

/// A descriptor function that returns a vector, with a name for each entry.
/// Python's wrapper adds a function per name to a namespace, each calling
/// [`call_desc`](Self::call_desc) with its index, and caches the vector on
/// the molecule; this computes it on each call. Mirrors
/// `DescriptorUtilities.VectorDescriptorWrapper`.
pub struct VectorDescriptorWrapper {
    func: Box<VectorFn>,
    names: Vec<String>,
    version: String,
}

impl fmt::Debug for VectorDescriptorWrapper {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("VectorDescriptorWrapper")
            .field("names", &self.names)
            .field("version", &self.version)
            .finish_non_exhaustive()
    }
}

impl VectorDescriptorWrapper {
    /// Wraps `func`, whose entries `names` name, as descriptors of version
    /// `version`. Python's `namespace`, the dictionary it adds a function
    /// per name to, has no counterpart. Mirrors
    /// `VectorDescriptorWrapper(func, names, version, namespace)`.
    pub fn new<F>(func: F, names: Vec<String>, version: &str) -> VectorDescriptorWrapper
    where
        F: Fn(&Mol) -> Result<Vec<f64>> + Send + Sync + 'static,
    {
        VectorDescriptorWrapper {
            func: Box::new(func),
            names,
            version: version.to_owned(),
        }
    }

    /// The names of the vector's entries. Python's `names` attribute.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// The descriptors' version, which Python gives each function as its
    /// `version` attribute.
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Entry `index` of the function's vector for `mol`, or NaN when the
    /// function fails, as in Python. `Error::Index` for an index past the
    /// vector's end, where Python raises `IndexError`. Mirrors
    /// `VectorDescriptorWrapper.call_desc`.
    pub fn call_desc(&self, mol: &Mol, index: u32) -> Result<f64> {
        let Ok(results) = (self.func)(mol) else {
            return Ok(f64::NAN);
        };
        results
            .get(index as usize)
            .copied()
            .ok_or(Error::index(index, results.len()))
    }
}
