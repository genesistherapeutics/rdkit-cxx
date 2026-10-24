// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Rust bindings for RDKit's C++ API, mirroring its Python wrappers.

pub mod chem;
pub mod data_structs;
pub mod error;
mod ffi;
mod flags;
pub mod force_field;
pub mod geometry;
pub mod rd_base;
mod util;

pub use error::{Error, Result};
