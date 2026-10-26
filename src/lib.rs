// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Rust bindings for RDKit's C++ API, mirroring its Python wrappers.

pub mod avalon;
pub mod chem;
pub mod data_manip;
pub mod data_structs;
pub mod distance_geometry;
pub mod error;
mod ffi;
mod flags;
pub mod force_field;
pub mod geometry;
pub mod ml;
pub mod numerics;
pub mod py_random;
pub mod rd_base;
pub mod sim_div_filters;
mod util;

pub use error::{Error, Result};
