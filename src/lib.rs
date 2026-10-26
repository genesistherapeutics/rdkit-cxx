// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Rust bindings for RDKit's C++ API, with Rust ports of much of the
//! functionality RDKit currently offers only in Python. Module paths follow
//! Python's:
//! `rdkit.Chem.rdmolfiles` is [`chem::rdmolfiles`] and `rdkit.Chem.QED` is
//! [`chem::qed`].

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
pub mod rd_base;
pub mod rd_logger;
pub mod rd_paths;
pub mod rd_random;
pub mod sim_div_filters;
mod util;
pub mod utils;

pub use error::{Error, Result};
