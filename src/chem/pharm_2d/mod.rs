// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm2D`, a Python package: 2D (topological)
//! pharmacophore fingerprints.
//!
//! `rdkit.Chem.Pharm2D.LazyGenerator` is unfinished in Python (importing it
//! raises `NotImplementedError`), so it is not ported.

pub mod generate;
pub mod gobbi_pharm_2d;
mod init;
pub mod matcher;
pub mod sig_factory;
pub mod utils;

pub use init::*;
