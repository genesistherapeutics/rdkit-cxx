// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.MolStandardize`, whose compiled module is
//! `rdMolStandardize` and whose own Python function (`ReorderTautomers`) is
//! ported here.

mod init;
pub mod rd_mol_standardize;

pub use init::*;
