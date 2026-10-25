// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.MolStandardize.rdMolStandardize`: molecular
//! standardization, after MolVS.
//!
//! RDKit's multi-molecule forms (`CleanupInPlace(mols, numThreads)` and the
//! like) are the `_mols` functions here. Their threads are the wrapper's,
//! not RDKit's: RDKit's stride by the requested count rather than the
//! threads they run (so 0 loops forever, and a negative count has two
//! threads standardize one molecule at once), and an exception in one ends
//! the process.

mod charge;
mod cleanup;
mod fragment;
mod metal;
mod normalize;
mod pipeline;
mod tautomer;
mod validate;

pub use charge::*;
pub use cleanup::*;
pub use fragment::*;
pub use metal::*;
pub use normalize::*;
pub use pipeline::*;
pub use tautomer::*;
pub use validate::*;
