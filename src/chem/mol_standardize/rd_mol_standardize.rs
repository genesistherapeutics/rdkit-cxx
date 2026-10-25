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
//!
//! Most of what is here matches substructures, and every substructure match
//! swaps the process's SIGINT handler through RDKit's unsynchronized
//! statics, so those calls hold one process-wide lock for the whole RDKit
//! call (CONVENTIONS.md U4): the standardizations and parents, normalizing,
//! reionizing, uncharging, removing fragments, disconnecting metals,
//! `FragmentValidation` and `MolVSValidation`, tautomer enumeration and
//! canonicalization, the substructure scores, and the pipeline. Calls on
//! several threads take turns, the `_mols` functions' threads among them, so
//! those threads gain little. Choosing the largest fragment, the ring and
//! hetero-H scores, picking a tautomer by a closure's score, and the other
//! validations run without the lock. Closures RDKit calls (validation
//! methods, tautomer enumerator callbacks, and score functions) run without
//! it too, so they may call anything here.

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
