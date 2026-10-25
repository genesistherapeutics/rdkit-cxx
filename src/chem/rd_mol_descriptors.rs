// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolDescriptors`: molecular descriptors.
//!
//! RDKit caches many descriptors, and their per-atom contributions, in a
//! molecule's computed properties, and Python's `force` arguments say
//! whether to read those caches back. The crate never keeps them between
//! calls (CONVENTIONS.md U8), so `force` is not bound and every result is
//! computed afresh.

mod dclv;
mod descriptors;
mod descriptors_3d;
mod fingerprints;
mod properties;
mod usr;

pub use dclv::*;
pub use descriptors::*;
pub use descriptors_3d::*;
pub use fingerprints::*;
pub use properties::*;
pub use usr::*;
