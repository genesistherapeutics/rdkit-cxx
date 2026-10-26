// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.AtomPairs`, a Python package: atom-pair,
//! topological-torsion, and Sheridan's physicochemical-property
//! fingerprints, and the helpers that explain their codes.

pub mod pairs;
pub mod sheridan;
pub mod torsions;
pub mod utils;
