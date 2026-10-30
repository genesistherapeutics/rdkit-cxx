// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Subshape`, a Python package: shape skeletons of
//! molecules and alignments between them.
//!
//! `rdkit.Chem.Subshape.demoCombined` is a demonstration script that reads
//! files from RDKit's source tree, so it is not ported.

pub mod builder_utils;
pub mod subshape_aligner;
pub mod subshape_builder;
pub mod subshape_objects;
