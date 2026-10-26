// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FeatMaps`, a Python package: feature maps, sets
//! of weighted feature points that score other features by family and
//! distance.

pub mod feat_map_parser;
pub mod feat_map_point;
pub mod feat_map_utils;
// Python's module has its package's name (CONVENTIONS.md N1).
#[allow(clippy::module_inception)]
pub mod feat_maps;
