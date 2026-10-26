// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.ML`, RDKit's Python machine-learning package, ported to
//! Rust: Butina and hierarchical clustering ([`cluster`]), data sets and
//! their files ([`data`]), molecular descriptor calculators
//! ([`descriptors`]), information theory and bit ranking
//! ([`info_theory`]), and enrichment and ROC scoring ([`scoring`]).
//!
//! ```no_run
//! use rdkit_cxx::ml::cluster::butina::cluster_data;
//!
//! let points = vec![vec![0.0, 0.0], vec![0.1, 0.0], vec![5.0, 5.0]];
//! let clusters = cluster_data(&points, points.len(), 0.5)?;
//! assert_eq!(clusters.len(), 2);
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

pub mod cluster;
pub mod data;
pub mod descriptors;
pub mod files;
pub mod info_theory;
pub mod scoring;
pub mod slt;
