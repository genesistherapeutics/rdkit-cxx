// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.DataStructs`: bit vectors ([`ExplicitBitVect`],
//! [`SparseBitVect`]), sparse integer vectors for count fingerprints,
//! discrete and real value vectors, FPB files, and the similarity metrics
//! on them, with the Python package's own functions ported.
//!
//! ```no_run
//! use rdkit_cxx::data_structs::{ExplicitBitVect, dice_similarity, tanimoto_similarity};
//!
//! let mut a = ExplicitBitVect::new(64);
//! let mut b = ExplicitBitVect::new(64);
//! for bit in [1, 5, 9] {
//!     a.set_bit(bit)?;
//! }
//! for bit in [5, 9, 12] {
//!     b.set_bit(bit)?;
//! }
//! println!("{} {}", tanimoto_similarity(&a, &b)?, dice_similarity(&a, &b)?);
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

pub mod bit_ensemble;
pub mod bit_utils;
pub mod c_data_structs;
pub(crate) mod init;
pub mod lazy_signature;
pub mod top_n_container;
pub mod vect_collection;

pub use c_data_structs::*;
pub use init::*;
