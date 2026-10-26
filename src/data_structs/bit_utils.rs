// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.DataStructs.BitUtils`, a pure-Python module: utilities for
//! bit vectors.

use crate::Result;
use crate::data_structs::BitVect;
use crate::data_structs::init::{get_bit, new_like, set_bit};

/// A vector of the same type holding just the bits `bits_to_keep`, in that
/// order: bit `i` of the result is bit `bits_to_keep[i]` of `bv`. Mirrors
/// `BitUtils.ConstructEnsembleBV`.
///
/// `Error::Index` for a bit past the end of `bv`, as Python raises;
/// `Error::Value` for more than `i32::MAX` bits.
pub fn construct_ensemble_bv<B: BitVect + 'static>(bv: &B, bits_to_keep: &[u32]) -> Result<B> {
    let final_size = u32::try_from(bits_to_keep.len()).unwrap_or(u32::MAX);
    let mut res: B = new_like(final_size)?;
    for (i, &bit) in bits_to_keep.iter().enumerate() {
        if get_bit(bv, bit)? {
            set_bit(&mut res, i as u32)?;
        }
    }
    Ok(res)
}
