// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.DataStructs.BitEnsemble`, a Python module: a
//! collection of bits to score bit vectors against.

use crate::Result;
use crate::data_structs::BitVect;
use crate::data_structs::init::{get_bit, on_bits};

/// A collection of bits, scoring bit vectors by how many of them they set.
/// Mirrors `BitEnsemble.BitEnsemble`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BitEnsemble {
    bits: Vec<u32>,
}

impl BitEnsemble {
    /// An ensemble of the bits `bits`. Mirrors `BitEnsemble(bits)`.
    pub fn new_with_bits(bits: Vec<u32>) -> BitEnsemble {
        BitEnsemble { bits }
    }

    /// Replaces the bits. Mirrors `BitEnsemble.SetBits`.
    pub fn set_bits(&mut self, bits: Vec<u32>) {
        self.bits = bits;
    }

    /// Adds a bit. Mirrors `BitEnsemble.AddBit`.
    pub fn add_bit(&mut self, bit: u32) {
        self.bits.push(bit);
    }

    /// The bits, in the order they were added. Mirrors
    /// `BitEnsemble.GetBits`.
    pub fn get_bits(&self) -> &[u32] {
        &self.bits
    }

    /// The number of bits. Mirrors `BitEnsemble.GetNumBits`.
    pub fn get_num_bits(&self) -> usize {
        self.bits.len()
    }

    /// How many of the ensemble's bits are among `other`'s on bits, a bit
    /// listed twice counting twice. Mirrors `BitEnsemble.ScoreWithOnBits`.
    pub fn score_with_on_bits<B: BitVect>(&self, other: &B) -> usize {
        let obl = on_bits(other);
        self.bits
            .iter()
            .filter(|bit| obl.binary_search(bit).is_ok())
            .count()
    }

    /// How many of the ensemble's bits `other` sets, looking each up.
    /// Mirrors `BitEnsemble.ScoreWithIndex`.
    ///
    /// `Error::Index` for a bit past the end of `other`, as Python raises.
    pub fn score_with_index<B: BitVect>(&self, other: &B) -> Result<usize> {
        let mut cnt = 0;
        for &bit in &self.bits {
            if get_bit(other, bit)? {
                cnt += 1;
            }
        }
        Ok(cnt)
    }
}
