// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.DataStructs.VectCollection`, a Python module: a
//! collection of bit vectors that acts as their union.

use crate::data_structs::ExplicitBitVect;
use crate::{Error, Result};

/// A collection of keyed bit vectors that reads as their union (the OR of
/// them all). Mirrors `VectCollection.VectCollection`.
///
/// Python's pickling hooks are not ported; `Clone` copies a collection, as
/// Python's `copy.copy` does.
#[derive(Clone, Debug)]
pub struct VectCollection {
    /// The vectors in the order their keys were first added, as Python's
    /// dict keeps them.
    vects: Vec<(u32, ExplicitBitVect)>,
    or_vect: Option<ExplicitBitVect>,
    num_bits: i64,
    need_reset: bool,
}

impl Default for VectCollection {
    /// An empty collection. Mirrors `VectCollection()`.
    fn default() -> VectCollection {
        VectCollection {
            vects: Vec::new(),
            or_vect: None,
            num_bits: -1,
            need_reset: true,
        }
    }
}

impl VectCollection {
    fn reset_if_needed(&mut self) -> Result<()> {
        if self.need_reset {
            self.reset()?;
        }
        Ok(())
    }

    fn union(&self) -> Result<&ExplicitBitVect> {
        self.or_vect.as_ref().ok_or_else(|| Error::Value {
            message: "the collection holds no vectors".to_owned(),
        })
    }

    /// The union of the vectors, `None` when there are none. Mirrors
    /// `VectCollection.GetOrVect`.
    ///
    /// `Error::Value` for vectors of different sizes, as Python raises.
    pub fn get_or_vect(&mut self) -> Result<Option<&ExplicitBitVect>> {
        self.reset_if_needed()?;
        Ok(self.or_vect.as_ref())
    }

    /// As [`get_or_vect`](Self::get_or_vect). Mirrors the
    /// `VectCollection.orVect` property.
    pub fn or_vect(&mut self) -> Result<Option<&ExplicitBitVect>> {
        self.get_or_vect()
    }

    /// Adds a vector under the key `idx`, replacing any vector already
    /// there (which keeps its place). Mirrors `VectCollection.AddVect`.
    pub fn add_vect(&mut self, idx: u32, vect: ExplicitBitVect) {
        match self.vects.iter_mut().find(|(k, _)| *k == idx) {
            Some(entry) => entry.1 = vect,
            None => self.vects.push((idx, vect)),
        }
        self.need_reset = true;
    }

    /// Recomputes the union, if a change needs it. Mirrors
    /// `VectCollection.Reset`.
    ///
    /// As in Python, with no vectors the union is `None` and the number of
    /// bits keeps its last value. `Error::Value` for vectors of different
    /// sizes, as Python raises.
    pub fn reset(&mut self) -> Result<()> {
        if !self.need_reset {
            return Ok(());
        }
        self.or_vect = None;
        let Some((_, first)) = self.vects.first() else {
            return Ok(());
        };
        let num_bits = first.get_num_bits();
        if self.vects.iter().any(|(_, v)| v.get_num_bits() != num_bits) {
            return Err(Error::Value {
                message: "the vectors have different sizes".to_owned(),
            });
        }
        let mut or_vect = first.clone();
        for (_, v) in &self.vects[1..] {
            or_vect = &or_vect | v;
        }
        self.or_vect = Some(or_vect);
        self.num_bits = i64::from(num_bits);
        self.need_reset = false;
        Ok(())
    }

    /// The number of vectors. Mirrors `VectCollection.NumChildren`.
    pub fn num_children(&self) -> usize {
        self.vects.len()
    }

    /// The keys and vectors, in the order the keys were first added.
    /// Mirrors `VectCollection.GetChildren`.
    pub fn get_children(&self) -> &[(u32, ExplicitBitVect)] {
        &self.vects
    }

    /// Whether any vector sets bit `idx`. Mirrors `VectCollection.GetBit`
    /// and `collection[idx]`.
    ///
    /// `Error::Value` with no vectors (Python's `AttributeError`) or as
    /// [`reset`](Self::reset); `Error::Index` for a bit past the end.
    pub fn get_bit(&mut self, idx: u32) -> Result<bool> {
        self.reset_if_needed()?;
        self.union()?.get_bit(idx)
    }

    /// The vectors' size: -1 before any vector is added, and, as in Python,
    /// the last size seen once all are removed. Mirrors
    /// `VectCollection.GetNumBits` and `len(collection)`.
    ///
    /// `Error::Value` as [`reset`](Self::reset).
    pub fn get_num_bits(&mut self) -> Result<i64> {
        self.reset_if_needed()?;
        Ok(self.num_bits)
    }

    /// The bits any vector sets. Mirrors `VectCollection.GetOnBits`.
    ///
    /// Errors as [`get_bit`](Self::get_bit).
    pub fn get_on_bits(&mut self) -> Result<Vec<u32>> {
        self.reset_if_needed()?;
        Ok(self.union()?.get_on_bits())
    }

    fn detach(&mut self, bit: u32, keep_matching: bool) -> Result<()> {
        // Python deletes as it goes, so a bad bit leaves the earlier
        // deletions done.
        let mut i = 0;
        while i < self.vects.len() {
            if self.vects[i].1.get_bit(bit)? == keep_matching {
                i += 1;
            } else {
                self.vects.remove(i);
                self.need_reset = true;
            }
        }
        Ok(())
    }

    /// Removes the vectors that do not set `bit`. Mirrors
    /// `VectCollection.DetachVectsNotMatchingBit`.
    ///
    /// `Error::Index` for a bit past a vector's end, as Python raises.
    pub fn detach_vects_not_matching_bit(&mut self, bit: u32) -> Result<()> {
        self.detach(bit, true)
    }

    /// Removes the vectors that set `bit`. Mirrors
    /// `VectCollection.DetachVectsMatchingBit`.
    ///
    /// `Error::Index` for a bit past a vector's end, as Python raises.
    pub fn detach_vects_matching_bit(&mut self, bit: u32) -> Result<()> {
        self.detach(bit, false)
    }

    /// Removes vectors with the same on bits as an earlier one. Mirrors
    /// `VectCollection.Uniquify`, whose `verbose` argument prints the
    /// counts before and after.
    pub fn uniquify(&mut self) {
        let obls: Vec<Vec<u32>> = self.vects.iter().map(|(_, v)| v.get_on_bits()).collect();
        let mut keep = vec![true; self.vects.len()];
        for i in 0..self.vects.len() {
            if keep[i] {
                for j in i + 1..self.vects.len() {
                    if keep[j] && obls[i] == obls[j] {
                        keep[j] = false;
                    }
                }
            }
        }
        let mut flags = keep.into_iter();
        self.vects.retain(|_| flags.next().unwrap_or(false));
    }
}
