// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.DataStructs.TopNContainer`, a Python module: a sorted
//! container of the largest values inserted, each with an extra item.

use crate::py::bisect_right;
use crate::{Error, Result};

/// Keeps the `size` largest values inserted, smallest first, each with an
/// optional extra item; with a negative size, keeps every value. Mirrors
/// `TopNContainer.TopNContainer`.
#[derive(Clone, Debug, PartialEq)]
pub struct TopNContainer<E = ()> {
    size: i64,
    best: Vec<f64>,
    extras: Vec<Option<E>>,
}

fn index_error(index: i64, len: usize) -> Error {
    Error::Index { index, len }
}

impl<E> TopNContainer<E> {
    /// A container for the `size` largest values, which start as `-1e99`.
    /// Mirrors `TopNContainer(size)`.
    pub fn new(size: i64) -> TopNContainer<E> {
        TopNContainer::new_with_most_neg(size, -1e99)
    }

    /// As [`new`](Self::new), the values starting as `most_neg`. Mirrors
    /// `TopNContainer(size, mostNeg)`.
    pub fn new_with_most_neg(size: i64, most_neg: f64) -> TopNContainer<E> {
        let n = usize::try_from(size).unwrap_or(0);
        TopNContainer {
            size,
            best: vec![most_neg; n],
            extras: std::iter::repeat_with(|| None).take(n).collect(),
        }
    }

    /// Python's `bisect.bisect`, whose probes give Python's index even once
    /// a NaN leaves `best` unsorted.
    fn bisect(&self, val: f64) -> usize {
        bisect_right(&self.best, val)
    }

    /// Inserts `val`, with no extra item, if it is larger than the smallest
    /// value kept. Mirrors `TopNContainer.Insert(val)`.
    ///
    /// `Error::Index` for a container of size zero, as Python raises.
    pub fn insert(&mut self, val: f64) -> Result<()> {
        self.insert_option(val, None)
    }

    /// As [`insert`](Self::insert), with the extra item `extra`. Mirrors
    /// `TopNContainer.Insert(val, extra)`.
    pub fn insert_with_extra(&mut self, val: f64, extra: E) -> Result<()> {
        self.insert_option(val, Some(extra))
    }

    fn insert_option(&mut self, val: f64, extra: Option<E>) -> Result<()> {
        if self.size >= 0 {
            let smallest = *self.best.first().ok_or(index_error(0, 0))?;
            if val > smallest {
                let idx = self.bisect(val);
                self.best.insert(idx, val);
                self.extras.insert(idx, extra);
                self.best.remove(0);
                self.extras.remove(0);
            }
        } else {
            let idx = self.bisect(val);
            self.best.insert(idx, val);
            self.extras.insert(idx, extra);
        }
        Ok(())
    }

    /// The values, smallest first (unless [`reverse`](Self::reverse)d).
    /// Mirrors `TopNContainer.GetPts`.
    pub fn get_pts(&self) -> &[f64] {
        &self.best
    }

    /// The extra items, in the values' order. Mirrors
    /// `TopNContainer.GetExtras`.
    pub fn get_extras(&self) -> &[Option<E>] {
        &self.extras
    }

    /// The size, or the number of values for a container that keeps every
    /// value. Mirrors `len(container)`.
    pub fn len(&self) -> usize {
        usize::try_from(self.size).unwrap_or(self.best.len())
    }

    /// Whether [`len`](Self::len) is zero.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The value and extra item at `which`, counting from the end when
    /// negative. Mirrors `container[which]`.
    ///
    /// `Error::Index` outside the container.
    pub fn get(&self, which: i64) -> Result<(f64, Option<&E>)> {
        let n = self.best.len();
        let idx = if which < 0 {
            i64::try_from(n).unwrap_or(i64::MAX) + which
        } else {
            which
        };
        let idx = usize::try_from(idx)
            .ok()
            .filter(|&i| i < n)
            .ok_or(index_error(which, n))?;
        Ok((self.best[idx], self.extras[idx].as_ref()))
    }

    /// Reverses the values and extra items, largest first. Mirrors
    /// `TopNContainer.reverse`.
    pub fn reverse(&mut self) {
        self.best.reverse();
        self.extras.reverse();
    }
}
