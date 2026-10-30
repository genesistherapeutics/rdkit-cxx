// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.SimDivFilters.SimilarityPickers`, a Python module:
//! picking the fingerprints of a data set most similar to a set of probes.
//!
//! Python's pickers hand back the data set's fingerprints themselves; these
//! hand back their indices into the data set. Python's `GenericPicker` base
//! class, whose `MakePicks` only raises, has no counterpart; each picker
//! has its own `len` and `get`.

use std::fmt;

use crate::data_structs::top_n_container::TopNContainer;
use crate::data_structs::{
    ExplicitBitVect, fingerprint_similarity_with_metric, tanimoto_similarity,
};
use crate::py::bisect_right_by;
use crate::{Error, Result};

fn index_error(index: i64, len: usize) -> Error {
    Error::Index { index, len }
}

/// Python's list indexing, negative indices counting from the end.
fn pick_at<T: Copy>(picks: &[(T, f64)], which: i64) -> Result<(T, f64)> {
    let n = picks.len();
    let idx = if which < 0 {
        i64::try_from(n).unwrap_or(i64::MAX) + which
    } else {
        which
    };
    usize::try_from(idx)
        .ok()
        .and_then(|i| picks.get(i).copied())
        .ok_or(index_error(which, n))
}

/// Parameters for [`TopNOverallPicker::new_with_params`] and
/// [`SpreadPicker::new_with_params`], with Python's defaults (A2), except
/// that Python's default probes and data set, `None`, are empty here.
#[derive(Clone, Copy)]
pub struct PickerParams<'a> {
    /// How many fingerprints to pick.
    pub num_to_pick: usize,
    /// The fingerprints to compare the data set with.
    pub probe_fps: &'a [ExplicitBitVect],
    /// The fingerprints to pick from.
    pub data_set: &'a [ExplicitBitVect],
    /// The similarity, taken with `fingerprint_similarity_with_metric`,
    /// which folds fingerprints of different sizes. Any function or closure
    /// will do, such as `&|a, b| tversky_similarity(a, b, 0.3, 0.7)`.
    pub sim_metric: &'a dyn Fn(&ExplicitBitVect, &ExplicitBitVect) -> Result<f64>,
}

impl fmt::Debug for PickerParams<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PickerParams")
            .field("num_to_pick", &self.num_to_pick)
            .field("probe_fps", &self.probe_fps)
            .field("data_set", &self.data_set)
            .finish_non_exhaustive()
    }
}

impl Default for PickerParams<'_> {
    fn default() -> Self {
        PickerParams {
            num_to_pick: 10,
            probe_fps: &[],
            data_set: &[],
            sim_metric: &tanimoto_similarity::<ExplicitBitVect>,
        }
    }
}

/// Picks the fingerprints with the highest similarity to any probe, best
/// first. Mirrors `SimilarityPickers.TopNOverallPicker`.
#[derive(Clone, Debug)]
pub struct TopNOverallPicker<'a> {
    params: PickerParams<'a>,
    picks: Option<Vec<(Option<usize>, f64)>>,
}

impl<'a> TopNOverallPicker<'a> {
    /// A picker. Mirrors `TopNOverallPicker(numToPick, probeFps, dataSet,
    /// simMetric)`.
    pub fn new_with_params(params: &PickerParams<'a>) -> TopNOverallPicker<'a> {
        TopNOverallPicker {
            params: *params,
            picks: None,
        }
    }

    /// Makes the picks, unless they are made. Mirrors
    /// `TopNOverallPicker.MakePicks()`.
    ///
    /// Errors from the similarity metric pass through.
    pub fn make_picks(&mut self) -> Result<()> {
        self.make_picks_with_force(false)
    }

    /// Makes the picks, again if `force`. Mirrors
    /// `TopNOverallPicker.MakePicks(force)`.
    pub fn make_picks_with_force(&mut self, force: bool) -> Result<()> {
        if self.picks.is_some() && !force {
            return Ok(());
        }
        let p = &self.params;
        let size = i64::try_from(p.num_to_pick).unwrap_or(i64::MAX);
        let mut picks: TopNContainer<usize> = TopNContainer::new(size);
        for (i, fp) in p.data_set.iter().enumerate() {
            let mut best_score: f64 = -1.0;
            for probe in p.probe_fps {
                let score = fingerprint_similarity_with_metric(fp, probe, p.sim_metric)?;
                // Python's `max(score, bestScore)`.
                if best_score.partial_cmp(&score) != Some(std::cmp::Ordering::Greater) {
                    best_score = score;
                }
            }
            picks.insert_with_extra(best_score, i)?;
        }
        // Python reads every slot of the container, including those no
        // fingerprint filled, largest first.
        let mut res: Vec<(Option<usize>, f64)> = picks
            .get_pts()
            .iter()
            .zip(picks.get_extras())
            .map(|(&score, &idx)| (idx, score))
            .collect();
        res.reverse();
        self.picks = Some(res);
        Ok(())
    }

    /// The number of picks, making them if needed. Mirrors `len(picker)`.
    ///
    /// Errors as [`make_picks`](Self::make_picks).
    pub fn len(&mut self) -> Result<usize> {
        self.make_picks()?;
        Ok(self.picks.as_ref().map_or(0, Vec::len))
    }

    /// Whether there are no picks, making them if needed.
    pub fn is_empty(&mut self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// The pick at `which` (counting from the end when negative), as the
    /// fingerprint's index in the data set and its score, making the picks
    /// if needed. Mirrors `picker[which]`. With fewer fingerprints than
    /// picks, the last picks are `(None, -1e99)`, as in Python.
    ///
    /// `Error::Index` past the picks; otherwise as
    /// [`make_picks`](Self::make_picks).
    pub fn get(&mut self, which: i64) -> Result<(Option<usize>, f64)> {
        self.make_picks()?;
        pick_at(self.picks.as_deref().unwrap_or(&[]), which)
    }
}

/// Picks fingerprints probe by probe in turn, each probe's most similar
/// fingerprint not yet picked. Mirrors `SimilarityPickers.SpreadPicker`.
///
/// Python's `expectPickles` argument is unused, and its `onlyNames` picks
/// names read from a database (`rdkit.Dbase`); neither is taken.
#[derive(Clone, Debug)]
pub struct SpreadPicker<'a> {
    params: PickerParams<'a>,
    picks: Option<Vec<(usize, f64)>>,
}

impl<'a> SpreadPicker<'a> {
    /// A picker. Mirrors `SpreadPicker(numToPick, probeFps, dataSet,
    /// simMetric)`.
    pub fn new_with_params(params: &PickerParams<'a>) -> SpreadPicker<'a> {
        SpreadPicker {
            params: *params,
            picks: None,
        }
    }

    /// Makes the picks, unless they are made. Mirrors
    /// `SpreadPicker.MakePicks()`, which also prints its progress every
    /// 1000 fingerprints unless `silent`.
    ///
    /// `Error::Value` with no probes (Python's `ZeroDivisionError`), and
    /// `Error::Index` when the probes run out of fingerprints to pick or
    /// there are more probes than picks (Python's `IndexError`s); errors
    /// from the similarity metric pass through.
    pub fn make_picks(&mut self) -> Result<()> {
        self.make_picks_with_force(false)
    }

    /// Makes the picks, again if `force`. Mirrors
    /// `SpreadPicker.MakePicks(force)`.
    pub fn make_picks_with_force(&mut self, force: bool) -> Result<()> {
        if self.picks.is_some() && !force {
            return Ok(());
        }
        let p = &self.params;
        let mut scores: Vec<Vec<(f64, usize)>> = vec![Vec::new(); p.probe_fps.len()];
        for (j, fp) in p.data_set.iter().enumerate() {
            for i in 0..p.probe_fps.len() {
                // Python truncates the list of rows, not the row, when a row
                // reaches `numToPick` scores, and then fails on a missing
                // row.
                let n_rows = scores.len();
                let row = scores.get_mut(i).ok_or(index_error(i as i64, n_rows))?;
                let score = fingerprint_similarity_with_metric(&p.probe_fps[i], fp, p.sim_metric)?;
                // Python's `bisect.insort`, comparing (score, index) tuples.
                let lt = |&(s, k): &(f64, usize)| if score != s { score < s } else { j < k };
                let at = bisect_right_by(row, lt);
                row.insert(at, (score, j));
                if row.len() >= p.num_to_pick {
                    scores.truncate(p.num_to_pick);
                }
            }
        }
        let mut picks = Vec::new();
        let mut taken = vec![false; p.data_set.len()];
        while picks.len() < p.num_to_pick {
            if scores.is_empty() {
                return Err(Error::Value {
                    message: "integer modulo by zero: there are no probes".to_owned(),
                });
            }
            let n_rows = scores.len();
            let row = &mut scores[picks.len() % n_rows];
            let mut next = row.pop().ok_or(index_error(-1, 0))?;
            while taken[next.1] && !row.is_empty() {
                next = row.pop().ok_or(index_error(-1, 0))?;
            }
            let (score, idx) = next;
            if !taken[idx] {
                picks.push((idx, score));
                taken[idx] = true;
            }
        }
        self.picks = Some(picks);
        Ok(())
    }

    /// The number of picks, making them if needed. Mirrors `len(picker)`.
    ///
    /// Errors as [`make_picks`](Self::make_picks).
    pub fn len(&mut self) -> Result<usize> {
        self.make_picks()?;
        Ok(self.picks.as_ref().map_or(0, Vec::len))
    }

    /// Whether there are no picks, making them if needed.
    pub fn is_empty(&mut self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// The pick at `which` (counting from the end when negative), as the
    /// fingerprint's index in the data set and its score, making the picks
    /// if needed. Mirrors `picker[which]`.
    ///
    /// `Error::Index` past the picks; otherwise as
    /// [`make_picks`](Self::make_picks).
    pub fn get(&mut self, which: i64) -> Result<(usize, f64)> {
        self.make_picks()?;
        pick_at(self.picks.as_deref().unwrap_or(&[]), which)
    }
}
