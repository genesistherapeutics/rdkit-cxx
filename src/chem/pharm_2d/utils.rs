// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm2D.Utils`, a Python module: the
//! combinatorics of 2D pharmacophores.
//!
//! An N-point pharmacophore has the distances between its points that
//! [`n_point_dist_dict`] lists, in that order; a scaffold is one distance
//! bin index per distance. Python's arguments "used in recursion" are not
//! bound.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::{Error, Result};

/// The point pairs whose distances describe an N-point pharmacophore, for
/// 2 to 10 points.
const N_POINT_DISTS: [&[(usize, usize)]; 9] = [
    &[(0, 1)],
    &[(0, 1), (0, 2), (1, 2)],
    &[(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)],
    &[(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)],
    &[
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (0, 5),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
    ],
    &[
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (0, 5),
        (0, 6),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
    ],
    &[
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (0, 5),
        (0, 6),
        (0, 7),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
    ],
    &[
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (0, 5),
        (0, 6),
        (0, 7),
        (0, 8),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
    ],
    &[
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (0, 5),
        (0, 6),
        (0, 7),
        (0, 8),
        (0, 9),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
    ],
];

/// The point pairs of an `n_pts`-point pharmacophore's distances, or
/// `Error::Key` outside 2 to 10 points, as Python's dictionary lookup.
pub(crate) fn point_dists(n_pts: usize) -> Result<&'static [(usize, usize)]> {
    n_pts
        .checked_sub(2)
        .and_then(|i| N_POINT_DISTS.get(i))
        .copied()
        .ok_or_else(|| Error::Key {
            key: n_pts.to_string(),
        })
}

/// The point pairs whose distances describe an N-point pharmacophore, by
/// N from 2 to 10. Mirrors `Utils.nPointDistDict`.
pub fn n_point_dist_dict() -> BTreeMap<usize, Vec<(usize, usize)>> {
    (2..=10)
        .map(|n| (n, point_dists(n).expect("2 to 10 points").to_vec()))
        .collect()
}

/// The number of points of a pharmacophore with N distances, by N. Mirrors
/// `Utils.nDistPointDict`.
pub fn n_dist_point_dict() -> BTreeMap<usize, usize> {
    (2..=10).map(|n| (2 * n - 3, n)).collect()
}

/// The distance indices of the triangles that make up an `n_pts`-point
/// pharmacophore; none below 3 points. Mirrors `Utils.GetTriangles`.
pub fn get_triangles(n_pts: usize) -> Vec<(usize, usize, usize)> {
    if n_pts < 3 {
        return Vec::new();
    }
    (0..n_pts - 2).map(|i| (i, i + 1, i + n_pts - 1)).collect()
}

/// Whether three distance bins can make a triangle: each bin's upper bound
/// plus another's is at least the third's lower bound. Mirrors
/// `Utils.BinsTriangleInequality`.
pub fn bins_triangle_inequality(d1: (u32, u32), d2: (u32, u32), d3: (u32, u32)) -> bool {
    let sum = |a: u32, b: u32| u64::from(a) + u64::from(b);
    !(sum(d1.1, d2.1) < u64::from(d3.0)
        || sum(d2.1, d3.1) < u64::from(d1.0)
        || sum(d3.1, d1.1) < u64::from(d2.0))
}

/// Whether every triangle of a scaffold (one bin index per distance) can
/// satisfy the triangle inequality. `Error::Key` for a scaffold whose
/// length no pharmacophore has, and `Error::Index` for a bin index past
/// `bins`. Mirrors `Utils.ScaffoldPasses`.
pub fn scaffold_passes(combo: &[usize], bins: &[(u32, u32)]) -> Result<bool> {
    let n_pts = *n_dist_point_dict()
        .get(&combo.len())
        .ok_or_else(|| Error::Key {
            key: combo.len().to_string(),
        })?;
    let bin = |i: usize| -> Result<(u32, u32)> {
        let b = combo[i];
        bins.get(b).copied().ok_or(Error::index(b, bins.len()))
    };
    for (a, b, c) in get_triangles(n_pts) {
        if !bins_triangle_inequality(bin(a)?, bin(b)?, bin(c)?) {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The binomial coefficient, as Python's `math.comb`.
fn comb(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut res: u128 = 1;
    for i in 0..k {
        res = res * u128::from(n - i) / u128::from(i + 1);
    }
    u64::try_from(res).expect("the number of combinations fits in 64 bits")
}

/// The number of ways to put `n_items` into `n_slots`, where order does
/// not matter and an item may repeat: (N + S − 1)! / ((N − 1)! S!).
///
/// # Panics
///
/// When the number does not fit in 64 bits (Python's integers grow).
///
/// Mirrors `Utils.NumCombinations`.
pub fn num_combinations(n_items: usize, n_slots: usize) -> u64 {
    if n_items + n_slots == 0 {
        // Python's comb(-1, 0) is a ValueError; comb(n, 0) is 1.
        return 1;
    }
    comb((n_items + n_slots - 1) as u64, n_slots as u64)
}

/// Where a combination of indices (sorted, as
/// [`get_index_combinations`] makes them) falls among the combinations of
/// `n_items` in `n_slots`. `Error::Index` when `vs` is shorter than
/// `n_slots`. Mirrors `Utils.CountUpTo`.
pub fn count_up_to(n_items: usize, n_slots: usize, vs: &[usize]) -> Result<u64> {
    if vs.len() < n_slots {
        return Err(Error::index(vs.len(), vs.len()));
    }
    let mut accum: i128 = 0;
    let mut start_at = 0usize;
    for (idx, &v) in vs[..n_slots].iter().enumerate() {
        if idx == n_slots - 1 {
            accum += v as i128 - start_at as i128;
        } else {
            for i in start_at..v {
                let n_vals_over = n_items.saturating_sub(i);
                accum += i128::from(num_combinations(n_vals_over, n_slots - idx - 1));
            }
            start_at = v;
        }
    }
    u64::try_from(accum).map_err(|_| Error::Value {
        message: format!("{vs:?} is not a sorted combination of {n_items} items"),
    })
}

/// Every way to put `n_items` into `n_slots` with each slot's item no
/// smaller than the last's, in order. Mirrors `Utils.GetIndexCombinations`.
pub fn get_index_combinations(n_items: usize, n_slots: usize) -> Vec<Vec<usize>> {
    fn recurse(n_items: usize, n_slots: usize, slot: usize, last: usize) -> Vec<Vec<usize>> {
        if slot >= n_slots {
            Vec::new()
        } else if slot == n_slots - 1 {
            (last..n_items).map(|x| vec![x]).collect()
        } else {
            let mut res = Vec::new();
            for x in last..n_items {
                for entry in recurse(n_items, n_slots, slot + 1, x) {
                    let mut combo = vec![x];
                    combo.extend(entry);
                    res.push(combo);
                }
            }
            res
        }
    }
    recurse(n_items, n_slots, 0, 0)
}

/// Every combination of one element from each of `choices`, in order,
/// leaving out those with an element twice. Mirrors
/// `Utils.GetAllCombinations`.
pub fn get_all_combinations<T: Clone + PartialEq>(choices: &[Vec<T>]) -> Vec<Vec<T>> {
    get_all_combinations_with_no_dups(choices, true)
}

/// Keeps combinations with repeated elements unless `no_dups`. Mirrors
/// `Utils.GetAllCombinations(choices, noDups)`.
pub fn get_all_combinations_with_no_dups<T: Clone + PartialEq>(
    choices: &[Vec<T>],
    no_dups: bool,
) -> Vec<Vec<T>> {
    let Some((first, rest)) = choices.split_first() else {
        return Vec::new();
    };
    if rest.is_empty() {
        return first.iter().map(|x| vec![x.clone()]).collect();
    }
    let tmp = get_all_combinations_with_no_dups(rest, no_dups);
    let mut res = Vec::new();
    for thing in first {
        for other in &tmp {
            if !no_dups || !other.contains(thing) {
                let mut combo = vec![thing.clone()];
                combo.extend(other.iter().cloned());
                res.push(combo);
            }
        }
    }
    res
}

/// The distinct combinations of one element from each of `choices`
/// without repeated elements, each as its elements paired with their
/// choice's class and sorted, in sorted order. Choices and classes of
/// different lengths are `Error::Value`. Mirrors
/// `Utils.GetUniqueCombinations`.
pub fn get_unique_combinations<T: Clone + Ord>(
    choices: &[Vec<T>],
    classes: &[usize],
) -> Result<Vec<Vec<(usize, T)>>> {
    if choices.len() != classes.len() {
        return Err(Error::Value {
            message: "choices and classes must have the same length".to_owned(),
        });
    }
    let mut combos = BTreeSet::new();
    if choices.iter().any(Vec::is_empty) {
        return Ok(Vec::new());
    }
    // Walk the product of the choices like an odometer, as itertools.product.
    let mut positions = vec![0usize; choices.len()];
    loop {
        let choice: Vec<&T> = positions.iter().zip(choices).map(|(&p, c)| &c[p]).collect();
        let distinct = choice.iter().collect::<BTreeSet<_>>().len() == choice.len();
        if distinct {
            let mut combo: Vec<(usize, T)> = classes
                .iter()
                .zip(&choice)
                .map(|(&cls, &ch)| (cls, ch.clone()))
                .collect();
            combo.sort();
            combos.insert(combo);
        }
        let mut i = choices.len();
        loop {
            if i == 0 {
                return Ok(combos.into_iter().collect());
            }
            i -= 1;
            positions[i] += 1;
            if positions[i] < choices[i].len() {
                break;
            }
            positions[i] = 0;
        }
    }
}

/// One of each set of combinations with the same elements in any order:
/// the last of them, where the first appeared. Mirrors
/// `Utils.UniquifyCombinations`.
pub fn uniquify_combinations<T: Clone + Ord + std::hash::Hash>(combos: &[Vec<T>]) -> Vec<Vec<T>> {
    let mut res: Vec<Vec<T>> = Vec::new();
    let mut seen: HashMap<Vec<T>, usize> = HashMap::new();
    for combo in combos {
        let mut key = combo.clone();
        key.sort();
        match seen.get(&key) {
            Some(&i) => res[i] = combo.clone(),
            None => {
                seen.insert(key, res.len());
                res.push(combo.clone());
            }
        }
    }
    res
}

/// Every scaffold of `n_pts` points over `bins` whose triangles satisfy
/// the triangle inequality, in order; none below 2 points, and
/// `Error::Key` above 10. Mirrors `Utils.GetPossibleScaffolds`.
pub fn get_possible_scaffolds(n_pts: usize, bins: &[(u32, u32)]) -> Result<Vec<Vec<usize>>> {
    get_possible_scaffolds_with_use_triangle_inequality(n_pts, bins, true)
}

/// Keeps every scaffold unless `use_triangle_inequality`. Mirrors
/// `Utils.GetPossibleScaffolds(nPts, bins, useTriangleInequality)`.
pub fn get_possible_scaffolds_with_use_triangle_inequality(
    n_pts: usize,
    bins: &[(u32, u32)],
    use_triangle_inequality: bool,
) -> Result<Vec<Vec<usize>>> {
    if n_pts < 2 {
        return Ok(Vec::new());
    }
    if n_pts == 2 {
        return Ok((0..bins.len()).map(|x| vec![x]).collect());
    }
    let n_dists = point_dists(n_pts)?.len();
    let all: Vec<usize> = (0..bins.len()).collect();
    let combos = get_all_combinations_with_no_dups(&vec![all; n_dists], false);
    let mut res = Vec::new();
    for combo in combos {
        if !use_triangle_inequality || scaffold_passes(&combo, bins)? {
            res.push(combo);
        }
    }
    Ok(res)
}

/// A triangle's feature indices and distances in canonical order, so the
/// same triangle gives the same bit however its features are numbered.
/// Lists of other lengths than 3 are `Error::Value`. Mirrors
/// `Utils.OrderTriangle`.
pub fn order_triangle(feat_indices: &[usize], dists: &[u32]) -> Result<(Vec<usize>, Vec<u32>)> {
    if feat_indices.len() != 3 {
        return Err(Error::Value {
            message: "bad indices".to_owned(),
        });
    }
    if dists.len() != 3 {
        return Err(Error::Value {
            message: "bad dists".to_owned(),
        });
    }
    let distinct = feat_indices.iter().collect::<BTreeSet<_>>().len();
    if distinct == 3 {
        return Ok((feat_indices.to_vec(), dists.to_vec()));
    }
    let d = |i: usize| u64::from(dists[i]);
    let d_sums = [d(0) + d(1), d(0) + d(2), d(1) + d(2)];
    let m_d = d_sums.iter().copied().max().expect("three sums");
    let (ireorder, dreorder) = if distinct == 1 {
        if d_sums[0] == m_d {
            if dists[0] > dists[1] {
                ([0, 1, 2], [0, 1, 2])
            } else {
                ([0, 2, 1], [1, 0, 2])
            }
        } else if d_sums[1] == m_d {
            if dists[0] > dists[2] {
                ([1, 0, 2], [0, 2, 1])
            } else {
                ([1, 2, 0], [2, 0, 1])
            }
        } else if dists[1] > dists[2] {
            ([2, 0, 1], [1, 2, 0])
        } else {
            ([2, 1, 0], [2, 1, 0])
        }
    } else if feat_indices[0] == feat_indices[1] {
        if dists[1] > dists[2] {
            ([0, 1, 2], [0, 1, 2])
        } else {
            ([1, 0, 2], [0, 2, 1])
        }
    } else if feat_indices[0] == feat_indices[2] {
        if dists[0] > dists[2] {
            ([0, 1, 2], [0, 1, 2])
        } else {
            ([2, 1, 0], [2, 1, 0])
        }
    } else if dists[0] > dists[1] {
        ([0, 1, 2], [0, 1, 2])
    } else {
        ([0, 2, 1], [1, 0, 2])
    };
    Ok((
        ireorder.iter().map(|&i| feat_indices[i]).collect(),
        dreorder.iter().map(|&i| dists[i]).collect(),
    ))
}
