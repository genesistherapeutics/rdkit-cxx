// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.AtomPairs.Utils`, a Python module: what an atom
//! code means, and similarities between sorted lists of bit ids.
//!
//! Python's `GetAtomCode` is `rdMolDescriptors.GetAtomPairAtomCode`,
//! [`get_atom_pair_atom_code`](crate::chem::rd_mol_descriptors::get_atom_pair_atom_code).

use crate::chem::get_periodic_table;
use crate::chem::rd_mol_descriptors::AtomPairsParameters;
use crate::{Error, Result};

/// An atom code's element symbol, number of heavy neighbors (less the
/// branch subtraction), and number of pi electrons, as Python's tuple.
pub type AtomCodeExplanation = (String, u32, u32);

/// An [`AtomCodeExplanation`] with the atom's CIP label: Python's 4-tuple.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct ChiralAtomCodeExplanation {
    /// The element symbol.
    pub symbol: String,
    /// The number of heavy neighbors, less the branch subtraction.
    pub n_branch: u32,
    /// The number of pi electrons.
    pub n_pi: u32,
    /// The CIP label: `""`, `"R"`, or `"S"`.
    pub chirality: String,
}

fn mask(bits: u32) -> u64 {
    (1u64 << bits) - 1
}

fn explain(code: u64) -> (AtomCodeExplanation, u64) {
    let n_branch = (code & mask(AtomPairsParameters::num_branch_bits())) as u32;
    let code = code >> AtomPairsParameters::num_branch_bits();
    let n_pi = (code & mask(AtomPairsParameters::num_pi_bits())) as u32;
    let code = code >> AtomPairsParameters::num_pi_bits();
    let type_idx = (code & mask(AtomPairsParameters::num_type_bits())) as usize;
    let symbol = match AtomPairsParameters::atom_types().get(type_idx) {
        Some(&atomic_num) => get_periodic_table()
            .get_element_symbol(atomic_num)
            .expect("the atom-pair types are elements"),
        None => "X".to_owned(),
    };
    let rest = code >> AtomPairsParameters::num_type_bits();
    ((symbol, n_branch, n_pi), rest)
}

/// What an atom code means: its element (`"X"` past the atom-pair types),
/// branches, and pi electrons. A code made with chirality explains the
/// same, without its CIP label. Python's `branchSubtract` argument is
/// unused, so it is not bound. Mirrors `Utils.ExplainAtomCode`.
pub fn explain_atom_code(code: u64) -> AtomCodeExplanation {
    explain(code).0
}

/// With the CIP label of a code made with `includeChirality`; a label code
/// of 3 is `Error::Key`, as Python's dictionary lookup is. Mirrors
/// `Utils.ExplainAtomCode(code, includeChirality=True)` (A4).
pub fn explain_atom_code_include_chirality(code: u64) -> Result<ChiralAtomCodeExplanation> {
    let ((symbol, n_branch, n_pi), rest) = explain(code);
    let chirality = match rest & mask(AtomPairsParameters::num_chiral_bits()) {
        0 => "",
        1 => "R",
        2 => "S",
        other => {
            return Err(Error::Key {
                key: other.to_string(),
            });
        }
    };
    Ok(ChiralAtomCodeExplanation {
        symbol,
        n_branch,
        n_pi,
        chirality: chirality.to_owned(),
    })
}

/// The number of bit ids two sorted lists share, counting a duplicate once
/// for each time both lists hold it. Mirrors `Utils.BitsInCommon`.
pub fn bits_in_common<T: PartialOrd>(v1: &[T], v2: &[T]) -> usize {
    let mut res = 0;
    let mut v2_pos = 0;
    for val in v1 {
        while v2_pos < v2.len() && v2[v2_pos] < *val {
            v2_pos += 1;
        }
        if v2_pos >= v2.len() {
            break;
        }
        if v2[v2_pos] == *val {
            res += 1;
            v2_pos += 1;
        }
    }
    res
}

/// Dice similarity between two sorted lists of bit ids; 0 for two empty
/// lists. Mirrors `Utils.DiceSimilarity`.
pub fn dice_similarity<T: PartialOrd>(v1: &[T], v2: &[T]) -> f64 {
    dice_similarity_with_bounds(v1, v2, None)
}

/// Returns 0 without comparing the lists when the shorter one's share of
/// both is below `bounds` (Python skips a bound of `None` or 0). Mirrors
/// `Utils.DiceSimilarity(v1, v2, bounds)`.
pub fn dice_similarity_with_bounds<T: PartialOrd>(v1: &[T], v2: &[T], bounds: Option<f64>) -> f64 {
    let denom = (v1.len() + v2.len()) as f64;
    if denom == 0.0 {
        return 0.0;
    }
    let numer = match bounds {
        Some(b) if b != 0.0 && (v1.len().min(v2.len()) as f64 / denom) < b => 0.0,
        _ => 2.0 * bits_in_common(v1, v2) as f64,
    };
    numer / denom
}

/// The dot product of two sorted lists of bit ids as count vectors, a
/// shared id adding the square of the smaller count. Mirrors `Utils.Dot`.
pub fn dot<T: PartialOrd>(v1: &[T], v2: &[T]) -> usize {
    let mut res = 0;
    let (mut i, mut j) = (0, 0);
    while i < v1.len() {
        let v1_val = &v1[i];
        let mut v1_count = 1;
        i += 1;
        while i < v1.len() && v1[i] == *v1_val {
            v1_count += 1;
            i += 1;
        }
        while j < v2.len() && v2[j] < *v1_val {
            j += 1;
        }
        if j < v2.len() && v2[j] == *v1_val {
            let mut v2_count = 1;
            j += 1;
            while j < v2.len() && v2[j] == *v1_val {
                v2_count += 1;
                j += 1;
            }
            let common = usize::min(v1_count, v2_count);
            res += common * common;
        } else if j >= v2.len() {
            break;
        }
    }
    res
}

/// Cosine similarity between two sorted lists of bit ids as count vectors;
/// 0 when either is empty. Mirrors `Utils.CosineSimilarity`.
pub fn cosine_similarity<T: PartialOrd>(v1: &[T], v2: &[T]) -> f64 {
    let denom = ((dot(v1, v1) * dot(v2, v2)) as f64).sqrt();
    if denom == 0.0 {
        0.0
    } else {
        dot(v1, v2) as f64 / denom
    }
}
