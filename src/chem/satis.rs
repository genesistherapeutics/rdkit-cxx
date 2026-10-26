// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.SATIS`, a Python module: SATIS atom types
//! (_J. Chem. Inf. Comput. Sci._ 39, 751-757, 1999).
//!
//! Python parses its patterns once, on import; each function here parses a
//! new molecule, since a `Mol` cannot be shared between threads
//! (CONVENTIONS.md S1).

use std::collections::BTreeSet;

use crate::chem::{Mol, mol_from_smarts};
use crate::{Error, Result};

fn pattern(smarts: &str) -> Mol {
    mol_from_smarts(smarts).expect("the SATIS SMARTS parse")
}

/// An aldehyde's carbonyl. Mirrors `SATIS.aldehydePatt`.
pub fn aldehyde_patt() -> Mol {
    pattern("[CD2]=[OD1]")
}

/// A ketone's carbonyl. Mirrors `SATIS.ketonePatt`.
pub fn ketone_patt() -> Mol {
    pattern("[CD3]=[OD1]")
}

/// An amide. Mirrors `SATIS.amidePatt`.
pub fn amide_patt() -> Mol {
    pattern("[CD3](=[OD1])-[#7]")
}

/// An ester. Mirrors `SATIS.esterPatt`.
pub fn ester_patt() -> Mol {
    pattern("C(=[OD1])-O-[#6]")
}

/// A carboxylate. Mirrors `SATIS.carboxylatePatt`.
pub fn carboxylate_patt() -> Mol {
    pattern("C(=[OD1])-[OX1]")
}

/// A carboxyl group. Mirrors `SATIS.carboxylPatt`.
pub fn carboxyl_patt() -> Mol {
    pattern("C(=[OD1])-[OX2]")
}

/// The special cases, in the order they are tried, with the code each gives
/// the last place of the SATIS code of a carbon or oxygen that has fewer
/// neighbors than it has places. Mirrors `SATIS.specialCases`.
pub fn special_cases() -> Vec<(Mol, u32)> {
    vec![
        (carboxylate_patt(), 97),
        (ester_patt(), 96),
        (carboxyl_patt(), 98),
        (amide_patt(), 95),
        (ketone_patt(), 94),
        (aldehyde_patt(), 93),
    ]
}

/// Each atom's SATIS code, with four neighbors: two-digit numbers for the
/// atom's atomic number and its neighbors', hydrogens included, largest
/// last. Mirrors `SATIS.SATISTypes(mol)`.
pub fn satis_types(mol: &Mol) -> Result<Vec<String>> {
    satis_types_with_neighbors_to_include(mol, 4)
}

/// As [`satis_types`], with `neighbors_to_include` neighbors. Zero, for a
/// molecule with a bond or a hydrogen, is `Error::Index`, as Python's
/// `IndexError`. Mirrors `SATIS.SATISTypes(mol, neighborsToInclude)`.
pub fn satis_types_with_neighbors_to_include(
    mol: &Mol,
    neighbors_to_include: u32,
) -> Result<Vec<String>> {
    let n = neighbors_to_include as usize;
    let atomic_nums: Vec<i32> = mol.get_atoms().map(|a| a.get_atomic_num()).collect();
    let mut special_case_matches = Vec::new();
    for (patt, special_case_idx) in special_cases() {
        let matches = mol.get_substruct_matches(&patt)?;
        if !matches.is_empty() {
            let matches: BTreeSet<u32> = matches.into_iter().flatten().collect();
            special_case_matches.push((special_case_idx, matches));
        }
    }
    let mut codes = Vec::with_capacity(atomic_nums.len());
    for atom in mol.get_atoms() {
        let mut code: Vec<i64> = vec![99; n + 1];
        code[0] = i64::from(atom.get_atomic_num()).min(99);
        let mut other_nums: Vec<i64> = atom
            .get_neighbors()
            .iter()
            .map(|x| i64::from(atomic_nums[x.get_idx() as usize]))
            .collect();
        other_nums.extend(std::iter::repeat_n(1, atom.get_total_num_hs()? as usize));
        other_nums.sort();
        // Python keeps the last n, and `[-0:]` keeps them all.
        if n > 0 && other_nums.len() > n {
            other_nums.drain(..other_nums.len() - n);
        }
        for (j, other_num) in other_nums.iter().enumerate() {
            let len = code.len();
            *code.get_mut(j + 1).ok_or(Error::index(j + 1, len))? = (*other_num).min(99);
        }
        if other_nums.len() < n && (code[0] == 6 || code[0] == 8) {
            let atom_idx = atom.get_idx();
            for (special_case_idx, matches) in &special_case_matches {
                if matches.contains(&atom_idx) {
                    code[n] = i64::from(*special_case_idx);
                    break;
                }
            }
        }
        codes.push(code.iter().map(|x| format!("{x:02}")).collect());
    }
    Ok(codes)
}
