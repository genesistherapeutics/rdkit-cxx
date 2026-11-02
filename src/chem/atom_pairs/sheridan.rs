// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.AtomPairs.Sheridan`, a Python module:
//! physicochemical-property fingerprints (Kearsley et al., J. Chem. Inf.
//! Model. 36, 118-127 (1996)), atom pairs and torsions over Patty's atom
//! types.
//!
//! Python reads the types from `RDConfig.RDDataDir`'s
//! `SmartsLib/patty_rules.txt`; the crate embeds RDKit's copy of that file
//! (`wrapper/third_party/rdkit`).

use std::collections::BTreeMap;

use crate::chem::rd_fingerprint_generator::{
    GetFingerprintParams, get_atom_pair_generator, get_topological_torsion_generator,
};
use crate::chem::rd_mol_descriptors::AtomPairsParameters;
use crate::chem::{Mol, SubstructMatchParameters, mol_from_smarts};
use crate::data_structs::ULongSparseIntVect;
use crate::{Error, Result};

const PATTY_RULES: &str =
    include_str!("../../../wrapper/third_party/rdkit/Data/SmartsLib/patty_rules.txt");

/// The bits for a pair's distance. Mirrors `Sheridan.numPathBits`.
pub fn num_path_bits() -> u32 {
    AtomPairsParameters::num_path_bits()
}

/// The bits in a pair's code. Mirrors `Sheridan.numFpBits`.
pub fn num_fp_bits() -> u32 {
    num_path_bits() + 2 * AtomPairsParameters::code_size()
}

/// The number of possible pair codes. Mirrors `Sheridan.fpLen`.
pub fn fp_len() -> u32 {
    1 << num_fp_bits()
}

/// Each Patty type's atom invariant. Mirrors `Sheridan.typMap`.
pub fn typ_map() -> BTreeMap<&'static str, u32> {
    BTreeMap::from([
        ("CAT", 1),
        ("ANI", 2),
        ("POL", 3),
        ("DON", 4),
        ("ACC", 5),
        ("HYD", 6),
        ("OTH", 7),
    ])
}

/// Patty's rules, as Python's `_readPattyDefs` reads them: a pattern and
/// its type on each line, skipping comments and patterns that do not parse.
fn patty_defs() -> Vec<(Mol, String)> {
    let mut matchers = Vec::new();
    for line in PATTY_RULES.lines() {
        let line = line.trim();
        let line = line.split("# ").next().unwrap_or_default().trim();
        if line.is_empty() {
            continue;
        }
        // Python splits on runs of spaces, as `re.split('[ ]+', x)`.
        let fields: Vec<&str> = line.split(' ').filter(|f| !f.is_empty()).collect();
        if fields.len() > 1
            && let Ok(mol) = mol_from_smarts(fields[0])
        {
            matchers.push((mol, fields[1].to_owned()));
        }
    }
    matchers
}

/// Each atom's Patty type (`CAT`, `ANI`, `POL`, `DON`, `ACC`, `HYD`, or
/// `OTH`), from the last rule that matches it first. Mirrors
/// `Sheridan.AssignPattyTypes`.
pub fn assign_patty_types(mol: &Mol) -> Result<Vec<String>> {
    assign_patty_types_with_defns(mol, &patty_defs())
}

/// With other rules: each pattern's matches give their first atom its
/// type, later rules overriding earlier ones; an atom no rule matches has
/// the empty type. Mirrors `Sheridan.AssignPattyTypes(mol, defns)`.
pub fn assign_patty_types_with_defns(mol: &Mol, defns: &[(Mol, String)]) -> Result<Vec<String>> {
    let mut res = vec![String::new(); mol.get_num_atoms() as usize];
    let params = SubstructMatchParameters {
        uniquify: false,
        ..Default::default()
    };
    for (matcher, name) in defns {
        for m in mol.get_substruct_matches_with_params(matcher, &params)? {
            if let Some(&first) = m.first() {
                res[first as usize] = name.clone();
            }
        }
    }
    Ok(res)
}

fn patty_invariants(mol: &Mol) -> Result<Vec<u32>> {
    let typ_map = typ_map();
    assign_patty_types(mol)?
        .into_iter()
        .map(|t| {
            typ_map
                .get(t.as_str())
                .copied()
                .ok_or(Error::Key { key: t })
        })
        .collect()
}

fn sparse_counts_with_invariants(
    mol: &Mol,
    invariants: &[u32],
    torsions: bool,
) -> Result<ULongSparseIntVect> {
    let generator = if torsions {
        get_topological_torsion_generator()
    } else {
        get_atom_pair_generator()
    };
    let params = GetFingerprintParams {
        custom_atom_invariants: Some(invariants.to_vec()),
        ..Default::default()
    };
    generator.get_sparse_count_fingerprint_with_params(mol, &params)
}

/// Sheridan's binding-property pairs: the unfolded atom-pair count
/// fingerprint over the atoms' Patty types. Mirrors
/// `Sheridan.GetBPFingerprint`.
pub fn get_bp_fingerprint(mol: &Mol) -> Result<ULongSparseIntVect> {
    get_bp_fingerprint_with_fpfn(mol, |mol, invariants| {
        sparse_counts_with_invariants(mol, invariants, false)
    })
}

/// Fingerprints the atoms' Patty-type invariants with `fpfn` instead.
/// Mirrors `Sheridan.GetBPFingerprint(mol, fpfn)`.
pub fn get_bp_fingerprint_with_fpfn<T>(
    mol: &Mol,
    fpfn: impl FnOnce(&Mol, &[u32]) -> Result<T>,
) -> Result<T> {
    let typs = patty_invariants(mol)?;
    fpfn(mol, &typs)
}

/// Sheridan's binding-property torsions: the unfolded topological-torsion
/// count fingerprint over the atoms' Patty types. Mirrors
/// `Sheridan.GetBTFingerprint`.
pub fn get_bt_fingerprint(mol: &Mol) -> Result<ULongSparseIntVect> {
    get_bt_fingerprint_with_fpfn(mol, |mol, invariants| {
        sparse_counts_with_invariants(mol, invariants, true)
    })
}

/// Fingerprints the atoms' Patty-type invariants with `fpfn` instead.
/// Mirrors `Sheridan.GetBTFingerprint(mol, fpfn)`.
pub fn get_bt_fingerprint_with_fpfn<T>(
    mol: &Mol,
    fpfn: impl FnOnce(&Mol, &[u32]) -> Result<T>,
) -> Result<T> {
    get_bp_fingerprint_with_fpfn(mol, fpfn)
}
