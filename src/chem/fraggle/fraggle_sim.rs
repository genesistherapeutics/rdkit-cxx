// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Fraggle.FraggleSim`, a pure-Python module: Fraggle
//! similarity, which compares molecules by the fragments of the query left
//! by cutting its acyclic and exocyclic ring bonds.
//!
//! Fragmentation cuts every pair of acyclic bonds, keeping terminal pieces
//! that hold more than 60% of the query's atoms; every pair of exocyclic
//! ring bonds, keeping a ring piece of more than 40%; and each such ring cut
//! with one acyclic cut, keeping two pieces of more than 60% together.
//!
//! Python caches a query's fragmentation on the molecule (`_fraggleDecomp`)
//! and the modified queries' fingerprints in the module
//! (`modified_query_fps`); the port computes them for each call, with the
//! same results. Python's `verbose` argument, which prints the matched
//! atoms, is not ported.

use std::collections::BTreeSet;

use crate::chem::rdqueries::atom_num_equals_query_atom;
use crate::chem::{
    Atom, FragmentOnBondsParams, Mol, RDKFingerprintParams, SanitizeFlags, SanitizeMolParams,
    fragment_on_bonds_with_params, get_mol_frags_as_mols_with_sanitize_frags, mol_from_smarts,
    mol_from_smiles, mol_to_smiles, rdk_fingerprint_with_outputs, rdk_fingerprint_with_params,
    sanitize_mol_with_params,
};
use crate::data_structs::{ExplicitBitVect, tanimoto_similarity, tversky_similarity};
use crate::{Error, Result};

/// The RDKit fingerprint's parameters for Fraggle: paths of up to 5 bonds,
/// 1024 bits, 2 bits per path. Mirrors `FraggleSim.rdkitFpParams`.
pub fn rdkit_fp_params() -> RDKFingerprintParams {
    RDKFingerprintParams {
        max_path: 5,
        fp_size: 1024,
        n_bits_per_hash: 2,
        ..Default::default()
    }
}

/// The type of fragmentation for acyclic cuts. Mirrors
/// `FraggleSim.FTYPE_ACYCLIC`.
pub fn ftype_acyclic() -> &'static str {
    "acyclic"
}

/// The type of fragmentation for exocyclic ring cuts. Mirrors
/// `FraggleSim.FTYPE_CYCLIC`.
pub fn ftype_cyclic() -> &'static str {
    "cyclic"
}

/// The type of fragmentation for a ring cut with an acyclic cut. Mirrors
/// `FraggleSim.FTYPE_CYCLIC_ACYCLIC`.
pub fn ftype_cyclic_acyclic() -> &'static str {
    "cyclic_and_acyclic"
}

fn pattern(smarts: &str) -> Mol {
    mol_from_smarts(smarts).expect("the Fraggle SMARTS parse")
}

/// The acyclic single bonds Fraggle cuts. Mirrors `FraggleSim.ACYC_SMARTS`.
pub fn acyc_smarts() -> Mol {
    pattern("*!@!=!#*")
}

/// The exocyclic ring bonds Fraggle cuts. Mirrors `FraggleSim.CYC_SMARTS`.
pub fn cyc_smarts() -> Mol {
    pattern("[R1,R2]@[r;!R1]")
}

/// A ring piece with attachment points on two ring atoms. Mirrors
/// `FraggleSim.cSma1`.
pub fn c_sma1() -> Mol {
    pattern("[#0][r].[r][#0]")
}

/// A ring piece with two attachment points on one ring atom. Mirrors
/// `FraggleSim.cSma2`.
pub fn c_sma2() -> Mol {
    pattern("[#0][r][#0]")
}

/// A query atom for attachment points, dummies. Mirrors
/// `FraggleSim.dummyAtomQuery`.
pub fn dummy_atom_query() -> Atom {
    atom_num_equals_query_atom(0).expect("0 is a valid atomic number query")
}

/// The molecule cut at the bonds, given as their atoms, each end capped
/// with an unlabeled dummy, as Fraggle's fragmentation SMILES for `ftype`
/// (see [`select_fragments`]); `None` if the cut gives none. A pair of
/// atoms that is not bonded is `Error::Value`. Mirrors
/// `FraggleSim.delete_bonds`.
pub fn delete_bonds(
    mol: &Mol,
    bonds: &[(u32, u32)],
    ftype: &str,
    hac: u32,
) -> Result<Option<String>> {
    let mut bond_idx = Vec::with_capacity(bonds.len());
    for &(a, b) in bonds {
        let bond = mol
            .get_bond_between_atoms(a, b)?
            .ok_or_else(|| Error::Value {
                message: format!("atoms {a} and {b} are not bonded"),
            })?;
        bond_idx.push(bond.get_idx());
    }
    let mut modified_mol = fragment_on_bonds_with_params(
        mol,
        &bond_idx,
        &FragmentOnBondsParams {
            dummy_labels: Some(vec![(0, 0); bond_idx.len()]),
            ..Default::default()
        },
    )?;
    // Python finds rings and valences without a full sanitization.
    sanitize_mol_with_params(
        &mut modified_mol,
        &SanitizeMolParams {
            sanitize_ops: SanitizeFlags::SANITIZE_PROPERTIES | SanitizeFlags::SANITIZE_SYMMRINGS,
            ..Default::default()
        },
    )?;
    let mut fragments = get_mol_frags_as_mols_with_sanitize_frags(&modified_mol, false)?;
    select_fragments(&mut fragments, ftype, hac)
}

fn n_attachments(mol: &Mol, query: &Atom) -> Result<usize> {
    Ok(mol.get_atoms_matching_query(query)?.len())
}

/// Fraggle's fragmentation SMILES from the pieces of a cut molecule of
/// `hac` atoms, or `None`: for acyclic cuts, the terminal pieces of more
/// than 3 atoms, when they hold more than 60% of the atoms; for cyclic
/// cuts, of two pieces, a valid ring cut of more than 3 atoms and 40%; for
/// both, two pieces of more than 3 atoms, each terminal or a valid ring cut,
/// holding more than 60%. Valid ring cuts are found with their rings, so
/// the pieces are changed. Another `ftype` is `Error::Value`, where Python
/// raises `NotImplementedError`. Mirrors `FraggleSim.select_fragments`.
pub fn select_fragments(fragments: &mut [Mol], ftype: &str, hac: u32) -> Result<Option<String>> {
    let query = dummy_atom_query();
    let hac = f64::from(hac);
    if ftype == ftype_acyclic() {
        let mut result = Vec::new();
        let mut result_hcount = 0;
        for f_mol in fragments.iter() {
            if n_attachments(f_mol, &query)? == 1 {
                let fhac = f_mol.get_num_atoms();
                if fhac > 3 {
                    result.push(mol_to_smiles(f_mol)?);
                    result_hcount += fhac;
                }
            }
        }
        return Ok(
            (!result.is_empty() && f64::from(result_hcount) > 0.6 * hac).then(|| result.join("."))
        );
    }
    if ftype == ftype_cyclic() {
        if fragments.len() != 2 {
            return Ok(None);
        }
        let mut result = None;
        for f_mol in fragments.iter_mut() {
            let f = mol_to_smiles(f_mol)?;
            if is_valid_ring_cut(f_mol)? {
                let result_hcount = f_mol.get_num_atoms();
                if result_hcount > 3 && f64::from(result_hcount) > 0.4 * hac {
                    result = Some(f);
                }
            }
        }
        return Ok(result);
    }
    if ftype == ftype_cyclic_acyclic() {
        let mut result = Vec::new();
        let mut result_hcount = 0;
        for f_mol in fragments.iter_mut() {
            let n = n_attachments(f_mol, &query)?;
            if n >= 3 {
                continue;
            }
            let fhac = f_mol.get_num_atoms();
            if fhac <= 3 {
                continue;
            }
            if n == 2 {
                if is_valid_ring_cut(f_mol)? {
                    result.push(mol_to_smiles(f_mol)?);
                    result_hcount += fhac;
                }
            } else if n == 1 {
                result.push(mol_to_smiles(f_mol)?);
                result_hcount += fhac;
            }
        }
        return Ok(
            (result.len() == 2 && f64::from(result_hcount) > 0.6 * hac).then(|| result.join("."))
        );
    }
    Err(Error::Value {
        message: format!("Invalid fragmentation type {ftype}"),
    })
}

/// Whether a piece is a valid ring cut, matching `c_sma1` or `c_sma2`; its
/// rings are found first. Mirrors `FraggleSim.isValidRingCut`.
pub fn is_valid_ring_cut(mol: &mut Mol) -> Result<bool> {
    sanitize_mol_with_params(
        mol,
        &SanitizeMolParams {
            sanitize_ops: SanitizeFlags::SANITIZE_SYMMRINGS,
            ..Default::default()
        },
    )?;
    Ok(mol.has_substruct_match(&c_sma1())? || mol.has_substruct_match(&c_sma2())?)
}

/// Every Fraggle fragmentation of the molecule, as fragmentation SMILES,
/// each once, sorted. Mirrors `FraggleSim.generate_fraggle_fragmentation`.
pub fn generate_fraggle_fragmentation(mol: &Mol) -> Result<Vec<String>> {
    let hac = mol.get_num_atoms();
    let pair = |m: &Vec<u32>| (m[0], m[1]);
    let acyclic: Vec<(u32, u32)> = mol
        .get_substruct_matches(&acyc_smarts())?
        .iter()
        .map(pair)
        .collect();
    let cyclic: Vec<(u32, u32)> = mol
        .get_substruct_matches(&cyc_smarts())?
        .iter()
        .map(pair)
        .collect();
    let mut out_fragments = BTreeSet::new();
    for (i, &bond1) in acyclic.iter().enumerate() {
        for &bond2 in &acyclic[i + 1..] {
            if let Some(f) = delete_bonds(mol, &[bond1, bond2], ftype_acyclic(), hac)? {
                out_fragments.insert(f);
            }
        }
    }
    for (i, &bond1) in cyclic.iter().enumerate() {
        for &bond2 in &cyclic[i + 1..] {
            let Some(f) = delete_bonds(mol, &[bond1, bond2], ftype_cyclic(), hac)? else {
                continue;
            };
            out_fragments.insert(f);
            for &abond in &acyclic {
                if let Some(f) =
                    delete_bonds(mol, &[bond1, bond2, abond], ftype_cyclic_acyclic(), hac)?
                {
                    out_fragments.insert(f);
                }
            }
        }
    }
    Ok(out_fragments.into_iter().collect())
}

/// The molecule with its atoms that are not like the substructure `subs`
/// (a SMILES) replaced: each atom whose own fingerprint bits are less than
/// 0.8 of a Tversky match to the substructure's, and every ring holding
/// one, become dummies if aromatic and scandium otherwise. Where RDKit
/// cannot sanitize the result, as Python, writes why to standard error and
/// returns the molecule unchanged. Mirrors `FraggleSim.atomContrib(subs,
/// mol)`.
pub fn atom_contrib(subs: &str, mol: &Mol) -> Result<Mol> {
    atom_contrib_with_tversky_thresh(subs, mol, 0.8)
}

/// As [`atom_contrib`], with the threshold `tversky_thresh`. Mirrors
/// `FraggleSim.atomContrib(subs, mol, tverskyThresh)`.
pub fn atom_contrib_with_tversky_thresh(subs: &str, mol: &Mol, tversky_thresh: f64) -> Result<Mol> {
    let mut p_mol = mol.clone();
    let (_, outputs) = rdk_fingerprint_with_outputs(&p_mol, &rdkit_fp_params())?;
    let a_bits = outputs.atom_bits;
    let qs_mol = mol_from_smiles(subs)?;
    let subs_fp = rdk_fingerprint_with_params(&qs_mol, &rdkit_fp_params())?;
    let mut marked = BTreeSet::new();
    for atom_idx in 0..p_mol.get_num_atoms() {
        let mut modified_fp = ExplicitBitVect::new(1024);
        let bits: Vec<u32> = a_bits[atom_idx as usize]
            .iter()
            .map(|&b| b as u32)
            .collect();
        modified_fp.set_bits_from_list(&bits)?;
        if tversky_similarity(&subs_fp, &modified_fp, 0.0, 1.0)? < tversky_thresh {
            marked.insert(atom_idx);
        }
    }
    let mut mark_ring_atoms = BTreeSet::new();
    for ring in p_mol.get_ring_info().atom_rings() {
        if ring.iter().any(|a| marked.contains(a)) {
            mark_ring_atoms.extend(ring);
        }
    }
    marked.extend(mark_ring_atoms);
    if !marked.is_empty() {
        for &idx in &marked {
            let mut atom = p_mol.get_atom_with_idx_mut(idx)?;
            if atom.get_is_aromatic() {
                atom.set_atomic_num(0);
                atom.set_no_implicit(true);
            } else {
                atom.set_atomic_num(21);
            }
        }
        let ops = SanitizeFlags::from_bits_retain(
            SanitizeFlags::SANITIZE_ALL.bits()
                ^ SanitizeFlags::SANITIZE_KEKULIZE.bits()
                ^ SanitizeFlags::SANITIZE_SETAROMATICITY.bits(),
        );
        if sanitize_mol_with_params(
            &mut p_mol,
            &SanitizeMolParams {
                sanitize_ops: ops,
                ..Default::default()
            },
        )
        .is_err()
        {
            eprintln!("Can't parse smiles: {}", mol_to_smiles(&p_mol)?);
            p_mol = mol.clone();
        }
    }
    Ok(p_mol)
}

/// The RDKit similarity of two molecules, and their Fraggle similarity for
/// one fragmentation `q_subs` of the query `q_mol`: the larger of that and
/// the Tanimoto similarity of the two molecules after `atom_contrib`. Where
/// RDKit cannot fingerprint the changed molecule, as Python, writes why to
/// standard error and gives 0. `q_smi`, the query's SMILES, only keyed
/// Python's cache. Mirrors `FraggleSim.compute_fraggle_similarity_for_subs(
/// inMol, qMol, qSmi, qSubs)`.
pub fn compute_fraggle_similarity_for_subs(
    in_mol: &Mol,
    q_mol: &Mol,
    q_smi: &str,
    q_subs: &str,
) -> Result<(f64, f64)> {
    compute_fraggle_similarity_for_subs_with_tversky_thresh(in_mol, q_mol, q_smi, q_subs, 0.8)
}

/// As [`compute_fraggle_similarity_for_subs`], with `tversky_thresh` for
/// `atom_contrib`. Mirrors `FraggleSim.compute_fraggle_similarity_for_subs(
/// inMol, qMol, qSmi, qSubs, tverskyThresh)`.
pub fn compute_fraggle_similarity_for_subs_with_tversky_thresh(
    in_mol: &Mol,
    q_mol: &Mol,
    _q_smi: &str,
    q_subs: &str,
    tversky_thresh: f64,
) -> Result<(f64, f64)> {
    let params = rdkit_fp_params();
    let q_fp = rdk_fingerprint_with_params(q_mol, &params)?;
    let i_fp = rdk_fingerprint_with_params(in_mol, &params)?;
    let rdkit_sim = tanimoto_similarity(&q_fp, &i_fp)?;
    let qm_mol = atom_contrib_with_tversky_thresh(q_subs, q_mol, tversky_thresh)?;
    let qm_mol_fp = rdk_fingerprint_with_params(&qm_mol, &params)?;
    let rm_mol = atom_contrib_with_tversky_thresh(q_subs, in_mol, tversky_thresh)?;
    let fraggle_sim = match rdk_fingerprint_with_params(&rm_mol, &params) {
        Ok(rm_mol_fp) => tanimoto_similarity(&qm_mol_fp, &rm_mol_fp)?.max(rdkit_sim),
        Err(_) => {
            eprintln!("Can't generate fp for: {}", mol_to_smiles(&rm_mol)?);
            0.0
        }
    };
    Ok((rdkit_sim, fraggle_sim))
}

/// The Fraggle similarity of two molecules, and the query's fragmentation
/// that gives it; `None` when no fragmentation gives more than 0. Mirrors
/// `FraggleSim.GetFraggleSimilarity(queryMol, refMol)`.
pub fn get_fraggle_similarity(query_mol: &Mol, ref_mol: &Mol) -> Result<(f64, Option<String>)> {
    get_fraggle_similarity_with_tversky_thresh(query_mol, ref_mol, 0.8)
}

/// As [`get_fraggle_similarity`], with `tversky_thresh` for `atom_contrib`.
/// Mirrors `FraggleSim.GetFraggleSimilarity(queryMol, refMol,
/// tverskyThresh)`.
pub fn get_fraggle_similarity_with_tversky_thresh(
    query_mol: &Mol,
    ref_mol: &Mol,
    tversky_thresh: f64,
) -> Result<(f64, Option<String>)> {
    let frags = generate_fraggle_fragmentation(query_mol)?;
    let q_smi = mol_to_smiles(query_mol)?;
    let mut result = 0.0;
    let mut best_match = None;
    for frag in frags {
        let (_, fragsim) = compute_fraggle_similarity_for_subs_with_tversky_thresh(
            ref_mol,
            query_mol,
            &q_smi,
            &frag,
            tversky_thresh,
        )?;
        if fragsim > result {
            result = fragsim;
            best_match = Some(frag);
        }
    }
    Ok((result, best_match))
}
