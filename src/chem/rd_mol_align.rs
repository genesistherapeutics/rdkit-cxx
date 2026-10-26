// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolAlign`: RMSD alignment of molecules and
//! conformers, and Open3DAlign (O3A).
//!
//! ```no_run
//! use rdkit_cxx::chem::rd_dist_geom::embed_multiple_confs;
//! use rdkit_cxx::chem::rd_mol_align::{align_mol, get_all_conformer_best_rms};
//! use rdkit_cxx::chem::{add_hs, mol_from_smiles};
//!
//! let mut mol = add_hs(&mol_from_smiles("CCOC(=O)c1ccccc1").unwrap()).unwrap();
//! embed_multiple_confs(&mut mol).unwrap();
//! // The RMSD between each pair of conformers, after aligning them.
//! let rms = get_all_conformer_best_rms(&mol).unwrap();
//! // Align a copy of the first conformer onto the molecule.
//! let mut probe = mol.clone();
//! align_mol(&mut probe, &mol).unwrap();
//! ```
//!
//! Atom maps are pairs of (probe atom, reference atom). The probe and the
//! reference are separate molecules: to compare conformers of one molecule,
//! pass a copy of it as the reference, or use the functions over all
//! conformers. Transforms are 4 x 4 matrices, as rows.
//!
//! RDKit's threads (`num_threads`) call `std::terminate` on any exception,
//! so the wrapper checks atom maps and weights before they start, and runs
//! at most as many threads as the hardware has. Finding the atom maps of a
//! molecule's symmetries writes a static RDKit variable, under a lock.

use std::marker::PhantomData;
use std::ptr;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::mol_align::ffi;
use crate::force_field::MMFFMolProperties;

fn to_pairs(pairs: &[(u32, u32)]) -> Vec<IndexPair> {
    pairs
        .iter()
        .map(|&(first, second)| IndexPair { first, second })
        .collect()
}

fn from_pairs(pairs: Vec<IndexPair>) -> Vec<(u32, u32)> {
    pairs.into_iter().map(|p| (p.first, p.second)).collect()
}

/// A best alignment: the RMSD, the transform (4 x 4, as rows), and the atom
/// map, as (probe, reference) pairs. A Rust addition, naming Python's tuple.
pub type BestAlignment = (f64, Vec<Vec<f64>>, Vec<(u32, u32)>);

fn rows(values: &[f64]) -> Vec<Vec<f64>> {
    values.chunks(4).map(<[f64]>::to_vec).collect()
}

/// Parameters for the best-alignment functions. Mirrors
/// `rdMolAlign.BestAlignmentParams`, plus the keyword arguments of
/// `GetBestAlignmentTransform` (`prbCid`, `refCid`, `reflect`, `maxIters`)
/// and `GetBestRMS` (`prbId` and `refId`, as `prb_cid` and `ref_cid`)
/// (A3). The functions over all conformers read neither.
#[derive(Clone, Debug, PartialEq)]
pub struct BestAlignmentParams {
    /// The most substructure matches to try when no map is given.
    pub max_matches: i32,
    /// Matches conjugated terminal groups (nitro, carboxylate) both ways.
    pub symmetrize_conjugated_terminal_groups: bool,
    /// Leaves the hydrogens out of the matches found (not of a given map).
    pub ignore_hs: bool,
    /// The threads to use; 0 or less counts back from the hardware's.
    pub num_threads: i32,
    /// The atom maps to try; empty finds them by substructure matching.
    pub map: Vec<Vec<(u32, u32)>>,
    /// A weight for each pair of a map; empty weighs them equally.
    pub weights: Vec<f64>,
    /// The probe's conformer; -1 for the first.
    pub prb_cid: i32,
    /// The reference's conformer; -1 for the first.
    pub ref_cid: i32,
    /// Reflects the probe.
    pub reflect: bool,
    /// The most iterations of the RMSD minimization.
    pub max_iters: u32,
}

impl Default for BestAlignmentParams {
    fn default() -> BestAlignmentParams {
        let data = expect_ok(ffi::best_alignment_params_default());
        BestAlignmentParams {
            max_matches: data.max_matches,
            symmetrize_conjugated_terminal_groups: data.symmetrize_conjugated_terminal_groups,
            ignore_hs: data.ignore_hs,
            num_threads: data.num_threads,
            map: Vec::new(),
            weights: Vec::new(),
            prb_cid: data.prb_cid,
            ref_cid: data.ref_cid,
            reflect: data.reflect,
            max_iters: data.max_iters,
        }
    }
}

impl BestAlignmentParams {
    fn to_ffi(&self) -> ffi::BestAlignmentParamsData {
        ffi::BestAlignmentParamsData {
            max_matches: self.max_matches,
            symmetrize_conjugated_terminal_groups: self.symmetrize_conjugated_terminal_groups,
            ignore_hs: self.ignore_hs,
            num_threads: self.num_threads,
            map_pairs: self.map.iter().flat_map(|m| to_pairs(m)).collect(),
            map_counts: self
                .map
                .iter()
                .map(|m| u32::try_from(m.len()).expect("an atom map that long"))
                .collect(),
            weights: self.weights.clone(),
            prb_cid: self.prb_cid,
            ref_cid: self.ref_cid,
            reflect: self.reflect,
            max_iters: self.max_iters,
        }
    }

    /// The keyword forms' defaults, which keep the hydrogens.
    fn keyword_defaults() -> BestAlignmentParams {
        BestAlignmentParams {
            ignore_hs: false,
            ..Default::default()
        }
    }
}

/// Parameters for `get_alignment_transform_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetAlignmentTransformParams {
    /// The probe's conformer; -1 for the first.
    pub prb_cid: i32,
    /// The reference's conformer; -1 for the first.
    pub ref_cid: i32,
    /// The atom map; empty finds one by substructure matching.
    pub atom_map: Vec<(u32, u32)>,
    /// A weight for each pair (each probe atom without a map); empty weighs
    /// them equally.
    pub weights: Vec<f64>,
    /// Reflects the probe.
    pub reflect: bool,
    /// The most iterations of the RMSD minimization.
    pub max_iters: u32,
}

impl Default for GetAlignmentTransformParams {
    fn default() -> GetAlignmentTransformParams {
        GetAlignmentTransformParams {
            prb_cid: -1,
            ref_cid: -1,
            atom_map: Vec::new(),
            weights: Vec::new(),
            reflect: false,
            max_iters: 50,
        }
    }
}

/// The RMSD and the transform that align the probe's first conformer to
/// the reference's. Mirrors `rdMolAlign.GetAlignmentTransform`.
pub fn get_alignment_transform(prb_mol: &Mol, ref_mol: &Mol) -> Result<(f64, Vec<Vec<f64>>)> {
    get_alignment_transform_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.GetAlignmentTransform(prbMol, refMol, prbCid,
/// refCid, atomMap, weights, reflect, maxIters)`. A map's atom out of range
/// is `Error::Index`, and the wrong number of weights `Error::Value`; no
/// substructure match is `Error::Runtime`.
pub fn get_alignment_transform_with_params(
    prb_mol: &Mol,
    ref_mol: &Mol,
    params: &GetAlignmentTransformParams,
) -> Result<(f64, Vec<Vec<f64>>)> {
    let mut trans = Vec::new();
    let rmsd = ffi::get_alignment_transform(
        prb_mol.as_ffi(),
        ref_mol.as_ffi(),
        params.prb_cid,
        params.ref_cid,
        &to_pairs(&params.atom_map),
        &params.weights,
        params.reflect,
        params.max_iters,
        &mut trans,
    )?;
    Ok((rmsd, rows(&trans)))
}

/// Aligns the probe's first conformer to the reference's, and returns the
/// RMSD. Mirrors `rdMolAlign.AlignMol`.
pub fn align_mol(prb_mol: &mut Mol, ref_mol: &Mol) -> Result<f64> {
    align_mol_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.AlignMol(prbMol, refMol, prbCid, refCid, atomMap,
/// weights, reflect, maxIters)`, with the errors of
/// `get_alignment_transform_with_params`.
pub fn align_mol_with_params(
    prb_mol: &mut Mol,
    ref_mol: &Mol,
    params: &GetAlignmentTransformParams,
) -> Result<f64> {
    Ok(ffi::align_mol(
        prb_mol.pin_mut(),
        ref_mol.as_ffi(),
        params.prb_cid,
        params.ref_cid,
        &to_pairs(&params.atom_map),
        &params.weights,
        params.reflect,
        params.max_iters,
    )?)
}

/// The best alignment of the probe to the reference over their atom maps:
/// the RMSD, the transform, and the map. Mirrors
/// `rdMolAlign.GetBestAlignmentTransform`, whose keyword form keeps the
/// hydrogens (`ignore_hs` false).
pub fn get_best_alignment_transform(prb_mol: &Mol, ref_mol: &Mol) -> Result<BestAlignment> {
    get_best_alignment_transform_with_params(
        prb_mol,
        ref_mol,
        &BestAlignmentParams::keyword_defaults(),
    )
}

/// Mirrors `rdMolAlign.GetBestAlignmentTransform(prbMol, refMol, params,
/// prbCid, refCid, reflect, maxIters)`, and its keyword form. A map's atom
/// out of range is `Error::Index`, and weights other than one per pair of
/// each map `Error::Value`; no substructure match is `Error::Runtime`.
pub fn get_best_alignment_transform_with_params(
    prb_mol: &Mol,
    ref_mol: &Mol,
    params: &BestAlignmentParams,
) -> Result<BestAlignment> {
    let mut trans = Vec::new();
    let mut best_match = Vec::new();
    let rmsd = ffi::get_best_alignment_transform(
        prb_mol.as_ffi(),
        ref_mol.as_ffi(),
        &params.to_ffi(),
        &mut trans,
        &mut best_match,
    )?;
    Ok((rmsd, rows(&trans), from_pairs(best_match)))
}

/// Aligns the probe's first conformer to the reference's by its best atom
/// map, and returns the RMSD. Mirrors `rdMolAlign.GetBestRMS`, whose
/// keyword form keeps the hydrogens (`ignore_hs` false).
pub fn get_best_rms(prb_mol: &mut Mol, ref_mol: &Mol) -> Result<f64> {
    get_best_rms_with_params(prb_mol, ref_mol, &BestAlignmentParams::keyword_defaults())
}

/// Mirrors `rdMolAlign.GetBestRMS(prbMol, refMol, params, prbId, refId)`,
/// and its keyword form, with the errors of
/// `get_best_alignment_transform_with_params`.
pub fn get_best_rms_with_params(
    prb_mol: &mut Mol,
    ref_mol: &Mol,
    params: &BestAlignmentParams,
) -> Result<f64> {
    let data = params.to_ffi();
    Ok(ffi::get_best_rms(
        prb_mol.pin_mut(),
        ref_mol.as_ffi(),
        &data,
    )?)
}

/// The best RMSD between each pair of the molecule's conformers, in the
/// order (1, 0), (2, 0), (2, 1), (3, 0), and so on. Mirrors
/// `rdMolAlign.GetAllConformerBestRMS`.
pub fn get_all_conformer_best_rms(mol: &Mol) -> Result<Vec<f64>> {
    get_all_conformer_best_rms_with_params(mol, &Default::default())
}

/// Mirrors `rdMolAlign.GetAllConformerBestRMS(mol, params)`, and its keyword
/// form, with the errors of `get_best_alignment_transform_with_params`.
pub fn get_all_conformer_best_rms_with_params(
    mol: &Mol,
    params: &BestAlignmentParams,
) -> Result<Vec<f64>> {
    Ok(ffi::get_all_conformer_best_rms(
        mol.as_ffi(),
        &params.to_ffi(),
    )?)
}

/// The best RMSD between each of the reference's conformers and each of
/// the probe's, by reference conformer. Mirrors
/// `rdMolAlign.GetAllConformerBestRMSToRef`.
pub fn get_all_conformer_best_rms_to_ref(prb_mol: &Mol, ref_mol: &Mol) -> Result<Vec<f64>> {
    get_all_conformer_best_rms_to_ref_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.GetAllConformerBestRMSToRef(prbMol, refMol, params)`,
/// with the errors of `get_best_alignment_transform_with_params`.
pub fn get_all_conformer_best_rms_to_ref_with_params(
    prb_mol: &Mol,
    ref_mol: &Mol,
    params: &BestAlignmentParams,
) -> Result<Vec<f64>> {
    Ok(ffi::get_all_conformer_best_rms_to_ref(
        prb_mol.as_ffi(),
        ref_mol.as_ffi(),
        &params.to_ffi(),
    )?)
}

/// Parameters for `calc_rms_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CalcRMSParams {
    /// The probe's conformer; -1 for the first.
    pub prb_id: i32,
    /// The reference's conformer; -1 for the first.
    pub ref_id: i32,
    /// The atom maps to try; empty finds them by substructure matching.
    pub map: Vec<Vec<(u32, u32)>>,
    /// The most substructure matches to try when no map is given.
    pub max_matches: i32,
    /// Matches conjugated terminal groups (nitro, carboxylate) both ways.
    pub symmetrize_conjugated_terminal_groups: bool,
    /// A weight for each pair of a map; empty weighs them equally.
    pub weights: Vec<f64>,
}

impl Default for CalcRMSParams {
    fn default() -> CalcRMSParams {
        CalcRMSParams {
            prb_id: -1,
            ref_id: -1,
            map: Vec::new(),
            max_matches: 1_000_000,
            symmetrize_conjugated_terminal_groups: true,
            weights: Vec::new(),
        }
    }
}

/// The best RMSD between the conformers over the atom maps, without
/// aligning them. Mirrors `rdMolAlign.CalcRMS`.
pub fn calc_rms(prb_mol: &Mol, ref_mol: &Mol) -> Result<f64> {
    calc_rms_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.CalcRMS(prbMol, refMol, prbId, refId, map,
/// maxMatches, symmetrizeConjugatedTerminalGroups, weights)`, with the
/// errors of `get_best_alignment_transform_with_params`.
pub fn calc_rms_with_params(prb_mol: &Mol, ref_mol: &Mol, params: &CalcRMSParams) -> Result<f64> {
    let best = BestAlignmentParams {
        max_matches: params.max_matches,
        symmetrize_conjugated_terminal_groups: params.symmetrize_conjugated_terminal_groups,
        map: params.map.clone(),
        weights: params.weights.clone(),
        prb_cid: params.prb_id,
        ref_cid: params.ref_id,
        num_threads: 1,
        ..BestAlignmentParams::keyword_defaults()
    };
    Ok(ffi::calc_rms(
        prb_mol.as_ffi(),
        ref_mol.as_ffi(),
        &best.to_ffi(),
    )?)
}

/// Parameters for `align_mol_conformers_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AlignMolConformersParams {
    /// The atoms to align; empty for all.
    pub atom_ids: Vec<u32>,
    /// The conformers, the first the reference; empty for all, the first
    /// conformer the reference.
    pub conf_ids: Vec<i32>,
    /// A weight for each atom aligned; empty weighs them equally.
    pub weights: Vec<f64>,
    /// Reflects the conformers.
    pub reflect: bool,
    /// The most iterations of the RMSD minimization.
    pub max_iters: u32,
}

impl Default for AlignMolConformersParams {
    fn default() -> AlignMolConformersParams {
        AlignMolConformersParams {
            atom_ids: Vec::new(),
            conf_ids: Vec::new(),
            weights: Vec::new(),
            reflect: false,
            max_iters: 50,
        }
    }
}

/// Aligns each of the molecule's conformers to its first, and returns the
/// RMSD of each, which Python's `RMSlist` collects. Mirrors
/// `rdMolAlign.AlignMolConformers`.
pub fn align_mol_conformers(mol: &mut Mol) -> Result<Vec<f64>> {
    align_mol_conformers_with_params(mol, &Default::default())
}

/// Mirrors `rdMolAlign.AlignMolConformers(mol, atomIds, confIds, weights,
/// reflect, maxIters, RMSlist)`. An atom out of range is `Error::Index`; a
/// missing conformer is `Error::Value`, and the wrong number of weights
/// `Error::Runtime`.
pub fn align_mol_conformers_with_params(
    mol: &mut Mol,
    params: &AlignMolConformersParams,
) -> Result<Vec<f64>> {
    // Python passes the ids as unsigned int, and RDKit takes them as int.
    let conf_ids: Vec<u32> = params.conf_ids.iter().map(|&id| id as u32).collect();
    Ok(ffi::align_mol_conformers(
        mol.pin_mut(),
        &params.atom_ids,
        &conf_ids,
        &params.weights,
        params.reflect,
        params.max_iters,
    )?)
}

/// Parameters for `random_transform_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct RandomTransformParams {
    /// The conformer; -1 for the first.
    pub cid: i32,
    /// A positive seed reseeds RDKit's random number generator first.
    pub seed: i32,
}

impl Default for RandomTransformParams {
    fn default() -> RandomTransformParams {
        RandomTransformParams { cid: -1, seed: -1 }
    }
}

/// Rotates and translates the molecule's first conformer at random. Mirrors
/// `rdMolAlign.RandomTransform`.
pub fn random_transform(mol: &mut Mol) -> Result<()> {
    random_transform_with_params(mol, &Default::default())
}

/// Mirrors `rdMolAlign.RandomTransform(mol, cid, seed)`. RDKit's random
/// number generator is shared by the process; a conformer of no atoms is
/// left as it is.
pub fn random_transform_with_params(mol: &mut Mol, params: &RandomTransformParams) -> Result<()> {
    Ok(ffi::random_transform(
        mol.pin_mut(),
        params.cid,
        params.seed,
    )?)
}

/// An Open3DAlign alignment of a probe conformer to a reference conformer.
/// Mirrors `rdMolAlign.O3A`.
///
/// RDKit's O3A keeps pointers to both molecules, and Python's does not keep
/// them alive: an O3A whose probe Python has freed reads freed memory. This
/// one borrows the probe exclusively and the reference for as long as it
/// lives:
///
/// ```compile_fail
/// use rdkit_cxx::chem::mol_from_smiles;
/// use rdkit_cxx::chem::rd_mol_align::get_crippen_o3a;
///
/// let mut prb = mol_from_smiles("CCO").unwrap();
/// let reference = prb.clone();
/// let mut o3a = get_crippen_o3a(&mut prb, &reference).unwrap();
/// drop(prb); // the O3A still borrows the probe
/// o3a.align().unwrap();
/// ```
pub struct O3A<'m> {
    ptr: UniquePtr<ffi::O3AHolder>,
    _mols: PhantomData<(&'m mut Mol, &'m Mol)>,
}

impl std::fmt::Debug for O3A<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("O3A").finish_non_exhaustive()
    }
}

impl O3A<'_> {
    /// Aligns the probe's conformer to the reference's by the matched
    /// atoms, and returns the RMSD. Mirrors `O3A.Align`.
    pub fn align(&mut self) -> Result<f64> {
        Ok(ffi::o3a_align(self.ptr.pin_mut())?)
    }

    /// The RMSD and the transform `align` applies, without applying it.
    /// Mirrors `O3A.Trans`.
    pub fn trans(&self) -> Result<(f64, Vec<Vec<f64>>)> {
        let mut trans = Vec::new();
        let rmsd = ffi::o3a_trans(self.raw(), &mut trans)?;
        Ok((rmsd, rows(&trans)))
    }

    /// The O3A score of the alignment. Mirrors `O3A.Score`.
    pub fn score(&self) -> f64 {
        expect_ok(ffi::o3a_score(self.raw()))
    }

    /// The matched atoms, as (probe, reference) pairs. Mirrors
    /// `O3A.Matches`.
    pub fn matches(&self) -> Vec<(u32, u32)> {
        from_pairs(expect_ok(ffi::o3a_matches(self.raw())))
    }

    /// The weight of each matched pair. Mirrors `O3A.Weights`.
    pub fn weights(&self) -> Vec<f64> {
        expect_ok(ffi::o3a_weights(self.raw()))
    }

    fn raw(&self) -> &ffi::O3AHolder {
        self.ptr.as_ref().expect("O3A is never null")
    }
}

fn o3a_list<'m>(mut list: UniquePtr<ffi::O3AList>) -> Result<Vec<O3A<'m>>> {
    let n = ffi::o3a_list_size(&list)?;
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        result.push(O3A {
            ptr: ffi::o3a_list_take(list.pin_mut(), i)?,
            _mols: PhantomData,
        });
    }
    Ok(result)
}

fn single<'m>(list: UniquePtr<ffi::O3AList>) -> Result<O3A<'m>> {
    Ok(o3a_list(list)?.pop().expect("one O3A"))
}

struct O3AOptions<'a> {
    prb_cid: i32,
    ref_cid: i32,
    reflect: bool,
    max_iters: u32,
    options: u32,
    num_threads: i32,
    constraint_map: &'a [(u32, u32)],
    constraint_weights: &'a [f64],
}

impl O3AOptions<'_> {
    fn to_ffi(&self) -> ffi::O3AOptionsData {
        ffi::O3AOptionsData {
            prb_cid: self.prb_cid,
            ref_cid: self.ref_cid,
            reflect: self.reflect,
            max_iters: self.max_iters,
            options: self.options,
            num_threads: self.num_threads,
            constraint_map: to_pairs(self.constraint_map),
            has_constraint_weights: !self.constraint_weights.is_empty(),
            constraint_weights: self.constraint_weights.to_vec(),
        }
    }
}

fn props_ptr(props: Option<&MMFFMolProperties>) -> *const ffi::MMFFMolProperties {
    props.map_or(ptr::null(), |p| p.raw() as *const _)
}

fn mmff_o3a(
    prb_mol: &mut Mol,
    ref_mol: &Mol,
    prb_props: Option<&MMFFMolProperties>,
    ref_props: Option<&MMFFMolProperties>,
    for_probe_confs: bool,
    options: &O3AOptions<'_>,
) -> Result<UniquePtr<ffi::O3AList>> {
    // SAFETY: the pointers are null or come from live references.
    Ok(unsafe {
        ffi::get_mmff_o3a(
            prb_mol.pin_mut(),
            ref_mol.as_ffi(),
            props_ptr(prb_props),
            props_ptr(ref_props),
            for_probe_confs,
            &options.to_ffi(),
        )
    }?)
}

fn logp(contribs: &Option<Vec<(f64, f64)>>) -> (bool, Vec<f64>) {
    match contribs {
        Some(c) => (true, c.iter().map(|&(logp, _)| logp).collect()),
        None => (false, Vec::new()),
    }
}

fn crippen_o3a(
    prb_mol: &mut Mol,
    ref_mol: &Mol,
    prb_contribs: &Option<Vec<(f64, f64)>>,
    ref_contribs: &Option<Vec<(f64, f64)>>,
    for_probe_confs: bool,
    options: &O3AOptions<'_>,
) -> Result<UniquePtr<ffi::O3AList>> {
    let (has_prb, prb) = logp(prb_contribs);
    let (has_ref, refs) = logp(ref_contribs);
    Ok(ffi::get_crippen_o3a(
        prb_mol.pin_mut(),
        ref_mol.as_ffi(),
        has_prb,
        &prb,
        has_ref,
        &refs,
        for_probe_confs,
        &options.to_ffi(),
    )?)
}

/// Parameters for `get_o3a_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug)]
pub struct GetO3AParams<'a> {
    /// The probe's MMFF properties; `None` computes them, from a copy.
    pub prb_mmff_mol_properties: Option<&'a MMFFMolProperties>,
    /// The reference's MMFF properties; `None` computes them, from a copy.
    pub ref_mmff_mol_properties: Option<&'a MMFFMolProperties>,
    /// The probe's conformer; -1 for the first.
    pub prb_cid: i32,
    /// The reference's conformer; -1 for the first.
    pub ref_cid: i32,
    /// Reflects the probe.
    pub reflect: bool,
    /// The most iterations of the RMSD minimization.
    pub max_iters: u32,
    /// O3A's option bits: the accuracy (0 to 3) and `O3_LOCAL_ONLY` (4).
    pub options: u32,
    /// Pairs of (probe, reference) heavy atoms to keep matched.
    pub constraint_map: Vec<(u32, u32)>,
    /// A weight for each constraint; empty for RDKit's default, 100.
    pub constraint_weights: Vec<f64>,
}

impl Default for GetO3AParams<'_> {
    fn default() -> Self {
        GetO3AParams {
            prb_mmff_mol_properties: None,
            ref_mmff_mol_properties: None,
            prb_cid: -1,
            ref_cid: -1,
            reflect: false,
            max_iters: 50,
            options: 0,
            constraint_map: Vec::new(),
            constraint_weights: Vec::new(),
        }
    }
}

/// Aligns the probe's first conformer to the reference's with O3A, scoring
/// atoms by MMFF properties. Mirrors `rdMolAlign.GetO3A`.
///
/// Properties that do not type every atom, or are for a molecule of
/// another size, and molecules with dummy atoms or coordinates beyond 1e9,
/// which RDKit would mishandle, are `Error::Value`; so are constraints on
/// hydrogens and the wrong number of constraint weights, and a
/// constraint's atom out of range is `Error::Index`.
pub fn get_o3a<'m>(prb_mol: &'m mut Mol, ref_mol: &'m Mol) -> Result<O3A<'m>> {
    get_o3a_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.GetO3A(prbMol, refMol, prbMMFFMolProperties,
/// refMMFFMolProperties, prbCid, refCid, reflect, maxIters, options,
/// constraintMap, constraintWeights)`, with the errors of `get_o3a`.
pub fn get_o3a_with_params<'m>(
    prb_mol: &'m mut Mol,
    ref_mol: &'m Mol,
    params: &GetO3AParams<'_>,
) -> Result<O3A<'m>> {
    let options = O3AOptions {
        prb_cid: params.prb_cid,
        ref_cid: params.ref_cid,
        reflect: params.reflect,
        max_iters: params.max_iters,
        options: params.options,
        num_threads: 1,
        constraint_map: &params.constraint_map,
        constraint_weights: &params.constraint_weights,
    };
    single(mmff_o3a(
        prb_mol,
        ref_mol,
        params.prb_mmff_mol_properties,
        params.ref_mmff_mol_properties,
        false,
        &options,
    )?)
}

/// Parameters for `get_crippen_o3a_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetCrippenO3AParams {
    /// The probe's Crippen (logP, MR) contributions, by atom; `None`, or
    /// other than one per atom, computes them, as in Python.
    pub prb_crippen_contribs: Option<Vec<(f64, f64)>>,
    /// The reference's Crippen contributions, likewise.
    pub ref_crippen_contribs: Option<Vec<(f64, f64)>>,
    /// The probe's conformer; -1 for the first.
    pub prb_cid: i32,
    /// The reference's conformer; -1 for the first.
    pub ref_cid: i32,
    /// Reflects the probe.
    pub reflect: bool,
    /// The most iterations of the RMSD minimization.
    pub max_iters: u32,
    /// O3A's option bits: the accuracy (0 to 3) and `O3_LOCAL_ONLY` (4).
    pub options: u32,
    /// Pairs of (probe, reference) heavy atoms to keep matched.
    pub constraint_map: Vec<(u32, u32)>,
    /// A weight for each constraint; empty for RDKit's default, 100.
    pub constraint_weights: Vec<f64>,
}

impl Default for GetCrippenO3AParams {
    fn default() -> GetCrippenO3AParams {
        GetCrippenO3AParams {
            prb_crippen_contribs: None,
            ref_crippen_contribs: None,
            prb_cid: -1,
            ref_cid: -1,
            reflect: false,
            max_iters: 50,
            options: 0,
            constraint_map: Vec::new(),
            constraint_weights: Vec::new(),
        }
    }
}

/// Aligns the probe's first conformer to the reference's with O3A, scoring
/// atoms by their Crippen logP contributions. Mirrors
/// `rdMolAlign.GetCrippenO3A`, with the errors of `get_o3a`; given
/// contributions that are not finite or beyond 100 are `Error::Value`, since
/// RDKit sums costs made from them in `int`.
pub fn get_crippen_o3a<'m>(prb_mol: &'m mut Mol, ref_mol: &'m Mol) -> Result<O3A<'m>> {
    get_crippen_o3a_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.GetCrippenO3A(prbMol, refMol, prbCrippenContribs,
/// refCrippenContribs, prbCid, refCid, reflect, maxIters, options,
/// constraintMap, constraintWeights)`.
pub fn get_crippen_o3a_with_params<'m>(
    prb_mol: &'m mut Mol,
    ref_mol: &'m Mol,
    params: &GetCrippenO3AParams,
) -> Result<O3A<'m>> {
    let options = O3AOptions {
        prb_cid: params.prb_cid,
        ref_cid: params.ref_cid,
        reflect: params.reflect,
        max_iters: params.max_iters,
        options: params.options,
        num_threads: 1,
        constraint_map: &params.constraint_map,
        constraint_weights: &params.constraint_weights,
    };
    single(crippen_o3a(
        prb_mol,
        ref_mol,
        &params.prb_crippen_contribs,
        &params.ref_crippen_contribs,
        false,
        &options,
    )?)
}

/// Parameters for `get_o3a_for_probe_confs_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug)]
pub struct GetO3AForProbeConfsParams<'a> {
    /// The threads to use; 0 or less counts back from the hardware's.
    pub num_threads: i32,
    /// The probe's MMFF properties; `None` computes them, from a copy.
    pub prb_mmff_mol_properties: Option<&'a MMFFMolProperties>,
    /// The reference's MMFF properties; `None` computes them, from a copy.
    pub ref_mmff_mol_properties: Option<&'a MMFFMolProperties>,
    /// The reference's conformer; -1 for the first.
    pub ref_cid: i32,
    /// Reflects the probe.
    pub reflect: bool,
    /// The most iterations of the RMSD minimization.
    pub max_iters: u32,
    /// O3A's option bits: the accuracy (0 to 3) and `O3_LOCAL_ONLY` (4).
    pub options: u32,
    /// Pairs of (probe, reference) heavy atoms to keep matched.
    pub constraint_map: Vec<(u32, u32)>,
    /// A weight for each constraint; empty for RDKit's default, 100.
    pub constraint_weights: Vec<f64>,
}

impl Default for GetO3AForProbeConfsParams<'_> {
    fn default() -> Self {
        GetO3AForProbeConfsParams {
            num_threads: 1,
            prb_mmff_mol_properties: None,
            ref_mmff_mol_properties: None,
            ref_cid: -1,
            reflect: false,
            max_iters: 50,
            options: 0,
            constraint_map: Vec::new(),
            constraint_weights: Vec::new(),
        }
    }
}

/// An O3A alignment of each of the probe's conformers to the reference's
/// first, scoring atoms by MMFF properties. Mirrors
/// `rdMolAlign.GetO3AForProbeConfs`, with the errors of `get_o3a`.
pub fn get_o3a_for_probe_confs<'m>(prb_mol: &'m mut Mol, ref_mol: &'m Mol) -> Result<Vec<O3A<'m>>> {
    get_o3a_for_probe_confs_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.GetO3AForProbeConfs(prbMol, refMol, numThreads,
/// prbMMFFMolProperties, refMMFFMolProperties, refCid, reflect, maxIters,
/// options, constraintMap, constraintWeights)`.
pub fn get_o3a_for_probe_confs_with_params<'m>(
    prb_mol: &'m mut Mol,
    ref_mol: &'m Mol,
    params: &GetO3AForProbeConfsParams<'_>,
) -> Result<Vec<O3A<'m>>> {
    let options = O3AOptions {
        prb_cid: -1,
        ref_cid: params.ref_cid,
        reflect: params.reflect,
        max_iters: params.max_iters,
        options: params.options,
        num_threads: params.num_threads,
        constraint_map: &params.constraint_map,
        constraint_weights: &params.constraint_weights,
    };
    o3a_list(mmff_o3a(
        prb_mol,
        ref_mol,
        params.prb_mmff_mol_properties,
        params.ref_mmff_mol_properties,
        true,
        &options,
    )?)
}

/// Parameters for `get_crippen_o3a_for_probe_confs_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetCrippenO3AForProbeConfsParams {
    /// The threads to use; 0 or less counts back from the hardware's.
    pub num_threads: i32,
    /// The probe's Crippen (logP, MR) contributions, by atom; `None`, or
    /// other than one per atom, computes them, as in Python.
    pub prb_crippen_contribs: Option<Vec<(f64, f64)>>,
    /// The reference's Crippen contributions, likewise.
    pub ref_crippen_contribs: Option<Vec<(f64, f64)>>,
    /// The reference's conformer; -1 for the first.
    pub ref_cid: i32,
    /// Reflects the probe.
    pub reflect: bool,
    /// The most iterations of the RMSD minimization.
    pub max_iters: u32,
    /// O3A's option bits: the accuracy (0 to 3) and `O3_LOCAL_ONLY` (4).
    pub options: u32,
    /// Pairs of (probe, reference) heavy atoms to keep matched.
    pub constraint_map: Vec<(u32, u32)>,
    /// A weight for each constraint; empty for RDKit's default, 100.
    pub constraint_weights: Vec<f64>,
}

impl Default for GetCrippenO3AForProbeConfsParams {
    fn default() -> GetCrippenO3AForProbeConfsParams {
        GetCrippenO3AForProbeConfsParams {
            num_threads: 1,
            prb_crippen_contribs: None,
            ref_crippen_contribs: None,
            ref_cid: -1,
            reflect: false,
            max_iters: 50,
            options: 0,
            constraint_map: Vec::new(),
            constraint_weights: Vec::new(),
        }
    }
}

/// An O3A alignment of each of the probe's conformers to the reference's
/// first, scoring atoms by their Crippen logP contributions. Mirrors
/// `rdMolAlign.GetCrippenO3AForProbeConfs`, with the errors of
/// `get_crippen_o3a`.
pub fn get_crippen_o3a_for_probe_confs<'m>(
    prb_mol: &'m mut Mol,
    ref_mol: &'m Mol,
) -> Result<Vec<O3A<'m>>> {
    get_crippen_o3a_for_probe_confs_with_params(prb_mol, ref_mol, &Default::default())
}

/// Mirrors `rdMolAlign.GetCrippenO3AForProbeConfs(prbMol, refMol,
/// numThreads, prbCrippenContribs, refCrippenContribs, refCid, reflect,
/// maxIters, options, constraintMap, constraintWeights)`.
pub fn get_crippen_o3a_for_probe_confs_with_params<'m>(
    prb_mol: &'m mut Mol,
    ref_mol: &'m Mol,
    params: &GetCrippenO3AForProbeConfsParams,
) -> Result<Vec<O3A<'m>>> {
    let options = O3AOptions {
        prb_cid: -1,
        ref_cid: params.ref_cid,
        reflect: params.reflect,
        max_iters: params.max_iters,
        options: params.options,
        num_threads: params.num_threads,
        constraint_map: &params.constraint_map,
        constraint_weights: &params.constraint_weights,
    };
    o3a_list(crippen_o3a(
        prb_mol,
        ref_mol,
        &params.prb_crippen_contribs,
        &params.ref_crippen_contribs,
        true,
        &options,
    )?)
}
