// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.TorsionFingerprints`, a Python module: the
//! torsion fingerprint deviation (TFD) between conformers (Schulz-Gasch et
//! al., JCIM 52, 1499-1512 (2012)).
//!
//! A torsion list holds, for each rotatable bond or ring, the atom
//! quartets whose dihedrals describe it and the deviation that normalizes
//! it; torsion angles hold the quartets' angles in degrees in its place.
//! The sums, means, and standard deviations are Python's and NumPy's, so
//! the values match Python's to the last bit.

use crate::chem::rd_fingerprint_generator::{
    AdditionalOutput, GetFingerprintParams, GetMorganGeneratorParams,
    get_morgan_generator_with_params,
};
use crate::chem::rd_mol_descriptors::get_connectivity_invariants;
use crate::chem::rdchem::{AtomRef, Conformer, Mol};
use crate::chem::{get_distance_matrix, get_symm_sssr, mol_from_smarts, mol_to_smiles};
use crate::geometry::rd_geometry::compute_signed_dihedral_angle;
use crate::py::{np_std, py_sum};
use crate::{Error, Result};

/// One torsion of a torsion list: the atom quartets that describe a bond or
/// ring, and the deviation that normalizes it.
/// A Rust addition, naming Python's tuple.
pub type Torsion = (Vec<[u32; 4]>, f64);

/// One torsion's angles in degrees (one for a ring, the mean of its
/// quartets' absolute angles), and the deviation that normalizes it.
/// A Rust addition, naming Python's tuple.
pub type TorsionAngles = (Vec<f64>, f64);

/// How [`calculate_torsion_lists_with_params`] normalizes the torsions of
/// bonds outside rings. Mirrors the `maxDev` argument's strings.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MaxDev {
    /// Every torsion by 180 degrees (`'equal'`).
    #[default]
    Equal,
    /// Each torsion by its own maximal deviation, as the paper gives it
    /// (`'spec'`).
    Spec,
}

fn value_error(message: &str) -> Error {
    Error::Value {
        message: message.to_owned(),
    }
}

/// Whether all the atoms have the same invariant.
fn do_match(inv: &[u64], atoms: &[AtomRef<'_>]) -> bool {
    for i in 0..atoms.len().saturating_sub(1) {
        for j in i + 1..atoms.len() {
            if inv[atoms[i].get_idx() as usize] != inv[atoms[j].get_idx() as usize] {
                return false;
            }
        }
    }
    true
}

/// Whether no two of the atoms have the same invariant.
fn do_not_match(inv: &[u64], atoms: &[AtomRef<'_>]) -> bool {
    for i in 0..atoms.len().saturating_sub(1) {
        for j in i + 1..atoms.len() {
            if inv[atoms[i].get_idx() as usize] == inv[atoms[j].get_idx() as usize] {
                return false;
            }
        }
    }
    true
}

/// Of three atoms, the one whose invariant the other two share without it.
fn do_match_except1<'m>(inv: &[u64], atoms: &[AtomRef<'m>]) -> Result<Option<AtomRef<'m>>> {
    if atoms.len() != 3 {
        return Err(value_error("Number of atoms must be three"));
    }
    let a1 = inv[atoms[0].get_idx() as usize];
    let a2 = inv[atoms[1].get_idx() as usize];
    let a3 = inv[atoms[2].get_idx() as usize];
    Ok(if a1 == a2 && a1 != a3 && a2 != a3 {
        Some(atoms[2])
    } else if a1 != a2 && a1 == a3 && a2 != a3 {
        Some(atoms[1])
    } else if a1 != a2 && a1 != a3 && a2 == a3 {
        Some(atoms[0])
    } else {
        None
    })
}

/// Each atom's invariant: the id of its Morgan environment of `radius`.
fn atom_invariants_with_radius(mol: &Mol, radius: u32) -> Result<Vec<u64>> {
    let generator = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius,
        include_redundant_environments: true,
        ..Default::default()
    })?;
    let mut ao = AdditionalOutput::default();
    ao.allocate_bit_info_map();
    let params = GetFingerprintParams {
        additional_output: Some(ao),
        ..Default::default()
    };
    let (_, ao) = generator.get_sparse_count_fingerprint_with_outputs(mol, &params)?;
    let mut inv = vec![0; mol.get_num_atoms() as usize];
    for (&k, vs) in ao.get_bit_info_map().into_iter().flatten() {
        for &(aid, r) in vs {
            if r == radius {
                inv[aid as usize] = k;
            }
        }
    }
    Ok(inv)
}

/// The atom's neighbors other than hydrogens and `exclude`.
fn heavy_atom_neighbors<'m>(atom: &AtomRef<'m>, exclude: Option<u32>) -> Vec<AtomRef<'m>> {
    atom.get_neighbors()
        .into_iter()
        .filter(|n| n.get_atomic_num() != 1 && Some(n.get_idx()) != exclude)
        .collect()
}

/// The neighbors whose torsions describe a bond's end.
fn index_for_torsion<'m>(neighbors: &[AtomRef<'m>], inv: &[u64]) -> Result<Vec<AtomRef<'m>>> {
    let first_by_invariant = || {
        let mut sorted = neighbors.to_vec();
        sorted.sort_by_key(|n| inv[n.get_idx() as usize]);
        vec![sorted[0]]
    };
    Ok(if neighbors.len() == 1 {
        vec![neighbors[0]]
    } else if do_match(inv, neighbors) {
        neighbors.to_vec()
    } else if do_not_match(inv, neighbors) || neighbors.len() != 3 {
        first_by_invariant()
    } else {
        match do_match_except1(inv, neighbors)? {
            Some(at) => vec![at],
            None => {
                return Err(value_error(
                    "Atom neighbors are either all the same or all different",
                ));
            }
        }
    })
}

/// A bond (or a pair of atoms standing for one) whose torsions count, with
/// its ends' heavy neighbors.
struct TorsionBond<'m> {
    a1: u32,
    a2: u32,
    nb1: Vec<AtomRef<'m>>,
    nb2: Vec<AtomRef<'m>>,
}

/// The bonds outside rings, between atoms with heavy neighbors, whose
/// torsions count, as Python's `_getBondsForTorsions` finds them.
fn bonds_for_torsions(mol: &Mol, ignore_colinear_bonds: bool) -> Result<Vec<TorsionBond<'_>>> {
    // Triple bonds and allenes cannot hold a torsion's middle atoms.
    let mut atom_flags = vec![false; mol.get_num_atoms() as usize];
    for smarts in ["*#*", "[$([C](=*)=*)]"] {
        let p = mol_from_smarts(smarts).expect("the torsion SMARTS parse");
        for m in mol.get_substruct_matches(&p)? {
            for a in m {
                atom_flags[a as usize] = true;
            }
        }
    }
    let heavy = |idx: u32, exclude: u32| -> Result<Vec<AtomRef<'_>>> {
        Ok(heavy_atom_neighbors(
            &mol.get_atom_with_idx(idx)?,
            Some(exclude),
        ))
    };
    let mut bonds = Vec::new();
    let mut done_bonds = vec![false; mol.get_num_bonds() as usize];
    for bond in mol.get_bonds() {
        if bond.is_in_ring() {
            continue;
        }
        let (mut a1, mut a2) = (bond.get_begin_atom_idx(), bond.get_end_atom_idx());
        let mut nb1 = heavy(a1, a2)?;
        let mut nb2 = heavy(a2, a1)?;
        if done_bonds[bond.get_idx() as usize] || nb1.is_empty() || nb2.is_empty() {
            continue;
        }
        done_bonds[bond.get_idx() as usize] = true;
        if !atom_flags[a1 as usize] && !atom_flags[a2 as usize] {
            bonds.push(TorsionBond { a1, a2, nb1, nb2 });
            continue;
        }
        if ignore_colinear_bonds {
            continue;
        }
        // Look past the linear atoms for others that are not bonded. Python
        // marks the bond it last followed, the first loop's or this one's,
        // done in the second loop, as this does.
        let mut b = bond.get_idx();
        while nb1.len() == 1 && atom_flags[a1 as usize] {
            let a1old = a1;
            a1 = nb1[0].get_idx();
            let next = mol
                .get_bond_between_atoms(a1old, a1)?
                .expect("neighbors are bonded");
            nb1 = heavy(a1, a1old)?;
            b = next.get_idx();
            done_bonds[b as usize] = true;
        }
        while nb2.len() == 1 && atom_flags[a2 as usize] {
            done_bonds[b as usize] = true;
            let a2old = a2;
            a2 = nb2[0].get_idx();
            let next = mol
                .get_bond_between_atoms(a2old, a2)?
                .expect("neighbors are bonded");
            nb2 = heavy(a2, a2old)?;
            b = next.get_idx();
            done_bonds[b as usize] = true;
        }
        if !nb1.is_empty() && !nb2.is_empty() {
            bonds.push(TorsionBond { a1, a2, nb1, nb2 });
        }
    }
    Ok(bonds)
}

/// Parameters for [`calculate_torsion_lists_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CalculateTorsionListsParams {
    /// How bonds' torsions are normalized.
    pub max_dev: MaxDev,
    /// The radius of the Morgan environments that tell symmetric atoms
    /// apart; 0 uses connectivity invariants.
    pub symm_radius: u32,
    /// Ignores single bonds next to triple bonds and allenes; otherwise
    /// atoms past the linear ones stand in for the bond's.
    pub ignore_colinear_bonds: bool,
}

impl Default for CalculateTorsionListsParams {
    fn default() -> CalculateTorsionListsParams {
        CalculateTorsionListsParams {
            max_dev: MaxDev::Equal,
            symm_radius: 2,
            ignore_colinear_bonds: true,
        }
    }
}

/// The molecule's torsions: those of bonds outside rings, and those of its
/// rings (symmetrized SSSR), each ring's quartets in ring order. Mirrors
/// `TorsionFingerprints.CalculateTorsionLists`.
pub fn calculate_torsion_lists(mol: &Mol) -> Result<(Vec<Torsion>, Vec<Torsion>)> {
    calculate_torsion_lists_with_params(mol, &CalculateTorsionListsParams::default())
}

/// Mirrors `TorsionFingerprints.CalculateTorsionLists(mol, maxDev,
/// symmRadius, ignoreColinearBonds)`.
pub fn calculate_torsion_lists_with_params(
    mol: &Mol,
    params: &CalculateTorsionListsParams,
) -> Result<(Vec<Torsion>, Vec<Torsion>)> {
    let bonds = bonds_for_torsions(mol, params.ignore_colinear_bonds)?;
    let inv = if params.symm_radius > 0 {
        atom_invariants_with_radius(mol, params.symm_radius)?
    } else {
        get_connectivity_invariants(mol)?
            .into_iter()
            .map(u64::from)
            .collect()
    };
    let mut tors_list = Vec::new();
    for TorsionBond { a1, a2, nb1, nb2 } in &bonds {
        let (a1, a2) = (*a1, *a2);
        let d1 = index_for_torsion(nb1, &inv)?;
        let d2 = index_for_torsion(nb2, &inv)?;
        let quartet = |n1: &AtomRef<'_>, n2: &AtomRef<'_>| [n1.get_idx(), a1, a2, n2.get_idx()];
        let (quartets, max_dev) = if d1.len() == 1 && d2.len() == 1 {
            (vec![quartet(&d1[0], &d2[0])], 180.0)
        } else if d1.len() == 1 {
            let quartets = d2.iter().map(|nb| quartet(&d1[0], nb)).collect();
            (quartets, if nb2.len() == 2 { 90.0 } else { 60.0 })
        } else if d2.len() == 1 {
            let quartets = d1.iter().map(|nb| quartet(nb, &d2[0])).collect();
            (quartets, if nb1.len() == 2 { 90.0 } else { 60.0 })
        } else {
            let mut quartets = Vec::new();
            for n1 in &d1 {
                for n2 in &d2 {
                    quartets.push(quartet(n1, n2));
                }
            }
            let max_dev = if nb1.len() == 2 && nb2.len() == 2 {
                90.0
            } else if nb1.len() == 3 && nb2.len() == 3 {
                60.0
            } else {
                30.0
            };
            (quartets, max_dev)
        };
        let max_dev = match params.max_dev {
            MaxDev::Equal => 180.0,
            MaxDev::Spec => max_dev,
        };
        tors_list.push((quartets, max_dev));
    }
    let mut copy = mol.clone();
    let mut tors_list_rings = Vec::new();
    for r in get_symm_sssr(&mut copy)? {
        let num = r.len();
        let max_dev = if num >= 14 {
            180.0
        } else {
            let d = num as f64 - 14.0;
            180.0 * (-0.025 * d * d).exp()
        };
        let quartets = (0..num)
            .map(|i| [r[i], r[(i + 1) % num], r[(i + 2) % num], r[(i + 3) % num]])
            .collect();
        tors_list_rings.push((quartets, max_dev));
    }
    Ok((tors_list, tors_list_rings))
}

/// The angles of the torsions in the default conformer: each quartet's
/// signed dihedral from 0 to 360 degrees for bonds, the mean absolute
/// dihedral for rings. Mirrors `TorsionFingerprints.CalculateTorsionAngles`.
pub fn calculate_torsion_angles(
    mol: &Mol,
    tors_list: &[Torsion],
    tors_list_rings: &[Torsion],
) -> Result<Vec<TorsionAngles>> {
    calculate_torsion_angles_with_conf_id(mol, tors_list, tors_list_rings, -1)
}

/// In the conformer with `conf_id`. Mirrors
/// `TorsionFingerprints.CalculateTorsionAngles(mol, tors_list,
/// tors_list_rings, confId)`.
pub fn calculate_torsion_angles_with_conf_id(
    mol: &Mol,
    tors_list: &[Torsion],
    tors_list_rings: &[Torsion],
    conf_id: i32,
) -> Result<Vec<TorsionAngles>> {
    let conf = mol.get_conformer_with_id(conf_id)?;
    let dihedral = |atoms: &[u32; 4]| -> Result<f64> {
        let p1 = conf.get_atom_position(atoms[0])?;
        let p2 = conf.get_atom_position(atoms[1])?;
        let p3 = conf.get_atom_position(atoms[2])?;
        let p4 = conf.get_atom_position(atoms[3])?;
        Ok(compute_signed_dihedral_angle(&p1, &p2, &p3, &p4) / std::f64::consts::PI * 180.0)
    };
    let mut torsions = Vec::new();
    for (quartets, max_dev) in tors_list {
        let mut tors = Vec::new();
        for atoms in quartets {
            let mut t = dihedral(atoms)?;
            if t < 0.0 {
                t += 360.0;
            }
            tors.push(t);
        }
        torsions.push((tors, *max_dev));
    }
    for (quartets, max_dev) in tors_list_rings {
        // Python adds from the integer 0, which is the same as from 0.0.
        let mut tors = 0.0;
        for atoms in quartets {
            tors += dihedral(atoms)?.abs();
        }
        tors /= quartets.len() as f64;
        torsions.push((vec![tors], *max_dev));
    }
    Ok(torsions)
}

/// The two most central atoms, bonded, the most central first: those whose
/// topological distances to the others vary least.
fn find_central_bond(mol: &Mol, distmat: &[Vec<f64>]) -> Result<(u32, u32)> {
    let mut stds = Vec::new();
    for atom in mol.get_atoms() {
        if heavy_atom_neighbors(&atom, None).len() < 2 {
            continue;
        }
        let i = atom.get_idx() as usize;
        let mut row = distmat[i].clone();
        row.remove(i);
        stds.push((np_std(&row), atom.get_idx()));
    }
    // Python sorts the (std, index) tuples; the deviations are never NaN.
    stds.sort_by(|a, b| a.partial_cmp(b).expect("the deviations are finite"));
    let len = stds.len();
    let aid1 = stds.first().ok_or(Error::Index { index: 0, len })?.1;
    let mut i = 1;
    loop {
        let &(_, aid2) = stds.get(i).ok_or(Error::index(i, len))?;
        if mol.get_bond_between_atoms(aid1, aid2)?.is_some() {
            return Ok((aid1, aid2));
        }
        i += 1;
    }
}

/// The weights' beta, which gives half the largest distance from `aid1` to
/// a bond whose end atom has heavy neighbors a weight of 0.1.
///
/// RDKit 2024.03.6's Python means to take the bonds between atoms with heavy
/// neighbors, but tests `len(nb2) > 1` twice, so only the end atom counts;
/// this computes what it does.
fn calculate_beta(mol: &Mol, distmat: &[Vec<f64>], aid1: u32) -> Result<f64> {
    let mut dmax = 0.0;
    for b in mol.get_bonds() {
        let begin = b.get_begin_atom();
        let end = b.get_end_atom();
        if heavy_atom_neighbors(&end, None).len() > 1 {
            let row = &distmat[aid1 as usize];
            let d = f64::max(row[begin.get_idx() as usize], row[end.get_idx() as usize]);
            if d > dmax {
                dmax = d;
            }
        }
    }
    let dmax2 = dmax / 2.0;
    if dmax2 == 0.0 {
        // Python's float division raises ZeroDivisionError.
        return Err(value_error("float division by zero"));
    }
    Ok(-(0.1f64.ln()) / (dmax2 * dmax2))
}

/// Parameters for [`calculate_torsion_weights_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CalculateTorsionWeightsParams {
    /// The first atom of the central bond; `None`, with `aid2`, for the
    /// most central atom.
    pub aid1: Option<u32>,
    /// The second atom of the central bond; `None`, with `aid1`, for the
    /// most central atom bonded to the first.
    pub aid2: Option<u32>,
    /// As in [`CalculateTorsionListsParams`].
    pub ignore_colinear_bonds: bool,
}

impl Default for CalculateTorsionWeightsParams {
    fn default() -> CalculateTorsionWeightsParams {
        CalculateTorsionWeightsParams {
            aid1: None,
            aid2: None,
            ignore_colinear_bonds: true,
        }
    }
}

/// The weights of the torsions, bonds' then rings', highest for the bond
/// between the two most central atoms and falling off with distance from
/// it. `Error::Index` when the molecule has no two bonded central atoms,
/// as Python's `IndexError`. Mirrors
/// `TorsionFingerprints.CalculateTorsionWeights`.
pub fn calculate_torsion_weights(mol: &Mol) -> Result<Vec<f64>> {
    calculate_torsion_weights_with_params(mol, &CalculateTorsionWeightsParams::default())
}

/// From the bond between `aid1` and `aid2`, which must be bonded
/// (`Error::Value` otherwise; giving one of them is the same). Mirrors
/// `TorsionFingerprints.CalculateTorsionWeights(mol, aid1, aid2,
/// ignoreColinearBonds)`.
pub fn calculate_torsion_weights_with_params(
    mol: &Mol,
    params: &CalculateTorsionWeightsParams,
) -> Result<Vec<f64>> {
    let distmat = get_distance_matrix(mol)?;
    let (aid1, aid2) = match (params.aid1, params.aid2) {
        (None, None) => find_central_bond(mol, &distmat)?,
        (Some(aid1), Some(aid2)) if mol.get_bond_between_atoms(aid1, aid2)?.is_some() => {
            (aid1, aid2)
        }
        _ => return Err(value_error("Specified atoms must be connected by a bond.")),
    };
    let beta = calculate_beta(mol, &distmat, aid1)?;
    let dist = |bid1: u32, bid2: u32| -> f64 {
        let (r1, r2) = (&distmat[aid1 as usize], &distmat[aid2 as usize]);
        let (b1, b2) = (bid1 as usize, bid2 as usize);
        r1[b1].min(r1[b2]).min(r2[b1]).min(r2[b2]) + 1.0
    };
    let mut weights = Vec::new();
    for TorsionBond { a1, a2, .. } in bonds_for_torsions(mol, params.ignore_colinear_bonds)? {
        let d = if (a1, a2) == (aid1, aid2) || (a2, a1) == (aid1, aid2) {
            0.0
        } else {
            dist(a1, a2)
        };
        weights.push((-beta * (d * d)).exp());
    }
    for r in mol.get_ring_info().bond_rings() {
        let mut tmp = Vec::new();
        for &bidx in &r {
            let b = mol.get_bond_with_idx(bidx)?;
            tmp.push(dist(b.get_begin_atom_idx(), b.get_end_atom_idx()));
        }
        let num = r.len() as f64;
        let w = py_sum(tmp) / num;
        let w = (-beta * (w * w)).exp();
        weights.push(w * (num / 2.0));
    }
    Ok(weights)
}

/// The torsion fingerprint deviation between two conformers' torsion
/// angles: the mean of each torsion's smallest angle difference over its
/// deviation. Lists of different lengths are `Error::Value`. Mirrors
/// `TorsionFingerprints.CalculateTFD`.
pub fn calculate_tfd(torsions1: &[TorsionAngles], torsions2: &[TorsionAngles]) -> Result<f64> {
    calculate_tfd_with_weights(torsions1, torsions2, None)
}

/// The weighted mean, with one weight per torsion (`Error::Value`
/// otherwise). Mirrors `TorsionFingerprints.CalculateTFD(torsions1,
/// torsions2, weights)`.
pub fn calculate_tfd_with_weights(
    torsions1: &[TorsionAngles],
    torsions2: &[TorsionAngles],
    weights: Option<&[f64]>,
) -> Result<f64> {
    if torsions1.len() != torsions2.len() {
        return Err(value_error(
            "List of torsions angles must have the same size.",
        ));
    }
    let mut deviations = Vec::with_capacity(torsions1.len());
    for (tors1, tors2) in torsions1.iter().zip(torsions2) {
        let mut mindiff = 180.0f64;
        for &t1 in &tors1.0 {
            for &t2 in &tors2.0 {
                let mut diff = (t1 - t2).abs();
                if 360.0 - diff < diff {
                    diff = 360.0 - diff;
                }
                if diff < mindiff {
                    mindiff = diff;
                }
            }
        }
        deviations.push(mindiff / tors1.1);
    }
    let sum_weights = match weights {
        Some(weights) => {
            if weights.len() != torsions1.len() {
                return Err(value_error(
                    "List of torsions angles and weights must have the same size.",
                ));
            }
            for (d, w) in deviations.iter_mut().zip(weights) {
                *d *= w;
            }
            py_sum(weights.iter().copied())
        }
        None => deviations.len() as f64,
    };
    let mut tfd = py_sum(deviations);
    if sum_weights != 0.0 {
        tfd /= sum_weights;
    }
    Ok(tfd)
}

/// `mol2`, with `mol1`'s atom order when they differ: a copy of `mol1` with
/// `mol2`'s conformers, their positions reordered by `mol2`'s match of
/// `mol1`.
fn same_atom_order(mol1: &Mol, mol2: &Mol) -> Result<Mol> {
    let m = mol2.get_substruct_match(mol1)?;
    if m.iter().copied().eq(0..mol1.get_num_atoms()) {
        return Ok(mol2.clone());
    }
    let mut mol3 = mol1.clone();
    mol3.remove_all_conformers();
    for conf2 in mol2.get_conformers() {
        let mut conf = Conformer::new(mol1.get_num_atoms());
        conf.set_id(conf2.get_id())?;
        for i in 0..mol1.get_num_atoms() {
            let j = *m.get(i as usize).ok_or(Error::index(i, m.len()))?;
            conf.set_atom_position(i, conf2.get_atom_position(j)?)?;
        }
        mol3.add_conformer(&conf)?;
    }
    Ok(mol3)
}

fn check_same_molecule(mol1: &Mol, mol2: &Mol) -> Result<()> {
    if mol_to_smiles(mol1)? != mol_to_smiles(mol2)? {
        return Err(value_error(
            "The two molecules must be instances of the same molecule!",
        ));
    }
    Ok(())
}

fn torsion_weights(
    mol: &Mol,
    use_weights: bool,
    ignore_colinear_bonds: bool,
) -> Result<Option<Vec<f64>>> {
    if !use_weights {
        return Ok(None);
    }
    let params = CalculateTorsionWeightsParams {
        ignore_colinear_bonds,
        ..Default::default()
    };
    calculate_torsion_weights_with_params(mol, &params).map(Some)
}

macro_rules! tfd_params {
    ($(#[$meta:meta])* $name:ident { $($(#[$fmeta:meta])* $field:ident: $ty:ty = $default:expr,)* }) => {
        $(#[$meta])*
        #[derive(Clone, Debug, PartialEq)]
        pub struct $name {
            $($(#[$fmeta])* pub $field: $ty,)*
            /// Weights the torsions by their distance from the central bond.
            pub use_weights: bool,
            /// How bonds' torsions are normalized.
            pub max_dev: MaxDev,
            /// As in [`CalculateTorsionListsParams`].
            pub symm_radius: u32,
            /// As in [`CalculateTorsionListsParams`].
            pub ignore_colinear_bonds: bool,
        }

        impl Default for $name {
            fn default() -> $name {
                $name {
                    $($field: $default,)*
                    use_weights: true,
                    max_dev: MaxDev::Equal,
                    symm_radius: 2,
                    ignore_colinear_bonds: true,
                }
            }
        }

        impl $name {
            fn lists_params(&self) -> CalculateTorsionListsParams {
                CalculateTorsionListsParams {
                    max_dev: self.max_dev,
                    symm_radius: self.symm_radius,
                    ignore_colinear_bonds: self.ignore_colinear_bonds,
                }
            }
        }
    };
}

tfd_params! {
    /// Parameters for [`get_tfd_between_conformers_with_params`], with
    /// Python's defaults (A2).
    GetTFDBetweenConformersParams {}
}

/// The TFD between each conformer of `conf_ids1` and each of `conf_ids2`,
/// row by row. Mirrors `TorsionFingerprints.GetTFDBetweenConformers`.
pub fn get_tfd_between_conformers(
    mol: &Mol,
    conf_ids1: &[i32],
    conf_ids2: &[i32],
) -> Result<Vec<f64>> {
    get_tfd_between_conformers_with_params(
        mol,
        conf_ids1,
        conf_ids2,
        &GetTFDBetweenConformersParams::default(),
    )
}

/// Mirrors `TorsionFingerprints.GetTFDBetweenConformers(mol, confIds1,
/// confIds2, useWeights, maxDev, symmRadius, ignoreColinearBonds)`.
pub fn get_tfd_between_conformers_with_params(
    mol: &Mol,
    conf_ids1: &[i32],
    conf_ids2: &[i32],
    params: &GetTFDBetweenConformersParams,
) -> Result<Vec<f64>> {
    let (tl, tlr) = calculate_torsion_lists_with_params(mol, &params.lists_params())?;
    let angles = |ids: &[i32]| -> Result<Vec<_>> {
        ids.iter()
            .map(|&cid| calculate_torsion_angles_with_conf_id(mol, &tl, &tlr, cid))
            .collect()
    };
    let torsions1 = angles(conf_ids1)?;
    let torsions2 = angles(conf_ids2)?;
    let weights = torsion_weights(mol, params.use_weights, params.ignore_colinear_bonds)?;
    let mut tfd = Vec::new();
    for t1 in &torsions1 {
        for t2 in &torsions2 {
            tfd.push(calculate_tfd_with_weights(t1, t2, weights.as_deref())?);
        }
    }
    Ok(tfd)
}

tfd_params! {
    /// Parameters for [`get_tfd_between_molecules_with_params`], with
    /// Python's defaults (A2).
    GetTFDBetweenMoleculesParams {
        /// The conformer of `mol1`; −1 for the default.
        conf_id1: i32 = -1,
        /// The conformer of `mol2`; −1 for the default.
        conf_id2: i32 = -1,
    }
}

/// The TFD between two instances of one molecule, which must have the same
/// SMILES (`Error::Value` otherwise), in either atom order. Mirrors
/// `TorsionFingerprints.GetTFDBetweenMolecules`.
pub fn get_tfd_between_molecules(mol1: &Mol, mol2: &Mol) -> Result<f64> {
    get_tfd_between_molecules_with_params(mol1, mol2, &GetTFDBetweenMoleculesParams::default())
}

/// Mirrors `TorsionFingerprints.GetTFDBetweenMolecules(mol1, mol2, confId1,
/// confId2, useWeights, maxDev, symmRadius, ignoreColinearBonds)`.
pub fn get_tfd_between_molecules_with_params(
    mol1: &Mol,
    mol2: &Mol,
    params: &GetTFDBetweenMoleculesParams,
) -> Result<f64> {
    check_same_molecule(mol1, mol2)?;
    let mol2 = same_atom_order(mol1, mol2)?;
    let (tl, tlr) = calculate_torsion_lists_with_params(mol1, &params.lists_params())?;
    let torsion1 = calculate_torsion_angles_with_conf_id(mol1, &tl, &tlr, params.conf_id1)?;
    let torsion2 = calculate_torsion_angles_with_conf_id(&mol2, &tl, &tlr, params.conf_id2)?;
    let weights = torsion_weights(mol1, params.use_weights, params.ignore_colinear_bonds)?;
    calculate_tfd_with_weights(&torsion1, &torsion2, weights.as_deref())
}

tfd_params! {
    /// Parameters for [`get_best_tfd_between_molecules_with_params`], with
    /// Python's defaults (A2).
    GetBestTFDBetweenMoleculesParams {
        /// The conformer of `mol1`; −1 for the default.
        conf_id1: i32 = -1,
    }
}

/// The smallest TFD between a conformer of `mol1` and each of `mol2`'s, or
/// 1e8 when `mol2` has none; the molecules as in
/// [`get_tfd_between_molecules`]. Mirrors
/// `TorsionFingerprints.GetBestTFDBetweenMolecules`.
pub fn get_best_tfd_between_molecules(mol1: &Mol, mol2: &Mol) -> Result<f64> {
    get_best_tfd_between_molecules_with_params(
        mol1,
        mol2,
        &GetBestTFDBetweenMoleculesParams::default(),
    )
}

/// Mirrors `TorsionFingerprints.GetBestTFDBetweenMolecules(mol1, mol2,
/// confId1, useWeights, maxDev, symmRadius, ignoreColinearBonds)`.
pub fn get_best_tfd_between_molecules_with_params(
    mol1: &Mol,
    mol2: &Mol,
    params: &GetBestTFDBetweenMoleculesParams,
) -> Result<f64> {
    check_same_molecule(mol1, mol2)?;
    let mol2 = same_atom_order(mol1, mol2)?;
    let (tl, tlr) = calculate_torsion_lists_with_params(mol1, &params.lists_params())?;
    let torsion1 = calculate_torsion_angles_with_conf_id(mol1, &tl, &tlr, params.conf_id1)?;
    let weights = torsion_weights(mol1, params.use_weights, params.ignore_colinear_bonds)?;
    let mut best = 1e8f64;
    for conf in mol2.get_conformers() {
        let torsion2 =
            calculate_torsion_angles_with_conf_id(&mol2, &tl, &tlr, conf.get_id() as i32)?;
        let tfd = calculate_tfd_with_weights(&torsion1, &torsion2, weights.as_deref())?;
        best = best.min(tfd);
    }
    Ok(best)
}

tfd_params! {
    /// Parameters for [`get_tfd_matrix_with_params`], with Python's defaults
    /// (A2).
    GetTFDMatrixParams {}
}

/// The TFD between each pair of the molecule's conformers, as the lower
/// triangle of the matrix, row by row: for conformers 0 to 3, (1, 0), (2,
/// 0), (2, 1), (3, 0), and so on. Mirrors `TorsionFingerprints.GetTFDMatrix`.
pub fn get_tfd_matrix(mol: &Mol) -> Result<Vec<f64>> {
    get_tfd_matrix_with_params(mol, &GetTFDMatrixParams::default())
}

/// Mirrors `TorsionFingerprints.GetTFDMatrix(mol, useWeights, maxDev,
/// symmRadius, ignoreColinearBonds)`.
pub fn get_tfd_matrix_with_params(mol: &Mol, params: &GetTFDMatrixParams) -> Result<Vec<f64>> {
    let (tl, tlr) = calculate_torsion_lists_with_params(mol, &params.lists_params())?;
    let torsions = mol
        .get_conformers()
        .map(|conf| calculate_torsion_angles_with_conf_id(mol, &tl, &tlr, conf.get_id() as i32))
        .collect::<Result<Vec<_>>>()?;
    let weights = torsion_weights(mol, params.use_weights, params.ignore_colinear_bonds)?;
    let mut tfdmat = Vec::new();
    for i in 0..torsions.len() {
        for j in 0..i {
            tfdmat.push(calculate_tfd_with_weights(
                &torsions[i],
                &torsions[j],
                weights.as_deref(),
            )?);
        }
    }
    Ok(tfdmat)
}
