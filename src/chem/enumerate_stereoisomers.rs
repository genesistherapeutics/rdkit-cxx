// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.EnumerateStereoisomers`, a Python module: the
//! stereoisomers of a molecule, from flipping its unassigned stereocenters,
//! double bonds, and stereo groups.
//!
//! The enumeration is an iterator, made as it is asked for, as Python's
//! generator is. Where Python samples isomers at random it draws from a
//! [`Random`], which draws what Python's `random.Random` would from the same
//! seed; without one it seeds from the molecule, with the hash Python
//! computes.

use std::collections::HashSet;

use crate::chem::rd_dist_geom::{EmbedParameters, embed_molecule_with_params};
use crate::chem::{
    AssignStereochemistryParams, BondDir, BondStereo, ChiralType, Conformer, Mol, RWMol,
    SmilesWriteParams, StereoGroupType, StereoInfo, StereoSpecified, StereoType, add_hs,
    assign_stereochemistry_with_params, find_potential_stereo, mol_to_smiles_with_params,
    set_double_bond_neighbor_directions,
};
use crate::rd_random::{Random, hash_tuple};
use crate::{Error, Result};

/// How stereoisomers are enumerated. Mirrors
/// `EnumerateStereoisomers.StereoEnumerationOptions`; `Default` mirrors
/// `StereoEnumerationOptions()`.
#[derive(Clone, Debug)]
pub struct StereoEnumerationOptions {
    /// Keeps only the isomers RDKit's distance geometry can embed, each
    /// with the conformer it found. Expensive, and only a heuristic.
    pub try_embedding: bool,
    /// Leaves specified stereocenters as they are, unless they are in a
    /// relative stereo group.
    pub only_unassigned: bool,
    /// The most isomers to yield, 0 for all; when there are more, a random
    /// sample of them.
    pub max_isomers: u32,
    /// The generator that samples isomers; `None` for one seeded from the
    /// molecule. Python also takes an `int` seed, here `Random::new(seed)`.
    /// Each enumeration draws from its own copy.
    pub rand: Option<Random>,
    /// Yields each isomer, by its canonical SMILES, once.
    pub unique: bool,
    /// Only flips the molecule's stereo groups.
    pub only_stereo_groups: bool,
    /// Prints the SMILES of each isomer that does not embed, as
    /// `EnumerateStereoisomers`' `verbose` argument does (A3).
    pub verbose: bool,
}

impl Default for StereoEnumerationOptions {
    fn default() -> StereoEnumerationOptions {
        StereoEnumerationOptions {
            try_embedding: false,
            only_unassigned: true,
            max_isomers: 1024,
            rand: None,
            unique: true,
            only_stereo_groups: false,
            verbose: false,
        }
    }
}

/// What flips one stereo element of the molecule being enumerated.
enum Flipper {
    /// Sets the atom clockwise or counterclockwise.
    Atom(u32),
    /// Sets the bond cis or trans.
    Bond(u32),
    /// Keeps the group's atoms' original tags or inverts them all.
    StereoGroup(Vec<(u32, ChiralType)>),
}

impl Flipper {
    fn flip(&self, mol: &mut Mol, flag: bool) -> Result<()> {
        match self {
            Flipper::Atom(idx) => mol.get_atom_with_idx_mut(*idx)?.set_chiral_tag(if flag {
                ChiralType::ChiTetrahedralCw
            } else {
                ChiralType::ChiTetrahedralCcw
            }),
            Flipper::Bond(idx) => mol.get_bond_with_idx_mut(*idx)?.set_stereo(if flag {
                BondStereo::Stereocis
            } else {
                BondStereo::Stereotrans
            }),
            Flipper::StereoGroup(original) => {
                for &(a, parity) in original {
                    let tag = match (flag, parity) {
                        (true, _) => parity,
                        (false, ChiralType::ChiTetrahedralCw) => ChiralType::ChiTetrahedralCcw,
                        (false, ChiralType::ChiTetrahedralCcw) => ChiralType::ChiTetrahedralCw,
                        _ => continue,
                    };
                    mol.get_atom_with_idx_mut(a)?.set_chiral_tag(tag);
                }
            }
        }
        Ok(())
    }
}

fn get_flippers(mol: &mut Mol, options: &StereoEnumerationOptions) -> Result<Vec<Flipper>> {
    let sinfo: Vec<StereoInfo> = find_potential_stereo(mol)?;
    let mut flippers = Vec::new();
    if !options.only_stereo_groups {
        for si in &sinfo {
            if options.only_unassigned
                && !matches!(
                    si.specified,
                    StereoSpecified::Unspecified | StereoSpecified::Unknown
                )
            {
                continue;
            }
            match si.type_ {
                StereoType::AtomTetrahedral => flippers.push(Flipper::Atom(si.centered_on)),
                StereoType::BondDouble => {
                    let mut bnd = mol.get_bond_with_idx_mut(si.centered_on)?;
                    if bnd.get_stereo_atoms().is_empty() {
                        let ca = &si.controlling_atoms;
                        if ca.len() < 3
                            || ca[0] == StereoInfo::NOATOM
                            || ca[2] == StereoInfo::NOATOM
                        {
                            continue;
                        }
                        bnd.set_stereo_atoms(ca[0], ca[2])?;
                    }
                    flippers.push(Flipper::Bond(si.centered_on));
                }
                _ => {}
            }
        }
    }
    if options.only_unassigned {
        for group in mol.get_stereo_groups() {
            if group.get_group_type() != StereoGroupType::StereoAbsolute {
                let original = group
                    .get_atoms()
                    .iter()
                    .map(|&a| Ok((a, mol.get_atom_with_idx(a)?.get_chiral_tag())))
                    .collect::<Result<_>>()?;
                flippers.push(Flipper::StereoGroup(original));
            }
        }
    }
    Ok(flippers)
}

/// An estimate, an upper bound, of the molecule's number of stereoisomers:
/// 2 to the number of elements the options would flip. More than 63 such
/// elements is `Error::Value`, where Python's `int` grows. Mirrors
/// `EnumerateStereoisomers.GetStereoisomerCount(m)`.
pub fn get_stereoisomer_count(m: &Mol) -> Result<u64> {
    get_stereoisomer_count_with_options(m, &StereoEnumerationOptions::default())
}

/// As [`get_stereoisomer_count`]. Mirrors
/// `EnumerateStereoisomers.GetStereoisomerCount(m, options)`.
pub fn get_stereoisomer_count_with_options(
    m: &Mol,
    options: &StereoEnumerationOptions,
) -> Result<u64> {
    let mut tm = m.clone();
    let n = get_flippers(&mut tm, options)?.len();
    1u64.checked_shl(n as u32)
        .filter(|_| n < 64)
        .ok_or_else(|| Error::Value {
            message: format!("2^{n} stereoisomers overflows u64"),
        })
}

/// The flags for each enumerated isomer: `n` bits, least significant word
/// first.
enum BitSource {
    /// Every combination in order, from 0 to `end` (2^n) when `n` is below
    /// 128, and without end otherwise.
    Range {
        next: u128,
        end: Option<u128>,
        n_words: usize,
    },
    /// Combinations at random, each once, until all have been drawn.
    UniqueRandom {
        n: u32,
        rand: Box<Random>,
        already_seen: HashSet<Vec<u32>>,
    },
}

impl BitSource {
    fn next(&mut self) -> Option<Vec<u32>> {
        match self {
            BitSource::Range { next, end, n_words } => {
                if end.is_some_and(|end| *next >= end) {
                    return None;
                }
                let bits = (0..*n_words)
                    .map(|w| if w < 4 { (*next >> (32 * w)) as u32 } else { 0 })
                    .collect();
                *next += 1;
                Some(bits)
            }
            BitSource::UniqueRandom {
                n,
                rand,
                already_seen,
            } => loop {
                if *n < 64 && already_seen.len() as u64 >= 1u64 << *n {
                    return None;
                }
                let bits = rand.getrandbits_words(*n);
                if already_seen.insert(bits.clone()) {
                    return Some(bits);
                }
            },
        }
    }
}

/// The stereoisomers of a molecule, made as they are asked for. Built by
/// [`enumerate_stereoisomers`].
/// A Rust addition: the iterator of Python's generator.
pub struct Stereoisomers {
    tm: Mol,
    flippers: Vec<Flipper>,
    options: StereoEnumerationOptions,
    /// `None` once the enumeration ends; for a molecule with nothing to
    /// flip, the molecule is yielded first.
    bitsource: Option<BitSource>,
    unflipped: Option<Mol>,
    isomers_seen: HashSet<String>,
    num_isomers: u32,
}

impl std::fmt::Debug for Stereoisomers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Stereoisomers").finish_non_exhaustive()
    }
}

impl Stereoisomers {
    fn isomer(&mut self, bitflag: &[u32]) -> Result<Option<Mol>> {
        for (i, flipper) in self.flippers.iter().enumerate() {
            let flag = bitflag
                .get(i / 32)
                .is_some_and(|w| w & (1 << (i % 32)) != 0);
            flipper.flip(&mut self.tm, flag)?;
        }
        let mut isomer = if self.tm.get_stereo_groups().is_empty() {
            self.tm.clone()
        } else {
            let mut isomer = RWMol::from(self.tm.clone());
            isomer.set_stereo_groups(&[])?;
            Mol::from(isomer)
        };
        set_double_bond_neighbor_directions(&mut isomer)?;
        isomer.clear_computed_props_with_include_rings(false);
        assign_stereochemistry_with_params(
            &mut isomer,
            &AssignStereochemistryParams {
                clean_it: true,
                force: true,
                flag_possible_stereo_centers: true,
            },
        )?;
        if self.options.unique {
            let cansmi = isomeric_smiles(&isomer)?;
            if !self.isomers_seen.insert(cansmi) {
                return Ok(None);
            }
        }
        if self.options.try_embedding {
            let mut ntm = add_hs(&isomer)?;
            let mut params = EmbedParameters {
                use_exp_torsion_angle_prefs: true,
                use_basic_knowledge: true,
                use_macrocycle_torsions: true,
                use_macrocycle14config: true,
                random_seed: (bitflag.first().copied().unwrap_or(0) & 0x7fff_ffff) as i32,
                ..EmbedParameters::default()
            };
            let cid = embed_molecule_with_params(&mut ntm, &mut params)?;
            if cid < 0 {
                if self.options.verbose {
                    println!("{}    failed to embed", isomeric_smiles(&isomer)?);
                }
                return Ok(None);
            }
            let mut conf = Conformer::new(isomer.get_num_atoms());
            let positions = ntm.get_conformer()?;
            for aid in 0..isomer.get_num_atoms() {
                conf.set_atom_position(aid, positions.get_atom_position(aid)?)?;
            }
            isomer.add_conformer(&conf)?;
        }
        Ok(Some(isomer))
    }
}

fn isomeric_smiles(mol: &Mol) -> Result<String> {
    mol_to_smiles_with_params(
        mol,
        &SmilesWriteParams {
            do_isomeric_smiles: true,
            ..Default::default()
        },
    )
}

/// The next stereoisomer.
impl Iterator for Stereoisomers {
    type Item = Result<Mol>;

    fn next(&mut self) -> Option<Result<Mol>> {
        if let Some(mol) = self.unflipped.take() {
            return Some(Ok(mol));
        }
        loop {
            let max = self.options.max_isomers;
            if max != 0 && self.num_isomers >= max {
                self.bitsource = None;
            }
            let bitflag = self.bitsource.as_mut()?.next()?;
            match self.isomer(&bitflag) {
                Ok(Some(isomer)) => {
                    self.num_isomers += 1;
                    return Some(Ok(isomer));
                }
                Ok(None) => continue,
                Err(e) => {
                    self.bitsource = None;
                    return Some(Err(e));
                }
            }
        }
    }
}

/// The molecule's stereoisomers, expanding its stereo groups: every
/// combination of its unassigned stereocenters and double bonds, at most
/// 1024 of them (a random sample of them, when there are more), each once.
/// A molecule with nothing to flip yields itself. Mirrors
/// `EnumerateStereoisomers.EnumerateStereoisomers(m)`.
pub fn enumerate_stereoisomers(m: &Mol) -> Result<Stereoisomers> {
    enumerate_stereoisomers_with_options(m, &StereoEnumerationOptions::default())
}

/// As [`enumerate_stereoisomers`]. Mirrors
/// `EnumerateStereoisomers.EnumerateStereoisomers(m, options, verbose)`,
/// whose `verbose` is `options.verbose`.
pub fn enumerate_stereoisomers_with_options(
    m: &Mol,
    options: &StereoEnumerationOptions,
) -> Result<Stereoisomers> {
    let mut tm = m.clone();
    for idx in 0..tm.get_num_atoms() {
        tm.get_atom_with_idx_mut(idx)?.clear_prop("_CIPCode");
    }
    for idx in 0..tm.get_num_bonds() {
        let mut bond = tm.get_bond_with_idx_mut(idx)?;
        // RDKit 2024.03's Python clears only EITHERDOUBLE (later releases
        // clear UNKNOWN too).
        if bond.get_bond_dir() == BondDir::Eitherdouble {
            bond.set_bond_dir(BondDir::None);
        }
    }
    let flippers = get_flippers(&mut tm, options)?;
    let n = flippers.len() as u32;
    let mut res = Stereoisomers {
        tm,
        flippers,
        options: options.clone(),
        bitsource: None,
        unflipped: None,
        isomers_seen: HashSet::new(),
        num_isomers: 0,
    };
    if n == 0 {
        res.unflipped = Some(res.tm.clone());
        return Ok(res);
    }
    let max = u64::from(options.max_isomers);
    res.bitsource = Some(if max == 0 || (n < 64 && 1u64 << n <= max) {
        BitSource::Range {
            next: 0,
            end: (n < 128).then(|| 1u128 << n),
            n_words: n.div_ceil(32) as usize,
        }
    } else {
        let rand = match &options.rand {
            Some(rand) => rand.clone(),
            None => {
                // Python seeds with hash(tuple(sorted((degree, atomic
                // number) for each atom))).
                let mut keys: Vec<(u64, u64)> = res
                    .tm
                    .get_atoms()
                    .map(|a| (u64::from(a.get_degree()), a.get_atomic_num() as u64))
                    .collect();
                keys.sort();
                let lanes: Vec<u64> = keys
                    .iter()
                    .map(|&(d, z)| hash_tuple(&[d, z]) as u64)
                    .collect();
                Random::new(hash_tuple(&lanes))
            }
        };
        BitSource::UniqueRandom {
            n,
            rand: Box::new(rand),
            already_seen: HashSet::new(),
        }
    });
    Ok(res)
}
