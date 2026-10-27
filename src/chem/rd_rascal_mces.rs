// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdRascalMCES`: maximum common edge substructures
//! with the RASCAL algorithm, and the clusterings built on its similarity.
//!
//! The wrapper checks what RDKit trusts (CONVENTIONS.md U1): a
//! `max_frag_separation` other than -1 and a NaN similarity threshold, which
//! RDKit reads past a matrix (or, in 2025.09, through a pruned fragment)
//! and casts, a timeout other than -1, on which RDKit 2025.09 reads an
//! emptied list of cliques (the clusterings run their searches without
//! one), and equivalent-atom SMARTS that do not parse, which RDKit
//! dereferences; `complete_smallest_rings` sets
//! `complete_aromatic_rings`, as RDKit documents but does not do (it reads
//! ring data only the latter makes). Molecules past the sizes RDKit's `int`
//! arithmetic allows are refused. Each MCES matches substructures, so the
//! clusterings run their pairs on the calling thread (CONVENTIONS.md U4),
//! where an error is returned rather than ending the process.

use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::error::{expect_ok, to_string};
use crate::ffi::chirality::ffi::IndexPair;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::rascal_mces::ffi;

/// Options for `find_mces_with_opts`. Mirrors `rdRascalMCES.RascalOptions`;
/// `Default` reads RDKit's defaults.
#[derive(Clone, Debug, PartialEq)]
pub struct RascalOptions {
    /// The similarity estimate below which no MCES is found, from 0 to 1;
    /// NaN is `Error::Value`.
    pub similarity_threshold: f64,
    /// Leaves out partial aromatic rings.
    pub complete_aromatic_rings: bool,
    /// Matches ring bonds only to ring bonds.
    pub ring_matches_ring_only: bool,
    /// Keeps only rings complete in both molecules' ring info; implies
    /// `complete_aromatic_rings` and `ring_matches_ring_only`.
    pub complete_smallest_rings: bool,
    /// Matches atoms only to atoms with as many explicit connections.
    pub exact_connections_match: bool,
    /// Returns only the largest fragment of the MCES.
    pub single_largest_frag: bool,
    /// The fewest atoms in a fragment of the MCES; -1 for no minimum.
    pub min_frag_size: i32,
    /// The most bonds between fragments of the MCES; -1 for no limit, and
    /// any other value is `Error::Value`: RDKit 2025.09 measures from a
    /// fragment it has just pruned.
    pub max_frag_separation: i32,
    /// Returns every MCES of the largest size, not just the first.
    pub all_best_mcess: bool,
    /// The most MCESs of one size to keep.
    pub max_best_mcess: u32,
    /// Seconds to spend on one MCES; -1 for no limit. RDKit 2025.09 reads
    /// an emptied list when a search times out, so other values are
    /// `Error::Value`, and the default is -1, not Python's 60.
    pub timeout: i32,
    /// Returns a result with the similarity estimates when they miss the
    /// threshold.
    pub return_empty_mces: bool,
    /// The most matching bond pairs to search; memory grows with its square.
    pub max_bond_match_pairs: u32,
    /// Space-separated SMARTS, each a class of atoms to treat as equivalent,
    /// at most 9; one that does not parse is `Error::Parse`.
    pub equivalent_atoms: String,
    /// Matches bonds whatever their order.
    pub ignore_bond_orders: bool,
    /// Matches atoms on their atomic numbers alone.
    pub ignore_atom_aromaticity: bool,
    /// The fewest bonds in the MCES; overrides `similarity_threshold` when
    /// above 0.
    pub min_clique_size: u32,
}

impl Default for RascalOptions {
    fn default() -> RascalOptions {
        let d = expect_ok(ffi::rascal_options_default());
        RascalOptions {
            similarity_threshold: d.similarity_threshold,
            complete_aromatic_rings: d.complete_aromatic_rings,
            ring_matches_ring_only: d.ring_matches_ring_only,
            complete_smallest_rings: d.complete_smallest_rings,
            exact_connections_match: d.exact_connections_match,
            single_largest_frag: d.single_largest_frag,
            min_frag_size: d.min_frag_size,
            max_frag_separation: d.max_frag_separation,
            all_best_mcess: d.all_best_mcess,
            max_best_mcess: d.max_best_mcess,
            timeout: d.timeout,
            return_empty_mces: d.return_empty_mces,
            max_bond_match_pairs: d.max_bond_match_pairs,
            equivalent_atoms: d.equivalent_atoms,
            ignore_bond_orders: d.ignore_bond_orders,
            ignore_atom_aromaticity: d.ignore_atom_aromaticity,
            min_clique_size: d.min_clique_size,
        }
    }
}

impl RascalOptions {
    pub(crate) fn to_ffi(&self) -> ffi::RascalOptionsData {
        ffi::RascalOptionsData {
            similarity_threshold: self.similarity_threshold,
            complete_aromatic_rings: self.complete_aromatic_rings,
            ring_matches_ring_only: self.ring_matches_ring_only,
            complete_smallest_rings: self.complete_smallest_rings,
            exact_connections_match: self.exact_connections_match,
            single_largest_frag: self.single_largest_frag,
            min_frag_size: self.min_frag_size,
            max_frag_separation: self.max_frag_separation,
            all_best_mcess: self.all_best_mcess,
            max_best_mcess: self.max_best_mcess,
            timeout: self.timeout,
            return_empty_mces: self.return_empty_mces,
            max_bond_match_pairs: self.max_bond_match_pairs,
            equivalent_atoms: self.equivalent_atoms.clone(),
            ignore_bond_orders: self.ignore_bond_orders,
            ignore_atom_aromaticity: self.ignore_atom_aromaticity,
            min_clique_size: self.min_clique_size,
        }
    }
}

/// Options for the clusterings. Mirrors `rdRascalMCES.RascalClusterOptions`;
/// `Default` reads RDKit's defaults. `rascal_butina_cluster` uses only
/// `similarity_cutoff` and `num_threads`.
#[derive(Clone, Debug, PartialEq)]
pub struct RascalClusterOptions {
    /// The similarity two molecules need to share a cluster, from 0 to 1;
    /// NaN is `Error::Value`.
    pub similarity_cutoff: f64,
    /// The penalty for each unconnected fragment of an MCES.
    pub a: f64,
    /// The weight of matched bonds over matched atoms.
    pub b: f64,
    /// The fewest atoms in a fragment counted in an MCES.
    pub min_frag_size: u32,
    /// The similarity of two pairs' MCESs for them to share a cluster.
    pub min_intra_cluster_sim: f64,
    /// The fraction of shared molecules above which two clusters merge.
    pub cluster_merge_sim: f64,
    /// The most fragments in an MCES.
    pub max_num_frags: u32,
    /// Threads RDKit would use: this many above 0, else all but that many
    /// of the cores. Not read: every MCES matches substructures, which take
    /// a lock (CONVENTIONS.md U4), so the pairs run on the calling thread.
    pub num_threads: i32,
}

impl Default for RascalClusterOptions {
    fn default() -> RascalClusterOptions {
        let d = expect_ok(ffi::rascal_cluster_options_default());
        RascalClusterOptions {
            similarity_cutoff: d.similarity_cutoff,
            a: d.a,
            b: d.b,
            min_frag_size: d.min_frag_size,
            min_intra_cluster_sim: d.min_intra_cluster_sim,
            cluster_merge_sim: d.cluster_merge_sim,
            max_num_frags: d.max_num_frags,
            num_threads: d.num_threads,
        }
    }
}

impl RascalClusterOptions {
    fn to_ffi(&self) -> ffi::RascalClusterOptionsData {
        ffi::RascalClusterOptionsData {
            similarity_cutoff: self.similarity_cutoff,
            a: self.a,
            b: self.b,
            min_frag_size: self.min_frag_size,
            min_intra_cluster_sim: self.min_intra_cluster_sim,
            cluster_merge_sim: self.cluster_merge_sim,
            max_num_frags: self.max_num_frags,
            num_threads: self.num_threads,
        }
    }
}

/// An MCES between two molecules. Mirrors `rdRascalMCES.RascalResult`.
///
/// `RascalResult` is `Send` but not `Sync`: its getters cache what they
/// compute in the result.
#[non_exhaustive]
pub struct RascalResult {
    /// Never null.
    ptr: UniquePtr<ffi::RascalResult>,
}

fn pairs(pairs: Vec<IndexPair>) -> Vec<(u32, u32)> {
    pairs.into_iter().map(|p| (p.first, p.second)).collect()
}

impl RascalResult {
    fn raw(&self) -> &ffi::RascalResult {
        self.ptr.as_ref().unwrap()
    }

    /// SMARTS for the MCES; empty for a result without one. Mirrors
    /// `RascalResult.smartsString`.
    pub fn smarts_string(&self) -> String {
        to_string(expect_ok(ffi::rascal_result_smarts_string(self.raw())))
            .expect("RDKit's SMARTS is ASCII")
    }

    /// The matching bonds, as (index in the first molecule, index in the
    /// second). Mirrors `RascalResult.bondMatches`.
    pub fn bond_matches(&self) -> Vec<(u32, u32)> {
        pairs(expect_ok(ffi::rascal_result_bond_matches(self.raw())))
    }

    /// The matching atoms, as `bond_matches`. Mirrors
    /// `RascalResult.atomMatches`.
    pub fn atom_matches(&self) -> Vec<(u32, u32)> {
        pairs(expect_ok(ffi::rascal_result_atom_matches(self.raw())))
    }

    /// Cuts the MCES down to its largest fragment, for good. A result
    /// without matches is unchanged (Python's crashes for one
    /// `return_empty_mces` made). Mirrors `RascalResult.largestFragmentOnly`.
    pub fn largest_fragment_only(&mut self) {
        expect_ok(ffi::rascal_result_largest_fragment_only(self.ptr.pin_mut()))
    }

    /// The Johnson similarity of the two molecules. Mirrors
    /// `RascalResult.similarity`.
    pub fn similarity(&self) -> f64 {
        expect_ok(ffi::rascal_result_similarity(self.raw()))
    }

    /// The number of fragments in the MCES. Mirrors
    /// `RascalResult.numFragments`.
    pub fn num_fragments(&self) -> i32 {
        expect_ok(ffi::rascal_result_num_fragments(self.raw()))
    }

    /// The number of atoms in the MCES's largest fragment. Mirrors
    /// `RascalResult.largestFragmentSize`.
    pub fn largest_fragment_size(&self) -> u32 {
        expect_ok(ffi::rascal_result_largest_fragment_size(self.raw()))
    }

    /// The first, quicker similarity estimate. Mirrors
    /// `RascalResult.tier1Sim`.
    pub fn tier1_sim(&self) -> f64 {
        expect_ok(ffi::rascal_result_tier1_sim(self.raw()))
    }

    /// The second similarity estimate. Mirrors `RascalResult.tier2Sim`.
    pub fn tier2_sim(&self) -> f64 {
        expect_ok(ffi::rascal_result_tier2_sim(self.raw()))
    }

    /// Whether the search timed out, leaving the best MCES found so far.
    /// Mirrors `RascalResult.timedOut`.
    pub fn timed_out(&self) -> bool {
        expect_ok(ffi::rascal_result_timed_out(self.raw()))
    }
}

/// A copy, with its own molecules.
impl Clone for RascalResult {
    fn clone(&self) -> RascalResult {
        RascalResult {
            ptr: expect_ok(ffi::rascal_result_clone(self.raw())),
        }
    }
}

impl fmt::Debug for RascalResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RascalResult")
            .field("smarts_string", &self.smarts_string())
            .field("similarity", &self.similarity())
            .finish_non_exhaustive()
    }
}

/// The MCESs of two molecules: the best, or every one of the best size with
/// `all_best_mcess`; none when the similarity estimates miss the threshold.
/// The molecules need ring info. Mirrors `rdRascalMCES.FindMCES`.
pub fn find_mces(mol1: &Mol, mol2: &Mol) -> Result<Vec<RascalResult>> {
    find_mces_with_opts(mol1, mol2, &RascalOptions::default())
}

/// `Error::Value` for the options the type documents, or for molecules
/// whose sizes RDKit's `int` arithmetic would overflow: one with more than
/// 46340 atoms (its distance matrix), or both with more than 46340 atoms and
/// bonds (the similarity estimates square a count of both). Mirrors
/// `rdRascalMCES.FindMCES(mol1, mol2, opts)`.
pub fn find_mces_with_opts(
    mol1: &Mol,
    mol2: &Mol,
    opts: &RascalOptions,
) -> Result<Vec<RascalResult>> {
    let mut list = ffi::find_mces(mol1.as_ffi(), mol2.as_ffi(), &opts.to_ffi())?;
    let size = expect_ok(ffi::rascal_result_list_size(&list));
    Ok((0..size)
        .map(|i| RascalResult {
            ptr: expect_ok(ffi::rascal_result_list_take(list.pin_mut(), i)),
        })
        .collect())
}

/// Fuzzy clusters of molecules by their MCES similarity, as indices into
/// `mols`; a molecule may be in several. The last cluster holds the
/// molecules in no other. Mirrors `rdRascalMCES.RascalCluster`.
pub fn rascal_cluster(mols: &[Mol]) -> Result<Vec<Vec<u32>>> {
    rascal_cluster_with_opts(mols, &RascalClusterOptions::default())
}

/// The molecules' sizes are checked as `find_mces_with_opts` checks a
/// pair's. The pairs run on the calling thread, and an error from any is
/// returned (Python's ends the process when it comes from another thread).
/// Mirrors `rdRascalMCES.RascalCluster(mols, opts)`.
pub fn rascal_cluster_with_opts(
    mols: &[Mol],
    opts: &RascalClusterOptions,
) -> Result<Vec<Vec<u32>>> {
    let clusters = with_borrowed_mols_shared(mols.iter().map(|m| ("", m)), |list| {
        ffi::rascal_cluster(list, &opts.to_ffi())
    })?;
    Ok(IndexList::into_nested(clusters))
}

/// Butina clusters of molecules by their MCES similarity, as indices into
/// `mols`, each molecule in one. Mirrors `rdRascalMCES.RascalButinaCluster`.
pub fn rascal_butina_cluster(mols: &[Mol]) -> Result<Vec<Vec<u32>>> {
    rascal_butina_cluster_with_opts(mols, &RascalClusterOptions::default())
}

/// As `rascal_cluster_with_opts`. Mirrors
/// `rdRascalMCES.RascalButinaCluster(mols, opts)`.
pub fn rascal_butina_cluster_with_opts(
    mols: &[Mol],
    opts: &RascalClusterOptions,
) -> Result<Vec<Vec<u32>>> {
    let clusters = with_borrowed_mols_shared(mols.iter().map(|m| ("", m)), |list| {
        ffi::rascal_butina_cluster(list, &opts.to_ffi())
    })?;
    Ok(IndexList::into_nested(clusters))
}
