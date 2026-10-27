// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdFMCS`: maximum common substructures of a set of
//! molecules (FMCS), with closures in place of the Python classes users
//! subclass to compare atoms and bonds, follow progress, and accept
//! matches.
//!
//! The wrapper checks what RDKit trusts (CONVENTIONS.md U1): an initial
//! seed that is not SMARTS, which RDKit dereferences, and a NaN threshold,
//! which it converts to `int`; the predefined checks the compare classes
//! offer take indices, which are checked. A closure never throws back into
//! RDKit's search: its first error (or panic) stops the search, which then
//! returns that error. Molecules without ring info get the empty one RDKit
//! gives them for the search and have none afterwards, even when RDKit
//! throws.
//!
//! The search matches substructures, so it holds the lock those take
//! (CONVENTIONS.md U4), and releases it while a closure runs, except a final
//! match check, which RDKit calls from inside a match: wrappers that match
//! substructures (or depict, or fragment on BRICS bonds) fail with
//! `Error::Runtime` when a final match check calls them, rather than
//! deadlock.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::Arc;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::rdmolops::{mol_list_into_map, with_borrowed_mols_shared};
use crate::error::{expect_ok, to_string};
use crate::ffi::fmcs::{McsCallbacks, ffi};

rdkit_enum! {
    /// RDKit's predefined atom comparisons. Mirrors `rdFMCS.AtomCompare`.
    AtomCompare: ffi::AtomComparator {
        CompareAny = AtomCompareAny,
        CompareElements = AtomCompareElements,
        CompareIsotopes = AtomCompareIsotopes,
        CompareAnyHeavyAtom = AtomCompareAnyHeavyAtom,
    }
}

rdkit_enum! {
    /// RDKit's predefined bond comparisons. Mirrors `rdFMCS.BondCompare`.
    BondCompare: ffi::BondComparator {
        CompareAny = BondCompareAny,
        CompareOrder = BondCompareOrder,
        CompareOrderExact = BondCompareOrderExact,
    }
}

/// How ring fusion constrains a match. Mirrors `rdFMCS.RingCompare`, which
/// only `FindMCS`'s keyword form takes; here
/// [`MCSBondCompareParameters::set_ring_compare`] sets the flags it stands
/// for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum RingCompare {
    /// Mirrors `RingCompare.IgnoreRingFusion`.
    IgnoreRingFusion,
    /// Mirrors `RingCompare.PermissiveRingFusion`.
    PermissiveRingFusion,
    /// Mirrors `RingCompare.StrictRingFusion`.
    StrictRingFusion,
}

/// How atoms are matched. Mirrors `rdFMCS.MCSAtomCompareParameters`;
/// `Default` reads RDKit's defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MCSAtomCompareParameters {
    /// Matches atoms only to atoms of the same total valence.
    pub match_valences: bool,
    /// Matches tetrahedral centers only to tetrahedral centers.
    pub match_chiral_tag: bool,
    /// Matches atoms only to atoms of the same formal charge.
    pub match_formal_charge: bool,
    /// Matches ring atoms only to ring atoms.
    pub ring_matches_ring_only: bool,
    /// Leaves out lone ring atoms; implies the bond parameters'
    /// `complete_rings_only`.
    pub complete_rings_only: bool,
    /// Writes isotope queries in the result's SMARTS.
    pub match_isotope: bool,
    /// Matches atoms only within this many angstroms of each other in the
    /// molecules' first conformers, when positive; a molecule without a
    /// conformer is then `Error::Value`.
    pub max_distance: f64,
}

impl Default for MCSAtomCompareParameters {
    fn default() -> MCSAtomCompareParameters {
        MCSAtomCompareParameters::from_ffi(&expect_ok(ffi::mcs_atom_compare_parameters_default()))
    }
}

impl MCSAtomCompareParameters {
    pub(crate) fn from_ffi(data: &ffi::MCSAtomCompareParametersData) -> MCSAtomCompareParameters {
        MCSAtomCompareParameters {
            match_valences: data.match_valences,
            match_chiral_tag: data.match_chiral_tag,
            match_formal_charge: data.match_formal_charge,
            ring_matches_ring_only: data.ring_matches_ring_only,
            complete_rings_only: data.complete_rings_only,
            match_isotope: data.match_isotope,
            max_distance: data.max_distance,
        }
    }

    fn to_ffi(self) -> ffi::MCSAtomCompareParametersData {
        ffi::MCSAtomCompareParametersData {
            match_valences: self.match_valences,
            match_chiral_tag: self.match_chiral_tag,
            match_formal_charge: self.match_formal_charge,
            ring_matches_ring_only: self.ring_matches_ring_only,
            complete_rings_only: self.complete_rings_only,
            match_isotope: self.match_isotope,
            max_distance: self.max_distance,
        }
    }
}

/// How bonds are matched. Mirrors `rdFMCS.MCSBondCompareParameters`;
/// `Default` reads RDKit's defaults.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MCSBondCompareParameters {
    /// Matches ring bonds only to ring bonds.
    pub ring_matches_ring_only: bool,
    /// Leaves out partial rings; implies `ring_matches_ring_only`.
    pub complete_rings_only: bool,
    /// Keeps ring fusion: alpha-methylnaphthalene does not match
    /// beta-methylnaphthalene, though decalin matches cyclodecane unless
    /// `match_fused_rings_strict` is set.
    pub match_fused_rings: bool,
    /// With `match_fused_rings`, requires the same fusion in both molecules,
    /// so decalin does not match cyclodecane.
    pub match_fused_rings_strict: bool,
    /// Matches double bonds only to double bonds of the same stereo.
    pub match_stereo: bool,
}

impl Default for MCSBondCompareParameters {
    fn default() -> MCSBondCompareParameters {
        MCSBondCompareParameters::from_ffi(&expect_ok(ffi::mcs_bond_compare_parameters_default()))
    }
}

impl MCSBondCompareParameters {
    pub(crate) fn from_ffi(data: &ffi::MCSBondCompareParametersData) -> MCSBondCompareParameters {
        MCSBondCompareParameters {
            ring_matches_ring_only: data.ring_matches_ring_only,
            complete_rings_only: data.complete_rings_only,
            match_fused_rings: data.match_fused_rings,
            match_fused_rings_strict: data.match_fused_rings_strict,
            match_stereo: data.match_stereo,
        }
    }

    fn to_ffi(self) -> ffi::MCSBondCompareParametersData {
        ffi::MCSBondCompareParametersData {
            ring_matches_ring_only: self.ring_matches_ring_only,
            complete_rings_only: self.complete_rings_only,
            match_fused_rings: self.match_fused_rings,
            match_fused_rings_strict: self.match_fused_rings_strict,
            match_stereo: self.match_stereo,
        }
    }

    /// Sets `match_fused_rings` and `match_fused_rings_strict` as
    /// `FindMCS`'s `ringCompare` keyword does: any fusion check but
    /// `IgnoreRingFusion` sets the first, and `StrictRingFusion` the second.
    pub fn set_ring_compare(&mut self, ring_compare: RingCompare) {
        self.match_fused_rings = ring_compare != RingCompare::IgnoreRingFusion;
        self.match_fused_rings_strict = ring_compare == RingCompare::StrictRingFusion;
    }
}

type AtomCompareFn =
    dyn Fn(&MCSAtomCompareParameters, &Mol, u32, &Mol, u32) -> Result<bool> + Send + Sync;
type BondCompareFn =
    dyn Fn(&MCSBondCompareParameters, &Mol, u32, &Mol, u32) -> Result<bool> + Send + Sync;
type ProgressFn = dyn Fn(&MCSProgressData, &MCSParameters) -> Result<bool> + Send + Sync;
type MatchCheckFn =
    dyn Fn(&Mol, &Mol, &[(u32, u32)], &[(u32, u32)], &MCSParameters) -> Result<bool> + Send + Sync;

/// An atom comparison written in Rust, for `MCSParameters::atom_typer`.
/// Mirrors `rdFMCS.MCSAtomCompare`: Python subclasses it and overrides
/// `__call__`, which here is the closure `new` takes. The closure gets the
/// atom compare parameters, then each molecule with an atom index, and
/// returns whether the atoms match; it is shared by clones of the value.
/// The predefined checks Python's subclasses call on `self` are associated
/// functions here, since a closure has no `self`.
#[derive(Clone)]
pub struct MCSAtomCompare(Arc<AtomCompareFn>);

impl fmt::Debug for MCSAtomCompare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MCSAtomCompare").finish_non_exhaustive()
    }
}

impl MCSAtomCompare {
    /// An atom comparison that runs `call`. An error it returns, or a panic
    /// in it, stops the search, and `find_mcs_with_params` returns the
    /// error (a panic as `Error::Runtime`). Mirrors `MCSAtomCompare()`.
    pub fn new<F>(call: F) -> MCSAtomCompare
    where
        F: Fn(&MCSAtomCompareParameters, &Mol, u32, &Mol, u32) -> Result<bool>
            + Send
            + Sync
            + 'static,
    {
        MCSAtomCompare(Arc::new(call))
    }

    /// Runs the comparison. Mirrors `MCSAtomCompare.__call__`.
    pub fn call(
        &self,
        parameters: &MCSAtomCompareParameters,
        mol1: &Mol,
        atom1: u32,
        mol2: &Mol,
        atom2: u32,
    ) -> Result<bool> {
        (self.0)(parameters, mol1, atom1, mol2, atom2)
    }

    /// Whether both atoms are in a ring or both are not, when
    /// `ring_matches_ring_only` is set (and otherwise true). A molecule
    /// without ring info is then `Error::Runtime`, and an atom index out of
    /// range `Error::Index`. Mirrors `MCSAtomCompare.CheckAtomRingMatch`.
    pub fn check_atom_ring_match(
        parameters: &MCSAtomCompareParameters,
        mol1: &Mol,
        atom1: u32,
        mol2: &Mol,
        atom2: u32,
    ) -> Result<bool> {
        Ok(ffi::mcs_atom_compare_check_atom_ring_match(
            &parameters.to_ffi(),
            mol1.as_ffi(),
            atom1,
            mol2.as_ffi(),
            atom2,
        )?)
    }

    /// Whether both atoms have the same formal charge; an atom index out of
    /// range is `Error::Index`. Mirrors `MCSAtomCompare.CheckAtomCharge`.
    pub fn check_atom_charge(
        parameters: &MCSAtomCompareParameters,
        mol1: &Mol,
        atom1: u32,
        mol2: &Mol,
        atom2: u32,
    ) -> Result<bool> {
        Ok(ffi::mcs_atom_compare_check_atom_charge(
            &parameters.to_ffi(),
            mol1.as_ffi(),
            atom1,
            mol2.as_ffi(),
            atom2,
        )?)
    }

    /// Whether the second atom is a tetrahedral center when the first is;
    /// an atom index out of range is `Error::Index`. Mirrors
    /// `MCSAtomCompare.CheckAtomChirality`.
    pub fn check_atom_chirality(
        parameters: &MCSAtomCompareParameters,
        mol1: &Mol,
        atom1: u32,
        mol2: &Mol,
        atom2: u32,
    ) -> Result<bool> {
        Ok(ffi::mcs_atom_compare_check_atom_chirality(
            &parameters.to_ffi(),
            mol1.as_ffi(),
            atom1,
            mol2.as_ffi(),
            atom2,
        )?)
    }
}

/// A bond comparison written in Rust, for `MCSParameters::bond_typer`.
/// Mirrors `rdFMCS.MCSBondCompare`, as [`MCSAtomCompare`] mirrors its atom
/// counterpart; the closure gets bond indices, and the predefined checks are
/// associated functions.
#[derive(Clone)]
pub struct MCSBondCompare(Arc<BondCompareFn>);

impl fmt::Debug for MCSBondCompare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MCSBondCompare").finish_non_exhaustive()
    }
}

impl MCSBondCompare {
    /// A bond comparison that runs `call`, whose errors and panics stop the
    /// search as [`MCSAtomCompare::new`]'s do. Mirrors `MCSBondCompare()`.
    pub fn new<F>(call: F) -> MCSBondCompare
    where
        F: Fn(&MCSBondCompareParameters, &Mol, u32, &Mol, u32) -> Result<bool>
            + Send
            + Sync
            + 'static,
    {
        MCSBondCompare(Arc::new(call))
    }

    /// Runs the comparison. Mirrors `MCSBondCompare.__call__`.
    pub fn call(
        &self,
        parameters: &MCSBondCompareParameters,
        mol1: &Mol,
        bond1: u32,
        mol2: &Mol,
        bond2: u32,
    ) -> Result<bool> {
        (self.0)(parameters, mol1, bond1, mol2, bond2)
    }

    /// Whether two double bonds have the same stereo (true for any other
    /// pair of bonds); a bond index out of range is `Error::Index`. Mirrors
    /// `MCSBondCompare.CheckBondStereo`.
    pub fn check_bond_stereo(
        parameters: &MCSBondCompareParameters,
        mol1: &Mol,
        bond1: u32,
        mol2: &Mol,
        bond2: u32,
    ) -> Result<bool> {
        Ok(ffi::mcs_bond_compare_check_bond_stereo(
            &parameters.to_ffi(),
            mol1.as_ffi(),
            bond1,
            mol2.as_ffi(),
            bond2,
        )?)
    }

    /// Whether both bonds are in a ring or both are not, and with
    /// `complete_rings_only` whether ring bonds share a ring that could
    /// match. A molecule without ring info is `Error::Runtime`, and a bond
    /// index out of range `Error::Index`. Mirrors
    /// `MCSBondCompare.CheckBondRingMatch`.
    pub fn check_bond_ring_match(
        parameters: &MCSBondCompareParameters,
        mol1: &Mol,
        bond1: u32,
        mol2: &Mol,
        bond2: u32,
    ) -> Result<bool> {
        Ok(ffi::mcs_bond_compare_check_bond_ring_match(
            &parameters.to_ffi(),
            mol1.as_ffi(),
            bond1,
            mol2.as_ffi(),
            bond2,
        )?)
    }
}

/// A progress callback written in Rust, for
/// `MCSParameters::progress_callback`. Mirrors `rdFMCS.MCSProgress`: the
/// closure gets the search's progress and parameters and returns whether to
/// go on; false cancels the search, whose result then has `canceled` set.
#[derive(Clone)]
pub struct MCSProgress(Arc<ProgressFn>);

impl fmt::Debug for MCSProgress {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MCSProgress").finish_non_exhaustive()
    }
}

impl MCSProgress {
    /// A progress callback that runs `call`, whose errors and panics stop
    /// the search as [`MCSAtomCompare::new`]'s do. Mirrors `MCSProgress()`.
    pub fn new<F>(call: F) -> MCSProgress
    where
        F: Fn(&MCSProgressData, &MCSParameters) -> Result<bool> + Send + Sync + 'static,
    {
        MCSProgress(Arc::new(call))
    }

    /// Runs the callback. Mirrors `MCSProgress.__call__`.
    pub fn call(&self, stat: &MCSProgressData, parameters: &MCSParameters) -> Result<bool> {
        (self.0)(stat, parameters)
    }
}

/// A check on each growing match, written in Rust, for
/// `MCSParameters::final_match_checker`. Mirrors `rdFMCS.MCSFinalMatchCheck`:
/// the closure gets the two molecules, the matched (first, second) atom and
/// bond index pairs, and the parameters, and returns whether to keep the
/// match. RDKit runs its own ring-fusion and chirality checks first. RDKit
/// calls it from inside a substructure match, whose lock it cannot release,
/// so wrappers that match substructures, depict, or fragment on BRICS bonds
/// return `Error::Runtime` when the closure calls them (the other callbacks
/// run without the lock).
#[derive(Clone)]
pub struct MCSFinalMatchCheck(Arc<MatchCheckFn>);

impl fmt::Debug for MCSFinalMatchCheck {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MCSFinalMatchCheck").finish_non_exhaustive()
    }
}

impl MCSFinalMatchCheck {
    /// A match check that runs `call`, whose errors and panics stop the
    /// search as [`MCSAtomCompare::new`]'s do. Mirrors
    /// `MCSFinalMatchCheck()`.
    pub fn new<F>(call: F) -> MCSFinalMatchCheck
    where
        F: Fn(&Mol, &Mol, &[(u32, u32)], &[(u32, u32)], &MCSParameters) -> Result<bool>
            + Send
            + Sync
            + 'static,
    {
        MCSFinalMatchCheck(Arc::new(call))
    }

    /// Runs the check. Mirrors `MCSFinalMatchCheck.__call__`.
    pub fn call(
        &self,
        mol1: &Mol,
        mol2: &Mol,
        atom_idx_match: &[(u32, u32)],
        bond_idx_match: &[(u32, u32)],
        parameters: &MCSParameters,
    ) -> Result<bool> {
        (self.0)(mol1, mol2, atom_idx_match, bond_idx_match, parameters)
    }
}

/// A check on each fully grown MCS candidate, written in Rust, for
/// `MCSParameters::should_accept_mcs`. Mirrors `rdFMCS.MCSAcceptance`: the
/// closure gets the query and a target molecule, the matched (query,
/// target) atom and bond index pairs, and the parameters, and returns
/// whether to accept the candidate.
#[derive(Clone)]
pub struct MCSAcceptance(Arc<MatchCheckFn>);

impl fmt::Debug for MCSAcceptance {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MCSAcceptance").finish_non_exhaustive()
    }
}

impl MCSAcceptance {
    /// An acceptance check that runs `call`, whose errors and panics stop
    /// the search as [`MCSAtomCompare::new`]'s do. Mirrors `MCSAcceptance()`.
    pub fn new<F>(call: F) -> MCSAcceptance
    where
        F: Fn(&Mol, &Mol, &[(u32, u32)], &[(u32, u32)], &MCSParameters) -> Result<bool>
            + Send
            + Sync
            + 'static,
    {
        MCSAcceptance(Arc::new(call))
    }

    /// Runs the check. Mirrors `MCSAcceptance.__call__`.
    pub fn call(
        &self,
        query: &Mol,
        target: &Mol,
        atom_idx_match: &[(u32, u32)],
        bond_idx_match: &[(u32, u32)],
        parameters: &MCSParameters,
    ) -> Result<bool> {
        (self.0)(query, target, atom_idx_match, bond_idx_match, parameters)
    }
}

/// The atom comparison `MCSParameters::atom_typer` holds: one of RDKit's or
/// a closure, as Python's `MCSParameters.AtomTyper` takes an `AtomCompare`
/// or an `MCSAtomCompare`.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum AtomTyper {
    /// One of RDKit's comparisons.
    Compare(AtomCompare),
    /// A comparison written in Rust.
    Custom(MCSAtomCompare),
}

impl From<AtomCompare> for AtomTyper {
    fn from(compare: AtomCompare) -> AtomTyper {
        AtomTyper::Compare(compare)
    }
}

impl From<MCSAtomCompare> for AtomTyper {
    fn from(compare: MCSAtomCompare) -> AtomTyper {
        AtomTyper::Custom(compare)
    }
}

/// The bond comparison `MCSParameters::bond_typer` holds: one of RDKit's or
/// a closure, as Python's `MCSParameters.BondTyper` takes a `BondCompare` or
/// an `MCSBondCompare`.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum BondTyper {
    /// One of RDKit's comparisons.
    Compare(BondCompare),
    /// A comparison written in Rust.
    Custom(MCSBondCompare),
}

impl From<BondCompare> for BondTyper {
    fn from(compare: BondCompare) -> BondTyper {
        BondTyper::Compare(compare)
    }
}

impl From<MCSBondCompare> for BondTyper {
    fn from(compare: MCSBondCompare) -> BondTyper {
        BondTyper::Custom(compare)
    }
}

/// How far a search has come, as a progress callback sees it. Mirrors
/// `rdFMCS.MCSProgressData`; `Default` is Python's `MCSProgressData()`, all
/// zeros.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MCSProgressData {
    num_atoms: u32,
    num_bonds: u32,
    seed_processed: u32,
}

impl MCSProgressData {
    pub(crate) fn new(num_atoms: u32, num_bonds: u32, seed_processed: u32) -> MCSProgressData {
        MCSProgressData {
            num_atoms,
            num_bonds,
            seed_processed,
        }
    }

    /// The atoms in the largest MCS so far. Mirrors
    /// `MCSProgressData.numAtoms`.
    pub fn num_atoms(&self) -> u32 {
        self.num_atoms
    }

    /// The bonds in the largest MCS so far. Mirrors
    /// `MCSProgressData.numBonds`.
    pub fn num_bonds(&self) -> u32 {
        self.num_bonds
    }

    /// The seeds processed (RDKit does not count them, so zero). Mirrors
    /// `MCSProgressData.seedProcessed`.
    pub fn seed_processed(&self) -> u32 {
        self.seed_processed
    }
}

/// How `find_mcs_with_params` searches. Mirrors `rdFMCS.MCSParameters`;
/// `Default` reads RDKit's defaults. The parameters a callback gets are
/// those of the search, which RDKit adjusts: `complete_rings_only` sets the
/// `ring_matches_ring_only` fields, the atom's sets the bond's, and
/// `match_chiral_tag` sets `match_stereo`.
#[derive(Clone, Debug)]
pub struct MCSParameters {
    /// Maximizes the bonds of the MCS, not its atoms.
    pub maximize_bonds: bool,
    /// The fraction of the molecules that must contain the MCS, at most 1
    /// (above is `Error::Runtime`; NaN is `Error::Value`).
    pub threshold: f64,
    /// Seconds to search before returning the best MCS so far, with
    /// `canceled` set; 0 for no limit. Not enforced when a progress
    /// callback is set, as in Python.
    pub timeout: u32,
    /// Prints statistics to standard output. With `store_all` it is
    /// `Error::Value`: RDKit 2024.03 then dereferences a null query
    /// molecule.
    pub verbose: bool,
    /// How atoms are matched.
    pub atom_compare_parameters: MCSAtomCompareParameters,
    /// How bonds are matched.
    pub bond_compare_parameters: MCSBondCompareParameters,
    /// The atom comparison.
    pub atom_typer: AtomTyper,
    /// The bond comparison.
    pub bond_typer: BondTyper,
    /// Called as the search goes on; returning false cancels it.
    pub progress_callback: Option<MCSProgress>,
    /// Checks each growing match.
    pub final_match_checker: Option<MCSFinalMatchCheck>,
    /// Checks each fully grown MCS candidate.
    pub should_accept_mcs: Option<MCSAcceptance>,
    /// SMARTS for the substructure the MCS grows from, or empty; one that
    /// does not parse is `Error::Parse`, and one that is not common to the
    /// molecules is ignored, with a warning.
    pub initial_seed: String,
    /// Keeps every MCS of the largest size, in the result's
    /// `degenerate_smarts_query_mol_dict`.
    pub store_all: bool,
}

impl Default for MCSParameters {
    fn default() -> MCSParameters {
        let data = expect_ok(ffi::mcs_parameters_default());
        MCSParameters {
            maximize_bonds: data.maximize_bonds,
            threshold: data.threshold,
            timeout: data.timeout,
            verbose: data.verbose,
            atom_compare_parameters: MCSAtomCompareParameters::from_ffi(
                &data.atom_compare_parameters,
            ),
            bond_compare_parameters: MCSBondCompareParameters::from_ffi(
                &data.bond_compare_parameters,
            ),
            atom_typer: AtomTyper::Compare(AtomCompare::from_ffi(data.atom_typer)),
            bond_typer: BondTyper::Compare(BondCompare::from_ffi(data.bond_typer)),
            progress_callback: None,
            final_match_checker: None,
            should_accept_mcs: None,
            initial_seed: data.initial_seed,
            store_all: data.store_all,
        }
    }
}

impl MCSParameters {
    /// The parameters of the search a callback is called from: RDKit's
    /// values, with the closures (which RDKit holds as trampolines) from
    /// the parameters the search was given.
    pub(crate) fn from_callback(
        data: &ffi::MCSParametersData,
        given: &MCSParameters,
    ) -> MCSParameters {
        MCSParameters {
            maximize_bonds: data.maximize_bonds,
            threshold: data.threshold,
            timeout: data.timeout,
            verbose: data.verbose,
            atom_compare_parameters: MCSAtomCompareParameters::from_ffi(
                &data.atom_compare_parameters,
            ),
            bond_compare_parameters: MCSBondCompareParameters::from_ffi(
                &data.bond_compare_parameters,
            ),
            atom_typer: given.atom_typer.clone(),
            bond_typer: given.bond_typer.clone(),
            progress_callback: given.progress_callback.clone(),
            final_match_checker: given.final_match_checker.clone(),
            should_accept_mcs: given.should_accept_mcs.clone(),
            initial_seed: given.initial_seed.clone(),
            store_all: data.store_all,
        }
    }

    fn to_ffi(&self) -> ffi::MCSParametersData {
        let (atom_typer, custom_atom_typer) = match &self.atom_typer {
            AtomTyper::Compare(compare) => (compare.to_ffi(), false),
            AtomTyper::Custom(_) => (AtomCompare::CompareElements.to_ffi(), true),
        };
        let (bond_typer, custom_bond_typer) = match &self.bond_typer {
            BondTyper::Compare(compare) => (compare.to_ffi(), false),
            BondTyper::Custom(_) => (BondCompare::CompareOrder.to_ffi(), true),
        };
        ffi::MCSParametersData {
            maximize_bonds: self.maximize_bonds,
            threshold: self.threshold,
            timeout: self.timeout,
            verbose: self.verbose,
            atom_compare_parameters: self.atom_compare_parameters.to_ffi(),
            bond_compare_parameters: self.bond_compare_parameters.to_ffi(),
            atom_typer,
            custom_atom_typer,
            bond_typer,
            custom_bond_typer,
            progress_callback: self.progress_callback.is_some(),
            final_match_checker: self.final_match_checker.is_some(),
            should_accept_mcs: self.should_accept_mcs.is_some(),
            initial_seed: self.initial_seed.clone(),
            store_all: self.store_all,
        }
    }
}

/// A maximum common substructure. Mirrors `rdFMCS.MCSResult`.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct MCSResult {
    num_atoms: u32,
    num_bonds: u32,
    smarts_string: String,
    canceled: bool,
    query_mol: Option<Mol>,
    degenerate_smarts_query_mol_dict: BTreeMap<String, Mol>,
}

impl MCSResult {
    fn from_ffi(data: ffi::MCSResultData) -> Result<MCSResult> {
        Ok(MCSResult {
            num_atoms: data.num_atoms,
            num_bonds: data.num_bonds,
            smarts_string: to_string(data.smarts_string)?,
            canceled: data.canceled,
            query_mol: (!data.query_mol.is_null()).then(|| Mol::from_ffi(data.query_mol)),
            degenerate_smarts_query_mol_dict: mol_list_into_map(
                data.degenerate_smarts_query_mol_dict,
            )?,
        })
    }

    /// The atoms in the MCS. Mirrors `MCSResult.numAtoms`.
    pub fn num_atoms(&self) -> u32 {
        self.num_atoms
    }

    /// The bonds in the MCS. Mirrors `MCSResult.numBonds`.
    pub fn num_bonds(&self) -> u32 {
        self.num_bonds
    }

    /// The MCS as SMARTS; empty when there is none or `store_all` was set.
    /// Mirrors `MCSResult.smartsString`.
    pub fn smarts_string(&self) -> &str {
        &self.smarts_string
    }

    /// Whether the search stopped early, at its timeout or when a progress
    /// callback returned false; the result is the best MCS found by then.
    /// Mirrors `MCSResult.canceled`.
    pub fn canceled(&self) -> bool {
        self.canceled
    }

    /// The MCS as a query molecule; `None` when there is none or
    /// `store_all` was set. Mirrors `MCSResult.queryMol`.
    pub fn query_mol(&self) -> Option<&Mol> {
        self.query_mol.as_ref()
    }

    /// With `store_all`, every MCS of the largest size, by SMARTS; empty
    /// otherwise. Mirrors `MCSResult.degenerateSmartsQueryMolDict`.
    pub fn degenerate_smarts_query_mol_dict(&self) -> &BTreeMap<String, Mol> {
        &self.degenerate_smarts_query_mol_dict
    }
}

/// The MCS of `mols`, with the defaults of `FindMCS`'s keyword form: those
/// of `MCSParameters::default()` with a `timeout` of 3600 seconds. Fewer
/// than two molecules are `Error::Runtime`. Mirrors `rdFMCS.FindMCS(mols)`.
pub fn find_mcs(mols: &[Mol]) -> Result<MCSResult> {
    find_mcs_with_params(
        mols,
        &MCSParameters {
            timeout: 3600,
            ..MCSParameters::default()
        },
    )
}

/// The MCS of `mols`, with `parameters`. Mirrors
/// `rdFMCS.FindMCS(mols, parameters)` and the keyword form, whose arguments
/// are fields here: `ringMatchesRingOnly` sets both compare parameters'
/// `ring_matches_ring_only`, `completeRingsOnly` and `matchValences` and
/// `matchChiralTag` set the fields of those names,
/// `atomCompare` and `bondCompare` are the typers, `ringCompare` is
/// [`MCSBondCompareParameters::set_ring_compare`], and `seedSmarts` is
/// `initial_seed`.
///
/// Fewer than two molecules, or a `threshold` above 1, are
/// `Error::Runtime`; a NaN `threshold` is `Error::Value`, and an
/// `initial_seed` that is not SMARTS `Error::Parse`. An error a closure
/// returns, or a panic in one (as `Error::Runtime`), stops the search and is
/// returned. RDKit's errors pass through: `Error::Value` for a positive
/// `max_distance` and a molecule without a conformer, `Error::Runtime` for
/// `match_valences` and a molecule whose valences are not computed.
///
/// `verbose` with `store_all` is `Error::Value`: RDKit 2024.03 then
/// matches the molecules against a query molecule it never set.
///
/// ```no_run
/// use rdkit_cxx::chem::mol_from_smiles;
/// use rdkit_cxx::chem::rd_fmcs::{AtomCompare, MCSParameters, find_mcs_with_params};
///
/// let mols = [mol_from_smiles("c1ccccc1CCO")?, mol_from_smiles("c1ccccc1CCN")?];
/// let params = MCSParameters {
///     atom_typer: AtomCompare::CompareAny.into(),
///     ..MCSParameters::default()
/// };
/// let mcs = find_mcs_with_params(&mols, &params)?;
/// assert_eq!(mcs.num_atoms(), 9);
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
pub fn find_mcs_with_params(mols: &[Mol], parameters: &MCSParameters) -> Result<MCSResult> {
    let data = parameters.to_ffi();
    let result = McsCallbacks::with_callbacks(mols, parameters, |callbacks| {
        with_borrowed_mols_shared(mols.iter().map(|mol| ("", mol)), |list| {
            ffi::find_mcs(list, &data, callbacks)
        })
    })?;
    MCSResult::from_ffi(result?)
}
