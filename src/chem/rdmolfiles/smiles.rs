// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! CXSMILES, fragments, random SMILES, atoms from SMILES, and canonical
//! ranking.

use crate::Result;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::{Atom, Bond, Mol, SmilesWriteParams};
use crate::error::to_string;
use crate::ffi::molfiles::ffi;

/// Which CXSMILES extensions to write, as bits. Mirrors
/// `rdmolfiles.CXSmilesFields` (CONVENTIONS.md F1); the constants keep
/// Python's spelling.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct CXSmilesFields(u32);

macro_rules! cx_smiles_fields {
    ($($name:ident),* $(,)?) => {
        impl CXSmilesFields {
            $(
                #[doc = concat!("Mirrors `CXSmilesFields.", stringify!($name), "`.")]
                pub const $name: CXSmilesFields = CXSmilesFields(ffi::CXSmilesFields::$name.repr);
            )*

            const NAMES: &[(&str, CXSmilesFields)] =
                &[$((stringify!($name), CXSmilesFields::$name)),*];
        }
    };
}

cx_smiles_fields!(
    CX_NONE,
    CX_ATOM_LABELS,
    CX_MOLFILE_VALUES,
    CX_COORDS,
    CX_RADICALS,
    CX_ATOM_PROPS,
    CX_LINKNODES,
    CX_ENHANCEDSTEREO,
    CX_SGROUPS,
    CX_POLYMER,
    CX_BOND_CFG,
    CX_BOND_ATROPISOMER,
    CX_COORDINATE_BONDS,
    CX_ZERO_BONDS,
    CX_ALL,
    CX_ALL_BUT_COORDS,
);

crate::flags::flag_ops!(CXSmilesFields);

rdkit_enum! {
    /// What CXSMILES writing does with bond directions. Mirrors
    /// `rdmolfiles.RestoreBondDirOption`.
    RestoreBondDirOption: ffi::RestoreBondDirOption {
        RestoreBondDirOptionTrue = RestoreBondDirOptionTrue,
        RestoreBondDirOptionClear = RestoreBondDirOptionClear,
    }
}

/// An atom from the SMILES for one atom. Mirrors `rdmolfiles.AtomFromSmiles`,
/// which returns `None` where this returns `Error::Parse`.
///
/// # Safety
///
/// RDKit writes the SMILES parser's debug level, a global that every SMILES
/// parse reads, without a lock (`yysmiles_debug` in `SmilesParse.cpp`). No
/// other thread may parse SMILES during the call, through this crate or inside
/// RDKit, including the parsing threads of a live
/// `MultithreadedSmilesMolSupplier` (CONVENTIONS.md U11).
pub unsafe fn atom_from_smiles(smiles: &str) -> Result<Atom> {
    Ok(Atom::from_ffi(ffi::atom_from_smiles(smiles)?))
}

/// A query atom from the SMARTS for one atom. Mirrors
/// `rdmolfiles.AtomFromSmarts`, which returns `None` where this returns
/// `Error::Parse`.
///
/// # Safety
///
/// RDKit writes the SMARTS parser's debug level, a global that every SMARTS
/// parse reads, without a lock (`yysmarts_debug` in `SmilesParse.cpp`). No
/// other thread may parse SMARTS during the call, through this crate or inside
/// RDKit (U11).
pub unsafe fn atom_from_smarts(smarts: &str) -> Result<Atom> {
    Ok(Atom::from_ffi(ffi::atom_from_smarts(smarts)?))
}

/// A bond from the SMILES for one bond. Mirrors `rdmolfiles.BondFromSmiles`,
/// which returns `None` where this returns `Error::Parse`.
///
/// # Safety
///
/// As `atom_from_smiles`: RDKit writes the SMILES parser's debug level.
pub unsafe fn bond_from_smiles(smiles: &str) -> Result<Bond> {
    Ok(Bond::from_ffi(ffi::bond_from_smiles(smiles)?))
}

/// A query bond from the SMARTS for one bond. Mirrors
/// `rdmolfiles.BondFromSmarts`, which returns `None` where this returns
/// `Error::Parse`.
///
/// # Safety
///
/// As `atom_from_smarts`: RDKit writes the SMARTS parser's debug level.
pub unsafe fn bond_from_smarts(smarts: &str) -> Result<Bond> {
    Ok(Bond::from_ffi(ffi::bond_from_smarts(smarts)?))
}

/// Parameters for `mol_to_cx_smiles_with_params`. Mirrors
/// `rdmolfiles.MolToCXSmiles(mol, params, flags, restoreBondDirs)`, whose
/// keyword overload takes `params`' fields.
#[derive(Clone, Debug, PartialEq)]
pub struct MolToCXSmilesParams {
    pub params: SmilesWriteParams,
    /// The extensions to write.
    pub flags: CXSmilesFields,
    pub restore_bond_dirs: RestoreBondDirOption,
}

impl Default for MolToCXSmilesParams {
    fn default() -> MolToCXSmilesParams {
        MolToCXSmilesParams {
            params: SmilesWriteParams::default(),
            flags: CXSmilesFields::CX_ALL,
            restore_bond_dirs: RestoreBondDirOption::RestoreBondDirOptionClear,
        }
    }
}

/// The canonical CXSMILES for a molecule. Mirrors `rdmolfiles.MolToCXSmiles`.
pub fn mol_to_cx_smiles(mol: &Mol) -> Result<String> {
    mol_to_cx_smiles_with_params(mol, &MolToCXSmilesParams::default())
}

/// Mirrors `rdmolfiles.MolToCXSmiles(mol, params, flags, restoreBondDirs)` and
/// its keyword form.
pub fn mol_to_cx_smiles_with_params(mol: &Mol, params: &MolToCXSmilesParams) -> Result<String> {
    to_string(ffi::mol_to_cx_smiles(
        mol.as_ffi(),
        &params.params.to_ffi(),
        params.flags.bits(),
        params.restore_bond_dirs.to_ffi(),
    )?)
}

/// A CXSMARTS for a molecule. Mirrors `rdmolfiles.MolToCXSmarts`.
pub fn mol_to_cx_smarts(mol: &Mol) -> Result<String> {
    mol_to_cx_smarts_with_isomeric_smiles(mol, true)
}

/// Mirrors `rdmolfiles.MolToCXSmarts(mol, isomericSmiles)`.
pub fn mol_to_cx_smarts_with_isomeric_smiles(mol: &Mol, isomeric_smiles: bool) -> Result<String> {
    to_string(ffi::mol_to_cx_smarts(mol.as_ffi(), isomeric_smiles)?)
}

/// Parameters for `mol_fragment_to_smiles_with_params` and
/// `mol_fragment_to_cx_smiles_with_params`, with Python's defaults. Mirrors
/// `rdmolfiles.MolFragmentToSmiles(mol, params, atomsToUse, bondsToUse,
/// atomSymbols, bondSymbols)`, whose keyword overload takes `params`' fields.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MolFragmentToSmilesParams {
    pub params: SmilesWriteParams,
    /// The bonds to include; `None` for every bond between the atoms.
    pub bonds_to_use: Option<Vec<u32>>,
    /// A symbol for each atom of the molecule, written in place of its SMILES.
    pub atom_symbols: Option<Vec<String>>,
    /// A symbol for each bond of the molecule.
    pub bond_symbols: Option<Vec<String>>,
}

fn fragment_to_smiles(
    mol: &Mol,
    atoms_to_use: &[u32],
    params: &MolFragmentToSmilesParams,
    cx: bool,
) -> Result<String> {
    let none = Vec::new();
    to_string(ffi::mol_fragment_to_smiles(
        mol.as_ffi(),
        &params.params.to_ffi(),
        atoms_to_use,
        params.bonds_to_use.is_some(),
        params.bonds_to_use.as_deref().unwrap_or_default(),
        params.atom_symbols.is_some(),
        params.atom_symbols.as_ref().unwrap_or(&none),
        params.bond_symbols.is_some(),
        params.bond_symbols.as_ref().unwrap_or(&none),
        cx,
    )?)
}

/// The SMILES for the fragment of a molecule made of `atoms_to_use`, which
/// must not be empty (`Error::Value`); an index out of range is
/// `Error::Index`. Mirrors `rdmolfiles.MolFragmentToSmiles`.
pub fn mol_fragment_to_smiles(mol: &Mol, atoms_to_use: &[u32]) -> Result<String> {
    mol_fragment_to_smiles_with_params(mol, atoms_to_use, &MolFragmentToSmilesParams::default())
}

/// Symbol lists must have one entry per atom or bond (`Error::Value`).
/// Mirrors `rdmolfiles.MolFragmentToSmiles` with its optional arguments.
pub fn mol_fragment_to_smiles_with_params(
    mol: &Mol,
    atoms_to_use: &[u32],
    params: &MolFragmentToSmilesParams,
) -> Result<String> {
    fragment_to_smiles(mol, atoms_to_use, params, false)
}

/// The CXSMILES for a fragment, as `mol_fragment_to_smiles`. Mirrors
/// `rdmolfiles.MolFragmentToCXSmiles`.
pub fn mol_fragment_to_cx_smiles(mol: &Mol, atoms_to_use: &[u32]) -> Result<String> {
    mol_fragment_to_cx_smiles_with_params(mol, atoms_to_use, &MolFragmentToSmilesParams::default())
}

/// Mirrors `rdmolfiles.MolFragmentToCXSmiles` with its optional arguments.
pub fn mol_fragment_to_cx_smiles_with_params(
    mol: &Mol,
    atoms_to_use: &[u32],
    params: &MolFragmentToSmilesParams,
) -> Result<String> {
    fragment_to_smiles(mol, atoms_to_use, params, true)
}

/// Parameters for `mol_fragment_to_smarts_with_params` and
/// `mol_fragment_to_cx_smarts_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolFragmentToSmartsParams {
    /// The bonds to include; `None` for every bond between the atoms.
    pub bonds_to_use: Option<Vec<u32>>,
    pub isomeric_smarts: bool,
}

impl Default for MolFragmentToSmartsParams {
    fn default() -> MolFragmentToSmartsParams {
        MolFragmentToSmartsParams {
            bonds_to_use: None,
            isomeric_smarts: true,
        }
    }
}

fn fragment_to_smarts(
    mol: &Mol,
    atoms_to_use: &[u32],
    params: &MolFragmentToSmartsParams,
    cx: bool,
) -> Result<String> {
    to_string(ffi::mol_fragment_to_smarts(
        mol.as_ffi(),
        atoms_to_use,
        params.bonds_to_use.is_some(),
        params.bonds_to_use.as_deref().unwrap_or_default(),
        params.isomeric_smarts,
        cx,
    )?)
}

/// The SMARTS for a fragment, as `mol_fragment_to_smiles`. Mirrors
/// `rdmolfiles.MolFragmentToSmarts`.
pub fn mol_fragment_to_smarts(mol: &Mol, atoms_to_use: &[u32]) -> Result<String> {
    mol_fragment_to_smarts_with_params(mol, atoms_to_use, &MolFragmentToSmartsParams::default())
}

/// Mirrors `rdmolfiles.MolFragmentToSmarts(mol, atomsToUse, bondsToUse,
/// isomericSmarts)`.
pub fn mol_fragment_to_smarts_with_params(
    mol: &Mol,
    atoms_to_use: &[u32],
    params: &MolFragmentToSmartsParams,
) -> Result<String> {
    fragment_to_smarts(mol, atoms_to_use, params, false)
}

/// The CXSMARTS for a fragment. Mirrors `rdmolfiles.MolFragmentToCXSmarts`.
pub fn mol_fragment_to_cx_smarts(mol: &Mol, atoms_to_use: &[u32]) -> Result<String> {
    mol_fragment_to_cx_smarts_with_params(mol, atoms_to_use, &MolFragmentToSmartsParams::default())
}

/// Mirrors `rdmolfiles.MolFragmentToCXSmarts(mol, atomsToUse, bondsToUse,
/// isomericSmarts)`.
pub fn mol_fragment_to_cx_smarts_with_params(
    mol: &Mol,
    atoms_to_use: &[u32],
    params: &MolFragmentToSmartsParams,
) -> Result<String> {
    fragment_to_smarts(mol, atoms_to_use, params, true)
}

/// Parameters for `mol_to_random_smiles_vect_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MolToRandomSmilesVectParams {
    /// Seeds RDKit's global generator when positive.
    pub random_seed: u32,
    pub isomeric_smiles: bool,
    pub kekule_smiles: bool,
    pub all_bonds_explicit: bool,
    pub all_hs_explicit: bool,
}

impl Default for MolToRandomSmilesVectParams {
    fn default() -> MolToRandomSmilesVectParams {
        MolToRandomSmilesVectParams {
            random_seed: 0,
            isomeric_smiles: true,
            kekule_smiles: false,
            all_bonds_explicit: false,
            all_hs_explicit: false,
        }
    }
}

/// `num_smiles` SMILES for a molecule from random traversals. RDKit draws
/// them from one process-wide generator, which the call holds for its
/// duration, so a seeded call repeats. Mirrors
/// `rdmolfiles.MolToRandomSmilesVect`.
pub fn mol_to_random_smiles_vect(mol: &Mol, num_smiles: u32) -> Result<Vec<String>> {
    mol_to_random_smiles_vect_with_params(mol, num_smiles, &MolToRandomSmilesVectParams::default())
}

/// RDKit converts the seed to `int`, so a seed above `i32::MAX`, like 0, does
/// not reseed the generator. Mirrors `rdmolfiles.MolToRandomSmilesVect` with
/// its optional arguments.
pub fn mol_to_random_smiles_vect_with_params(
    mol: &Mol,
    num_smiles: u32,
    params: &MolToRandomSmilesVectParams,
) -> Result<Vec<String>> {
    ffi::mol_to_random_smiles_vect(
        mol.as_ffi(),
        num_smiles,
        params.random_seed,
        params.isomeric_smiles,
        params.kekule_smiles,
        params.all_bonds_explicit,
        params.all_hs_explicit,
    )?
    .iter()
    .map(|smiles| {
        smiles
            .to_str()
            .map(str::to_owned)
            .map_err(|err| crate::Error::Value {
                message: format!("invalid UTF-8 in a string from RDKit: {err}"),
            })
    })
    .collect()
}

/// Parameters for `canonical_rank_atoms_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CanonicalRankAtomsParams {
    /// Breaks ties; otherwise equal ranks mark symmetry classes.
    pub break_ties: bool,
    pub include_chirality: bool,
    pub include_isotopes: bool,
    pub include_atom_maps: bool,
    pub include_chiral_presence: bool,
}

impl Default for CanonicalRankAtomsParams {
    fn default() -> CanonicalRankAtomsParams {
        CanonicalRankAtomsParams {
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
        }
    }
}

/// Each atom's canonical rank. A molecule without a property cache is
/// `Error::Runtime`. Mirrors `rdmolfiles.CanonicalRankAtoms`.
pub fn canonical_rank_atoms(mol: &Mol) -> Result<Vec<u32>> {
    canonical_rank_atoms_with_params(mol, &CanonicalRankAtomsParams::default())
}

/// Mirrors `rdmolfiles.CanonicalRankAtoms(mol, breakTies, includeChirality,
/// includeIsotopes, includeAtomMaps, includeChiralPresence)`.
pub fn canonical_rank_atoms_with_params(
    mol: &Mol,
    params: &CanonicalRankAtomsParams,
) -> Result<Vec<u32>> {
    Ok(ffi::canonical_rank_atoms(
        mol.as_ffi(),
        params.break_ties,
        params.include_chirality,
        params.include_isotopes,
        params.include_atom_maps,
        params.include_chiral_presence,
    )?)
}

/// Parameters for `canonical_rank_atoms_in_fragment_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CanonicalRankAtomsInFragmentParams {
    /// The bonds to use; `None` for none.
    pub bonds_to_use: Option<Vec<u32>>,
    /// A symbol for each atom of the molecule.
    pub atom_symbols: Option<Vec<String>>,
    pub break_ties: bool,
    pub include_chirality: bool,
    pub include_isotopes: bool,
    pub include_atom_maps: bool,
    pub include_chiral_presence: bool,
}

impl Default for CanonicalRankAtomsInFragmentParams {
    fn default() -> CanonicalRankAtomsInFragmentParams {
        CanonicalRankAtomsInFragmentParams {
            bonds_to_use: None,
            atom_symbols: None,
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
        }
    }
}

/// Each atom's canonical rank within the fragment made of `atoms_to_use`,
/// and -1 for atoms outside it. Mirrors
/// `rdmolfiles.CanonicalRankAtomsInFragment`.
pub fn canonical_rank_atoms_in_fragment(mol: &Mol, atoms_to_use: &[u32]) -> Result<Vec<i32>> {
    canonical_rank_atoms_in_fragment_with_params(
        mol,
        atoms_to_use,
        &CanonicalRankAtomsInFragmentParams::default(),
    )
}

/// Mirrors `rdmolfiles.CanonicalRankAtomsInFragment` with its optional
/// arguments.
pub fn canonical_rank_atoms_in_fragment_with_params(
    mol: &Mol,
    atoms_to_use: &[u32],
    params: &CanonicalRankAtomsInFragmentParams,
) -> Result<Vec<i32>> {
    let none = Vec::new();
    Ok(ffi::canonical_rank_atoms_in_fragment(
        mol.as_ffi(),
        atoms_to_use,
        params.bonds_to_use.is_some(),
        params.bonds_to_use.as_deref().unwrap_or_default(),
        params.atom_symbols.is_some(),
        params.atom_symbols.as_ref().unwrap_or(&none),
        params.break_ties,
        params.include_chirality,
        params.include_isotopes,
        params.include_atom_maps,
        params.include_chiral_presence,
    )?)
}

/// Orders the atoms of each enhanced stereo group canonically. Mirrors
/// `rdmolfiles.CanonicalizeEnhancedStereo`.
pub fn canonicalize_enhanced_stereo(mol: &mut Mol) -> Result<()> {
    Ok(ffi::canonicalize_enhanced_stereo(mol.pin_mut())?)
}
