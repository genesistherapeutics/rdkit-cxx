// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdmolops`: hydrogens, sanitization, kekulization,
//! aromaticity, fragments, rings, stereochemistry, query adjustment, and
//! the classic RDKit, layered, and pattern fingerprints.
//!
//! Operations that change a molecule in place take `&mut Mol`, as
//! [`sanitize_mol`] does; those that return a new molecule take `&Mol`, as
//! [`add_hs`] does.
//!
//! ```no_run
//! use rdkit_cxx::chem::{add_hs, get_mol_frags, mol_from_smiles, remove_hs};
//!
//! let mol = mol_from_smiles("CC(=O)[O-].[Na+]")?;
//! assert_eq!(get_mol_frags(&mol), [vec![0, 1, 2, 3], vec![4]]);
//! let with_hs = add_hs(&mol)?;
//! assert_eq!(remove_hs(&with_hs)?.get_num_atoms(), mol.get_num_atoms());
//! # Ok::<(), rdkit_cxx::Error>(())
//! ```

mod adjust_query;
mod fingerprints;
mod graph;
mod hs;
mod perception;
mod stereo;
mod transforms;

use std::collections::BTreeMap;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::chem::rdchem::StereoInfo;
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::molops::ffi;
use crate::{Error, Result};

pub use adjust_query::*;
pub use fingerprints::*;
pub use graph::*;
pub use hs::*;
pub use perception::*;
pub use stereo::*;
pub use transforms::*;

/// Sanitization operations, as bits. Mirrors `rdmolops.SanitizeFlags`
/// (CONVENTIONS.md F1).
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct SanitizeFlags(u32);

macro_rules! sanitize_flags {
    ($($name:ident),* $(,)?) => {
        impl SanitizeFlags {
            $(
                #[doc = concat!("Mirrors `SanitizeFlags.", stringify!($name), "`.")]
                pub const $name: SanitizeFlags = SanitizeFlags(ffi::SanitizeFlags::$name.repr);
            )*

            const NAMES: &[(&str, SanitizeFlags)] = &[$((stringify!($name), SanitizeFlags::$name)),*];
        }
    };
}

sanitize_flags!(
    SANITIZE_NONE,
    SANITIZE_CLEANUP,
    SANITIZE_PROPERTIES,
    SANITIZE_SYMMRINGS,
    SANITIZE_KEKULIZE,
    SANITIZE_FINDRADICALS,
    SANITIZE_SETAROMATICITY,
    SANITIZE_SETCONJUGATION,
    SANITIZE_SETHYBRIDIZATION,
    SANITIZE_CLEANUPCHIRALITY,
    SANITIZE_ADJUSTHS,
    SANITIZE_CLEANUP_ORGANOMETALLICS,
    SANITIZE_CLEANUPATROPISOMERS,
    SANITIZE_ALL,
);

crate::flags::flag_ops!(SanitizeFlags);

/// Parameters for `sanitize_mol_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SanitizeMolParams {
    /// The operations to run.
    pub sanitize_ops: SanitizeFlags,
    /// Returns the operation that failed instead of an error.
    pub catch_errors: bool,
}

impl Default for SanitizeMolParams {
    fn default() -> SanitizeMolParams {
        SanitizeMolParams {
            sanitize_ops: SanitizeFlags::SANITIZE_ALL,
            catch_errors: false,
        }
    }
}

/// Parameters for `add_hs_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AddHsParams {
    /// Adds only the explicit hydrogens.
    pub explicit_only: bool,
    /// Gives the hydrogens coordinates.
    pub add_coords: bool,
    /// The atoms to add hydrogens to; `None` for all, and an empty list for
    /// none (Python treats an empty list as `None`).
    pub only_on_atoms: Option<Vec<u32>>,
    /// Gives the hydrogens PDB residue information.
    pub add_residue_info: bool,
}

impl AddHsParams {
    fn to_ffi(&self) -> ffi::AddHsParametersData {
        ffi::AddHsParametersData {
            explicit_only: self.explicit_only,
            add_coords: self.add_coords,
            add_residue_info: self.add_residue_info,
            has_only_on_atoms: self.only_on_atoms.is_some(),
            only_on_atoms: self.only_on_atoms.clone().unwrap_or_default(),
        }
    }
}

/// Sanitizes a molecule in place, returning `SANITIZE_NONE`. Mirrors
/// `rdmolops.SanitizeMol`.
pub fn sanitize_mol(mol: &mut Mol) -> Result<SanitizeFlags> {
    sanitize_mol_with_params(mol, &SanitizeMolParams::default())
}

/// Runs the chosen operations. With `catch_errors`, a failure returns the
/// operation that failed rather than an error. Mirrors
/// `rdmolops.SanitizeMol(mol, sanitizeOps, catchErrors)`.
pub fn sanitize_mol_with_params(
    mol: &mut Mol,
    params: &SanitizeMolParams,
) -> Result<SanitizeFlags> {
    let failed = ffi::sanitize_mol(
        mol.pin_mut(),
        params.sanitize_ops.bits(),
        params.catch_errors,
    )?;
    Ok(SanitizeFlags(failed))
}

/// Kekulizes a molecule in place. Mirrors `rdmolops.Kekulize`.
pub fn kekulize(mol: &mut Mol) -> Result<()> {
    kekulize_with_clear_aromatic_flags(mol, false)
}

/// Mirrors `rdmolops.Kekulize(mol, clearAromaticFlags)`: also clears the
/// aromatic flags on atoms and bonds.
pub fn kekulize_with_clear_aromatic_flags(mol: &mut Mol, clear_aromatic_flags: bool) -> Result<()> {
    Ok(ffi::kekulize(mol.pin_mut(), clear_aromatic_flags)?)
}

/// A copy of a molecule with its hydrogens made explicit. Mirrors
/// `rdmolops.AddHs`.
pub fn add_hs(mol: &Mol) -> Result<Mol> {
    add_hs_with_params(mol, &AddHsParams::default())
}

/// Mirrors `rdmolops.AddHs(mol, explicitOnly, addCoords, onlyOnAtoms,
/// addResidueInfo)`.
pub fn add_hs_with_params(mol: &Mol, params: &AddHsParams) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::add_hs(mol.as_ffi(), &params.to_ffi())?))
}

/// The atom indices of each disconnected fragment. Mirrors
/// `rdmolops.GetMolFrags`. Python's `frags` and `fragsMolAtomMapping` output
/// arguments are not bound (A4).
pub fn get_mol_frags(mol: &Mol) -> Vec<Vec<u32>> {
    IndexList::into_nested(expect_ok(ffi::get_mol_frags(mol.as_ffi())))
}

/// Each disconnected fragment as a sanitized molecule. Mirrors
/// `rdmolops.GetMolFrags(mol, asMols=True)` (A4).
pub fn get_mol_frags_as_mols(mol: &Mol) -> Result<Vec<Mol>> {
    get_mol_frags_as_mols_with_sanitize_frags(mol, true)
}

/// Mirrors `rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags)`.
pub fn get_mol_frags_as_mols_with_sanitize_frags(
    mol: &Mol,
    sanitize_frags: bool,
) -> Result<Vec<Mol>> {
    Ok(mol_list_into_vec(ffi::get_mol_frags_as_mols(
        mol.as_ffi(),
        sanitize_frags,
    )?))
}

/// Parameters for `find_potential_stereo_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct FindPotentialStereoParams {
    /// Removes stereo specifications that cannot hold.
    pub clean_it: bool,
    /// Includes potential stereo elements that are not specified.
    pub flag_possible: bool,
}

impl Default for FindPotentialStereoParams {
    fn default() -> FindPotentialStereoParams {
        FindPotentialStereoParams {
            clean_it: false,
            flag_possible: true,
        }
    }
}

/// The molecule's potential stereo elements. Mirrors
/// `rdmolops.FindPotentialStereo`; RDKit calls it experimental.
pub fn find_potential_stereo(mol: &mut Mol) -> Result<Vec<StereoInfo>> {
    find_potential_stereo_with_params(mol, &FindPotentialStereoParams::default())
}

/// Mirrors `rdmolops.FindPotentialStereo(mol, cleanIt, flagPossible)`.
pub fn find_potential_stereo_with_params(
    mol: &mut Mol,
    params: &FindPotentialStereoParams,
) -> Result<Vec<StereoInfo>> {
    Ok(crate::ffi::stereo::ffi::find_potential_stereo(
        mol.pin_mut(),
        params.clean_it,
        params.flag_possible,
    )?
    .into_iter()
    .map(StereoInfo::from_ffi)
    .collect())
}

/// Moves each molecule out of a C++ list (W6).
pub(crate) fn mol_list_into_vec(mut list: UniquePtr<ffi::MolList>) -> Vec<Mol> {
    let size = expect_ok(ffi::mol_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    (0..size)
        .map(|i| Mol::from_ffi(expect_ok(ffi::mol_list_take(list.pin_mut(), i))))
        .collect()
}

/// Moves each named molecule out of a C++ list (W6); invalid UTF-8 in a name
/// is `Error::Value` (T2).
pub(crate) fn mol_list_into_map(
    mut list: UniquePtr<ffi::MolList>,
) -> Result<BTreeMap<String, Mol>> {
    let size = expect_ok(ffi::mol_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    let mut map = BTreeMap::new();
    for i in 0..size {
        let name = expect_ok(ffi::mol_list_name(&list, i))
            .to_str()
            .map_err(|err| Error::Value {
                message: format!("invalid UTF-8 in a string from RDKit: {err}"),
            })?
            .to_owned();
        let mol = Mol::from_ffi(expect_ok(ffi::mol_list_take(list.pin_mut(), i)));
        map.insert(name, mol);
    }
    Ok(map)
}

/// Calls `f` with a C++ list that borrows the molecules, which RDKit may
/// modify (W8).
pub(crate) fn with_borrowed_mols<'a, T>(
    mols: impl IntoIterator<Item = (&'a str, &'a mut Mol)>,
    f: impl FnOnce(Pin<&mut ffi::BorrowedMols>) -> T,
) -> T {
    let mut list = expect_ok(ffi::borrowed_mols_new());
    for (name, mol) in mols {
        expect_ok(ffi::borrowed_mols_push(list.pin_mut(), name, mol.pin_mut()));
    }
    // The list drops here, before the borrows of the molecules end.
    f(list.pin_mut())
}

/// Calls `f` with a C++ list that borrows the molecules, for calls that only
/// read them (W8).
pub(crate) fn with_borrowed_mols_shared<'a, T>(
    mols: impl IntoIterator<Item = (&'a str, &'a Mol)>,
    f: impl FnOnce(Pin<&mut ffi::BorrowedMols>) -> T,
) -> T {
    let mut list = expect_ok(ffi::borrowed_mols_new());
    for (name, mol) in mols {
        expect_ok(ffi::borrowed_mols_push_shared(
            list.pin_mut(),
            name,
            mol.as_ffi(),
        ));
    }
    // The list drops here, before the borrows of the molecules end.
    f(list.pin_mut())
}
