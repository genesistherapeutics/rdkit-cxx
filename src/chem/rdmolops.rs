//! Mirrors `rdkit.Chem.rdmolops`.

use std::fmt;
use std::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, Not};

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::molops::ffi;

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

impl SanitizeFlags {
    /// The raw bits.
    pub const fn bits(self) -> u32 {
        self.0
    }

    /// Flags from raw bits, keeping bits no constant names.
    pub const fn from_bits_retain(bits: u32) -> SanitizeFlags {
        SanitizeFlags(bits)
    }
}

impl fmt::Debug for SanitizeFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match SanitizeFlags::NAMES.iter().find(|(_, flags)| flags == self) {
            Some((name, _)) => write!(f, "SanitizeFlags::{name}"),
            None => write!(f, "SanitizeFlags({:#x})", self.0),
        }
    }
}

impl BitOr for SanitizeFlags {
    type Output = SanitizeFlags;

    fn bitor(self, rhs: SanitizeFlags) -> SanitizeFlags {
        SanitizeFlags(self.0 | rhs.0)
    }
}

impl BitOrAssign for SanitizeFlags {
    fn bitor_assign(&mut self, rhs: SanitizeFlags) {
        self.0 |= rhs.0;
    }
}

impl BitAnd for SanitizeFlags {
    type Output = SanitizeFlags;

    fn bitand(self, rhs: SanitizeFlags) -> SanitizeFlags {
        SanitizeFlags(self.0 & rhs.0)
    }
}

impl BitAndAssign for SanitizeFlags {
    fn bitand_assign(&mut self, rhs: SanitizeFlags) {
        self.0 &= rhs.0;
    }
}

impl Not for SanitizeFlags {
    type Output = SanitizeFlags;

    fn not(self) -> SanitizeFlags {
        SanitizeFlags(!self.0)
    }
}

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

/// Parameters for `kekulize_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct KekulizeParams {
    /// Clears the aromatic flags on atoms and bonds.
    pub clear_aromatic_flags: bool,
    /// Picks the canonical Kekulé form.
    pub canonical: bool,
}

impl Default for KekulizeParams {
    fn default() -> KekulizeParams {
        KekulizeParams {
            clear_aromatic_flags: false,
            canonical: true,
        }
    }
}

/// Parameters for `add_hs_with_params`. Mirrors `rdmolops.AddHsParameters`,
/// plus `AddHs`'s `onlyOnAtoms` (A3).
#[derive(Clone, Debug, PartialEq)]
pub struct AddHsParameters {
    /// Adds only the explicit hydrogens.
    pub explicit_only: bool,
    /// Gives the hydrogens coordinates.
    pub add_coords: bool,
    /// Gives the hydrogens PDB residue information.
    pub add_residue_info: bool,
    /// Skips query atoms and atoms with query bonds.
    pub skip_queries: bool,
    /// The atoms to add hydrogens to; `None` for all.
    pub only_on_atoms: Option<Vec<u32>>,
}

impl Default for AddHsParameters {
    fn default() -> AddHsParameters {
        let data = expect_ok(ffi::add_hs_parameters_default());
        AddHsParameters {
            explicit_only: data.explicit_only,
            add_coords: data.add_coords,
            add_residue_info: data.add_residue_info,
            skip_queries: data.skip_queries,
            only_on_atoms: None,
        }
    }
}

impl AddHsParameters {
    fn to_ffi(&self) -> ffi::AddHsParametersData {
        ffi::AddHsParametersData {
            explicit_only: self.explicit_only,
            add_coords: self.add_coords,
            add_residue_info: self.add_residue_info,
            skip_queries: self.skip_queries,
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
    kekulize_with_params(mol, &KekulizeParams::default())
}

/// Mirrors `rdmolops.Kekulize(mol, clearAromaticFlags, canonical)`.
pub fn kekulize_with_params(mol: &mut Mol, params: &KekulizeParams) -> Result<()> {
    Ok(ffi::kekulize(
        mol.pin_mut(),
        params.clear_aromatic_flags,
        params.canonical,
    )?)
}

/// A copy of a molecule with its hydrogens made explicit. Mirrors
/// `rdmolops.AddHs`.
pub fn add_hs(mol: &Mol) -> Result<Mol> {
    add_hs_with_params(mol, &AddHsParameters::default())
}

/// Mirrors `rdmolops.AddHs(mol, params, onlyOnAtoms)` and its keyword form.
pub fn add_hs_with_params(mol: &Mol, params: &AddHsParameters) -> Result<Mol> {
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

/// Moves each molecule out of a C++ list (W6).
pub(crate) fn mol_list_into_vec(mut list: UniquePtr<ffi::MolList>) -> Vec<Mol> {
    let size = expect_ok(ffi::mol_list_size(
        list.as_ref().expect("RDKit returned a null list"),
    ));
    (0..size)
        .map(|i| Mol::from_ffi(expect_ok(ffi::mol_list_take(list.pin_mut(), i))))
        .collect()
}
