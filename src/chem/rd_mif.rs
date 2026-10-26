// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMIF`: molecular interaction fields, the energies
//! of a probe with a molecule at the points of a grid.
//!
//! ```no_run
//! use rdkit_cxx::chem::rd_mif::{Coulomb, calculate_descriptors, construct_grid};
//! use rdkit_cxx::chem::rd_partial_charges::compute_gasteiger_charges;
//! use rdkit_cxx::chem::mol_from_mol_file;
//!
//! let mut mol = mol_from_mol_file("ligand.mol").unwrap();
//! compute_gasteiger_charges(&mut mol).unwrap();
//! let mut grid = construct_grid(&mol).unwrap();
//! calculate_descriptors(&mut grid, &Coulomb::new(&mol).unwrap());
//! ```
//!
//! Each functor computes from a copy of the conformer's coordinates, so it
//! does not borrow the molecule. RDKit trusts some inputs that the wrappers
//! check first, refusing with `Error::Value`: probe atom types without
//! force-field parameters (`MMFFVdWaals`, `UFFVdWaals`), charges and
//! positions of different lengths, nitrogens and oxygens with more than 12
//! neighbors (`HBond`, `Hydrophilic`), grids too large for `i32`, and cube
//! files with values that do not parse or are missing, or a negative atom
//! count. `CoulombDielectric` for a molecule without atoms, which never
//! returns in Python, gives a field of 0.
//!
//! `CoulombDielectric`, `HBond`, and `Hydrophilic` compute in buffers they
//! own while they are called, so they are `Send` but not `Sync`.

use std::path::Path;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::mif::ffi;
use crate::geometry::{Point3D, UniformRealValueGrid3D};
use crate::util::path_bytes;

/// Parameters for `construct_grid_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct ConstructGridParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// The least distance from an atom to the grid's faces, in Å.
    pub margin: f64,
    /// The grid's spacing, in Å.
    pub spacing: f64,
}

impl Default for ConstructGridParams {
    fn default() -> ConstructGridParams {
        ConstructGridParams {
            conf_id: -1,
            margin: 5.0,
            spacing: 0.5,
        }
    }
}

/// A grid around the molecule's first conformer, 5 Å from every atom, of
/// spacing 0.5 Å. Mirrors `rdMIF.ConstructGrid`.
///
/// A molecule without conformers is `Error::Runtime`, a missing conformer
/// `Error::Value`, and so is a grid of 2^31 points or more, or with none
/// along an axis.
pub fn construct_grid(mol: &Mol) -> Result<UniformRealValueGrid3D> {
    construct_grid_with_params(mol, &ConstructGridParams::default())
}

/// Mirrors `rdMIF.ConstructGrid(mol, confId, margin, spacing)`.
pub fn construct_grid_with_params(
    mol: &Mol,
    params: &ConstructGridParams,
) -> Result<UniformRealValueGrid3D> {
    Ok(UniformRealValueGrid3D::from_ffi(ffi::construct_grid(
        mol.as_ffi(),
        params.conf_id,
        params.margin,
        params.spacing,
    )?))
}

/// Defines a functor's `call` and its `calculate_descriptors` functions.
macro_rules! functor {
    ($name:ident, $python:literal, $call:ident, $calc:ident, $calc_threshold:ident, $ffi_calc:ident) => {
        impl $name {
            /// The energy of the probe at (x, y, z), in kJ/mol, from the
            /// atoms whose squared distance from it, in Å², is less than
            #[doc = concat!(
                "`threshold`, which RDKit compares with squared distances. Mirrors `",
                $python, ".__call__`."
            )]
            pub fn call(&self, x: f64, y: f64, z: f64, threshold: f64) -> f64 {
                expect_ok(ffi::$call(&self.ptr, x, y, z, threshold))
            }
        }

        #[doc = concat!(
            "Sets every point of the grid to the descriptor's energy there. Mirrors `rdMIF.CalculateDescriptors(grid, ",
            $python, ")`."
        )]
        pub fn $calc(grid: &mut UniformRealValueGrid3D, descriptor: &$name) {
            $calc_threshold(grid, descriptor, -1.0)
        }

        /// Mirrors `rdMIF.CalculateDescriptors(grid, descriptor,
        /// threshold)`: atoms `threshold` Å or more from a point add nothing
        /// to it; a negative threshold is the grid's spacing times its
        /// number of points.
        pub fn $calc_threshold(grid: &mut UniformRealValueGrid3D, descriptor: &$name, threshold: f64) {
            expect_ok(ffi::$ffi_calc(grid.raw_mut(), &descriptor.ptr, threshold))
        }
    };
}

/// The Coulomb energy of a probe charge with a molecule's partial charges,
/// in vacuum. Mirrors `rdMIF.Coulomb`.
pub struct Coulomb {
    ptr: UniquePtr<ffi::Coulomb>,
}

impl std::fmt::Debug for Coulomb {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Coulomb").finish_non_exhaustive()
    }
}

/// Parameters for `Coulomb::new_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CoulombParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// The probe's charge, in e.
    pub probe_charge: f64,
    /// Gives each energy as minus its absolute value.
    pub abs_val: bool,
    /// The atom property that holds each atom's partial charge.
    pub charge_key: String,
    /// The softcore parameter, in Å²; 0 for a least distance of `cutoff`.
    pub softcore_param: f64,
    /// The least distance, in Å.
    pub cutoff: f64,
}

impl Default for CoulombParams {
    fn default() -> CoulombParams {
        CoulombParams {
            conf_id: -1,
            probe_charge: 1.0,
            abs_val: false,
            charge_key: "_GasteigerCharge".to_string(),
            softcore_param: 0.0,
            cutoff: 1.0,
        }
    }
}

/// Parameters for `Coulomb::from_charges_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CoulombFromChargesParams {
    /// The probe's charge, in e.
    pub probe_charge: f64,
    /// Gives each energy as minus its absolute value.
    pub abs_val: bool,
    /// The softcore parameter, in Å²; 0 for a least distance of `cutoff`.
    pub softcore_param: f64,
    /// The least distance, in Å.
    pub cutoff: f64,
}

impl Default for CoulombFromChargesParams {
    fn default() -> CoulombFromChargesParams {
        CoulombFromChargesParams {
            probe_charge: 1.0,
            abs_val: false,
            softcore_param: 0.0,
            cutoff: 1.0,
        }
    }
}

impl Coulomb {
    /// The field of the molecule's first conformer, with the partial charges
    /// in its atoms' `_GasteigerCharge` properties. Mirrors
    /// `Coulomb(mol)`.
    ///
    /// A missing conformer is `Error::Value`, and a missing charge
    /// `Error::Key`.
    pub fn new(mol: &Mol) -> Result<Coulomb> {
        Coulomb::new_with_params(mol, &CoulombParams::default())
    }

    /// Mirrors `Coulomb(mol, confId, probeCharge, absVal, chargeKey,
    /// softcoreParam, cutoff)`.
    pub fn new_with_params(mol: &Mol, params: &CoulombParams) -> Result<Coulomb> {
        Ok(Coulomb {
            ptr: ffi::coulomb_new(
                mol.as_ffi(),
                params.conf_id,
                params.probe_charge,
                params.abs_val,
                &params.charge_key,
                params.softcore_param,
                params.cutoff,
            )?,
        })
    }

    /// The field of the given charges at the given positions; different
    /// lengths are `Error::Value`. Mirrors `Coulomb(charges, positions)`.
    pub fn from_charges(charges: &[f64], positions: &[Point3D]) -> Result<Coulomb> {
        Coulomb::from_charges_with_params(charges, positions, &CoulombFromChargesParams::default())
    }

    /// Mirrors `Coulomb(charges, positions, probeCharge, absVal,
    /// softcoreParam, cutoff)`.
    pub fn from_charges_with_params(
        charges: &[f64],
        positions: &[Point3D],
        params: &CoulombFromChargesParams,
    ) -> Result<Coulomb> {
        Ok(Coulomb {
            ptr: ffi::coulomb_from_charges(
                charges,
                positions,
                params.probe_charge,
                params.abs_val,
                params.softcore_param,
                params.cutoff,
            )?,
        })
    }
}

functor!(
    Coulomb,
    "Coulomb",
    coulomb_call,
    calculate_descriptors,
    calculate_descriptors_with_threshold,
    calculate_descriptors_coulomb
);

/// The Coulomb energy of a probe charge with a molecule's partial charges,
/// with a distance-dependent dielectric, as GRID computes it. Mirrors
/// `rdMIF.CoulombDielectric`.
pub struct CoulombDielectric {
    ptr: UniquePtr<ffi::CoulombDielectric>,
}

impl std::fmt::Debug for CoulombDielectric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CoulombDielectric").finish_non_exhaustive()
    }
}

/// Parameters for `CoulombDielectric::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CoulombDielectricParams {
    /// The shared parameters.
    pub coulomb: CoulombParams,
    /// The solvent's relative permittivity.
    pub epsilon: f64,
    /// The solute's relative permittivity.
    pub xi: f64,
}

impl Default for CoulombDielectricParams {
    fn default() -> CoulombDielectricParams {
        CoulombDielectricParams {
            coulomb: CoulombParams::default(),
            epsilon: 80.0,
            xi: 4.0,
        }
    }
}

/// Parameters for `CoulombDielectric::from_charges_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CoulombDielectricFromChargesParams {
    /// The shared parameters.
    pub coulomb: CoulombFromChargesParams,
    /// The solvent's relative permittivity.
    pub epsilon: f64,
    /// The solute's relative permittivity.
    pub xi: f64,
}

impl Default for CoulombDielectricFromChargesParams {
    fn default() -> CoulombDielectricFromChargesParams {
        CoulombDielectricFromChargesParams {
            coulomb: CoulombFromChargesParams::default(),
            epsilon: 80.0,
            xi: 4.0,
        }
    }
}

impl CoulombDielectric {
    /// The field of the molecule's first conformer, with the partial charges
    /// in its atoms' `_GasteigerCharge` properties. Mirrors
    /// `CoulombDielectric(mol)`.
    ///
    /// A molecule without conformers is `Error::Runtime`, a missing
    /// conformer `Error::Value`, and a missing charge `Error::Key`.
    pub fn new(mol: &Mol) -> Result<CoulombDielectric> {
        CoulombDielectric::new_with_params(mol, &CoulombDielectricParams::default())
    }

    /// Mirrors `CoulombDielectric(mol, confId, probeCharge, absVal,
    /// chargeKey, softcoreParam, cutoff, epsilon, xi)`.
    pub fn new_with_params(
        mol: &Mol,
        params: &CoulombDielectricParams,
    ) -> Result<CoulombDielectric> {
        let c = &params.coulomb;
        Ok(CoulombDielectric {
            ptr: ffi::coulomb_dielectric_new(
                mol.as_ffi(),
                c.conf_id,
                c.probe_charge,
                c.abs_val,
                &c.charge_key,
                c.softcore_param,
                c.cutoff,
                params.epsilon,
                params.xi,
            )?,
        })
    }

    /// The field of the given charges at the given positions; different
    /// lengths are `Error::Value`. Mirrors `CoulombDielectric(charges,
    /// positions)`.
    pub fn from_charges(charges: &[f64], positions: &[Point3D]) -> Result<CoulombDielectric> {
        CoulombDielectric::from_charges_with_params(
            charges,
            positions,
            &CoulombDielectricFromChargesParams::default(),
        )
    }

    /// Mirrors `CoulombDielectric(charges, positions, probeCharge, absVal,
    /// softcoreParam, cutoff, epsilon, xi)`.
    pub fn from_charges_with_params(
        charges: &[f64],
        positions: &[Point3D],
        params: &CoulombDielectricFromChargesParams,
    ) -> Result<CoulombDielectric> {
        let c = &params.coulomb;
        Ok(CoulombDielectric {
            ptr: ffi::coulomb_dielectric_from_charges(
                charges,
                positions,
                c.probe_charge,
                c.abs_val,
                c.softcore_param,
                c.cutoff,
                params.epsilon,
                params.xi,
            )?,
        })
    }
}

functor!(
    CoulombDielectric,
    "CoulombDielectric",
    coulomb_dielectric_call,
    calculate_descriptors_coulomb_dielectric,
    calculate_descriptors_coulomb_dielectric_with_threshold,
    calculate_descriptors_coulomb_dielectric
);

/// The van der Waals energy of a probe atom with a molecule, by MMFF94.
/// Mirrors `rdMIF.MMFFVdWaals`.
pub struct MMFFVdWaals {
    ptr: UniquePtr<ffi::MMFFVdWaals>,
}

impl std::fmt::Debug for MMFFVdWaals {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MMFFVdWaals").finish_non_exhaustive()
    }
}

/// Parameters for `MMFFVdWaals::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MMFFVdWaalsParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// The probe's MMFF94 atom type; 6 is an sp3 oxygen.
    pub probe_atom_type: u32,
    /// Scales the parameters to account for hydrogen bonds.
    pub scaling: bool,
    /// The least distance, in Å.
    pub cutoff: f64,
}

impl Default for MMFFVdWaalsParams {
    fn default() -> MMFFVdWaalsParams {
        MMFFVdWaalsParams {
            conf_id: -1,
            probe_atom_type: 6,
            scaling: false,
            cutoff: 1.0,
        }
    }
}

impl MMFFVdWaals {
    /// The field of the molecule's first conformer for an sp3 oxygen probe.
    /// Mirrors `MMFFVdWaals(mol)`.
    ///
    /// A missing conformer, an atom MMFF cannot type, and a probe type
    /// without MMFF parameters (which RDKit dereferences as null) are
    /// `Error::Value`.
    pub fn new(mol: &Mol) -> Result<MMFFVdWaals> {
        MMFFVdWaals::new_with_params(mol, &MMFFVdWaalsParams::default())
    }

    /// Mirrors `MMFFVdWaals(mol, confId, probeAtomType, scaling, cutoff)`.
    pub fn new_with_params(mol: &Mol, params: &MMFFVdWaalsParams) -> Result<MMFFVdWaals> {
        Ok(MMFFVdWaals {
            ptr: ffi::mmff_vd_waals_new(
                mol.as_ffi(),
                params.conf_id,
                params.probe_atom_type,
                params.scaling,
                params.cutoff,
            )?,
        })
    }
}

functor!(
    MMFFVdWaals,
    "MMFFVdWaals",
    mmff_vd_waals_call,
    calculate_descriptors_mmff_vd_waals,
    calculate_descriptors_mmff_vd_waals_with_threshold,
    calculate_descriptors_mmff_vd_waals
);

/// The van der Waals energy of a probe atom with a molecule, by UFF.
/// Mirrors `rdMIF.UFFVdWaals`.
pub struct UFFVdWaals {
    ptr: UniquePtr<ffi::UFFVdWaals>,
}

impl std::fmt::Debug for UFFVdWaals {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UFFVdWaals").finish_non_exhaustive()
    }
}

/// Parameters for `UFFVdWaals::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UFFVdWaalsParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// The probe's UFF atom type; `O_3` is an sp3 oxygen.
    pub probe_atom_type: String,
    /// The least distance, in Å.
    pub cutoff: f64,
}

impl Default for UFFVdWaalsParams {
    fn default() -> UFFVdWaalsParams {
        UFFVdWaalsParams {
            conf_id: -1,
            probe_atom_type: "O_3".to_string(),
            cutoff: 1.0,
        }
    }
}

impl UFFVdWaals {
    /// The field of the molecule's first conformer for an sp3 oxygen probe.
    /// Mirrors `UFFVdWaals(mol)`.
    ///
    /// A missing conformer, an atom UFF has no parameters for, and a probe
    /// type UFF lacks (which RDKit dereferences as null) are
    /// `Error::Value`.
    pub fn new(mol: &Mol) -> Result<UFFVdWaals> {
        UFFVdWaals::new_with_params(mol, &UFFVdWaalsParams::default())
    }

    /// Mirrors `UFFVdWaals(mol, confId, probeAtomType, cutoff)`.
    pub fn new_with_params(mol: &Mol, params: &UFFVdWaalsParams) -> Result<UFFVdWaals> {
        Ok(UFFVdWaals {
            ptr: ffi::uff_vd_waals_new(
                mol.as_ffi(),
                params.conf_id,
                &params.probe_atom_type,
                params.cutoff,
            )?,
        })
    }
}

functor!(
    UFFVdWaals,
    "UFFVdWaals",
    uff_vd_waals_call,
    calculate_descriptors_uff_vd_waals,
    calculate_descriptors_uff_vd_waals_with_threshold,
    calculate_descriptors_uff_vd_waals
);

/// The hydrogen-bonding energy of a probe with a molecule, as GRID computes
/// it. Mirrors `rdMIF.HBond`.
pub struct HBond {
    ptr: UniquePtr<ffi::HBond>,
}

impl std::fmt::Debug for HBond {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HBond").finish_non_exhaustive()
    }
}

/// Parameters for `HBond::new_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct HBondParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// The probe: `OH` or `NH`, which donate to the molecule's acceptors,
    /// or `O` or `N`, which accept from its donors.
    pub probe_atom_type: String,
    /// Uses GRID's angular terms for rigid groups, rather than those that
    /// let hydrogens and lone pairs rotate.
    pub fixed: bool,
    /// The least distance, in Å.
    pub cutoff: f64,
}

impl Default for HBondParams {
    fn default() -> HBondParams {
        HBondParams {
            conf_id: -1,
            probe_atom_type: "OH".to_string(),
            fixed: true,
            cutoff: 1.0,
        }
    }
}

impl HBond {
    /// The field of the molecule's first conformer for an OH probe. Mirrors
    /// `HBond(mol)`.
    ///
    /// A missing conformer, another probe type, and a nitrogen or oxygen
    /// with more than 12 neighbors (whose bonds RDKit writes past an array
    /// of 12) are `Error::Value`.
    pub fn new(mol: &Mol) -> Result<HBond> {
        HBond::new_with_params(mol, &HBondParams::default())
    }

    /// Mirrors `HBond(mol, confId, probeAtomType, fixed, cutoff)`.
    pub fn new_with_params(mol: &Mol, params: &HBondParams) -> Result<HBond> {
        Ok(HBond {
            ptr: ffi::h_bond_new(
                mol.as_ffi(),
                params.conf_id,
                &params.probe_atom_type,
                params.fixed,
                params.cutoff,
            )?,
        })
    }
}

functor!(
    HBond,
    "HBond",
    h_bond_call,
    calculate_descriptors_h_bond,
    calculate_descriptors_h_bond_with_threshold,
    calculate_descriptors_h_bond
);

/// The more favorable of the hydrogen-bonding energies of water's hydrogen
/// (an OH probe) and oxygen (an O probe) with a molecule. Mirrors
/// `rdMIF.Hydrophilic`.
pub struct Hydrophilic {
    ptr: UniquePtr<ffi::Hydrophilic>,
}

impl std::fmt::Debug for Hydrophilic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Hydrophilic").finish_non_exhaustive()
    }
}

/// Parameters for `Hydrophilic::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct HydrophilicParams {
    /// The conformer; -1 for the first.
    pub conf_id: i32,
    /// Uses GRID's angular terms for rigid groups.
    pub fixed: bool,
    /// The least distance, in Å.
    pub cutoff: f64,
}

impl Default for HydrophilicParams {
    fn default() -> HydrophilicParams {
        HydrophilicParams {
            conf_id: -1,
            fixed: true,
            cutoff: 1.0,
        }
    }
}

impl Hydrophilic {
    /// The field of the molecule's first conformer. Mirrors
    /// `Hydrophilic(mol)`, with the errors of `HBond::new`.
    pub fn new(mol: &Mol) -> Result<Hydrophilic> {
        Hydrophilic::new_with_params(mol, &HydrophilicParams::default())
    }

    /// Mirrors `Hydrophilic(mol, confId, fixed, cutoff)`.
    pub fn new_with_params(mol: &Mol, params: &HydrophilicParams) -> Result<Hydrophilic> {
        Ok(Hydrophilic {
            ptr: ffi::hydrophilic_new(mol.as_ffi(), params.conf_id, params.fixed, params.cutoff)?,
        })
    }
}

functor!(
    Hydrophilic,
    "Hydrophilic",
    hydrophilic_call,
    calculate_descriptors_hydrophilic,
    calculate_descriptors_hydrophilic_with_threshold,
    calculate_descriptors_hydrophilic
);

/// Parameters for `write_to_cube_file_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug)]
pub struct WriteToCubeFileParams<'a> {
    /// A molecule whose atoms the file records; `None` for none.
    pub mol: Option<&'a Mol>,
    /// The molecule's conformer; -1 for the first.
    pub conf_id: i32,
}

impl Default for WriteToCubeFileParams<'_> {
    fn default() -> Self {
        WriteToCubeFileParams {
            mol: None,
            conf_id: -1,
        }
    }
}

/// Writes the grid to a file in Gaussian's cube format. Mirrors
/// `rdMIF.WriteToCubeFile`.
///
/// A file that cannot be written is `Error::Os`.
pub fn write_to_cube_file(grid: &UniformRealValueGrid3D, filename: impl AsRef<Path>) -> Result<()> {
    write_to_cube_file_with_params(grid, filename, &WriteToCubeFileParams::default())
}

/// Mirrors `rdMIF.WriteToCubeFile(grid, filename, mol, confId)`; a missing
/// conformer is `Error::Value`.
pub fn write_to_cube_file_with_params(
    grid: &UniformRealValueGrid3D,
    filename: impl AsRef<Path>,
    params: &WriteToCubeFileParams<'_>,
) -> Result<()> {
    let mol = params
        .mol
        .map_or(std::ptr::null(), |mol| mol.as_ffi() as *const _);
    // SAFETY: the molecule is null or borrowed for the call.
    Ok(unsafe {
        ffi::write_to_cube_file(
            grid.raw(),
            path_bytes(filename.as_ref()),
            mol,
            params.conf_id,
        )
    }?)
}

/// Reads a grid from a file in Gaussian's cube format, with the molecule it
/// records, if any: its atoms and coordinates, without bonds. Mirrors
/// `rdMIF.ReadFromCubeFile`.
///
/// A file that cannot be opened is `Error::Os`. Spacings that differ are
/// `Error::Parse`; a value that does not parse, too few values, a negative
/// atom count, and a grid of 2^31 points or more are `Error::Value`; a grid
/// with no points along an axis is `Error::Runtime`.
pub fn read_from_cube_file(
    filename: impl AsRef<Path>,
) -> Result<(UniformRealValueGrid3D, Option<Mol>)> {
    let mut grid = UniformRealValueGrid3D::default();
    let mol = ffi::read_from_cube_file(path_bytes(filename.as_ref()), grid.raw_mut())?;
    Ok((grid, (!mol.is_null()).then(|| Mol::from_ffi(mol))))
}
