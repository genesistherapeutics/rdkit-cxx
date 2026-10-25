// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdForceFieldHelpers`: UFF and MMFF force fields for
//! molecules, and minimization of their conformers.
//!
//! The `..._confs` functions take a thread count as Python's do: 0 uses
//! every hardware thread, and a negative count that many fewer.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdmolops::SanitizeFlags;
use crate::ffi::force_field::ffi::TermParams;
use crate::ffi::force_field_helpers::ffi;
use crate::force_field::{ForceField, MMFFMolProperties, MMFFVerbosity};

fn conf_results(results: Vec<ffi::ConfOptimization>) -> Vec<(i32, f64)> {
    results
        .into_iter()
        .map(|r| (r.not_converged, r.energy))
        .collect()
}

fn values(params: TermParams) -> Option<Vec<f64>> {
    params.found.then_some(params.values)
}

/// Parameters for `uff_optimize_molecule_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UFFOptimizeMoleculeParams {
    /// The most minimizer iterations.
    pub max_iters: u32,
    /// Leaves out van der Waals contacts farther apart than this times their
    /// minimum.
    pub vdw_thresh: f64,
    /// The conformer; −1 is the default.
    pub conf_id: i32,
    /// Leaves out nonbonded terms between fragments.
    pub ignore_interfrag_interactions: bool,
}

impl Default for UFFOptimizeMoleculeParams {
    fn default() -> UFFOptimizeMoleculeParams {
        UFFOptimizeMoleculeParams {
            max_iters: 200,
            vdw_thresh: 10.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    }
}

/// Minimizes a conformer with UFF; 0 when it converged, 1 when it needs more
/// iterations. Mirrors `rdForceFieldHelpers.UFFOptimizeMolecule`.
/// `Error::Value` for a missing conformer, and above 21845 atoms, where
/// RDKit's minimizer overflows (U3).
pub fn uff_optimize_molecule(mol: &mut Mol) -> Result<i32> {
    uff_optimize_molecule_with_params(mol, &UFFOptimizeMoleculeParams::default())
}

/// Mirrors `rdForceFieldHelpers.UFFOptimizeMolecule(mol, maxIters,
/// vdwThresh, confId, ignoreInterfragInteractions)`.
pub fn uff_optimize_molecule_with_params(
    mol: &mut Mol,
    params: &UFFOptimizeMoleculeParams,
) -> Result<i32> {
    Ok(ffi::uff_optimize_molecule(
        mol.pin_mut(),
        params.max_iters,
        params.vdw_thresh,
        params.conf_id,
        params.ignore_interfrag_interactions,
    )?)
}

/// Parameters for `uff_optimize_molecule_confs_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UFFOptimizeMoleculeConfsParams {
    /// Threads to use; 0 for all, negative for that many fewer.
    pub num_threads: i32,
    /// The most minimizer iterations per conformer.
    pub max_iters: u32,
    /// Leaves out van der Waals contacts farther apart than this times their
    /// minimum.
    pub vdw_thresh: f64,
    /// Leaves out nonbonded terms between fragments.
    pub ignore_interfrag_interactions: bool,
}

impl Default for UFFOptimizeMoleculeConfsParams {
    fn default() -> UFFOptimizeMoleculeConfsParams {
        UFFOptimizeMoleculeConfsParams {
            num_threads: 1,
            max_iters: 200,
            vdw_thresh: 10.0,
            ignore_interfrag_interactions: true,
        }
    }
}

/// Minimizes every conformer with UFF, returning `(not_converged, energy)`
/// for each. Mirrors `rdForceFieldHelpers.UFFOptimizeMoleculeConfs`.
pub fn uff_optimize_molecule_confs(mol: &mut Mol) -> Result<Vec<(i32, f64)>> {
    uff_optimize_molecule_confs_with_params(mol, &UFFOptimizeMoleculeConfsParams::default())
}

/// Mirrors `rdForceFieldHelpers.UFFOptimizeMoleculeConfs(mol, numThreads,
/// maxIters, vdwThresh, ignoreInterfragInteractions)`.
pub fn uff_optimize_molecule_confs_with_params(
    mol: &mut Mol,
    params: &UFFOptimizeMoleculeConfsParams,
) -> Result<Vec<(i32, f64)>> {
    Ok(conf_results(ffi::uff_optimize_molecule_confs(
        mol.pin_mut(),
        params.num_threads,
        params.max_iters,
        params.vdw_thresh,
        params.ignore_interfrag_interactions,
    )?))
}

/// Parameters for `uff_get_molecule_force_field_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct UFFGetMoleculeForceFieldParams {
    /// Leaves out van der Waals contacts farther apart than this times their
    /// minimum.
    pub vdw_thresh: f64,
    /// The conformer; −1 is the default.
    pub conf_id: i32,
    /// Leaves out nonbonded terms between fragments.
    pub ignore_interfrag_interactions: bool,
}

impl Default for UFFGetMoleculeForceFieldParams {
    fn default() -> UFFGetMoleculeForceFieldParams {
        UFFGetMoleculeForceFieldParams {
            vdw_thresh: 10.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    }
}

/// A UFF force field for a conformer, which borrows the molecule. Mirrors
/// `rdForceFieldHelpers.UFFGetMoleculeForceField`. `Error::Value` for a
/// missing conformer and above 46340 atoms (U3); a molecule without atoms is
/// `Error::Runtime`.
pub fn uff_get_molecule_force_field(mol: &mut Mol) -> Result<ForceField<'_>> {
    uff_get_molecule_force_field_with_params(mol, &UFFGetMoleculeForceFieldParams::default())
}

/// Mirrors `rdForceFieldHelpers.UFFGetMoleculeForceField(mol, vdwThresh,
/// confId, ignoreInterfragInteractions)`.
pub fn uff_get_molecule_force_field_with_params<'m>(
    mol: &'m mut Mol,
    params: &UFFGetMoleculeForceFieldParams,
) -> Result<ForceField<'m>> {
    let ptr = ffi::uff_get_molecule_force_field(
        mol.pin_mut(),
        params.vdw_thresh,
        params.conf_id,
        params.ignore_interfrag_interactions,
    )?;
    Ok(ForceField::from_ffi(ptr, mol))
}

/// Whether UFF has parameters for every atom. Mirrors
/// `rdForceFieldHelpers.UFFHasAllMoleculeParams`.
pub fn uff_has_all_molecule_params(mol: &Mol) -> Result<bool> {
    Ok(ffi::uff_has_all_molecule_params(mol.as_ffi())?)
}

/// Parameters for `mmff_optimize_molecule_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MMFFOptimizeMoleculeParams {
    /// `"MMFF94"` or `"MMFF94s"`.
    pub mmff_variant: String,
    /// The most minimizer iterations.
    pub max_iters: u32,
    /// Leaves out nonbonded contacts farther apart than this times their
    /// minimum.
    pub non_bonded_thresh: f64,
    /// The conformer; −1 is the default.
    pub conf_id: i32,
    /// Leaves out nonbonded terms between fragments.
    pub ignore_interfrag_interactions: bool,
}

impl Default for MMFFOptimizeMoleculeParams {
    fn default() -> MMFFOptimizeMoleculeParams {
        MMFFOptimizeMoleculeParams {
            mmff_variant: "MMFF94".to_owned(),
            max_iters: 200,
            non_bonded_thresh: 100.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    }
}

/// Minimizes a conformer with MMFF; 0 when it converged, 1 when it needs
/// more iterations, −1 when MMFF cannot type the molecule. Mirrors
/// `rdForceFieldHelpers.MMFFOptimizeMolecule`. `Error::Value` above 21845
/// atoms, where RDKit's minimizer overflows (U3).
pub fn mmff_optimize_molecule(mol: &mut Mol) -> Result<i32> {
    mmff_optimize_molecule_with_params(mol, &MMFFOptimizeMoleculeParams::default())
}

/// Mirrors `rdForceFieldHelpers.MMFFOptimizeMolecule(mol, mmffVariant,
/// maxIters, nonBondedThresh, confId, ignoreInterfragInteractions)`.
pub fn mmff_optimize_molecule_with_params(
    mol: &mut Mol,
    params: &MMFFOptimizeMoleculeParams,
) -> Result<i32> {
    Ok(ffi::mmff_optimize_molecule(
        mol.pin_mut(),
        &params.mmff_variant,
        params.max_iters,
        params.non_bonded_thresh,
        params.conf_id,
        params.ignore_interfrag_interactions,
    )?)
}

/// Sanitizes a molecule as MMFF needs, returning the operation that failed
/// or `SANITIZE_NONE`. Mirrors `rdForceFieldHelpers.MMFFSanitizeMolecule`.
pub fn mmff_sanitize_molecule(mol: &mut Mol) -> Result<SanitizeFlags> {
    Ok(SanitizeFlags::from_bits_retain(
        ffi::mmff_sanitize_molecule(mol.pin_mut())?,
    ))
}

/// Parameters for `mmff_get_molecule_properties_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MMFFGetMoleculePropertiesParams {
    /// `"MMFF94"` or `"MMFF94s"`.
    pub mmff_variant: String,
    /// What the setup prints.
    pub mmff_verbosity: MMFFVerbosity,
}

impl Default for MMFFGetMoleculePropertiesParams {
    fn default() -> MMFFGetMoleculePropertiesParams {
        MMFFGetMoleculePropertiesParams {
            mmff_variant: "MMFF94".to_owned(),
            mmff_verbosity: MMFFVerbosity::MmffVerbosityNone,
        }
    }
}

/// MMFF properties for a molecule, or `None` when MMFF cannot type every
/// atom. Mirrors `rdForceFieldHelpers.MMFFGetMoleculeProperties`; like
/// `MMFFMolProperties::new`, it kekulizes the molecule.
pub fn mmff_get_molecule_properties(mol: &mut Mol) -> Result<Option<MMFFMolProperties>> {
    mmff_get_molecule_properties_with_params(mol, &MMFFGetMoleculePropertiesParams::default())
}

/// Mirrors `rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant,
/// mmffVerbosity)`.
pub fn mmff_get_molecule_properties_with_params(
    mol: &mut Mol,
    params: &MMFFGetMoleculePropertiesParams,
) -> Result<Option<MMFFMolProperties>> {
    let ptr = ffi::mmff_get_molecule_properties(
        mol.pin_mut(),
        &params.mmff_variant,
        params.mmff_verbosity.to_ffi(),
    )?;
    Ok((!ptr.is_null()).then(|| MMFFMolProperties::from_ffi(ptr)))
}

/// Parameters for `mmff_get_molecule_force_field_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MMFFGetMoleculeForceFieldParams {
    /// Leaves out nonbonded contacts farther apart than this times their
    /// minimum.
    pub non_bonded_thresh: f64,
    /// The conformer; −1 is the default.
    pub conf_id: i32,
    /// Leaves out nonbonded terms between fragments.
    pub ignore_interfrag_interactions: bool,
}

impl Default for MMFFGetMoleculeForceFieldParams {
    fn default() -> MMFFGetMoleculeForceFieldParams {
        MMFFGetMoleculeForceFieldParams {
            non_bonded_thresh: 100.0,
            conf_id: -1,
            ignore_interfrag_interactions: true,
        }
    }
}

/// An MMFF force field for a conformer, which borrows the molecule. Mirrors
/// `rdForceFieldHelpers.MMFFGetMoleculeForceField`. Properties for a
/// molecule with another number of atoms are `Error::Value`, invalid ones
/// `Error::Runtime`; Python's `None` for the properties, which returns
/// `None`, is not bound.
pub fn mmff_get_molecule_force_field<'m>(
    mol: &'m mut Mol,
    mmff_mol_properties: &MMFFMolProperties,
) -> Result<ForceField<'m>> {
    mmff_get_molecule_force_field_with_params(
        mol,
        mmff_mol_properties,
        &MMFFGetMoleculeForceFieldParams::default(),
    )
}

/// Mirrors `rdForceFieldHelpers.MMFFGetMoleculeForceField(mol,
/// MMFFMolProperties, nonBondedThresh, confId,
/// ignoreInterfragInteractions)`.
pub fn mmff_get_molecule_force_field_with_params<'m>(
    mol: &'m mut Mol,
    mmff_mol_properties: &MMFFMolProperties,
    params: &MMFFGetMoleculeForceFieldParams,
) -> Result<ForceField<'m>> {
    let ptr = ffi::mmff_get_molecule_force_field(
        mol.pin_mut(),
        mmff_mol_properties.raw(),
        params.non_bonded_thresh,
        params.conf_id,
        params.ignore_interfrag_interactions,
    )?;
    Ok(ForceField::from_ffi(ptr, mol))
}

/// A force field with a conformer's positions and no contributions, which
/// borrows the molecule. Mirrors
/// `rdForceFieldHelpers.CreateEmptyForceFieldForMol`.
pub fn create_empty_force_field_for_mol(mol: &mut Mol) -> Result<ForceField<'_>> {
    create_empty_force_field_for_mol_with_conf_id(mol, -1)
}

/// Mirrors `rdForceFieldHelpers.CreateEmptyForceFieldForMol(mol, confId)`.
pub fn create_empty_force_field_for_mol_with_conf_id(
    mol: &mut Mol,
    conf_id: i32,
) -> Result<ForceField<'_>> {
    let ptr = ffi::create_empty_force_field_for_mol(mol.pin_mut(), conf_id)?;
    Ok(ForceField::from_ffi(ptr, mol))
}

/// Whether MMFF can type every atom; the molecule is left unchanged.
/// Mirrors `rdForceFieldHelpers.MMFFHasAllMoleculeParams`.
pub fn mmff_has_all_molecule_params(mol: &Mol) -> Result<bool> {
    Ok(ffi::mmff_has_all_molecule_params(mol.as_ffi())?)
}

/// Parameters for `mmff_optimize_molecule_confs_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MMFFOptimizeMoleculeConfsParams {
    /// Threads to use; 0 for all, negative for that many fewer.
    pub num_threads: i32,
    /// The most minimizer iterations per conformer.
    pub max_iters: u32,
    /// `"MMFF94"` or `"MMFF94s"`.
    pub mmff_variant: String,
    /// Leaves out nonbonded contacts farther apart than this times their
    /// minimum.
    pub non_bonded_thresh: f64,
    /// Leaves out nonbonded terms between fragments.
    pub ignore_interfrag_interactions: bool,
}

impl Default for MMFFOptimizeMoleculeConfsParams {
    fn default() -> MMFFOptimizeMoleculeConfsParams {
        MMFFOptimizeMoleculeConfsParams {
            num_threads: 1,
            max_iters: 200,
            mmff_variant: "MMFF94".to_owned(),
            non_bonded_thresh: 100.0,
            ignore_interfrag_interactions: true,
        }
    }
}

/// Minimizes every conformer with MMFF, returning `(not_converged, energy)`
/// for each, or `(-1, -1.0)` for each when MMFF cannot type the molecule.
/// Mirrors `rdForceFieldHelpers.MMFFOptimizeMoleculeConfs`.
pub fn mmff_optimize_molecule_confs(mol: &mut Mol) -> Result<Vec<(i32, f64)>> {
    mmff_optimize_molecule_confs_with_params(mol, &MMFFOptimizeMoleculeConfsParams::default())
}

/// Mirrors `rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol, numThreads,
/// maxIters, mmffVariant, nonBondedThresh, ignoreInterfragInteractions)`.
pub fn mmff_optimize_molecule_confs_with_params(
    mol: &mut Mol,
    params: &MMFFOptimizeMoleculeConfsParams,
) -> Result<Vec<(i32, f64)>> {
    Ok(conf_results(ffi::mmff_optimize_molecule_confs(
        mol.pin_mut(),
        params.num_threads,
        params.max_iters,
        &params.mmff_variant,
        params.non_bonded_thresh,
        params.ignore_interfrag_interactions,
    )?))
}

/// UFF's `(kb, r0)` for the bond between two atoms, or `None`. Mirrors
/// `rdForceFieldHelpers.GetUFFBondStretchParams`.
pub fn get_uff_bond_stretch_params(mol: &Mol, idx1: u32, idx2: u32) -> Result<Option<(f64, f64)>> {
    let params = ffi::get_uff_bond_stretch_params(mol.as_ffi(), idx1, idx2)?;
    Ok(values(params).map(|v| (v[0], v[1])))
}

/// UFF's `(ka, theta0)` for an angle, or `None`. Mirrors
/// `rdForceFieldHelpers.GetUFFAngleBendParams`.
pub fn get_uff_angle_bend_params(
    mol: &Mol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
) -> Result<Option<(f64, f64)>> {
    let params = ffi::get_uff_angle_bend_params(mol.as_ffi(), idx1, idx2, idx3)?;
    Ok(values(params).map(|v| (v[0], v[1])))
}

/// UFF's torsion barrier `V`, or `None`. Mirrors
/// `rdForceFieldHelpers.GetUFFTorsionParams`.
pub fn get_uff_torsion_params(
    mol: &Mol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
) -> Result<Option<f64>> {
    let params = ffi::get_uff_torsion_params(mol.as_ffi(), idx1, idx2, idx3, idx4)?;
    Ok(values(params).map(|v| v[0]))
}

/// UFF's inversion force constant `K`, or `None`. Mirrors
/// `rdForceFieldHelpers.GetUFFInversionParams`.
pub fn get_uff_inversion_params(
    mol: &Mol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
) -> Result<Option<f64>> {
    let params = ffi::get_uff_inversion_params(mol.as_ffi(), idx1, idx2, idx3, idx4)?;
    Ok(values(params).map(|v| v[0]))
}

/// UFF's van der Waals `(x_ij, D_ij)` for two atoms, or `None`. Mirrors
/// `rdForceFieldHelpers.GetUFFVdWParams`.
pub fn get_uff_vdw_params(mol: &Mol, idx1: u32, idx2: u32) -> Result<Option<(f64, f64)>> {
    let params = ffi::get_uff_vdw_params(mol.as_ffi(), idx1, idx2)?;
    Ok(values(params).map(|v| (v[0], v[1])))
}

/// Initializes and minimizes a force field; 0 when it converged, 1 when it
/// needs more iterations. Mirrors `rdForceFieldHelpers.OptimizeMolecule`.
pub fn optimize_molecule(ff: &mut ForceField<'_>) -> Result<i32> {
    optimize_molecule_with_max_iters(ff, 200)
}

/// Mirrors `rdForceFieldHelpers.OptimizeMolecule(ff, maxIters)`.
pub fn optimize_molecule_with_max_iters(ff: &mut ForceField<'_>, max_iters: u32) -> Result<i32> {
    Ok(ffi::optimize_molecule(ff.raw_mut(), max_iters)?)
}

/// Parameters for `optimize_molecule_confs_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct OptimizeMoleculeConfsParams {
    /// Threads to use; 0 for all, negative for that many fewer.
    pub num_threads: i32,
    /// The most minimizer iterations per conformer.
    pub max_iters: u32,
}

impl Default for OptimizeMoleculeConfsParams {
    fn default() -> OptimizeMoleculeConfsParams {
        OptimizeMoleculeConfsParams {
            num_threads: 1,
            max_iters: 200,
        }
    }
}

/// Minimizes each of `mol`'s conformers with a force field built for
/// another molecule with the same atoms, returning `(not_converged, energy)`
/// for each; the force field's own points are restored afterwards. Mirrors
/// `rdForceFieldHelpers.OptimizeMoleculeConfs`.
///
/// `Error::Value` for a force field with fewer points than `mol` has atoms,
/// which RDKit would write past, and, with more than one thread, for one
/// with extra points, which RDKit's per-thread copies drop. As in RDKit,
/// those copies also drop the fixed points.
pub fn optimize_molecule_confs(mol: &mut Mol, ff: &mut ForceField<'_>) -> Result<Vec<(i32, f64)>> {
    optimize_molecule_confs_with_params(mol, ff, &OptimizeMoleculeConfsParams::default())
}

/// Mirrors `rdForceFieldHelpers.OptimizeMoleculeConfs(mol, ff, numThreads,
/// maxIters)`.
pub fn optimize_molecule_confs_with_params(
    mol: &mut Mol,
    ff: &mut ForceField<'_>,
    params: &OptimizeMoleculeConfsParams,
) -> Result<Vec<(i32, f64)>> {
    Ok(conf_results(ffi::optimize_molecule_confs(
        mol.pin_mut(),
        ff.raw_mut(),
        params.num_threads,
        params.max_iters,
    )?))
}
