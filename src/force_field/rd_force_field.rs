// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.ForceField.rdForceField`: force fields and MMFF molecule
//! properties.
//!
//! A [`ForceField`] moves the atoms of the conformer it was built for, so it
//! borrows that molecule mutably for as long as it lives; read the molecule
//! again once the force field is dropped.

use std::fmt;
use std::marker::PhantomData;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::rdtrajectory::Snapshot;
use crate::error::expect_ok;
use crate::ffi::force_field::ffi;
use crate::geometry::Point3D;

rdkit_enum! {
    /// How much an MMFF setup prints. Mirrors `rdForceField.MMFFVerbosity`.
    MMFFVerbosity: ffi::MMFFVerbosity {
        MmffVerbosityNone = MMFF_VERBOSITY_NONE,
        MmffVerbosityLow = MMFF_VERBOSITY_LOW,
        MmffVerbosityHigh = MMFF_VERBOSITY_HIGH,
    }
}

/// A force field over a conformer's atoms, plus any extra points. Mirrors
/// `rdForceField.ForceField`; the `rdForceFieldHelpers` functions build one.
///
/// The force field holds pointers to the conformer's positions, which
/// [`minimize`](ForceField::minimize) moves, so it borrows the molecule
/// mutably for `'m`.
pub struct ForceField<'m> {
    /// Never null.
    ptr: UniquePtr<ffi::ForceField>,
    _mol: PhantomData<&'m mut Mol>,
}

/// Parameters for `ForceField::minimize_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MinimizeParams {
    /// The most iterations to run.
    pub max_its: u32,
    /// The convergence threshold on the gradient.
    pub force_tol: f64,
    /// The convergence threshold on the energy.
    pub energy_tol: f64,
}

impl Default for MinimizeParams {
    fn default() -> MinimizeParams {
        MinimizeParams {
            max_its: 200,
            force_tol: 1e-4,
            energy_tol: 1e-6,
        }
    }
}

/// Parameters for `ForceField::minimize_trajectory_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MinimizeTrajectoryParams {
    /// The most iterations to run.
    pub max_its: u32,
    /// The convergence threshold on the gradient.
    pub force_tol: f64,
    /// The convergence threshold on the energy.
    pub energy_tol: f64,
}

impl Default for MinimizeTrajectoryParams {
    fn default() -> MinimizeTrajectoryParams {
        MinimizeTrajectoryParams {
            max_its: 200,
            force_tol: 1e-4,
            energy_tol: 1e-6,
        }
    }
}

impl<'m> ForceField<'m> {
    /// Takes a force field whose positions point into `_mol`.
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::ForceField>, _mol: &'m mut Mol) -> ForceField<'m> {
        assert!(!ptr.is_null(), "RDKit returned a null force field");
        ForceField {
            ptr,
            _mol: PhantomData,
        }
    }

    pub(crate) fn raw(&self) -> &ffi::ForceField {
        self.ptr.as_ref().unwrap()
    }

    pub(crate) fn raw_mut(&mut self) -> Pin<&mut ffi::ForceField> {
        self.ptr.pin_mut()
    }

    /// The energy, in kcal/mol, of the current positions. Mirrors
    /// `ForceField.CalcEnergy()`.
    pub fn calc_energy(&self) -> Result<f64> {
        Ok(ffi::force_field_calc_energy(self.raw())?)
    }

    /// The energy of the coordinates `pos`, `dimension() * num_points()` of
    /// them. Mirrors `ForceField.CalcEnergy(pos)`. `Error::Value` for another
    /// length, and for points added since the last `initialize`, which
    /// RDKit would read past `pos` for.
    pub fn calc_energy_with_pos(&mut self, pos: &[f64]) -> Result<f64> {
        Ok(ffi::force_field_calc_energy_with_pos(self.raw_mut(), pos)?)
    }

    /// The gradient at the current positions, one value per coordinate.
    /// Mirrors `ForceField.CalcGrad()`. `Error::Value` for points added since
    /// the last `initialize`, which RDKit would write past the gradient for.
    pub fn calc_grad(&self) -> Result<Vec<f64>> {
        Ok(ffi::force_field_calc_grad(self.raw())?)
    }

    /// The gradient at the coordinates `pos`. Mirrors
    /// `ForceField.CalcGrad(pos)`; `Error::Value` as for
    /// `calc_energy_with_pos`.
    pub fn calc_grad_with_pos(&mut self, pos: &[f64]) -> Result<Vec<f64>> {
        Ok(ffi::force_field_calc_grad_with_pos(self.raw_mut(), pos)?)
    }

    /// The coordinates of every point, extra points included, flattened.
    /// Mirrors `ForceField.Positions`.
    pub fn positions(&self) -> Vec<f64> {
        expect_ok(ffi::force_field_positions(self.raw()))
    }

    /// Mirrors `ForceField.Dimension`.
    pub fn dimension(&self) -> u32 {
        expect_ok(ffi::force_field_dimension(self.raw()))
    }

    /// The number of points as of the last `initialize`. Mirrors
    /// `ForceField.NumPoints`.
    pub fn num_points(&self) -> u32 {
        expect_ok(ffi::force_field_num_points(self.raw()))
    }

    /// Minimizes the energy, moving the conformer's atoms; 0 when it
    /// converged, 1 when it needs more iterations. Mirrors
    /// `ForceField.Minimize`. `Error::Value` above 65535 coordinates, where
    /// RDKit's minimizer overflows (U3).
    pub fn minimize(&mut self) -> Result<i32> {
        self.minimize_with_params(&MinimizeParams::default())
    }

    /// Mirrors `ForceField.Minimize(maxIts, forceTol, energyTol)`.
    pub fn minimize_with_params(&mut self, params: &MinimizeParams) -> Result<i32> {
        Ok(ffi::force_field_minimize(
            self.raw_mut(),
            params.max_its,
            params.force_tol,
            params.energy_tol,
        )?)
    }

    /// As `minimize`, also returning a snapshot every `snapshot_freq`
    /// iterations. Mirrors `ForceField.MinimizeTrajectory`.
    ///
    /// The snapshots read their points as a trajectory of the force field's
    /// dimension (at most 3) and point count would. Python's belong to no
    /// trajectory, so their point getters always raise.
    pub fn minimize_trajectory(&mut self, snapshot_freq: u32) -> Result<(i32, Vec<Snapshot>)> {
        self.minimize_trajectory_with_params(snapshot_freq, &MinimizeTrajectoryParams::default())
    }

    /// Mirrors `ForceField.MinimizeTrajectory(snapshotFreq, maxIts,
    /// forceTol, energyTol)`.
    pub fn minimize_trajectory_with_params(
        &mut self,
        snapshot_freq: u32,
        params: &MinimizeTrajectoryParams,
    ) -> Result<(i32, Vec<Snapshot>)> {
        let mut energies = Vec::new();
        let mut positions = Vec::new();
        let res = ffi::force_field_minimize_trajectory(
            self.raw_mut(),
            snapshot_freq,
            params.max_its,
            params.force_tol,
            params.energy_tol,
            &mut energies,
            &mut positions,
        )?;
        let dimension = self.dimension().min(3);
        let size = (dimension * self.num_points()) as usize;
        let snapshots = energies
            .into_iter()
            .enumerate()
            .map(|(i, energy)| {
                Snapshot::with_layout(
                    positions[i * size..(i + 1) * size].to_vec(),
                    energy,
                    dimension,
                    self.num_points(),
                )
            })
            .collect();
        Ok((res, snapshots))
    }

    /// Adds a flat-bottomed distance constraint. Mirrors
    /// `ForceField.AddDistanceConstraint`, which RDKit deprecates for
    /// `uff_add_distance_constraint`.
    pub fn add_distance_constraint(
        &mut self,
        idx1: u32,
        idx2: u32,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) -> Result<()> {
        Ok(ffi::force_field_add_distance_constraint(
            self.raw_mut(),
            idx1,
            idx2,
            min_len,
            max_len,
            force_constant,
        )?)
    }

    /// Keeps a point from moving. Mirrors `ForceField.AddFixedPoint`, which
    /// accepts any index; here an index past the points is `Error::Index`.
    pub fn add_fixed_point(&mut self, idx: u32) -> Result<()> {
        Ok(ffi::force_field_add_fixed_point(self.raw_mut(), idx)?)
    }

    /// Adds a distance constraint; with `relative`, the bounds are offsets
    /// from the current distance. Mirrors
    /// `ForceField.UFFAddDistanceConstraint`.
    pub fn uff_add_distance_constraint(
        &mut self,
        idx1: u32,
        idx2: u32,
        relative: bool,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) -> Result<()> {
        Ok(ffi::force_field_add_relative_distance_constraint(
            self.raw_mut(),
            idx1,
            idx2,
            relative,
            min_len,
            max_len,
            force_constant,
        )?)
    }

    /// Adds an angle constraint, in degrees; with `relative`, the bounds are
    /// offsets from the current angle. Mirrors
    /// `ForceField.UFFAddAngleConstraint`.
    #[allow(clippy::too_many_arguments)]
    pub fn uff_add_angle_constraint(
        &mut self,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        relative: bool,
        min_angle_deg: f64,
        max_angle_deg: f64,
        force_constant: f64,
    ) -> Result<()> {
        Ok(ffi::force_field_add_angle_constraint(
            self.raw_mut(),
            idx1,
            idx2,
            idx3,
            relative,
            min_angle_deg,
            max_angle_deg,
            force_constant,
        )?)
    }

    /// Adds a dihedral constraint, in degrees; with `relative`, the bounds
    /// are offsets from the current dihedral. Mirrors
    /// `ForceField.UFFAddTorsionConstraint`.
    #[allow(clippy::too_many_arguments)]
    pub fn uff_add_torsion_constraint(
        &mut self,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        relative: bool,
        min_dihedral_deg: f64,
        max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Result<()> {
        Ok(ffi::force_field_add_torsion_constraint(
            self.raw_mut(),
            idx1,
            idx2,
            idx3,
            idx4,
            relative,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        )?)
    }

    /// Keeps a point within `max_displ` of where it is now. Mirrors
    /// `ForceField.UFFAddPositionConstraint`.
    pub fn uff_add_position_constraint(
        &mut self,
        idx: u32,
        max_displ: f64,
        force_constant: f64,
    ) -> Result<()> {
        Ok(ffi::force_field_add_position_constraint(
            self.raw_mut(),
            idx,
            max_displ,
            force_constant,
        )?)
    }

    /// Mirrors `ForceField.MMFFAddDistanceConstraint`; RDKit's MMFF and UFF
    /// constraints are the same.
    pub fn mmff_add_distance_constraint(
        &mut self,
        idx1: u32,
        idx2: u32,
        relative: bool,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) -> Result<()> {
        self.uff_add_distance_constraint(idx1, idx2, relative, min_len, max_len, force_constant)
    }

    /// Mirrors `ForceField.MMFFAddAngleConstraint`.
    #[allow(clippy::too_many_arguments)]
    pub fn mmff_add_angle_constraint(
        &mut self,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        relative: bool,
        min_angle_deg: f64,
        max_angle_deg: f64,
        force_constant: f64,
    ) -> Result<()> {
        self.uff_add_angle_constraint(
            idx1,
            idx2,
            idx3,
            relative,
            min_angle_deg,
            max_angle_deg,
            force_constant,
        )
    }

    /// Mirrors `ForceField.MMFFAddTorsionConstraint`.
    #[allow(clippy::too_many_arguments)]
    pub fn mmff_add_torsion_constraint(
        &mut self,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        relative: bool,
        min_dihedral_deg: f64,
        max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Result<()> {
        self.uff_add_torsion_constraint(
            idx1,
            idx2,
            idx3,
            idx4,
            relative,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        )
    }

    /// Mirrors `ForceField.MMFFAddPositionConstraint`.
    pub fn mmff_add_position_constraint(
        &mut self,
        idx: u32,
        max_displ: f64,
        force_constant: f64,
    ) -> Result<()> {
        self.uff_add_position_constraint(idx, max_displ, force_constant)
    }

    /// Resets the force field after points are added. A force field with no
    /// points is `Error::Runtime`. Mirrors `ForceField.Initialize`.
    pub fn initialize(&mut self) -> Result<()> {
        Ok(ffi::force_field_initialize(self.raw_mut())?)
    }

    /// Adds a fixed point, returning the new number of points. Mirrors
    /// `ForceField.AddExtraPoint`.
    pub fn add_extra_point(&mut self, x: f64, y: f64, z: f64) -> u32 {
        self.add_extra_point_with_fixed(x, y, z, true)
    }

    /// Mirrors `ForceField.AddExtraPoint(x, y, z, fixed)`.
    pub fn add_extra_point_with_fixed(&mut self, x: f64, y: f64, z: f64, fixed: bool) -> u32 {
        expect_ok(ffi::force_field_add_extra_point(
            self.raw_mut(),
            x,
            y,
            z,
            fixed,
        ))
    }

    /// The position of the extra point `idx`, counting from 0. Mirrors
    /// `ForceField.GetExtraPointPos`.
    pub fn get_extra_point_pos(&self, idx: u32) -> Result<Point3D> {
        Ok(ffi::force_field_get_extra_point_pos(self.raw(), idx)?)
    }
}

impl fmt::Debug for ForceField<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ForceField")
            .field("dimension", &self.dimension())
            .field("num_points", &self.num_points())
            .finish()
    }
}

/// MMFF torsion parameters: Python's `(torsionType, V1, V2, V3)`, from
/// `MMFFMolProperties::get_mmff_torsion_params`.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct MMFFTorsionParams {
    /// The torsion type.
    pub torsion_type: u32,
    /// The onefold term, `V1`.
    pub v1: f64,
    /// The twofold term, `V2`.
    pub v2: f64,
    /// The threefold term, `V3`.
    pub v3: f64,
}

/// MMFF van der Waals parameters: Python's `(R_ij_starUnscaled,
/// epsilonUnscaled, R_ij_star, epsilon)`, from
/// `MMFFMolProperties::get_mmff_vdw_params`.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct MMFFVdWParams {
    /// `R_ij_starUnscaled`.
    pub r_ij_star_unscaled: f64,
    /// `epsilonUnscaled`.
    pub epsilon_unscaled: f64,
    /// `R_ij_star`.
    pub r_ij_star: f64,
    /// `epsilon`.
    pub epsilon: f64,
}

/// MMFF atom types, charges, and settings for a molecule. Mirrors
/// `rdForceField.MMFFMolProperties`.
pub struct MMFFMolProperties {
    /// Never null.
    ptr: UniquePtr<ffi::MMFFMolProperties>,
}

/// Parameters for `MMFFMolProperties::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct MMFFMolPropertiesParams {
    /// `"MMFF94"` or `"MMFF94s"`; RDKit takes any other name as MMFF94.
    pub mmff_variant: String,
    /// What the setup prints.
    pub verbosity: MMFFVerbosity,
}

impl Default for MMFFMolPropertiesParams {
    fn default() -> MMFFMolPropertiesParams {
        MMFFMolPropertiesParams {
            mmff_variant: "MMFF94".to_owned(),
            verbosity: MMFFVerbosity::MmffVerbosityNone,
        }
    }
}

fn term_params(params: ffi::TermParams) -> Option<(u32, Vec<f64>)> {
    params.found.then_some((params.term_type, params.values))
}

impl MMFFMolProperties {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::MMFFMolProperties>) -> MMFFMolProperties {
        assert!(!ptr.is_null(), "RDKit returned null MMFF properties");
        MMFFMolProperties { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::MMFFMolProperties {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::MMFFMolProperties> {
        self.ptr.pin_mut()
    }

    /// Types the molecule's atoms for MMFF94, kekulizing it and setting
    /// MMFF's aromaticity on it. Mirrors `MMFFMolProperties(mol)`; atoms MMFF
    /// cannot type leave the properties invalid rather than failing.
    pub fn new(mol: &mut Mol) -> Result<MMFFMolProperties> {
        MMFFMolProperties::new_with_params(mol, &MMFFMolPropertiesParams::default())
    }

    /// Mirrors `MMFFMolProperties(mol, mmffVariant, verbosity)`.
    pub fn new_with_params(
        mol: &mut Mol,
        params: &MMFFMolPropertiesParams,
    ) -> Result<MMFFMolProperties> {
        Ok(MMFFMolProperties::from_ffi(ffi::mmff_mol_properties_new(
            mol.pin_mut(),
            &params.mmff_variant,
            params.verbosity.to_ffi(),
        )?))
    }

    /// Whether MMFF typed every atom.
    pub(crate) fn is_valid(&self) -> bool {
        expect_ok(ffi::mmff_mol_properties_is_valid(self.raw()))
    }

    /// Mirrors `MMFFMolProperties.GetMMFFAtomType`.
    pub fn get_mmff_atom_type(&self, idx: u32) -> Result<u32> {
        Ok(ffi::mmff_mol_properties_get_atom_type(self.raw(), idx)?)
    }

    /// Mirrors `MMFFMolProperties.GetMMFFFormalCharge`.
    pub fn get_mmff_formal_charge(&self, idx: u32) -> Result<f64> {
        Ok(ffi::mmff_mol_properties_get_formal_charge(self.raw(), idx)?)
    }

    /// Mirrors `MMFFMolProperties.GetMMFFPartialCharge`.
    pub fn get_mmff_partial_charge(&self, idx: u32) -> Result<f64> {
        Ok(ffi::mmff_mol_properties_get_partial_charge(
            self.raw(),
            idx,
        )?)
    }

    /// `(bondType, kb, r0)` for the bond between two atoms of `mol`, or
    /// `None`. Mirrors `MMFFMolProperties.GetMMFFBondStretchParams`; an
    /// index past `mol` or these properties is `Error::Index`.
    pub fn get_mmff_bond_stretch_params(
        &self,
        mol: &Mol,
        idx1: u32,
        idx2: u32,
    ) -> Result<Option<(u32, f64, f64)>> {
        let params =
            ffi::mmff_mol_properties_get_bond_stretch_params(self.raw(), mol.as_ffi(), idx1, idx2)?;
        Ok(term_params(params).map(|(t, v)| (t, v[0], v[1])))
    }

    /// `(angleType, ka, theta0)`, or `None`. Mirrors
    /// `MMFFMolProperties.GetMMFFAngleBendParams`.
    pub fn get_mmff_angle_bend_params(
        &self,
        mol: &Mol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
    ) -> Result<Option<(u32, f64, f64)>> {
        let params = ffi::mmff_mol_properties_get_angle_bend_params(
            self.raw(),
            mol.as_ffi(),
            idx1,
            idx2,
            idx3,
        )?;
        Ok(term_params(params).map(|(t, v)| (t, v[0], v[1])))
    }

    /// `(stretchBendType, kbaIJK, kbaKJI)`, or `None`. Mirrors
    /// `MMFFMolProperties.GetMMFFStretchBendParams`.
    pub fn get_mmff_stretch_bend_params(
        &self,
        mol: &Mol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
    ) -> Result<Option<(u32, f64, f64)>> {
        let params = ffi::mmff_mol_properties_get_stretch_bend_params(
            self.raw(),
            mol.as_ffi(),
            idx1,
            idx2,
            idx3,
        )?;
        Ok(term_params(params).map(|(t, v)| (t, v[0], v[1])))
    }

    /// The torsion's parameters, or `None`. Mirrors
    /// `MMFFMolProperties.GetMMFFTorsionParams`.
    pub fn get_mmff_torsion_params(
        &self,
        mol: &Mol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
    ) -> Result<Option<MMFFTorsionParams>> {
        let params = ffi::mmff_mol_properties_get_torsion_params(
            self.raw(),
            mol.as_ffi(),
            idx1,
            idx2,
            idx3,
            idx4,
        )?;
        Ok(term_params(params).map(|(t, v)| MMFFTorsionParams {
            torsion_type: t,
            v1: v[0],
            v2: v[1],
            v3: v[2],
        }))
    }

    /// The out-of-plane bending force constant `koop`, or `None`. Mirrors
    /// `MMFFMolProperties.GetMMFFOopBendParams`.
    pub fn get_mmff_oop_bend_params(
        &self,
        mol: &Mol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
    ) -> Result<Option<f64>> {
        let params = ffi::mmff_mol_properties_get_oop_bend_params(
            self.raw(),
            mol.as_ffi(),
            idx1,
            idx2,
            idx3,
            idx4,
        )?;
        Ok(term_params(params).map(|(_, v)| v[0]))
    }

    /// The pair's van der Waals parameters, or `None`. Mirrors
    /// `MMFFMolProperties.GetMMFFVdWParams`.
    pub fn get_mmff_vdw_params(&self, idx1: u32, idx2: u32) -> Result<Option<MMFFVdWParams>> {
        let params = ffi::mmff_mol_properties_get_vdw_params(self.raw(), idx1, idx2)?;
        Ok(term_params(params).map(|(_, v)| MMFFVdWParams {
            r_ij_star_unscaled: v[0],
            epsilon_unscaled: v[1],
            r_ij_star: v[2],
            epsilon: v[3],
        }))
    }

    /// Sets the constant dielectric model. Mirrors
    /// `MMFFMolProperties.SetMMFFDielectricModel()`.
    pub fn set_mmff_dielectric_model(&mut self) {
        self.set_mmff_dielectric_model_with_diel_model(1);
    }

    /// 1 for a constant dielectric, 2 for a distance-dependent one. Mirrors
    /// `MMFFMolProperties.SetMMFFDielectricModel(dielModel)`.
    pub fn set_mmff_dielectric_model_with_diel_model(&mut self, diel_model: u8) {
        expect_ok(ffi::mmff_mol_properties_set_dielectric_model(
            self.raw_mut(),
            diel_model,
        ))
    }

    /// Sets the dielectric constant to 1. Mirrors
    /// `MMFFMolProperties.SetMMFFDielectricConstant()`.
    pub fn set_mmff_dielectric_constant(&mut self) {
        expect_ok(ffi::mmff_mol_properties_set_dielectric_constant(
            self.raw_mut(),
            1.0,
        ))
    }

    /// A constant that is not positive is `Error::Runtime`. Mirrors
    /// `MMFFMolProperties.SetMMFFDielectricConstant(dielConst)`.
    pub fn set_mmff_dielectric_constant_with_diel_const(&mut self, diel_const: f64) -> Result<()> {
        Ok(ffi::mmff_mol_properties_set_dielectric_constant(
            self.raw_mut(),
            diel_const,
        )?)
    }

    fn set_term(&mut self, term: ffi::MMFFTerm, state: bool) {
        expect_ok(ffi::mmff_mol_properties_set_term(
            self.raw_mut(),
            term,
            state,
        ))
    }

    /// Sets the MMFF variant. Mirrors `MMFFMolProperties.SetMMFFVariant()`,
    /// which sets MMFF94.
    pub fn set_mmff_variant(&mut self) {
        expect_ok(ffi::mmff_mol_properties_set_variant(
            self.raw_mut(),
            "MMFF94",
        ))
    }

    /// `"MMFF94"` or `"MMFF94s"`; another name is `Error::Runtime`. Mirrors
    /// `MMFFMolProperties.SetMMFFVariant(mmffVariant)`.
    pub fn set_mmff_variant_with_mmff_variant(&mut self, mmff_variant: &str) -> Result<()> {
        Ok(ffi::mmff_mol_properties_set_variant(
            self.raw_mut(),
            mmff_variant,
        )?)
    }

    /// Turns printing off. Mirrors `MMFFMolProperties.SetMMFFVerbosity()`.
    pub fn set_mmff_verbosity(&mut self) {
        self.set_mmff_verbosity_with_verbosity(MMFFVerbosity::MmffVerbosityNone);
    }

    /// Mirrors `MMFFMolProperties.SetMMFFVerbosity(verbosity)`.
    pub fn set_mmff_verbosity_with_verbosity(&mut self, verbosity: MMFFVerbosity) {
        expect_ok(ffi::mmff_mol_properties_set_verbosity(
            self.raw_mut(),
            verbosity.to_ffi(),
        ))
    }
}

/// Defines the setters for one MMFF term (A2).
macro_rules! mmff_terms {
    ($($term:ident: $set:ident, $set_with:ident, $python:literal, $what:literal;)*) => {
        impl MMFFMolProperties {
            $(
                #[doc = concat!("Includes the ", $what, " term. Mirrors `MMFFMolProperties.SetMMFF", $python, "Term()`.")]
                pub fn $set(&mut self) {
                    self.set_term(ffi::MMFFTerm::$term, true);
                }

                #[doc = concat!("Mirrors `MMFFMolProperties.SetMMFF", $python, "Term(state)`.")]
                pub fn $set_with(&mut self, state: bool) {
                    self.set_term(ffi::MMFFTerm::$term, state);
                }
            )*
        }
    };
}

mmff_terms! {
    Bond: set_mmff_bond_term, set_mmff_bond_term_with_state, "Bond", "bond stretching";
    Angle: set_mmff_angle_term, set_mmff_angle_term_with_state, "Angle", "angle bending";
    StretchBend: set_mmff_stretch_bend_term, set_mmff_stretch_bend_term_with_state, "StretchBend", "stretch-bend";
    Oop: set_mmff_oop_term, set_mmff_oop_term_with_state, "Oop", "out-of-plane bending";
    Torsion: set_mmff_torsion_term, set_mmff_torsion_term_with_state, "Torsion", "torsion";
    VdW: set_mmff_vdw_term, set_mmff_vdw_term_with_state, "VdW", "van der Waals";
    Ele: set_mmff_ele_term, set_mmff_ele_term_with_state, "Ele", "electrostatic";
}

impl fmt::Debug for MMFFMolProperties {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MMFFMolProperties")
            .field("valid", &self.is_valid())
            .finish()
    }
}
