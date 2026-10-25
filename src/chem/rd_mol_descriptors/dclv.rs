// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's double cubic lattice method for surface areas and
//! volumes.

use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::dclv::ffi;
use crate::geometry::Point3D;

/// Parameters for `DoubleCubicLatticeVolume::new_with_params` and
/// `DoubleCubicLatticeVolume::from_radii_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleCubicLatticeVolumeParams {
    /// Treats the molecule as a protein: solvent residues (HOH, SO4, and
    /// the like) are left out, by their PDB residue info.
    pub is_protein: bool,
    /// For a protein, keeps its HETATM ligands.
    pub include_ligand: bool,
    /// The solvent probe's radius.
    pub probe_radius: f64,
    /// The conformer; −1 is the default conformer.
    pub conf_id: i32,
}

impl Default for DoubleCubicLatticeVolumeParams {
    fn default() -> DoubleCubicLatticeVolumeParams {
        DoubleCubicLatticeVolumeParams {
            is_protein: false,
            include_ligand: true,
            probe_radius: 1.4,
            conf_id: -1,
        }
    }
}

/// Parameters for `DoubleCubicLatticeVolume::get_polar_surface_area_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetPolarSurfaceAreaParams {
    /// Counts sulfur and phosphorus as polar, with nitrogen and oxygen.
    pub include_s_and_p: bool,
    /// Counts the hydrogens on polar atoms.
    pub include_hs: bool,
}

/// Parameters for `DoubleCubicLatticeVolume::get_polar_volume_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetPolarVolumeParams {
    /// Counts sulfur and phosphorus as polar, with nitrogen and oxygen.
    pub include_s_and_p: bool,
    /// Counts the hydrogens on polar atoms.
    pub include_hs: bool,
}

/// Surface areas and volumes of a conformer by Eisenhaber's double cubic
/// lattice method. Mirrors `rdMolDescriptors.DoubleCubicLatticeVolume`.
///
/// RDKit keeps a reference to the molecule, so this borrows it. Python's
/// object does not keep its molecule alive, and dangles once the molecule
/// is gone. The results are cached as they are computed, so the methods
/// take `&mut self`.
///
/// ```no_run
/// use rdkit_cxx::chem::rd_mol_descriptors::DoubleCubicLatticeVolume;
/// # fn embedded() -> rdkit_cxx::chem::Mol { unimplemented!() }
/// let mol = embedded();
/// let mut dclv = DoubleCubicLatticeVolume::new(&mol)?;
/// let area = dclv.get_surface_area();
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
///
/// The molecule cannot change while it is borrowed:
///
/// ```compile_fail
/// use rdkit_cxx::chem::rd_mol_descriptors::DoubleCubicLatticeVolume;
/// let mut mol = rdkit_cxx::chem::mol_from_smiles("CCO")?;
/// let mut dclv = DoubleCubicLatticeVolume::new(&mol)?;
/// mol.update_property_cache()?;
/// dclv.get_surface_area();
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
///
/// RDKit trusts its inputs, so the constructors refuse what it would mishandle,
/// with `Error::Value`:
///
/// - a molecule without atoms, whose largest radius RDKit reads anyway;
/// - radii other than one per atom, and radii or a probe radius beyond 1e6
///   or not finite, and coordinates beyond 1e8 or not finite: RDKit casts
///   them, scaled, to integer indices into its voxel grid;
/// - with `is_protein`, monomer info that is not PDB residue info, which
///   RDKit casts to it regardless;
/// - an atom with radius 0, which masks it, or an atom left out as solvent,
///   before an atom with a nonzero radius: RDKit lists neighbors only for
///   the atoms it keeps but looks them up by atom index, so the later atoms
///   read another atom's list, and the last past the end (Python crashes).
pub struct DoubleCubicLatticeVolume<'m> {
    /// Never null.
    ptr: UniquePtr<ffi::DoubleCubicLatticeVolume>,
    mol: PhantomData<&'m Mol>,
}

impl<'m> DoubleCubicLatticeVolume<'m> {
    fn make(mol: &'m Mol, radii: &[f64], params: &DoubleCubicLatticeVolumeParams) -> Result<Self> {
        let ptr = ffi::dclv_new(
            mol.as_ffi(),
            radii,
            params.is_protein,
            params.include_ligand,
            params.probe_radius,
            params.conf_id,
        )?;
        Ok(DoubleCubicLatticeVolume {
            ptr,
            mol: PhantomData,
        })
    }

    fn raw(&mut self) -> Pin<&mut ffi::DoubleCubicLatticeVolume> {
        self.ptr.pin_mut()
    }

    /// With van der Waals radii, for the default conformer. Mirrors
    /// `DoubleCubicLatticeVolume(mol)`.
    pub fn new(mol: &'m Mol) -> Result<Self> {
        DoubleCubicLatticeVolume::new_with_params(mol, &DoubleCubicLatticeVolumeParams::default())
    }

    /// Mirrors `DoubleCubicLatticeVolume(mol, isProtein, includeLigand,
    /// probeRadius, confId)`.
    pub fn new_with_params(mol: &'m Mol, params: &DoubleCubicLatticeVolumeParams) -> Result<Self> {
        DoubleCubicLatticeVolume::make(mol, &[], params)
    }

    /// With one radius per atom; empty means van der Waals radii, as
    /// Python's empty list does. Mirrors `DoubleCubicLatticeVolume(mol,
    /// radii)`.
    pub fn from_radii(mol: &'m Mol, radii: &[f64]) -> Result<Self> {
        DoubleCubicLatticeVolume::from_radii_with_params(
            mol,
            radii,
            &DoubleCubicLatticeVolumeParams::default(),
        )
    }

    /// Mirrors `DoubleCubicLatticeVolume(mol, radii, isProtein,
    /// includeLigand, probeRadius, confId)`.
    pub fn from_radii_with_params(
        mol: &'m Mol,
        radii: &[f64],
        params: &DoubleCubicLatticeVolumeParams,
    ) -> Result<Self> {
        DoubleCubicLatticeVolume::make(mol, radii, params)
    }

    /// The solvent-accessible surface area. Mirrors `GetSurfaceArea`.
    pub fn get_surface_area(&mut self) -> f64 {
        expect_ok(ffi::dclv_get_surface_area(self.raw()))
    }

    /// One atom's solvent-accessible surface area; `Error::Index` for an
    /// atom out of range. Mirrors `GetAtomSurfaceArea`.
    pub fn get_atom_surface_area(&mut self, atom_idx: u32) -> Result<f64> {
        Ok(ffi::dclv_get_atom_surface_area(self.raw(), atom_idx)?)
    }

    /// The surface area of the nitrogens and oxygens. Mirrors
    /// `GetPolarSurfaceArea`.
    pub fn get_polar_surface_area(&mut self) -> f64 {
        self.get_polar_surface_area_with_params(&GetPolarSurfaceAreaParams::default())
    }

    /// Mirrors `GetPolarSurfaceArea(includeSandP, includeHs)`.
    pub fn get_polar_surface_area_with_params(
        &mut self,
        params: &GetPolarSurfaceAreaParams,
    ) -> f64 {
        expect_ok(ffi::dclv_get_polar_surface_area(
            self.raw(),
            params.include_s_and_p,
            params.include_hs,
        ))
    }

    /// The surface area of the given atoms; `Error::Index` for an atom out
    /// of range, which Python skips. Mirrors `GetPartialSurfaceArea`.
    pub fn get_partial_surface_area(&mut self, atom_indices: &[u32]) -> Result<f64> {
        Ok(ffi::dclv_get_partial_surface_area(
            self.raw(),
            atom_indices,
        )?)
    }

    /// The points of the solvent-accessible surface, by atom. Mirrors
    /// `GetSurfacePoints`.
    pub fn get_surface_points(&mut self) -> BTreeMap<u32, Vec<Point3D>> {
        self.get_surface_points_with_all_points(false)
    }

    /// With `all_points`, every lattice point of every atom's probe sphere,
    /// buried or not. Mirrors `GetSurfacePoints(allPoints)`.
    pub fn get_surface_points_with_all_points(
        &mut self,
        all_points: bool,
    ) -> BTreeMap<u32, Vec<Point3D>> {
        let (mut atoms, mut counts, mut points) = (Vec::new(), Vec::new(), Vec::new());
        expect_ok(ffi::dclv_get_surface_points(
            self.raw(),
            all_points,
            &mut atoms,
            &mut counts,
            &mut points,
        ));
        let mut points = points.into_iter();
        atoms
            .into_iter()
            .zip(counts)
            .map(|(atom, count)| (atom, points.by_ref().take(count as usize).collect()))
            .collect()
    }

    /// The volume inside the probe's surface. Mirrors `GetVolume`.
    pub fn get_volume(&mut self) -> f64 {
        expect_ok(ffi::dclv_get_volume(self.raw()))
    }

    /// The van der Waals volume. Mirrors `GetVDWVolume`.
    pub fn get_vdw_volume(&mut self) -> f64 {
        expect_ok(ffi::dclv_get_vdw_volume(self.raw()))
    }

    /// One atom's volume with the given solvent radius; `Error::Index` for
    /// an atom out of range. Mirrors `GetAtomVolume`.
    pub fn get_atom_volume(&mut self, atom_idx: u32, solvent_radius: f64) -> Result<f64> {
        Ok(ffi::dclv_get_atom_volume(
            self.raw(),
            atom_idx,
            solvent_radius,
        )?)
    }

    /// The volume of the nitrogens and oxygens. Mirrors `GetPolarVolume`.
    pub fn get_polar_volume(&mut self) -> f64 {
        self.get_polar_volume_with_params(&GetPolarVolumeParams::default())
    }

    /// Mirrors `GetPolarVolume(includeSandP, includeHs)`.
    pub fn get_polar_volume_with_params(&mut self, params: &GetPolarVolumeParams) -> f64 {
        expect_ok(ffi::dclv_get_polar_volume(
            self.raw(),
            params.include_s_and_p,
            params.include_hs,
        ))
    }

    /// The volume of the given atoms; `Error::Index` for an atom out of
    /// range, which Python skips. Mirrors `GetPartialVolume`.
    pub fn get_partial_volume(&mut self, atom_indices: &[u32]) -> Result<f64> {
        Ok(ffi::dclv_get_partial_volume(self.raw(), atom_indices)?)
    }

    /// The surface area over that of a sphere of the same volume. Mirrors
    /// `GetCompactness`.
    pub fn get_compactness(&mut self) -> f64 {
        expect_ok(ffi::dclv_get_compactness(self.raw()))
    }

    /// The van der Waals volume over the total volume. Mirrors
    /// `GetPackingDensity`.
    pub fn get_packing_density(&mut self) -> f64 {
        expect_ok(ffi::dclv_get_packing_density(self.raw()))
    }
}
