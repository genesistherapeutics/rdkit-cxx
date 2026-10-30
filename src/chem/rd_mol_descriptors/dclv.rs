// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's double cubic lattice method for surface areas and
//! volumes.

use crate::Result;
use crate::chem::Mol;
use crate::ffi::dclv::ffi;

/// Parameters for `DoubleCubicLatticeVolume::new_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleCubicLatticeVolumeParams {
    /// Treats the molecule as a protein: solvent residues (HOH, SO4, and
    /// the like) are left out, by their PDB residue info, and only backbone
    /// atoms are typed. With the default probe radius and depth, RDKit then
    /// uses a probe radius of 1.4 and a depth of 2.
    pub is_protein: bool,
    /// For a protein, keeps its HETATM ligands.
    pub include_ligand: bool,
    /// The solvent probe's radius; a negative one is RDKit's default (1.2,
    /// or 1.4 for a protein).
    pub probe_radius: f64,
    /// How finely each atom's sphere is divided: 20 × 4^depth dots, up to a
    /// depth of 8 (more is `Error::Runtime`); a negative depth is 0.
    pub depth: i32,
    /// Dots per unit of area instead of the sphere's divisions; only 0 is
    /// accepted (see `DoubleCubicLatticeVolume`).
    pub dot_density: i32,
}

impl Default for DoubleCubicLatticeVolumeParams {
    fn default() -> DoubleCubicLatticeVolumeParams {
        DoubleCubicLatticeVolumeParams {
            is_protein: false,
            include_ligand: true,
            probe_radius: 1.2,
            depth: 4,
            dot_density: 0,
        }
    }
}

/// Surface areas and volumes of a molecule's default conformer, without its
/// hydrogens, by Eisenhaber's double cubic lattice method. Mirrors
/// `rdMolDescriptors.DoubleCubicLatticeVolume`, which computes them all
/// when it is made.
///
/// ```no_run
/// use rdkit_cxx::chem::rd_mol_descriptors::DoubleCubicLatticeVolume;
/// # fn embedded() -> rdkit_cxx::chem::Mol { unimplemented!() }
/// let mol = embedded();
/// let dclv = DoubleCubicLatticeVolume::new(&mol)?;
/// let area = dclv.get_surface_area();
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
///
/// RDKit 2025.03 trusts its inputs, so the constructors refuse what it would
/// mishandle, with `Error::Value`:
///
/// - a dot density other than 0: RDKit then divides by a count it never
///   sets;
/// - a probe radius beyond 1e6 or not finite, and coordinates beyond 1e6 or
///   not finite: RDKit casts them, scaled, to integer indices into its
///   voxel grid;
/// - monomer info that is not PDB residue info, which RDKit casts to it
///   regardless, or whose atom name is empty (or, with `is_protein`, one
///   character), which RDKit reads past.
///
/// A molecule without a conformer is `Error::Value`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleCubicLatticeVolume {
    surface_area: f64,
    volume: f64,
    vdw_volume: f64,
    compactness: f64,
    packing_density: f64,
}

impl DoubleCubicLatticeVolume {
    /// For a molecule that is not a protein. Mirrors
    /// `DoubleCubicLatticeVolume(mol)`.
    pub fn new(mol: &Mol) -> Result<Self> {
        DoubleCubicLatticeVolume::new_with_params(mol, &DoubleCubicLatticeVolumeParams::default())
    }

    /// Mirrors `DoubleCubicLatticeVolume(mol, isProtein, includeLigand,
    /// probeRadius, depth, dotDensity)`.
    pub fn new_with_params(mol: &Mol, params: &DoubleCubicLatticeVolumeParams) -> Result<Self> {
        let r = ffi::dclv_compute(
            mol.as_ffi(),
            params.is_protein,
            params.include_ligand,
            params.probe_radius,
            params.depth,
            params.dot_density,
        )?;
        Ok(DoubleCubicLatticeVolume {
            surface_area: r.surface_area,
            volume: r.volume,
            vdw_volume: r.vdw_volume,
            compactness: r.compactness,
            packing_density: r.packing_density,
        })
    }

    /// The solvent-accessible surface area. Mirrors `GetSurfaceArea`.
    pub fn get_surface_area(&self) -> f64 {
        self.surface_area
    }

    /// The volume inside the probe's surface. Mirrors `GetVolume`.
    pub fn get_volume(&self) -> f64 {
        self.volume
    }

    /// The van der Waals volume. Mirrors `GetVDWVolume`.
    pub fn get_vdw_volume(&self) -> f64 {
        self.vdw_volume
    }

    /// The surface area over that of a sphere of the same volume. Mirrors
    /// `GetCompactness`.
    pub fn get_compactness(&self) -> f64 {
        self.compactness
    }

    /// The van der Waals volume over the total volume. Mirrors
    /// `GetPackingDensity`.
    pub fn get_packing_density(&self) -> f64 {
        self.packing_density
    }
}
