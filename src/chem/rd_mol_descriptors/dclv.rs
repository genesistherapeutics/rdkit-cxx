// Copyright (c) 2026 Steven Kearnes
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
    /// the like) are left out, by their PDB residue info. With the default
    /// probe radius and depth, a protein is computed with a probe radius of
    /// 1.4 and a depth of 2.
    pub is_protein: bool,
    /// For a protein, keeps its HETATM ligands.
    pub include_ligand: bool,
    /// The solvent probe's radius; negative for 1.2, or 1.4 for a protein.
    pub probe_radius: f64,
    /// How finely each atom's sphere is divided: 20 × 4^`depth` dots, from 0
    /// (below is 0) to 8 (above is `Error::Runtime`).
    pub depth: i32,
    /// The dots per unit area, instead of `depth`'s; only 0, which uses
    /// `depth`, is accepted (see `DoubleCubicLatticeVolume`).
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

/// Surface areas and volumes of a molecule's default conformer by
/// Eisenhaber's double cubic lattice method, without its hydrogens.
/// Mirrors `rdMolDescriptors.DoubleCubicLatticeVolume`, which computes them
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
/// A molecule without a conformer is `Error::Conformer`. RDKit trusts its
/// inputs, so the constructors refuse what it would mishandle, with
/// `Error::Value`:
///
/// - a nonzero `dot_density`: RDKit 2024.09 never keeps the dots it makes for
///   one, and divides the volumes by a count it never set;
/// - a probe radius beyond 1e6, and coordinates beyond 1e6, or either not
///   finite: RDKit casts them, scaled, to integer indices into its voxel
///   grid;
/// - monomer info that is not PDB residue info, which RDKit casts to it
///   regardless, and PDB atom names shorter than RDKit reads (one character,
///   or two for a protein);
/// - for a protein without its ligand, no atom of a standard residue (and
///   some not solvent): RDKit's grid is then empty, and it reads the grid's
///   uninitialized scales.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleCubicLatticeVolume {
    results: ffi::DclvResults,
}

impl DoubleCubicLatticeVolume {
    /// With Python's defaults. Mirrors `DoubleCubicLatticeVolume(mol)`.
    pub fn new(mol: &Mol) -> Result<DoubleCubicLatticeVolume> {
        DoubleCubicLatticeVolume::new_with_params(mol, &DoubleCubicLatticeVolumeParams::default())
    }

    /// Mirrors `DoubleCubicLatticeVolume(mol, isProtein, includeLigand,
    /// probeRadius, depth, dotDensity)`.
    pub fn new_with_params(
        mol: &Mol,
        params: &DoubleCubicLatticeVolumeParams,
    ) -> Result<DoubleCubicLatticeVolume> {
        let results = ffi::dclv_compute(
            mol.as_ffi(),
            params.is_protein,
            params.include_ligand,
            params.probe_radius,
            params.depth,
            params.dot_density,
        )?;
        Ok(DoubleCubicLatticeVolume { results })
    }

    /// The solvent-accessible surface area. Mirrors
    /// `DoubleCubicLatticeVolume.GetSurfaceArea`.
    pub fn get_surface_area(&self) -> f64 {
        self.results.surface_area
    }

    /// The volume the probe sphere bounds. Mirrors
    /// `DoubleCubicLatticeVolume.GetVolume`.
    pub fn get_volume(&self) -> f64 {
        self.results.volume
    }

    /// The van der Waals volume. Mirrors
    /// `DoubleCubicLatticeVolume.GetVDWVolume`.
    pub fn get_vdw_volume(&self) -> f64 {
        self.results.vdw_volume
    }

    /// The compactness, for a protein. Mirrors
    /// `DoubleCubicLatticeVolume.GetCompactness`.
    pub fn get_compactness(&self) -> f64 {
        self.results.compactness
    }

    /// The packing density, for a protein. Mirrors
    /// `DoubleCubicLatticeVolume.GetPackingDensity`.
    pub fn get_packing_density(&self) -> f64 {
        self.results.packing_density
    }
}
