// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's double cubic lattice method for surface areas and
//! volumes.

use crate::Result;
use crate::chem::Mol;
use crate::ffi::dclv::ffi;

/// Parameters for `DoubleCubicLatticeVolume::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleCubicLatticeVolumeParams {
    /// Treats the molecule as a protein: solvent residues (HOH, SO4, and
    /// the like) are left out, by their PDB residue info, and the surface
    /// area is the residues'.
    pub is_protein: bool,
    /// For a protein, keeps its HETATM ligands.
    pub include_ligand: bool,
    /// The solvent probe's radius; a negative one is RDKit's default (1.4
    /// for a protein, else 1.2). With `is_protein`, the defaults of this and
    /// `depth` together become 1.4 and 2.
    pub probe_radius: f64,
    /// The tessellation depth, which sets the number of dots on each atom's
    /// sphere: 20 × 4^depth. At most 8; a negative depth is 0.
    pub depth: i32,
    /// Dots per unit area instead of a tessellation. RDKit 2024.03 reads an
    /// uninitialized count for any value but 0, so this refuses others with
    /// `Error::Value`.
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

/// Surface area and volumes of a molecule's first conformer, without its
/// hydrogens, by Eisenhaber's double cubic lattice method. Mirrors
/// `rdMolDescriptors.DoubleCubicLatticeVolume`, which computes them all when
/// it is made.
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
/// RDKit trusts its inputs, so the constructors refuse what it would mishandle,
/// with `Error::Value`:
///
/// - coordinates beyond 1e6 or not finite, and a probe radius beyond 1e5 or
///   not finite: RDKit converts them, scaled, to integer indices into its
///   voxel grid;
/// - monomer info that is not PDB residue info, which RDKit casts to it
///   regardless, and PDB atom names shorter than RDKit reads (two
///   characters for a protein, else one);
/// - with `is_protein` and not `include_ligand`, a molecule without an atom
///   of the 20 standard residues, for which RDKit reads its grid's
///   uninitialized scale;
/// - a nonzero `dot_density`.
///
/// A molecule without a conformer is `Error::Value`. A depth over 8 is
/// `Error::Runtime`, as Python's raises.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleCubicLatticeVolume {
    data: Results,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Results {
    surface_area: f64,
    volume: f64,
    vdw_volume: f64,
    compactness: f64,
    packing_density: f64,
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
        let data = ffi::dclv_compute(
            mol.as_ffi(),
            params.is_protein,
            params.include_ligand,
            params.probe_radius,
            params.depth,
            params.dot_density,
        )?;
        Ok(DoubleCubicLatticeVolume {
            data: Results {
                surface_area: data.surface_area,
                volume: data.volume,
                vdw_volume: data.vdw_volume,
                compactness: data.compactness,
                packing_density: data.packing_density,
            },
        })
    }

    /// The solvent-accessible surface area. Mirrors `GetSurfaceArea`.
    pub fn get_surface_area(&self) -> f64 {
        self.data.surface_area
    }

    /// The volume inside the probe's surface. Mirrors `GetVolume`.
    pub fn get_volume(&self) -> f64 {
        self.data.volume
    }

    /// The van der Waals volume. Mirrors `GetVDWVolume`.
    pub fn get_vdw_volume(&self) -> f64 {
        self.data.vdw_volume
    }

    /// The surface area over that of a sphere of the same volume. Mirrors
    /// `GetCompactness`.
    pub fn get_compactness(&self) -> f64 {
        self.data.compactness
    }

    /// The van der Waals volume over the total volume. Mirrors
    /// `GetPackingDensity`.
    pub fn get_packing_density(&self) -> f64 {
        self.data.packing_density
    }
}
