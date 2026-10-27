// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdDetermineBonds`: bonds, and their orders, from a
//! molecule's 3D coordinates (xyz2mol).
//!
//! ```no_run
//! use rdkit_cxx::chem::rd_determine_bonds::determine_bonds;
//! use rdkit_cxx::chem::rd_dist_geom::embed_molecule;
//! use rdkit_cxx::chem::{add_hs, mol_from_smiles};
//!
//! let mut mol = add_hs(&mol_from_smiles("CC=O").unwrap()).unwrap();
//! embed_molecule(&mut mol).unwrap();
//! // The bonds, and their orders, come back from the coordinates alone.
//! determine_bonds(&mut mol).unwrap();
//! ```
//!
//! The functions refuse, with `Error::Value`, what RDKit would mishandle:
//! for connect-the-dots, coordinates that are not finite or beyond 1e7, and
//! for the bond orders, charges near `i32::MIN` or `i32::MAX`. They hide
//! from connect-the-dots any monomer info that is not PDB residue info,
//! which RDKit would read as one.
//!
//! The Hückel option runs YAeHMOP, which needs `RDBASE` set to RDKit's
//! `share/RDKit` directory, or `BIND_PARM_FILE` to an EHT parameter file:
//! without them it reads the parameter file it was built with, and frees
//! memory twice when that is missing, as conda-forge's is, so the functions
//! refuse it. YAeHMOP exits the process for input it cannot handle, and the
//! functions check that input first: no dummy atoms, parameters for every
//! element, at least one electron and at most two per orbital, at most 46340
//! orbitals, and a parameter file it reads safely.
//!
//! Finding bond orders swaps the process's SIGINT handler for RDKit's, so
//! one such call runs at a time in the process. RDKit's handler stops the
//! search at a Ctrl-C, which then leaves the bond orders unassigned without
//! an error.

use crate::Result;
use crate::chem::Mol;
use crate::error::expect_ok;
use crate::ffi::determine_bonds::ffi;

/// Parameters for `determine_connectivity_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DetermineConnectivityParams {
    /// Uses extended Hückel theory (YAeHMOP) for the connectivity.
    pub use_hueckel: bool,
    /// The molecule's charge, which the Hückel method needs.
    pub charge: i32,
    /// The factor that scales covalent radii for the van der Waals method.
    pub cov_factor: f64,
    /// Uses the van der Waals method instead of connect-the-dots.
    pub use_vdw: bool,
}

impl Default for DetermineConnectivityParams {
    fn default() -> DetermineConnectivityParams {
        DetermineConnectivityParams {
            use_hueckel: false,
            charge: 0,
            cov_factor: 1.3,
            use_vdw: false,
        }
    }
}

/// Replaces the molecule's bonds with single bonds between atoms close
/// enough in its first conformer. Mirrors
/// `rdDetermineBonds.DetermineConnectivity`.
pub fn determine_connectivity(mol: &mut Mol) -> Result<()> {
    determine_connectivity_with_params(mol, &DetermineConnectivityParams::default())
}

/// Mirrors `rdDetermineBonds.DetermineConnectivity(mol, useHueckel, charge,
/// covFactor, useVdw)`. The Hückel method sets the first atom's formal
/// charge to `charge`; without YAeHMOP in RDKit it is `Error::Value`.
pub fn determine_connectivity_with_params(
    mol: &mut Mol,
    params: &DetermineConnectivityParams,
) -> Result<()> {
    Ok(ffi::determine_connectivity(
        mol.pin_mut(),
        params.use_hueckel,
        params.charge,
        params.cov_factor,
        params.use_vdw,
    )?)
}

/// Parameters for `determine_bond_orders_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DetermineBondOrdersParams {
    /// The molecule's charge.
    pub charge: i32,
    /// Places formal charges on atoms by their valence; otherwise radical
    /// electrons.
    pub allow_charged_fragments: bool,
    /// Sanitizes the molecule and assigns its stereochemistry from the
    /// conformer.
    pub embed_chiral: bool,
    /// Numbers the atoms' map numbers from 1.
    pub use_atom_map: bool,
}

impl Default for DetermineBondOrdersParams {
    fn default() -> DetermineBondOrdersParams {
        DetermineBondOrdersParams {
            charge: 0,
            allow_charged_fragments: true,
            embed_chiral: true,
            use_atom_map: false,
        }
    }
}

/// Assigns orders to the molecule's bonds, which give its connectivity.
/// Mirrors `rdDetermineBonds.DetermineBondOrders`.
pub fn determine_bond_orders(mol: &mut Mol) -> Result<()> {
    determine_bond_orders_with_params(mol, &DetermineBondOrdersParams::default())
}

/// Mirrors `rdDetermineBonds.DetermineBondOrders(mol, charge,
/// allowChargedFragments, embedChiral, useAtomMap)`. No ordering that gives
/// the charge is `Error::Value`. RDKit 2024.03 tries every combination of
/// valences, without a limit.
pub fn determine_bond_orders_with_params(
    mol: &mut Mol,
    params: &DetermineBondOrdersParams,
) -> Result<()> {
    Ok(ffi::determine_bond_orders(
        mol.pin_mut(),
        params.charge,
        params.allow_charged_fragments,
        params.embed_chiral,
        params.use_atom_map,
    )?)
}

/// Parameters for `determine_bonds_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DetermineBondsParams {
    /// Uses extended Hückel theory (YAeHMOP) for the connectivity.
    pub use_hueckel: bool,
    /// The molecule's charge.
    pub charge: i32,
    /// The factor that scales covalent radii for the van der Waals method.
    pub cov_factor: f64,
    /// Places formal charges on atoms by their valence; otherwise radical
    /// electrons.
    pub allow_charged_fragments: bool,
    /// Sanitizes the molecule and assigns its stereochemistry from the
    /// conformer.
    pub embed_chiral: bool,
    /// Numbers the atoms' map numbers from 1.
    pub use_atom_map: bool,
    /// Uses the van der Waals method instead of connect-the-dots.
    pub use_vdw: bool,
}

impl Default for DetermineBondsParams {
    fn default() -> DetermineBondsParams {
        DetermineBondsParams {
            use_hueckel: false,
            charge: 0,
            cov_factor: 1.3,
            allow_charged_fragments: true,
            embed_chiral: true,
            use_atom_map: false,
            use_vdw: false,
        }
    }
}

/// Determines the molecule's connectivity, then its bond orders. Mirrors
/// `rdDetermineBonds.DetermineBonds`; a molecule of one atom or none is
/// left as it is.
pub fn determine_bonds(mol: &mut Mol) -> Result<()> {
    determine_bonds_with_params(mol, &DetermineBondsParams::default())
}

/// Mirrors `rdDetermineBonds.DetermineBonds(mol, useHueckel, charge,
/// covFactor, allowChargedFragments, embedChiral, useAtomMap, useVdw)`, with
/// the errors of `determine_connectivity_with_params`
/// and `determine_bond_orders_with_params`.
pub fn determine_bonds_with_params(mol: &mut Mol, params: &DetermineBondsParams) -> Result<()> {
    Ok(ffi::determine_bonds(
        mol.pin_mut(),
        params.use_hueckel,
        params.charge,
        params.cov_factor,
        params.allow_charged_fragments,
        params.embed_chiral,
        params.use_atom_map,
        params.use_vdw,
    )?)
}

/// Whether RDKit was built with YAeHMOP, which the Hückel option needs.
/// Mirrors `rdDetermineBonds.hueckelEnabled`.
pub fn hueckel_enabled() -> bool {
    expect_ok(ffi::hueckel_enabled())
}
