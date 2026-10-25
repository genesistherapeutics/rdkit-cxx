// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Lipinski`, a pure-Python module: Lipinski's
//! parameters and the SMARTS patterns behind them.
//!
//! Python builds the patterns once, on import; each pattern function here
//! parses a new molecule, since a `Mol` cannot be shared between threads
//! (CONVENTIONS.md S1). The descriptor functions call
//! `rdkit_cxx::chem::rd_mol_descriptors`, as Python's do. Python's private
//! `_HDonors`-style matchers and the functions' `version` attributes are not
//! ported.

use crate::Result;
use crate::chem::rd_mol_descriptors as descriptors;
use crate::chem::{Mol, mol_from_smarts};

fn pattern(smarts: &str) -> Mol {
    mol_from_smarts(smarts).expect("the Lipinski SMARTS parse")
}

/// Hydrogen-bond donors, adapted from the Gobbi paper. Mirrors
/// `Lipinski.HDonorSmarts`.
pub fn h_donor_smarts() -> Mol {
    pattern("[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]")
}

/// Hydrogen-bond acceptors, adapted from the Gobbi paper. For traditional
/// Lipinski numbers use [`no_count`] instead. Mirrors
/// `Lipinski.HAcceptorSmarts`.
pub fn h_acceptor_smarts() -> Mol {
    pattern(concat!(
        "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),",
        "$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),",
        "$([nH0,o,s;+0])]"
    ))
}

/// Atoms other than carbon and hydrogen. Mirrors `Lipinski.HeteroatomSmarts`.
pub fn heteroatom_smarts() -> Mol {
    pattern("[!#6;!#1]")
}

/// Rotatable bonds. Deuteriums, which stay in the graph, count toward a
/// carbon's degree where hydrogens do not, so the bond in
/// `[2H]C([2H])([2H])C([2H])([2H])[2H]` is rotatable. Mirrors
/// `Lipinski.RotatableBondSmarts`.
pub fn rotatable_bond_smarts() -> Mol {
    pattern("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")
}

/// NH and OH groups. Mirrors `Lipinski.NHOHSmarts`.
pub fn nhoh_smarts() -> Mol {
    pattern("[#8H1,#7H1,#7H2,#7H3]")
}

/// Nitrogens and oxygens. Mirrors `Lipinski.NOCountSmarts`.
pub fn no_count_smarts() -> Mol {
    pattern("[#7,#8]")
}

/// The number of hydrogen-bond donors. Mirrors `Lipinski.NumHDonors`.
pub fn num_h_donors(mol: &Mol) -> Result<u32> {
    descriptors::calc_num_hbd(mol)
}

/// The number of hydrogen-bond acceptors. Mirrors `Lipinski.NumHAcceptors`.
pub fn num_h_acceptors(mol: &Mol) -> Result<u32> {
    descriptors::calc_num_hba(mol)
}

/// The number of heteroatoms. Mirrors `Lipinski.NumHeteroatoms`.
pub fn num_heteroatoms(mol: &Mol) -> Result<u32> {
    descriptors::calc_num_heteroatoms(mol)
}

/// The number of rotatable bonds. Mirrors `Lipinski.NumRotatableBonds`.
pub fn num_rotatable_bonds(mol: &Mol) -> Result<u32> {
    descriptors::calc_num_rotatable_bonds(mol)
}

/// The number of nitrogens and oxygens. Mirrors `Lipinski.NOCount`.
pub fn no_count(mol: &Mol) -> u32 {
    descriptors::calc_num_lipinski_hba(mol)
}

/// The number of NHs and OHs. Mirrors `Lipinski.NHOHCount`.
pub fn nhoh_count(mol: &Mol) -> Result<u32> {
    descriptors::calc_num_lipinski_hbd(mol)
}

/// The number of rings; `Error::Runtime` before rings are found. Mirrors
/// `Lipinski.RingCount`.
pub fn ring_count(mol: &Mol) -> Result<u32> {
    descriptors::calc_num_rings(mol)
}

/// The number of heavy atoms. Mirrors `Lipinski.HeavyAtomCount`.
pub fn heavy_atom_count(mol: &Mol) -> u32 {
    mol.get_num_heavy_atoms()
}

/// Defines the functions Python's `Lipinski` makes from
/// `rdMolDescriptors.Calc*` functions by dropping `Calc` from their names.
macro_rules! bulk_convert {
    ($($name:ident = $calc:ident -> $out:ty, $python:literal;)*) => {
        $(
            #[doc = concat!("Mirrors `Lipinski.", $python, "`, which calls `rdMolDescriptors.Calc", $python, "`.")]
            pub fn $name(mol: &Mol) -> $out {
                descriptors::$calc(mol)
            }
        )*
    };
}

bulk_convert! {
    fraction_csp3 = calc_fraction_csp3 -> Result<f64>, "FractionCSP3";
    num_aromatic_rings = calc_num_aromatic_rings -> u32, "NumAromaticRings";
    num_saturated_rings = calc_num_saturated_rings -> u32, "NumSaturatedRings";
    num_aromatic_heterocycles = calc_num_aromatic_heterocycles -> u32, "NumAromaticHeterocycles";
    num_aromatic_carbocycles = calc_num_aromatic_carbocycles -> u32, "NumAromaticCarbocycles";
    num_saturated_heterocycles = calc_num_saturated_heterocycles -> u32, "NumSaturatedHeterocycles";
    num_saturated_carbocycles = calc_num_saturated_carbocycles -> u32, "NumSaturatedCarbocycles";
    num_aliphatic_rings = calc_num_aliphatic_rings -> u32, "NumAliphaticRings";
    num_aliphatic_heterocycles = calc_num_aliphatic_heterocycles -> u32, "NumAliphaticHeterocycles";
    num_aliphatic_carbocycles = calc_num_aliphatic_carbocycles -> u32, "NumAliphaticCarbocycles";
    num_heterocycles = calc_num_heterocycles -> u32, "NumHeterocycles";
    num_bridgehead_atoms = calc_num_bridgehead_atoms -> Result<u32>, "NumBridgeheadAtoms";
    num_amide_bonds = calc_num_amide_bonds -> Result<u32>, "NumAmideBonds";
    num_atom_stereo_centers = calc_num_atom_stereo_centers -> Result<u32>, "NumAtomStereoCenters";
    num_unspecified_atom_stereo_centers = calc_num_unspecified_atom_stereo_centers -> Result<u32>, "NumUnspecifiedAtomStereoCenters";
    num_spiro_atoms = calc_num_spiro_atoms -> Result<u32>, "NumSpiroAtoms";
    phi = calc_phi -> Result<f64>, "Phi";
}
