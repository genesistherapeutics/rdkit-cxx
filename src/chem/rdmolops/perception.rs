// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ring finding, aromaticity, cleanup, and the other perception steps.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::AsAtom;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::perception::ffi;

rdkit_enum! {
    /// An aromaticity model. Mirrors `rdmolops.AromaticityModel`.
    AromaticityModel: ffi::AromaticityModel {
        AromaticityDefault = AROMATICITY_DEFAULT,
        AromaticityRdkit = AROMATICITY_RDKIT,
        AromaticitySimple = AROMATICITY_SIMPLE,
        AromaticityMdl = AROMATICITY_MDL,
        AromaticityMmff94 = AROMATICITY_MMFF94,
        AromaticityCustom = AROMATICITY_CUSTOM,
    }
}

/// The smallest set of smallest rings, as atom indices; also stores them in
/// the molecule's ring info. Mirrors `rdmolops.GetSSSR`.
pub fn get_sssr(mol: &mut Mol) -> Result<Vec<Vec<u32>>> {
    get_sssr_with_include_dative_bonds(mol, false)
}

/// Mirrors `rdmolops.GetSSSR(mol, includeDativeBonds)`: also uses dative
/// bonds in ring finding.
pub fn get_sssr_with_include_dative_bonds(
    mol: &mut Mol,
    include_dative_bonds: bool,
) -> Result<Vec<Vec<u32>>> {
    Ok(IndexList::into_nested(ffi::get_sssr(
        mol.pin_mut(),
        include_dative_bonds,
    )?))
}

/// The symmetrized SSSR, as atom indices: the SSSR plus the rings that make
/// it symmetric, as for cubane's sixth ring. Also stores them in the ring
/// info. Mirrors `rdmolops.GetSymmSSSR`.
pub fn get_symm_sssr(mol: &mut Mol) -> Result<Vec<Vec<u32>>> {
    get_symm_sssr_with_include_dative_bonds(mol, false)
}

/// Mirrors `rdmolops.GetSymmSSSR(mol, includeDativeBonds)`: also uses
/// dative bonds in ring finding.
pub fn get_symm_sssr_with_include_dative_bonds(
    mol: &mut Mol,
    include_dative_bonds: bool,
) -> Result<Vec<Vec<u32>>> {
    Ok(IndexList::into_nested(ffi::get_symm_sssr(
        mol.pin_mut(),
        include_dative_bonds,
    )?))
}

/// Finds rings quickly, for ring membership only: the ring info answers
/// whether an atom or bond is in a ring, and little else. Mirrors
/// `rdmolops.FastFindRings`.
pub fn fast_find_rings(mol: &mut Mol) -> Result<()> {
    Ok(ffi::fast_find_rings(mol.pin_mut())?)
}

/// Finds the molecule's ring families and stores them in its ring info.
/// Mirrors `rdmolops.FindRingFamilies`.
pub fn find_ring_families(mol: &mut Mol) -> Result<()> {
    Ok(ffi::find_ring_families(mol.pin_mut())?)
}

/// Removes common mistakes in drawn structures: charge-separated nitro
/// groups, five-valent nitrogens, and the like. Mirrors `rdmolops.Cleanup`.
pub fn cleanup(mol: &mut Mol) -> Result<()> {
    Ok(ffi::cleanup(mol.pin_mut())?)
}

/// Removes chirality from atoms that can't be stereocenters. Mirrors
/// `rdmolops.CleanupChirality`.
pub fn cleanup_chirality(mol: &mut Mol) -> Result<()> {
    Ok(ffi::cleanup_chirality(mol.pin_mut())?)
}

/// Removes atropisomer stereo from bonds that can't hold it. Mirrors
/// `rdmolops.CleanupAtropisomers`.
pub fn cleanup_atropisomers(mol: &mut Mol) -> Result<()> {
    Ok(ffi::cleanup_atropisomers(mol.pin_mut())?)
}

/// Converts bonds to metals that would break valence rules into dative
/// bonds. Mirrors `rdmolops.CleanupOrganometallics`.
pub fn cleanup_organometallics(mol: &mut Mol) -> Result<()> {
    Ok(ffi::cleanup_organometallics(mol.pin_mut())?)
}

/// Perceives aromaticity with RDKit's default model, on a molecule already
/// kekulized. Mirrors `rdmolops.SetAromaticity`.
pub fn set_aromaticity(mol: &mut Mol) -> Result<()> {
    set_aromaticity_with_model(mol, AromaticityModel::AromaticityDefault)
}

/// Perceives aromaticity with `model`. `AromaticityCustom` needs a C++
/// function and is `Error::Runtime`, as in Python. Mirrors
/// `rdmolops.SetAromaticity(mol, model)`.
pub fn set_aromaticity_with_model(mol: &mut Mol, model: AromaticityModel) -> Result<()> {
    Ok(ffi::set_aromaticity(mol.pin_mut(), model.to_ffi())?)
}

/// Flags the molecule's conjugated bonds. Mirrors `rdmolops.SetConjugation`.
pub fn set_conjugation(mol: &mut Mol) -> Result<()> {
    Ok(ffi::set_conjugation(mol.pin_mut())?)
}

/// Assigns each atom's hybridization. Mirrors `rdmolops.SetHybridization`.
pub fn set_hybridization(mol: &mut Mol) -> Result<()> {
    Ok(ffi::set_hybridization(mol.pin_mut())?)
}

/// Assigns radical electrons to atoms with too few bonds. Mirrors
/// `rdmolops.AssignRadicals`.
pub fn assign_radicals(mol: &mut Mol) -> Result<()> {
    Ok(ffi::assign_radicals(mol.pin_mut())?)
}

/// Kekulizes the molecule, leaving it as it was if that fails. Mirrors
/// `rdmolops.KekulizeIfPossible`.
pub fn kekulize_if_possible(mol: &mut Mol) -> Result<()> {
    kekulize_if_possible_with_clear_aromatic_flags(mol, false)
}

/// Mirrors `rdmolops.KekulizeIfPossible(mol, clearAromaticFlags)`: also
/// clears the aromatic flags on atoms and bonds.
pub fn kekulize_if_possible_with_clear_aromatic_flags(
    mol: &mut Mol,
    clear_aromatic_flags: bool,
) -> Result<()> {
    Ok(ffi::kekulize_if_possible(
        mol.pin_mut(),
        clear_aromatic_flags,
    )?)
}

/// A copy with each haptic bond, drawn as a dummy atom with a dative bond,
/// replaced by dative bonds from the atoms it stands for. Mirrors
/// `rdmolops.HapticBondsToDative`.
pub fn haptic_bonds_to_dative(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::haptic_bonds_to_dative(mol.as_ffi())?))
}

/// The reverse of `haptic_bonds_to_dative`. Mirrors
/// `rdmolops.DativeBondsToHaptic`.
pub fn dative_bonds_to_haptic(mol: &Mol) -> Result<Mol> {
    Ok(Mol::from_ffi(ffi::dative_bonds_to_haptic(mol.as_ffi())?))
}

/// Whether the molecule has implicit hydrogens to add. Mirrors
/// `rdmolops.NeedsHs`.
pub fn needs_hs(mol: &Mol) -> Result<bool> {
    Ok(ffi::needs_hs(mol.as_ffi())?)
}

/// The electrons an atom in a molecule can give to an aromatic system.
/// Mirrors `rdmolops.CountAtomElec`.
pub fn count_atom_elec(atom: &impl AsAtom) -> Result<i32> {
    Ok(ffi::count_atom_elec(atom.atom())?)
}

/// Whether an atom in a molecule has a conjugated bond. Mirrors
/// `rdmolops.AtomHasConjugatedBond`.
pub fn atom_has_conjugated_bond(atom: &impl AsAtom) -> Result<bool> {
    Ok(ffi::atom_has_conjugated_bond(atom.atom())?)
}

/// The sum of the atoms' formal charges. RDKit sums them in `int`, which
/// overflows only for molecules of more than 2^24 atoms (each charge is 8
/// bits); a molecule whose running sum would leave `i32`'s range is
/// `Error::Value`. Mirrors `rdmolops.GetFormalCharge`.
pub fn get_formal_charge(mol: &Mol) -> Result<i32> {
    Ok(ffi::get_formal_charge(mol.as_ffi())?)
}

/// The atoms on a shortest path between two atoms, both included; empty
/// when they are not connected. Mirrors `rdmolops.GetShortestPath`.
pub fn get_shortest_path(mol: &Mol, aid1: u32, aid2: u32) -> Result<Vec<u32>> {
    Ok(ffi::get_shortest_path(mol.as_ffi(), aid1, aid2)?)
}

/// Places terminal atom `idx` next to its neighbor `other_idx` in every
/// conformer, by the neighbor's geometry. Mirrors
/// `rdmolops.SetTerminalAtomCoords`.
pub fn set_terminal_atom_coords(mol: &mut Mol, idx: u32, other_idx: u32) -> Result<()> {
    Ok(ffi::set_terminal_atom_coords(
        mol.pin_mut(),
        idx,
        other_idx,
    )?)
}
