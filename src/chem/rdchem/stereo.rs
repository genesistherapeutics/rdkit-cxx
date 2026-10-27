// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Result;
use crate::chem::rdchem::{
    Atom, Mol, RWMol, StereoDescriptor, StereoGroupType, StereoSpecified, StereoType,
};
use crate::error::expect_ok;
use crate::ffi::stereo::ffi;

/// A set of atoms and bonds with a defined stereochemical relationship.
/// Mirrors `rdchem.StereoGroup`.
///
/// RDKit's stereo groups point into their molecule; this one holds indices,
/// so it is a plain value that `RWMol::set_stereo_groups` checks against the
/// molecule it is given. Python's `GetAtoms` returns atoms; `get_atoms`
/// returns their indices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoGroup {
    group_type: StereoGroupType,
    atoms: Vec<u32>,
    bonds: Vec<u32>,
    read_id: u32,
    write_id: u32,
}

impl StereoGroup {
    fn from_ffi(data: ffi::StereoGroupData) -> StereoGroup {
        StereoGroup {
            group_type: StereoGroupType::from_ffi(data.group_type),
            atoms: data.atoms,
            bonds: data.bonds,
            read_id: data.read_id,
            write_id: data.write_id,
        }
    }

    fn to_ffi(&self) -> ffi::StereoGroupData {
        ffi::StereoGroupData {
            group_type: self.group_type.to_ffi(),
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
            read_id: self.read_id,
            write_id: self.write_id,
        }
    }

    /// Mirrors `StereoGroup.GetGroupType`.
    pub fn get_group_type(&self) -> StereoGroupType {
        self.group_type
    }

    /// The group's atom indices. Mirrors `StereoGroup.GetAtoms`.
    pub fn get_atoms(&self) -> &[u32] {
        &self.atoms
    }

    /// The group's bond indices. Mirrors `StereoGroup.GetBonds`.
    pub fn get_bonds(&self) -> &[u32] {
        &self.bonds
    }

    /// The id the group was read with; meaningful for AND and OR groups.
    /// Mirrors `StereoGroup.GetReadId`.
    pub fn get_read_id(&self) -> u32 {
        self.read_id
    }

    /// The id the group will be written with. Mirrors `StereoGroup.GetWriteId`.
    pub fn get_write_id(&self) -> u32 {
        self.write_id
    }

    /// Mirrors `StereoGroup.SetWriteId`.
    pub fn set_write_id(&mut self, id: u32) {
        self.write_id = id;
    }
}

/// Parameters for `create_stereo_group_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CreateStereoGroupParams {
    /// The group's atoms.
    pub atom_ids: Vec<u32>,
    /// The group's bonds.
    pub bond_ids: Vec<u32>,
    /// The group's id as read from a file; 0 for none.
    pub read_id: u32,
}

/// A stereo group of `atom_ids` in `mol`. Mirrors `rdchem.CreateStereoGroup`,
/// whose atoms are optional; a group needs at least one atom or bond, so they
/// are required here.
pub fn create_stereo_group(
    stereo_group_type: StereoGroupType,
    mol: &Mol,
    atom_ids: &[u32],
) -> Result<StereoGroup> {
    let params = CreateStereoGroupParams {
        atom_ids: atom_ids.to_vec(),
        ..Default::default()
    };
    create_stereo_group_with_params(stereo_group_type, mol, &params)
}

/// Mirrors `rdchem.CreateStereoGroup(stereoGroupType, mol, atomIds, bondIds,
/// readId)`.
pub fn create_stereo_group_with_params(
    stereo_group_type: StereoGroupType,
    mol: &Mol,
    params: &CreateStereoGroupParams,
) -> Result<StereoGroup> {
    Ok(StereoGroup::from_ffi(ffi::create_stereo_group(
        stereo_group_type.to_ffi(),
        mol.as_ffi(),
        &params.atom_ids,
        &params.bond_ids,
        params.read_id,
    )?))
}

/// Uses each stereo group's read id as its write id. Mirrors
/// `rdchem.ForwardStereoGroupIds`.
pub fn forward_stereo_group_ids(mol: &mut Mol) {
    expect_ok(ffi::forward_stereo_group_ids(mol.pin_mut()))
}

impl Mol {
    /// Copies of the molecule's stereo groups. Mirrors `Mol.GetStereoGroups`.
    pub fn get_stereo_groups(&self) -> Vec<StereoGroup> {
        expect_ok(ffi::mol_get_stereo_groups(self.as_ffi()))
            .into_iter()
            .map(StereoGroup::from_ffi)
            .collect()
    }
}

impl RWMol {
    /// Replaces the stereo groups. Mirrors `RWMol.SetStereoGroups`.
    pub fn set_stereo_groups(&mut self, stereo_groups: &[StereoGroup]) -> Result<()> {
        let data: Vec<_> = stereo_groups.iter().map(StereoGroup::to_ffi).collect();
        Ok(ffi::rwmol_set_stereo_groups(self.pin_mut(), &data)?)
    }
}

/// A potential stereo element. Mirrors `rdchem.StereoInfo`, which
/// `rdmolops.FindPotentialStereo` returns.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct StereoInfo {
    /// The kind of stereo element. Mirrors `StereoInfo.type`.
    pub type_: StereoType,
    /// Whether the element's stereo is specified.
    pub specified: StereoSpecified,
    /// The atom or bond index, or `StereoInfo::NOATOM`.
    pub centered_on: u32,
    /// The element's stereo descriptor.
    pub descriptor: StereoDescriptor,
    /// For the non-tetrahedral cases.
    pub permutation: u32,
    /// The atoms around the atom or bond.
    pub controlling_atoms: Vec<u32>,
}

/// An unspecified element centered on nothing. Mirrors `StereoInfo()`.
impl Default for StereoInfo {
    fn default() -> StereoInfo {
        StereoInfo {
            type_: StereoType::Unspecified,
            specified: StereoSpecified::Unspecified,
            centered_on: StereoInfo::NOATOM,
            descriptor: StereoDescriptor::NoValue,
            permutation: 0,
            controlling_atoms: Vec::new(),
        }
    }
}

impl StereoInfo {
    /// The marker for unset indices. Mirrors `StereoInfo.NOATOM`.
    pub const NOATOM: u32 = u32::MAX;

    pub(crate) fn from_ffi(data: ffi::StereoInfoData) -> StereoInfo {
        StereoInfo {
            type_: StereoType::from_ffi(data.type_),
            specified: StereoSpecified::from_ffi(data.specified),
            centered_on: data.centered_on,
            descriptor: StereoDescriptor::from_ffi(data.descriptor),
            permutation: data.permutation,
            controlling_atoms: data.controlling_atoms,
        }
    }
}
