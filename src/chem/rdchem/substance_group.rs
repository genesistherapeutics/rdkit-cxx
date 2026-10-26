// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::pin::Pin;

use crate::Result;
use crate::chem::rdchem::Mol;
use crate::chem::rdchem::atom::expect_ok_index;
use crate::chem::rdchem::props::{prop_getters, prop_setters, strings};
use crate::error::{expect_ok, to_string};
use crate::ffi::props::ffi::{RDProps, props_get_string_vect_prop, props_get_unsigned_vect_prop};
use crate::ffi::substance_group::ffi;
use crate::geometry::Point3D;

/// A crossing bond's connection state. Mirrors `rdchem.SubstanceGroupCState`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SubstanceGroupCState {
    /// The crossing bond.
    pub bond_idx: u32,
    /// The bond's display vector.
    pub vector: Point3D,
}

/// Mirrors `SubstanceGroupCState()`.
impl Default for SubstanceGroupCState {
    fn default() -> SubstanceGroupCState {
        SubstanceGroupCState::from_ffi(expect_ok(ffi::substance_group_cstate_default()))
    }
}

impl SubstanceGroupCState {
    fn from_ffi(data: ffi::SubstanceGroupCStateData) -> SubstanceGroupCState {
        SubstanceGroupCState {
            bond_idx: data.bond_idx,
            vector: data.vector,
        }
    }
}

/// An attachment point. Mirrors `rdchem.SubstanceGroupAttach`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubstanceGroupAttach {
    /// The attached atom.
    pub a_idx: u32,
    /// The leaving atom, or -1 for none.
    pub lv_idx: i32,
    /// The attachment point's id.
    pub id: String,
}

/// Mirrors `SubstanceGroupAttach()`.
impl Default for SubstanceGroupAttach {
    fn default() -> SubstanceGroupAttach {
        SubstanceGroupAttach::from_ffi(expect_ok(ffi::substance_group_attach_default()))
            .expect("RDKit's default attachment id is empty")
    }
}

impl SubstanceGroupAttach {
    fn from_ffi(data: ffi::SubstanceGroupAttachData) -> Result<SubstanceGroupAttach> {
        Ok(SubstanceGroupAttach {
            a_idx: data.a_idx,
            lv_idx: data.lv_idx,
            id: to_string(data.id)?,
        })
    }
}

/// Stamps the substance group getters onto a type with
/// `fn raw(&self) -> &ffi::SubstanceGroup` and `fn props(&self) -> &RDProps`
/// (CONVENTIONS.md W3).
macro_rules! substance_group_getters {
    ($type:ident) => {
        impl $type<'_> {
            /// The group's index in its molecule. Mirrors
            /// `SubstanceGroup.GetIndexInMol`.
            pub fn get_index_in_mol(&self) -> u32 {
                expect_ok(ffi::substance_group_get_index_in_mol(self.raw()))
            }

            /// Mirrors `SubstanceGroup.GetAtoms`.
            pub fn get_atoms(&self) -> Vec<u32> {
                expect_ok(ffi::substance_group_get_atoms(self.raw()))
            }

            /// Mirrors `SubstanceGroup.GetParentAtoms`.
            pub fn get_parent_atoms(&self) -> Vec<u32> {
                expect_ok(ffi::substance_group_get_parent_atoms(self.raw()))
            }

            /// Mirrors `SubstanceGroup.GetBonds`.
            pub fn get_bonds(&self) -> Vec<u32> {
                expect_ok(ffi::substance_group_get_bonds(self.raw()))
            }

            /// Mirrors `SubstanceGroup.GetCStates`.
            pub fn get_c_states(&self) -> Vec<SubstanceGroupCState> {
                expect_ok(ffi::substance_group_get_cstates(self.raw()))
                    .into_iter()
                    .map(SubstanceGroupCState::from_ffi)
                    .collect()
            }

            /// Mirrors `SubstanceGroup.GetAttachPoints`.
            pub fn get_attach_points(&self) -> Result<Vec<SubstanceGroupAttach>> {
                expect_ok(ffi::substance_group_get_attach_points(self.raw()))
                    .into_iter()
                    .map(SubstanceGroupAttach::from_ffi)
                    .collect()
            }

            /// Each bracket's three points; a bracket added with two ends
            /// with the origin. Mirrors `SubstanceGroup.GetBrackets`.
            pub fn get_brackets(&self) -> Vec<[Point3D; 3]> {
                expect_ok(ffi::substance_group_get_brackets(self.raw()))
                    .chunks_exact(3)
                    .map(|points| [points[0], points[1], points[2]])
                    .collect()
            }

            /// Mirrors `SubstanceGroup.GetUnsignedVectProp`.
            pub fn get_unsigned_vect_prop(&self, key: &str) -> Result<Vec<u32>> {
                Ok(props_get_unsigned_vect_prop(self.props(), key)?)
            }

            /// Mirrors `SubstanceGroup.GetStringVectProp`.
            pub fn get_string_vect_prop(&self, key: &str) -> Result<Vec<String>> {
                strings(props_get_string_vect_prop(self.props(), key)?)
            }
        }

        prop_getters!(
            $type<'m>,
            dict_defaults = crate::chem::rdchem::GetPropsAsDictParams {
                include_private: true,
                include_computed: true,
                auto_convert_strings: true,
            },
            dict_doc = "All properties, private and computed ones included"
        );
    };
}

/// A substance group, borrowed from its molecule. Mirrors
/// `rdchem.SubstanceGroup`.
///
/// RDKit keeps a molecule's substance groups in the molecule, so there is no
/// owned substance group: `create_mol_substance_group` makes one in a
/// molecule, and `add_mol_substance_group` copies one into another.
#[derive(Clone, Copy)]
pub struct SubstanceGroupRef<'m> {
    mol: &'m Mol,
    sgroup: &'m ffi::SubstanceGroup,
}

impl<'m> SubstanceGroupRef<'m> {
    fn raw(&self) -> &'m ffi::SubstanceGroup {
        self.sgroup
    }

    fn props(&self) -> &'m RDProps {
        expect_ok(crate::ffi::props::ffi::substance_group_as_props(
            self.sgroup,
        ))
    }

    /// Mirrors `SubstanceGroup.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &'m Mol {
        self.mol
    }
}

substance_group_getters!(SubstanceGroupRef);

impl fmt::Debug for SubstanceGroupRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        debug(
            f,
            "SubstanceGroupRef",
            self.get_index_in_mol(),
            self.get_prop("TYPE").ok(),
        )
    }
}

/// A substance group, borrowed exclusively from its molecule to edit it.
/// Mirrors `rdchem.SubstanceGroup` (W2).
pub struct SubstanceGroupMut<'m> {
    mol: &'m mut Mol,
    idx: u32,
}

impl SubstanceGroupMut<'_> {
    fn raw(&self) -> &ffi::SubstanceGroup {
        expect_ok_index(
            ffi::mol_get_substance_group(self.mol.as_ffi(), self.idx).map_err(Into::into),
        )
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::SubstanceGroup> {
        expect_ok_index(
            ffi::mol_get_substance_group_mut(self.mol.pin_mut(), self.idx).map_err(Into::into),
        )
    }

    fn props(&self) -> &RDProps {
        expect_ok(crate::ffi::props::ffi::substance_group_as_props(self.raw()))
    }

    /// Mirrors `SubstanceGroup.GetOwningMol`.
    pub fn get_owning_mol(&self) -> &Mol {
        self.mol
    }

    /// Replaces the atoms; properties, connection states, and attachment
    /// points are left as they are. Mirrors `SubstanceGroup.SetAtoms`.
    pub fn set_atoms(&mut self, atoms: &[u32]) -> Result<()> {
        Ok(ffi::substance_group_set_atoms(self.raw_mut(), atoms)?)
    }

    /// Replaces the parent atoms. Mirrors `SubstanceGroup.SetParentAtoms`.
    pub fn set_parent_atoms(&mut self, patoms: &[u32]) -> Result<()> {
        Ok(ffi::substance_group_set_parent_atoms(
            self.raw_mut(),
            patoms,
        )?)
    }

    /// Replaces the bonds. Mirrors `SubstanceGroup.SetBonds`.
    pub fn set_bonds(&mut self, bonds: &[u32]) -> Result<()> {
        Ok(ffi::substance_group_set_bonds(self.raw_mut(), bonds)?)
    }

    /// Mirrors `SubstanceGroup.AddAtomWithIdx`.
    pub fn add_atom_with_idx(&mut self, idx: u32) -> Result<()> {
        Ok(ffi::substance_group_add_atom_with_idx(self.raw_mut(), idx)?)
    }

    /// Adds a parent atom, which must be one of the group's atoms. Mirrors
    /// `SubstanceGroup.AddParentAtomWithIdx`.
    pub fn add_parent_atom_with_idx(&mut self, idx: u32) -> Result<()> {
        Ok(ffi::substance_group_add_parent_atom_with_idx(
            self.raw_mut(),
            idx,
        )?)
    }

    /// Mirrors `SubstanceGroup.AddBondWithIdx`.
    pub fn add_bond_with_idx(&mut self, idx: u32) -> Result<()> {
        Ok(ffi::substance_group_add_bond_with_idx(self.raw_mut(), idx)?)
    }

    /// Adds the atom with an RDKit bookmark; a missing bookmark is
    /// `Error::Runtime`. Mirrors `SubstanceGroup.AddAtomWithBookmark`.
    pub fn add_atom_with_bookmark(&mut self, mark: i32) -> Result<()> {
        Ok(ffi::substance_group_add_atom_with_bookmark(
            self.raw_mut(),
            mark,
        )?)
    }

    /// Mirrors `SubstanceGroup.AddParentAtomWithBookmark`.
    pub fn add_parent_atom_with_bookmark(&mut self, mark: i32) -> Result<()> {
        Ok(ffi::substance_group_add_parent_atom_with_bookmark(
            self.raw_mut(),
            mark,
        )?)
    }

    /// Mirrors `SubstanceGroup.AddBondWithBookmark`.
    pub fn add_bond_with_bookmark(&mut self, mark: i32) -> Result<()> {
        Ok(ffi::substance_group_add_bond_with_bookmark(
            self.raw_mut(),
            mark,
        )?)
    }

    /// Adds a connection state to one of the group's crossing bonds. Mirrors
    /// `SubstanceGroup.AddCState`.
    pub fn add_c_state(&mut self, bond_idx: u32, vector: Point3D) -> Result<()> {
        Ok(ffi::substance_group_add_cstate(
            self.raw_mut(),
            bond_idx,
            &vector,
        )?)
    }

    /// Adds an attachment point; `lv_idx` is the leaving atom, or -1 for
    /// none. Mirrors `SubstanceGroup.AddAttachPoint`.
    pub fn add_attach_point(&mut self, a_idx: u32, lv_idx: i32, id_str: &str) -> Result<()> {
        Ok(ffi::substance_group_add_attach_point(
            self.raw_mut(),
            a_idx,
            lv_idx,
            id_str,
        )?)
    }

    /// Adds a bracket of two or three points; any other number is
    /// `Error::Value`. Mirrors `SubstanceGroup.AddBracket`.
    pub fn add_bracket(&mut self, pts: &[Point3D]) -> Result<()> {
        Ok(ffi::substance_group_add_bracket(self.raw_mut(), pts)?)
    }

    /// Mirrors `SubstanceGroup.ClearBrackets`.
    pub fn clear_brackets(&mut self) {
        expect_ok(ffi::substance_group_clear_brackets(self.raw_mut()))
    }

    /// Mirrors `SubstanceGroup.ClearCStates`.
    pub fn clear_c_states(&mut self) {
        expect_ok(ffi::substance_group_clear_cstates(self.raw_mut()))
    }

    /// Mirrors `SubstanceGroup.ClearAttachPoints`.
    pub fn clear_attach_points(&mut self) {
        expect_ok(ffi::substance_group_clear_attach_points(self.raw_mut()))
    }
}

substance_group_getters!(SubstanceGroupMut);
prop_setters!(SubstanceGroupMut<'m>);

impl fmt::Debug for SubstanceGroupMut<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        debug(f, "SubstanceGroupMut", self.idx, self.get_prop("TYPE").ok())
    }
}

fn debug(f: &mut fmt::Formatter<'_>, name: &str, idx: u32, type_: Option<String>) -> fmt::Result {
    f.debug_struct(name)
        .field("idx", &idx)
        .field("type", &type_)
        .finish()
}

/// The molecule's substance groups. Mirrors `rdchem.GetMolSubstanceGroups`,
/// which returns copies.
pub fn get_mol_substance_groups(mol: &Mol) -> Vec<SubstanceGroupRef<'_>> {
    let num = expect_ok(ffi::mol_get_num_substance_groups(mol.as_ffi()));
    (0..num)
        .map(|idx| expect_ok_index(get_mol_substance_group_with_idx(mol, idx)))
        .collect()
}

/// Mirrors `rdchem.GetMolSubstanceGroupWithIdx`.
pub fn get_mol_substance_group_with_idx(mol: &Mol, idx: u32) -> Result<SubstanceGroupRef<'_>> {
    let sgroup = ffi::mol_get_substance_group(mol.as_ffi(), idx)?;
    Ok(SubstanceGroupRef { mol, sgroup })
}

/// As [`get_mol_substance_group_with_idx`], to edit the group.
/// Mirrors `rdchem.GetMolSubstanceGroupWithIdx` (A6).
pub fn get_mol_substance_group_with_idx_mut(
    mol: &mut Mol,
    idx: u32,
) -> Result<SubstanceGroupMut<'_>> {
    ffi::mol_get_substance_group(mol.as_ffi(), idx)?;
    Ok(SubstanceGroupMut { mol, idx })
}

/// Removes the molecule's substance groups. Mirrors
/// `rdchem.ClearMolSubstanceGroups`.
pub fn clear_mol_substance_groups(mol: &mut Mol) {
    expect_ok(ffi::mol_clear_substance_groups(mol.pin_mut()))
}

/// Adds an empty substance group of type `type_` (`"SRU"`, `"DAT"`, ...).
/// Mirrors `rdchem.CreateMolSubstanceGroup`.
pub fn create_mol_substance_group<'m>(mol: &'m mut Mol, type_: &str) -> SubstanceGroupMut<'m> {
    let idx = expect_ok(ffi::mol_create_substance_group(mol.pin_mut(), type_));
    SubstanceGroupMut { mol, idx }
}

/// Adds a data substance group holding one field. Mirrors
/// `rdchem.CreateMolDataSubstanceGroup`.
pub fn create_mol_data_substance_group<'m>(
    mol: &'m mut Mol,
    field_name: &str,
    value: &str,
) -> SubstanceGroupMut<'m> {
    let idx = expect_ok(ffi::mol_create_data_substance_group(
        mol.pin_mut(),
        field_name,
        value,
    ));
    SubstanceGroupMut { mol, idx }
}

/// Adds a copy of another molecule's substance group, whose atom and bond
/// indices must exist in `mol`. Mirrors `rdchem.AddMolSubstanceGroup`.
///
/// To copy a group within one molecule, copy it from a clone of the molecule.
pub fn add_mol_substance_group<'m>(
    mol: &'m mut Mol,
    sgroup: &SubstanceGroupRef<'_>,
) -> Result<SubstanceGroupMut<'m>> {
    let idx = ffi::mol_add_substance_group(mol.pin_mut(), sgroup.raw())?;
    Ok(SubstanceGroupMut { mol, idx })
}
