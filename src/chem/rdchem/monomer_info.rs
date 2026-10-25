// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::error::{expect_ok, to_string};
use crate::ffi::monomer_info::ffi;

/// The kind of monomer information. Mirrors `rdchem.AtomMonomerType`.
///
/// RDKit's enum has no fixed underlying type, so cxx cannot check it; it
/// crosses as `i32`, and a test sends every variant through C++.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AtomMonomerType {
    /// Mirrors `UNKNOWN`.
    Unknown,
    /// Mirrors `PDBRESIDUE`.
    Pdbresidue,
    /// Mirrors `OTHER`.
    Other,
}

impl AtomMonomerType {
    fn to_ffi(self) -> i32 {
        match self {
            AtomMonomerType::Unknown => 0,
            AtomMonomerType::Pdbresidue => 1,
            AtomMonomerType::Other => 2,
        }
    }

    fn from_ffi(value: i32) -> AtomMonomerType {
        match value {
            0 => AtomMonomerType::Unknown,
            1 => AtomMonomerType::Pdbresidue,
            2 => AtomMonomerType::Other,
            _ => panic!("unknown AtomMonomerType value {value}"),
        }
    }
}

/// Stamps the monomer getters onto a type with
/// `fn monomer(&self) -> &ffi::AtomMonomerInfo` (CONVENTIONS.md W3).
macro_rules! monomer_getters {
    ($type:ident) => {
        impl $type {
            /// Mirrors `GetName`.
            pub fn get_name(&self) -> Result<String> {
                to_string(ffi::monomer_info_get_name(self.monomer())?)
            }

            /// Mirrors `GetMonomerType`.
            pub fn get_monomer_type(&self) -> AtomMonomerType {
                AtomMonomerType::from_ffi(expect_ok(ffi::monomer_info_get_monomer_type(
                    self.monomer(),
                )))
            }

            /// Mirrors `GetResidueName`.
            pub fn get_residue_name(&self) -> Result<String> {
                to_string(ffi::monomer_info_get_residue_name(self.monomer())?)
            }

            /// Mirrors `GetResidueNumber`.
            pub fn get_residue_number(&self) -> i32 {
                expect_ok(ffi::monomer_info_get_residue_number(self.monomer()))
            }

            /// Mirrors `GetChainId`.
            pub fn get_chain_id(&self) -> Result<String> {
                to_string(ffi::monomer_info_get_chain_id(self.monomer())?)
            }

            /// Mirrors `GetMonomerClass`.
            pub fn get_monomer_class(&self) -> Result<String> {
                to_string(ffi::monomer_info_get_monomer_class(self.monomer())?)
            }
        }
    };
}

/// Stamps the monomer setters onto a type with
/// `fn monomer_mut(&mut self) -> Pin<&mut ffi::AtomMonomerInfo>` (W3).
macro_rules! monomer_setters {
    ($type:ident) => {
        impl $type {
            /// Mirrors `SetName`.
            pub fn set_name(&mut self, nm: &str) {
                expect_ok(ffi::monomer_info_set_name(self.monomer_mut(), nm))
            }

            /// Mirrors `SetResidueName`.
            pub fn set_residue_name(&mut self, val: &str) {
                expect_ok(ffi::monomer_info_set_residue_name(self.monomer_mut(), val))
            }

            /// Mirrors `SetResidueNumber`.
            pub fn set_residue_number(&mut self, val: i32) {
                expect_ok(ffi::monomer_info_set_residue_number(
                    self.monomer_mut(),
                    val,
                ))
            }

            /// Mirrors `SetChainId`.
            pub fn set_chain_id(&mut self, val: &str) {
                expect_ok(ffi::monomer_info_set_chain_id(self.monomer_mut(), val))
            }

            /// Mirrors `SetMonomerClass`.
            pub fn set_monomer_class(&mut self, val: &str) {
                expect_ok(ffi::monomer_info_set_monomer_class(self.monomer_mut(), val))
            }
        }
    };
}

/// Parameters for `AtomMonomerInfo::new_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AtomMonomerInfoParams {
    pub name: String,
    pub residue_name: String,
    pub res_num: i32,
    pub chain_id: String,
    pub monomer_class: String,
}

/// An atom's monomer information. Mirrors `rdchem.AtomMonomerInfo`.
///
/// Atoms hand out copies (`get_monomer_info`) and take copies
/// (`set_monomer_info`); edit a copy and set it back, where Python edits the
/// atom's info in place. A copy of a PDB residue's info keeps its PDB fields
/// for `set_monomer_info`, though only `get_pdb_residue_info` reads them.
///
/// Only such a copy, or an `AtomPDBResidueInfo`, has type
/// `AtomMonomerType::Pdbresidue`: RDKit casts any info of that type to PDB
/// residue info, reading past the end of a plain one, so `new` and
/// `set_monomer_type` refuse it, where Python allows it (CONVENTIONS.md U10).
pub struct AtomMonomerInfo {
    /// Never null.
    ptr: UniquePtr<ffi::AtomMonomerInfo>,
}

impl AtomMonomerInfo {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::AtomMonomerInfo>) -> AtomMonomerInfo {
        assert!(!ptr.is_null(), "RDKit returned null monomer info");
        AtomMonomerInfo { ptr }
    }

    pub(crate) fn monomer(&self) -> &ffi::AtomMonomerInfo {
        self.ptr.as_ref().unwrap()
    }

    fn monomer_mut(&mut self) -> Pin<&mut ffi::AtomMonomerInfo> {
        self.ptr.pin_mut()
    }

    /// Info of type `typ`, with Python's defaults for the rest. Mirrors
    /// `AtomMonomerInfo(type)`; `Error::Value` for `Pdbresidue`, for which
    /// `AtomPDBResidueInfo::new` makes PDB residue info.
    pub fn new(typ: AtomMonomerType) -> Result<AtomMonomerInfo> {
        AtomMonomerInfo::new_with_params(typ, &AtomMonomerInfoParams::default())
    }

    /// Mirrors `AtomMonomerInfo(type, name, residueName, resNum, chainId,
    /// monomerClass)`; `Error::Value` for `Pdbresidue`.
    pub fn new_with_params(
        typ: AtomMonomerType,
        params: &AtomMonomerInfoParams,
    ) -> Result<AtomMonomerInfo> {
        Ok(AtomMonomerInfo::from_ffi(ffi::monomer_info_with(
            typ.to_ffi(),
            &params.name,
            &params.residue_name,
            params.res_num,
            &params.chain_id,
            &params.monomer_class,
        )?))
    }

    /// Mirrors `SetMonomerType`; `Error::Value` for `Pdbresidue` unless the
    /// info is a copy of a PDB residue's.
    pub fn set_monomer_type(&mut self, typ: AtomMonomerType) -> Result<()> {
        Ok(ffi::monomer_info_set_monomer_type(
            self.monomer_mut(),
            typ.to_ffi(),
        )?)
    }
}

/// Info of type `Unknown`. Mirrors `AtomMonomerInfo()`.
impl Default for AtomMonomerInfo {
    fn default() -> AtomMonomerInfo {
        AtomMonomerInfo::from_ffi(expect_ok(ffi::monomer_info_new()))
    }
}

impl Clone for AtomMonomerInfo {
    fn clone(&self) -> AtomMonomerInfo {
        AtomMonomerInfo::from_ffi(expect_ok(ffi::monomer_info_clone(self.monomer())))
    }
}

impl fmt::Debug for AtomMonomerInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AtomMonomerInfo")
            .field("monomer_type", &self.get_monomer_type())
            .field("name", &self.get_name().ok())
            .finish()
    }
}

monomer_getters!(AtomMonomerInfo);
monomer_setters!(AtomMonomerInfo);

/// Parameters for `AtomPDBResidueInfo::new_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct AtomPDBResidueInfoParams {
    pub serial_number: i32,
    pub alt_loc: String,
    pub residue_name: String,
    pub residue_number: i32,
    pub chain_id: String,
    pub insertion_code: String,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub is_hetero_atom: bool,
    pub secondary_structure: u32,
    pub segment_number: u32,
    pub monomer_class: String,
}

impl Default for AtomPDBResidueInfoParams {
    fn default() -> AtomPDBResidueInfoParams {
        AtomPDBResidueInfoParams {
            serial_number: 1,
            alt_loc: String::new(),
            residue_name: String::new(),
            residue_number: 0,
            chain_id: String::new(),
            insertion_code: String::new(),
            occupancy: 1.0,
            temp_factor: 0.0,
            is_hetero_atom: false,
            secondary_structure: 0,
            segment_number: 0,
            monomer_class: String::new(),
        }
    }
}

/// An atom's PDB residue information. Mirrors `rdchem.AtomPDBResidueInfo`;
/// copied in and out of atoms like `AtomMonomerInfo`.
pub struct AtomPDBResidueInfo {
    /// Never null.
    ptr: UniquePtr<ffi::AtomPDBResidueInfo>,
}

impl AtomPDBResidueInfo {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::AtomPDBResidueInfo>) -> AtomPDBResidueInfo {
        assert!(!ptr.is_null(), "RDKit returned null residue info");
        AtomPDBResidueInfo { ptr }
    }

    fn raw(&self) -> &ffi::AtomPDBResidueInfo {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::AtomPDBResidueInfo> {
        self.ptr.pin_mut()
    }

    pub(crate) fn monomer(&self) -> &ffi::AtomMonomerInfo {
        expect_ok(ffi::pdb_residue_info_as_monomer_info(self.raw()))
    }

    fn monomer_mut(&mut self) -> Pin<&mut ffi::AtomMonomerInfo> {
        expect_ok(ffi::pdb_residue_info_as_monomer_info_mut(self.raw_mut()))
    }

    /// Info for the atom named `atom_name`, with Python's defaults for the
    /// rest. Mirrors `AtomPDBResidueInfo(atomName)`.
    pub fn new(atom_name: &str) -> AtomPDBResidueInfo {
        AtomPDBResidueInfo::new_with_params(atom_name, &AtomPDBResidueInfoParams::default())
    }

    /// Mirrors `AtomPDBResidueInfo(atomName, serialNumber, ...)`.
    pub fn new_with_params(atom_name: &str, p: &AtomPDBResidueInfoParams) -> AtomPDBResidueInfo {
        AtomPDBResidueInfo::from_ffi(expect_ok(ffi::pdb_residue_info_with(
            atom_name,
            p.serial_number,
            &p.alt_loc,
            &p.residue_name,
            p.residue_number,
            &p.chain_id,
            &p.insertion_code,
            p.occupancy,
            p.temp_factor,
            p.is_hetero_atom,
            p.secondary_structure,
            p.segment_number,
            &p.monomer_class,
        )))
    }

    /// Mirrors `AtomPDBResidueInfo.SetMonomerType`. Any type is safe here:
    /// RDKit casts only on `Pdbresidue`, and this info is PDB residue info.
    pub fn set_monomer_type(&mut self, typ: AtomMonomerType) {
        expect_ok(ffi::monomer_info_set_monomer_type(
            self.monomer_mut(),
            typ.to_ffi(),
        ))
    }

    /// Mirrors `AtomPDBResidueInfo.GetSerialNumber`.
    pub fn get_serial_number(&self) -> i32 {
        expect_ok(ffi::pdb_residue_info_get_serial_number(self.raw()))
    }

    /// Mirrors `AtomPDBResidueInfo.GetAltLoc`.
    pub fn get_alt_loc(&self) -> Result<String> {
        to_string(ffi::pdb_residue_info_get_alt_loc(self.raw())?)
    }

    /// Mirrors `AtomPDBResidueInfo.GetInsertionCode`.
    pub fn get_insertion_code(&self) -> Result<String> {
        to_string(ffi::pdb_residue_info_get_insertion_code(self.raw())?)
    }

    /// Mirrors `AtomPDBResidueInfo.GetOccupancy`.
    pub fn get_occupancy(&self) -> f64 {
        expect_ok(ffi::pdb_residue_info_get_occupancy(self.raw()))
    }

    /// Mirrors `AtomPDBResidueInfo.GetTempFactor`.
    pub fn get_temp_factor(&self) -> f64 {
        expect_ok(ffi::pdb_residue_info_get_temp_factor(self.raw()))
    }

    /// Mirrors `AtomPDBResidueInfo.GetIsHeteroAtom`.
    pub fn get_is_hetero_atom(&self) -> bool {
        expect_ok(ffi::pdb_residue_info_get_is_hetero_atom(self.raw()))
    }

    /// Mirrors `AtomPDBResidueInfo.GetSecondaryStructure`.
    pub fn get_secondary_structure(&self) -> u32 {
        expect_ok(ffi::pdb_residue_info_get_secondary_structure(self.raw()))
    }

    /// Mirrors `AtomPDBResidueInfo.GetSegmentNumber`.
    pub fn get_segment_number(&self) -> u32 {
        expect_ok(ffi::pdb_residue_info_get_segment_number(self.raw()))
    }

    /// Mirrors `AtomPDBResidueInfo.SetSerialNumber`.
    pub fn set_serial_number(&mut self, val: i32) {
        expect_ok(ffi::pdb_residue_info_set_serial_number(self.raw_mut(), val))
    }

    /// Mirrors `AtomPDBResidueInfo.SetAltLoc`.
    pub fn set_alt_loc(&mut self, val: &str) {
        expect_ok(ffi::pdb_residue_info_set_alt_loc(self.raw_mut(), val))
    }

    /// Mirrors `AtomPDBResidueInfo.SetInsertionCode`.
    pub fn set_insertion_code(&mut self, val: &str) {
        expect_ok(ffi::pdb_residue_info_set_insertion_code(
            self.raw_mut(),
            val,
        ))
    }

    /// Mirrors `AtomPDBResidueInfo.SetOccupancy`.
    pub fn set_occupancy(&mut self, val: f64) {
        expect_ok(ffi::pdb_residue_info_set_occupancy(self.raw_mut(), val))
    }

    /// Mirrors `AtomPDBResidueInfo.SetTempFactor`.
    pub fn set_temp_factor(&mut self, val: f64) {
        expect_ok(ffi::pdb_residue_info_set_temp_factor(self.raw_mut(), val))
    }

    /// Mirrors `AtomPDBResidueInfo.SetIsHeteroAtom`.
    pub fn set_is_hetero_atom(&mut self, val: bool) {
        expect_ok(ffi::pdb_residue_info_set_is_hetero_atom(
            self.raw_mut(),
            val,
        ))
    }

    /// Mirrors `AtomPDBResidueInfo.SetSecondaryStructure`.
    pub fn set_secondary_structure(&mut self, val: u32) {
        expect_ok(ffi::pdb_residue_info_set_secondary_structure(
            self.raw_mut(),
            val,
        ))
    }

    /// Mirrors `AtomPDBResidueInfo.SetSegmentNumber`.
    pub fn set_segment_number(&mut self, val: u32) {
        expect_ok(ffi::pdb_residue_info_set_segment_number(
            self.raw_mut(),
            val,
        ))
    }
}

/// Info for an unnamed atom. Mirrors `AtomPDBResidueInfo()`.
impl Default for AtomPDBResidueInfo {
    fn default() -> AtomPDBResidueInfo {
        AtomPDBResidueInfo::from_ffi(expect_ok(ffi::pdb_residue_info_new()))
    }
}

impl Clone for AtomPDBResidueInfo {
    fn clone(&self) -> AtomPDBResidueInfo {
        AtomPDBResidueInfo::from_ffi(expect_ok(ffi::pdb_residue_info_clone(self.raw())))
    }
}

impl fmt::Debug for AtomPDBResidueInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AtomPDBResidueInfo")
            .field("name", &self.get_name().ok())
            .field("residue_name", &self.get_residue_name().ok())
            .field("residue_number", &self.get_residue_number())
            .finish()
    }
}

monomer_getters!(AtomPDBResidueInfo);
monomer_setters!(AtomPDBResidueInfo);

/// `AtomMonomerType` values as C++ sees them, for a round-trip test. PDB
/// residue info takes every type.
#[doc(hidden)]
pub fn __monomer_type_round_trip_for_tests(typ: AtomMonomerType) -> AtomMonomerType {
    let mut info = AtomPDBResidueInfo::default();
    info.set_monomer_type(typ);
    info.get_monomer_type()
}
