// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdChemicalFeatures`: chemical features that belong to
//! no molecule.

use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::error::{expect_ok, to_string};
use crate::ffi::chemical_features::ffi;
use crate::geometry::rd_geometry::Point3D;

/// A chemical feature with a family, a type, an id, and a position, not tied
/// to a molecule. Mirrors `rdChemicalFeatures.FreeChemicalFeature`.
///
/// `Default` mirrors `FreeChemicalFeature()`: empty family and type, id -1,
/// and the origin. `Clone` copies it.
pub struct FreeChemicalFeature {
    /// Never null.
    ptr: UniquePtr<ffi::FreeChemicalFeature>,
}

impl Default for FreeChemicalFeature {
    fn default() -> FreeChemicalFeature {
        FreeChemicalFeature::from_ffi(expect_ok(ffi::free_chemical_feature_new()))
    }
}

impl Clone for FreeChemicalFeature {
    fn clone(&self) -> FreeChemicalFeature {
        FreeChemicalFeature::from_ffi(expect_ok(ffi::free_chemical_feature_clone(&self.ptr)))
    }
}

impl fmt::Debug for FreeChemicalFeature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FreeChemicalFeature")
            .field("id", &self.get_id())
            .field("family", &self.get_family())
            .field("type", &self.get_type())
            .field("pos", &self.get_pos())
            .finish()
    }
}

// Every string is UTF-8: Rust gave it, or from_pickle checked it.
fn text(value: std::result::Result<UniquePtr<cxx::CxxString>, cxx::Exception>) -> String {
    to_string(expect_ok(value)).expect("a FreeChemicalFeature string is not UTF-8")
}

impl FreeChemicalFeature {
    fn from_ffi(ptr: UniquePtr<ffi::FreeChemicalFeature>) -> FreeChemicalFeature {
        assert!(!ptr.is_null(), "RDKit returned a null feature");
        FreeChemicalFeature { ptr }
    }

    /// A feature with a family, a type, and a position, and id -1. Mirrors
    /// `FreeChemicalFeature(family, type, loc)`.
    pub fn new(family: &str, type_: &str, loc: Point3D) -> FreeChemicalFeature {
        FreeChemicalFeature::new_with_id(family, type_, loc, -1)
    }

    /// Mirrors `FreeChemicalFeature(family, type, loc, id)`.
    pub fn new_with_id(family: &str, type_: &str, loc: Point3D, id: i32) -> FreeChemicalFeature {
        FreeChemicalFeature::from_ffi(expect_ok(ffi::free_chemical_feature_from_parts(
            family, type_, &loc, id,
        )))
    }

    /// A feature with a family and a position, an empty type, and id -1.
    /// Mirrors `FreeChemicalFeature(family, loc)`.
    pub fn from_family(family: &str, loc: Point3D) -> FreeChemicalFeature {
        FreeChemicalFeature::new(family, "", loc)
    }

    /// A feature from a pickle that `to_pickle` (or Python's
    /// `__getstate__`) wrote. Each string is read up to its first NUL.
    /// `Error::Value` for a pickle RDKit would read past the end of (an
    /// unknown version, a string without its NUL, or bytes missing) or with a
    /// string that is not UTF-8. Mirrors `FreeChemicalFeature(pickle)`.
    pub fn from_pickle(pickle: &[u8]) -> Result<FreeChemicalFeature> {
        let feat = FreeChemicalFeature::from_ffi(ffi::free_chemical_feature_from_pickle(pickle)?);
        to_string(expect_ok(ffi::free_chemical_feature_get_family(&feat.ptr)))?;
        to_string(expect_ok(ffi::free_chemical_feature_get_type(&feat.ptr)))?;
        Ok(feat)
    }

    /// The feature as a pickle `from_pickle` reads: what Python's
    /// `__getstate__` returns, without the tuple.
    pub fn to_pickle(&self) -> Vec<u8> {
        expect_ok(ffi::free_chemical_feature_to_pickle(&self.ptr))
    }

    /// The feature's id. Mirrors `FreeChemicalFeature.GetId`.
    pub fn get_id(&self) -> i32 {
        expect_ok(ffi::free_chemical_feature_get_id(&self.ptr))
    }

    /// Mirrors `FreeChemicalFeature.SetId`.
    pub fn set_id(&mut self, id: i32) {
        expect_ok(ffi::free_chemical_feature_set_id(self.ptr.pin_mut(), id));
    }

    /// The feature's family (donor, acceptor, and so on). Mirrors
    /// `FreeChemicalFeature.GetFamily`.
    pub fn get_family(&self) -> String {
        text(ffi::free_chemical_feature_get_family(&self.ptr))
    }

    /// Mirrors `FreeChemicalFeature.SetFamily`.
    pub fn set_family(&mut self, family: &str) {
        expect_ok(ffi::free_chemical_feature_set_family(
            self.ptr.pin_mut(),
            family,
        ));
    }

    /// The feature's specific type. Mirrors `FreeChemicalFeature.GetType`.
    pub fn get_type(&self) -> String {
        text(ffi::free_chemical_feature_get_type(&self.ptr))
    }

    /// Mirrors `FreeChemicalFeature.SetType`.
    pub fn set_type(&mut self, type_: &str) {
        expect_ok(ffi::free_chemical_feature_set_type(
            self.ptr.pin_mut(),
            type_,
        ));
    }

    /// The feature's position. Mirrors `FreeChemicalFeature.GetPos`.
    pub fn get_pos(&self) -> Point3D {
        expect_ok(ffi::free_chemical_feature_get_pos(&self.ptr))
    }

    /// Mirrors `FreeChemicalFeature.SetPos`.
    pub fn set_pos(&mut self, loc: Point3D) {
        expect_ok(ffi::free_chemical_feature_set_pos(self.ptr.pin_mut(), &loc));
    }
}
