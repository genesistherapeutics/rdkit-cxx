// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Result;
use crate::chem::rdchem::{Mol, RWMol};
use crate::error::expect_ok;
use crate::ffi::mol::ffi;

/// What a pickle includes, as bits. Mirrors `rdchem.PropertyPickleOptions`
/// (CONVENTIONS.md F1); the constants keep Python's spelling.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct PropertyPickleOptions(u32);

macro_rules! pickle_options {
    ($($name:ident = $value:expr),* $(,)?) => {
        #[allow(non_upper_case_globals)]
        impl PropertyPickleOptions {
            $(
                #[doc = concat!("Mirrors `PropertyPickleOptions.", stringify!($name), "`.")]
                pub const $name: PropertyPickleOptions = PropertyPickleOptions($value);
            )*

            const NAMES: &[(&str, PropertyPickleOptions)] =
                &[$((stringify!($name), PropertyPickleOptions::$name)),*];
        }
    };
}

// In the order property_pickle_option_values returns them; a test checks the
// values against MolPickler.h.
pickle_options! {
    NoProps = 0x0,
    MolProps = 0x1,
    AtomProps = 0x2,
    BondProps = 0x4,
    QueryAtomData = 0x2,
    PrivateProps = 0x10,
    ComputedProps = 0x20,
    AllProps = 0xFFFF,
    CoordsAsDouble = 0x10000,
    NoConformers = 0x20000,
}

crate::flags::flag_ops!(PropertyPickleOptions);

/// What pickles include when not told. Mirrors
/// `rdchem.GetDefaultPickleProperties`.
pub fn get_default_pickle_properties() -> PropertyPickleOptions {
    PropertyPickleOptions(expect_ok(ffi::get_default_pickle_properties()))
}

/// Sets what pickles include when not told, for the whole process. Mirrors
/// `rdchem.SetDefaultPickleProperties`.
pub fn set_default_pickle_properties(arg1: PropertyPickleOptions) {
    expect_ok(ffi::set_default_pickle_properties(arg1.bits()))
}

#[doc(hidden)]
pub fn __property_pickle_option_values_for_tests() -> (Vec<u32>, Vec<u32>) {
    let rust = PropertyPickleOptions::NAMES
        .iter()
        .map(|(_, options)| options.bits())
        .collect();
    (rust, expect_ok(ffi::property_pickle_option_values()))
}

impl Mol {
    /// A molecule from a pickle (`to_binary`). Mirrors `Mol(pklString)`.
    pub fn from_pkl_string(pkl_string: &[u8]) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::mol_from_pkl_string(
            pkl_string, false, 0,
        )?))
    }

    /// Reads only what `property_flags` selects. Mirrors
    /// `Mol(pklString, propertyFlags)`.
    pub fn from_pkl_string_with_property_flags(
        pkl_string: &[u8],
        property_flags: PropertyPickleOptions,
    ) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::mol_from_pkl_string(
            pkl_string,
            true,
            property_flags.bits(),
        )?))
    }

    /// The molecule as a pickle, with the default properties. Mirrors
    /// `Mol.ToBinary`.
    pub fn to_binary(&self) -> Vec<u8> {
        expect_ok(ffi::mol_to_binary(self.as_ffi(), false, 0))
    }

    /// Mirrors `Mol.ToBinary(propertyFlags)`.
    pub fn to_binary_with_property_flags(&self, property_flags: PropertyPickleOptions) -> Vec<u8> {
        expect_ok(ffi::mol_to_binary(
            self.as_ffi(),
            true,
            property_flags.bits(),
        ))
    }
}

impl RWMol {
    /// A molecule from a pickle. Mirrors `RWMol(pklString)`.
    pub fn from_pkl_string(pkl_string: &[u8]) -> Result<RWMol> {
        Ok(Mol::from_pkl_string(pkl_string)?.into())
    }

    /// Mirrors `RWMol(pklString, propertyFlags)`.
    pub fn from_pkl_string_with_property_flags(
        pkl_string: &[u8],
        property_flags: PropertyPickleOptions,
    ) -> Result<RWMol> {
        Ok(Mol::from_pkl_string_with_property_flags(pkl_string, property_flags)?.into())
    }
}
