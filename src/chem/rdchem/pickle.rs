// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Result;
use crate::chem::rdchem::{Mol, RWMol};
use crate::error::{expect_ok, picklable};
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

#[cfg(test)]
pub(crate) fn property_pickle_option_values_for_tests() -> (Vec<u32>, Vec<u32>) {
    let rust = PropertyPickleOptions::NAMES
        .iter()
        .map(|(_, options)| options.bits())
        .collect();
    (rust, expect_ok(ffi::property_pickle_option_values()))
}

impl Mol {
    /// A molecule from a pickle (`to_binary`). Mirrors `Mol(pklString)`.
    ///
    /// # Safety
    ///
    /// `pkl_string` must be a molecule pickle RDKit wrote (this crate's
    /// `to_binary` or RDKit's), unmodified. RDKit's pickle reader trusts the
    /// counts, indices, enum values, and query trees it reads, and other
    /// bytes may make it, or code that later reads the molecule, read out of
    /// bounds (CONVENTIONS.md U2).
    pub unsafe fn from_pkl_string(pkl_string: &[u8]) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::mol_from_pkl_string(
            pkl_string, false, 0,
        )?))
    }

    /// Reads only what `property_flags` selects. Mirrors
    /// `Mol(pklString, propertyFlags)`.
    ///
    /// # Safety
    ///
    /// As for `from_pkl_string`.
    pub unsafe fn from_pkl_string_with_property_flags(
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
    ///
    /// # Panics
    ///
    /// For a molecule RDKit would pickle truncated, which Python's writes
    /// and cannot read back: a query of more than 255 children, or, in a
    /// molecule of at most 255 atoms, a substance group or stereo group list
    /// of more than 255 entries, more than 255 stereo groups, or a group
    /// naming a bond past 255; or a molecule, conformer, atom, bond, or
    /// substance group of more than 65,535 properties, which RDKit counts in
    /// two bytes.
    pub fn to_binary(&self) -> Vec<u8> {
        picklable(ffi::mol_to_binary(self.as_ffi(), false, 0))
    }

    /// Mirrors `Mol.ToBinary(propertyFlags)`.
    ///
    /// # Panics
    ///
    /// As `to_binary`.
    pub fn to_binary_with_property_flags(&self, property_flags: PropertyPickleOptions) -> Vec<u8> {
        picklable(ffi::mol_to_binary(
            self.as_ffi(),
            true,
            property_flags.bits(),
        ))
    }
}

impl RWMol {
    /// A molecule from a pickle. Mirrors `RWMol(pklString)`.
    ///
    /// # Safety
    ///
    /// As for `Mol::from_pkl_string`.
    pub unsafe fn from_pkl_string(pkl_string: &[u8]) -> Result<RWMol> {
        // SAFETY: the caller's contract is `Mol::from_pkl_string`'s.
        Ok(unsafe { Mol::from_pkl_string(pkl_string) }?.into())
    }

    /// Mirrors `RWMol(pklString, propertyFlags)`.
    ///
    /// # Safety
    ///
    /// As for `Mol::from_pkl_string`.
    pub unsafe fn from_pkl_string_with_property_flags(
        pkl_string: &[u8],
        property_flags: PropertyPickleOptions,
    ) -> Result<RWMol> {
        // SAFETY: the caller's contract is `Mol::from_pkl_string`'s.
        Ok(unsafe { Mol::from_pkl_string_with_property_flags(pkl_string, property_flags) }?.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn option_values_match_rdkit() {
        let (rust, cxx) = property_pickle_option_values_for_tests();
        assert_eq!(rust, cxx);
        assert_eq!(
            PropertyPickleOptions::QueryAtomData,
            PropertyPickleOptions::AtomProps
        );
        let combined = PropertyPickleOptions::MolProps | PropertyPickleOptions::AtomProps;
        assert_eq!(combined.bits(), 0x3);
        assert_eq!(
            format!("{:?}", PropertyPickleOptions::AllProps),
            "PropertyPickleOptions::AllProps"
        );
    }
}
