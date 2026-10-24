// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Result;
use crate::error::{expect_ok, to_string};
use crate::ffi::periodic_table::ffi;

/// RDKit's periodic table. Mirrors `rdchem.PeriodicTable`.
///
/// It is a handle to a process-wide table that never changes, so it is `Copy`,
/// `Send`, and `Sync`. Atomic numbers outside the table are `Error::Index`.
/// Methods that take an element symbol end in `_str` (CONVENTIONS.md A5);
/// an unknown symbol is `Error::Runtime`.
#[derive(Clone, Copy, Debug, Default)]
pub struct PeriodicTable {
    _private: (),
}

/// Mirrors `rdchem.GetPeriodicTable`.
pub fn get_periodic_table() -> PeriodicTable {
    PeriodicTable::default()
}

/// Defines a lookup by atomic number and its `_str` form by element symbol.
macro_rules! lookups {
    ($($(#[$meta:meta])* $name:ident, $name_str:ident, $ffi:ident -> $out:ty;)*) => {
        impl PeriodicTable {
            $(
                $(#[$meta])*
                pub fn $name(&self, atomic_number: u32) -> Result<$out> {
                    Ok(ffi::$ffi(atomic_number)?)
                }

                #[doc = concat!("As [`", stringify!($name), "`](Self::", stringify!($name), "), by element symbol.")]
                pub fn $name_str(&self, element_symbol: &str) -> Result<$out> {
                    self.$name(self.get_atomic_number(element_symbol)?)
                }
            )*
        }
    };
}

lookups! {
    /// The average atomic weight. Mirrors `PeriodicTable.GetAtomicWeight`.
    get_atomic_weight, get_atomic_weight_str, periodic_table_get_atomic_weight -> f64;
    /// The radius for bond length estimates. Mirrors `PeriodicTable.GetRb0`.
    get_rb0, get_rb0_str, periodic_table_get_rb0 -> f64;
    /// The covalent radius. Mirrors `PeriodicTable.GetRcovalent`.
    get_rcovalent, get_rcovalent_str, periodic_table_get_rcovalent -> f64;
    /// The van der Waals radius. Mirrors `PeriodicTable.GetRvdw`.
    get_rvdw, get_rvdw_str, periodic_table_get_rvdw -> f64;
    /// The default valence, or -1 for any. Mirrors
    /// `PeriodicTable.GetDefaultValence`.
    get_default_valence, get_default_valence_str, periodic_table_get_default_valence -> i32;
    /// The allowed valences; -1 means any. Mirrors
    /// `PeriodicTable.GetValenceList`.
    get_valence_list, get_valence_list_str, periodic_table_get_valence_list -> Vec<i32>;
    /// The number of outer-shell electrons. Mirrors
    /// `PeriodicTable.GetNOuterElecs`.
    get_n_outer_elecs, get_n_outer_elecs_str, periodic_table_get_n_outer_elecs -> i32;
    /// Mirrors `PeriodicTable.GetMostCommonIsotope`.
    get_most_common_isotope, get_most_common_isotope_str, periodic_table_get_most_common_isotope -> i32;
    /// Mirrors `PeriodicTable.GetMostCommonIsotopeMass`.
    get_most_common_isotope_mass, get_most_common_isotope_mass_str, periodic_table_get_most_common_isotope_mass -> f64;
    /// The row of the periodic table. Mirrors `PeriodicTable.GetRow`.
    get_row, get_row_str, periodic_table_get_row -> u32;
}

impl PeriodicTable {
    /// Mirrors `PeriodicTable.GetAtomicNumber`.
    pub fn get_atomic_number(&self, element_symbol: &str) -> Result<u32> {
        Ok(ffi::periodic_table_get_atomic_number(element_symbol)?)
    }

    /// Mirrors `PeriodicTable.GetMaxAtomicNumber`.
    pub fn get_max_atomic_number(&self) -> u32 {
        expect_ok(ffi::periodic_table_get_max_atomic_number())
    }

    /// Mirrors `PeriodicTable.GetElementSymbol`.
    pub fn get_element_symbol(&self, atomic_number: u32) -> Result<String> {
        to_string(ffi::periodic_table_get_element_symbol(atomic_number)?)
    }

    /// Mirrors `PeriodicTable.GetElementName`.
    pub fn get_element_name(&self, atomic_number: u32) -> Result<String> {
        to_string(ffi::periodic_table_get_element_name(atomic_number)?)
    }

    /// An isotope's mass, or 0 for an unknown isotope. Mirrors
    /// `PeriodicTable.GetMassForIsotope`.
    pub fn get_mass_for_isotope(&self, atomic_number: u32, isotope: u32) -> Result<f64> {
        Ok(ffi::periodic_table_get_mass_for_isotope(
            atomic_number,
            isotope,
        )?)
    }

    /// As [`get_mass_for_isotope`](Self::get_mass_for_isotope), by element
    /// symbol.
    pub fn get_mass_for_isotope_str(&self, element_symbol: &str, isotope: u32) -> Result<f64> {
        self.get_mass_for_isotope(self.get_atomic_number(element_symbol)?, isotope)
    }

    /// An isotope's natural abundance, in percent. Mirrors
    /// `PeriodicTable.GetAbundanceForIsotope`.
    pub fn get_abundance_for_isotope(&self, atomic_number: u32, isotope: u32) -> Result<f64> {
        Ok(ffi::periodic_table_get_abundance_for_isotope(
            atomic_number,
            isotope,
        )?)
    }

    /// As [`get_abundance_for_isotope`](Self::get_abundance_for_isotope), by
    /// element symbol.
    pub fn get_abundance_for_isotope_str(&self, element_symbol: &str, isotope: u32) -> Result<f64> {
        self.get_abundance_for_isotope(self.get_atomic_number(element_symbol)?, isotope)
    }
}
