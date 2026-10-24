// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/periodic_table.h");

        fn periodic_table_get_atomic_number(element_symbol: &str) -> Result<u32>;
        fn periodic_table_get_max_atomic_number() -> Result<u32>;
        fn periodic_table_get_element_symbol(atomic_number: u32) -> Result<UniquePtr<CxxString>>;
        fn periodic_table_get_element_name(atomic_number: u32) -> Result<UniquePtr<CxxString>>;
        fn periodic_table_get_atomic_weight(atomic_number: u32) -> Result<f64>;
        fn periodic_table_get_rb0(atomic_number: u32) -> Result<f64>;
        fn periodic_table_get_rcovalent(atomic_number: u32) -> Result<f64>;
        fn periodic_table_get_rvdw(atomic_number: u32) -> Result<f64>;
        fn periodic_table_get_default_valence(atomic_number: u32) -> Result<i32>;
        fn periodic_table_get_valence_list(atomic_number: u32) -> Result<Vec<i32>>;
        fn periodic_table_get_n_outer_elecs(atomic_number: u32) -> Result<i32>;
        fn periodic_table_get_most_common_isotope(atomic_number: u32) -> Result<i32>;
        fn periodic_table_get_most_common_isotope_mass(atomic_number: u32) -> Result<f64>;
        fn periodic_table_get_row(atomic_number: u32) -> Result<u32>;
        fn periodic_table_get_mass_for_isotope(atomic_number: u32, isotope: u32) -> Result<f64>;
        fn periodic_table_get_abundance_for_isotope(
            atomic_number: u32,
            isotope: u32,
        ) -> Result<f64>;
    }
}
