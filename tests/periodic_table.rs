// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::get_periodic_table;

#[test]
fn lookups() {
    let table = get_periodic_table();
    assert_eq!(table.get_atomic_number("C").unwrap(), 6);
    assert_eq!(table.get_element_symbol(8).unwrap(), "O");
    assert_eq!(table.get_element_name(8).unwrap(), "Oxygen");
    assert!((table.get_atomic_weight(6).unwrap() - 12.011).abs() < 1e-3);
    assert_eq!(
        table.get_atomic_weight_str("C").unwrap(),
        table.get_atomic_weight(6).unwrap()
    );
    assert_eq!(table.get_default_valence(6).unwrap(), 4);
    assert_eq!(table.get_default_valence_str("N").unwrap(), 3);
    assert_eq!(table.get_valence_list(8).unwrap(), [2]);
    assert_eq!(table.get_n_outer_elecs(6).unwrap(), 4);
    assert_eq!(table.get_most_common_isotope(6).unwrap(), 12);
    assert!((table.get_most_common_isotope_mass(6).unwrap() - 12.0).abs() < 1e-6);
    assert_eq!(table.get_row(11).unwrap(), 3);
    assert_eq!(table.get_row_str("Na").unwrap(), 3);
    assert!((table.get_mass_for_isotope(6, 13).unwrap() - 13.003).abs() < 1e-3);
    assert!(table.get_abundance_for_isotope_str("C", 12).unwrap() > 98.0);
    assert!(table.get_rvdw(6).unwrap() > table.get_rcovalent(6).unwrap());
    assert!(table.get_rb0(6).unwrap() > 0.0);
    assert!(table.get_max_atomic_number() >= 118);
}

#[test]
fn bad_inputs() {
    let table = get_periodic_table();
    let max = table.get_max_atomic_number();
    assert_eq!(
        table.get_atomic_weight(500).unwrap_err(),
        Error::Index {
            index: 500,
            len: max as usize + 1
        }
    );
    assert!(matches!(
        table.get_atomic_number("Xx"),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        table.get_row_str("Xx"),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn shared_across_threads() {
    fn assert_send_sync<T: Send + Sync + Copy>() {}
    assert_send_sync::<rdkit_cxx::chem::PeriodicTable>();
    let table = get_periodic_table();
    let weight = std::thread::spawn(move || table.get_atomic_weight(1).unwrap())
        .join()
        .unwrap();
    assert!((weight - 1.008).abs() < 1e-3);
}
