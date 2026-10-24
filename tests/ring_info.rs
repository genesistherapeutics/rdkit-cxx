// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeSet;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn ring_queries() {
    let mol = mol_from_smiles("c1ccccc1C1CC1").unwrap();
    let info = mol.get_ring_info();
    assert_eq!(info.num_rings().unwrap(), 2);
    let rings: BTreeSet<BTreeSet<u32>> = info
        .atom_rings()
        .into_iter()
        .map(|r| r.into_iter().collect())
        .collect();
    assert_eq!(
        rings,
        BTreeSet::from([
            BTreeSet::from([0, 1, 2, 3, 4, 5]),
            BTreeSet::from([6, 7, 8])
        ])
    );
    assert_eq!(info.bond_rings().iter().map(Vec::len).sum::<usize>(), 9);
    assert_eq!(info.num_atom_rings(6).unwrap(), 1);
    assert_eq!(info.num_atom_rings(5).unwrap(), 1);
    assert_eq!(info.num_bond_rings(0).unwrap(), 1);
    assert!(info.is_atom_in_ring_of_size(6, 3).unwrap());
    assert!(info.is_atom_in_ring_of_size(0, 6).unwrap());
    let exocyclic = mol.get_bond_between_atoms(5, 6).unwrap().unwrap().get_idx();
    assert!(!info.is_bond_in_ring_of_size(exocyclic, 3).unwrap());
    assert_eq!(info.num_bond_rings(exocyclic).unwrap(), 0);
    assert_eq!(
        info.num_atom_rings(99).unwrap_err(),
        Error::Index { index: 99, len: 9 }
    );
    assert_eq!(
        info.num_bond_rings(99).unwrap_err(),
        Error::Index { index: 99, len: 10 }
    );
}

#[test]
fn uninitialized_ring_info_is_an_error_not_an_abort() {
    let mut rw = RWMol::default();
    for _ in 0..3 {
        rw.add_atom(&Atom::from_num(6));
    }
    for (begin, end) in [(0, 1), (1, 2), (2, 0)] {
        rw.add_bond_with_order(begin, end, BondType::Single)
            .unwrap();
    }
    assert!(matches!(rw.get_ring_info().num_rings(),
                     Err(Error::Runtime { message }) if message.contains("RingInfo not initialized")));
    // As in Python, the ring lists are empty rather than an error.
    assert!(rw.get_ring_info().atom_rings().is_empty());
    assert!(matches!(
        rw.get_ring_info().num_atom_rings(0),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn unsanitized_smiles_still_finds_rings() {
    // RDKit 2026.09's SMILES parser runs fast ring finding without
    // sanitization.
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let mol = mol_from_smiles_with_params("C1CC1", &params).unwrap();
    assert_eq!(mol.get_ring_info().num_rings().unwrap(), 1);
}

#[test]
fn empty_molecule_has_no_rings() {
    let mol = mol_from_smiles("").unwrap();
    assert_eq!(mol.get_ring_info().num_rings().unwrap(), 0);
    assert!(mol.get_ring_info().atom_rings().is_empty());
}
