// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
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
    // RDKit 2025.03's SMILES parser runs fast ring finding without
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

#[test]
fn memberships_and_sizes() {
    // Naphthalene's rings share the bond between atoms 4 and 9.
    let mol = mol_from_smiles("c1ccc2ccccc2c1").unwrap();
    let info = mol.get_ring_info();
    let fusion = mol.get_bond_between_atoms(3, 8).unwrap().unwrap().get_idx();
    assert_eq!(info.atom_members(3).unwrap().len(), 2);
    assert_eq!(info.atom_members(0).unwrap().len(), 1);
    assert_eq!(info.bond_members(fusion).unwrap().len(), 2);
    assert_eq!(info.atom_ring_sizes(3).unwrap(), [6, 6]);
    assert_eq!(info.bond_ring_sizes(fusion).unwrap(), [6, 6]);
    assert_eq!(info.min_atom_ring_size(0).unwrap(), 6);
    assert_eq!(info.min_bond_ring_size(fusion).unwrap(), 6);
    assert!(info.are_atoms_in_same_ring(0, 1).unwrap());
    assert!(!info.are_atoms_in_same_ring(0, 6).unwrap());
    assert!(info.are_atoms_in_same_ring_of_size(0, 1, 6).unwrap());
    assert!(!info.are_atoms_in_same_ring_of_size(0, 1, 5).unwrap());
    assert!(info.are_bonds_in_same_ring(0, 1).unwrap());
    assert!(info.are_bonds_in_same_ring_of_size(0, 1, 6).unwrap());
    assert!(info.is_ring_fused(0).unwrap() && info.are_rings_fused(0, 1).unwrap());
    assert_eq!(info.num_fused_bonds(0).unwrap(), 1);
    assert_eq!(
        info.is_ring_fused(2).unwrap_err(),
        Error::Index { index: 2, len: 2 }
    );
    assert_eq!(
        info.atom_members(10).unwrap_err(),
        Error::Index { index: 10, len: 10 }
    );
    assert_eq!(
        info.min_bond_ring_size(11).unwrap_err(),
        Error::Index { index: 11, len: 11 }
    );

    let chain = mol_from_smiles("CCC").unwrap();
    assert!(chain.get_ring_info().atom_members(0).unwrap().is_empty());
    assert_eq!(chain.get_ring_info().min_atom_ring_size(0).unwrap(), 0);
}

#[test]
fn ring_families() {
    // Parsing finds the rings, but not their families.
    let mut mol = mol_from_smiles("C1CC1").unwrap();
    let info = mol.get_ring_info();
    assert!(!info.are_ring_families_initialized());
    assert_eq!(info.num_ring_families().unwrap(), 0);
    assert!(info.atom_ring_families().is_empty());
    assert!(matches!(
        info.num_relevant_cycles(),
        Err(Error::Runtime { .. })
    ));
    find_ring_families(&mut mol).unwrap();
    let info = mol.get_ring_info();
    assert!(info.are_ring_families_initialized());
    assert_eq!(
        (
            info.num_ring_families().unwrap(),
            info.num_relevant_cycles().unwrap()
        ),
        (1, 1)
    );
    assert_eq!(info.atom_ring_families(), [vec![0, 1, 2]]);
    assert_eq!(info.bond_ring_families().len(), 1);

    let mut rw = RWMol::default();
    rw.add_atom(&Atom::from_num(6));
    let info = rw.get_ring_info();
    assert!(!info.are_ring_families_initialized());
    assert!(info.atom_ring_families().is_empty());
    assert!(matches!(
        info.num_ring_families(),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        info.num_relevant_cycles(),
        Err(Error::Runtime { .. })
    ));
}

/// A ring of `n` diamonds: `n` hub atoms, each joined to the next by two
/// parallel two-bond paths, with 2^n + n relevant cycles. Its rings are found
/// without sanitization, whose symmetrized SSSR enumerates the relevant
/// cycles.
fn diamonds(n: u32) -> RWMol {
    let mut rw = RWMol::default();
    let hubs: Vec<u32> = (0..n).map(|_| rw.add_atom(&Atom::from_num(6))).collect();
    for (i, &hub) in hubs.iter().enumerate() {
        let next = hubs[(i + 1) % hubs.len()];
        for _ in 0..2 {
            let middle = rw.add_atom(&Atom::from_num(6));
            rw.add_bond_with_order(hub, middle, BondType::Single)
                .unwrap();
            rw.add_bond_with_order(middle, next, BondType::Single)
                .unwrap();
        }
    }
    rw.update_property_cache().unwrap();
    fast_find_rings(&mut rw).unwrap();
    rw
}

#[test]
fn relevant_cycles_beyond_u32_are_an_error() {
    let mut rw = diamonds(10);
    find_ring_families(&mut rw).unwrap();
    assert_eq!(rw.get_ring_info().num_relevant_cycles().unwrap(), 1034);
    // 2^33 + 33, which RDKit 2025.03's numRelevantCycles converts to
    // unsigned int unchecked: Python returns 33.
    let mut rw = diamonds(33);
    find_ring_families(&mut rw).unwrap();
    assert!(matches!(
        rw.get_ring_info().num_relevant_cycles(),
        Err(Error::Value { .. })
    ));
}

#[test]
fn add_ring() {
    let mut rw = RWMol::default();
    for _ in 0..3 {
        rw.add_atom(&Atom::from_num(6));
    }
    for (begin, end) in [(0, 1), (1, 2), (2, 0)] {
        rw.add_bond_with_order(begin, end, BondType::Single)
            .unwrap();
    }
    let mut info = rw.get_ring_info_mut();
    assert!(matches!(
        info.add_ring(&[0, 1, 2], &[0, 1]),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        info.add_ring(&[0, 1, 5], &[0, 1, 2]).unwrap_err(),
        Error::Index { index: 5, len: 3 }
    );
    info.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
    assert_eq!(info.num_rings().unwrap(), 1);
    assert!(rw.get_ring_info().is_atom_in_ring_of_size(2, 3).unwrap());
}
