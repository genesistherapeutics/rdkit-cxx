// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::{BondType, ChiralType, HybridizationType, mol_from_smiles};

#[test]
fn atom_getters() {
    let mol = mol_from_smiles("OCC").unwrap();
    let symbols: Vec<_> = mol.get_atoms().map(|a| a.get_symbol().unwrap()).collect();
    assert_eq!(symbols, ["O", "C", "C"]);
    let degrees: Vec<_> = mol.get_atoms().map(|a| a.get_degree()).collect();
    assert_eq!(degrees, [1, 2, 1]);
    let hs: Vec<_> = mol
        .get_atoms()
        .map(|a| a.get_total_num_hs().unwrap())
        .collect();
    assert_eq!(hs, [1, 2, 3]);
    let atom = mol.get_atom_with_idx(1).unwrap();
    let neighbors: Vec<_> = atom.get_neighbors().iter().map(|a| a.get_idx()).collect();
    assert_eq!(neighbors, [0, 2]);
    let bonds: Vec<_> = atom.get_bonds().iter().map(|b| b.get_idx()).collect();
    assert_eq!(bonds, [0, 1]);
    assert_eq!(atom.get_total_degree().unwrap(), 4);
    assert_eq!(atom.get_total_valence().unwrap(), 4);
    assert_eq!(atom.get_num_implicit_hs().unwrap(), 2);
    assert!(!atom.is_in_ring());
    assert!(std::ptr::eq(atom.get_owning_mol(), &mol));
}

#[test]
fn ring_membership_is_found_on_demand() {
    let mol = mol_from_smiles("C1CC1C").unwrap();
    let in_ring: Vec<_> = mol.get_atoms().map(|a| a.is_in_ring()).collect();
    assert_eq!(in_ring, [true, true, true, false]);
}

#[test]
fn enums() {
    let mol = mol_from_smiles("C=CO").unwrap();
    let types: Vec<_> = mol.get_bonds().map(|b| b.get_bond_type()).collect();
    assert_eq!(types, [BondType::Double, BondType::Single]);
    let hybridizations: Vec<_> = mol.get_atoms().map(|a| a.get_hybridization()).collect();
    assert_eq!(hybridizations, [HybridizationType::Sp2; 3]);
    let chiral = mol_from_smiles("F[C@H](Cl)Br").unwrap();
    assert_eq!(
        chiral.get_atom_with_idx(1).unwrap().get_chiral_tag(),
        ChiralType::ChiTetrahedralCcw
    );
    let benzene = mol_from_smiles("c1ccccc1").unwrap();
    let bond = benzene.get_bond_with_idx(0).unwrap();
    assert_eq!(bond.get_bond_type(), BondType::Aromatic);
    assert_eq!(bond.get_bond_type_as_double().unwrap(), 1.5);
    assert!(bond.get_is_aromatic() && bond.get_is_conjugated());
}

#[test]
fn bonds_and_lookups() {
    let mol = mol_from_smiles("C=CO").unwrap();
    let ends: Vec<_> = mol
        .get_bonds()
        .map(|b| (b.get_begin_atom_idx(), b.get_end_atom_idx()))
        .collect();
    assert_eq!(ends, [(0, 1), (1, 2)]);
    let bond = mol.get_bond_between_atoms(0, 1).unwrap().unwrap();
    assert_eq!(bond.get_idx(), 0);
    assert_eq!(bond.get_other_atom_idx(1).unwrap(), 0);
    assert!(matches!(
        bond.get_other_atom_idx(2),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(
        (
            bond.get_begin_atom().get_idx(),
            bond.get_end_atom().get_idx()
        ),
        (0, 1)
    );
    assert!(mol.get_bond_between_atoms(0, 2).unwrap().is_none());
    assert_eq!(
        mol.get_bond_between_atoms(0, 99).unwrap_err(),
        Error::Index { index: 99, len: 3 }
    );
    assert_eq!(
        mol.get_atom_with_idx(99).unwrap_err(),
        Error::Index { index: 99, len: 3 }
    );
    assert_eq!(
        mol.get_bond_with_idx(99).unwrap_err(),
        Error::Index { index: 99, len: 2 }
    );
}

#[test]
fn iterators() {
    let mol = mol_from_smiles("OCC").unwrap();
    assert_eq!(mol.get_atoms().len(), 3);
    assert_eq!(
        mol.get_atoms()
            .rev()
            .map(|a| a.get_idx())
            .collect::<Vec<_>>(),
        [2, 1, 0]
    );
    assert_eq!(mol.get_bonds().len(), 2);
    let empty = mol_from_smiles("").unwrap();
    assert_eq!(
        (empty.get_atoms().count(), empty.get_bonds().count()),
        (0, 0)
    );
}

#[test]
fn views_debug_as_their_index() {
    let mol = mol_from_smiles("OCC").unwrap();
    assert_eq!(
        format!("{:?}", mol.get_atom_with_idx(2).unwrap()),
        "AtomRef { idx: 2 }"
    );
    assert_eq!(
        format!("{:?}", mol.get_bond_with_idx(1).unwrap()),
        "BondRef { idx: 1 }"
    );
}
