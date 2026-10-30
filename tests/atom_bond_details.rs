// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn atom_details() {
    let mol = mol_from_smiles("[CH3:7][O]").unwrap();
    let carbon = mol.get_atom_with_idx(0).unwrap();
    assert_eq!(carbon.get_atom_map_num(), 7);
    assert!((carbon.get_mass().unwrap() - 12.011).abs() < 1e-3);
    assert!(carbon.get_no_implicit());
    assert!(carbon.has_owning_mol() && !carbon.has_query());
    assert!(!carbon.needs_update_property_cache());
    assert_eq!(carbon.get_explicit_valence().unwrap(), 4);
    assert_eq!(carbon.get_implicit_valence().unwrap(), 0);
    assert!(!carbon.has_valence_violation().unwrap());
    assert_eq!(carbon.get_smarts().unwrap(), "[CH3:7]");
    let oxygen = mol.get_atom_with_idx(1).unwrap();
    assert_eq!(oxygen.get_num_radical_electrons(), 1);
    assert_eq!(
        GetSmartsParams::default(),
        GetSmartsParams {
            do_kekule: false,
            all_hs_explicit: false,
            isomeric_smiles: true,
        }
    );
    let iso = mol_from_smiles("[13CH4]").unwrap();
    let all_hs = GetSmartsParams {
        isomeric_smiles: false,
        ..Default::default()
    };
    assert_eq!(
        iso.get_atom_with_idx(0)
            .unwrap()
            .get_smarts_with_params(&all_hs)
            .unwrap(),
        "C"
    );
    assert!(matches!(
        Atom::from_num(500).get_mass(),
        Err(Error::Runtime { .. })
    ));
    assert!(!Atom::from_num(6).has_owning_mol());
    assert!(matches!(
        Atom::from_num(6).get_explicit_valence(),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn query_atoms_write_smarts() {
    let query = mol_from_smarts("[#6;R]").unwrap();
    assert!(query.has_query());
    let atom = query.get_atom_with_idx(0).unwrap();
    assert!(atom.has_query());
    assert_eq!(atom.get_smarts().unwrap(), "[#6&R]");
    assert!(!mol_from_smiles("C").unwrap().has_query());
}

#[test]
fn atom_setters() {
    let mut mol = mol_from_smiles("CC").unwrap();
    let mut atom = mol.get_atom_with_idx_mut(0).unwrap();
    atom.set_atom_map_num(3);
    assert_eq!(atom.get_atom_map_num(), 3);
    atom.set_atom_map_num(0);
    assert_eq!(atom.get_atom_map_num(), 0);
    atom.set_atom_map_num(5000);
    assert!(matches!(
        atom.set_atom_map_num_with_strict(5000, true),
        Err(Error::Runtime { .. })
    ));
    atom.set_atom_map_num_with_strict(12, true).unwrap();
    atom.set_num_radical_electrons(1);
    assert_eq!(atom.get_num_radical_electrons(), 1);
    assert!(!atom.invert_chirality());
    let mut chiral = mol_from_smiles("F[C@H](Cl)Br").unwrap();
    assert!(chiral.get_atom_with_idx_mut(1).unwrap().invert_chirality());
    assert_eq!(mol_to_smiles(&chiral).unwrap(), "F[C@@H](Cl)Br");
    let mut owned = Atom::from_num(8);
    owned.set_num_radical_electrons(2);
    assert_eq!(owned.get_num_radical_electrons(), 2);
}

#[test]
fn ring_sizes() {
    let mol = mol_from_smiles("C1CC1C1CCCC1").unwrap();
    let atom = mol.get_atom_with_idx(0).unwrap();
    assert!(atom.is_in_ring_size(3) && !atom.is_in_ring_size(5));
    let bridge = mol.get_bond_between_atoms(2, 3).unwrap().unwrap();
    assert!(!bridge.is_in_ring());
    let ring_bond = mol.get_bond_between_atoms(0, 1).unwrap().unwrap();
    assert!(ring_bond.is_in_ring() && ring_bond.is_in_ring_size(3));
    assert!(!ring_bond.is_in_ring_size(5));
}

#[test]
fn bond_details() {
    let mol = mol_from_smiles("C/C=C/C").unwrap();
    let double = mol.get_bond_with_idx(1).unwrap();
    assert_eq!(double.get_stereo(), BondStereo::Stereoe);
    assert_eq!(double.get_stereo_atoms(), [0, 3]);
    assert_eq!(
        mol.get_bond_with_idx(0).unwrap().get_bond_dir(),
        BondDir::Endupright
    );
    assert!(double.has_owning_mol() && !double.has_query());
    assert_eq!(double.get_smarts().unwrap(), "=");
    let single = mol.get_bond_with_idx(2).unwrap();
    assert_eq!(
        single.get_smarts_with_all_bonds_explicit(true).unwrap(),
        "/"
    );
    let begin = double.get_begin_atom();
    assert_eq!(double.get_other_atom(&begin).unwrap().get_idx(), 2);
    assert!(matches!(
        double.get_other_atom(&mol.get_atom_with_idx(0).unwrap()),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(double.get_valence_contrib(&begin), 2.0);
    assert_eq!(
        double.get_valence_contrib(&mol.get_atom_with_idx(0).unwrap()),
        0.0
    );
    assert_eq!(
        mol_from_smarts("C-,=C")
            .unwrap()
            .get_bond_with_idx(0)
            .unwrap()
            .get_smarts()
            .unwrap(),
        "-,="
    );
}

#[test]
fn bond_setters() {
    let mut mol = mol_from_smiles("CC=CC").unwrap();
    let mut double = mol.get_bond_with_idx_mut(1).unwrap();
    double.set_stereo_atoms(0, 3).unwrap();
    double.set_stereo(BondStereo::Stereocis);
    assert_eq!(double.get_stereo_atoms(), [0, 3]);
    assert!(matches!(
        double.set_stereo_atoms(3, 0),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(
        double.set_stereo_atoms(0, 9).unwrap_err(),
        Error::Index { index: 9, len: 4 }
    );
    assert!(!double.invert_chirality());
    double.set_stereo(BondStereo::Stereoatropcw);
    assert!(double.invert_chirality());
    assert_eq!(double.get_stereo(), BondStereo::Stereoatropccw);
    double.set_stereo(BondStereo::Stereocis);
    mol.get_bond_with_idx_mut(0)
        .unwrap()
        .set_bond_dir(BondDir::Beginwedge);
    assert_eq!(
        mol.get_bond_with_idx(0).unwrap().get_bond_dir(),
        BondDir::Beginwedge
    );
    assert_eq!(mol_to_smiles(&mol).unwrap(), "C/C=C\\C");
}

#[test]
fn mol_details() {
    let mut mol = mol_from_smiles("c1ccccc1CO").unwrap();
    assert_eq!(
        mol.get_aromatic_atoms()
            .iter()
            .map(|a| a.get_idx())
            .collect::<Vec<_>>(),
        [0, 1, 2, 3, 4, 5]
    );
    assert!(!mol.needs_update_property_cache());
    mol.update_property_cache().unwrap();
    mol.debug_with_use_stdout(false);
    assert!(mol_from_smiles("").unwrap().get_aromatic_atoms().is_empty());
}

#[test]
fn rwmol_details() {
    let mut rw = RWMol::from(mol_from_smiles("CC").unwrap());
    rw.insert_mol(&mol_from_smiles("O").unwrap());
    assert_eq!(mol_to_smiles(&rw).unwrap(), "CC.O");

    let donor = mol_from_smiles("C=C").unwrap();
    let mut target = RWMol::from(mol_from_smiles("CC").unwrap());
    target
        .replace_bond(0, &donor.get_bond_with_idx(0).unwrap())
        .unwrap();
    target.update_property_cache().unwrap();
    assert_eq!(mol_to_smiles(&target).unwrap(), "C=C");
    assert_eq!(
        target
            .replace_bond(3, &donor.get_bond_with_idx(0).unwrap())
            .unwrap_err(),
        Error::Index { index: 3, len: 1 }
    );
    let params = ReplaceBondParams {
        preserve_props: true,
        ..Default::default()
    };
    target
        .replace_bond_with_params(0, &donor.get_bond_with_idx(0).unwrap(), &params)
        .unwrap();
    assert_eq!(
        ReplaceBondParams::default(),
        ReplaceBondParams {
            preserve_props: false,
            keep_s_groups: true
        }
    );

    let mut batch = RWMol::from(mol_from_smiles("OCCN").unwrap());
    batch.begin_batch_edit().unwrap();
    assert!(matches!(batch.begin_batch_edit(), Err(Error::Value { .. })));
    batch.remove_atom(0).unwrap();
    batch.remove_atom(3).unwrap();
    // Removals wait for the commit, so indices stay valid.
    assert_eq!(batch.get_num_atoms(), 4);
    batch.commit_batch_edit().unwrap();
    assert_eq!(mol_to_smiles(&batch).unwrap(), "CC");
    batch.begin_batch_edit().unwrap();
    batch.remove_atom(0).unwrap();
    batch.rollback_batch_edit();
    assert_eq!(batch.get_num_atoms(), 2);
}
