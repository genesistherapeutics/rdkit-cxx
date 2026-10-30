// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn stereo_groups_from_cxsmiles() {
    let mol = mol_from_smiles("C[C@H](F)C[C@H](O)C |&1:1,4|").unwrap();
    let groups = mol.get_stereo_groups();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].get_group_type(), StereoGroupType::StereoAnd);
    assert_eq!(groups[0].get_atoms(), [1, 4]);
    assert!(groups[0].get_bonds().is_empty());
    assert_eq!(groups[0].get_read_id(), 1);
    assert!(
        mol_from_smiles("CC")
            .unwrap()
            .get_stereo_groups()
            .is_empty()
    );
}

#[test]
fn create_and_set_stereo_groups() {
    let mut rw = RWMol::from(mol_from_smiles("C[C@H](F)C[C@H](O)C").unwrap());
    let mut group = create_stereo_group(StereoGroupType::StereoOr, &rw, &[1, 4]).unwrap();
    group.set_write_id(3);
    assert_eq!(group.get_write_id(), 3);
    rw.set_stereo_groups(&[group.clone()]).unwrap();
    assert_eq!(rw.get_stereo_groups(), [group.clone()]);

    assert!(matches!(
        create_stereo_group(StereoGroupType::StereoAnd, &rw, &[]),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        create_stereo_group(StereoGroupType::StereoAnd, &rw, &[1, 40]).unwrap_err(),
        Error::Index { index: 40, len: 7 }
    );
    let params = CreateStereoGroupParams {
        atom_ids: vec![1],
        read_id: 5,
        ..Default::default()
    };
    let with_id =
        create_stereo_group_with_params(StereoGroupType::StereoAnd, &rw, &params).unwrap();
    assert_eq!(with_id.get_read_id(), 5);

    // Groups hold indices, so a group checked against one molecule is checked
    // again against the one it is set on.
    let mut small = RWMol::from(mol_from_smiles("CC").unwrap());
    small.set_stereo_groups(&[with_id]).unwrap();
    assert_eq!(
        small.set_stereo_groups(&[group]).unwrap_err(),
        Error::Index { index: 4, len: 2 }
    );
}

#[test]
fn forwarding_ids() {
    let mut mol = mol_from_smiles("C[C@H](F)C[C@H](O)C |&7:1,4|").unwrap();
    forward_stereo_group_ids(&mut mol);
    assert_eq!(mol.get_stereo_groups()[0].get_write_id(), 7);
}

#[test]
fn potential_stereo() {
    let mut mol = mol_from_smiles("C[C@H](F)Cl").unwrap();
    let found = find_potential_stereo(&mut mol).unwrap();
    assert_eq!(found.len(), 1);
    let info = &found[0];
    assert_eq!(info.type_, StereoType::AtomTetrahedral);
    assert_eq!(info.specified, StereoSpecified::Specified);
    assert_eq!(info.centered_on, 1);
    assert_eq!(info.descriptor, StereoDescriptor::TetCcw);
    assert_eq!(info.controlling_atoms, [0, 2, 3]);
    let mut double = mol_from_smiles("CC=CC").unwrap();
    let params = FindPotentialStereoParams {
        flag_possible: true,
        ..Default::default()
    };
    let bonds = find_potential_stereo_with_params(&mut double, &params).unwrap();
    assert_eq!(bonds[0].type_, StereoType::BondDouble);
    assert_eq!(bonds[0].specified, StereoSpecified::Unspecified);
    assert_eq!(StereoInfo::NOATOM, Atom::NOATOM);
}

#[test]
fn stereo_info_default() {
    // Python's StereoInfo().
    let info = StereoInfo::default();
    assert_eq!(info.type_, StereoType::Unspecified);
    assert_eq!(info.specified, StereoSpecified::Unspecified);
    assert_eq!(info.centered_on, 4294967295);
    assert_eq!(info.descriptor, StereoDescriptor::NoValue);
    assert_eq!(info.permutation, 0);
    assert!(info.controlling_atoms.is_empty());
}

// RDKit 2025.09 accepts an atom listed twice in a stereo group, from its
// constructor and its CXSMILES reader, and removing the atom erases one
// entry and leaves a pointer to the deleted atom in the other. Both are
// refused (check_cx_extensions for the reader), as RDKit 2026.03 does.
#[test]
fn stereo_groups_hold_each_atom_once() {
    let mol = mol_from_smiles("C[C@H](F)Cl").unwrap();
    assert!(matches!(
        create_stereo_group(StereoGroupType::StereoOr, &mol, &[1, 1]),
        Err(Error::Value { .. })
    ));
    for text in ["C[C@H](F)Cl |o1:1,1|", "[H]C[C@H](F)Cl |o1:2,2|"] {
        assert!(
            matches!(mol_from_smiles(text), Err(Error::Value { .. })),
            "{text}"
        );
    }
    let parsed = mol_from_smiles("C[C@H](F)Cl |o1:1|").unwrap();
    assert_eq!(parsed.get_stereo_groups().len(), 1);
    let mut rw = RWMol::from(parsed);
    rw.remove_atom(1).unwrap();
    assert!(rw.get_stereo_groups().is_empty());
}

// RDKit 2025.09's setStereoGroups merges ABS groups into one, which holds an
// atom of two of them twice (rdkit-2024.03's "two ABS stereo groups").
#[test]
fn abs_groups_share_no_atom() {
    let mut rw = RWMol::from(mol_from_smiles("[*:1]C").unwrap());
    let abs = create_stereo_group(StereoGroupType::StereoAbsolute, &rw, &[0]).unwrap();
    match rw.set_stereo_groups(&[abs.clone(), abs.clone()]) {
        Err(Error::Value { message }) => {
            assert!(message.contains("two ABS stereo groups"), "{message}")
        }
        other => panic!("{other:?}"),
    }
    assert!(rw.get_stereo_groups().is_empty());
    rw.set_stereo_groups(&[abs]).unwrap();
    let a = Mol::from(rw);
    let mut partner = RWMol::from(mol_from_smiles("[*:1][C@H](F)Cl").unwrap());
    let or = create_stereo_group(StereoGroupType::StereoOr, &partner, &[1]).unwrap();
    partner.set_stereo_groups(&[or]).unwrap();
    let zipped = rdmolops::molzip(&a, &Mol::from(partner)).unwrap().unwrap();
    // The ABS group's one atom is the dummy molzip removes.
    let groups = zipped.get_stereo_groups();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].get_group_type(), StereoGroupType::StereoOr);
}
