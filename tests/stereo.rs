// Copyright (c) 2026 Steven Kearnes All rights reserved.
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
