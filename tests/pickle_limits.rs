// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// RDKit 2024.09's pickles: the property counts it writes in two bytes, and
// stereo groups that list an atom twice, which 2024.09 keeps and its
// removals then leave pointing at the atom (mol_checks.cpp, stereo.cpp).

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn hex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

// RDKit 2024.09.6's Mol.ToBinary() of "C[C@H](F)Cl |o1:1,1|": its OR group
// lists atom 1 twice, which 2024.09's StereoGroup keeps (2026.03 refuses).
const DUPLICATE_GROUP_ATOM: &str = "efbeadde0000000010000000020000000000000004000000030000008001060060000000010306203400000002010409002000000001110020000000010b00010001020001030042000000003e0101010201010017040000000000000016";
// The same molecule with the group listing atom 1 once.
const GROUP_ATOM_ONCE: &str = "efbeadde0000000010000000020000000000000004000000030000008001060060000000010306203400000002010409002000000001110020000000010b00010001020001030042000000003e01010101010017040000000000000016";

#[test]
fn a_pickled_stereo_group_may_not_list_an_atom_twice() {
    // SAFETY: a pickle RDKit 2024.09 wrote.
    let once = unsafe { Mol::from_pkl_string(&hex(GROUP_ATOM_ONCE)) }.unwrap();
    assert_eq!(once.get_stereo_groups()[0].get_atoms(), [1]);
    // SAFETY: a pickle RDKit 2024.09 wrote.
    match unsafe { Mol::from_pkl_string(&hex(DUPLICATE_GROUP_ATOM)) } {
        Err(Error::Value { message }) => assert!(message.contains("Duplicate atom"), "{message}"),
        other => panic!("expected Error::Value, got {other:?}"),
    }
}

#[test]
fn a_new_stereo_group_may_not_list_an_atom_twice() {
    let rw = RWMol::from(mol_from_smiles("C[C@H](F)Cl").unwrap());
    assert!(create_stereo_group(StereoGroupType::StereoOr, &rw, &[1]).is_ok());
    match create_stereo_group(StereoGroupType::StereoOr, &rw, &[1, 1]) {
        Err(Error::Value { message }) => assert!(message.contains("Duplicate atom"), "{message}"),
        other => panic!("expected Error::Value, got {other:?}"),
    }
}

// MolPickler counts each object's properties in a u16: past 65535 the
// count wraps and Python's ToBinary writes a pickle it cannot read back
// (2024.09.6: "Bad pickle format: ENDPROPS tag not found" for a molecule,
// "Bond index out of range" for a substance group).
#[test]
#[should_panic(expected = "more than a molecule pickle counts")]
fn to_binary_panics_for_more_than_65535_molecule_properties() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    for i in 0..65536 {
        mol.set_prop(&format!("p{i}"), "x");
    }
    mol.to_binary_with_property_flags(PropertyPickleOptions::AllProps);
}

#[test]
#[should_panic(expected = "more than a molecule pickle counts")]
fn to_binary_panics_for_more_than_65535_substance_group_properties() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    let mut group = create_mol_substance_group(&mut mol, "DAT");
    group.add_atom_with_idx(0).unwrap();
    for i in 0..65536 {
        group.set_prop(&format!("p{i}"), "x");
    }
    mol.to_binary();
}

#[test]
fn many_properties_within_the_count_round_trip() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    for i in 0..1000 {
        mol.set_prop(&format!("p{i}"), "x");
    }
    let pkl = mol.to_binary_with_property_flags(PropertyPickleOptions::AllProps);
    // SAFETY: a pickle RDKit wrote.
    let back = unsafe { Mol::from_pkl_string(&pkl) }.unwrap();
    assert_eq!(back.get_prop("p999").unwrap(), "x");
}
