// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// RDKit 2026.03's pickles: the property counts it writes in two bytes
// (mol_checks.cpp).

use rdkit_cxx::chem::*;

// MolPickler counts each object's properties in a u16: past 65535 the
// count wraps and Python's ToBinary writes a pickle it cannot read back.
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
