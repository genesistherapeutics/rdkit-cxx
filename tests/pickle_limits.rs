// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// MolPickler counts each object's properties in a u16 (streamWriteProps in
// RDGeneral/StreamOps.h): past 65535 the count wraps and Python's ToBinary
// writes a pickle it cannot read back (2026.09.1: "Bond index out of range"
// for 65536 properties; 65535 user properties beside a computed one are
// dropped), so check_picklable refuses them (mol_checks.cpp).

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chem_reactions::ChemicalReaction;
use rdkit_cxx::chem::rd_substruct_library::CachedMolHolder;
use rdkit_cxx::chem::*;

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
fn to_binary_panics_for_more_than_65535_atom_properties() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    let mut atom = mol.get_atom_with_idx_mut(1).unwrap();
    for i in 0..65536 {
        atom.set_prop(&format!("p{i}"), "x");
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
    // Substance groups' properties are pickled whatever the flags.
    mol.to_binary();
}

#[test]
#[should_panic(expected = "more than a molecule pickle counts")]
fn reaction_to_binary_panics_for_more_than_65535_template_properties() {
    let mut template = mol_from_smarts("[C:1]").unwrap();
    for i in 0..65536 {
        template.set_prop(&format!("p{i}"), "x");
    }
    let mut rxn = ChemicalReaction::default();
    rxn.add_reactant_template(&template);
    rxn.add_product_template(&mol_from_smarts("[C:1]O").unwrap());
    rxn.to_binary();
}

#[test]
fn pickling_wrappers_refuse_too_many_properties() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    let mut group = create_mol_substance_group(&mut mol, "DAT");
    group.add_atom_with_idx(0).unwrap();
    for i in 0..65536 {
        group.set_prop(&format!("p{i}"), "x");
    }
    let mut holder = CachedMolHolder::default();
    assert!(matches!(holder.add_mol(&mol), Err(Error::Value { .. })));
}

// 65535 user properties beside the parser's computed ones: RDKit writes all
// of them with AllProps, and its count wraps.
#[test]
#[should_panic(expected = "more than a molecule pickle counts")]
fn to_binary_panics_for_65535_user_properties_beside_computed_ones() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    for i in 0..65535 {
        mol.set_prop(&format!("p{i}"), "x");
    }
    mol.to_binary_with_property_flags(PropertyPickleOptions::AllProps);
}

#[test]
fn many_properties_within_the_count_round_trip() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    for i in 0..65000 {
        mol.set_prop(&format!("p{i}"), "x");
    }
    let pkl = mol.to_binary_with_property_flags(PropertyPickleOptions::AllProps);
    // SAFETY: a pickle RDKit wrote.
    let back = unsafe { Mol::from_pkl_string(&pkl) }.unwrap();
    assert_eq!(back.get_prop("p64999").unwrap(), "x");
}
