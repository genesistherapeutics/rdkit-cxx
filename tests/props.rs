// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn typed_getters_and_errors() {
    let mut mol = mol_from_smiles("OCC").unwrap();
    mol.set_prop("s", "abc");
    mol.set_int_prop("i", 7);
    assert!(mol.has_prop("s") && !mol.has_prop("nope"));
    assert_eq!(
        mol.get_prop("nope").unwrap_err(),
        Error::Key { key: "nope".into() }
    );
    assert_eq!(
        mol.get_int_prop("nope").unwrap_err(),
        Error::Key { key: "nope".into() }
    );
    assert!(
        matches!(mol.get_int_prop("s"), Err(Error::Value { message }) if message.contains("integer"))
    );
    assert_eq!(mol.get_prop("i").unwrap(), "7");
    assert_eq!(mol.get_int_prop("i").unwrap(), 7);
    assert!(matches!(mol.get_double_prop("i"), Err(Error::Value { .. })));
    mol.clear_prop("nope");
    assert_eq!(mol.get_prop_names().unwrap(), ["s", "i"]);
    mol.clear_prop("s");
    assert_eq!(mol.get_prop_names().unwrap(), ["i"]);
}

#[test]
fn every_type_round_trips() {
    let mut mol = Mol::default();
    mol.set_unsigned_prop("u", 3);
    mol.set_double_prop("d", 0.25);
    mol.set_bool_prop("b", true);
    assert_eq!(mol.get_unsigned_prop("u").unwrap(), 3);
    assert_eq!(mol.get_double_prop("d").unwrap(), 0.25);
    assert!(mol.get_bool_prop("b").unwrap());
    assert_eq!(
        mol.get_prop_auto_convert("d").unwrap(),
        PropValue::Double(0.25)
    );
    assert_eq!(
        mol.get_prop_auto_convert("u").unwrap(),
        PropValue::Unsigned(3)
    );
    assert_eq!(
        mol.get_prop_auto_convert("b").unwrap(),
        PropValue::Bool(true)
    );
    mol.set_prop("n", "12");
    mol.set_prop("x", "1.5");
    mol.set_prop("t", "text");
    assert_eq!(mol.get_prop_auto_convert("n").unwrap(), PropValue::Int(12));
    assert_eq!(
        mol.get_prop_auto_convert("x").unwrap(),
        PropValue::Double(1.5)
    );
    assert_eq!(
        mol.get_prop_auto_convert("t").unwrap(),
        PropValue::String("text".into())
    );
    assert_eq!(
        mol.get_prop_auto_convert("nope").unwrap_err(),
        Error::Key { key: "nope".into() }
    );
}

#[test]
fn computed_properties() {
    let mut mol = mol_from_smiles("C").unwrap();
    mol.set_int_prop_with_computed("c", 1, true);
    assert!(mol.get_prop_names().unwrap().is_empty());
    let computed = GetPropNamesParams {
        include_computed: true,
        ..Default::default()
    };
    // Sanitization also sets computed properties, such as numArom.
    assert!(
        mol.get_prop_names_with_params(&computed)
            .unwrap()
            .contains(&"c".to_string())
    );
    mol.clear_computed_props();
    assert!(!mol.has_prop("c"));
    mol.set_prop_with_computed("c", "x", true);
    mol.clear_computed_props_with_include_rings(false);
    assert!(!mol.has_prop("c"));
}

#[test]
fn writing_smiles_sets_computed_properties() {
    let mol = mol_from_smiles("OCC").unwrap();
    mol_to_smiles(&mol).unwrap();
    let all = GetPropNamesParams {
        include_private: true,
        include_computed: true,
    };
    let names = mol.get_prop_names_with_params(&all).unwrap();
    assert!(names.contains(&"_smilesAtomOutputOrder".to_string()));
    let dict = mol
        .get_props_as_dict_with_params(&GetPropsAsDictParams {
            include_private: true,
            include_computed: true,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        dict["_smilesAtomOutputOrder"],
        PropValue::UnsignedVec(vec![2, 1, 0])
    );
}

#[test]
fn props_as_dict_converts_strings() {
    let mut mol = mol_from_smiles("C").unwrap();
    mol.set_prop("s", "abc");
    mol.set_int_prop("i", 7);
    mol.set_prop("n", "5");
    let expected = BTreeMap::from([
        ("i".to_string(), PropValue::Int(7)),
        ("n".to_string(), PropValue::Int(5)),
        ("s".to_string(), PropValue::String("abc".into())),
    ]);
    assert_eq!(mol.get_props_as_dict().unwrap(), expected);
    assert_eq!(
        GetPropsAsDictParams::default(),
        GetPropsAsDictParams {
            include_private: false,
            include_computed: false,
            auto_convert_strings: true,
        }
    );
    let raw = GetPropsAsDictParams {
        auto_convert_strings: false,
        ..Default::default()
    };
    assert_eq!(
        mol.get_props_as_dict_with_params(&raw).unwrap()["n"],
        PropValue::String("5".into())
    );
}

#[test]
fn atom_and_bond_props() {
    let mut atom = Atom::from_num(6);
    atom.set_prop("x", "y");
    assert_eq!(atom.get_prop("x").unwrap(), "y");
    let mut rw = RWMol::default();
    rw.add_atom(&atom);
    assert_eq!(rw.get_atom_with_idx(0).unwrap().get_prop("x").unwrap(), "y");
    let mut view = rw.get_atom_with_idx_mut(0).unwrap();
    view.set_bool_prop("flag", true);
    assert!(view.get_bool_prop("flag").unwrap());

    let mut mol = mol_from_smiles("CC").unwrap();
    mol.get_bond_with_idx_mut(0)
        .unwrap()
        .set_double_prop("w", 0.5);
    assert_eq!(
        mol.get_bond_with_idx(0)
            .unwrap()
            .get_double_prop("w")
            .unwrap(),
        0.5
    );
    assert_eq!(
        mol.get_bond_with_idx(0).unwrap().get_prop_names().unwrap(),
        ["w"]
    );
}

#[test]
fn latin1_name_is_a_value_error() {
    let mol = mol_from_mol_file(common::latin1_mol_file()).unwrap();
    assert!(mol.has_prop("_Name"));
    assert!(matches!(mol.get_prop("_Name"), Err(Error::Value { .. })));
}
