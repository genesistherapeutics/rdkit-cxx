// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn option_values_match_rdkit() {
    let (rust, cxx) = __property_pickle_option_values_for_tests();
    assert_eq!(rust, cxx);
    assert_eq!(
        PropertyPickleOptions::QueryAtomData,
        PropertyPickleOptions::AtomProps
    );
    let combined = PropertyPickleOptions::MolProps | PropertyPickleOptions::AtomProps;
    assert_eq!(combined.bits(), 0x3);
    assert_eq!(
        format!("{:?}", PropertyPickleOptions::AllProps),
        "PropertyPickleOptions::AllProps"
    );
}

// The default is process-wide, so everything that depends on it is in one test.
#[test]
fn pickles_round_trip() {
    let mut mol = mol_from_smiles("c1ccccc1O").unwrap();
    mol.set_prop("note", "kept");
    let pickle = mol.to_binary();
    let back = Mol::from_pkl_string(&pickle).unwrap();
    assert_eq!(mol_to_smiles(&back).unwrap(), "Oc1ccccc1");
    assert!(!back.has_prop("note"));

    let all = mol.to_binary_with_property_flags(PropertyPickleOptions::AllProps);
    assert_eq!(
        Mol::from_pkl_string(&all)
            .unwrap()
            .get_prop("note")
            .unwrap(),
        "kept"
    );
    let no_props =
        Mol::from_pkl_string_with_property_flags(&all, PropertyPickleOptions::NoProps).unwrap();
    assert!(!no_props.has_prop("note"));
    let rw = RWMol::from_pkl_string(&all).unwrap();
    assert_eq!(rw.get_num_atoms(), 7);
    assert_eq!(
        RWMol::from_pkl_string_with_property_flags(&all, PropertyPickleOptions::AllProps)
            .unwrap()
            .get_num_atoms(),
        7
    );

    let original = get_default_pickle_properties();
    set_default_pickle_properties(PropertyPickleOptions::MolProps);
    assert_eq!(
        get_default_pickle_properties(),
        PropertyPickleOptions::MolProps
    );
    assert_eq!(
        Mol::from_pkl_string(&mol.to_binary())
            .unwrap()
            .get_prop("note")
            .unwrap(),
        "kept"
    );
    set_default_pickle_properties(original);

    assert!(matches!(
        Mol::from_pkl_string(b"not a pickle"),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn atom_labels() {
    let mut mol = mol_from_smiles("C*").unwrap();
    let mut dummy = mol.get_atom_with_idx_mut(1).unwrap();
    set_atom_r_label(&mut dummy, 2).unwrap();
    set_atom_alias(&mut dummy, "R2");
    set_atom_value(&mut dummy, "v");
    set_supplemental_smiles_label(&mut dummy, "xxx");
    let dummy = mol.get_atom_with_idx(1).unwrap();
    assert_eq!(get_atom_r_label(&dummy), 2);
    assert_eq!(get_atom_alias(&dummy).unwrap(), "R2");
    assert_eq!(get_atom_value(&dummy).unwrap(), "v");
    assert_eq!(get_supplemental_smiles_label(&dummy).unwrap(), "xxx");
    assert!(mol_to_smiles(&mol).unwrap().contains("xxx"));

    let mut owned = Atom::from_num(0);
    set_atom_alias(&mut owned, "Q");
    assert_eq!(get_atom_alias(&owned).unwrap(), "Q");
    assert_eq!(get_atom_r_label(&owned), 0);
    assert!(matches!(
        set_atom_r_label(&mut owned, 500),
        Err(Error::Runtime { .. })
    ));

    let benzene = mol_from_smiles("c1ccccc1").unwrap();
    assert_eq!(
        get_num_pi_electrons(&benzene.get_atom_with_idx(0).unwrap()).unwrap(),
        1
    );
    assert!(matches!(
        get_num_pi_electrons(&Atom::from_num(6)),
        Err(Error::Runtime { .. })
    ));
}
