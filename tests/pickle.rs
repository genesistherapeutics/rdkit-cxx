// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
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

// The one byte at which two pickles differ, with their values there.
fn differing_byte(a: &[u8], b: &[u8]) -> (usize, u8, u8) {
    assert_eq!(a.len(), b.len());
    let diffs: Vec<usize> = (0..a.len()).filter(|&i| a[i] != b[i]).collect();
    assert_eq!(diffs.len(), 1, "{diffs:?}");
    (diffs[0], a[diffs[0]], b[diffs[0]])
}

// MolPickler casts the bytes it reads for enums without checking them, and
// the Rust enums have no variant for a value RDKit's lack, so unpickling
// refuses one. A pickle leaves out a default value, so each byte is found
// by pickling two other values.
#[test]
fn pickles_with_unknown_enum_values_are_refused() {
    let flags = PropertyPickleOptions::NoProps;
    type Edit = fn(&mut Mol, bool);
    let edits: [(&str, Edit); 4] = [
        ("BondType", |m, other| {
            m.get_bond_with_idx_mut(0).unwrap().set_bond_type(if other {
                BondType::Triple
            } else {
                BondType::Double
            })
        }),
        ("BondDir", |m, other| {
            m.get_bond_with_idx_mut(0).unwrap().set_bond_dir(if other {
                BondDir::Unknown
            } else {
                BondDir::Beginwedge
            })
        }),
        ("ChiralType", |m, other| {
            m.get_atom_with_idx_mut(0)
                .unwrap()
                .set_chiral_tag(if other {
                    ChiralType::ChiOther
                } else {
                    ChiralType::ChiTetrahedralCw
                })
        }),
        ("HybridizationType", |m, other| {
            m.get_atom_with_idx_mut(1)
                .unwrap()
                .set_hybridization(if other {
                    HybridizationType::Sp2d
                } else {
                    HybridizationType::Sp
                })
        }),
    ];
    for (what, edit) in edits {
        let pickles: Vec<Vec<u8>> = [false, true]
            .into_iter()
            .map(|other| {
                let mut mol = mol_from_smiles("CC").unwrap();
                edit(&mut mol, other);
                mol.to_binary_with_property_flags(flags)
            })
            .collect();
        let (at, _, _) = differing_byte(&pickles[0], &pickles[1]);
        let mut pickle = pickles[1].clone();
        // The edited value itself reads back.
        Mol::from_pkl_string(&pickle).unwrap();
        // RDKit refuses a negative byte, but not one past its enum.
        pickle[at] = 100;
        match Mol::from_pkl_string(&pickle) {
            Err(Error::Value { message }) => assert!(message.contains(what), "{message}"),
            other => panic!("{what}: {other:?}"),
        }
    }

    // A stereo group's type, the same way.
    let pickles: Vec<Vec<u8>> = [StereoGroupType::StereoOr, StereoGroupType::StereoAnd]
        .into_iter()
        .map(|kind| {
            let mut rw = RWMol::from(mol_from_smiles("C[C@H](F)Cl").unwrap());
            let group = create_stereo_group(kind, &rw, &[1]).unwrap();
            rw.set_stereo_groups(&[group]).unwrap();
            rw.to_binary_with_property_flags(flags)
        })
        .collect();
    let (at, _, _) = differing_byte(&pickles[0], &pickles[1]);
    let mut pickle = pickles[1].clone();
    pickle[at] = 7;
    match Mol::from_pkl_string(&pickle) {
        Err(Error::Value { message }) => assert!(message.contains("StereoGroupType"), "{message}"),
        other => panic!("StereoGroupType: {other:?}"),
    }
}
