// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_substruct_library::CachedMolHolder;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::rd_geometry::Point3D;

// The default is process-wide, so everything that depends on it is in one test.
#[test]
fn pickles_round_trip() {
    let mut mol = mol_from_smiles("c1ccccc1O").unwrap();
    mol.set_prop("note", "kept");
    let pickle = mol.to_binary();
    // SAFETY: a pickle RDKit wrote.
    let back = unsafe { Mol::from_pkl_string(&pickle) }.unwrap();
    assert_eq!(mol_to_smiles(&back).unwrap(), "Oc1ccccc1");
    assert!(!back.has_prop("note"));

    let all = mol.to_binary_with_property_flags(PropertyPickleOptions::AllProps);
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { Mol::from_pkl_string(&all) }
            .unwrap()
            .get_prop("note")
            .unwrap(),
        "kept"
    );
    let no_props =
        // SAFETY: a pickle RDKit wrote.
        unsafe { Mol::from_pkl_string_with_property_flags(&all, PropertyPickleOptions::NoProps) }.unwrap();
    assert!(!no_props.has_prop("note"));
    // SAFETY: a pickle RDKit wrote.
    let rw = unsafe { RWMol::from_pkl_string(&all) }.unwrap();
    assert_eq!(rw.get_num_atoms(), 7);
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        unsafe {
            RWMol::from_pkl_string_with_property_flags(&all, PropertyPickleOptions::AllProps)
        }
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
        // SAFETY: a pickle RDKit wrote.
        unsafe { Mol::from_pkl_string(&mol.to_binary()) }
            .unwrap()
            .get_prop("note")
            .unwrap(),
        "kept"
    );
    set_default_pickle_properties(original);

    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { Mol::from_pkl_string(b"not a pickle") },
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
        // SAFETY: a pickle RDKit wrote.
        unsafe { Mol::from_pkl_string(&pickle) }.unwrap();
        // RDKit refuses a negative byte, but not one past its enum.
        pickle[at] = 100;
        // SAFETY: a pickle RDKit wrote but for an enum byte, which check_pickled_mol checks.
        match unsafe { Mol::from_pkl_string(&pickle) } {
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
    // SAFETY: a pickle RDKit wrote but for an enum byte, which check_pickled_mol checks.
    match unsafe { Mol::from_pkl_string(&pickle) } {
        Err(Error::Value { message }) => assert!(message.contains("StereoGroupType"), "{message}"),
        other => panic!("StereoGroupType: {other:?}"),
    }
}

// RDKit 2024.03's StereoGroup keeps an atom listed twice, which a removal
// leaves dangling, so unpickling refuses a group that lists one twice. The
// atom's byte is found by pickling two other groups.
#[test]
fn pickles_with_stereo_groups_listing_an_atom_twice_are_refused() {
    let pickles: Vec<Vec<u8>> = [4, 3]
        .into_iter()
        .map(|second| {
            let mut rw = RWMol::from(mol_from_smiles("C[C@H](F)C[C@H](O)C").unwrap());
            let group = create_stereo_group(StereoGroupType::StereoOr, &rw, &[1, second]).unwrap();
            rw.set_stereo_groups(&[group]).unwrap();
            rw.to_binary_with_property_flags(PropertyPickleOptions::NoProps)
        })
        .collect();
    let (at, four, _) = differing_byte(&pickles[0], &pickles[1]);
    assert_eq!(four, 4);
    // SAFETY: a pickle RDKit wrote.
    unsafe { Mol::from_pkl_string(&pickles[0]) }.unwrap();
    let mut pickle = pickles[0].clone();
    pickle[at] = 1;
    // SAFETY: a pickle RDKit wrote but for an atom index, which
    // check_pickled_mol checks.
    assert_eq!(
        unsafe { Mol::from_pkl_string(&pickle) }.unwrap_err(),
        Error::Value {
            message: "Duplicate atom in StereoGroup".into()
        }
    );
}

// MolPickler writes these values in one byte, and its reader takes the rest
// of a longer list for the pickle's next fields: Python's ToBinary writes
// pickles it cannot read back ("Bad pickle format"). RDKit 2024.03 also
// writes bond indices in one byte for a molecule of at most 255 atoms.

fn brackets(n: usize) -> Mol {
    let mut mol = mol_from_smiles("CCO").unwrap();
    let mut group = create_mol_substance_group(&mut mol, "GEN");
    group.add_atom_with_idx(0).unwrap();
    for _ in 0..n {
        group
            .add_bracket(&[Point3D::new(0.0, 0.0, 0.0), Point3D::new(1.0, 0.0, 0.0)])
            .unwrap();
    }
    mol
}

fn stereo_groups(n: usize) -> Mol {
    let mut rw = RWMol::from(mol_from_smiles("C[C@H](F)Cl").unwrap());
    let group = create_stereo_group(StereoGroupType::StereoOr, &rw, &[1]).unwrap();
    rw.set_stereo_groups(&vec![group; n]).unwrap();
    rw.into()
}

/// An atom query of `n` children, which RDKit's JSON reader makes as given.
fn query_children(n: usize) -> Mol {
    let child = r#"{"descr": "AtomType", "tag": 33, "val": 6}"#;
    let children = vec![child; n].join(",");
    let json = format!(
        r#"{{"commonchem": {{"version": 10}}, "defaults": {{"atom": {{"z": 6, "impHs": 0, "chg": 0, "nRad": 0, "isotope": 0, "stereo": "unspecified"}}, "bond": {{"bo": 1, "stereo": "unspecified"}}}}, "molecules": [{{"atoms": [{{"z": 6}}, {{"z": 6}}], "bonds": [{{"atoms": [0, 1]}}], "extensions": [{{"name": "rdkitQueries", "formatVersion": 10, "atomQueries": [{{"descr": "AtomOr", "tag": 31, "children": [{children}]}}, {child}], "bondQueries": [{{"descr": "SingleOrAromaticBond", "type": "BondOrder", "tag": 33, "val": 1}}]}}]}}]}}"#
    );
    json_to_mols(&json).unwrap().remove(0)
}

/// The complete graph on 30 atoms (435 bonds), with a substance group of
/// `bond`.
fn complete_graph_with_group_bond(bond: u32) -> Mol {
    let mut rw = RWMol::default();
    for _ in 0..30 {
        rw.add_atom(&Atom::from_num(6));
    }
    for i in 0..30 {
        for j in i + 1..30 {
            rw.add_bond_with_order(i, j, BondType::Single).unwrap();
        }
    }
    let mut mol: Mol = rw.into();
    create_mol_substance_group(&mut mol, "GEN")
        .add_bond_with_idx(bond)
        .unwrap();
    mol
}

/// Makes a molecule with `n` of what a case counts.
type MakeMol = fn(usize) -> Mol;

#[test]
fn molecules_a_pickle_cannot_hold() {
    let cases: [(&str, MakeMol); 3] = [
        ("number of brackets", brackets),
        ("number of stereo groups", stereo_groups),
        ("number of children", query_children),
    ];
    for (what, make) in cases {
        let fits = make(255);
        // SAFETY: a pickle RDKit wrote.
        let back = unsafe { Mol::from_pkl_string(&fits.to_binary()) }.unwrap();
        assert_eq!(back.get_num_atoms(), fits.get_num_atoms(), "{what}");
        let mut holder = CachedMolHolder::default();
        holder.add_mol(&fits).unwrap();
        assert_eq!(
            holder.get_mol(0).unwrap().get_num_atoms(),
            fits.get_num_atoms()
        );
        match holder.add_mol(&make(256)) {
            Err(Error::Value { message }) => assert!(message.contains(what), "{message}"),
            other => panic!("{what}: {other:?}"),
        }
    }

    let fits = complete_graph_with_group_bond(255);
    // SAFETY: a pickle RDKit wrote.
    let back = unsafe { Mol::from_pkl_string(&fits.to_binary()) }.unwrap();
    assert_eq!(get_mol_substance_groups(&back)[0].get_bonds(), [255]);
    // Python's reads back bond 44.
    let mut holder = CachedMolHolder::default();
    assert!(matches!(
        holder.add_mol(&complete_graph_with_group_bond(300)),
        Err(Error::Value { .. })
    ));
}

#[test]
#[should_panic(expected = "number of brackets")]
fn to_binary_panics_for_what_a_pickle_cannot_hold() {
    brackets(256).to_binary();
}
