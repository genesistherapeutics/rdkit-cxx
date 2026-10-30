// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

#[test]
fn substance_groups_from_cxsmiles() {
    let mol = mol_from_smiles("*-CC(-*)C |$star_e;;;star_e;$,Sg:n:1,2::ht:::|").unwrap();
    let groups = get_mol_substance_groups(&mol);
    assert_eq!(groups.len(), 1);
    let group = groups[0];
    assert_eq!(group.get_index_in_mol(), 0);
    assert_eq!(group.get_prop("TYPE").unwrap(), "SRU");
    assert_eq!(group.get_atoms(), [1, 2]);
    assert_eq!(group.get_bonds(), [0, 2, 3]);
    assert!(std::ptr::eq(group.get_owning_mol(), &mol));
    assert_eq!(group.get_unsigned_vect_prop("XBCORR").unwrap(), [0, 2]);

    // Unlike Mol's, SubstanceGroup's GetPropsAsDict includes private and
    // computed properties by default.
    let props = group.get_props_as_dict().unwrap();
    assert_eq!(props["CONNECT"], PropValue::String("HT".into()));
    assert_eq!(props["_cxsmilesindex"], PropValue::Unsigned(0));
    assert_eq!(
        group.get_prop_names().unwrap(),
        ["TYPE", "CONNECT", "XBHEAD", "XBCORR", "index"]
    );

    let same = get_mol_substance_group_with_idx(&mol, 0).unwrap();
    assert_eq!(same.get_atoms(), group.get_atoms());
    assert_eq!(
        get_mol_substance_group_with_idx(&mol, 1).unwrap_err(),
        Error::Index { index: 1, len: 1 }
    );
    assert!(matches!(
        group.get_int_prop("TYPE"),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        group.get_string_vect_prop("missing"),
        Err(Error::Key { .. })
    ));
}

#[test]
fn create_and_edit_substance_groups() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    let data = create_mol_data_substance_group(&mut mol, "pKa", "4.5");
    assert_eq!(data.get_prop("TYPE").unwrap(), "DAT");
    assert_eq!(data.get_prop("FIELDNAME").unwrap(), "pKa");
    assert_eq!(data.get_string_vect_prop("DATAFIELDS").unwrap(), ["4.5"]);

    let mut sru = create_mol_substance_group(&mut mol, "SRU");
    assert_eq!(sru.get_index_in_mol(), 1);
    sru.add_atom_with_idx(1).unwrap();
    sru.add_atom_with_idx(2).unwrap();
    sru.add_bond_with_idx(0).unwrap();
    sru.add_parent_atom_with_idx(1).unwrap();
    sru.add_attach_point(1, 0, "1").unwrap();
    sru.add_c_state(0, Point3D::new(1.0, 0.0, 0.0)).unwrap();
    sru.add_bracket(&[Point3D::new(0.0, 0.0, 0.0), Point3D::new(1.0, 1.0, 0.0)])
        .unwrap();
    sru.set_int_prop("count", 3);
    assert_eq!(sru.get_int_prop("count").unwrap(), 3);

    assert_eq!(sru.get_atoms(), [1, 2]);
    assert_eq!(sru.get_parent_atoms(), [1]);
    assert_eq!(
        sru.get_attach_points().unwrap(),
        [SubstanceGroupAttach {
            a_idx: 1,
            lv_idx: 0,
            id: "1".into()
        }]
    );
    assert_eq!(
        sru.get_c_states(),
        [SubstanceGroupCState {
            bond_idx: 0,
            vector: Point3D::new(1.0, 0.0, 0.0)
        }]
    );
    assert_eq!(
        sru.get_brackets(),
        [[
            Point3D::new(0.0, 0.0, 0.0),
            Point3D::new(1.0, 1.0, 0.0),
            Point3D::default()
        ]]
    );

    // RDKit's own checks.
    assert!(matches!(
        sru.add_bracket(&[Point3D::default()]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        sru.add_c_state(1, Point3D::default()),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        sru.add_parent_atom_with_idx(0),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        sru.add_atom_with_bookmark(3),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        sru.add_bond_with_bookmark(3),
        Err(Error::Runtime { .. })
    ));
    // The wrappers' index checks (E3).
    assert_eq!(
        sru.add_atom_with_idx(3).unwrap_err(),
        Error::Index { index: 3, len: 3 }
    );
    assert_eq!(
        sru.add_bond_with_idx(2).unwrap_err(),
        Error::Index { index: 2, len: 2 }
    );
    assert_eq!(
        sru.set_atoms(&[0, 5]).unwrap_err(),
        Error::Index { index: 5, len: 3 }
    );
    assert_eq!(
        sru.add_attach_point(1, 7, "2").unwrap_err(),
        Error::Index { index: 7, len: 3 }
    );
    sru.add_attach_point(2, -1, "2").unwrap();

    sru.set_atoms(&[0, 1]).unwrap();
    sru.set_parent_atoms(&[0]).unwrap();
    sru.set_bonds(&[0, 1]).unwrap();
    assert_eq!(sru.get_atoms(), [0, 1]);
    assert_eq!(sru.get_parent_atoms(), [0]);
    assert_eq!(sru.get_bonds(), [0, 1]);
    sru.clear_brackets();
    sru.clear_c_states();
    sru.clear_attach_points();
    assert!(sru.get_brackets().is_empty());
    assert!(sru.get_c_states().is_empty());
    assert!(sru.get_attach_points().unwrap().is_empty());

    let block = mol_to_mol_block(&mol).unwrap();
    assert!(block.contains("M  STY  2   1 DAT   2 SRU"));
    let read = mol_from_mol_block(&block).unwrap();
    assert_eq!(get_mol_substance_groups(&read).len(), 2);

    let mut edit = get_mol_substance_group_with_idx_mut(&mut mol, 0).unwrap();
    edit.set_prop("FIELDNAME", "pKb");
    edit.clear_prop("missing");
    assert_eq!(
        get_mol_substance_group_with_idx(&mol, 0)
            .unwrap()
            .get_prop("FIELDNAME")
            .unwrap(),
        "pKb"
    );
    assert!(get_mol_substance_group_with_idx_mut(&mut mol, 2).is_err());

    clear_mol_substance_groups(&mut mol);
    assert!(get_mol_substance_groups(&mol).is_empty());
}

#[test]
fn removing_an_atom_drops_its_substance_groups() {
    let mut rw = RWMol::from(mol_from_smiles("CCO").unwrap());
    create_mol_data_substance_group(&mut rw, "pKa", "4.5");
    let mut sru = create_mol_substance_group(&mut rw, "SRU");
    sru.set_atoms(&[1, 2]).unwrap();
    rw.remove_atom(2).unwrap();
    let groups = get_mol_substance_groups(&rw);
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].get_prop("TYPE").unwrap(), "DAT");
}

#[test]
fn add_a_copy_from_another_molecule() {
    let mut source = mol_from_smiles("CCO").unwrap();
    let mut sru = create_mol_substance_group(&mut source, "SRU");
    sru.set_atoms(&[1, 2]).unwrap();
    sru.set_bonds(&[1]).unwrap();
    let sgroup = get_mol_substance_group_with_idx(&source, 0).unwrap();

    let mut target = mol_from_smiles("CCOC").unwrap();
    let copy = add_mol_substance_group(&mut target, &sgroup).unwrap();
    assert_eq!(copy.get_index_in_mol(), 0);
    assert_eq!(copy.get_atoms(), [1, 2]);

    // A copy's indices must exist in the molecule it joins.
    let mut small = mol_from_smiles("C").unwrap();
    assert_eq!(
        add_mol_substance_group(&mut small, &sgroup).unwrap_err(),
        Error::Index { index: 1, len: 1 }
    );
    assert!(get_mol_substance_groups(&small).is_empty());

    // Clones carry their substance groups, owned by the clone.
    let clone = target.clone();
    let cloned = get_mol_substance_groups(&clone);
    assert_eq!(cloned.len(), 1);
    assert!(std::ptr::eq(cloned[0].get_owning_mol(), &clone));
}

#[test]
fn defaults() {
    assert_eq!(
        SubstanceGroupCState::default(),
        SubstanceGroupCState {
            bond_idx: 0,
            vector: Point3D::default()
        }
    );
    assert_eq!(
        SubstanceGroupAttach::default(),
        SubstanceGroupAttach {
            a_idx: 0,
            lv_idx: 0,
            id: String::new()
        }
    );
}
