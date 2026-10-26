// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::rd_geometry::Point3D;

#[test]
fn missing_and_added_conformers() {
    let mut mol = mol_from_smiles("CO").unwrap();
    assert_eq!(mol.get_num_conformers(), 0);
    assert!(matches!(mol.get_conformer(),
                     Err(Error::Value { message }) if message.contains("No conformations")));
    let mut conf = Conformer::new(2);
    conf.set_atom_position(1, Point3D::new(1.5, 0.0, 0.0))
        .unwrap();
    assert!(!conf.has_owning_mol());
    assert_eq!(mol.add_conformer_with_assign_id(&conf, true).unwrap(), 0);
    let got = mol.get_conformer().unwrap();
    assert_eq!(
        got.get_atom_position(1).unwrap(),
        Point3D::new(1.5, 0.0, 0.0)
    );
    assert!(got.is_3d() && got.has_owning_mol());
    assert!(std::ptr::eq(got.get_owning_mol(), &mol));
    // A copy has no owner: RDKit's copy constructor would keep this
    // molecule, which the copy can outlive.
    let mut copy = Conformer::from(got);
    assert!(!copy.has_owning_mol() && !copy.clone().has_owning_mol());
    assert_eq!(copy.get_positions(), got.get_positions());
    assert!(copy.is_3d() && copy.get_id() == 0);
    copy.set_id(3).unwrap();
    copy.set_int_prop("n", 1);
    let again = copy.clone();
    assert_eq!(again.get_id(), 3);
    assert_eq!(again.get_int_prop("n").unwrap(), 1);
    assert_eq!(got.get_positions().len(), 2);
    assert_eq!(
        got.get_atom_position(5).unwrap_err(),
        Error::Index { index: 5, len: 2 }
    );
    assert_eq!(mol.get_conformer_with_id(0).unwrap().get_id(), 0);
    assert!(matches!(
        mol.get_conformer_with_id(7),
        Err(Error::Value { .. })
    ));
}

#[test]
fn conformer_must_match_atom_count() {
    let mut mol = mol_from_smiles("CO").unwrap();
    assert!(matches!(mol.add_conformer(&Conformer::new(3)),
                     Err(Error::Runtime { message }) if message.contains("mismatch")));
}

#[test]
fn edit_through_conformer_mut() {
    let block = mol_to_mol_block(&mol_from_smiles("CO").unwrap()).unwrap();
    let mut mol = mol_from_mol_block(&block).unwrap();
    assert_eq!(mol.get_num_conformers(), 1);
    let mut conf = mol.get_conformer_mut().unwrap();
    conf.set_atom_position(0, Point3D::new(0.0, 2.0, 0.0))
        .unwrap();
    conf.set_3d(true);
    conf.set_int_prop("n", 1);
    conf.set_int_prop_with_computed("c", 1, true);
    conf.clear_computed_props();
    assert!(!conf.has_prop("c"));
    assert_eq!(
        conf.get_atom_position(9).unwrap_err(),
        Error::Index { index: 9, len: 2 }
    );
    let conf = mol.get_conformer().unwrap();
    assert_eq!(conf.get_atom_position(0).unwrap().y, 2.0);
    assert!(conf.is_3d());
    assert_eq!(conf.get_int_prop("n").unwrap(), 1);
}

#[test]
fn positions_and_ids() {
    let mut conf = Conformer::new(2);
    let points = [Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)];
    conf.set_positions(&points).unwrap();
    assert_eq!(conf.get_positions(), points);
    assert!(matches!(
        conf.set_positions(&points[..1]),
        Err(Error::Value { .. })
    ));
    conf.set_id(4).unwrap();
    conf.set_3d(false);
    let copy = conf.clone();
    assert_eq!(
        (copy.get_id(), copy.get_num_atoms(), copy.is_3d()),
        (4, 2, false)
    );
    assert_eq!(Conformer::default().get_num_atoms(), 0);
    assert_eq!(Point3D::default(), Point3D::new(0.0, 0.0, 0.0));
}

#[test]
fn many_conformers() {
    let mut mol = mol_from_smiles("CO").unwrap();
    for i in 0..3 {
        let mut conf = Conformer::new(2);
        conf.set_atom_position(0, Point3D::new(i as f64, 0.0, 0.0))
            .unwrap();
        mol.add_conformer_with_assign_id(&conf, true).unwrap();
    }
    let ids: Vec<_> = mol.get_conformers().map(|c| c.get_id()).collect();
    assert_eq!(ids, [0, 1, 2]);
    assert_eq!(mol.get_conformers().rev().len(), 3);
    mol.get_conformer_with_id_mut(2)
        .unwrap()
        .set_atom_position(1, Point3D::new(9.0, 9.0, 9.0))
        .unwrap();
    assert_eq!(
        mol.get_conformer_with_id(2)
            .unwrap()
            .get_atom_position(1)
            .unwrap()
            .x,
        9.0
    );
    mol.remove_conformer(1);
    mol.remove_conformer(42);
    assert_eq!(
        mol.get_conformers().map(|c| c.get_id()).collect::<Vec<_>>(),
        [0, 2]
    );
    mol.remove_all_conformers();
    assert_eq!(mol.get_conformers().count(), 0);
}

/// A view finds its conformer by id, so it follows the id it sets, and
/// refuses an id that another conformer has, which it could not reach.
#[test]
fn conformer_mut_follows_its_new_id() {
    let mut mol = mol_from_smiles("CC").unwrap();
    for x in [1.0, 2.0] {
        let mut conf = Conformer::new(2);
        conf.set_atom_position(0, Point3D::new(x, 0.0, 0.0))
            .unwrap();
        mol.add_conformer_with_assign_id(&conf, true).unwrap();
    }
    let mut view = mol.get_conformer_with_id_mut(0).unwrap();
    view.set_id(5).unwrap();
    assert_eq!(view.get_id(), 5);
    assert_eq!(view.get_atom_position(0).unwrap().x, 1.0);
    view.set_atom_position(0, Point3D::new(3.0, 0.0, 0.0))
        .unwrap();
    assert!(matches!(view.set_id(1), Err(Error::Value { .. })));
    view.set_id(5).unwrap();
    assert_eq!(view.get_id(), 5);
    let xs: Vec<_> = mol
        .get_conformers()
        .map(|c| (c.get_id(), c.get_atom_position(0).unwrap().x))
        .collect();
    assert_eq!(xs, [(5, 3.0), (1, 2.0)]);
}

/// Lookups take an i32, where -1 is the default conformer, so larger ids are
/// refused rather than made unreachable.
#[test]
fn conformer_ids_fit_in_i32() {
    let max = i32::MAX as u32;
    let mut conf = Conformer::new(2);
    assert!(matches!(
        conf.set_id(3_000_000_000),
        Err(Error::Value { .. })
    ));
    assert!(matches!(conf.set_id(max + 1), Err(Error::Value { .. })));
    conf.set_id(max).unwrap();
    let mut mol = mol_from_smiles("CC").unwrap();
    assert_eq!(mol.add_conformer(&conf).unwrap(), max);
    assert_eq!(mol.get_conformer_with_id(i32::MAX).unwrap().get_id(), max);
    // RDKit would assign i32::MAX + 1, which overflows its int.
    assert!(matches!(
        mol.add_conformer_with_assign_id(&conf, true),
        Err(Error::Value { .. })
    ));
    let mut view = mol.get_conformer_mut().unwrap();
    assert!(matches!(view.set_id(max + 1), Err(Error::Value { .. })));
    assert_eq!(view.get_id(), max);
    assert_eq!(mol.get_num_conformers(), 1);
}
