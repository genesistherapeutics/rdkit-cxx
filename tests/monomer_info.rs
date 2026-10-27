// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn monomer_info() {
    assert_eq!(
        AtomMonomerInfo::default().get_monomer_type(),
        AtomMonomerType::Unknown
    );
    let params = AtomMonomerInfoParams { name: "N1".into() };
    let mut info = AtomMonomerInfo::new_with_params(AtomMonomerType::Other, &params).unwrap();
    assert_eq!(info.get_name().unwrap(), "N1");
    info.set_name("C1");
    info.set_monomer_type(AtomMonomerType::Unknown).unwrap();
    let copy = info.clone();
    assert_eq!(copy.get_name().unwrap(), "C1");
    assert_eq!(copy.get_monomer_type(), AtomMonomerType::Unknown);
    assert_eq!(
        AtomMonomerInfo::new(AtomMonomerType::Other)
            .unwrap()
            .get_name()
            .unwrap(),
        ""
    );
}

#[test]
fn pdb_residue_info() {
    let defaults = AtomPDBResidueInfoParams::default();
    assert_eq!(
        (
            defaults.serial_number,
            defaults.occupancy,
            defaults.temp_factor
        ),
        (1, 1.0, 0.0)
    );
    let mut info = AtomPDBResidueInfo::new(" CA ");
    assert_eq!(info.get_monomer_type(), AtomMonomerType::Pdbresidue);
    assert_eq!(
        (info.get_name().unwrap(), info.get_serial_number()),
        (" CA ".to_string(), 1)
    );
    info.set_residue_name("ALA");
    info.set_residue_number(7);
    info.set_chain_id("A");
    info.set_serial_number(12);
    info.set_alt_loc("B");
    info.set_insertion_code("C");
    info.set_occupancy(0.5);
    info.set_temp_factor(20.0);
    info.set_is_hetero_atom(true);
    info.set_secondary_structure(2);
    info.set_segment_number(3);
    let copy = info.clone();
    assert_eq!(copy.get_residue_name().unwrap(), "ALA");
    assert_eq!(
        (copy.get_residue_number(), copy.get_serial_number()),
        (7, 12)
    );
    assert_eq!(
        (
            copy.get_alt_loc().unwrap(),
            copy.get_insertion_code().unwrap()
        ),
        ("B".into(), "C".into())
    );
    assert_eq!((copy.get_occupancy(), copy.get_temp_factor()), (0.5, 20.0));
    assert!(copy.get_is_hetero_atom());
    assert_eq!(
        (copy.get_secondary_structure(), copy.get_segment_number()),
        (2, 3)
    );
    let p = AtomPDBResidueInfoParams {
        residue_name: "GLY".into(),
        chain_id: "Z".into(),
        ..Default::default()
    };
    let made = AtomPDBResidueInfo::new_with_params(" N  ", &p);
    assert_eq!(
        (
            made.get_residue_name().unwrap(),
            made.get_chain_id().unwrap()
        ),
        ("GLY".into(), "Z".into())
    );
    assert_eq!(AtomPDBResidueInfo::default().get_name().unwrap(), "");
}

#[test]
fn atoms_copy_info_in_and_out() {
    let mut mol = mol_from_smiles("CO").unwrap();
    assert!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_monomer_info()
            .is_none()
    );
    assert!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_pdb_residue_info()
            .unwrap()
            .is_none()
    );

    let mut residue = AtomPDBResidueInfo::new(" C1 ");
    residue.set_residue_name("LIG");
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_pdb_residue_info(&residue);
    let atom = mol.get_atom_with_idx(0).unwrap();
    let back = atom.get_pdb_residue_info().unwrap().unwrap();
    assert_eq!(back.get_residue_name().unwrap(), "LIG");
    let generic = atom.get_monomer_info().unwrap();
    assert_eq!(generic.get_monomer_type(), AtomMonomerType::Pdbresidue);
    assert_eq!(generic.get_name().unwrap(), " C1 ");

    // Editing a copy leaves the atom alone until the copy is set back.
    let mut edited = back.clone();
    edited.set_residue_name("ALT");
    assert_eq!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_pdb_residue_info()
            .unwrap()
            .unwrap()
            .get_residue_name()
            .unwrap(),
        "LIG"
    );

    let mut atom = Atom::from_num(8);
    atom.set_monomer_info(&AtomMonomerInfo::new(AtomMonomerType::Other).unwrap());
    assert!(matches!(
        atom.get_pdb_residue_info(),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        atom.get_monomer_info().unwrap().get_monomer_type(),
        AtomMonomerType::Other
    );
    let mut rw = RWMol::default();
    rw.add_atom(&atom);
    assert_eq!(
        rw.get_atom_with_idx(0)
            .unwrap()
            .get_monomer_info()
            .unwrap()
            .get_monomer_type(),
        AtomMonomerType::Other
    );
}

fn is_runtime<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(rdkit_cxx::Error::Runtime { .. }))
}

fn is_value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

// RDKit casts any info of type PDBRESIDUE to AtomPDBResidueInfo, which reads
// past the end of a plain AtomMonomerInfo (CONVENTIONS.md U10). Python allows
// one.
#[test]
fn plain_info_cannot_have_the_pdb_type() {
    assert!(is_value(AtomMonomerInfo::new(AtomMonomerType::Pdbresidue)));
    assert!(is_value(AtomMonomerInfo::new_with_params(
        AtomMonomerType::Pdbresidue,
        &AtomMonomerInfoParams::default()
    )));
    let mut info = AtomMonomerInfo::new(AtomMonomerType::Other).unwrap();
    assert!(is_value(info.set_monomer_type(AtomMonomerType::Pdbresidue)));
    assert_eq!(info.get_monomer_type(), AtomMonomerType::Other);
    let mut info = AtomMonomerInfo::default();
    assert!(is_value(info.set_monomer_type(AtomMonomerType::Pdbresidue)));
    assert_eq!(info.get_monomer_type(), AtomMonomerType::Unknown);
}

#[test]
fn pdb_residue_info_changes_type_freely() {
    let mut residue = AtomPDBResidueInfo::new(" CA ");
    residue.set_residue_name("ALA");
    residue.set_monomer_type(AtomMonomerType::Other);
    assert_eq!(residue.get_monomer_type(), AtomMonomerType::Other);
    let mut atom = Atom::from_num(6);
    atom.set_pdb_residue_info(&residue);
    // As in Python, PDB residue info of another type is not read as such.
    assert!(is_value(atom.get_pdb_residue_info()));

    // A plain copy of PDB residue info is still one, so it may have its type
    // back, and its clones too.
    let mut copy = atom.get_monomer_info().unwrap();
    copy.set_monomer_type(AtomMonomerType::Unknown).unwrap();
    let mut clone = copy.clone();
    clone.set_monomer_type(AtomMonomerType::Pdbresidue).unwrap();
    atom.set_monomer_info(&clone);
    let back = atom.get_pdb_residue_info().unwrap().unwrap();
    assert_eq!(
        (back.get_name().unwrap(), back.get_residue_name().unwrap()),
        (" CA ".into(), "ALA".into())
    );
    residue.set_monomer_type(AtomMonomerType::Pdbresidue);
    assert_eq!(residue.get_monomer_type(), AtomMonomerType::Pdbresidue);
}

/// A pickle of methane whose atom has plain info named "MARK", and where its
/// type byte is.
fn marked_pickle() -> (Vec<u8>, usize) {
    let mut mol = mol_from_smiles("C").unwrap();
    let params = AtomMonomerInfoParams {
        name: "MARK".into(),
        ..Default::default()
    };
    mol.get_atom_with_idx_mut(0).unwrap().set_monomer_info(
        &AtomMonomerInfo::new_with_params(AtomMonomerType::Other, &params).unwrap(),
    );
    let pkl = mol.to_binary();
    // The info's name (a u32 length, then its bytes), then its type as a byte.
    let name = pkl
        .windows(9)
        .position(|w| w == b"\x04\0\0\0MARK\x02")
        .unwrap();
    (pkl, name + 8)
}

#[test]
fn pickles_check_monomer_types() {
    let (pkl, at) = marked_pickle();
    let typ = |pkl: &[u8]| {
        // SAFETY: a pickle RDKit wrote.
        unsafe { Mol::from_pkl_string(pkl) }.map(|mol| {
            mol.get_atom_with_idx(0)
                .unwrap()
                .get_monomer_info()
                .unwrap()
                .get_monomer_type()
        })
    };
    assert_eq!(typ(&pkl).unwrap(), AtomMonomerType::Other);
    let mut unknown = pkl.clone();
    unknown[at] = 0;
    assert_eq!(typ(&unknown).unwrap(), AtomMonomerType::Unknown);
    // RDKit 2025.03 reads the PDBRESIDUE type as PDB residue info, from the
    // fields that follow it.
    let mut pdb = pkl.clone();
    pdb[at] = 1;
    assert_eq!(typ(&pdb).unwrap(), AtomMonomerType::Pdbresidue);
    // Other types are RDKit's "unrecognized MonomerType", a RuntimeError in
    // Python.
    for bad in [3, 255] {
        let mut changed = pkl.clone();
        changed[at] = bad;
        assert!(is_runtime(typ(&changed)), "type {bad}");
        // SAFETY: a pickle RDKit wrote but for its monomer type byte, which RDKit checks.
        assert!(is_runtime(unsafe {
            Mol::from_pkl_string_with_property_flags(&changed, PropertyPickleOptions::AllProps)
        }));
        // SAFETY: as above.
        assert!(is_runtime(unsafe { RWMol::from_pkl_string(&changed) }));
    }

    // MolBundle pickles hold molecule pickles, byte for byte.
    let mut bundle = MolBundle::default();
    // SAFETY: a pickle RDKit wrote.
    bundle.add_mol(&unsafe { Mol::from_pkl_string(&pkl) }.unwrap());
    let bundle_pkl = bundle.to_binary().unwrap();
    let name = bundle_pkl
        .windows(9)
        .position(|w| w == b"\x04\0\0\0MARK\x02")
        .unwrap();
    // SAFETY: a pickle RDKit wrote.
    assert_eq!(
        unsafe { MolBundle::from_pkl_string(&bundle_pkl) }
            .unwrap()
            .size(),
        1
    );
    let mut changed = bundle_pkl;
    changed[name + 8] = 3;
    // SAFETY: a bundle pickle RDKit wrote but for a monomer type byte, which RDKit checks.
    assert!(is_runtime(unsafe { MolBundle::from_pkl_string(&changed) }));
}

// MolPickler pickles a recursive query's molecule as a molecule, monomer
// infos included, and reads it back the same way.
#[test]
fn recursive_query_molecules_are_checked() {
    let (pkl, _) = marked_pickle();
    // SAFETY: a pickle RDKit wrote.
    let marked = unsafe { Mol::from_pkl_string(&pkl) }.unwrap();
    let mut query = mol_from_smarts("[#6]").unwrap();
    add_recursive_query(&mut query, &marked, 0).unwrap();
    // Twice nested, under an AND with the atom's own query.
    let mut outer = mol_from_smarts("[#6]").unwrap();
    add_recursive_query(&mut outer, &query, 0).unwrap();
    let mut pkl = outer.to_binary();
    let at = pkl
        .windows(9)
        .position(|w| w == b"\x04\0\0\0MARK\x02")
        .unwrap()
        + 8;
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { Mol::from_pkl_string(&pkl) }.is_ok());
    pkl[at] = 1;
    // SAFETY: a pickle RDKit wrote but for a monomer type byte, which RDKit
    // 2025.03 reads as PDB residue info.
    assert!(unsafe { Mol::from_pkl_string(&pkl) }.is_ok());
    pkl[at] = 9;
    // SAFETY: a pickle RDKit wrote but for a monomer type byte, which RDKit checks.
    assert!(is_runtime(unsafe { Mol::from_pkl_string(&pkl) }));
}

// RDKit pickles PDB residue info by its type, so a real one round-trips.
#[test]
fn pickles_keep_pdb_residue_info() {
    let mut mol = mol_from_smiles("C").unwrap();
    let mut residue = AtomPDBResidueInfo::new(" C1 ");
    residue.set_residue_name("LIG");
    residue.set_serial_number(7);
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_pdb_residue_info(&residue);
    // SAFETY: a pickle RDKit wrote.
    let back = unsafe { Mol::from_pkl_string(&mol.to_binary()) }.unwrap();
    let info = back
        .get_atom_with_idx(0)
        .unwrap()
        .get_pdb_residue_info()
        .unwrap()
        .unwrap();
    assert_eq!(
        (info.get_residue_name().unwrap(), info.get_serial_number()),
        ("LIG".into(), 7)
    );
}
