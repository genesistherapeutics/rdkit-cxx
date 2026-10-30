// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// What RDKit 2024.03's MolPickler writes but cannot read back: the counts it
// writes in one byte, and the property counts it writes in two
// (check_picklable in mol_checks.cpp).

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chem_reactions::ChemicalReaction;
use rdkit_cxx::chem::rd_substruct_library::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

/// Whether `to_binary` refuses the molecule (a panic naming `what`), or else
/// reads its pickle back.
fn refused_for(mol: &Mol, what: &str) -> bool {
    match std::panic::catch_unwind(|| mol.to_binary()) {
        Err(payload) => {
            let message = payload.downcast_ref::<String>().unwrap();
            assert!(message.contains(what), "{message}");
            true
        }
        Ok(pickle) => {
            // SAFETY: a pickle RDKit wrote.
            let back = unsafe { Mol::from_pkl_string(&pickle) }.unwrap();
            assert_eq!(back.get_num_atoms(), mol.get_num_atoms());
            false
        }
    }
}

/// `smiles` with a substance group of atom 0 and bond 0 that `edit` fills
/// `n` times.
fn with_substance_group(smiles: &str, n: usize, edit: impl Fn(&mut SubstanceGroupMut<'_>)) -> Mol {
    let mut mol = mol_from_smiles(smiles).unwrap();
    let mut sgroup = create_mol_substance_group(&mut mol, "GEN");
    sgroup.add_atom_with_idx(0).unwrap();
    sgroup.add_bond_with_idx(0).unwrap();
    for _ in 0..n {
        edit(&mut sgroup);
    }
    mol
}

/// A molecule whose first atom has a recursive query with `n` children.
fn with_query_children(smiles: &str, n: usize) -> Mol {
    let mut mol = mol_from_smiles(smiles).unwrap();
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("query", &vec!["a"; n].join(","));
    let queries = BTreeMap::from([("a".to_owned(), mol_from_smiles("C").unwrap())]);
    mol_add_recursive_queries(&mut mol, &queries, "query").unwrap();
    mol
}

// MolPickler writes these counts in one byte (in a molecule of at most 255
// atoms, whatever its bond count, in RDKit 2024.03; query children in any
// molecule), truncating larger ones, and then cannot read its own pickle
// back.
#[test]
fn one_byte_counts_past_255_are_refused() {
    let bracket = [Point3D::default(), Point3D::new(1.0, 0.0, 0.0)];
    for n in [255, 256] {
        let refused = n > 255;
        let mols = [
            (
                with_substance_group("CCO", n, |s| s.add_bracket(&bracket).unwrap()),
                "number of brackets",
            ),
            (
                with_substance_group("CCO", n, |s| s.add_attach_point(0, 1, "1").unwrap()),
                "number of attachment points",
            ),
            (
                with_substance_group("CCO", n - 1, |s| s.add_atom_with_idx(1).unwrap()),
                "number of atoms",
            ),
            (
                with_substance_group("CCO", n, |s| s.add_parent_atom_with_idx(0).unwrap()),
                "number of parent atoms",
            ),
            (
                with_substance_group("CCO", n - 1, |s| s.add_bond_with_idx(1).unwrap()),
                "number of bonds",
            ),
            (
                with_substance_group("CCO", n, |s| s.add_c_state(0, Point3D::default()).unwrap()),
                "number of CStates",
            ),
            (with_query_children("CC", n), "number of children"),
        ];
        for (mol, what) in &mols {
            assert_eq!(refused_for(mol, what), refused, "{what} {n}");
        }

        let mut rw = RWMol::from(mol_from_smiles("C[C@H](F)Cl").unwrap());
        let groups: Vec<_> = (0..n)
            .map(|_| create_stereo_group(StereoGroupType::StereoOr, &rw, &[1]).unwrap())
            .collect();
        rw.set_stereo_groups(&groups).unwrap();
        assert_eq!(
            refused_for(&Mol::from(rw), "number of stereo groups"),
            refused,
            "stereo groups {n}"
        );
    }

    // A molecule of more than 255 atoms is pickled with 32-bit counts, but
    // query children are counted in one byte whatever its size, and so are
    // a recursive query's molecule's counts, by its own size.
    let large = "C".repeat(300);
    let many_brackets = with_substance_group(&large, 300, |s| s.add_bracket(&bracket).unwrap());
    assert!(!refused_for(&many_brackets, ""));
    assert!(refused_for(
        &with_query_children(&large, 256),
        "number of children"
    ));
    let mut mol = mol_from_smiles(&large).unwrap();
    mol.get_atom_with_idx_mut(0).unwrap().set_prop("query", "a");
    let inner = with_substance_group("CCO", 256, |s| s.add_bracket(&bracket).unwrap());
    mol_add_recursive_queries(
        &mut mol,
        &BTreeMap::from([("a".to_owned(), inner)]),
        "query",
    )
    .unwrap();
    assert!(refused_for(
        &mol,
        "a recursive query's molecule's substance group"
    ));

    // Formats that hold molecule pickles refuse them too.
    let bad = with_substance_group("CCO", 256, |s| s.add_bracket(&bracket).unwrap());
    let good = with_substance_group("CCO", 255, |s| s.add_bracket(&bracket).unwrap());
    let mut holder = CachedMolHolder::default();
    assert!(matches!(holder.add_mol(&bad), Err(Error::Value { .. })));
    assert_eq!(holder.add_mol(&good).unwrap(), 0);
    assert_eq!(holder.len(), 1);
    assert_eq!(holder.get_mol(0).unwrap().get_num_atoms(), 3);
    let mut library =
        SubstructLibrary::new(MolHolderBase::from(CachedMolHolder::default()), None, None);
    assert!(matches!(library.add_mol(&bad), Err(Error::Value { .. })));
    assert_eq!(library.add_mol(&good).unwrap(), 0);
    let mut bundle = MolBundle::default();
    bundle.add_mol(&good);
    assert!(bundle.to_binary().is_ok());
    bundle.add_mol(&bad);
    assert!(matches!(bundle.to_binary(), Err(Error::Value { .. })));
    let mut rxn = ChemicalReaction::default();
    rxn.add_reactant_template(&good);
    rxn.to_binary();
    rxn.add_product_template(&with_query_children("CC", 256));
    assert!(std::panic::catch_unwind(|| rxn.to_binary()).is_err());
}

// MolPickler counts each object's properties in a u16: past 65535 the
// count wraps and Python's ToBinary writes a pickle it cannot read back
// (2024.03.6: "Bad pickle format: ENDPROPS tag not found" for a molecule).
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
fn more_than_65535_atom_properties_are_refused() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    {
        let mut atom = mol.get_atom_with_idx_mut(1).unwrap();
        for i in 0..65536 {
            atom.set_prop(&format!("p{i}"), "x");
        }
    }
    let mut holder = CachedMolHolder::default();
    match holder.add_mol(&mol) {
        Err(Error::Value { message }) => {
            assert!(message.contains("atom 1 has"), "{message}")
        }
        other => panic!("{other:?}"),
    }
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
