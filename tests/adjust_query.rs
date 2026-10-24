// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::chem::*;

fn smiles(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

#[test]
fn adjust_queries() {
    let query = smiles("c1ccccc1*");
    let adjusted = adjust_query_properties(&query).unwrap();
    assert_eq!(
        mol_to_smarts(&adjusted).unwrap(),
        "[#6&D2]1:[#6&D2]:[#6&D2]:[#6&D2]:[#6&D2]:[#6&D3]:1-*"
    );
    let mut none = AdjustQueryParameters::no_adjustments();
    assert_eq!(
        mol_to_smarts(&adjust_query_properties_with_params(&query, &none).unwrap()).unwrap(),
        "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#0]"
    );
    none.make_atoms_generic = true;
    assert_eq!(
        mol_to_smarts(&adjust_query_properties_with_params(&query, &none).unwrap()).unwrap(),
        "*1:*:*:*:*:*:1-*"
    );
}

#[test]
fn parameter_defaults() {
    let defaults = AdjustQueryParameters::default();
    assert!(
        defaults.adjust_degree && defaults.make_dummies_queries && defaults.aromatize_if_possible
    );
    assert!(
        !defaults.adjust_ring_count && !defaults.make_atoms_generic && !defaults.make_bonds_generic
    );
    let ignore =
        AdjustQueryWhichFlags::ADJUST_IGNOREDUMMIES | AdjustQueryWhichFlags::ADJUST_IGNORECHAINS;
    assert_eq!(ignore.bits(), 3);
    assert_eq!(defaults.adjust_degree_flags, ignore);
    assert_eq!(
        defaults.make_atoms_generic_flags,
        AdjustQueryWhichFlags::ADJUST_IGNORENONE
    );
    let none = AdjustQueryParameters::no_adjustments();
    assert!(!none.adjust_degree && !none.make_dummies_queries && !none.aromatize_if_possible);
    assert_eq!(
        format!("{:?}", AdjustQueryWhichFlags::ADJUST_IGNOREALL),
        "AdjustQueryWhichFlags::ADJUST_IGNOREALL"
    );
    let generic = SetGenericQueriesFromPropertiesParams::default();
    assert!(generic.use_atom_labels && generic.use_s_groups);
}

#[test]
fn generic_groups() {
    let mut mol = smiles("*c1ccccc1 |$ARY_p;;;;;;$|");
    set_generic_queries_from_properties(&mut mol).unwrap();
    assert_eq!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_prop("_QueryAtomGenericLabel")
            .unwrap(),
        "ARY"
    );
    assert_eq!(
        adjust_query_properties_with_generic_groups(&mol)
            .unwrap()
            .get_num_atoms(),
        7
    );
    convert_generic_queries_to_substance_groups(&mut mol).unwrap();
    assert_eq!(get_mol_substance_groups(&mol).len(), 1);
}
