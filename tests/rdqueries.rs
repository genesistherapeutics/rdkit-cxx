// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rdqueries::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::data_structs::ExplicitBitVect;

fn mol() -> Mol {
    mol_from_smiles("C[13CH2]C(=O)[O-]").unwrap()
}

fn hits(mol: &Mol, query: &Atom) -> Vec<u32> {
    mol.get_atoms_matching_query(query)
        .unwrap()
        .iter()
        .map(|a| a.get_idx())
        .collect()
}

// Expected values come from Python's rdqueries.

#[test]
fn value_queries() {
    let m = mol();
    assert_eq!(
        hits(&m, &atom_num_equals_query_atom(6).unwrap()),
        vec![0, 1, 2]
    );
    assert_eq!(
        hits(
            &m,
            &atom_num_equals_query_atom_with_negate(6, true).unwrap()
        ),
        vec![3, 4]
    );
    // As in Python, "less" matches atoms whose value is less than val.
    assert_eq!(
        hits(&m, &atom_num_less_query_atom(7).unwrap()),
        vec![0, 1, 2]
    );
    assert_eq!(
        hits(&m, &atom_num_greater_query_atom(6).unwrap()),
        vec![3, 4]
    );
    assert_eq!(
        atom_num_less_query_atom(7)
            .unwrap()
            .describe_query()
            .unwrap(),
        "AtomNumLess 7 > \n"
    );
    assert_eq!(hits(&m, &isotope_equals_query_atom(13).unwrap()), vec![1]);
    assert!(hits(&m, &mass_equals_query_atom(12).unwrap()).is_empty());
    assert_eq!(
        hits(&m, &formal_charge_equals_query_atom(-1).unwrap()),
        vec![4]
    );
    assert_eq!(
        hits(&m, &formal_charge_less_query_atom(0).unwrap()),
        vec![4]
    );
    assert_eq!(hits(&m, &h_count_equals_query_atom(3).unwrap()), vec![0]);
    assert_eq!(
        hits(&m, &h_count_greater_query_atom(1).unwrap()),
        vec![0, 1]
    );
    assert_eq!(
        hits(&m, &min_ring_size_equals_query_atom(0).unwrap()),
        vec![0, 1, 2, 3, 4]
    );
}

#[test]
fn value_queries_refuse_values_that_overflow() {
    // RDKit subtracts the atom's value from the query's in int.
    let limit = 1 << 30;
    assert!(atom_num_equals_query_atom(limit).is_ok());
    assert!(atom_num_equals_query_atom(-limit).is_ok());
    for val in [limit + 1, -limit - 1, i32::MAX, i32::MIN] {
        assert!(matches!(
            formal_charge_less_query_atom(val),
            Err(Error::Value { .. })
        ));
    }
    // Masses are scaled by 1000 first.
    assert!(mass_equals_query_atom(limit / 1000).is_ok());
    assert!(matches!(
        mass_greater_query_atom_with_negate(limit / 1000 + 1, true),
        Err(Error::Value { .. })
    ));
}

#[test]
fn flag_queries() {
    let m = mol();
    assert_eq!(hits(&m, &is_unsaturated_query_atom()), vec![2, 3]);
    assert!(hits(&m, &x_atom_query_atom()).is_empty());
    assert_eq!(hits(&m, &q_atom_query_atom()), vec![3, 4]);
    assert_eq!(hits(&m, &a_atom_query_atom()), vec![0, 1, 2, 3, 4]);
    assert!(hits(&m, &is_in_ring_query_atom()).is_empty());
    assert_eq!(
        hits(&m, &is_in_ring_query_atom_with_negate(true)),
        vec![0, 1, 2, 3, 4]
    );
    let bicycle = mol_from_smiles("C1CC2CCC1C2").unwrap();
    assert_eq!(hits(&bicycle, &is_bridgehead_query_atom()), vec![2, 5]);
}

#[test]
fn property_queries() {
    let mut m = mol();
    m.get_atom_with_idx_mut(0).unwrap().set_int_prop("n", 4);
    m.get_atom_with_idx_mut(1).unwrap().set_prop("s", "x");
    m.get_atom_with_idx_mut(2)
        .unwrap()
        .set_double_prop("d", 1.5);
    m.get_atom_with_idx_mut(3).unwrap().set_bool_prop("b", true);
    let mut bv = ExplicitBitVect::new(8);
    bv.set_bit(1).unwrap();
    m.get_atom_with_idx_mut(4)
        .unwrap()
        .set_explicit_bit_vect_prop("bv", &bv);

    let near = HasIntPropWithValueQueryAtomParams {
        tolerance: 1,
        ..Default::default()
    };
    assert_eq!(
        hits(
            &m,
            &has_int_prop_with_value_query_atom_with_params("n", 3, &near)
        ),
        vec![0]
    );
    assert!(hits(&m, &has_int_prop_with_value_query_atom("n", 3)).is_empty());
    assert_eq!(
        hits(&m, &has_string_prop_with_value_query_atom("s", "x")),
        vec![1]
    );
    assert_eq!(hits(&m, &has_prop_query_atom("s")), vec![1]);
    assert_eq!(
        hits(&m, &has_prop_query_atom_with_negate("s", true)),
        vec![0, 2, 3, 4]
    );
    let close = HasDoublePropWithValueQueryAtomParams {
        tolerance: 0.6,
        ..Default::default()
    };
    assert_eq!(
        hits(
            &m,
            &has_double_prop_with_value_query_atom_with_params("d", 1.0, &close)
        ),
        vec![2]
    );
    assert!(hits(&m, &has_double_prop_with_value_query_atom("d", 1.0)).is_empty());
    assert_eq!(
        hits(&m, &has_bool_prop_with_value_query_atom("b", true)),
        vec![3]
    );
    assert!(hits(&m, &has_bool_prop_with_value_query_atom("b", false)).is_empty());
    assert_eq!(
        hits(&m, &has_bit_vect_prop_with_value_query_atom("bv", &bv)),
        vec![4]
    );
    assert_eq!(
        has_int_prop_with_value_query_atom("n", 3)
            .describe_query()
            .unwrap(),
        "HasPropWithValue 0 = val\n"
    );
}

#[test]
fn int_property_queries_at_the_limits() {
    // RDKit's queryCmp would subtract these in int; the wrapper's query
    // compares in 64 bits.
    let mut m = mol();
    m.get_atom_with_idx_mut(0)
        .unwrap()
        .set_int_prop("n", i32::MIN);
    m.get_atom_with_idx_mut(1)
        .unwrap()
        .set_int_prop("n", i32::MAX);
    assert!(hits(&m, &has_int_prop_with_value_query_atom("n", 1)).is_empty());
    assert_eq!(
        hits(&m, &has_int_prop_with_value_query_atom("n", i32::MIN)),
        vec![0]
    );
    // i32::MIN is one further from 0 than i32::MAX reaches.
    let wide = HasIntPropWithValueQueryAtomParams {
        tolerance: i32::MAX,
        ..Default::default()
    };
    assert_eq!(
        hits(
            &m,
            &has_int_prop_with_value_query_atom_with_params("n", 0, &wide)
        ),
        vec![1]
    );
    assert_eq!(
        hits(
            &m,
            &has_int_prop_with_value_query_atom_with_params("n", -1, &wide)
        ),
        vec![0]
    );
    let none = HasIntPropWithValueQueryAtomParams {
        tolerance: i32::MIN,
        ..Default::default()
    };
    assert!(
        hits(
            &m,
            &has_int_prop_with_value_query_atom_with_params("n", 0, &none)
        )
        .is_empty()
    );
    // A copy keeps the 64-bit comparison.
    let copy = has_int_prop_with_value_query_atom("n", 1).clone();
    assert!(hits(&m, &copy).is_empty());
    let negated = HasIntPropWithValueQueryAtomParams {
        negate: true,
        ..Default::default()
    };
    assert_eq!(
        hits(
            &m,
            &has_int_prop_with_value_query_atom_with_params("n", i32::MAX, &negated)
        ),
        vec![0, 2, 3, 4]
    );
}

#[test]
fn bond_queries() {
    let mut m = RWMol::from(mol());
    m.get_bond_with_idx_mut(0).unwrap().set_prop("bp", "1");
    let query = has_prop_query_bond("bp");
    assert!(query.has_query());
    assert!(query.match_(&m.get_bond_with_idx(0).unwrap()).unwrap());
    assert!(!query.match_(&m.get_bond_with_idx(1).unwrap()).unwrap());
    assert!(
        has_prop_query_bond_with_negate("bp", true)
            .match_(&m.get_bond_with_idx(1).unwrap())
            .unwrap()
    );
    let int_query = has_int_prop_with_value_query_bond("i", 2);
    m.replace_bond(1, &int_query).unwrap();
    let replaced = m.get_bond_with_idx(1).unwrap();
    assert!(replaced.has_query());
    assert_eq!(
        replaced.describe_query().unwrap(),
        "HasPropWithValue 0 = val\n"
    );

    m.get_bond_with_idx_mut(2)
        .unwrap()
        .set_int_prop("i", i32::MIN);
    let far = has_int_prop_with_value_query_bond("i", 2);
    assert!(!far.match_(&m.get_bond_with_idx(2).unwrap()).unwrap());
    let params = HasIntPropWithValueQueryBondParams {
        tolerance: i32::MAX,
        ..Default::default()
    };
    let wide = has_int_prop_with_value_query_bond_with_params("i", -1, &params);
    assert!(wide.match_(&m.get_bond_with_idx(2).unwrap()).unwrap());

    m.get_bond_with_idx_mut(3).unwrap().set_bool_prop("b", true);
    assert!(
        has_bool_prop_with_value_query_bond("b", true)
            .match_(&m.get_bond_with_idx(3).unwrap())
            .unwrap()
    );
    m.get_bond_with_idx_mut(3).unwrap().set_prop("s", "v");
    assert!(
        has_string_prop_with_value_query_bond_with_negate("s", "v", true)
            .match_(&m.get_bond_with_idx(0).unwrap())
            .unwrap()
    );
    m.get_bond_with_idx_mut(3)
        .unwrap()
        .set_double_prop("d", 2.0);
    assert!(
        has_double_prop_with_value_query_bond("d", 2.0)
            .match_(&m.get_bond_with_idx(3).unwrap())
            .unwrap()
    );
}

#[test]
fn owned_bonds() {
    let mut query = has_prop_query_bond("a");
    let other = has_prop_query_bond("b");
    query.expand_query(&other).unwrap();
    let copy = query.clone();
    assert_eq!(
        copy.describe_query().unwrap(),
        query.describe_query().unwrap()
    );
    assert!(copy.has_query());
    assert!(!copy.has_owning_mol());
    query.set_int_prop("x", 3);
    assert_eq!(query.get_int_prop("x").unwrap(), 3);
    query.set_bond_type(BondType::Double);
    assert_eq!(query.get_bond_type(), BondType::Double);
    assert!(query.set_stereo_atoms(0, 1).is_err());
    // RDKit 2025.03 throws a bare C string for a query it cannot write, a
    // RuntimeError in Python.
    assert!(matches!(query.get_smarts(), Err(Error::Runtime { .. })));
    let mut m = RWMol::from(mol());
    m.get_bond_with_idx_mut(0)
        .unwrap()
        .set_query(&other)
        .unwrap_err();
}

#[test]
fn replace_with_query_atom() {
    let mut m = mol();
    let mut atom = replace_atom_with_query_atom(&mut m, 0).unwrap();
    assert!(atom.has_query());
    assert_eq!(atom.describe_query().unwrap(), "AtomAtomicNum 6 = val\n");
    atom.expand_query(&h_count_equals_query_atom(3).unwrap())
        .unwrap();
    assert_eq!(
        m.get_atom_with_idx(0).unwrap().get_smarts().unwrap(),
        "[#6&H3]"
    );
    // An atom with a query already is left alone.
    replace_atom_with_query_atom(&mut m, 0).unwrap();
    assert_eq!(
        m.get_atom_with_idx(0).unwrap().get_smarts().unwrap(),
        "[#6&H3]"
    );
    assert!(matches!(
        replace_atom_with_query_atom(&mut m, 5),
        Err(Error::Index { .. })
    ));
}

#[test]
fn query_atoms_in_molecules() {
    let mut m = RWMol::from(mol_from_smiles("C").unwrap());
    let idx = m.add_atom(&q_atom_query_atom());
    m.add_bond_with_order(0, idx, BondType::Single).unwrap();
    assert_eq!(
        mol_from_smiles("CO")
            .unwrap()
            .get_substruct_matches(&m)
            .unwrap(),
        vec![vec![0, 1]]
    );
    assert!(
        mol_from_smiles("CC")
            .unwrap()
            .get_substruct_matches(&m)
            .unwrap()
            .is_empty()
    );
}
