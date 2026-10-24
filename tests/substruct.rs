use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn matches() {
    let catechol = mol_from_smiles("Oc1ccccc1O").unwrap();
    let hydroxyl = mol_from_smarts("[OX2H]").unwrap();
    assert!(catechol.has_substruct_match(&hydroxyl).unwrap());
    assert_eq!(catechol.get_substruct_match(&hydroxyl).unwrap(), [0]);
    assert_eq!(
        catechol.get_substruct_matches(&hydroxyl).unwrap(),
        [vec![0], vec![7]]
    );
    let amine = mol_from_smarts("N").unwrap();
    assert!(!catechol.has_substruct_match(&amine).unwrap());
    assert!(catechol.get_substruct_match(&amine).unwrap().is_empty());
    assert!(catechol.get_substruct_matches(&amine).unwrap().is_empty());
    let one = SubstructMatchParameters {
        max_matches: 1,
        ..Default::default()
    };
    assert_eq!(
        catechol
            .get_substruct_matches_with_params(&hydroxyl, &one)
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn matches_are_ordered_by_query_atom() {
    let mol = mol_from_smiles("CCO").unwrap();
    let query = mol_from_smarts("OC").unwrap();
    assert_eq!(mol.get_substruct_match(&query).unwrap(), [2, 1]);
    let all = SubstructMatchParameters {
        uniquify: false,
        ..Default::default()
    };
    let ethane = mol_from_smiles("CC").unwrap();
    let cc = mol_from_smarts("CC").unwrap();
    assert_eq!(ethane.get_substruct_matches(&cc).unwrap().len(), 1);
    assert_eq!(
        ethane
            .get_substruct_matches_with_params(&cc, &all)
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn chirality_is_opt_in() {
    let target = mol_from_smiles("C[C@H](F)Cl").unwrap();
    let query = mol_from_smiles("C[C@@H](F)Cl").unwrap();
    assert!(target.has_substruct_match(&query).unwrap());
    let chiral = SubstructMatchParameters {
        use_chirality: true,
        ..Default::default()
    };
    assert!(
        !target
            .has_substruct_match_with_params(&query, &chiral)
            .unwrap()
    );
    assert!(
        target
            .get_substruct_match_with_params(&query, &chiral)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn defaults_come_from_rdkit() {
    let p = SubstructMatchParameters::default();
    assert!(p.uniquify && p.recursion_possible);
    assert!(!(p.use_chirality || p.use_query_query_matches || p.use_enhanced_stereo));
    assert_eq!(
        (p.max_matches, p.max_recursive_matches, p.num_threads),
        (1000, 1000, 1)
    );
    assert!(p.atom_properties.is_empty() && p.bond_properties.is_empty());
}

#[test]
fn property_matching() {
    let mut target = mol_from_smiles("CC").unwrap();
    target
        .get_atom_with_idx_mut(1)
        .unwrap()
        .set_prop("tag", "x");
    let mut query = mol_from_smiles("C").unwrap();
    query.get_atom_with_idx_mut(0).unwrap().set_prop("tag", "x");
    let params = SubstructMatchParameters {
        atom_properties: vec!["tag".into()],
        ..Default::default()
    };
    assert_eq!(
        target
            .get_substruct_matches_with_params(&query, &params)
            .unwrap(),
        [vec![1]]
    );
}

#[test]
fn ring_query_on_an_unsanitized_target() {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let target = mol_from_smiles_with_params("C1CC1", &params).unwrap();
    assert!(
        target
            .has_substruct_match(&mol_from_smarts("[R]").unwrap())
            .unwrap()
    );
}

#[test]
fn ring_query_without_rings_is_an_error_not_an_abort() {
    let mut rw = RWMol::default();
    for _ in 0..3 {
        rw.add_atom(&Atom::from_num(6));
    }
    for (begin, end) in [(0, 1), (1, 2), (2, 0)] {
        rw.add_bond_with_order(begin, end, BondType::Single)
            .unwrap();
    }
    let ring = mol_from_smarts("[R]").unwrap();
    assert!(matches!(rw.has_substruct_match(&ring),
                     Err(Error::Runtime { message }) if message.contains("RingInfo not initialized")));
    // Matching without ring queries works on the same molecule.
    assert!(
        rw.has_substruct_match(&mol_from_smarts("C").unwrap())
            .unwrap()
    );
}
