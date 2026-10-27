// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn smarts(text: &str) -> Mol {
    mol_from_smarts(text).unwrap()
}

#[test]
fn describe_and_match() {
    let query = smarts("[#6;R]-[O,N]");
    let atom = query.get_atom_with_idx(0).unwrap();
    assert_eq!(
        atom.describe_query().unwrap(),
        "AtomAnd\n  AtomAtomicNum 6 = val\n  AtomInNRings -1 = val\n"
    );
    assert_eq!(atom.get_query_type().unwrap(), "");
    assert_eq!(
        query
            .get_bond_with_idx(0)
            .unwrap()
            .describe_query()
            .unwrap(),
        "BondOrder 1 = val\n"
    );

    let mol = mol_from_smiles("C1CC1CO").unwrap();
    let plain = mol.get_atom_with_idx(0).unwrap();
    assert_eq!(plain.describe_query().unwrap(), "");
    assert!(atom.match_(&plain).unwrap());
    assert!(!atom.match_(&mol.get_atom_with_idx(3).unwrap()).unwrap());
    assert!(plain.match_(&mol.get_atom_with_idx(1).unwrap()).unwrap());
    // A ring query needs the other atom's molecule.
    assert!(matches!(
        atom.match_(&Atom::from_num(6)),
        Err(Error::Runtime { .. })
    ));

    let double = smarts("C=C");
    let bond = double.get_bond_with_idx(0).unwrap();
    let ethene = mol_from_smiles("C=C").unwrap();
    let ethane = mol_from_smiles("CC").unwrap();
    assert!(bond.match_(&ethene.get_bond_with_idx(0).unwrap()).unwrap());
    assert!(!bond.match_(&ethane.get_bond_with_idx(0).unwrap()).unwrap());
}

#[test]
fn atoms_matching_query() {
    let mol = mol_from_smiles("C1CC1CO").unwrap();
    let ring_carbon = smarts("[#6;R]");
    let atoms = mol
        .get_atoms_matching_query(&ring_carbon.get_atom_with_idx(0).unwrap())
        .unwrap();
    let indices: Vec<u32> = atoms.iter().map(|atom| atom.get_idx()).collect();
    assert_eq!(indices, [0, 1, 2]);

    // Recursive queries are matched first.
    let recursive = smarts("[$(CO)]");
    let atoms = mol
        .get_atoms_matching_query(&recursive.get_atom_with_idx(0).unwrap())
        .unwrap();
    assert_eq!(atoms.len(), 1);
    assert_eq!(atoms[0].get_idx(), 3);
    // Twice, since the recursive query's matches are cleared after each call.
    let other = mol_from_smiles("OCCC").unwrap();
    let atoms = other
        .get_atoms_matching_query(&recursive.get_atom_with_idx(0).unwrap())
        .unwrap();
    assert_eq!(atoms.len(), 1);
    assert_eq!(atoms[0].get_idx(), 1);

    assert!(matches!(
        mol.get_atoms_matching_query(&mol.get_atom_with_idx(0).unwrap()),
        Err(Error::Value { .. })
    ));
}

#[test]
fn expand_and_set_queries() {
    let nitrogen = smarts("[#7]");
    let other = nitrogen.get_atom_with_idx(0).unwrap();
    for (params, expected) in [
        (ExpandQueryParams::default(), "[#6&#7]"),
        (
            ExpandQueryParams {
                how: CompositeQueryType::CompositeOr,
                ..Default::default()
            },
            "[#6,#7]",
        ),
    ] {
        let mut query = smarts("[#6]");
        query
            .get_atom_with_idx_mut(0)
            .unwrap()
            .expand_query_with_params(&other, &params)
            .unwrap();
        assert_eq!(mol_to_smarts(&query).unwrap(), expected);
    }

    let mut xor = smarts("[#6]");
    let mut atom = xor.get_atom_with_idx_mut(0).unwrap();
    let params = ExpandQueryParams {
        how: CompositeQueryType::CompositeXor,
        maintain_order: false,
    };
    atom.expand_query_with_params(&other, &params).unwrap();
    assert_eq!(
        atom.describe_query().unwrap(),
        "AtomXor\n  AtomAtomicNum 7 = val\n  AtomAtomicNum 6 = val\n"
    );

    let mut set = smarts("[#6]");
    set.get_atom_with_idx_mut(0)
        .unwrap()
        .set_query(&other)
        .unwrap();
    assert_eq!(mol_to_smarts(&set).unwrap(), "[#7]");

    // Both atoms must be query atoms.
    let mut plain = mol_from_smiles("C").unwrap();
    assert!(matches!(
        plain.get_atom_with_idx_mut(0).unwrap().set_query(&other),
        Err(Error::Value { .. })
    ));
    let mut query = smarts("[#6]");
    assert!(matches!(
        query
            .get_atom_with_idx_mut(0)
            .unwrap()
            .expand_query(&Atom::from_num(7)),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        Atom::from_num(6).expand_query(&other),
        Err(Error::Value { .. })
    ));

    let double = smarts("C=C");
    let bond = double.get_bond_with_idx(0).unwrap();
    let mut single = smarts("C-C");
    let params = ExpandQueryParams {
        how: CompositeQueryType::CompositeOr,
        ..Default::default()
    };
    single
        .get_bond_with_idx_mut(0)
        .unwrap()
        .expand_query_with_params(&bond, &params)
        .unwrap();
    assert_eq!(mol_to_smarts(&single).unwrap(), "C-,=C");
    single
        .get_bond_with_idx_mut(0)
        .unwrap()
        .set_query(&bond)
        .unwrap();
    assert_eq!(mol_to_smarts(&single).unwrap(), "C=C");
    let ethane = mol_from_smiles("CC").unwrap();
    assert!(matches!(
        single
            .get_bond_with_idx_mut(0)
            .unwrap()
            .set_query(&ethane.get_bond_with_idx(0).unwrap()),
        Err(Error::Value { .. })
    ));
}
