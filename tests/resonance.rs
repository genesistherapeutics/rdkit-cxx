// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

const AMIDINIUM_BENZOATE: &str = "NC(=[NH2+])c1ccc(cc1)C(=O)[O-]";

fn supplier(smiles: &str) -> ResonanceMolSupplier {
    ResonanceMolSupplier::new(&mol_from_smiles(smiles).unwrap()).unwrap()
}

#[test]
fn structures() {
    let mut suppl = supplier("[O-]C=CC=O");
    assert!(!suppl.get_is_enumerated());
    assert_eq!(suppl.len(), 3);
    assert!(suppl.get_is_enumerated());
    assert_eq!(suppl.get_num_conj_grps(), 1);

    let mut smiles: Vec<String> = suppl
        .structures()
        .map(|mol| mol_to_smiles(&mol.unwrap()).unwrap())
        .collect();
    smiles.sort();
    assert_eq!(smiles, ["O=CC=C[O-]", "O=CC=C[O-]", "O=C[CH-]C=O"]);
    assert!(suppl.at_end());
    assert!(suppl.structures_from_here().next().is_none());
    // structures starts over, as Python's for loop does.
    assert_eq!(suppl.structures().count(), 3);

    suppl.reset();
    assert!(!suppl.at_end());
    let first = mol_to_smiles(&suppl.get(0).unwrap()).unwrap();
    let mut structures = suppl.structures_from_here();
    assert_eq!(
        mol_to_smiles(&structures.next().unwrap().unwrap()).unwrap(),
        first
    );
    assert_eq!(structures.count(), 2);
    assert_eq!(suppl.get(3).unwrap_err(), Error::Index { index: 3, len: 3 });
    assert!(!suppl.was_canceled());
}

#[test]
fn flags_and_groups() {
    let mol = mol_from_smiles(AMIDINIUM_BENZOATE).unwrap();
    let mut suppl = ResonanceMolSupplier::new(&mol).unwrap();
    suppl.set_num_threads(2);
    suppl.enumerate();
    assert_eq!(suppl.len(), 4);
    // RDKit negates a negative count in int, undefined for i32::MIN.
    let mut counted = ResonanceMolSupplier::new(&mol).unwrap();
    counted.set_num_threads(i32::MIN);
    counted.enumerate();
    assert_eq!(counted.len(), 4);
    assert_eq!(suppl.get_atom_conj_grp_idx(0).unwrap(), 0);
    assert_eq!(suppl.get_bond_conj_grp_idx(0).unwrap(), 0);
    assert_eq!(
        suppl.get_atom_conj_grp_idx(99).unwrap_err(),
        Error::Index { index: 99, len: 12 }
    );
    assert_eq!(
        suppl.get_bond_conj_grp_idx(12).unwrap_err(),
        Error::Index { index: 12, len: 12 }
    );

    let params = ResonanceMolSupplierParams {
        flags: ResonanceFlags::KEKULE_ALL,
        ..Default::default()
    };
    assert_eq!(
        ResonanceMolSupplier::new_with_params(&mol, &params)
            .unwrap()
            .len(),
        8
    );
    let params = ResonanceMolSupplierParams {
        flags: ResonanceFlags::ALLOW_CHARGE_SEPARATION | ResonanceFlags::ALLOW_INCOMPLETE_OCTETS,
        max_structs: 2,
    };
    assert_eq!(
        ResonanceMolSupplier::new_with_params(&mol, &params)
            .unwrap()
            .len(),
        2
    );

    assert_eq!(ResonanceFlags::ALLOW_INCOMPLETE_OCTETS.bits(), 1);
    assert_eq!(ResonanceFlags::UNCONSTRAINED_ANIONS.bits(), 16);
    assert_eq!(
        format!("{:?}", ResonanceFlags::KEKULE_ALL),
        "ResonanceFlags::KEKULE_ALL"
    );
    assert_eq!(ResonanceMolSupplierParams::default().max_structs, 1000);
}

#[test]
fn substructure_matches() {
    let mut suppl = supplier(AMIDINIUM_BENZOATE);
    let query = mol_from_smarts("NC(=N)").unwrap();
    // Across structures, the amidinium's two ends both match the double bond.
    assert_eq!(suppl.get_substruct_match(&query).unwrap(), [2, 1, 0]);
    assert_eq!(
        suppl.get_substruct_matches(&query).unwrap(),
        [vec![0, 1, 2], vec![2, 1, 0]]
    );
    let params = GetSubstructMatchesParams {
        uniquify: true,
        ..Default::default()
    };
    assert_eq!(
        suppl
            .get_substruct_matches_with_params(&query, &params)
            .unwrap(),
        [vec![0, 1, 2]]
    );
    assert!(
        suppl
            .get_substruct_match_with_params(&mol_from_smarts("S").unwrap(), &Default::default())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn unkekulizable_molecules() {
    let mol = mol_from_smiles_with_params(
        "c1cccc1",
        &SmilesParserParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(matches!(
        ResonanceMolSupplier::new(&mol),
        Err(Error::Kekulize { .. })
    ));
}
