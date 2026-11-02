// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mmpa::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdMMPA.

fn smiles_pairs(pairs: &[(Option<Mol>, Mol)]) -> Vec<(Option<String>, String)> {
    pairs
        .iter()
        .map(|(core, side)| {
            (
                core.as_ref().map(|core| mol_to_smiles(core).unwrap()),
                mol_to_smiles(side).unwrap(),
            )
        })
        .collect()
}

fn some(core: &str, side: &str) -> (Option<String>, String) {
    (Some(core.to_string()), side.to_string())
}

#[test]
fn fragmentations() {
    let mol = mol_from_smiles("CC(C)Oc1ccccc1N").unwrap();
    let all = fragment_mol(&mol).unwrap();
    assert_eq!(all.len(), 14);
    assert_eq!(
        smiles_pairs(&all[..4]),
        [
            (None, "CC(Oc1ccccc1N)[*:1].C[*:1]".to_string()),
            some("Nc1ccccc1OC([*:1])[*:2]", "C[*:1].C[*:2]"),
            some("C([*:1])([*:2])[*:3]", "C[*:1].C[*:2].Nc1ccccc1O[*:3]"),
            some("O(C([*:1])[*:2])[*:3]", "C[*:1].C[*:2].Nc1ccccc1[*:3]"),
        ]
    );
    assert_eq!(
        fragment_mol_as_smiles(&mol).unwrap()[..2],
        [
            (String::new(), "CC(Oc1ccccc1N)[*:1].C[*:1]".to_string()),
            (
                "Nc1ccccc1OC([*:1])[*:2]".to_string(),
                "C[*:1].C[*:2]".to_string()
            ),
        ]
    );
    let two = fragment_mol_min_cuts(&mol, 2, 2, 20).unwrap();
    assert_eq!(two.len(), 7);
    assert_eq!(
        smiles_pairs(&two[..2]),
        [
            some("Nc1ccccc1OC([*:1])[*:2]", "C[*:1].C[*:2]"),
            some("CC([*:1])[*:2]", "C[*:1].Nc1ccccc1O[*:2]"),
        ]
    );
    assert_eq!(
        smiles_pairs(&fragment_mol_bonds(&mol, &[3, 9]).unwrap()),
        [
            (None, "CC(C)O[*:1].Nc1ccccc1[*:1]".to_string()),
            some("c1ccc([*:2])c([*:1])c1", "CC(C)O[*:2].N[*:1]"),
            (None, "CC(C)Oc1ccccc1[*:1].N[*:1]".to_string()),
        ]
    );
    // No match, or more matched bonds than maxCutBonds, gives nothing.
    assert!(
        fragment_mol(&mol_from_smiles("C").unwrap())
            .unwrap()
            .is_empty()
    );
    let few = FragmentMolParams {
        max_cut_bonds: 1,
        ..Default::default()
    };
    assert!(fragment_mol_with_params(&mol, &few).unwrap().is_empty());
}

#[test]
fn arguments_that_rdkit_trusts() {
    let mol = mol_from_smiles("CC(C)Oc1ccccc1N").unwrap();
    // RDKit dereferences the parsed pattern and reads two atoms of each match.
    let bad = FragmentMolParams {
        pattern: "[C".to_string(),
        ..Default::default()
    };
    assert!(matches!(
        fragment_mol_with_params(&mol, &bad),
        Err(Error::Parse { .. })
    ));
    let one_atom = FragmentMolParams {
        pattern: "[#6]".to_string(),
        ..Default::default()
    };
    assert!(matches!(
        fragment_mol_with_params(&mol, &one_atom),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        fragment_mol_bonds(&mol, &[]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        fragment_mol_bonds(&mol, &[mol.get_num_bonds()]),
        Err(Error::Index { .. })
    ));
    assert!(matches!(
        fragment_mol_min_cuts(&mol, 0, 2, 20),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        fragment_mol_min_cuts(&mol, 3, 2, 20),
        Err(Error::Value { .. })
    ));
}

#[test]
fn params_defaults() {
    assert_eq!(FragmentMolParams::default().max_cuts, 3);
    assert_eq!(FragmentMolParams::default().max_cut_bonds, 20);
    assert_eq!(
        FragmentMolParams::default().pattern,
        "[#6+0;!$(*=,#[!#6])]!@!=!#[*]"
    );
    assert_eq!(
        FragmentMolBondsParams::default(),
        FragmentMolBondsParams {
            min_cuts: 1,
            max_cuts: 3,
        }
    );
}
