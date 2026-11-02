// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_enumerator::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdMolEnumerator.

fn smiles(smi: &str) -> Mol {
    mol_from_smiles(smi).unwrap()
}

fn sorted(bundle: &MolBundle) -> Vec<String> {
    let mut smiles: Vec<String> = (0..bundle.size())
        .map(|i| mol_to_smiles(bundle.get_mol(i).unwrap()).unwrap())
        .collect();
    smiles.sort();
    smiles
}

fn with_link_nodes(smi: &str, nodes: &str) -> Mol {
    let mut mol = smiles(smi);
    mol.set_prop("_molLinkNodes", nodes);
    mol
}

fn assert_value(res: Result<MolBundle, Error>) {
    assert_value_with(res, "");
}

fn assert_value_with(res: Result<MolBundle, Error>, needle: &str) {
    match res.unwrap_err() {
        Error::Value { message } => assert!(message.contains(needle), "{message}"),
        err => panic!("{err:?}"),
    }
}

const PHENOL: &str = "*Cl.Oc1ccccc1 |m:0:3.4.5|";
const LINK: &str = "COC |LN:1:1.3|";
const SRU: &str = "FCN(CC)-* |Sg:n:2,5:2-3:ht|";

#[test]
fn params() {
    assert_eq!(
        MolEnumeratorParams::default(),
        MolEnumeratorParams {
            sanitize: false,
            max_to_enumerate: 1000,
            do_random: false,
            random_seed: -1,
            operation: None,
        }
    );
    let mut params = MolEnumeratorParams::new(EnumeratorType::LinkNode);
    assert_eq!(params.operation, Some(EnumeratorType::LinkNode));
    params.set_enumeration_operator(EnumeratorType::RepeatUnit);
    assert_eq!(params.operation, Some(EnumeratorType::RepeatUnit));
}

#[test]
fn position_variation() {
    let mol = smiles(PHENOL);
    let expected = ["Oc1(Cl)ccccc1", "Oc1cccc(Cl)c1", "Oc1ccccc1Cl"];
    assert_eq!(sorted(&enumerate(&mol).unwrap()), expected);
    let params = MolEnumeratorParams::new(EnumeratorType::PositionVariation);
    assert_eq!(
        sorted(&enumerate_with_params(&mol, &params).unwrap()),
        expected
    );
    let params = MolEnumeratorParams::new(EnumeratorType::LinkNode);
    assert_eq!(enumerate_with_params(&mol, &params).unwrap().size(), 0);
}

#[test]
fn link_nodes() {
    let mol = smiles(LINK);
    assert_eq!(sorted(&enumerate(&mol).unwrap()), ["COC", "COOC", "COOOC"]);
    assert_eq!(
        sorted(&enumerate_with_max_per_operation(&mol, 2).unwrap()),
        ["COC", "COOC"]
    );
    let params = MolEnumeratorParams {
        max_to_enumerate: 2,
        ..MolEnumeratorParams::new(EnumeratorType::LinkNode)
    };
    assert_eq!(
        sorted(&enumerate_with_params(&mol, &params).unwrap()),
        ["COC", "COOC"]
    );
    assert_value(enumerate(&with_link_nodes("CCC.O", "1 3 2 2 1 2 3")));
}

#[test]
fn repeat_units() {
    let mol = smiles(SRU);
    let expected = ["*N(CC)N(*)CF", "*N(CC)N(*)N(*)CF"];
    assert_eq!(sorted(&enumerate(&mol).unwrap()), expected);
    assert_eq!(
        sorted(&enumerate_with_max_per_operation(&mol, 1).unwrap()),
        ["*N(CC)N(*)CF"]
    );
    let mut params = MolEnumeratorParams::default();
    params.set_enumeration_operator(EnumeratorType::RepeatUnit);
    assert_eq!(
        sorted(&enumerate_with_params(&mol, &params).unwrap()),
        expected
    );
    let mol = smiles("*-CO-* |$star_e;;;star_e$,Sg:n:2,1::ht|");
    assert_eq!(
        sorted(&enumerate(&mol).unwrap()),
        ["**", "*CO*", "*COCO*", "*COCOCO*"]
    );
}

#[test]
fn nothing_to_enumerate_and_errors() {
    assert_eq!(enumerate(&smiles("CCO")).unwrap().size(), 0);
    let mol = smiles(LINK);
    assert_value(enumerate_with_params(&mol, &MolEnumeratorParams::default()));
    let random = MolEnumeratorParams {
        do_random: true,
        ..MolEnumeratorParams::new(EnumeratorType::LinkNode)
    };
    let err = enumerate_with_params(&mol, &random).unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
}

// Position variation removes the dummy atom before the link nodes are
// labeled, which brings the node's atoms to indices 62 to 64.
#[test]
fn position_variation_then_link_nodes() {
    let mut mol = smiles("*Cl.c1ccccc1OCN |m:0:2.3|");
    mol.set_prop("_molLinkNodes", "1 2 2 10 9 10 11");
    assert_eq!(
        sorted(&enumerate(&mol).unwrap()),
        [
            "NCCOc1cccc(Cl)c1",
            "NCCOc1ccccc1Cl",
            "NCOc1cccc(Cl)c1",
            "NCOc1ccccc1Cl"
        ]
    );
    let smi = format!("*Cl.c1ccccc1{}OCN |m:0:2.3|", "C".repeat(55));
    let mut mol = smiles(&smi);
    mol.set_prop("_molLinkNodes", "1 2 2 65 64 65 66");
    assert_eq!(enumerate(&mol).unwrap().size(), 4);
}

// RDKit takes the first product of each link node's reaction without
// checking that there is one; Python crashes on each of these.
#[test]
fn refuses_link_nodes_reactions_cannot_find() {
    // A double bond, which the reaction's single bond does not match.
    assert_value_with(
        enumerate(&with_link_nodes("C=CO", "1 3 2 2 1 2 3")),
        "single",
    );
    // Atoms past 64, whose isotope labels do not fit 16 bits.
    let long = format!("{}O", "C".repeat(70));
    assert_value_with(
        enumerate(&with_link_nodes(&long, "1 2 2 68 67 68 69")),
        "isotopes",
    );
    assert_eq!(
        enumerate(&with_link_nodes(&long, "1 2 2 5 4 5 6"))
            .unwrap()
            .size(),
        2
    );
    // Two nodes on one atom, and both bonds to one atom.
    let two = "two link nodes";
    assert_value_with(
        enumerate(&with_link_nodes("CCO", "1 2 2 2 1 2 3|1 2 2 2 1 2 3")),
        two,
    );
    assert_value_with(
        enumerate(&with_link_nodes("CC(F)O", "1 2 2 2 1 2 4|1 2 2 2 1 2 3")),
        two,
    );
    assert_value_with(
        enumerate(&with_link_nodes("CCO", "1 2 2 2 1 2 1")),
        "one atom",
    );
    let params = MolEnumeratorParams::new(EnumeratorType::LinkNode);
    assert_value(enumerate_with_params(
        &with_link_nodes("C=CO", "1 3 2 2 1 2 3"),
        &params,
    ));
    // Nodes on neighboring atoms are found.
    assert_eq!(
        sorted(&enumerate(&with_link_nodes("CCCO", "1 2 2 2 1 2 3|1 2 2 3 2 3 4")).unwrap()),
        ["CCCCCO", "CCCCO", "CCCCO", "CCCO"]
    );
    // RDKit itself rejects malformed nodes.
    assert_value(enumerate(&with_link_nodes("CCO", "1 2 2 2 1 x 3")));
    assert_value(enumerate(&with_link_nodes("CCO", "1 2 2 9 1 9 3")));
}

// RDKit recurses once per variation point.
#[test]
fn refuses_too_many_variation_points() {
    let mut mol = smiles(&"C".repeat(4002));
    for idx in 0..4001 {
        let mut bond = mol.get_bond_with_idx_mut(idx).unwrap();
        bond.set_prop("_MolFileBondEndPts", "(1 1)");
        bond.set_prop("_MolFileBondAttach", "ANY");
    }
    let deeper = "recursion deeper";
    assert_value_with(enumerate(&mol), deeper);
    let params = MolEnumeratorParams::new(EnumeratorType::PositionVariation);
    assert_value_with(enumerate_with_params(&mol, &params), deeper);
    let nodes = vec!["1 2 2 2 1 2 3"; 4001].join("|");
    assert_value_with(enumerate(&with_link_nodes("CCO", &nodes)), deeper);
}

// Removing atoms rewrites every bond's _MolFileBondEndPts and follows the
// substance groups' PARENT chains.
#[test]
fn refuses_what_removing_atoms_misreads() {
    let mut mol = smiles(PHENOL);
    mol.get_bond_with_idx_mut(3)
        .unwrap()
        .set_prop("_MolFileBondEndPts", "()");
    assert_value_with(enumerate(&mol), "_MolFileBondEndPts");

    let mut mol = smiles(PHENOL);
    for (i, (atom, parent)) in [(1, 2), (3, 1)].into_iter().enumerate() {
        let mut group = create_mol_substance_group(&mut mol, "DAT");
        group.add_atom_with_idx(atom).unwrap();
        group.set_unsigned_prop("index", i as u32 + 1);
        group.set_unsigned_prop("PARENT", parent);
    }
    let mut group = create_mol_substance_group(&mut mol, "DAT");
    group.add_atom_with_idx(0).unwrap();
    assert_value_with(enumerate(&mol), "PARENT");
}

#[test]
fn threads() {
    std::thread::scope(|s| {
        for _ in 0..4 {
            s.spawn(|| {
                let mol = smiles(LINK);
                for _ in 0..20 {
                    assert_eq!(sorted(&enumerate(&mol).unwrap()), ["COC", "COOC", "COOOC"]);
                }
            });
        }
    });
}
