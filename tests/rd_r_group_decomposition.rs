// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_r_group_decomposition::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdRGroupDecomposition.

const MOLS: [&str; 4] = ["c1ccccc1Cl", "c1ccccc1CCO", "Cc1ccccc1Br", "CCCC"];
const CORE: &str = "c1ccc([*:2])c([*:1])c1";

fn mols(smiles: &[&str]) -> Vec<Mol> {
    smiles.iter().map(|s| mol_from_smiles(s).unwrap()).collect()
}

fn benzene() -> Vec<Mol> {
    mols(&["c1ccccc1"])
}

fn row(entries: &[(&str, &str)]) -> BTreeMap<String, String> {
    entries
        .iter()
        .map(|(label, smiles)| (label.to_string(), smiles.to_string()))
        .collect()
}

fn expected_rows() -> Vec<BTreeMap<String, String>> {
    vec![
        row(&[("Core", CORE), ("R1", "[H][*:1]"), ("R2", "Cl[*:2]")]),
        row(&[("Core", CORE), ("R1", "[H][*:1]"), ("R2", "OCC[*:2]")]),
        row(&[("Core", CORE), ("R1", "C[*:1]"), ("R2", "Br[*:2]")]),
    ]
}

#[test]
fn params_defaults() {
    let p = RGroupDecompositionParameters::default();
    assert_eq!(p.labels, RGroupLabels::AutoDetect);
    assert_eq!(p.labels.bits(), 0xFF);
    assert_eq!(p.matching_strategy, RGroupMatching::GreedyChunks);
    assert_eq!(p.score_method, RGroupScore::Match);
    assert_eq!(
        p.rgroup_labelling,
        RGroupLabelling::AtomMap | RGroupLabelling::MDLRGroup
    );
    assert_eq!(p.alignment, RGroupCoreAlignment::MCS);
    assert_eq!(p.chunk_size, 5);
    assert!(!p.only_match_at_r_groups);
    assert!(p.remove_all_hydrogen_r_groups);
    assert!(p.remove_all_hydrogen_r_groups_and_labels);
    assert!(p.remove_hydrogens_post_match);
    assert!(!p.allow_non_terminal_r_groups);
    assert!(!p.allow_multiple_r_groups_on_unlabelled);
    assert!(!p.allow_multiple_cores_in_same_mol);
    assert!(!p.do_tautomers);
    assert!(!p.do_enumeration);
    assert!(!p.include_target_mol_in_results);
    assert_eq!(p.timeout, -1.0);
    assert_eq!(p.ga_population_size, -1);
    assert_eq!(p.ga_maximum_operations, -1);
    assert_eq!(p.ga_number_operations_without_improvement, -1);
    assert_eq!(p.ga_random_seed, -1);
    assert_eq!(p.ga_number_runs, 1);
    assert!(p.ga_parallel_runs);
    assert!(p.substruct_match_params.use_chirality);
    assert_eq!(format!("{:?}", RGroupMatching::GA), "RGroupMatching::GA");
    assert_eq!(RGroupLabels::IsotopeLabels.bits(), 1);
    assert_eq!(RGroupCoreAlignment::NoAlignment.bits(), 0);
    assert_eq!(RGroupScore::FingerprintVariance.bits(), 4);
}

#[test]
fn decompose() {
    let (rows, unmatched) = r_group_decompose_as_smiles(&benzene(), &mols(&MOLS)).unwrap();
    assert_eq!(rows, expected_rows());
    assert_eq!(unmatched, [3]);

    let (rows, unmatched) = r_group_decompose(&benzene(), &mols(&MOLS)).unwrap();
    assert_eq!(unmatched, [3]);
    assert_eq!(mol_to_smiles(&rows[2]["R2"]).unwrap(), "Br[*:2]");

    let (columns, unmatched) =
        r_group_decompose_as_columns_as_smiles(&benzene(), &mols(&MOLS)).unwrap();
    assert_eq!(unmatched, [3]);
    assert_eq!(columns["Core"], [CORE; 3]);
    assert_eq!(columns["R1"], ["[H][*:1]", "[H][*:1]", "C[*:1]"]);
    assert_eq!(columns["R2"], ["Cl[*:2]", "OCC[*:2]", "Br[*:2]"]);
    let (columns, _) = r_group_decompose_as_columns(&benzene(), &mols(&MOLS)).unwrap();
    assert_eq!(columns.len(), 3);
    assert_eq!(columns["R1"].len(), 3);

    // Two cores, from most to least specific.
    let (rows, unmatched) = r_group_decompose_as_smiles(
        &mols(&["c1ccncc1", "c1ccccc1"]),
        &mols(&["Clc1ccncc1", "Brc1ccccc1"]),
    )
    .unwrap();
    assert!(unmatched.is_empty());
    assert_eq!(
        rows,
        [
            row(&[("Core", "c1cc([*:1])ccn1"), ("R1", "Cl[*:1]")]),
            row(&[("Core", "c1ccc([*:1])cc1"), ("R1", "Br[*:1]")]),
        ]
    );

    // A labeled core, matched anywhere or only at its labels.
    let labelled = mols(&["[*:1]c1ccc([*:2])cc1"]);
    let expected = [
        row(&[
            ("Core", "c1cc([*:2])ccc1[*:1]"),
            ("R1", "O[*:1]"),
            ("R2", "Cl[*:2]"),
        ]),
        row(&[
            ("Core", "c1cc([*:2])ccc1[*:1]"),
            ("R1", "Br[*:1]"),
            ("R2", "[H][*:2]"),
        ]),
    ];
    let (rows, unmatched) = r_group_decompose_as_smiles(
        &labelled,
        &mols(&["Clc1ccc(O)cc1", "Brc1ccccc1", "c1ccncc1"]),
    )
    .unwrap();
    assert_eq!(rows, expected);
    assert_eq!(unmatched, [2]);
    let only = RGroupDecompositionParameters {
        only_match_at_r_groups: true,
        ..RGroupDecompositionParameters::default()
    };
    let (rows, unmatched) = r_group_decompose_as_smiles_with_options(
        &labelled,
        &mols(&["Clc1ccc(O)cc1", "Brc1ccccc1", "Clc1ccc(O)cc1C"]),
        &only,
    )
    .unwrap();
    assert_eq!(rows, expected);
    assert_eq!(unmatched, [2]);
}

#[test]
fn options() {
    for options in [
        RGroupDecompositionParameters {
            score_method: RGroupScore::FingerprintVariance,
            ..RGroupDecompositionParameters::default()
        },
        // Too few permutations for the GA: RDKit searches them all.
        RGroupDecompositionParameters {
            matching_strategy: RGroupMatching::GA,
            ga_number_runs: 2,
            ga_random_seed: 42,
            ..RGroupDecompositionParameters::default()
        },
    ] {
        let (rows, unmatched) =
            r_group_decompose_as_smiles_with_options(&benzene(), &mols(&MOLS), &options).unwrap();
        assert_eq!(rows, expected_rows());
        assert_eq!(unmatched, [3]);
    }

    let target = RGroupDecompositionParameters {
        include_target_mol_in_results: true,
        ..RGroupDecompositionParameters::default()
    };
    let (rows, _) =
        r_group_decompose_as_smiles_with_options(&benzene(), &mols(&MOLS[..1]), &target).unwrap();
    assert_eq!(
        rows,
        [row(&[
            ("Core", "c1ccc([*:1])cc1"),
            ("Mol", "Clc1ccccc1"),
            ("R1", "Cl[*:1]")
        ])]
    );

    let timeout = RGroupDecompositionParameters {
        timeout: 1e-9,
        ..RGroupDecompositionParameters::default()
    };
    assert!(matches!(
        r_group_decompose_with_options(&benzene(), &mols(&MOLS), &timeout),
        Err(Error::Runtime { message }) if message.contains("timed out")
    ));
}

#[test]
fn decomposition() {
    let mut decomp = RGroupDecomposition::new(&benzene()).unwrap();
    let added: Vec<i32> = mols(&MOLS)
        .iter()
        .map(|mol| decomp.add(mol).unwrap())
        .collect();
    assert_eq!(added, [0, 1, 2, -1]);
    // Before processing, RDKit's labels are not final.
    assert!(matches!(
        decomp.get_r_groups_as_rows(),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        decomp.get_r_groups_as_columns(),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(
        decomp.add(&mol_from_smiles("c1ccccc1").unwrap()).unwrap(),
        3
    );

    let (idx, outputs) = decomp
        .get_matching_core_idx_with_outputs(&mol_from_smiles("Clc1ccccc1Br").unwrap())
        .unwrap();
    assert_eq!(idx, 0);
    assert_eq!(outputs.matches.len(), 12);
    assert_eq!(
        outputs.matches[0],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]
    );
    assert_eq!(
        decomp
            .get_matching_core_idx(&mol_from_smiles("CC").unwrap())
            .unwrap(),
        -1
    );

    assert!(decomp.process());
    assert_eq!(decomp.get_r_group_labels().unwrap(), ["Core", "R1", "R2"]);
    let mut expected = expected_rows();
    expected.push(row(&[
        ("Core", CORE),
        ("R1", "[H][*:1]"),
        ("R2", "[H][*:2]"),
    ]));
    assert_eq!(decomp.get_r_groups_as_rows_as_smiles().unwrap(), expected);
    // RDKit rewrites its matched cores when it returns them; again gives
    // the same.
    assert_eq!(decomp.get_r_groups_as_rows_as_smiles().unwrap(), expected);
    let columns = decomp.get_r_groups_as_columns_as_smiles().unwrap();
    assert_eq!(
        columns["R1"],
        ["[H][*:1]", "[H][*:1]", "C[*:1]", "[H][*:1]"]
    );
    assert_eq!(decomp.get_r_groups_as_columns().unwrap()["R2"].len(), 4);
    let (success, score) = decomp.process_and_score();
    assert!(success);
    assert!(score.is_finite());

    // Processing nothing fails, and leaves nothing.
    let mut empty = RGroupDecomposition::new(&benzene()).unwrap();
    assert!(!empty.process());
    assert!(empty.get_r_groups_as_rows().unwrap().is_empty());
}

// RDKit takes the number of molecules modulo chunkSize; Python dies of
// SIGFPE with chunkSize 0.
#[test]
fn chunk_size_zero() {
    let params = RGroupDecompositionParameters {
        chunk_size: 0,
        ..RGroupDecompositionParameters::default()
    };
    assert!(matches!(
        RGroupDecomposition::new_with_params(&benzene(), &params),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        r_group_decompose_with_options(&benzene(), &mols(&MOLS), &params),
        Err(Error::Value { .. })
    ));
    // Greedy never divides by it.
    let greedy = RGroupDecompositionParameters {
        matching_strategy: RGroupMatching::Greedy,
        ..params
    };
    let (rows, _) =
        r_group_decompose_as_smiles_with_options(&benzene(), &mols(&MOLS), &greedy).unwrap();
    assert_eq!(rows.len(), 3);
}

// prepareCore adds one to labels in int.
#[test]
fn core_labels() {
    let setters: [fn(&mut AtomMut<'_>); 3] = [
        |atom: &mut AtomMut| atom.set_int_prop("molAtomMapNumber", i32::MAX),
        |atom: &mut AtomMut| atom.set_int_prop("tempRlabel", i32::MAX),
        |atom: &mut AtomMut| atom.set_unsigned_prop("_MolFileRLabel", 1 << 30),
    ];
    for (i, set) in setters.into_iter().enumerate() {
        let mut core = mol_from_smiles("c1ccccc1[1*]").unwrap();
        set(&mut core.get_atom_with_idx_mut(6).unwrap());
        let result = RGroupDecomposition::new(&[core]);
        assert!(
            matches!(result, Err(Error::Value { .. })),
            "{i}: {result:?}"
        );
    }
    let mut core = mol_from_smiles("c1ccccc1[1*]").unwrap();
    core.get_atom_with_idx_mut(6)
        .unwrap()
        .set_int_prop("molAtomMapNumber", (1 << 30) - 1);
    RGroupDecomposition::new(&[core]).unwrap();
}

// MCS alignment's bond comparison indexes a table by bond types, which a
// pickle could set past BondType::Zero; unpickling refuses such a type
// (check_pickled_mol), so no such core reaches the decomposition.
#[test]
fn core_bond_types() {
    let mut pickle = mol_from_smiles("CC=O").unwrap().to_binary();
    assert_eq!((pickle[55], pickle[62]), (11, 2));
    pickle[62] = 30;
    assert!(matches!(
        Mol::from_pkl_string(&pickle),
        Err(Error::Value { .. })
    ));
}

// The core with its first atom replaced by the SMARTS.
fn with_first_atom(smarts: &str) -> Mol {
    let mut mol = RWMol::from(mol_from_smarts(smarts).unwrap());
    mol.insert_mol(&mol_from_smiles("CCC(=O)C[*:1]").unwrap());
    mol.add_bond_with_order(0, 1, BondType::Single).unwrap();
    mol.update_property_cache_with_strict(false).unwrap();
    mol.into()
}

// doTautomers matches with RDKit's tautomer queries, whose check reads a
// recursive SMARTS's own matches at the template's modified atoms: Python
// fails a precondition on a bond between atoms it read past the match.
#[test]
fn recursive_cores_with_tautomers() {
    let mols = mols(&["CCCC(=O)CCl", "CCCC(O)=CBr", "CC(C)CC(=O)CO"]);
    let tautomers = RGroupDecompositionParameters {
        do_tautomers: true,
        ..RGroupDecompositionParameters::default()
    };
    let (rows, unmatched) =
        r_group_decompose_as_smiles_with_options(&[with_first_atom("[CH3]")], &mols, &tautomers)
            .unwrap();
    assert_eq!(rows.len(), 3);
    assert!(unmatched.is_empty());
    assert_eq!(rows[1]["Core"], "CC(CC(O)=C[*:1])[*:2]");

    let recursive = [with_first_atom("[$([CH3])]")];
    assert!(matches!(
        r_group_decompose_with_options(&recursive, &mols, &tautomers),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        RGroupDecomposition::new_with_params(&recursive, &tautomers),
        Err(Error::Value { .. })
    ));
    let (rows, unmatched) = r_group_decompose_as_smiles_with_options(
        &recursive,
        &mols,
        &RGroupDecompositionParameters::default(),
    )
    .unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0]["Core"], "O=C(CC[*:2])C[*:1]");
    assert_eq!(unmatched, [1]);
}

// doEnumeration enumerates the cores' link nodes, which RDKit trusts:
// Python segfaults on a link node over a double bond.
#[test]
fn link_node_cores() {
    let core = |text: &str, nodes: &str| {
        let mut mol = mol_from_smiles(text).unwrap();
        mol.set_prop("_molLinkNodes", nodes);
        mol
    };
    let enumeration = RGroupDecompositionParameters {
        do_enumeration: true,
        ..RGroupDecompositionParameters::default()
    };
    let target = mols(&["C1=CCCC1CO"]);
    let (rows, unmatched) = r_group_decompose_as_smiles_with_options(
        &[core("C1CCCC1", "1 3 2 2 1 2 3")],
        &target,
        &enumeration,
    )
    .unwrap();
    assert!(rows.is_empty());
    assert_eq!(unmatched, [0]);
    let bad = [core("C1=CCC1", "1 3 2 1 2 1 4")];
    assert!(matches!(
        r_group_decompose_with_options(&bad, &target, &enumeration),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        RGroupDecomposition::new_with_params(&bad, &enumeration),
        Err(Error::Value { .. })
    ));
    // Without enumeration the link nodes are not read.
    RGroupDecomposition::new(&bad).unwrap();
}

#[test]
fn relabel() {
    let mut mol = mol_from_smiles("[*:1]c1ccccc1[1*]").unwrap();
    relabel_mapped_dummies(&mut mol).unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "*c1ccccc1*");
    let mut mol = mol_from_smiles("[*:1]c1ccccc1[1*]").unwrap();
    relabel_mapped_dummies_with_params(
        &mut mol,
        &RelabelMappedDummiesParams {
            output_labels: RGroupLabelling::AtomMap | RGroupLabelling::Isotope,
            ..RelabelMappedDummiesParams::default()
        },
    )
    .unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "c1ccc([1*:1])c([1*:1])c1");

    // RDKit takes abs() of the atom map number.
    let mut mol = mol_from_smiles("*c1ccccc1").unwrap();
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_int_prop("molAtomMapNumber", i32::MIN);
    assert!(matches!(
        relabel_mapped_dummies(&mut mol),
        Err(Error::Value { .. })
    ));
}

// Decompositions match under control_c_mutex and process under
// rgroup_mutex; several threads at once agree.
#[test]
fn threads() {
    let handles: Vec<_> = (0..4)
        .map(|i| {
            std::thread::spawn(move || {
                let options = RGroupDecompositionParameters {
                    score_method: if i % 2 == 0 {
                        RGroupScore::Match
                    } else {
                        RGroupScore::FingerprintVariance
                    },
                    ..RGroupDecompositionParameters::default()
                };
                (0..3)
                    .map(|_| {
                        r_group_decompose_as_smiles_with_options(&benzene(), &mols(&MOLS), &options)
                            .unwrap()
                    })
                    .collect::<Vec<_>>()
            })
        })
        .collect();
    for handle in handles {
        for (rows, unmatched) in handle.join().unwrap() {
            assert_eq!(rows, expected_rows());
            assert_eq!(unmatched, [3]);
        }
    }
}

// Enumerating the cores takes RDKit's enumerator a recursion level per
// variation point, as rdMolEnumerator's own checks say.
#[test]
fn enumeration_refuses_too_many_variation_points() {
    let mut core = mol_from_smiles(&"C".repeat(4002)).unwrap();
    for idx in 0..4001 {
        let mut bond = core.get_bond_with_idx_mut(idx).unwrap();
        bond.set_prop("_MolFileBondEndPts", "(1 1)");
        bond.set_prop("_MolFileBondAttach", "ANY");
    }
    let params = RGroupDecompositionParameters {
        do_enumeration: true,
        ..Default::default()
    };
    match RGroupDecomposition::new_with_params(&[core], &params).unwrap_err() {
        Error::Value { message } => assert!(message.contains("recursion deeper"), "{message}"),
        err => panic!("{err:?}"),
    }
}
