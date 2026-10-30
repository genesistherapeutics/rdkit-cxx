// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_rascal_mces::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdRascalMCES. Searches run without a
// timeout, the only one the bindings accept.

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-5
}

fn sorted(mut pairs: Vec<(u32, u32)>) -> Vec<(u32, u32)> {
    pairs.sort();
    pairs
}

const AD1: &str = "CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934";
const AD2: &str = "N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336";

#[test]
fn defaults() {
    let opts = RascalOptions::default();
    assert_eq!(opts.similarity_threshold, 0.7);
    assert!(opts.complete_aromatic_rings);
    assert!(!opts.ring_matches_ring_only);
    assert!(!opts.exact_connections_match);
    assert!(!opts.single_largest_frag);
    assert_eq!(opts.min_frag_size, -1);
    assert_eq!(opts.max_frag_separation, -1);
    assert!(!opts.all_best_mcess);
    assert_eq!(opts.timeout, -1);
    assert!(!opts.return_empty_mces);
    assert_eq!(opts.max_bond_match_pairs, 1000);
    assert_eq!(opts.equivalent_atoms, "");
    assert!(!opts.ignore_bond_orders);
    assert!(opts.ignore_atom_aromaticity);
    assert_eq!(opts.min_clique_size, 0);

    let clus = RascalClusterOptions::default();
    assert_eq!(clus.similarity_cutoff, 0.7);
    assert_eq!(clus.a, 0.05);
    assert_eq!(clus.b, 2.0);
    assert_eq!(clus.min_frag_size, 3);
    assert_eq!(clus.min_intra_cluster_sim, 0.9);
    assert_eq!(clus.cluster_merge_sim, 0.6);
    assert_eq!(clus.max_num_frags, 2);
    assert_eq!(clus.num_threads, -1);
}

#[test]
fn find_mces_results() {
    let results = find_mces(&smiles("c1ccccc1Cl"), &smiles("c1ccccc1F")).unwrap();
    assert_eq!(results.len(), 1);
    let r = &results[0];
    assert_eq!(r.smarts_string(), "c1:c:c:c:c:c:1");
    assert_eq!(
        sorted(r.bond_matches()),
        [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (6, 6)]
    );
    assert_eq!(
        sorted(r.atom_matches()),
        [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    );
    assert!(close(r.similarity(), 0.734694));
    assert_eq!(r.num_fragments(), 1);
    assert_eq!(r.largest_fragment_size(), 6);
    assert!(close(r.tier1_sim(), 0.734694));
    assert!(close(r.tier2_sim(), 0.734694));
    assert!(!r.timed_out());

    let mut results = find_mces(&smiles(AD1), &smiles(AD2)).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(
        results[0].smarts_string(),
        "N(-C)-c1:c:c:c(-CC(=O)-NCCCCCCCCCC):c:c:1.NC12CC3CC(-C1)-CC(-C2)-C3"
    );
    assert_eq!(results[0].num_fragments(), 2);
    assert!(close(results[0].similarity(), 0.892147));
    let copy = results[0].clone();
    results[0].largest_fragment_only();
    assert_eq!(
        results[0].smarts_string(),
        "N(-C)-c1:c:c:c(-CC(=O)-NCCCCCCCCCC):c:c:1"
    );
    assert_eq!(results[0].num_fragments(), 1);
    assert_eq!(results[0].largest_fragment_size(), 22);
    assert!(close(results[0].similarity(), 0.373529));
    // The copy kept its fragments.
    assert_eq!(copy.num_fragments(), 2);
}

#[test]
fn options() {
    let (ad1, ad2) = (smiles(AD1), smiles(AD2));
    let opts = RascalOptions {
        single_largest_frag: true,
        ..Default::default()
    };
    let results = find_mces_with_opts(&ad1, &ad2, &opts).unwrap();
    assert_eq!(
        results[0].smarts_string(),
        "CCCCCCCCCCNC12CC3CC(-C1)-CC(-C2)-C3"
    );

    let mut opts = RascalOptions {
        similarity_threshold: 0.95,
        ..Default::default()
    };
    assert!(find_mces_with_opts(&ad1, &ad2, &opts).unwrap().is_empty());
    opts.return_empty_mces = true;
    let mut results = find_mces_with_opts(&ad1, &ad2, &opts).unwrap();
    assert_eq!(results.len(), 1);
    let empty = &mut results[0];
    assert_eq!(empty.smarts_string(), "");
    assert!(empty.bond_matches().is_empty());
    assert!(close(empty.tier1_sim(), 0.945398));
    assert!(close(empty.tier2_sim(), -1.0));
    assert_eq!(empty.num_fragments(), 0);
    assert_eq!(empty.similarity(), 0.0);
    // Python's dereferences the molecules this result lacks.
    empty.largest_fragment_only();
    assert_eq!(empty.largest_fragment_size(), 0);

    let a = smiles("Oc1cccc2C(=O)C=CC(=O)c12");
    let b = smiles("O1C(=O)C=Cc2cc(OC)c(O)cc12");
    assert!(find_mces(&a, &b).unwrap().is_empty());
    let opts = RascalOptions {
        similarity_threshold: 0.5,
        ..Default::default()
    };
    let results = find_mces_with_opts(&a, &b, &opts).unwrap();
    assert_eq!(results[0].smarts_string(), "Oc1:c:c:c:c:c:1.[#6]=O");
    assert_eq!(results[0].num_fragments(), 2);

    let opts = RascalOptions {
        equivalent_atoms: "[F,Cl,Br,I]".to_string(),
        ..Default::default()
    };
    let results = find_mces_with_opts(&smiles("c1ccccc1Cl"), &smiles("c1ccccc1F"), &opts).unwrap();
    assert_eq!(results[0].smarts_string(), "c1:c:c:c:c:c:1-[F,Cl,Br,I]");
    assert_eq!(results[0].similarity(), 1.0);

    let opts = RascalOptions {
        similarity_threshold: 0.1,
        all_best_mcess: true,
        ..Default::default()
    };
    let results = find_mces_with_opts(&smiles("CC(C)C"), &smiles("CCCC"), &opts).unwrap();
    assert_eq!(results.len(), 12);
    assert_eq!(results[0].smarts_string(), "C(-C)-C");
    assert_eq!(results[11].smarts_string(), "CCC");

    // RDKit 2024.09 reads a fragment it dropped through a null pointer when
    // it applies a maximum fragment separation, so only -1 is accepted.
    let opts = RascalOptions {
        similarity_threshold: 0.1,
        max_frag_separation: 2,
        ..Default::default()
    };
    assert!(matches!(
        find_mces_with_opts(
            &smiles("c1ccccc1CCCCC(=O)O"),
            &smiles("c1ccccc1CC(=O)O"),
            &opts,
        ),
        Err(Error::Value { .. })
    ));

    let opts = RascalOptions {
        similarity_threshold: 0.5,
        ignore_bond_orders: true,
        ..Default::default()
    };
    let results = find_mces_with_opts(&smiles("CC=CC"), &smiles("CCCC"), &opts).unwrap();
    assert_eq!(results[0].smarts_string(), "C~C~C~C");
}

#[test]
fn checked_options() {
    let (a, b) = (smiles("c1ccccc1CC(=O)O"), smiles("c1ccccc1CCC(=O)O"));
    let find = |opts: RascalOptions| find_mces_with_opts(&a, &b, &opts);
    // RDKit computes the distance matrices only for -1 and up, and reads
    // them for anything but -1.
    assert!(matches!(
        find(RascalOptions {
            max_frag_separation: -2,
            similarity_threshold: 0.1,
            ..Default::default()
        }),
        Err(Error::Value { .. })
    ));
    // RDKit casts a NaN threshold to unsigned int.
    assert!(matches!(
        find(RascalOptions {
            similarity_threshold: f64::NAN,
            ..Default::default()
        }),
        Err(Error::Value { .. })
    ));
    // RDKit dereferences the null molecule an unparsable class makes.
    assert!(matches!(
        find(RascalOptions {
            equivalent_atoms: "[F,Cl] [C".to_string(),
            ..Default::default()
        }),
        Err(Error::Parse { .. })
    ));
    // RDKit refuses more than nine classes itself.
    assert!(matches!(
        find(RascalOptions {
            equivalent_atoms: "C C C C C C C C C C".to_string(),
            ..Default::default()
        }),
        Err(Error::Value { .. })
    ));
    assert!(
        find(RascalOptions {
            equivalent_atoms: " [F,Cl]  [N,O] ".to_string(),
            ..Default::default()
        })
        .is_ok()
    );
}

#[test]
fn molecule_sizes_are_checked() {
    // The similarity estimates square a count of matched atoms and bonds in
    // int: 23171 atoms and 23170 bonds each is one past.
    let chain = smiles(&"C".repeat(23171));
    assert!(matches!(
        find_mces(&chain, &chain),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        rascal_cluster(&[chain.clone(), chain.clone()]),
        Err(Error::Value { .. })
    ));
    // One such molecule bounds every count by the other's size.
    let small = smiles("CCCC");
    assert!(find_mces(&chain, &small).unwrap().is_empty());
}

const CDK2: [&str; 14] = [
    "CC(C)C(=O)COc1nc(N)nc2[nH]cnc12",
    "Nc1nc(OCC2CCCO2)c2nc[nH]c2n1",
    "Nc1nc(OCC2CCC(=O)N2)c2nc[nH]c2n1",
    "Nc1nc(OCC2CCCCC2)c2nc[nH]c2n1",
    "Nc1nc(OCC2CC=CCC2)c2nc[nH]c2n1",
    "Cn1cnc2c(NCc3ccccc3)nc(NCCO)nc21",
    "CCC(CO)Nc1nc(NCc2ccccc2)c2ncn(C(C)C)c2n1",
    "COc1ccc(CNc2nc(N(CCO)CCO)nc3c2ncn3C(C)C)cc1",
    "Nc1nc(N)c(N=O)c(OCC2CCCCC2)n1",
    "COc1ccc2c(c1)C(=Cc1cnc[nH]1)C(=O)N2",
    "COc1cc[nH]c1C=C1C(=O)Nc2ccc([N+](=O)[O-])cc21",
    "CCc1cnc(CSc2cnc(NC(=O)C(C)C)s2)o1",
    "CCCCOc1c(C(=O)c2nccs2)cnc2[nH]ncc12",
    "COc1cc(-c2ccc[nH]2)c2c3c(ccc(F)c13)NC2=O",
];

#[test]
fn clusters() {
    let mols: Vec<Mol> = CDK2.iter().map(|s| smiles(s)).collect();
    let fuzzy = vec![
        vec![3, 4, 1, 2, 0],
        vec![6, 5, 7],
        vec![9, 10, 8, 11, 12, 13],
    ];
    let butina = vec![
        vec![4, 3, 1, 2, 0],
        vec![7, 6, 5],
        vec![13],
        vec![12],
        vec![11],
        vec![10],
        vec![9],
        vec![8],
    ];
    assert_eq!(rascal_cluster(&mols).unwrap(), fuzzy);
    assert_eq!(rascal_butina_cluster(&mols).unwrap(), butina);
    // num_threads changes nothing: the pairs run on the calling thread.
    for num_threads in [1, 2, 4, 64, 0, i32::MIN] {
        let opts = RascalClusterOptions {
            num_threads,
            ..Default::default()
        };
        assert_eq!(rascal_cluster_with_opts(&mols, &opts).unwrap(), fuzzy);
        assert_eq!(
            rascal_butina_cluster_with_opts(&mols, &opts).unwrap(),
            butina
        );
    }

    let opts = RascalClusterOptions {
        similarity_cutoff: 0.6,
        ..Default::default()
    };
    assert_eq!(
        rascal_cluster_with_opts(&mols, &opts).unwrap(),
        [
            vec![3, 4, 1, 2, 0, 8],
            vec![6, 5, 7],
            vec![10, 13, 9, 11, 12]
        ]
    );
    assert_eq!(
        rascal_butina_cluster_with_opts(&mols, &opts).unwrap(),
        [
            vec![4, 3, 1, 2, 0, 8],
            vec![7, 6, 5],
            vec![13],
            vec![12],
            vec![11],
            vec![10],
            vec![9]
        ]
    );

    let empty: Vec<Vec<u32>> = vec![vec![]];
    assert_eq!(rascal_cluster(&[]).unwrap(), empty);
    assert_eq!(rascal_cluster(&mols[..1]).unwrap(), empty);
    assert!(rascal_butina_cluster(&[]).unwrap().is_empty());
    assert!(rascal_butina_cluster(&mols[..1]).unwrap().is_empty());
    assert!(matches!(
        rascal_cluster_with_opts(
            &mols,
            &RascalClusterOptions {
                similarity_cutoff: f64::NAN,
                ..Default::default()
            }
        ),
        Err(Error::Value { .. })
    ));
}

#[test]
fn errors_on_cluster_threads_are_returned() {
    // A molecule without ring info makes RDKit throw. Python's threads
    // catch nothing, so the process ends; here the pairs run on the calling
    // thread, and the error is returned.
    let without_rings = |text: &str| {
        let mut mol = smiles(text);
        mol.clear_computed_props();
        mol
    };
    let mut mols: Vec<Mol> = CDK2[..5].iter().map(|s| smiles(s)).collect();
    mols.push(without_rings("c1ccccc1CCO"));
    mols.push(without_rings("c1ccccc1CCN"));
    for num_threads in [1, 4] {
        let opts = RascalClusterOptions {
            num_threads,
            ..Default::default()
        };
        assert!(matches!(
            rascal_cluster_with_opts(&mols, &opts),
            Err(Error::Runtime { .. })
        ));
        assert!(matches!(
            rascal_butina_cluster_with_opts(&mols, &opts),
            Err(Error::Runtime { .. })
        ));
    }
    assert!(matches!(
        find_mces(&mols[5], &mols[6]),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn results_move_between_threads() {
    let results = find_mces(&smiles("c1ccccc1Cl"), &smiles("c1ccccc1F")).unwrap();
    let handle = std::thread::spawn(move || results[0].smarts_string());
    assert_eq!(handle.join().unwrap(), "c1:c:c:c:c:c:1");
}

/// RDKit 2024.09's checkTimeout reads past the best cliques it cleared once
/// the search runs past the timeout, which depends on the clock.
#[test]
fn rascal_timeout_must_be_none() {
    for timeout in [-2, 0, 60] {
        let opts = RascalOptions {
            timeout,
            ..Default::default()
        };
        let err =
            find_mces_with_opts(&smiles("c1ccccc1Cl"), &smiles("c1ccccc1F"), &opts).unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{timeout}: {err:?}");
    }
    find_mces_with_opts(
        &smiles("c1ccccc1Cl"),
        &smiles("c1ccccc1F"),
        &RascalOptions::default(),
    )
    .unwrap();
}
