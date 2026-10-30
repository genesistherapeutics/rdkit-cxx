// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_fmcs::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdFMCS (testFMCS.py where it covers a
// case). The searches are small enough to finish well within any timeout.

fn mol(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn mols(texts: &[&str]) -> Vec<Mol> {
    texts.iter().map(|text| mol(text)).collect()
}

fn summary(result: &MCSResult) -> (u32, u32, &str, bool) {
    (
        result.num_atoms(),
        result.num_bonds(),
        result.smarts_string(),
        result.canceled(),
    )
}

fn with<F: FnOnce(&mut MCSParameters)>(f: F) -> MCSParameters {
    let mut params = MCSParameters::default();
    f(&mut params);
    params
}

fn mcs(texts: &[&str], params: &MCSParameters) -> MCSResult {
    find_mcs_with_params(&mols(texts), params).unwrap()
}

#[test]
fn parameter_defaults() {
    let params = MCSParameters::default();
    assert!(params.maximize_bonds);
    assert_eq!(params.threshold, 1.0);
    assert_eq!(params.timeout, 0);
    assert!(!params.verbose);
    assert!(matches!(
        params.atom_typer,
        AtomTyper::Compare(AtomCompare::CompareElements)
    ));
    assert!(matches!(
        params.bond_typer,
        BondTyper::Compare(BondCompare::CompareOrder)
    ));
    assert!(params.progress_callback.is_none());
    assert!(params.final_match_checker.is_none());
    assert!(params.should_accept_mcs.is_none());
    assert_eq!(params.initial_seed, "");
    assert!(!params.store_all);
    let atoms = MCSAtomCompareParameters {
        match_valences: false,
        match_chiral_tag: false,
        match_formal_charge: false,
        ring_matches_ring_only: false,
        complete_rings_only: false,
        match_isotope: false,
        max_distance: -1.0,
    };
    assert_eq!(params.atom_compare_parameters, atoms);
    assert_eq!(MCSAtomCompareParameters::default(), atoms);
    let bonds = MCSBondCompareParameters {
        ring_matches_ring_only: false,
        complete_rings_only: false,
        match_fused_rings: false,
        match_fused_rings_strict: false,
        match_stereo: false,
    };
    assert_eq!(params.bond_compare_parameters, bonds);
    assert_eq!(MCSBondCompareParameters::default(), bonds);
    let stat = MCSProgressData::default();
    assert_eq!(
        (stat.num_atoms(), stat.num_bonds(), stat.seed_processed()),
        (0, 0, 0)
    );
}

#[test]
fn finds_the_mcs() {
    let pair = ["Nc1ccc(cc1)C-Cc1c(N)cccc1", "Nc1ccc(cc1)C=Cc1c(N)cccc1"];
    let result = find_mcs(&mols(&pair)).unwrap();
    assert_eq!(
        summary(&result),
        (8, 8, "[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]", false)
    );
    assert_eq!(result.query_mol().unwrap().get_num_atoms(), 8);
    assert!(result.degenerate_smarts_query_mol_dict().is_empty());

    // testFMCS.py test18GitHub3693.
    let result = mcs(
        &[
            "Nc1ccc(O)cc1c1ccc2ccccc2c1",
            "Oc1cnc(NC2CCC2)c(c1)c1ccc2ccccc2c1",
        ],
        &MCSParameters::default(),
    );
    assert_eq!(
        summary(&result),
        (
            17,
            18,
            "[#7]-,:[#6]:[#6](:[#6]:[#6](:[#6])-[#8])-[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1):[#6]:[#6]:[#6]:[#6]:2",
            false
        )
    );

    // StoreAll keeps both MCSs of the largest size, and no single one.
    let result = mcs(&pair, &with(|p| p.store_all = true));
    assert_eq!(summary(&result), (8, 8, "", false));
    assert!(result.query_mol().is_none());
    let dict = result.degenerate_smarts_query_mol_dict();
    assert_eq!(
        dict.keys().collect::<Vec<_>>(),
        [
            "[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7]",
            "[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]"
        ]
    );
    for query in dict.values() {
        assert_eq!((query.get_num_atoms(), query.get_num_bonds()), (8, 8));
    }
}

#[test]
fn threshold() {
    let texts = ["c1ccccc1CCO", "c1ccccc1CCN", "CCCCN"];
    let all = mcs(&texts, &MCSParameters::default());
    assert_eq!(summary(&all), (4, 3, "[#6]-,:[#6]-,:[#6]-,:[#6]", false));
    let ring = (8, 8, "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]", false);
    for threshold in [0.5, 0.0, -5.0, f64::NEG_INFINITY] {
        let result = mcs(&texts, &with(|p| p.threshold = threshold));
        assert_eq!(summary(&result), ring, "{threshold}");
    }
}

#[test]
fn atom_typers() {
    let pair = ["c1ccccc1CCO", "c1ccccc1CCN"];
    let any = mcs(
        &pair,
        &with(|p| p.atom_typer = AtomCompare::CompareAny.into()),
    );
    assert_eq!(
        summary(&any),
        (
            9,
            9,
            "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]-[#8,#7]",
            false
        )
    );
    let isotopes = mcs(
        &["[13CH3]CC[2H]", "[13CH3]C[13CH3]"],
        &with(|p| p.atom_typer = AtomCompare::CompareIsotopes.into()),
    );
    assert_eq!(summary(&isotopes), (2, 1, "[13*]-[0*]", false));

    let keep_hs = SmilesParserParams {
        remove_hs: false,
        ..SmilesParserParams::default()
    };
    let with_hs = [
        mol_from_smiles_with_params("[H]OCC", &keep_hs).unwrap(),
        mol_from_smiles_with_params("[H]NCCC", &keep_hs).unwrap(),
    ];
    let heavy = find_mcs_with_params(
        &with_hs,
        &with(|p| p.atom_typer = AtomCompare::CompareAnyHeavyAtom.into()),
    )
    .unwrap();
    assert_eq!(summary(&heavy), (4, 3, "[#1]-[#8,#7]-[#6]-[#6]", false));

    let atoms = mcs(&["C1CCC1C", "CCCCC"], &with(|p| p.maximize_bonds = false));
    assert_eq!(summary(&atoms), (5, 4, "[#6]-[#6]-[#6]-[#6]-[#6]", false));
}

#[test]
fn bond_typers() {
    let pair = ["c1ccccc1C", "C1CCCCC1C"];
    let order = mcs(&pair, &MCSParameters::default());
    assert_eq!(
        summary(&order),
        (
            7,
            7,
            "[#6]1:,-[#6]:,-[#6]:,-[#6]:,-[#6]:,-[#6]:,-1-[#6]",
            false
        )
    );
    let exact = mcs(
        &pair,
        &with(|p| p.bond_typer = BondCompare::CompareOrderExact.into()),
    );
    assert_eq!(summary(&exact), (2, 1, "[#6]-[#6]", false));
    let any = mcs(
        &["CC=CC", "CCCC"],
        &with(|p| p.bond_typer = BondCompare::CompareAny.into()),
    );
    assert_eq!(summary(&any), (4, 3, "[#6]-[#6]=,-[#6]-[#6]", false));
}

#[test]
fn ring_compare() {
    let naphthalenes = ["Cc1cccc2ccccc12", "Cc1ccc2ccccc2c1"];
    let decalin = ["C1CCC2CCCCC2C1", "C1CCCCCCCCC1"];
    let ignore = "[#6]-,:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6]:[#6]:[#6]:,-[#6]";
    let ring = "[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-1";
    for (compare, naphthalene, decalin_result) in [
        (
            RingCompare::IgnoreRingFusion,
            (11, 11, ignore),
            (10, 10, ring),
        ),
        (
            RingCompare::PermissiveRingFusion,
            (11, 11, ignore),
            (10, 10, ring),
        ),
        (
            RingCompare::StrictRingFusion,
            (
                11,
                11,
                "[#6]-,:;!@[#6;r6,!R1]1:&@[#6;r6,!R1]:&@[#6;r6,!R1]:&@[#6;r6,!R1]:&@[#6]:&@[#6]:&@1:&@[#6]:&@[#6]:&@[#6]:,-;@[#6]",
            ),
            (
                10,
                9,
                "[#6](-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@[#6])-&@[#6]-&@[#6]-&@[#6]-&@[#6]",
            ),
        ),
    ] {
        let params = with(|p| p.bond_compare_parameters.set_ring_compare(compare));
        let result = mcs(&naphthalenes, &params);
        assert_eq!(
            (
                result.num_atoms(),
                result.num_bonds(),
                result.smarts_string()
            ),
            naphthalene,
            "{compare:?}"
        );
        let result = mcs(&decalin, &params);
        assert_eq!(
            (
                result.num_atoms(),
                result.num_bonds(),
                result.smarts_string()
            ),
            decalin_result,
            "{compare:?}"
        );
    }
}

#[test]
fn initial_seed() {
    // testFMCS.py test7Seed.
    let texts = [
        "C1CCC1CC1CC1",
        "C1CCC1OC1CC1",
        "C1CCC1NC1CC1",
        "C1CCC1SC1CC1",
    ];
    let four = (4, 4, "[#6]1-[#6]-[#6]-[#6]-1", false);
    assert_eq!(summary(&mcs(&texts, &MCSParameters::default())), four);
    let three = mcs(&texts, &with(|p| p.initial_seed = "C1CC1".to_string()));
    assert_eq!(summary(&three), (3, 3, "[#6]1-[#6]-[#6]-1", false));
    // A seed that is not common to the molecules is ignored.
    let ignored = mcs(&texts, &with(|p| p.initial_seed = "C1OC1".to_string()));
    assert_eq!(summary(&ignored), four);
    let rings = (4, 4, "[#6]1-&@[#6]-&@[#6]-&@[#6]-&@1", false);
    let ring_only = with(|p| {
        p.initial_seed = "C1OC1".to_string();
        p.atom_compare_parameters.ring_matches_ring_only = true;
        p.bond_compare_parameters.ring_matches_ring_only = true;
    });
    assert_eq!(summary(&mcs(&texts, &ring_only)), rings);
    let complete = with(|p| {
        p.initial_seed = "C1OC1".to_string();
        p.bond_compare_parameters.complete_rings_only = true;
    });
    assert_eq!(summary(&mcs(&texts, &complete)), rings);

    // RDKit dereferences the parsed seed without a check (Python crashes).
    let bad = with(|p| p.initial_seed = "C((".to_string());
    assert!(matches!(
        find_mcs_with_params(&mols(&texts), &bad),
        Err(Error::Parse { .. })
    ));
}

#[test]
fn atom_compare_parameters() {
    let valences = mcs(
        &["NC1OC1", "C1OC1[N+](=O)[O-]"],
        &with(|p| p.atom_compare_parameters.match_valences = true),
    );
    assert_eq!(summary(&valences), (3, 3, "[#6]1-[#8]-[#6]-1", false));
    let pair = ["CCC[NH3+]", "CCCN"];
    let charge = mcs(
        &pair,
        &with(|p| p.atom_compare_parameters.match_formal_charge = true),
    );
    assert_eq!(summary(&charge), (3, 2, "[#6]-[#6]-[#6]", false));
    assert_eq!(
        summary(&mcs(&pair, &MCSParameters::default())),
        (4, 3, "[#6]-[#6]-[#6]-[#7]", false)
    );
    let alanines = ["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"];
    let chiral = mcs(
        &alanines,
        &with(|p| p.atom_compare_parameters.match_chiral_tag = true),
    );
    assert_eq!(
        summary(&chiral),
        (5, 4, "[#6]-[#6]-[#6](=[#8])-[#8]", false)
    );
    assert_eq!(
        summary(&mcs(&alanines, &MCSParameters::default())),
        (6, 5, "[#6]-[#6](-[#7])-[#6](=[#8])-[#8]", false)
    );
}

#[test]
fn custom_typers() {
    // As testFMCS.py's CompareElements, through the predefined checks.
    let elements = MCSAtomCompare::new(|p, mol1, atom1, mol2, atom2| {
        if mol1.get_atom_with_idx(atom1)?.get_atomic_num()
            != mol2.get_atom_with_idx(atom2)?.get_atomic_num()
        {
            return Ok(false);
        }
        if p.match_formal_charge && !MCSAtomCompare::check_atom_charge(p, mol1, atom1, mol2, atom2)?
        {
            return Ok(false);
        }
        if p.ring_matches_ring_only {
            return MCSAtomCompare::check_atom_ring_match(p, mol1, atom1, mol2, atom2);
        }
        Ok(true)
    });
    let result = mcs(
        &["c1ccccc1CCO", "c1ccccc1CCN"],
        &with(|p| p.atom_typer = elements.clone().into()),
    );
    assert_eq!(
        summary(&result),
        (8, 8, "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]", false)
    );
    let any = MCSBondCompare::new(|_, _, _, _, _| Ok(true));
    let result = mcs(
        &["CC=CC", "CCCC"],
        &with(|p| p.bond_typer = any.clone().into()),
    );
    assert_eq!(summary(&result), (4, 3, "[#6]-[#6]=,-[#6]-[#6]", false));
    let (a, b) = (mol("CO"), mol("CN"));
    let params = MCSAtomCompareParameters::default();
    assert!(elements.call(&params, &a, 0, &b, 0).unwrap());
    assert!(!elements.call(&params, &a, 1, &b, 1).unwrap());
    assert!(matches!(
        elements.call(&params, &a, 2, &b, 0),
        Err(Error::Index { index: 2, len: 2 })
    ));
}

#[test]
fn progress_callback_cancels() {
    // testFMCS.py test17MCSProgressCallbackCancel, with and without a
    // custom atom typer.
    let texts = ["CCCC(C)CC(CC)CC", "OC(N)CC(C)CC(CC)CC"];
    for custom in [false, true] {
        let seen = Arc::new(Mutex::new(Vec::new()));
        let progress = {
            let seen = Arc::clone(&seen);
            MCSProgress::new(move |stat, params| {
                let mut seen = seen.lock().unwrap();
                seen.push((stat.num_atoms(), stat.num_bonds(), stat.seed_processed()));
                assert_eq!(params.timeout, 0);
                Ok(seen.len() < 3)
            })
        };
        let params = with(|p| {
            p.progress_callback = Some(progress);
            if custom {
                p.atom_typer = MCSAtomCompare::new(|_, mol1, atom1, mol2, atom2| {
                    Ok(mol1.get_atom_with_idx(atom1)?.get_atomic_num()
                        == mol2.get_atom_with_idx(atom2)?.get_atomic_num())
                })
                .into();
            }
        });
        let result = mcs(&texts, &params);
        assert_eq!(
            summary(&result),
            (6, 5, "[#6]-[#6]-[#6]-[#6](-[#6])-[#6]", true)
        );
        assert_eq!(*seen.lock().unwrap(), [(3, 2, 0), (4, 3, 0), (6, 5, 0)]);
    }
}

#[test]
fn callbacks_see_the_search_parameters() {
    // RDKit's adjustments show: complete_rings_only implies the rest, and
    // match_chiral_tag implies match_stereo. The given parameters keep
    // their values.
    let seen = Arc::new(Mutex::new(None));
    let progress = {
        let seen = Arc::clone(&seen);
        MCSProgress::new(move |_, params| {
            *seen.lock().unwrap() = Some((
                params.bond_compare_parameters.ring_matches_ring_only,
                params.bond_compare_parameters.complete_rings_only,
                params.atom_compare_parameters.ring_matches_ring_only,
                params.bond_compare_parameters.match_stereo,
                params.threshold,
                params.timeout,
                params.initial_seed.clone(),
                params.progress_callback.is_some(),
            ));
            Ok(true)
        })
    };
    let params = with(|p| {
        p.progress_callback = Some(progress);
        p.atom_compare_parameters.complete_rings_only = true;
        p.atom_compare_parameters.match_chiral_tag = true;
        p.timeout = 7;
    });
    let result = mcs(
        &[
            "Nc1ccc(O)cc1c1ccc2ccccc2c1",
            "Oc1cnc(NC2CCC2)c(c1)c1ccc2ccccc2c1",
        ],
        &params,
    );
    assert_eq!(
        summary(&result),
        (
            10,
            11,
            "[#6]1:&@[#6]:&@[#6]:&@[#6]2:&@[#6](:&@[#6]:&@1):&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2",
            false
        )
    );
    assert_eq!(
        seen.lock().unwrap().clone().unwrap(),
        (true, true, true, true, 1.0, 7, String::new(), true)
    );
    assert!(!params.bond_compare_parameters.complete_rings_only);
}

#[test]
fn match_checks() {
    // testFMCS.py test21ShouldAcceptMCS: a dummy atom matches a substituent.
    let dummy_or_element = MCSAtomCompare::new(|p, mol1, atom1, mol2, atom2| {
        let (a1, a2) = (
            mol1.get_atom_with_idx(atom1)?,
            mol2.get_atom_with_idx(atom2)?,
        );
        let (d1, d2) = (a1.get_atomic_num() == 0, a2.get_atomic_num() == 0);
        if d1 != d2 {
            let (dummy, other) = if d1 { (a1, a2) } else { (a2, a1) };
            return Ok(dummy.get_degree() == 1 && other.get_degree() > 1);
        }
        if a1.get_atomic_num() != a2.get_atomic_num() {
            return Ok(false);
        }
        MCSAtomCompare::check_atom_ring_match(p, mol1, atom1, mol2, atom2)
    });
    let accepted = Arc::new(AtomicUsize::new(0));
    let accept = {
        let accepted = Arc::clone(&accepted);
        MCSAcceptance::new(move |query, target, atom_idx_match, bond_idx_match, _| {
            accepted.fetch_add(1, Ordering::Relaxed);
            assert!(!atom_idx_match.is_empty());
            for &(q, t) in bond_idx_match {
                query.get_bond_with_idx(q)?;
                target.get_bond_with_idx(t)?;
            }
            for &(q, t) in atom_idx_match {
                if query.get_atom_with_idx(q)?.get_atomic_num() == 0
                    || target.get_atom_with_idx(t)?.get_atomic_num() == 0
                {
                    return Ok(true);
                }
            }
            Ok(false)
        })
    };
    let checked = Arc::new(AtomicUsize::new(0));
    let check = {
        let checked = Arc::clone(&checked);
        MCSFinalMatchCheck::new(move |mol1, mol2, atom_idx_match, bond_idx_match, params| {
            checked.fetch_add(1, Ordering::Relaxed);
            assert!(!atom_idx_match.is_empty() && !bond_idx_match.is_empty());
            assert!(params.final_match_checker.is_some());
            for &(q, t) in atom_idx_match {
                mol1.get_atom_with_idx(q)?;
                mol2.get_atom_with_idx(t)?;
            }
            for &(q, t) in bond_idx_match {
                mol1.get_bond_with_idx(q)?;
                mol2.get_bond_with_idx(t)?;
            }
            Ok(true)
        })
    };
    let mut params = with(|p| {
        p.atom_typer = dummy_or_element.into();
        p.bond_typer = BondCompare::CompareAny.into();
        p.should_accept_mcs = Some(accept);
        p.final_match_checker = Some(check);
        p.bond_compare_parameters.ring_matches_ring_only = true;
    });
    let result = mcs(
        &["OC(=O)CCc1c(C)ccc2ccccc12", "c1(C)cc(cccc2)c2c(*)c1"],
        &params,
    );
    assert_eq!(
        summary(&result),
        (
            11,
            12,
            "[#6]1:&@[#6]:&@[#6]2:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2:&@[#6](:&@[#6]:&@1)-&!@[#0,#6]",
            false
        )
    );
    assert!(accepted.load(Ordering::Relaxed) > 0);
    assert!(checked.load(Ordering::Relaxed) > 0);
    params.atom_compare_parameters.ring_matches_ring_only = true;
    let result = mcs(
        &[
            "OC(=O)C(C1CC1)Cc1c(C)ccc2ccccc12",
            "c1(C)cc(cccc2)c2c(C2CC2*)c1",
        ],
        &params,
    );
    assert_eq!(
        summary(&result),
        (4, 4, "[#6]1-&@[#6]-&@[#6]-&@1-&!@[#0,#6;!R]", false)
    );

    let pair = ["c1ccccc1CCO", "c1ccccc1CCN"];
    let reject = with(|p| {
        p.should_accept_mcs = Some(MCSAcceptance::new(|_, _, _, _, _| Ok(false)));
    });
    assert_eq!(summary(&mcs(&pair, &reject)), (0, 0, "", false));
    let reject = with(|p| {
        p.final_match_checker = Some(MCSFinalMatchCheck::new(|_, _, _, _, _| Ok(false)));
    });
    assert_eq!(summary(&mcs(&pair, &reject)), (2, 1, "[#6]:[#6]", false));
}

#[test]
fn callback_errors_are_returned() {
    let pair = mols(&["c1ccccc1CCO", "c1ccccc1CCN"]);
    // The first error stops the search: the closure is not called again.
    let calls = Arc::new(AtomicUsize::new(0));
    let failing = {
        let calls = Arc::clone(&calls);
        MCSAtomCompare::new(move |_, _, _, _, _| {
            calls.fetch_add(1, Ordering::Relaxed);
            Err(Error::Value {
                message: "no".to_string(),
            })
        })
    };
    let result = find_mcs_with_params(&pair, &with(|p| p.atom_typer = failing.into()));
    assert_eq!(
        result.unwrap_err(),
        Error::Value {
            message: "no".to_string()
        }
    );
    assert_eq!(calls.load(Ordering::Relaxed), 1);

    // An error from a closure called later, mid-search.
    let calls = Arc::new(AtomicUsize::new(0));
    let fails_later = {
        let calls = Arc::clone(&calls);
        MCSFinalMatchCheck::new(move |mol1, _, _, _, _| {
            if calls.fetch_add(1, Ordering::Relaxed) == 2 {
                mol1.get_atom_with_idx(100)?;
            }
            Ok(true)
        })
    };
    let result = find_mcs_with_params(&pair, &with(|p| p.final_match_checker = Some(fails_later)));
    assert!(matches!(result, Err(Error::Index { index: 100, .. })));
    assert_eq!(calls.load(Ordering::Relaxed), 3);

    // A panic is caught at the bridge.
    let panics = MCSBondCompare::new(|_, _, _, _, _| panic!("bond compare broke"));
    let result = find_mcs_with_params(&pair, &with(|p| p.bond_typer = panics.into()));
    match result {
        Err(Error::Runtime { message }) => {
            assert!(message.contains("bond compare broke"), "{message}");
        }
        other => panic!("unexpected {other:?}"),
    }
    let panics = MCSProgress::new(|_, _| panic!("progress broke"));
    let result = find_mcs_with_params(&pair, &with(|p| p.progress_callback = Some(panics)));
    assert!(
        matches!(result, Err(Error::Runtime { message }) if message.contains("progress broke"))
    );
    let panics = MCSAcceptance::new(|_, _, _, _, _| panic!("acceptance broke"));
    let result = find_mcs_with_params(&pair, &with(|p| p.should_accept_mcs = Some(panics)));
    assert!(
        matches!(result, Err(Error::Runtime { message }) if message.contains("acceptance broke"))
    );

    // The molecules still work afterwards.
    assert_eq!(find_mcs(&pair).unwrap().num_atoms(), 8);
}

#[test]
fn invalid_inputs() {
    let pair = mols(&["c1ccccc1CCO", "c1ccccc1CCN"]);
    assert!(matches!(find_mcs(&pair[..1]), Err(Error::Runtime { .. })));
    assert!(matches!(find_mcs(&[]), Err(Error::Runtime { .. })));
    assert!(matches!(
        find_mcs_with_params(&pair, &with(|p| p.threshold = 1.5)),
        Err(Error::Runtime { .. })
    ));
    // RDKit converts the NaN to int.
    assert!(matches!(
        find_mcs_with_params(&pair, &with(|p| p.threshold = f64::NAN)),
        Err(Error::Value { .. })
    ));
    // MaxDistance needs conformers.
    assert!(matches!(
        find_mcs_with_params(
            &pair,
            &with(|p| p.atom_compare_parameters.max_distance = 1.0)
        ),
        Err(Error::Value { .. })
    ));
    // Molecules without atoms have no MCS.
    let empty = [Mol::default(), Mol::default()];
    assert_eq!(summary(&find_mcs(&empty).unwrap()), (0, 0, "", false));
}

#[test]
fn ring_info_is_restored() {
    let mut pair = mols(&["c1ccccc1CCO", "c1ccccc1CCN"]);
    for mol in &mut pair {
        mol.clear_computed_props();
        assert!(mol.get_ring_info().num_rings().is_err());
    }
    // RDKit searches with empty ring info and resets it.
    let result = find_mcs(&pair).unwrap();
    assert_eq!(
        summary(&result),
        (8, 8, "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]", false)
    );
    for mol in &pair {
        assert!(mol.get_ring_info().num_rings().is_err());
    }
    // RDKit throws for the missing conformers after it set the empty ring
    // info, which Python then keeps (reporting no rings); it is reset here.
    let result = find_mcs_with_params(
        &pair,
        &with(|p| p.atom_compare_parameters.max_distance = 1.0),
    );
    assert!(matches!(result, Err(Error::Value { .. })));
    for mol in &pair {
        assert!(mol.get_ring_info().num_rings().is_err());
    }
    // So it is after a callback's error.
    let failing = MCSAtomCompare::new(|_, _, _, _, _| {
        Err(Error::Value {
            message: "no".to_string(),
        })
    });
    assert!(find_mcs_with_params(&pair, &with(|p| p.atom_typer = failing.into())).is_err());
    for mol in &pair {
        assert!(mol.get_ring_info().num_rings().is_err());
    }
}

#[test]
fn predefined_checks() {
    let (m1, m2) = (mol("C1CC1C[NH3+]"), mol("CCCC[C@H](F)Cl"));
    let mut atoms = MCSAtomCompareParameters::default();
    assert!(MCSAtomCompare::check_atom_ring_match(&atoms, &m1, 0, &m2, 0).unwrap());
    atoms.ring_matches_ring_only = true;
    assert!(!MCSAtomCompare::check_atom_ring_match(&atoms, &m1, 0, &m2, 0).unwrap());
    assert!(MCSAtomCompare::check_atom_ring_match(&atoms, &m1, 3, &m2, 0).unwrap());
    assert!(!MCSAtomCompare::check_atom_charge(&atoms, &m1, 4, &m2, 0).unwrap());
    assert!(MCSAtomCompare::check_atom_charge(&atoms, &m1, 3, &m2, 0).unwrap());
    assert!(!MCSAtomCompare::check_atom_chirality(&atoms, &m2, 4, &m1, 0).unwrap());
    assert!(MCSAtomCompare::check_atom_chirality(&atoms, &m1, 0, &m2, 4).unwrap());
    assert!(MCSAtomCompare::check_atom_chirality(&atoms, &m2, 4, &m2, 4).unwrap());
    // Python raises RuntimeError from RDKit's range check (or reads ring info
    // past the end, for the ring match); these are checked.
    assert!(matches!(
        MCSAtomCompare::check_atom_charge(&atoms, &m1, 40, &m2, 0),
        Err(Error::Index { index: 40, len: 5 })
    ));
    assert!(matches!(
        MCSAtomCompare::check_atom_ring_match(&atoms, &m1, 0, &m2, 7),
        Err(Error::Index { index: 7, len: 7 })
    ));
    assert!(matches!(
        MCSAtomCompare::check_atom_chirality(&atoms, &m1, 5, &m2, 0),
        Err(Error::Index { index: 5, len: 5 })
    ));

    let (trans, cis, plain) = (mol("C/C=C/C"), mol("C/C=C\\C"), mol("CC=CC"));
    let mut bonds = MCSBondCompareParameters::default();
    assert!(!MCSBondCompare::check_bond_stereo(&bonds, &trans, 1, &cis, 1).unwrap());
    assert!(MCSBondCompare::check_bond_stereo(&bonds, &trans, 1, &trans, 1).unwrap());
    assert!(!MCSBondCompare::check_bond_stereo(&bonds, &trans, 1, &plain, 1).unwrap());
    assert!(!MCSBondCompare::check_bond_stereo(&bonds, &plain, 1, &trans, 1).unwrap());
    assert!(MCSBondCompare::check_bond_stereo(&bonds, &trans, 0, &cis, 1).unwrap());
    let (three, four) = (mol("C1CC1C"), mol("C1CCC1C"));
    // Bond 0 is in each ring, and bonds 2 and 3 are the methyls'.
    assert!(MCSBondCompare::check_bond_ring_match(&bonds, &three, 0, &four, 0).unwrap());
    assert!(!MCSBondCompare::check_bond_ring_match(&bonds, &three, 0, &four, 3).unwrap());
    assert!(MCSBondCompare::check_bond_ring_match(&bonds, &three, 2, &four, 3).unwrap());
    bonds.complete_rings_only = true;
    assert!(!MCSBondCompare::check_bond_ring_match(&bonds, &three, 0, &four, 0).unwrap());
    assert!(MCSBondCompare::check_bond_ring_match(&bonds, &three, 2, &four, 3).unwrap());
    assert!(MCSBondCompare::check_bond_ring_match(&bonds, &three, 0, &three, 1).unwrap());
    assert!(matches!(
        MCSBondCompare::check_bond_stereo(&bonds, &trans, 3, &cis, 1),
        Err(Error::Index { index: 3, len: 3 })
    ));
    assert!(matches!(
        MCSBondCompare::check_bond_ring_match(&bonds, &three, 0, &four, 5),
        Err(Error::Index { index: 5, len: 5 })
    ));

    // Ring checks need ring info.
    let mut bare = mol("C1CC1C");
    bare.clear_computed_props();
    assert!(matches!(
        MCSAtomCompare::check_atom_ring_match(&atoms, &bare, 0, &m2, 0),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        MCSBondCompare::check_bond_ring_match(&bonds, &bare, 0, &four, 0),
        Err(Error::Runtime { .. })
    ));
}

fn unhex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

#[test]
fn unknown_bond_types_are_refused() {
    // Python's pickle of C=C with the bond's type byte set to 30, past
    // Bond::ZERO, which RDKit's order comparisons would read past their
    // table. The pickle itself is refused (check_pickled_mol), so no such
    // bond reaches FindMCS.
    let pkl = unhex(concat!(
        "efbeadde00000000100000000400000000000000020000000100000080010600",
        "680000000302020600680000000302020b0001081e4200000000170400000000",
        "00000016"
    ));
    // SAFETY: a pickle RDKit wrote.
    match unsafe { Mol::from_pkl_string(&pkl) } {
        Err(Error::Value { message }) => assert!(message.contains("BondType"), "{message}"),
        other => panic!("{other:?}"),
    }
}

// A final match check runs inside a substructure match, holding
// control_c_mutex. R-group decomposition takes its own lock before that one
// (CONVENTIONS.md U4), so while another thread decomposes, holding the
// R-group lock and waiting for the match's, a check that decomposed would
// deadlock with it. The check's call fails instead.
#[test]
fn r_group_decomposition_in_final_match_checks_fails() {
    use rdkit_cxx::chem::rd_r_group_decomposition::r_group_decompose;
    use std::sync::mpsc;
    use std::time::Duration;
    let (in_check, check_entered) = mpsc::channel();
    let (done, finished) = mpsc::channel();
    let matcher = {
        let done = done.clone();
        std::thread::spawn(move || {
            let check = MCSFinalMatchCheck::new(move |_, _, _, _, _| {
                let _ = in_check.send(());
                // Let the other thread take the R-group lock first.
                std::thread::sleep(Duration::from_millis(300));
                let result = r_group_decompose(&[mol("c1ccccc1")], &[mol("Cc1ccccc1")]);
                assert!(matches!(result, Err(Error::Runtime { .. })), "{result:?}");
                Ok(true)
            });
            let params = MCSParameters {
                final_match_checker: Some(check),
                ..MCSParameters::default()
            };
            find_mcs_with_params(&mols(&["CCO", "CCN"]), &params).unwrap();
            done.send("matcher").unwrap();
        })
    };
    check_entered.recv_timeout(Duration::from_secs(60)).unwrap();
    let decomposer = std::thread::spawn(move || {
        r_group_decompose(&[mol("c1ccccc1")], &[mol("Cc1ccccc1")]).unwrap();
        done.send("decomposer").unwrap();
    });
    for _ in 0..2 {
        finished
            .recv_timeout(Duration::from_secs(60))
            .expect("the threads deadlocked");
    }
    matcher.join().unwrap();
    decomposer.join().unwrap();
}
