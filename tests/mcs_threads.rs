// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The substructure matches of rdFMCS, rdRascalMCES, rdTautomerQuery, and
//! rdGeneralizedSubstruct on many threads at once. Each match constructs a
//! ControlCHandler, so they share control_c_mutex (CONVENTIONS.md U4); FMCS
//! releases it for its closures, except a final match check, which runs
//! inside a match.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_depictor::compute_2d_coords;
use rdkit_cxx::chem::rd_fmcs::*;
use rdkit_cxx::chem::rd_generalized_substruct::*;
use rdkit_cxx::chem::rd_rascal_mces::*;
use rdkit_cxx::chem::rd_tautomer_query::*;
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

const PAIR: [&str; 2] = ["Nc1ccc(cc1)C-Cc1c(N)cccc1", "Nc1ccc(cc1)C=Cc1c(N)cccc1"];
const PAIR_MCS: &str = "[#7]-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]";

/// FMCS whose atom typer and progress callback match patterns themselves,
/// as closures may: the search releases the lock for them. (Molecules and
/// queries are not Sync, so the closures hold them behind mutexes.)
fn fmcs_with_matching_typer() -> MCSResult {
    let nitrogen = Mutex::new(smiles("N"));
    let amine = Mutex::new(create_extended_query_mol(&smiles("N")).unwrap());
    let typer = MCSAtomCompare::new(move |_, mol1, atom1, mol2, atom2| {
        assert!(mol1.has_substruct_match(&nitrogen.lock().unwrap())?);
        assert!(mol_has_substruct_match(mol2, &amine.lock().unwrap())?);
        Ok(mol1.get_atom_with_idx(atom1)?.get_atomic_num()
            == mol2.get_atom_with_idx(atom2)?.get_atomic_num())
    });
    let tautomers = Mutex::new(TautomerQuery::new(&smiles("O=C1CCCCC1")).unwrap());
    let progress = MCSProgress::new(move |_, _| {
        tautomers
            .lock()
            .unwrap()
            .is_substruct_of(&smiles("OC1=CCCC(CC)C1"))
    });
    let params = MCSParameters {
        atom_typer: typer.into(),
        progress_callback: Some(progress),
        ..MCSParameters::default()
    };
    let mols: Vec<Mol> = PAIR.iter().map(|text| smiles(text)).collect();
    find_mcs_with_params(&mols, &params).unwrap()
}

#[test]
fn matching_on_threads() {
    std::thread::scope(|scope| {
        for t in 0..8 {
            scope.spawn(move || {
                for _ in 0..3 {
                    match t % 4 {
                        0 => {
                            let result = fmcs_with_matching_typer();
                            assert_eq!(result.smarts_string(), PAIR_MCS);
                        }
                        1 => {
                            let results =
                                find_mces(&smiles("c1ccccc1Cl"), &smiles("c1ccccc1F")).unwrap();
                            assert_eq!(results[0].smarts_string(), "c1:c:c:c:c:c:1");
                            let mols: Vec<Mol> = ["c1ccccc1CCO", "c1ccccc1CCN", "CCCC"]
                                .iter()
                                .map(|text| smiles(text))
                                .collect();
                            rascal_butina_cluster(&mols).unwrap();
                        }
                        2 => {
                            let query = TautomerQuery::new(&smiles("O=C1CCCCC1")).unwrap();
                            let target = smiles("OC1=CCCC(CC)C1");
                            assert!(query.is_substruct_of(&target).unwrap());
                            assert_eq!(query.get_substruct_matches(&target).unwrap().len(), 1);
                            query.pattern_fingerprint_template().unwrap();
                        }
                        _ => {
                            let query =
                                create_extended_query_mol(&smiles("COC |LN:1:1.3|")).unwrap();
                            assert!(mol_has_substruct_match(&smiles("COOC"), &query).unwrap());
                            assert!(!mol_has_substruct_match(&smiles("CC"), &query).unwrap());
                            pattern_fingerprint_target(&smiles("COOC")).unwrap();
                        }
                    }
                }
            });
        }
    });
}

#[test]
fn final_match_checks_cannot_take_the_lock_again() {
    // The check runs inside vf2, with a ControlCHandler alive, so the search
    // keeps the lock: wrappers that would take it, or a lock before it,
    // fail instead of deadlocking.
    let failures = Arc::new(AtomicUsize::new(0));
    let check = {
        let failures = Arc::clone(&failures);
        MCSFinalMatchCheck::new(move |mol1, mol2, _, _, _| {
            let nested = find_mcs(&[mol1.clone(), mol2.clone()]);
            assert!(matches!(nested, Err(Error::Runtime { .. })), "{nested:?}");
            let query = TautomerQuery::new(mol1);
            assert!(matches!(query, Err(Error::Runtime { .. })));
            let mut copy = mol1.clone();
            assert!(matches!(
                compute_2d_coords(&mut copy),
                Err(Error::Runtime { .. })
            ));
            assert!(matches!(
                fragment_on_brics_bonds(mol2),
                Err(Error::Runtime { .. })
            ));
            failures.fetch_add(1, Ordering::Relaxed);
            Ok(true)
        })
    };
    let params = MCSParameters {
        final_match_checker: Some(check),
        ..MCSParameters::default()
    };
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                let mols: Vec<Mol> = PAIR.iter().map(|text| smiles(text)).collect();
                let result = find_mcs_with_params(&mols, &params).unwrap();
                assert_eq!(result.smarts_string(), PAIR_MCS);
            });
        }
        scope.spawn(|| {
            for _ in 0..3 {
                assert_eq!(fmcs_with_matching_typer().smarts_string(), PAIR_MCS);
            }
        });
    });
    assert!(failures.load(Ordering::Relaxed) > 0);
    // Outside a check, the same calls work.
    let mols: Vec<Mol> = PAIR.iter().map(|text| smiles(text)).collect();
    assert_eq!(find_mcs(&mols).unwrap().smarts_string(), PAIR_MCS);
    let mut copy = mols[0].clone();
    compute_2d_coords(&mut copy).unwrap();
    fragment_on_brics_bonds(&mols[1]).unwrap();
}

#[test]
fn an_error_from_a_nested_call_is_returned() {
    let check = MCSFinalMatchCheck::new(|mol1, mol2, _, _, _| {
        find_mcs(&[mol1.clone(), mol2.clone()])?;
        Ok(true)
    });
    let params = MCSParameters {
        final_match_checker: Some(check),
        ..MCSParameters::default()
    };
    let mols: Vec<Mol> = PAIR.iter().map(|text| smiles(text)).collect();
    let result = find_mcs_with_params(&mols, &params);
    assert!(matches!(result, Err(Error::Runtime { .. })), "{result:?}");
}
