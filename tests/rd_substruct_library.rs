// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_generalized_substruct::{
    CreateExtendedQueryMolParams, ExtendedQueryMol, create_extended_query_mol,
    create_extended_query_mol_with_params,
};
use rdkit_cxx::chem::rd_substruct_library::*;
use rdkit_cxx::chem::rd_tautomer_query::TautomerQuery;
use rdkit_cxx::chem::*;
use rdkit_cxx::data_structs::c_data_structs::ExplicitBitVect;

// Expected values come from Python's rdSubstructLibrary, with numThreads=1
// where the number of threads changes them (bundles, repeated search order
// entries).

const SMILES: [&str; 8] = [
    "CCO",
    "c1ccccc1O",
    "CC(=O)O",
    "c1ccccc1N",
    "CCN",
    "C1CCCCC1O",
    "c1ccncc1",
    "OCCO",
];

fn mols() -> Vec<Mol> {
    SMILES
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut mol = mol_from_smiles(s).unwrap();
            mol.set_prop("_Name", &format!("m{i}"));
            mol
        })
        .collect()
}

fn library(molecules: MolHolderBase) -> SubstructLibrary {
    let mut lib = SubstructLibrary::new(
        molecules,
        Some(PatternHolder::default().into()),
        Some(KeyFromPropHolder::default().into()),
    );
    for mol in mols() {
        lib.add_mol(&mol).unwrap();
    }
    lib
}

fn hydroxyl() -> Mol {
    mol_from_smarts("[OX2H]").unwrap()
}

fn bundle(smarts: &[&str]) -> MolBundle {
    let mut bundle = MolBundle::default();
    for s in smarts {
        bundle.add_mol(&mol_from_smarts(s).unwrap());
    }
    bundle
}

#[test]
fn params_defaults() {
    let params = GetMatchesParams::default();
    assert!(params.parameters.use_chirality);
    assert!(params.parameters.recursion_possible);
    assert!(!params.parameters.use_query_query_matches);
    assert_eq!(params.num_threads, -1);
    assert_eq!(params.max_results, 1000);
    let params = CountMatchesParams::default();
    assert!(params.parameters.use_chirality);
    assert_eq!(params.num_threads, -1);
    let other = SubstructMatchParameters::default();
    assert_eq!(
        SubstructMatchParameters {
            use_chirality: other.use_chirality,
            ..GetMatchesParams::default().parameters
        },
        other
    );
}

#[test]
fn searches() {
    let lib = library(MolHolder::default().into());
    let q = hydroxyl();
    assert_eq!(lib.len(), 8);
    assert_eq!(lib.get_matches(&q).unwrap(), [0, 1, 2, 5, 7]);
    assert_eq!(lib.count_matches(&q).unwrap(), 5);
    assert!(lib.has_match(&q).unwrap());
    let max2 = GetMatchesParams {
        max_results: 2,
        ..GetMatchesParams::default()
    };
    assert_eq!(lib.get_matches_with_params(&q, &max2).unwrap(), [0, 1]);
    assert_eq!(lib.get_matches_range(&q, 2, 6).unwrap(), [2, 5]);
    assert_eq!(lib.count_matches_range(&q, 2, 6).unwrap(), 2);
    assert!(lib.has_match_range(&q, 2, 3).unwrap());
    assert!(lib.has_match_range(&q, 1, 2).unwrap());
    assert_eq!(lib.get_matches_range(&q, 2, 100).unwrap(), [2, 5, 7]);
    assert_eq!(
        lib.get_matches(&mol_from_smarts("[R]").unwrap()).unwrap(),
        [1, 3, 5, 6]
    );

    // Python's SubstructMatchParameters overloads.
    let plain = GetMatchesParams {
        parameters: SubstructMatchParameters::default(),
        ..GetMatchesParams::default()
    };
    assert_eq!(
        lib.get_matches_with_params(&q, &plain).unwrap(),
        [0, 1, 2, 5, 7]
    );
    let one = GetMatchesParams {
        max_results: 1,
        ..plain.clone()
    };
    assert_eq!(
        lib.get_matches_range_with_params(&q, 1, 4, &one).unwrap(),
        [1]
    );
    let mut no_chirality = GetMatchesParams {
        max_results: 3,
        ..GetMatchesParams::default()
    };
    no_chirality.parameters.use_chirality = false;
    assert_eq!(
        lib.get_matches_with_params(&q, &no_chirality).unwrap(),
        [0, 1, 2]
    );

    // The keyword defaults match chirality; SubstructMatchParameters' do not.
    let mut chiral = SubstructLibrary::default();
    chiral
        .add_mol(&mol_from_smiles("C[C@H](N)O").unwrap())
        .unwrap();
    chiral
        .add_mol(&mol_from_smiles("C[C@@H](N)O").unwrap())
        .unwrap();
    let query = mol_from_smiles("C[C@H](N)O").unwrap();
    assert_eq!(chiral.get_matches(&query).unwrap(), [0]);
    assert_eq!(
        chiral.get_matches_with_params(&query, &plain).unwrap(),
        [0, 1]
    );

    // RDKit's preconditions: an empty library, and ranges.
    let empty = SubstructLibrary::default();
    assert!(empty.is_empty());
    assert!(matches!(empty.get_matches(&q), Err(Error::Runtime { .. })));
    assert!(matches!(
        empty.count_matches(&q),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(empty.has_match(&q), Err(Error::Runtime { .. })));
    assert!(matches!(
        lib.get_matches_range(&q, 8, 9),
        Err(Error::Index { index: 8, len: 8 })
    ));
    assert!(matches!(
        lib.get_matches_range(&q, 3, 3),
        Err(Error::Runtime { .. })
    ));

    assert_eq!(
        mol_to_smiles(&lib.get_mol(1).unwrap()).unwrap(),
        "Oc1ccccc1"
    );
    assert!(matches!(lib.get_mol(8), Err(Error::Index { .. })));
}

#[test]
fn search_order() {
    let mut lib = library(MolHolder::default().into());
    let q = hydroxyl();
    lib.set_search_order(&[7, 6, 5, 4, 3, 2, 1, 0]).unwrap();
    assert_eq!(lib.get_search_order(), [7, 6, 5, 4, 3, 2, 1, 0]);
    assert_eq!(lib.get_matches(&q).unwrap(), [7, 5, 2, 1, 0]);
    let max2 = GetMatchesParams {
        max_results: 2,
        ..GetMatchesParams::default()
    };
    assert_eq!(lib.get_matches_with_params(&q, &max2).unwrap(), [7, 5]);
    assert_eq!(lib.get_matches_range(&q, 0, 3).unwrap(), [7, 5]);
    assert_eq!(lib.count_matches(&q).unwrap(), 5);
    // A range starts in the order.
    assert!(matches!(
        lib.get_matches_range(&q, 8, 9),
        Err(Error::Index { .. })
    ));

    // A repeated molecule is matched once.
    lib.set_search_order(&[1, 1, 5]).unwrap();
    assert_eq!(lib.get_matches(&q).unwrap(), [1, 5]);
    assert_eq!(lib.count_matches(&q).unwrap(), 2);

    assert!(matches!(
        lib.set_search_order(&[8]),
        Err(Error::Index { index: 8, len: 8 })
    ));
    assert_eq!(lib.get_search_order(), [1, 1, 5]);
    lib.set_search_order(&[]).unwrap();
    assert!(lib.get_search_order().is_empty());

    // An order the molecule holder no longer fits: RDKit would index its
    // molecules and its bitset of found molecules with it.
    lib.set_search_order(&[7]).unwrap();
    *lib.get_mol_holder_mut() = MolHolder::default().into();
    assert!(matches!(lib.get_matches(&q), Err(Error::Index { .. })));
}

#[test]
fn bundles() {
    let mut lib = library(MolHolder::default().into());
    let two = bundle(&["[NX3H2]", "[OX2H]"]);
    assert_eq!(
        lib.get_matches_mol_bundle(&two).unwrap(),
        [3, 4, 0, 1, 2, 5, 7]
    );
    let max = |max_results| GetMatchesParams {
        max_results,
        ..GetMatchesParams::default()
    };
    assert_eq!(
        lib.get_matches_mol_bundle_with_params(&two, &max(3))
            .unwrap(),
        [3, 4, 0]
    );
    assert_eq!(lib.count_matches_mol_bundle(&two).unwrap(), 7);
    assert!(lib.has_match_mol_bundle(&two).unwrap());

    // RDKit subtracts its running total from maxResults before each member:
    // once a member reaches it, the next finds nothing and the one after
    // finds everything left.
    let three = bundle(&["[NX3H2]", "[OX2H]", "[#7]"]);
    assert_eq!(
        lib.get_matches_mol_bundle(&three).unwrap(),
        [3, 4, 0, 1, 2, 5, 7, 6]
    );
    assert_eq!(
        lib.get_matches_mol_bundle_with_params(&three, &max(2))
            .unwrap(),
        [3, 4, 6]
    );
    // Any negative maxResults is every match; RDKit's subtraction from
    // i32::MIN would overflow.
    assert_eq!(
        lib.get_matches_mol_bundle_with_params(&three, &max(i32::MIN))
            .unwrap(),
        [3, 4, 0, 1, 2, 5, 7, 6]
    );

    lib.set_search_order(&[7, 6, 5, 4, 3, 2, 1, 0]).unwrap();
    assert_eq!(
        lib.get_matches_mol_bundle(&three).unwrap(),
        [4, 3, 7, 5, 2, 1, 0, 6]
    );
    assert_eq!(
        lib.get_matches_mol_bundle_with_params(&three, &max(2))
            .unwrap(),
        [4, 3, 6]
    );
    assert_eq!(
        lib.get_matches_mol_bundle_range(&three, 2, 6).unwrap(),
        [4, 3, 5, 2]
    );
    assert_eq!(lib.count_matches_mol_bundle(&three).unwrap(), 8);
    assert!(lib.has_match_mol_bundle_range(&three, 2, 3).unwrap());
    assert!(
        !lib.has_match_mol_bundle_range_with_params(
            &bundle(&["[Cl]"]),
            0,
            8,
            &HasMatchParams::default()
        )
        .unwrap()
    );
    assert_eq!(
        lib.count_matches_mol_bundle_range_with_params(
            &three,
            0,
            2,
            &CountMatchesParams::default()
        )
        .unwrap(),
        2
    );

    // An empty bundle searches nothing, as RDKit's loop does.
    assert!(
        SubstructLibrary::default()
            .get_matches_mol_bundle(&MolBundle::default())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn keys() {
    let lib = library(MolHolder::default().into());
    let keys = lib.get_key_holder().unwrap();
    let matches = lib.get_matches(&hydroxyl()).unwrap();
    assert_eq!(
        keys.get_keys(&matches).unwrap(),
        ["m0", "m1", "m2", "m5", "m7"]
    );
    assert_eq!(keys.get_key(0).unwrap(), "m0");
    assert_eq!(keys.len(), 8);
    let props = keys.downcast_ref::<KeyFromPropHolder>().unwrap();
    assert_eq!(props.get_prop_name().unwrap(), "_Name");
    assert!(matches!(keys.get_keys(&[99]), Err(Error::Index { .. })));
    assert!(matches!(keys.get_key(8), Err(Error::Index { .. })));

    let mut holder = KeyFromPropHolder::new("foo");
    assert_eq!(holder.get_prop_name().unwrap(), "foo");
    assert_eq!(holder.add_mol(&mol_from_smiles("C").unwrap()).unwrap(), 0);
    assert_eq!(holder.get_key(0).unwrap(), "LIBIDX-0");
    assert_eq!(holder.add_key("x").unwrap(), 1);
    assert_eq!(holder.get_keys(&[1, 0]).unwrap(), ["x", "LIBIDX-0"]);
}

#[test]
fn holders() {
    let q = hydroxyl();
    let ring = mol_from_smarts("[R]").unwrap();
    for holder in [
        MolHolderBase::from(CachedMolHolder::default()),
        CachedSmilesMolHolder::default().into(),
        CachedTrustedSmilesMolHolder::default().into(),
    ] {
        let lib = library(holder);
        assert_eq!(lib.get_matches(&q).unwrap(), [0, 1, 2, 5, 7]);
        assert_eq!(lib.get_matches(&ring).unwrap(), [1, 3, 5, 6]);
        assert_eq!(
            mol_to_smiles(&lib.get_mol(1).unwrap()).unwrap(),
            "Oc1ccccc1"
        );
        assert_eq!(lib.get_mol_holder().len(), 8);
    }

    let mut lib = SubstructLibrary::new(CachedSmilesMolHolder::default().into(), None, None);
    for mol in mols() {
        lib.add_mol(&mol).unwrap();
    }
    assert!(lib.get_mol_holder().downcast_ref::<MolHolder>().is_none());
    let smiles = lib
        .get_mol_holder_mut()
        .downcast_mut::<CachedSmilesMolHolder>()
        .unwrap();
    assert_eq!(smiles.add_smiles("garbage").unwrap(), 8);
    // A stored SMILES that does not parse is skipped by searches; Python's
    // GetMol returns None.
    assert!(matches!(lib.get_mol(8), Err(Error::Parse { .. })));
    assert_eq!(lib.get_matches(&q).unwrap(), [0, 1, 2, 5, 7]);
    assert!(format!("{:?}", lib.get_mol_holder()).contains("CachedSmilesMolHolder"));

    let mut cached = CachedMolHolder::default();
    // SAFETY: a pickle RDKit wrote.
    assert_eq!(
        unsafe { cached.add_binary(&mols()[0].to_binary()) }.unwrap(),
        0
    );
    assert_eq!(mol_to_smiles(&cached.get_mol(0).unwrap()).unwrap(), "CCO");
    // SAFETY: a pickle RDKit wrote.
    assert_eq!(unsafe { cached.add_binary(b"not a pickle") }.unwrap(), 1);
    assert!(cached.get_mol(1).is_err());

    let mut trusted = CachedTrustedSmilesMolHolder::default();
    trusted.add_smiles("C1CC1").unwrap();
    assert_eq!(trusted.get_mol(0).unwrap().get_num_atoms(), 3);
}

#[test]
fn fingerprints() {
    let lib = library(MolHolder::default().into());
    let fps = lib.get_fp_holder().unwrap();
    assert_eq!(fps.len(), 8);
    assert_eq!(fps.get_num_bits(), 2048);
    let fp = fps.get_fingerprint(0).unwrap();
    assert_eq!(fp.get_num_bits(), 2048);
    assert_eq!(fp.get_num_on_bits(), 16);
    assert!(
        fps.passes_filter(0, &fps.make_fingerprint(&hydroxyl()).unwrap())
            .unwrap()
    );
    let ethanol = fps
        .make_fingerprint(&mol_from_smiles("CCO").unwrap())
        .unwrap();
    assert!(!fps.passes_filter(6, &ethanol).unwrap());
    assert!(matches!(fps.get_fingerprint(8), Err(Error::Index { .. })));
    assert!(matches!(
        fps.passes_filter(8, &ethanol),
        Err(Error::Index { .. })
    ));
    assert!(fps.downcast_ref::<PatternHolder>().is_some());
    assert!(fps.downcast_ref::<TautomerPatternHolder>().is_none());
    assert_eq!(
        PatternHolder::new(1024)
            .make_fingerprint(&mols()[1])
            .unwrap()
            .get_num_bits(),
        1024
    );

    let mut tautomers = SubstructLibrary::new(
        MolHolder::default().into(),
        Some(TautomerPatternHolder::default().into()),
        None,
    );
    for s in ["Oc1ccccn1", "O=c1cccc[nH]1", "CCO"] {
        tautomers.add_mol(&mol_from_smiles(s).unwrap()).unwrap();
    }
    let query = mol_from_smiles("O=c1cccc[nH]1").unwrap();
    assert_eq!(tautomers.get_matches(&query).unwrap(), [1]);
    let fps = tautomers.get_fp_holder().unwrap();
    assert_eq!(fps.get_fingerprint(0).unwrap().get_num_on_bits(), 117);
}

// AllProbeBitsMatch reads the query's blocks from the stored fingerprint,
// past the end of a shorter one, so sizes must agree.
#[test]
fn fingerprint_sizes() {
    let mut holder = PatternHolder::new(1 << 22);
    assert!(matches!(
        holder.add_fingerprint(&ExplicitBitVect::new(1)),
        Err(Error::Value { .. })
    ));
    let fp = holder.make_fingerprint(&mols()[0]).unwrap();
    assert_eq!(holder.add_fingerprint(&fp).unwrap(), 0);
    assert!(holder.passes_filter(0, &fp).unwrap());
    assert!(matches!(
        holder.passes_filter(0, &ExplicitBitVect::new(1)),
        Err(Error::Value { .. })
    ));
    // A holder of no bits makes no fingerprints (RDKit's precondition).
    let mut empty = PatternHolder::new(0);
    assert!(matches!(
        empty.add_mol(&mols()[0]),
        Err(Error::Runtime { .. })
    ));
}

#[test]
#[should_panic(expected = "at most i32::MAX bits")]
fn pattern_holders_are_int_sized() {
    let _ = PatternHolder::new(u32::MAX);
}

// Python's searches with several threads abort with heap corruption when
// one thread throws, as reading a stored SMILES that does not sanitize
// does, and segfault when the library has fewer fingerprints than
// molecules; these search on the calling thread and fail.
#[test]
fn errors_in_searches() {
    let mut holder = CachedSmilesMolHolder::default();
    for i in 0..200 {
        holder
            .add_smiles(if i % 7 == 0 { "c1cccc1" } else { "c1ccccc1CC" })
            .unwrap();
    }
    let lib = SubstructLibrary::new(holder.into(), None, None);
    let params = GetMatchesParams {
        num_threads: 8,
        max_results: -1,
        ..GetMatchesParams::default()
    };
    let query = mol_from_smarts("cC").unwrap();
    assert!(matches!(
        lib.get_matches_with_params(&query, &params),
        Err(Error::Kekulize { .. })
    ));

    let mut fps = PatternHolder::default();
    let mut molecules = MolHolder::default();
    let ethanol = mol_from_smiles("CCO").unwrap();
    for i in 0..50 {
        molecules.add_mol(&ethanol).unwrap();
        if i < 25 {
            fps.add_mol(&ethanol).unwrap();
        }
    }
    let lib = SubstructLibrary::new(molecules.into(), Some(fps.into()), None);
    assert!(
        lib.get_matches_with_params(&mol_from_smiles("CC").unwrap(), &params)
            .is_err()
    );
    // A library whose holders differ in size adds to the molecules first.
    let mut lib = lib;
    assert!(matches!(lib.add_mol(&ethanol), Err(Error::Runtime { .. })));
    assert_eq!(lib.len(), 51);

    // RDKit negates numThreads in int.
    let mut small = SubstructLibrary::default();
    small.add_mol(&ethanol).unwrap();
    let min = GetMatchesParams {
        num_threads: i32::MIN,
        ..GetMatchesParams::default()
    };
    assert_eq!(
        small
            .get_matches_with_params(&mol_from_smiles("C").unwrap(), &min)
            .unwrap(),
        [0]
    );
}

#[test]
fn patterns() {
    let mut lib = SubstructLibrary::default();
    for mol in mols() {
        lib.add_mol(&mol).unwrap();
    }
    add_patterns(&mut lib).unwrap();
    assert_eq!(lib.get_fp_holder().unwrap().len(), 8);
    assert_eq!(lib.get_matches(&hydroxyl()).unwrap(), [0, 1, 2, 5, 7]);
    assert!(matches!(add_patterns(&mut lib), Err(Error::Runtime { .. })));

    let mut lib = SubstructLibrary::default();
    for mol in mols() {
        lib.add_mol(&mol).unwrap();
    }
    let mut holder = TautomerPatternHolder::new(512);
    holder.add_mol(&mols()[0]).unwrap();
    add_patterns_holder_with_num_threads(&mut lib, holder.into(), 4).unwrap();
    let fps = lib.get_fp_holder().unwrap();
    assert!(fps.downcast_ref::<TautomerPatternHolder>().is_some());
    // The holder's fingerprint is replaced.
    assert_eq!(fps.len(), 8);
    assert_eq!(fps.get_fingerprint(0).unwrap().get_num_bits(), 512);

    // RDKit's AddPatterns leaves null fingerprints behind when a molecule
    // fails (Python segfaults reading one); here the library is unchanged.
    let mut smiles = CachedSmilesMolHolder::default();
    smiles.add_smiles("c1cccc1").unwrap();
    let mut lib = SubstructLibrary::new(smiles.into(), None, None);
    assert!(matches!(
        add_patterns_with_num_threads(&mut lib, 2),
        Err(Error::Kekulize { .. })
    ));
    assert!(lib.get_fp_holder().is_none());
}

#[test]
fn pickles() {
    assert!(substruct_library_can_serialize());
    let mut lib = library(MolHolder::default().into());
    lib.set_search_order(&[3, 1]).unwrap();
    let pickle = lib.serialize().unwrap();
    assert!(pickle.starts_with(b"22 serialization::archive "));
    // SAFETY: RDKit wrote it.
    let copy = unsafe { SubstructLibrary::from_pickle(&pickle) }.unwrap();
    assert_eq!(copy.len(), 8);
    assert_eq!(copy.get_search_order(), [3, 1]);
    assert_eq!(copy.get_matches(&hydroxyl()).unwrap(), [1]);
    assert!(copy.get_mol_holder().downcast_ref::<MolHolder>().is_some());
    assert!(
        copy.get_fp_holder()
            .unwrap()
            .downcast_ref::<PatternHolder>()
            .is_some()
    );
    assert_eq!(copy.get_key_holder().unwrap().get_key(3).unwrap(), "m3");

    let mut stream = Vec::new();
    lib.to_stream(&mut stream).unwrap();
    assert_eq!(stream, pickle);
    let mut other = SubstructLibrary::default();
    // SAFETY: RDKit wrote it.
    unsafe { other.init_from_stream(&mut stream.as_slice()) }.unwrap();
    assert_eq!(other.len(), 8);

    for holder in [
        MolHolderBase::from(CachedMolHolder::default()),
        CachedSmilesMolHolder::default().into(),
        CachedTrustedSmilesMolHolder::default().into(),
    ] {
        let lib = SubstructLibrary::new(holder, None, None);
        let mut lib = lib;
        for mol in mols() {
            lib.add_mol(&mol).unwrap();
        }
        // SAFETY: RDKit wrote it.
        let copy = unsafe { SubstructLibrary::from_pickle(&lib.serialize().unwrap()) }.unwrap();
        assert_eq!(copy.get_matches(&hydroxyl()).unwrap(), [0, 1, 2, 5, 7]);
        assert!(copy.get_fp_holder().is_none());
    }
}

// A search order past the molecules would index them, and RDKit's bitset of
// found molecules, out of range.
#[test]
fn pickled_search_orders_are_checked() {
    let mut lib = SubstructLibrary::default();
    lib.add_mol(&mols()[0]).unwrap();
    lib.set_search_order(&[0]).unwrap();
    let mut pickle = lib.serialize().unwrap();
    // The order is the archive's first vector, after its header and the
    // library's: its count, its item version, and its entry.
    let at = pickle.windows(7).position(|w| w == b" 1 0 0 ").unwrap();
    pickle[at + 5] = b'7';
    // SAFETY: the archive is the one RDKit wrote with another index, so
    // Boost reads the same objects.
    let result = unsafe { SubstructLibrary::from_pickle(&pickle) };
    assert!(matches!(result, Err(Error::Index { index: 7, len: 1 })));
}

fn tautomer_library(fingerprints: Option<FPHolderBase>, smiles: &[&str]) -> SubstructLibrary {
    let mut lib = SubstructLibrary::new(MolHolder::default().into(), fingerprints, None);
    for s in smiles {
        lib.add_mol(&mol_from_smiles(s).unwrap()).unwrap();
    }
    lib
}

const PYRIDONES: [&str; 5] = ["Oc1ccccn1", "O=c1cccc[nH]1", "CCO", "Oc1ccccc1", "CC(=O)C"];

fn tautomer_query(smiles: &str) -> TautomerQuery {
    TautomerQuery::new(&mol_from_smiles(smiles).unwrap()).unwrap()
}

fn extended_query(smiles: &str) -> ExtendedQueryMol {
    create_extended_query_mol(&mol_from_smiles(smiles).unwrap()).unwrap()
}

#[test]
fn tautomer_and_extended_queries() {
    let lib = tautomer_library(Some(TautomerPatternHolder::default().into()), &PYRIDONES);
    let pyridone = tautomer_query("O=c1cccc[nH]1");
    assert_eq!(lib.get_matches_tautomer_query(&pyridone).unwrap(), [0, 1]);
    assert_eq!(lib.count_matches_tautomer_query(&pyridone).unwrap(), 2);
    assert!(lib.has_match_tautomer_query(&pyridone).unwrap());
    assert_eq!(
        lib.get_matches_tautomer_query_range(&pyridone, 1, 3)
            .unwrap(),
        [1]
    );
    let one = GetMatchesParams {
        max_results: 1,
        ..Default::default()
    };
    assert_eq!(
        lib.get_matches_tautomer_query_with_params(&pyridone, &one)
            .unwrap(),
        [0]
    );
    let acetone = tautomer_query("CC(=O)C");
    assert_eq!(lib.get_matches_tautomer_query(&acetone).unwrap(), [4]);
    assert_eq!(lib.count_matches_tautomer_query(&acetone).unwrap(), 1);

    let extended = extended_query("O=c1cccc[nH]1");
    assert_eq!(
        lib.get_matches_extended_query_mol(&extended).unwrap(),
        [0, 1]
    );
    assert_eq!(lib.count_matches_extended_query_mol(&extended).unwrap(), 2);
    assert!(
        !lib.has_match_extended_query_mol_range(&extended, 2, 3)
            .unwrap()
    );
    let ethanol = create_extended_query_mol_with_params(
        &mol_from_smiles("CCO").unwrap(),
        &CreateExtendedQueryMolParams {
            do_tautomers: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(lib.get_matches_extended_query_mol(&ethanol).unwrap(), [2]);
    let benzene = create_extended_query_mol(&mol_from_smarts("c1ccccc1").unwrap()).unwrap();
    assert_eq!(lib.get_matches_extended_query_mol(&benzene).unwrap(), [3]);

    // Without tautomer fingerprints RDKit searches every molecule.
    for fingerprints in [None, Some(PatternHolder::default().into())] {
        let lib = tautomer_library(fingerprints, &PYRIDONES);
        assert_eq!(lib.get_matches_tautomer_query(&pyridone).unwrap(), [0, 1]);
        assert_eq!(
            lib.get_matches_extended_query_mol(&extended).unwrap(),
            [0, 1]
        );
    }

    // A recursive query in a template whose query modifies no atoms or
    // bonds is never read past.
    for (smarts, matches) in [("[$(O=c)]=c1cccc[nH]1", [1]), ("[$([OH])]c1ccccn1", [0])] {
        let query = TautomerQuery::new(&mol_from_smarts(smarts).unwrap()).unwrap();
        assert!(query.get_modified_atoms().is_empty());
        assert_eq!(lib.get_matches_tautomer_query(&query).unwrap(), matches);
        assert_eq!(lib.count_matches_tautomer_query(&query).unwrap(), 1);
    }

    let handles: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                let lib =
                    tautomer_library(Some(TautomerPatternHolder::default().into()), &PYRIDONES);
                let query = tautomer_query("O=c1cccc[nH]1");
                let extended = extended_query("O=c1cccc[nH]1");
                (0..10)
                    .map(|_| {
                        (
                            lib.get_matches_tautomer_query(&query).unwrap(),
                            lib.count_matches_extended_query_mol(&extended).unwrap(),
                        )
                    })
                    .collect::<Vec<_>>()
            })
        })
        .collect();
    for handle in handles {
        for (matches, count) in handle.join().unwrap() {
            assert_eq!(matches, [0, 1]);
            assert_eq!(count, 2);
        }
    }
}

// The template with its first atom replaced by the SMARTS.
fn with_first_atom(smarts: &str) -> Mol {
    let mut mol = RWMol::from(mol_from_smarts(smarts).unwrap());
    mol.insert_mol(&mol_from_smiles("CCC(=O)C").unwrap());
    mol.add_bond_with_order(0, 1, BondType::Single).unwrap();
    mol.update_property_cache_with_strict(false).unwrap();
    mol.into()
}

#[test]
fn recursive_tautomer_templates() {
    let lib = tautomer_library(
        Some(TautomerPatternHolder::default().into()),
        &[
            "CCCC(=O)C",
            "CCCC(O)=C",
            "OCCC(=O)C",
            "CCC=C(O)C",
            "CC(C)CC(=O)C",
        ],
    );
    let plain = with_first_atom("[CH3]");
    let query = TautomerQuery::new(&plain).unwrap();
    assert_eq!(
        lib.get_matches_tautomer_query(&query).unwrap(),
        [0, 1, 3, 4]
    );
    let extended = create_extended_query_mol(&plain).unwrap();
    assert_eq!(
        lib.get_matches_extended_query_mol(&extended).unwrap(),
        [0, 1, 3, 4]
    );

    // RDKit's tautomer check reads the recursive SMARTS's one-atom matches
    // at the modified atoms 2, 4, and 5: Python's GetMatches reads atom
    // 2532319792 of six for the query, and finds nothing for the extended
    // query.
    let recursive = with_first_atom("[$([CH3])]");
    let query = TautomerQuery::new(&recursive).unwrap();
    assert_eq!(query.get_modified_atoms(), [2, 4, 5]);
    assert!(matches!(
        lib.get_matches_tautomer_query(&query),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        lib.count_matches_tautomer_query_range(&query, 0, 2),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        lib.has_match_tautomer_query(&query),
        Err(Error::Value { .. })
    ));
    let extended = create_extended_query_mol(&recursive).unwrap();
    assert!(matches!(
        lib.get_matches_extended_query_mol(&extended),
        Err(Error::Value { .. })
    ));
    // Without tautomers the extended query is a plain one.
    let extended = create_extended_query_mol_with_params(
        &recursive,
        &CreateExtendedQueryMolParams {
            do_tautomers: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        lib.get_matches_extended_query_mol(&extended).unwrap(),
        [0, 4]
    );
}

// Searches match under control_c_mutex; libraries searched on several
// threads at once agree.
#[test]
fn threads() {
    let handles: Vec<_> = (0..4)
        .map(|i| {
            std::thread::spawn(move || {
                let holder: MolHolderBase = if i % 2 == 0 {
                    MolHolder::default().into()
                } else {
                    CachedSmilesMolHolder::default().into()
                };
                let lib = library(holder);
                (0..10)
                    .map(|_| {
                        (
                            lib.get_matches(&hydroxyl()).unwrap(),
                            lib.count_matches_mol_bundle(&bundle(&["[NX3H2]", "[OX2H]"]))
                                .unwrap(),
                        )
                    })
                    .collect::<Vec<_>>()
            })
        })
        .collect();
    for handle in handles {
        for (matches, count) in handle.join().unwrap() {
            assert_eq!(matches, [0, 1, 2, 5, 7]);
            assert_eq!(count, 7);
        }
    }
}
