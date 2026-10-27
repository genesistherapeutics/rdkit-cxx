// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdChemReactions' EnumerateLibrary and enumeration strategies: the
//! expected values come from RDKit's Python bindings.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chem_reactions::*;
use rdkit_cxx::chem::{
    AtomMonomerInfo, AtomMonomerInfoParams, AtomMonomerType, Mol, mol_from_smiles, mol_to_smiles,
};

const AMIDE: &str = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]";

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn amide() -> ChemicalReaction {
    reaction_from_smarts(AMIDE).unwrap()
}

fn group(texts: &[&str]) -> Vec<Mol> {
    texts.iter().map(|text| smiles(text)).collect()
}

/// Three acids and an alcohol the acid template does not match, and two
/// amines.
fn building_blocks() -> Vec<Vec<Mol>> {
    vec![
        group(&["CC(=O)O", "CCC(=O)O", "OC(=O)c1ccccc1", "CCO"]),
        group(&["CN", "CCN"]),
    ]
}

/// The matching building blocks only.
fn matching() -> Vec<Vec<Mol>> {
    vec![
        group(&["CC(=O)O", "CCC(=O)O", "OC(=O)c1ccccc1"]),
        group(&["CN", "CCN"]),
    ]
}

fn group_smiles(groups: &[Vec<Mol>]) -> Vec<Vec<String>> {
    groups
        .iter()
        .map(|set| set.iter().map(|mol| mol_to_smiles(mol).unwrap()).collect())
        .collect()
}

fn walk(strategy: &mut EnumerationStrategyBase, n: usize) -> Vec<Vec<u64>> {
    (0..n).map(|_| strategy.next().unwrap()).collect()
}

fn is_value(result: Result<impl std::fmt::Debug, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

#[test]
fn library_enumerates_every_combination() {
    assert!(enumerate_library_can_serialize());
    let mut lib = EnumerateLibrary::new(&amide(), &building_blocks()).unwrap();
    assert_eq!(
        group_smiles(&lib.get_reagents()),
        [
            vec!["CC(=O)O", "CCC(=O)O", "O=C(O)c1ccccc1"],
            vec!["CN", "CCN"]
        ]
    );
    assert_eq!(lib.get_position().unwrap(), [0, 0]);
    assert!(lib.has_next().unwrap());
    let mut all = Vec::new();
    let mut positions = Vec::new();
    while lib.has_next().unwrap() {
        all.push(lib.next_smiles().unwrap());
        positions.push(lib.get_position().unwrap());
    }
    let expected = [
        "CNC(C)=O",
        "CCC(=O)NC",
        "CNC(=O)c1ccccc1",
        "CCNC(C)=O",
        "CCNC(=O)CC",
        "CCNC(=O)c1ccccc1",
    ];
    assert_eq!(
        all,
        expected
            .iter()
            .map(|s| vec![vec![s.to_string()]])
            .collect::<Vec<_>>()
    );
    assert_eq!(positions, [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]);
    // Past the last: RDKit's PRECONDITION.
    assert!(matches!(lib.next_smiles(), Err(Error::Runtime { .. })));
    assert!(lib.next().is_none());

    lib.reset_state().unwrap();
    assert_eq!(lib.get_position().unwrap(), [0, 0]);
    assert!(lib.has_next().unwrap());
    let products: Vec<_> = lib.by_ref().map(Result::unwrap).collect();
    assert_eq!(products.len(), 6);
    assert_eq!(mol_to_smiles(&products[0][0][0]).unwrap(), "CNC(C)=O");

    let enumerator = lib.get_enumerator().unwrap();
    assert_eq!(enumerator.type_(), "CartesianProductStrategy");
    assert_eq!(enumerator.get_num_permutations(), 6);
    assert_eq!(enumerator.get_position(), [2, 1]);
    assert_eq!(lib.get_reaction().get_num_reactant_templates(), 2);
}

#[test]
fn library_owns_its_building_blocks() {
    let lib = {
        let blocks = building_blocks();
        EnumerateLibrary::new(&amide(), &blocks).unwrap()
    };
    let mut lib = lib;
    assert_eq!(lib.next_smiles().unwrap(), [vec!["CNC(C)=O".to_string()]]);
}

#[test]
fn library_params_filter_building_blocks() {
    let defaults = EnumerationParams::default();
    assert_eq!(defaults.reagent_max_match_count, i32::MAX);
    assert!(!defaults.sane_partial_products);

    let mut blocks = building_blocks();
    blocks[0].push(smiles("OC(=O)CC(=O)O"));
    let lib = EnumerateLibrary::new(&amide(), &blocks).unwrap();
    assert_eq!(
        group_smiles(&lib.get_reagents())[0],
        ["CC(=O)O", "CCC(=O)O", "O=C(O)c1ccccc1", "O=C(O)CC(=O)O"]
    );
    let params = EnumerationParams {
        reagent_max_match_count: 1,
        ..EnumerationParams::default()
    };
    let lib = EnumerateLibrary::new_with_params(&amide(), &blocks, &params).unwrap();
    assert_eq!(
        group_smiles(&lib.get_reagents())[0],
        ["CC(=O)O", "CCC(=O)O", "O=C(O)c1ccccc1"]
    );
}

#[test]
fn library_with_strategies() {
    let mut random = EnumerateLibrary::new_with_enumerator(
        &amide(),
        &building_blocks(),
        &RandomSampleStrategy::default(),
    )
    .unwrap();
    let first: Vec<_> = (0..3).map(|_| random.next_smiles().unwrap()).collect();
    assert_eq!(
        first,
        [
            vec![vec!["CNC(C)=O".to_string()]],
            vec![vec!["CCNC(=O)CC".to_string()]],
            vec![vec!["CNC(=O)c1ccccc1".to_string()]],
        ]
    );
    assert_eq!(
        random.get_enumerator().unwrap().type_(),
        "RandomSampleStrategy"
    );

    let mut even = EnumerateLibrary::new_with_enumerator_and_params(
        &amide(),
        &building_blocks(),
        &EvenSamplePairsStrategy::default(),
        &EnumerationParams::default(),
    )
    .unwrap();
    assert_eq!(even.next_smiles().unwrap(), [vec!["CCC(=O)NC".to_string()]]);
    assert_eq!(even.get_position().unwrap(), [1, 0]);
}

#[test]
fn cartesian_strategy() {
    let mut strategy = CartesianProductStrategy::default();
    assert_eq!(strategy.type_(), "CartesianProductStrategy");
    strategy.initialize(&amide(), &matching()).unwrap();
    assert_eq!(strategy.get_num_permutations(), 6);
    assert_eq!(strategy.get_position(), [0, 0]);
    let mut seq = Vec::new();
    while strategy.has_next() {
        seq.push(strategy.next().unwrap());
    }
    assert_eq!(seq, [[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]]);

    let mut skipped = CartesianProductStrategy::default();
    skipped.initialize(&amide(), &matching()).unwrap();
    assert!(skipped.skip(3).unwrap());
    assert_eq!(skipped.get_position(), [2, 0]);
    let copy = skipped.clone();
    assert_eq!(copy.get_position(), [2, 0]);
    let base: EnumerationStrategyBase = copy.into();
    assert_eq!(base.type_(), "CartesianProductStrategy");
}

#[test]
fn random_strategies() {
    let mut random = RandomSampleStrategy::default();
    random.initialize(&amide(), &matching()).unwrap();
    assert!(random.has_next());
    assert_eq!(
        walk(&mut random, 8),
        [
            [0, 0],
            [1, 1],
            [2, 0],
            [1, 0],
            [0, 1],
            [0, 1],
            [1, 1],
            [1, 1]
        ]
    );
    let mut all_bbs = RandomSampleAllBBsStrategy::default();
    all_bbs.initialize(&amide(), &matching()).unwrap();
    assert_eq!(all_bbs.type_(), "RandomSampleAllBBsStrategy");
    assert_eq!(
        walk(&mut all_bbs, 8),
        [
            [1, 1],
            [2, 0],
            [0, 1],
            [0, 0],
            [1, 1],
            [2, 0],
            [0, 1],
            [1, 1]
        ]
    );
}

#[test]
fn even_sample_pairs_strategy() {
    let mut even = EvenSamplePairsStrategy::default();
    even.initialize(&amide(), &matching()).unwrap();
    assert_eq!(even.type_(), "EvenSamplePairsStrategy");
    assert_eq!(even.get_num_permutations(), 6);
    assert_eq!(
        walk(&mut even, 6),
        [[1, 0], [0, 1], [2, 0], [1, 1], [0, 0], [2, 1]]
    );
    assert_eq!(
        even.stats(),
        "#BEGIN# BBSTAT\n1\t3\t  2.00\t0|1\t0|2\t0|3\n2\t2\t  3.00\t0|1\t0|2\n\
         #END# BBSTAT\n#BEGIN# PAIRSTAT\n1\t2\t3\t2\t  1.00\t0|6\n#END# PAIRSTAT\n\
         Rejected Period: 1\nRejected (dupes): 0\nRejected Slack Conditions: 0\n\
         Rejected Pair Sampling: 0\n"
    );
    // Its next throws past the last permutation, and before any.
    assert!(matches!(even.next(), Err(Error::Runtime { .. })));
    let mut fresh = EvenSamplePairsStrategy::default();
    assert!(matches!(fresh.next(), Err(Error::Runtime { .. })));
}

#[test]
fn state_round_trip() {
    let mut lib = EnumerateLibrary::new(&amide(), &building_blocks()).unwrap();
    lib.next_smiles().unwrap();
    lib.next_smiles().unwrap();
    let state = lib.get_state().unwrap();
    let next = lib.next_smiles().unwrap();
    // SAFETY: the state came from get_state.
    unsafe { lib.set_state(&state) }.unwrap();
    assert_eq!(lib.next_smiles().unwrap(), next);
    assert!(state.contains("RDKit::CartesianProductStrategy"));
}

#[test]
fn serialize_round_trip() {
    for strategy in [
        EnumerationStrategyBase::from(CartesianProductStrategy::default()),
        RandomSampleStrategy::default().into(),
        RandomSampleAllBBsStrategy::default().into(),
        EvenSamplePairsStrategy::default().into(),
    ] {
        let mut lib =
            EnumerateLibrary::new_with_enumerator(&amide(), &building_blocks(), &strategy).unwrap();
        lib.next_smiles().unwrap();
        let bytes = lib.serialize().unwrap();
        let mut loaded = EnumerateLibrary::default();
        // SAFETY: the bytes came from serialize.
        unsafe { loaded.init_from_string(&bytes) }.unwrap();
        assert_eq!(loaded.get_position().unwrap(), lib.get_position().unwrap());
        assert_eq!(loaded.get_enumerator().unwrap().type_(), strategy.type_());
        assert_eq!(
            group_smiles(&loaded.get_reagents()),
            group_smiles(&lib.get_reagents())
        );
        assert_eq!(loaded.next_smiles().unwrap(), lib.next_smiles().unwrap());
        loaded.reset_state().unwrap();
        lib.reset_state().unwrap();
        assert_eq!(loaded.next_smiles().unwrap(), lib.next_smiles().unwrap());

        // A library that has a reaction: as Python's.
        // SAFETY: the bytes came from serialize.
        assert!(is_value(unsafe { loaded.init_from_string(&bytes) }));
    }
}

/// RDKit dereferences an EnumerateLibrary()'s missing strategy (Python
/// segfaults on GetPosition); the wrappers refuse.
#[test]
fn empty_library() {
    let mut lib = EnumerateLibrary::default();
    assert!(is_value(lib.get_position()));
    assert!(is_value(lib.has_next()));
    assert!(is_value(lib.get_enumerator()));
    assert!(is_value(lib.get_state()));
    assert!(is_value(lib.next_smiles()));
    assert!(lib.next().is_none());
    assert!(lib.get_reagents().is_empty());
    assert_eq!(lib.get_reaction().get_num_reactant_templates(), 0);
    // RDKit's PRECONDITION, as Python's RuntimeError.
    assert!(matches!(lib.reset_state(), Err(Error::Runtime { .. })));
    assert!(format!("{lib:?}").contains("None"));

    // An empty library serializes, and loads as one.
    let bytes = lib.serialize().unwrap();
    let mut loaded = EnumerateLibrary::default();
    // SAFETY: the bytes came from serialize.
    unsafe { loaded.init_from_string(&bytes) }.unwrap();
    assert!(is_value(loaded.get_position()));
}

/// RDKit's next indexes the building blocks with the strategy's position
/// unchecked, and Python's SetState installs a state from any library: one
/// from a larger library reads past the building blocks.
#[test]
fn state_from_another_library_is_refused() {
    let mut big = EnumerateLibrary::new(&amide(), &building_blocks()).unwrap();
    for _ in 0..5 {
        big.next_smiles().unwrap();
    }
    let state = big.get_state().unwrap();
    let small_blocks = vec![group(&["CC(=O)O"]), group(&["CN"])];
    let mut small = EnumerateLibrary::new(&amide(), &small_blocks).unwrap();
    // SAFETY: the state came from get_state.
    assert!(is_value(unsafe { small.set_state(&state) }));
    // The library keeps its own state.
    assert_eq!(small.get_position().unwrap(), [0, 0]);
    assert_eq!(small.next_smiles().unwrap(), [vec!["CNC(C)=O".to_string()]]);
}

/// The random strategies draw from uniform_int_distribution(0, size - 1),
/// which an empty group makes empty (boost does not allow it), and
/// RandomSampleAllBBsStrategy's initialize dereferences the largest group's
/// size, which no groups do not have.
#[test]
fn random_strategies_need_building_blocks() {
    let empty_group = vec![group(&["CC(=O)O"]), vec![]];
    let mut random = RandomSampleStrategy::default();
    assert!(is_value(random.initialize(&amide(), &empty_group)));
    let mut all_bbs = RandomSampleAllBBsStrategy::default();
    assert!(is_value(all_bbs.initialize(&amide(), &empty_group)));
    assert!(is_value(all_bbs.initialize(&amide(), &[])));
    assert!(is_value(all_bbs.initialize(&amide(), &[vec![]])));
    // The other strategies take them.
    let mut cartesian = CartesianProductStrategy::default();
    cartesian.initialize(&amide(), &empty_group).unwrap();
    assert_eq!(cartesian.get_num_permutations(), 0);
    assert!(!cartesian.has_next());
    let mut even = EvenSamplePairsStrategy::default();
    even.initialize(&amide(), &empty_group).unwrap();
    assert!(matches!(even.next(), Err(Error::Runtime { .. })));

    // A library checks the groups left once the templates have filtered
    // them.
    let unmatched = vec![group(&["CC(=O)O"]), group(&["CCO"])];
    assert!(is_value(EnumerateLibrary::new_with_enumerator(
        &amide(),
        &unmatched,
        &RandomSampleStrategy::default()
    )));
    assert!(is_value(EnumerateLibrary::new_with_enumerator(
        &amide(),
        &unmatched,
        &RandomSampleAllBBsStrategy::default()
    )));
    let mut lib = EnumerateLibrary::new(&amide(), &unmatched).unwrap();
    assert!(!lib.has_next().unwrap());
    assert!(lib.next().is_none());
}

/// EvenSamplePairsStrategy doubles a uint64 modulus until it reaches the
/// number of permutations, which never ends past 2^63; past 2^64, RDKit's
/// PRECONDITION throws.
#[test]
fn even_sample_pairs_needs_at_most_2_pow_63_permutations() {
    let pair = || group(&["C", "N"]);
    let mut groups: Vec<Vec<Mol>> = (0..62).map(|_| pair()).collect();
    groups.push(group(&["C", "N", "O"]));
    let mut even = EvenSamplePairsStrategy::default();
    assert!(is_value(even.initialize(&amide(), &groups)));

    let groups: Vec<Vec<Mol>> = (0..65).map(|_| pair()).collect();
    assert!(matches!(
        even.initialize(&amide(), &groups),
        Err(Error::Runtime { .. })
    ));

    let groups: Vec<Vec<Mol>> = (0..63).map(|_| pair()).collect();
    even.initialize(&amide(), &groups).unwrap();
    assert_eq!(even.get_num_permutations(), 1 << 63);
}

#[test]
fn libraries_and_strategies_move_between_threads() {
    let lib = EnumerateLibrary::new(&amide(), &building_blocks()).unwrap();
    let strategy = EvenSamplePairsStrategy::default();
    let (mut lib, strategy) = std::thread::spawn(move || (lib, strategy)).join().unwrap();
    assert_eq!(lib.next_smiles().unwrap(), [vec!["CNC(C)=O".to_string()]]);
    assert_eq!(strategy.type_(), "EvenSamplePairsStrategy");
}

/// A library's reaction and building blocks are molecule pickles, whose
/// monomer infos of type PDBRESIDUE RDKit 2024.09's pickler reads as
/// AtomPDBResidueInfo, with its residue fields, so they cannot give plain
/// info that type (U10): the building block's info comes back as PDB residue
/// info, as Python's InitFromString gives it.
#[test]
fn loaded_libraries_read_the_pdb_type_as_pdb_residue_info() {
    let mut amine = smiles("CN");
    let params = AtomMonomerInfoParams {
        name: "MARK".into(),
        ..Default::default()
    };
    amine.get_atom_with_idx_mut(0).unwrap().set_monomer_info(
        &AtomMonomerInfo::new_with_params(AtomMonomerType::Other, &params).unwrap(),
    );
    let blocks = vec![group(&["CC(=O)O"]), vec![amine]];
    let lib = EnumerateLibrary::new(&amide(), &blocks).unwrap();
    let bytes = lib.serialize().unwrap();
    let mut loaded = EnumerateLibrary::default();
    // SAFETY: the bytes came from serialize.
    unsafe { loaded.init_from_string(&bytes) }.unwrap();
    // The info's name (a u32 length, then its bytes), then its type.
    let at = bytes
        .windows(9)
        .position(|w| w == b"\x04\0\0\0MARK\x02")
        .unwrap()
        + 8;
    let mut changed = bytes.clone();
    changed[at] = 1;
    let mut loaded = EnumerateLibrary::default();
    // SAFETY: serialize's bytes with one byte changed inside a molecule
    // pickle, which RDKit reads with its own checks; Boost's framing is
    // serialize's.
    unsafe { loaded.init_from_string(&changed) }.unwrap();
    let reagents = loaded.get_reagents();
    let atom = reagents[1][0].get_atom_with_idx(0).unwrap();
    let info = atom.get_pdb_residue_info().unwrap().unwrap();
    assert_eq!(info.get_name().unwrap(), "MARK");
}
