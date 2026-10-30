// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports the BuildFragmentCatalog tests of rdkit/Chem/UnitTestCatalog.py and
// tests rdkit.Chem.FragmentCatalog; values from Python's BuildFragmentCatalog
// and FragmentCatalog in RDKit 2024.09.6 (ProcessGainsFile's from a copy that
// reads lines as its Python 2 xreadlines meant to).

use std::collections::BTreeSet;

use rdkit_cxx::chem::build_fragment_catalog::*;
use rdkit_cxx::chem::fragment_catalog::*;
use rdkit_cxx::chem::{Mol, mol_from_smiles, mol_to_smiles};
use rdkit_cxx::{Error, Result};

const SMI_LIST2: [&str; 16] = [
    "OCCC", "CCC", "C=CC", "OC=CC", "CC(O)C", "C=C(O)C", "OCCCC", "CC(O)CC", "C=CCC", "CC=CC",
    "OC=CCC", "CC=C(O)C", "OCC=CC", "C=C(O)CC", "C=CC(O)C", "C=CCCO",
];
const LIST2_ACTS: [i64; 16] = [1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1];

fn mols() -> Vec<Mol> {
    SMI_LIST2
        .iter()
        .map(|s| mol_from_smiles(s).unwrap())
        .collect()
}

fn catalog(mols: &[Mol]) -> Result<FragCatalog> {
    build_catalog_with_params(
        mols,
        &BuildCatalogParams {
            min_path: 1,
            ..Default::default()
        },
    )
}

fn acts() -> Activities {
    Activities::Values(LIST2_ACTS.to_vec())
}

fn table(rows: &[[[i32; 2]; 2]]) -> ScoreTable {
    rows.iter().map(|[a, b]| [a.to_vec(), b.to_vec()]).collect()
}

fn obls() -> Vec<Vec<u32>> {
    [
        &[0, 1, 2][..],
        &[1, 3],
        &[1, 4, 5],
        &[1, 6, 7],
        &[0, 8],
        &[0, 6, 9],
        &[0, 1, 2, 3, 10],
        &[0, 1, 2, 8, 11],
        &[1, 3, 4, 5, 12],
        &[1, 4, 5, 13],
        &[1, 3, 6, 7, 14],
        &[0, 1, 6, 7, 9, 15],
        &[0, 1, 4, 5, 16, 17],
        &[0, 1, 2, 6, 9, 18],
        &[0, 4, 8, 16, 19],
        &[0, 1, 2, 4, 5, 20],
    ]
    .iter()
    .map(|b| b.to_vec())
    .collect()
}

/// Gains as `get_top_n` gives them, compared with Python's: in Python's order
/// of gains, but in any order among equal gains, which RDKit's ranker orders
/// by its C++ library's heap (libc++ and libstdc++ differ, as Python's own
/// test allows). When `top_n` cuts through the last tie, its bits need only
/// be among `last_tie`'s.
fn assert_gains(actual: &[Vec<f64>], expected: &[[f64; 4]], last_tie: Option<&[[f64; 4]]>) {
    assert_eq!(actual.len(), expected.len());
    let key = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };
    let mut start = 0;
    while start < expected.len() {
        let gain = expected[start][1];
        let end = start
            + expected[start..]
                .iter()
                .take_while(|e| e[1] == gain)
                .count();
        let mut want: Vec<Vec<u64>> = expected[start..end].iter().map(|e| key(e)).collect();
        let mut got: Vec<Vec<u64>> = actual[start..end]
            .iter()
            .map(|a| {
                assert!((a[1] - gain).abs() < 1e-12, "{a:?} != gain {gain}");
                let mut a = a.clone();
                a[1] = gain;
                key(&a)
            })
            .collect();
        want.sort();
        got.sort();
        match last_tie {
            Some(tie) if end == expected.len() => {
                let tie: Vec<Vec<u64>> = tie.iter().map(|e| key(e)).collect();
                assert!(
                    got.iter().all(|g| tie.contains(g)),
                    "{got:?} not in {tie:?}"
                );
            }
            _ => assert_eq!(got, want),
        }
        start = end;
    }
}

const TIE_0269: [[f64; 4]; 8] = [
    [18.0, 0.02692742888928845, 0.0, 1.0],
    [17.0, 0.02692742888928845, 0.0, 1.0],
    [20.0, 0.02692742888928845, 0.0, 1.0],
    [19.0, 0.02692742888928845, 0.0, 1.0],
    [15.0, 0.02692742888928845, 0.0, 1.0],
    [14.0, 0.02692742888928845, 0.0, 1.0],
    [11.0, 0.02692742888928845, 0.0, 1.0],
    [10.0, 0.02692742888928845, 0.0, 1.0],
];

/// Python's gains for `CalcGains(topN=10, biasList=(1,))`.
const BIASED_TOP10: [[f64; 4]; 10] = [
    [0.0, 0.4669171866886993, 0.0, 10.0],
    [6.0, 0.16113666504619273, 0.0, 5.0],
    [2.0, 0.16113666504619273, 0.0, 5.0],
    [7.0, 0.0877536667807961, 0.0, 3.0],
    [8.0, 0.0877536667807961, 0.0, 3.0],
    [9.0, 0.0877536667807961, 0.0, 3.0],
    [16.0, 0.056047626963330655, 0.0, 2.0],
    [15.0, 0.02692742888928845, 0.0, 1.0],
    [14.0, 0.02692742888928845, 0.0, 1.0],
    [11.0, 0.02692742888928845, 0.0, 1.0],
];

/// Python's gains for `CalcGains()`, all 21 bits.
const ALL_GAINS: [[f64; 4]; 21] = [
    [0.0, 0.4669171866886993, 0.0, 10.0],
    [5.0, 0.2057018851835084, 3.0, 2.0],
    [6.0, 0.16113666504619273, 0.0, 5.0],
    [2.0, 0.16113666504619273, 0.0, 5.0],
    [4.0, 0.1431558784658321, 3.0, 3.0],
    [12.0, 0.13447053550223065, 1.0, 0.0],
    [13.0, 0.13447053550223065, 1.0, 0.0],
    [7.0, 0.0877536667807961, 0.0, 3.0],
    [8.0, 0.0877536667807961, 0.0, 3.0],
    [9.0, 0.0877536667807961, 0.0, 3.0],
    [1.0, 0.0877536667807961, 4.0, 9.0],
    [3.0, 0.0737613082228673, 2.0, 2.0],
    [16.0, 0.056047626963330655, 0.0, 2.0],
    [18.0, 0.02692742888928845, 0.0, 1.0],
    [17.0, 0.02692742888928845, 0.0, 1.0],
    [20.0, 0.02692742888928845, 0.0, 1.0],
    [19.0, 0.02692742888928845, 0.0, 1.0],
    [15.0, 0.02692742888928845, 0.0, 1.0],
    [14.0, 0.02692742888928845, 0.0, 1.0],
    [11.0, 0.02692742888928845, 0.0, 1.0],
    [10.0, 0.02692742888928845, 0.0, 1.0],
];

// test6Builder (with test1CatGen's bits)
#[test]
fn builder() -> Result<()> {
    let mols = mols();
    let cat = catalog(&mols)?;
    assert_eq!(cat.get_num_entries(), 21);
    assert_eq!(cat.get_fp_length(), 21);
    let fpgen = FragFPGenerator;
    let obits = [3, 2, 3, 3, 2, 3, 5, 5, 5, 4, 5, 6];
    for (i, mol) in mols.iter().enumerate() {
        let bits = fpgen.get_fp_for_mol(mol, &cat)?.get_on_bits();
        if i < obits.len() {
            assert_eq!(bits.len(), obits[i]);
            assert_eq!(bits, obls()[i]);
        }
    }
    // max_pts stops after the first molecules.
    let cat = build_catalog_with_params(
        &mols,
        &BuildCatalogParams {
            min_path: 1,
            max_pts: 2,
            ..Default::default()
        },
    )?;
    assert_eq!(cat.get_fp_length(), 4);
    Ok(())
}

// test7ScoreMolecules
#[test]
fn score_molecules_() -> Result<()> {
    let mols = mols();
    let cat = catalog(&mols)?;
    let params = ScoreMoleculesParams {
        activities: acts(),
        ..Default::default()
    };
    let opt: Vec<Option<Mol>> = mols.iter().cloned().map(Some).collect();
    let (scores, bit_lists) = score_molecules_with_params(&opt, &cat, &params)?;
    assert_eq!(bit_lists, obls());
    #[rustfmt::skip]
    let expected = table(&[
        [[0, 3], [4, 9]], [[4, 7], [0, 5]], [[2, 10], [2, 2]], [[1, 9], [3, 3]],
        [[1, 10], [3, 2]], [[4, 7], [0, 5]], [[4, 9], [0, 3]], [[4, 9], [0, 3]],
        [[4, 9], [0, 3]], [[4, 11], [0, 1]], [[4, 11], [0, 1]], [[3, 12], [1, 0]],
        [[3, 12], [1, 0]], [[4, 11], [0, 1]], [[4, 11], [0, 1]], [[4, 10], [0, 2]],
        [[4, 11], [0, 1]], [[4, 11], [0, 1]], [[4, 11], [0, 1]], [[4, 11], [0, 1]],
        [[4, 2], [0, 10]],
    ]);
    assert_eq!(scores, expected);

    let scores2 = score_from_lists_with_params(
        &bit_lists,
        &mols,
        &cat,
        &ScoreFromListsParams {
            activities: acts(),
            ..Default::default()
        },
    )?;
    assert_eq!(scores2, scores);

    // Activities read from a property score the same.
    let mut named = mols.clone();
    for (m, a) in named.iter_mut().zip(LIST2_ACTS) {
        m.set_prop("act", &a.to_string());
    }
    let by_prop = ScoreFromListsParams {
        activities: Activities::Prop("act".into()),
        ..Default::default()
    };
    assert_eq!(
        score_from_lists_with_params(&bit_lists, &named, &cat, &by_prop)?,
        scores
    );
    let last = ScoreFromListsParams::default();
    assert_eq!(
        score_from_lists_with_params(&bit_lists, &named, &cat, &last)?,
        scores
    );
    // Without properties, there is no activity to read.
    assert!(matches!(
        score_from_lists(&bit_lists, &mols, &cat),
        Err(Error::Value { .. })
    ));

    // max_pts limits the molecules scored.
    let scores5 = score_from_lists_with_params(
        &bit_lists,
        &mols,
        &cat,
        &ScoreFromListsParams {
            max_pts: 5,
            activities: acts(),
            ..Default::default()
        },
    )?;
    #[rustfmt::skip]
    let expected5 = table(&[
        [[0, 1], [2, 2]], [[2, 2], [0, 1]], [[1, 3], [1, 0]], [[1, 3], [1, 0]],
        [[1, 3], [1, 0]], [[2, 2], [0, 1]], [[2, 2], [0, 1]], [[2, 2], [0, 1]],
        [[2, 3], [0, 0]], [[2, 3], [0, 0]], [[2, 3], [0, 0]], [[2, 3], [0, 0]],
        [[2, 3], [0, 0]], [[2, 3], [0, 0]], [[2, 3], [0, 0]], [[2, 3], [0, 0]],
        [[2, 3], [0, 0]], [[2, 3], [0, 0]], [[2, 3], [0, 0]], [[2, 3], [0, 0]],
        [[2, 1], [0, 2]],
    ]);
    assert_eq!(scores5, expected5);

    // A molecule that didn't parse has no bits and isn't counted.
    let mut opt = opt;
    opt[1] = None;
    let (scores_none, bit_lists_none) = score_molecules_with_params(&opt, &cat, &params)?;
    assert!(bit_lists_none[1].is_empty());
    assert_eq!(scores_none[0][0], vec![0, 3]);
    assert_eq!(scores_none[1][0], vec![3, 7]);
    Ok(())
}

// test8MolRanks
#[test]
fn mol_ranks() -> Result<()> {
    let mols = mols();
    let cat = catalog(&mols)?;
    let (gains, fps) = calc_gains_with_params(
        &mols,
        &cat,
        &CalcGainsParams {
            top_n: 10,
            activities: acts(),
            bias_list: Some(vec![1]),
            collect_fps: true,
            ..Default::default()
        },
    )?;
    assert_gains(&gains, &BIASED_TOP10, Some(&TIE_0269));
    assert_eq!(cat.get_bit_description(0)?, "C<-O>C");
    assert_eq!(cat.get_bit_description(16)?, "C=CC<-O>");
    assert_eq!(fps.len(), 16);

    let (gains, no_fps) = calc_gains_with_params(
        &mols,
        &cat,
        &CalcGainsParams {
            activities: acts(),
            ..Default::default()
        },
    )?;
    assert!(no_fps.is_empty());
    assert_eq!(gains.len(), 21);
    assert_gains(&gains, &ALL_GAINS, None);
    assert_eq!(cat.get_bit_description(5)?, "C=CC");

    let from_fps = calc_gains_from_fps_with_params(
        &mols,
        &fps,
        &CalcGainsFromFpsParams {
            top_n: 5,
            activities: acts(),
            ..Default::default()
        },
    )?;
    assert_eq!(from_fps, gains[..5]);
    assert!(matches!(
        calc_gains_from_fps(&mols, &[]),
        Err(Error::Index { .. })
    ));
    Ok(())
}

// test9Issue116
#[test]
fn issue116() -> Result<()> {
    let cat = build_catalog_with_params(
        &[mol_from_smiles("Cc1ccccc1")?],
        &BuildCatalogParams {
            min_path: 2,
            max_path: 2,
            ..Default::default()
        },
    )?;
    assert_eq!(cat.get_fp_length(), 2);
    assert_eq!(cat.get_bit_description(0)?, "ccC");
    let fpgen = FragFPGenerator;
    let fp = fpgen.get_fp_for_mol(&mol_from_smiles("Cc1ccccc1")?, &cat)?;
    assert!(fp.get_bit(0)?);
    assert!(fp.get_bit(1)?);
    let fp = fpgen.get_fp_for_mol(&mol_from_smiles("c1ccccc1-c1ccccc1")?, &cat)?;
    assert!(!fp.get_bit(0)?);
    assert!(fp.get_bit(1)?);
    Ok(())
}

const GAINS_CSV: &str = "id,Description,Gain,Act-0,Act-1\n0,C<-O>C,0.466917,0,10\n\
6,C<-O>=C,0.161137,0,5\n2,C<-O>CC,0.161137,0,5\n7,C<-O>=CC,0.087754,0,3\n\
8,CC<-O>C,0.087754,0,3\n9,C=C<-O>C,0.087754,0,3\n16,C=CC<-O>,0.056048,0,2\n\
15,CC=C<-O>C,0.026927,0,1\n14,C<-O>=CCC,0.026927,0,1\n11,CC<-O>CC,0.026927,0,1\n";

#[test]
fn gains_data() -> Result<()> {
    // Python's gains, in Python's order, so the text is Python's whatever
    // order this platform's ranker breaks ties in.
    let cat = catalog(&mols())?;
    let to_vecs = |g: &[[f64; 4]]| -> Vec<Vec<f64>> { g.iter().map(|e| e.to_vec()).collect() };
    let mut out = Vec::new();
    output_gains_data(&mut out, &to_vecs(&BIASED_TOP10), Some(&cat))?;
    assert_eq!(String::from_utf8(out).unwrap(), GAINS_CSV);

    let gains5 = to_vecs(&ALL_GAINS[..5]);
    let mut out = Vec::new();
    output_gains_data_with_n_acts(&mut out, &gains5, None, 2)?;
    let no_cat = String::from_utf8(out).unwrap();
    assert_eq!(
        no_cat,
        "id,Gain,Act-0,Act-1\n0,0.466917,0,10\n5,0.205702,3,2\n6,0.161137,0,5\n\
         2,0.161137,0,5\n4,0.143156,3,3\n"
    );

    let read = process_gains_data(no_cat.as_bytes())?;
    let expected: Vec<(String, f64)> = [
        ("0", 0.466917),
        ("5", 0.205702),
        ("6", 0.161137),
        ("2", 0.161137),
        ("4", 0.143156),
    ]
    .iter()
    .map(|(i, g)| (i.to_string(), *g))
    .collect();
    assert_eq!(read, expected);
    let read = process_gains_data_with_params(
        GAINS_CSV.as_bytes(),
        &ProcessGainsDataParams {
            gain_col: 2,
            ..Default::default()
        },
    )?;
    assert_eq!(read.len(), 10);
    assert_eq!(read[6], ("16".to_owned(), 0.056048));
    // The description isn't a gain.
    assert!(matches!(
        process_gains_data(GAINS_CSV.as_bytes()),
        Err(Error::Value { .. })
    ));
    Ok(())
}

fn gains_file() -> std::path::PathBuf {
    let path =
        std::env::temp_dir().join(format!("rdkit-cxx-test-gains-{}.csv", std::process::id()));
    std::fs::write(&path, GAINS_CSV).unwrap();
    path
}

#[test]
fn fragment_catalog() -> Result<()> {
    let path = gains_file();
    let bits = process_gains_file(&path)?;
    let first = process_gains_file_with_params(
        &path,
        &ProcessGainsFileParams {
            n_to_do: Some(3),
            ..Default::default()
        },
    )?;
    let no_descr = process_gains_file_with_params(
        &path,
        &ProcessGainsFileParams {
            have_descriptions: false,
            ..Default::default()
        },
    );
    std::fs::remove_file(&path).unwrap();
    assert_eq!(bits.len(), 10);
    assert_eq!(
        bits[0],
        BitGainsInfo {
            id: 0,
            description: "C<-O>C".into(),
            gain: 0.466917,
            n_per_class: Some(vec![0, 10]),
        }
    );
    assert_eq!(
        bits[6],
        BitGainsInfo {
            id: 16,
            description: "C=CC<-O>".into(),
            gain: 0.056048,
            n_per_class: Some(vec![0, 2]),
        }
    );
    assert_eq!(first, bits[..3]);
    assert!(matches!(no_descr, Err(Error::Value { .. })));
    assert_eq!(BitGainsInfo::default().id, -1);

    let mols = mols();
    let cat = catalog(&mols)?;
    let set = |ids: &[u32]| ids.iter().copied().collect::<BTreeSet<u32>>();
    let adj = build_adjacency_list(&cat, &bits)?;
    assert_eq!(
        adj.adjs,
        vec![
            (0, vec![2, 8, 9, 16]),
            (6, vec![7, 9]),
            (2, vec![11]),
            (7, vec![14, 15]),
            (8, vec![11]),
            (9, vec![15]),
            (16, vec![]),
            (15, vec![]),
            (14, vec![]),
            (11, vec![]),
        ]
    );
    assert_eq!(
        adj.levels,
        vec![
            (1, set(&[0, 6])),
            (2, set(&[2, 7, 8, 9, 16])),
            (3, set(&[11, 14, 15])),
        ]
    );
    let adj = build_adjacency_list_with_limit_inclusion(&cat, &bits, false)?;
    assert_eq!(adj.adjs[2], (2, vec![10, 11, 18, 20]));
    assert_eq!(adj.adjs[6], (16, vec![17, 19]));
    assert_eq!(adj.levels[2], (3, set(&[10, 11, 14, 15, 17, 18, 19, 20])));

    let (_, fps) = calc_gains_with_params(
        &mols,
        &cat,
        &CalcGainsParams {
            activities: acts(),
            collect_fps: true,
            ..Default::default()
        },
    )?;
    let smiles =
        |ms: Vec<&Mol>| -> Vec<String> { ms.iter().map(|m| mol_to_smiles(m).unwrap()).collect() };
    assert_eq!(
        smiles(get_mols_matching_bit(&mols, bits[1].id as u32, &fps)?),
        ["CC=CO", "C=C(C)O", "CCC=CO", "CC=C(C)O", "C=C(O)CC"]
    );
    assert_eq!(
        smiles(get_mols_matching_bit(&mols, 16, &fps)?),
        ["CC=CCO", "C=CC(C)O"]
    );
    assert!(get_mols_matching_bit(&mols, 16, &fps[..2]).is_err());
    Ok(())
}
