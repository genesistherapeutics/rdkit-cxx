// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/SimDivFilters/UnitTestSimilarityPickers.py and the doctests
// of SimilarityPickers, with the molecules of RDKit's test database
// (simple_mols1) and values from Python's
// rdkit.SimDivFilters.SimilarityPickers in RDKit 2026.09.1.

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rdmolops::rdk_fingerprint;
use rdkit_cxx::data_structs::{ExplicitBitVect, dice_similarity};
use rdkit_cxx::sim_div_filters::similarity_pickers::*;

const DATA: [(&str, &str); 12] = [
    ("COC", "ether-1"),
    ("CC(=O)O", "acid-1"),
    ("c1ccccc1", "arom-1"),
    ("CCOC", "ether-2"),
    ("CCC(=O)O", "acid-2"),
    ("Cc1ccccc1", "arom-2"),
    ("CCOCC", "ether-3"),
    ("CC(C)C(=O)O", "acid-3"),
    ("Cc1cc(C)ccc1", "arom-3"),
    ("CCCOCC", "ether-4"),
    ("CCC(C)C(=O)O", "acid-4"),
    ("Cc1cc(C)cc(C)c1", "arom-4"),
];

fn fp(smiles: &str) -> ExplicitBitVect {
    rdk_fingerprint(&mol_from_smiles(smiles).unwrap()).unwrap()
}

fn data_set() -> Vec<ExplicitBitVect> {
    DATA.iter().map(|(smi, _)| fp(smi)).collect()
}

fn id(idx: usize) -> &'static str {
    DATA[idx].1
}

#[test]
fn top_n_overall_picker() {
    let fps = data_set();
    let probe = [fp("COC")];
    let mut picker = TopNOverallPicker::new_with_params(PickerParams {
        num_to_pick: 2,
        probe_fps: &probe,
        data_set: &fps,
        ..Default::default()
    });
    assert_eq!(picker.len().unwrap(), 2);
    let (idx, score) = picker.get(0).unwrap();
    assert_eq!((id(idx.unwrap()), score), ("ether-1", 1.0));
    assert_eq!(id(picker.get(1).unwrap().0.unwrap()), "ether-2");

    let probes = [fp("COC"), fp("CC(=O)O")];
    let expected = [
        ("acid-1", 1.0),
        ("ether-1", 1.0),
        ("acid-2", 0.6363636363636364),
        ("acid-3", 0.4666666666666667),
        ("ether-2", 0.4),
        ("ether-3", 0.3333333333333333),
        ("acid-4", 0.30434782608695654),
        ("ether-4", 0.21428571428571427),
        ("arom-1", 0.04),
        ("arom-2", 0.025),
        ("arom-3", 0.017857142857142856),
        ("arom-4", 0.017241379310344827),
    ];
    for n in [2, 3, 5, 12, 14] {
        let mut picker = TopNOverallPicker::new_with_params(PickerParams {
            num_to_pick: n,
            probe_fps: &probes,
            data_set: &fps,
            ..Default::default()
        });
        picker.make_picks().unwrap();
        assert_eq!(picker.len().unwrap(), n);
        for i in 0..n {
            let (idx, score) = picker.get(i as i64).unwrap();
            match idx {
                Some(idx) => assert_eq!((id(idx), score), expected[i]),
                // Python's empty slots, past the data set.
                None => assert_eq!(score, -1e99),
            }
        }
        assert!(matches!(picker.get(n as i64), Err(Error::Index { .. })));
    }
    picker.make_picks_with_force(true).unwrap();
}

#[test]
fn spread_picker() {
    let fps = data_set();
    let probe = [fp("COC")];
    let mut picker = SpreadPicker::new_with_params(PickerParams {
        num_to_pick: 2,
        probe_fps: &probe,
        data_set: &fps,
        ..Default::default()
    });
    assert_eq!(picker.len().unwrap(), 2);
    assert_eq!(picker.get(0).unwrap(), (0, 1.0));
    assert_eq!(id(picker.get(1).unwrap().0), "ether-2");

    let probes = [fp("COC"), fp("CC(=O)O")];
    let mut picker = SpreadPicker::new_with_params(PickerParams {
        num_to_pick: 3,
        probe_fps: &probes,
        data_set: &fps,
        ..Default::default()
    });
    let picks: Vec<(&str, f64)> = (0..3)
        .map(|i| {
            let (idx, score) = picker.get(i).unwrap();
            (id(idx), score)
        })
        .collect();
    assert_eq!(picks[..2], [("ether-1", 1.0), ("acid-1", 1.0)]);
    assert_eq!(picks[2].0, "ether-2");

    let expected = [
        ("ether-1", 1.0),
        ("acid-1", 1.0),
        ("ether-2", 0.5714285714285714),
        ("acid-2", 0.7777777777777778),
        ("ether-3", 0.5),
        ("acid-3", 0.6363636363636364),
        ("ether-4", 0.3333333333333333),
        ("acid-4", 0.4666666666666667),
        ("arom-2", 0.04878048780487805),
        ("arom-1", 0.07692307692307693),
        ("arom-3", 0.03508771929824561),
        ("arom-4", 0.028985507246376812),
    ];
    for n in [2, 3, 5, 12] {
        let mut picker = SpreadPicker::new_with_params(PickerParams {
            num_to_pick: n,
            probe_fps: &probes,
            data_set: &fps,
            sim_metric: dice_similarity,
        });
        let picks: Vec<(&str, f64)> = (0..n as i64)
            .map(|i| {
                let (idx, score) = picker.get(i).unwrap();
                (id(idx), score)
            })
            .collect();
        assert_eq!(picks, expected[..n]);
    }
    // More probes than picks: Python fails with an IndexError.
    let mut picker = SpreadPicker::new_with_params(PickerParams {
        num_to_pick: 1,
        probe_fps: &probes,
        data_set: &fps,
        ..Default::default()
    });
    assert!(matches!(picker.len(), Err(Error::Index { .. })));
    let mut no_probes = SpreadPicker::new_with_params(PickerParams {
        num_to_pick: 1,
        data_set: &fps,
        ..Default::default()
    });
    assert!(matches!(no_probes.make_picks(), Err(Error::Value { .. })));
}
