// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The ported rdkit.Chem.Draw.SimilarityMaps against Python's in RDKit
// 2026.09.1, whose SVGs it reproduces byte for byte (by SHA-1), after
// UnitTestSimilarityMaps.py. Python's tests of the atom-pair, torsion, and
// Morgan fingerprint functions, which are not ported, are ported with the
// RDKit fingerprint and a Morgan generator.

use rdkit_cxx::Error;
use rdkit_cxx::chem::draw::rd_mol_draw_2d::{DrawColour, MolDraw2DSVG};
use rdkit_cxx::chem::draw::similarity_maps::*;
use rdkit_cxx::chem::rd_fingerprint_generator::{
    FingerprintGenerator64, GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::chem::{Mol, mol_from_smiles};
use rdkit_cxx::data_structs::{
    ExplicitBitVect, dice_similarity, dice_similarity_sparse_int_vect, tanimoto_similarity,
};

fn sha1(text: &str) -> String {
    let mut h = sha1_smol::Sha1::new();
    h.update(text.as_bytes());
    h.digest().to_string()
}

/// Checks a similarity map against Python's. Its contours are floating-point
/// output (exponentials and the contour lines through them), whose last
/// digits differ across platforms, so the SVG is pinned (by SHA-1) only on
/// Linux x86-64, where Python's was made; elsewhere, the map must draw each
/// of the molecule's bonds over its contours.
fn assert_map(svg: &str, python_sha1: &str, mol: &Mol) {
    if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
        assert_eq!(sha1(svg), python_sha1);
    }
    assert!(svg.trim_end().ends_with("</svg>"));
    let first_bond = svg
        .find("class='bond-0 ")
        .expect("the map draws the molecule");
    assert!(
        svg[..first_bond].contains("<path"),
        "the contours come first"
    );
    for bond in 0..mol.get_num_bonds() {
        assert!(svg.contains(&format!("class='bond-{bond} ")));
    }
}

fn benzene() -> Mol {
    mol_from_smiles("c1ccccc1").unwrap()
}

fn pyridine() -> Mol {
    mol_from_smiles("c1ccncc1").unwrap()
}

fn rdk(m: &Mol, atom_id: Option<u32>) -> rdkit_cxx::Result<ExplicitBitVect> {
    get_rdk_fingerprint_with_params(
        m,
        &GetRDKFingerprintParams {
            atom_id,
            ..Default::default()
        },
    )
}

fn morgan2() -> FingerprintGenerator64 {
    get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 2,
        ..Default::default()
    })
    .unwrap()
}

fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() < 1e-12, "{actual:?} != {expected:?}");
    }
}

// UnitTestSimilarityMaps.testSimilarityMap, the RDKit fingerprint.
#[test]
fn atomic_weights_for_fingerprint_match_python() {
    let weights = get_atomic_weights_for_fingerprint(&benzene(), &pyridine(), |m, i| {
        get_rdk_fingerprint_with_params(
            m,
            &GetRDKFingerprintParams {
                atom_id: i,
                n_bits: 1024,
                n_bits_per_hash: 1,
                ..Default::default()
            },
        )
    })
    .unwrap();
    let (a, b) = (0.42105263157894735, -0.32894736842105265);
    assert_close(&weights, &[a, a, a, b, a, a]);
    let weights = get_atomic_weights_for_fingerprint(&benzene(), &pyridine(), rdk).unwrap();
    let (a, b) = (0.41025641025641024, -0.4132730015082956);
    assert_close(&weights, &[a, a, a, b, a, a]);
    assert!(matches!(rdk(&benzene(), Some(6)), Err(Error::Value { .. })));
    assert_eq!(
        rdk(&benzene(), None).unwrap(),
        get_rdk_fingerprint(&benzene()).unwrap()
    );
}

#[test]
fn atomic_weights_for_fingerprint_generator_match_python() {
    let fpg = morgan2();
    let weights =
        get_atomic_weights_for_fingerprint_generator(&benzene(), &pyridine(), &fpg).unwrap();
    assert_close(&weights, &[0.5, 0.5, 0.5, -0.5, 0.5, 0.5]);
    let weights =
        get_atomic_weights_for_fingerprint_generator_use_counts(&benzene(), &pyridine(), &fpg)
            .unwrap();
    assert_close(
        &weights,
        &[
            0.20370370370370372,
            0.12962962962962965,
            0.05555555555555558,
            -0.16666666666666663,
            0.05555555555555558,
            0.12962962962962965,
        ],
    );
    let weights = get_atomic_weights_for_fingerprint_generator_with_metric(
        &benzene(),
        &pyridine(),
        &fpg,
        tanimoto_similarity,
    )
    .unwrap();
    let (a, b) = (1.0 / 3.0, -0.6666666666666667);
    assert_close(&weights, &[a, a, a, b, a, a]);
}

#[test]
fn standardized_weights_match_python() {
    // Morgan weights of 0.5 standardize to 1 (testSimilarityMap).
    let (weights, max) = get_standardized_weights(&[0.5, 0.5, 0.5, -0.5, 0.5, 0.5]).unwrap();
    assert_eq!(weights, [1.0, 1.0, 1.0, -1.0, 1.0, 1.0]);
    assert_eq!(max, 0.5);
    assert_eq!(
        get_standardized_weights(&[0.5, -1.0, 0.25]).unwrap(),
        (vec![0.5, -1.0, 0.25], 1.0)
    );
    assert_eq!(
        get_standardized_weights(&[0.0, 0.0]).unwrap(),
        (vec![0.0, 0.0], 0.0)
    );
    assert!(get_standardized_weights(&[]).is_err());
}

fn refmol() -> Mol {
    mol_from_smiles("CCCN(CCCCN1CCN(c2ccccc2OC)CC1)Cc1ccc2ccccc2c1").unwrap()
}

fn probe() -> Mol {
    mol_from_smiles("COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21").unwrap()
}

// UnitTestSimilarityMaps.testSimilarityMapsMolDraw2D.
#[test]
fn similarity_maps_match_python() {
    let mut d = MolDraw2DSVG::new(400, 400).unwrap();
    d.clear_drawing().unwrap();
    let max_weight = get_similarity_map_for_fingerprint(&refmol(), &probe(), rdk, &mut d).unwrap();
    d.finish_drawing();
    assert_eq!(max_weight, 0.18541345659989728);
    assert_map(
        &d.get_drawing_text().unwrap(),
        "848e90f20af14f4d3afc8f1aa04632bdd2908d84",
        &probe(),
    );
    assert!(!d.draw_options().clear_background);

    let fpg = morgan2();
    let mut d = MolDraw2DSVG::new(400, 400).unwrap();
    let max_weight =
        get_similarity_map_for_fingerprint_generator(&refmol(), &probe(), &fpg, &mut d).unwrap();
    d.finish_drawing();
    assert_eq!(max_weight, 0.12255947497949138);
    assert_map(
        &d.get_drawing_text().unwrap(),
        "09d0631bc598e8d48209261d0966d20c4bce3ff8",
        &probe(),
    );

    // Github #2904: a color map of three colors.
    let mut d = MolDraw2DSVG::new(400, 400).unwrap();
    let params = GetSimilarityMapForFingerprintGeneratorParams {
        map: GetSimilarityMapFromWeightsParams {
            color_map: Some([
                DrawColour {
                    r: 0.0,
                    g: 1.0,
                    b: 0.0,
                    a: 0.5,
                },
                DrawColour::new(1.0, 1.0, 1.0),
                DrawColour {
                    r: 0.0,
                    g: 0.0,
                    b: 1.0,
                    a: 0.5,
                },
            ]),
            contour_lines: 5,
            sigma: Some(0.4),
            ..Default::default()
        },
        ..Default::default()
    };
    let max_weight = get_similarity_map_for_fingerprint_generator_use_counts_with_params(
        &refmol(),
        &probe(),
        &fpg,
        &mut d,
        &params,
    )
    .unwrap();
    d.finish_drawing();
    assert_eq!(max_weight, 0.08955961037179305);
    assert_map(
        &d.get_drawing_text().unwrap(),
        "5c925d0355e36038b2efe56a7e80ef72407981bc",
        &probe(),
    );
}

#[test]
fn similarity_map_from_weights_matches_python() {
    let mut d = MolDraw2DSVG::new(300, 300).unwrap();
    let params = GetSimilarityMapFromWeightsParams {
        use_fill_threshold: true,
        grid_resolution: 0.05,
        extra_grid_padding: 1.0,
        ..Default::default()
    };
    get_similarity_map_from_weights_with_params(
        &pyridine(),
        &[0.1, -0.2, 0.3, -0.4, 0.5, -0.6],
        &mut d,
        &params,
    )
    .unwrap();
    d.finish_drawing();
    assert_map(
        &d.get_drawing_text().unwrap(),
        "4b38cf3a5fb9cd6f5b5ef2f67e1666dadce9c8b5",
        &pyridine(),
    );
    let mut d = MolDraw2DSVG::new(300, 300).unwrap();
    let methane = mol_from_smiles("C").unwrap();
    let Err(Error::Value { message }) = get_similarity_map_from_weights(&methane, &[1.0], &mut d)
    else {
        panic!("one atom maps");
    };
    assert_eq!(message, "too few atoms");
}

#[test]
fn similarity_map_for_model_matches_python() {
    let predict = |fp: &ExplicitBitVect| Ok(f64::from(fp.get_num_on_bits()) / 100.0);
    let weights = get_atomic_weights_for_model(&pyridine(), rdk, predict).unwrap();
    assert_close(
        &weights,
        &[
            0.19999999999999998,
            0.26,
            0.29,
            0.21999999999999997,
            0.29,
            0.26,
        ],
    );
    let mut d = MolDraw2DSVG::new(300, 300).unwrap();
    let max_weight = get_similarity_map_for_model(&pyridine(), rdk, predict, &mut d).unwrap();
    d.finish_drawing();
    assert_eq!(max_weight, 0.29);
    assert_map(
        &d.get_drawing_text().unwrap(),
        "19e4fe241fc791873a000592fbcbf9c170879d83",
        &pyridine(),
    );
}

#[test]
fn similarity_map_params_defaults_match_python() {
    let map = GetSimilarityMapFromWeightsParams::default();
    assert_eq!(map.color_map, None);
    assert_eq!(map.sigma, None);
    assert_eq!(map.contour_lines, 10);
    assert_eq!(map.grid_resolution, 0.1);
    assert_eq!(map.extra_grid_padding, 0.5);
    assert!(!map.use_fill_threshold);
    assert_eq!(map.fill_threshold, 0.01);
    assert!(map.fill_threshold_is_fraction);
    let rdk = GetRDKFingerprintParams::default();
    assert_eq!(
        (
            rdk.atom_id,
            rdk.n_bits,
            rdk.min_path,
            rdk.max_path,
            rdk.n_bits_per_hash
        ),
        (None, 2048, 1, 5, 2)
    );
}

fn weights<F>(
    a: &str,
    b: &str,
    fp: impl FnMut(&Mol, Option<u32>) -> rdkit_cxx::Result<F>,
    metric: impl Fn(&F, &F) -> rdkit_cxx::Result<f64>,
) -> Vec<f64> {
    get_atomic_weights_for_fingerprint_with_metric(
        &mol_from_smiles(a).unwrap(),
        &mol_from_smiles(b).unwrap(),
        fp,
        metric,
    )
    .unwrap()
}

fn bz_py<F>(
    fp: impl FnMut(&Mol, Option<u32>) -> rdkit_cxx::Result<F>,
    metric: impl Fn(&F, &F) -> rdkit_cxx::Result<f64>,
) -> Vec<f64> {
    weights("c1ccccc1", "c1ccncc1", fp, metric)
}

// UnitTestSimilarityMaps.testSimilarityMap, against Python's values.
#[test]
fn deprecated_fingerprint_weights_match_python() {
    let morgan = |m: &Mol, i| {
        get_morgan_fingerprint_with_params(
            m,
            &GetMorganFingerprintParams {
                atom_id: i,
                ..Default::default()
            },
        )
    };
    assert_close(
        &bz_py(morgan, dice_similarity),
        &[0.5, 0.5, 0.5, -0.5, 0.5, 0.5],
    );
    let morgan_count = |m: &Mol, i| {
        get_morgan_fingerprint_count_with_params(
            m,
            &GetMorganFingerprintParams {
                atom_id: i,
                ..Default::default()
            },
        )
    };
    assert_close(
        &bz_py(morgan_count, dice_similarity_sparse_int_vect),
        &[
            0.20370370370370372,
            0.12962962962962965,
            0.05555555555555558,
            -0.16666666666666663,
            0.05555555555555558,
            0.12962962962962965,
        ],
    );
    let morgan_features = |m: &Mol, i| {
        get_morgan_fingerprint_with_params(
            m,
            &GetMorganFingerprintParams {
                atom_id: i,
                use_features: true,
                ..Default::default()
            },
        )
    };
    assert_close(
        &bz_py(morgan_features, dice_similarity),
        &[0.5, 0.5, 0.5, -0.5, 0.5, 0.5],
    );

    let ap_bv = |m: &Mol, i| {
        get_ap_fingerprint_bv_with_params(
            m,
            &GetAPFingerprintParams {
                atom_id: i,
                n_bits: 1024,
                ..Default::default()
            },
        )
    };
    assert_close(
        &bz_py(ap_bv, dice_similarity),
        &[
            0.09523809523809523,
            0.17366946778711478,
            0.17366946778711478,
            -0.23809523809523814,
            0.17366946778711478,
            0.17366946778711478,
        ],
    );
    let ap_params = |i| GetAPFingerprintParams {
        atom_id: i,
        ..Default::default()
    };
    let (a, b) = (0.18666666666666665, -0.13333333333333341);
    assert_close(
        &bz_py(
            |m: &Mol, i| get_ap_fingerprint_with_params(m, &ap_params(i)),
            dice_similarity_sparse_int_vect,
        ),
        &[a, a, a, b, a, a],
    );
    assert_close(
        &bz_py(
            |m: &Mol, i| get_ap_fingerprint_hashed_with_params(m, &ap_params(i)),
            dice_similarity_sparse_int_vect,
        ),
        &[a, a, a, b, a, a],
    );

    let tt_bv = |m: &Mol, i| {
        get_tt_fingerprint_bv_with_params(
            m,
            &GetTTFingerprintParams {
                atom_id: i,
                n_bits: 1024,
                n_bits_per_entry: 1,
                ..Default::default()
            },
        )
    };
    let c = -0.16666666666666663;
    assert_close(&bz_py(tt_bv, dice_similarity), &[0.5, 0.5, c, -0.5, c, 0.5]);
    let tt_params = |i| GetTTFingerprintParams {
        atom_id: i,
        ..Default::default()
    };
    let expected = [
        1.0 / 3.0,
        1.0 / 3.0,
        0.08333333333333331,
        -0.16666666666666669,
        0.08333333333333331,
        1.0 / 3.0,
    ];
    assert_close(
        &bz_py(
            |m: &Mol, i| get_tt_fingerprint_with_params(m, &tt_params(i)),
            dice_similarity_sparse_int_vect,
        ),
        &expected,
    );
    assert_close(
        &bz_py(
            |m: &Mol, i| get_tt_fingerprint_hashed_with_params(m, &tt_params(i)),
            dice_similarity_sparse_int_vect,
        ),
        &expected,
    );
    for err in [
        get_ap_fingerprint_with_params(&benzene(), &ap_params(Some(6))).err(),
        get_tt_fingerprint_with_params(&benzene(), &tt_params(Some(6))).err(),
        get_morgan_fingerprint_with_params(
            &benzene(),
            &GetMorganFingerprintParams {
                atom_id: Some(6),
                ..Default::default()
            },
        )
        .err(),
    ] {
        assert!(matches!(err, Some(Error::Value { .. })));
    }
}

// UnitTestSimilarityMaps.testSimilarityMapKWArgs, against Python's values.
#[test]
fn deprecated_fingerprint_keywords_match_python() {
    let (m1, m2) = ("CC[C@](F)(Cl)c1ccccc1", "CC[C@@](F)(Cl)c1ccccc1");
    let ap = |include_chirality| {
        move |m: &Mol, i| {
            get_ap_fingerprint_with_params(
                m,
                &GetAPFingerprintParams {
                    atom_id: i,
                    include_chirality,
                    ..Default::default()
                },
            )
        }
    };
    let w = 0.09999999999999998;
    assert_close(
        &weights(m1, m2, ap(false), dice_similarity_sparse_int_vect),
        &[w; 11],
    );
    let (w, c) = (0.09818181818181826, -0.08181818181818179);
    assert_close(
        &weights(m1, m2, ap(true), dice_similarity_sparse_int_vect),
        &[w, w, c, w, w, w, w, w, w, w, w],
    );
    let tt = |m: &Mol, i| {
        get_tt_fingerprint_with_params(
            m,
            &GetTTFingerprintParams {
                atom_id: i,
                ..Default::default()
            },
        )
    };
    assert_close(
        &weights(m1, m2, tt, dice_similarity_sparse_int_vect),
        &[
            0.09677419354838712,
            0.1724137931034483,
            0.4782608695652174,
            0.09677419354838712,
            0.09677419354838712,
            0.6190476190476191,
            0.3076923076923077,
            0.1724137931034483,
            0.1333333333333333,
            0.1724137931034483,
            0.3076923076923077,
        ],
    );
    let tt_chiral = |m: &Mol, i| {
        get_tt_fingerprint_hashed_with_params(
            m,
            &GetTTFingerprintParams {
                atom_id: i,
                include_chirality: true,
                ..Default::default()
            },
        )
    };
    assert_close(
        &weights(m1, m2, tt_chiral, dice_similarity_sparse_int_vect),
        &[
            -0.03415559772296012,
            -0.0608519269776876,
            -0.16879795396419434,
            -0.03415559772296012,
            -0.03415559772296012,
            0.1624649859943978,
            0.1990950226244344,
            0.21501014198782964,
            0.21960784313725493,
            0.21501014198782964,
            0.1990950226244344,
        ],
    );
    let morgan = |use_chirality| {
        move |m: &Mol, i| {
            get_morgan_fingerprint_with_params(
                m,
                &GetMorganFingerprintParams {
                    atom_id: i,
                    radius: 1,
                    use_chirality,
                    ..Default::default()
                },
            )
        }
    };
    let weights1 = weights(m1, m2, morgan(false), dice_similarity);
    assert_close(
        &weights1,
        &[
            0.11111111111111116,
            0.15384615384615385,
            0.25,
            0.11111111111111116,
            0.11111111111111116,
            0.15384615384615385,
            0.15384615384615385,
            0.11111111111111116,
            0.0714285714285714,
            0.11111111111111116,
            0.15384615384615385,
        ],
    );
    let weights2 = weights(m1, m2, morgan(true), dice_similarity);
    assert_close(
        &weights2,
        &[
            0.11851851851851858,
            0.0871794871794872,
            0.18333333333333335,
            0.04444444444444451,
            0.04444444444444451,
            0.0871794871794872,
            0.16410256410256407,
            0.11851851851851858,
            0.07619047619047625,
            0.11851851851851858,
            0.16410256410256407,
        ],
    );
    // The chiral center's contribution drops.
    assert!(weights1[2] > weights2[2]);
}

// The common recipe: a map of Morgan bit-vector weights.
#[test]
fn morgan_similarity_map_matches_python() {
    let mut d = MolDraw2DSVG::new(400, 400).unwrap();
    d.clear_drawing().unwrap();
    let morgan = |m: &Mol, i| {
        get_morgan_fingerprint_with_params(
            m,
            &GetMorganFingerprintParams {
                atom_id: i,
                ..Default::default()
            },
        )
    };
    let max_weight =
        get_similarity_map_for_fingerprint(&refmol(), &probe(), morgan, &mut d).unwrap();
    d.finish_drawing();
    assert_eq!(max_weight, 0.12255947497949138);
    assert_map(
        &d.get_drawing_text().unwrap(),
        "3093893d6b352ae3a9e64085ed92acb384146395",
        &probe(),
    );
}

/// Python's bindings take colors as tuples of C floats, so a color map's
/// components are rounded to `f32`, which shows in the fill colors.
#[test]
fn color_maps_are_rounded_as_python_rounds_them() {
    let c = DrawColour::new(0.1, 0.2, 0.3);
    assert_eq!(
        (c.r, c.g, c.b, c.a),
        (
            0.10000000149011612,
            0.20000000298023224,
            0.30000001192092896,
            1.0
        )
    );
    let mut d = MolDraw2DSVG::new(250, 250).unwrap();
    let params = GetSimilarityMapFromWeightsParams {
        color_map: Some([
            DrawColour::new(0.1, 0.2, 0.9),
            DrawColour::new(1.0, 1.0, 1.0),
            DrawColour::new(0.9, 0.1, 0.1),
        ]),
        ..Default::default()
    };
    get_similarity_map_from_weights_with_params(
        &mol_from_smiles("CC").unwrap(),
        &[-1.0, 0.0],
        &mut d,
        &params,
    )
    .unwrap();
    d.finish_drawing();
    let svg = d.get_drawing_text().unwrap();
    // Python's fill; unrounded components give #1932E5.
    assert!(svg.contains("#1933E5"));
    assert!(!svg.contains("#1932E5"));
}
