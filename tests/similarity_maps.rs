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
use rdkit_cxx::data_structs::{ExplicitBitVect, tanimoto_similarity};

fn sha1(text: &str) -> String {
    let mut h = sha1_smol::Sha1::new();
    h.update(text.as_bytes());
    h.digest().to_string()
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
    assert_eq!(
        sha1(&d.get_drawing_text().unwrap()),
        "848e90f20af14f4d3afc8f1aa04632bdd2908d84"
    );
    assert!(!d.draw_options().clear_background);

    let fpg = morgan2();
    let mut d = MolDraw2DSVG::new(400, 400).unwrap();
    let max_weight =
        get_similarity_map_for_fingerprint_generator(&refmol(), &probe(), &fpg, &mut d).unwrap();
    d.finish_drawing();
    assert_eq!(max_weight, 0.12255947497949138);
    assert_eq!(
        sha1(&d.get_drawing_text().unwrap()),
        "09d0631bc598e8d48209261d0966d20c4bce3ff8"
    );

    // Github #2904: a colour map of three colours.
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
    assert_eq!(
        sha1(&d.get_drawing_text().unwrap()),
        "5c925d0355e36038b2efe56a7e80ef72407981bc"
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
    assert_eq!(
        sha1(&d.get_drawing_text().unwrap()),
        "4b38cf3a5fb9cd6f5b5ef2f67e1666dadce9c8b5"
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
    assert_eq!(
        sha1(&d.get_drawing_text().unwrap()),
        "19e4fe241fc791873a000592fbcbf9c170879d83"
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
