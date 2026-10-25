// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mhfp_fingerprint::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdMHFPFingerprint.

fn encoder() -> MHFPEncoder {
    MHFPEncoder::new_with_params(&MHFPEncoderParams {
        n_permutations: 128,
        seed: 42,
    })
}

fn secfp64() -> EncodeSECFPParams {
    EncodeSECFPParams {
        length: 64,
        ..Default::default()
    }
}

#[test]
fn min_hashes() {
    let enc = encoder();
    let mh = enc.encode_smiles("CCO").unwrap();
    assert_eq!(mh.len(), 128);
    assert_eq!(
        mh[..5],
        [998603975, 113829445, 99459723, 323656167, 2238542606]
    );
    assert_eq!(
        enc.encode_mol(&mol_from_smiles("CCO").unwrap()).unwrap(),
        mh
    );
    assert_eq!(
        enc.from_string_array(&["a", "b"])[..3],
        [2212273767, 865464383, 3220240811]
    );
    assert_eq!(
        enc.from_array(&[1, 2, 3])[..3],
        [734796314, 152805157, 1046592797]
    );
    let other = enc.encode_smiles("CCN").unwrap();
    assert_eq!(MHFPEncoder::distance(&mh, &other).unwrap(), 0.8984375);
    let bulk = enc.encode_smiles_bulk(&["CCO", "CCN"]).unwrap();
    assert_eq!(bulk, [mh.clone(), other]);
    assert_eq!(
        enc.encode_mols_bulk(&[mol_from_smiles("CCO").unwrap()])
            .unwrap(),
        [mh]
    );

    let default = MHFPEncoder::default();
    let mh = default.encode_smiles("CCO").unwrap();
    assert_eq!(mh.len(), 2048);
    assert_eq!(mh[..3], [998603975, 113829445, 99459723]);
}

#[test]
fn shinglings() {
    let enc = encoder();
    assert_eq!(
        enc.create_shingling_from_smiles("CCO").unwrap(),
        ["CC", "CCO", "C(C)O", "OC", "OCC"]
    );
    assert_eq!(
        enc.create_shingling_from_mol_with_params(
            &mol_from_smiles("c1ccccc1O").unwrap(),
            &EncodeParams {
                radius: 1,
                ..Default::default()
            },
        )
        .unwrap(),
        [
            "C1=CC=CC=C1",
            "C(=C)C",
            "C(=C)C",
            "C(=C)C",
            "C(=C)C",
            "C(=C)C",
            "C(=C)(C)O",
            "OC"
        ]
    );
    assert_eq!(
        enc.create_shingling_from_smiles_with_params(
            "CCO",
            &EncodeParams {
                min_radius: 0,
                rings: false,
                ..Default::default()
            },
        )
        .unwrap(),
        ["C", "C", "O", "CC", "CCO", "C(C)O", "OC", "OCC"]
    );
}

#[test]
fn secfp() {
    let enc = encoder();
    let fp = enc
        .encode_secfp_smiles_with_params("CCO", &secfp64())
        .unwrap();
    assert_eq!(fp.get_num_bits(), 64);
    assert_eq!(fp.get_on_bits(), [4, 7, 12, 15, 35]);
    let mol = mol_from_smiles("CCO").unwrap();
    assert_eq!(
        enc.encode_secfp_mol_with_params(&mol, &secfp64())
            .unwrap()
            .get_on_bits(),
        [4, 7, 12, 15, 35]
    );
    let bulk = enc
        .encode_secfp_smiles_bulk_with_params(&["CCO", "CCN"], &secfp64())
        .unwrap();
    assert_eq!(bulk[1].get_on_bits(), [15, 22, 23, 35, 60]);
    let bulk = enc
        .encode_secfp_mols_bulk_with_params(&[mol_from_smiles("CCN").unwrap()], &secfp64())
        .unwrap();
    assert_eq!(bulk[0].get_on_bits(), [15, 22, 23, 35, 60]);
    assert_eq!(
        MHFPEncoder::default()
            .encode_secfp_smiles("CCO")
            .unwrap()
            .get_on_bits(),
        [419, 972, 1223, 1668, 1999]
    );
}

#[test]
fn arguments_that_rdkit_trusts() {
    let enc = encoder();
    // CreateShingling asserts that the SMILES parsed.
    assert!(matches!(enc.encode_smiles("C(("), Err(Error::Parse { .. })));
    assert!(matches!(
        enc.encode_secfp_smiles("C(("),
        Err(Error::Parse { .. })
    ));
    // Distance reads the second MinHash at every index of the first.
    assert!(matches!(
        MHFPEncoder::distance(&[1, 2, 3], &[1]),
        Err(Error::Value { .. })
    ));
    assert!(MHFPEncoder::distance(&[], &[]).unwrap().is_nan());
    // SECFP folds modulo the length, taken modulo 2^32.
    for length in [0, 1 << 32, i32::MAX as usize + 1] {
        let params = EncodeSECFPParams {
            length,
            ..Default::default()
        };
        assert!(matches!(
            enc.encode_secfp_smiles_with_params("CCO", &params),
            Err(Error::Value { .. })
        ));
    }
    assert!(matches!(
        enc.encode_mol(
            &mol_from_smiles_with_params(
                "c1cccc1",
                &SmilesParserParams {
                    sanitize: false,
                    ..Default::default()
                },
            )
            .unwrap()
        ),
        Err(Error::Kekulize { .. })
    ));
}

#[test]
fn params_defaults() {
    assert_eq!(
        EncodeParams::default(),
        EncodeParams {
            radius: 3,
            rings: true,
            isomeric: false,
            kekulize: true,
            min_radius: 1,
        }
    );
    assert_eq!(EncodeSECFPParams::default().length, 2048);
    assert_eq!(
        MHFPEncoderParams::default(),
        MHFPEncoderParams {
            n_permutations: 2048,
            seed: 42,
        }
    );
}
