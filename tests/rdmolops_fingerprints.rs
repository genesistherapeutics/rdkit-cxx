// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;
use rdkit_cxx::data_structs::ExplicitBitVect;

fn smiles(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

#[test]
fn rdk_fingerprints() {
    let mol = smiles("c1ccccc1CCO");
    let fp = rdk_fingerprint(&mol).unwrap();
    assert_eq!(fp.get_num_bits(), 2048);
    assert_eq!(fp.get_num_on_bits(), 85);
    assert_eq!(fp.get_on_bits()[..5], [103, 105, 161, 194, 195]);

    let short = RDKFingerprintParams {
        max_path: 2,
        ..Default::default()
    };
    let (fp, outputs) = rdk_fingerprint_with_outputs(&mol, &short).unwrap();
    assert_eq!(fp.get_num_on_bits(), 16);
    assert_eq!(outputs.atom_bits.len(), 9);
    assert_eq!(outputs.atom_bits[0][..4], [842, 1026, 294, 330]);
    assert_eq!(outputs.bit_info.len(), 16);
    assert_eq!(
        outputs.bit_info[&294],
        [[0, 1], [0, 8], [1, 2], [2, 3], [3, 4], [4, 8]]
    );
    assert_eq!(rdk_fingerprint_with_params(&mol, &short).unwrap(), fp);

    let from = RDKFingerprintParams {
        from_atoms: Some(vec![7]),
        ..Default::default()
    };
    assert_eq!(
        rdk_fingerprint_with_params(&mol, &from)
            .unwrap()
            .get_num_on_bits(),
        48
    );
    let bad = RDKFingerprintParams {
        from_atoms: Some(vec![70]),
        ..Default::default()
    };
    assert_eq!(
        rdk_fingerprint_with_params(&mol, &bad).unwrap_err(),
        Error::Index { index: 70, len: 9 }
    );
    // RDKit's own preconditions check these.
    let short_invariants = RDKFingerprintParams {
        atom_invariants: Some(vec![1, 2]),
        ..Default::default()
    };
    assert!(matches!(
        rdk_fingerprint_with_params(&mol, &short_invariants),
        Err(Error::Runtime { .. })
    ));
    let huge = RDKFingerprintParams {
        fp_size: u32::MAX,
        ..Default::default()
    };
    assert!(matches!(
        rdk_fingerprint_with_params(&mol, &huge),
        Err(Error::Value { .. })
    ));

    let unfolded = unfolded_rdk_fingerprint_count_based(&mol).unwrap();
    assert_eq!(unfolded.get_nonzero_elements().len(), 43);
    assert_eq!(unfolded.get_total_val().unwrap(), 92);
    let (_, outputs) =
        unfolded_rdk_fingerprint_count_based_with_outputs(&mol, &Default::default()).unwrap();
    assert_eq!(outputs.atom_bits.len(), 9);
}

#[test]
fn layered_and_pattern_fingerprints() {
    let mol = smiles("c1ccccc1CCO");
    assert_eq!(layered_fingerprint(&mol).unwrap().get_num_on_bits(), 144);
    let counts = LayeredFingerprintParams {
        atom_counts: Some(vec![0; 9]),
        ..Default::default()
    };
    let (_, outputs) = layered_fingerprint_with_outputs(&mol, &counts).unwrap();
    assert_eq!(outputs.atom_counts, [64, 55, 52, 55, 64, 79, 60, 39, 17]);
    assert_eq!(layered_fingerprint_substruct_layers(), 7);
    let short_counts = LayeredFingerprintParams {
        atom_counts: Some(vec![0; 2]),
        ..Default::default()
    };
    assert!(layered_fingerprint_with_params(&mol, &short_counts).is_err());

    let pattern = pattern_fingerprint(&mol).unwrap();
    assert_eq!(pattern.get_num_on_bits(), 86);
    let mut only = ExplicitBitVect::new(2048);
    only.set_bits_from_list(&pattern.get_on_bits()[..10])
        .unwrap();
    let with_only = PatternFingerprintParams {
        set_only_bits: Some(only),
        ..Default::default()
    };
    pattern_fingerprint_with_params(&mol, &with_only).unwrap();
    // The vector must have fp_size bits.
    let wrong_size = PatternFingerprintParams {
        set_only_bits: Some(ExplicitBitVect::new(16)),
        ..Default::default()
    };
    assert!(matches!(
        pattern_fingerprint_with_params(&mol, &wrong_size),
        Err(Error::Runtime { .. })
    ));
    let tautomers = PatternFingerprintParams {
        tautomer_fingerprints: true,
        ..Default::default()
    };
    assert_eq!(
        pattern_fingerprint_with_params(&smiles("Oc1ccccn1"), &tautomers)
            .unwrap()
            .get_num_on_bits(),
        117
    );

    let mut bundle = MolBundle::default();
    // An empty bundle has no fingerprint (Python's None).
    assert!(pattern_fingerprint_mol_bundle(&bundle).unwrap().is_none());
    bundle.add_mol(&smiles("CCO"));
    bundle.add_mol(&smiles("CCN"));
    assert_eq!(
        pattern_fingerprint_mol_bundle(&bundle)
            .unwrap()
            .unwrap()
            .get_num_on_bits(),
        11
    );
}

#[test]
fn param_defaults() {
    let rdk = RDKFingerprintParams::default();
    assert_eq!(
        (
            rdk.min_path,
            rdk.max_path,
            rdk.fp_size,
            rdk.n_bits_per_hash,
            rdk.min_size
        ),
        (1, 7, 2048, 2, 128)
    );
    assert!(rdk.use_hs && rdk.branched_paths && rdk.use_bond_order && rdk.tgt_density == 0.0);
    let layered = LayeredFingerprintParams::default();
    assert_eq!(layered.layer_flags, 0xFFFFFFFF);
    assert_eq!(
        (layered.min_path, layered.max_path, layered.fp_size),
        (1, 7, 2048)
    );
    let pattern = PatternFingerprintParams::default();
    assert_eq!(pattern.fp_size, 2048);
    assert!(!pattern.tautomer_fingerprints);
}
