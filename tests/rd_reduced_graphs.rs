// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_reduced_graphs::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdReducedGraphs in RDKit 2024.03.6.

fn nonzero(fp: &[f64]) -> Vec<(usize, f64)> {
    fp.iter()
        .enumerate()
        .filter(|&(_, &v)| v != 0.0)
        .map(|(i, &v)| (i, (v * 1e6).round() / 1e6))
        .collect()
}

const ERG: [(usize, f64); 28] = [
    (3, 0.3),
    (4, 1.0),
    (5, 0.3),
    (17, 0.6),
    (18, 2.0),
    (19, 0.6),
    (31, 0.3),
    (32, 1.0),
    (33, 0.3),
    (70, 0.3),
    (71, 1.0),
    (72, 0.6),
    (73, 1.0),
    (74, 0.3),
    (87, 0.3),
    (88, 1.0),
    (89, 0.3),
    (101, 0.3),
    (102, 1.0),
    (103, 0.3),
    (140, 0.3),
    (141, 1.0),
    (142, 0.6),
    (143, 1.0),
    (144, 0.3),
    (196, 0.3),
    (197, 1.0),
    (198, 0.3),
];

#[test]
fn reduced_graphs_and_erg() {
    let mol = mol_from_smiles("OCCc1ccccc1N").unwrap();
    let graph = generate_mol_extended_reduced_graph(&mol).unwrap();
    assert_eq!(graph.get_num_atoms(), 7);
    assert_eq!(mol_to_smiles(&graph).unwrap(), "Nc1c(CCO)*1");
    let fp = generate_erg_fingerprint_for_reduced_graph(&graph).unwrap();
    assert_eq!(fp.len(), 315);
    assert_eq!(nonzero(&fp), ERG);
    assert_eq!(nonzero(&get_erg_fingerprint(&mol).unwrap()), ERG);
    let fp = get_erg_fingerprint_with_params(
        &mol,
        &GetErGFingerprintParams {
            fuzz_increment: 0.0,
            min_path: 2,
            max_path: 4,
        },
    )
    .unwrap();
    assert_eq!(fp.len(), 63);
    assert_eq!(
        nonzero(&fp),
        [(10, 1.0), (12, 1.0), (20, 1.0), (22, 1.0), (28, 1.0)]
    );
}

#[test]
fn arguments_that_rdkit_trusts() {
    let mol = mol_from_smiles("OCCc1ccccc1N").unwrap();
    // RDKit asserts minPath < maxPath and counts both as unsigned.
    for (min_path, max_path) in [(4, 4), (-1, 4), (0, i32::MAX)] {
        let params = GetErGFingerprintParams {
            min_path,
            max_path,
            ..Default::default()
        };
        assert!(
            matches!(
                get_erg_fingerprint_with_params(&mol, &params),
                Err(Error::Value { .. })
            ),
            "{min_path} {max_path}"
        );
    }
    // A molecule that is not a reduced graph has no types.
    assert!(matches!(
        generate_erg_fingerprint_for_reduced_graph(&mol),
        Err(Error::Key { .. })
    ));
}
