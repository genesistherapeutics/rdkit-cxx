// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_reduced_graphs::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdReducedGraphs.

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
    (18, 0.6),
    (19, 2.0),
    (20, 0.6),
    (33, 0.3),
    (34, 1.0),
    (35, 0.3),
    (75, 0.3),
    (76, 1.0),
    (77, 0.6),
    (78, 1.0),
    (79, 0.3),
    (93, 0.3),
    (94, 1.0),
    (95, 0.3),
    (108, 0.3),
    (109, 1.0),
    (110, 0.3),
    (150, 0.3),
    (151, 1.0),
    (152, 0.6),
    (153, 1.0),
    (154, 0.3),
    (210, 0.3),
    (211, 1.0),
    (212, 0.3),
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
        [(15, 1.0), (17, 1.0), (30, 1.0), (32, 1.0), (42, 1.0)]
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
