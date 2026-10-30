// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_partial_charges::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdPartialCharges.

fn charges(mol: &Mol, key: &str) -> Vec<f64> {
    mol.get_atoms()
        .map(|atom| (atom.get_double_prop(key).unwrap() * 1e6).round() / 1e6)
        .collect()
}

#[test]
fn gasteiger_charges() {
    let mut mol = mol_from_smiles("CC(=O)O").unwrap();
    compute_gasteiger_charges(&mut mol).unwrap();
    assert_eq!(
        charges(&mol, "_GasteigerCharge"),
        [0.033768, 0.299685, -0.25282, -0.481433]
    );
    assert_eq!(
        charges(&mol, "_GasteigerHCharge"),
        [0.104545, 0.0, 0.0, 0.296257]
    );
    compute_gasteiger_charges_with_params(
        &mut mol,
        &ComputeGasteigerChargesParams {
            n_iter: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        charges(&mol, "_GasteigerCharge"),
        [0.015978, 0.33588, -0.229022, -0.470062]
    );
}

#[test]
fn missing_parameters() {
    let mut mol = mol_from_smiles("[Xe]C").unwrap();
    compute_gasteiger_charges(&mut mol).unwrap();
    assert!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_double_prop("_GasteigerCharge")
            .unwrap()
            .is_nan()
    );
    let strict = ComputeGasteigerChargesParams {
        throw_on_param_failure: true,
        ..Default::default()
    };
    assert!(matches!(
        compute_gasteiger_charges_with_params(&mut mol, &strict),
        Err(Error::Value { .. })
    ));
}

#[test]
fn arguments_that_rdkit_trusts() {
    let mut mol = mol_from_smiles("CCO").unwrap();
    // RDKit counts the iterations as unsigned.
    let negative = ComputeGasteigerChargesParams {
        n_iter: -1,
        ..Default::default()
    };
    assert!(matches!(
        compute_gasteiger_charges_with_params(&mut mol, &negative),
        Err(Error::Value { .. })
    ));
    let mut unsanitized = mol_from_smiles_with_params(
        "CCO",
        &SmilesParserParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(matches!(
        compute_gasteiger_charges(&mut unsanitized),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn params_defaults() {
    assert_eq!(
        ComputeGasteigerChargesParams::default(),
        ComputeGasteigerChargesParams {
            n_iter: 12,
            throw_on_param_failure: false,
        }
    );
}
