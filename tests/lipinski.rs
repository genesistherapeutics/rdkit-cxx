// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.Lipinski in RDKit 2025.09.6.

use rdkit_cxx::Result;
use rdkit_cxx::chem::lipinski::*;
use rdkit_cxx::chem::{Mol, mol_from_smiles};

/// Every Lipinski function, in the order of the expected values, as `f64`.
fn descriptors(m: &Mol) -> Result<Vec<f64>> {
    Ok(vec![
        num_h_donors(m)?.into(),
        num_h_acceptors(m)?.into(),
        num_heteroatoms(m)?.into(),
        num_rotatable_bonds(m)?.into(),
        no_count(m).into(),
        nhoh_count(m)?.into(),
        ring_count(m)?.into(),
        heavy_atom_count(m).into(),
        fraction_csp3(m)?,
        num_aromatic_rings(m).into(),
        num_saturated_rings(m).into(),
        num_aromatic_heterocycles(m).into(),
        num_aromatic_carbocycles(m).into(),
        num_saturated_heterocycles(m).into(),
        num_saturated_carbocycles(m).into(),
        num_aliphatic_rings(m).into(),
        num_aliphatic_heterocycles(m).into(),
        num_aliphatic_carbocycles(m).into(),
        num_heterocycles(m).into(),
        num_bridgehead_atoms(m)?.into(),
        num_amide_bonds(m)?.into(),
        num_atom_stereo_centers(m)?.into(),
        num_unspecified_atom_stereo_centers(m)?.into(),
        num_spiro_atoms(m)?.into(),
        phi(m)?,
    ])
}

#[test]
fn descriptors_match_python() {
    for (smiles, expected) in [
        (
            "CC(=O)Nc1ccc(O)cc1",
            [
                2.0,
                2.0,
                3.0,
                1.0,
                3.0,
                2.0,
                1.0,
                11.0,
                0.125,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                2.1081569077654403,
            ],
        ),
        (
            "C[C@H](N)C(=O)NC1CC2(CCOC2)C1",
            [
                2.0,
                3.0,
                4.0,
                2.0,
                4.0,
                3.0,
                2.0,
                14.0,
                0.9,
                0.0,
                2.0,
                0.0,
                0.0,
                1.0,
                1.0,
                2.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                2.490982930884356,
            ],
        ),
        (
            "C1CC2CCC1C2",
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                2.0,
                7.0,
                1.0,
                0.0,
                2.0,
                0.0,
                0.0,
                0.0,
                2.0,
                2.0,
                0.0,
                2.0,
                0.0,
                2.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.6973140495867769,
            ],
        ),
    ] {
        let m = mol_from_smiles(smiles).unwrap();
        assert_eq!(descriptors(&m).unwrap(), expected, "{smiles}");
    }
}

#[test]
fn patterns_match_python() {
    let m = mol_from_smiles("CCC(CC(C)CC1CCC1)C(CC(=O)O)N").unwrap();
    for (pattern, expected) in [
        (h_donor_smarts(), 2),
        (h_acceptor_smarts(), 2),
        (heteroatom_smarts(), 3),
        (rotatable_bond_smarts(), 8),
        (nhoh_smarts(), 2),
        (no_count_smarts(), 3),
    ] {
        assert_eq!(m.get_substruct_matches(&pattern).unwrap().len(), expected);
    }
}

#[test]
fn issue_2183420() {
    // A problem with the acceptor definition.
    for (smiles, expected) in [
        ("NC", 1),
        ("CNC", 1),
        ("CN(C)C", 1),
        ("NC(=O)", 1),
        ("NC(=O)C", 1),
        ("CNC(=O)", 1),
        ("CNC(=O)C", 1),
        ("O=CNC(=O)C", 2),
        ("O=C(C)NC(=O)C", 2),
    ] {
        let m = mol_from_smiles(smiles).unwrap();
        assert_eq!(num_h_acceptors(&m).unwrap(), expected, "{smiles}");
    }
}
