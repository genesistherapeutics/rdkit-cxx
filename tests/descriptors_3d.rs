// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports testGet3DMolDescriptors of RDKit's rdkit/Chem/UnitTestDescriptors.py,
// with its reference file tests/data/descriptors_multiconf.sdf
// (rdkit/Chem/test_data), and pins values from Python's
// rdkit.Chem.Descriptors3D in RDKit 2025.03.6.

use std::collections::BTreeMap;

use rdkit_cxx::chem::descriptors_3d::*;
use rdkit_cxx::chem::rd_mol_descriptors::{CalcNPR1Params, CalcPMI1Params};
use rdkit_cxx::chem::{
    Conformer, Mol, MolFromMolBlockParams, mol_from_mol_block_with_params, mol_from_smiles,
};

/// Agrees with Python's value to the last bits RDKit's eigensolvers may
/// change across platforms.
fn assert_close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{what}: {actual} != {expected}"
    );
}

/// Python's `Chem.MultiConfMolFromSDF(..., removeHs=False)`: the first
/// record, with every record's conformer.
fn multi_conf_mol() -> Mol {
    let params = MolFromMolBlockParams {
        remove_hs: false,
        ..MolFromMolBlockParams::default()
    };
    let mut records = include_str!("data/descriptors_multiconf.sdf")
        .split("$$$$\n")
        .filter(|record| !record.trim().is_empty())
        .map(|record| mol_from_mol_block_with_params(record, &params).unwrap());
    let mut mol = records.next().unwrap();
    for other in records {
        let conf = Conformer::from(other.get_conformer().unwrap());
        mol.add_conformer_with_assign_id(&conf, true).unwrap();
    }
    mol
}

const MULTICONF_DEFAULT: [(&str, f64); 11] = [
    ("PMI1", 20.959498779728786),
    ("PMI2", 69.57141747216187),
    ("PMI3", 85.07741076486279),
    ("NPR1", 0.24635797670966636),
    ("NPR2", 0.8177425340839717),
    ("RadiusOfGyration", 1.299039315208758),
    ("InertialShapeFactor", 0.039015366859576836),
    ("Eccentricity", 0.9691789036661496),
    ("Asphericity", 0.43547565691919127),
    ("SpherocityIndex", 0.1000854885995871),
    ("PBF", 0.2223884029596845),
];

const PROPANOL: [(&str, f64); 11] = [
    ("PMI1", 20.9582722049221),
    ("PMI2", 69.5741688498054),
    ("PMI3", 85.07823500818624),
    ("NPR1", 0.24634117295575644),
    ("NPR2", 0.8177669511257605),
    ("RadiusOfGyration", 1.2990480035584586),
    ("InertialShapeFactor", 0.03901881525012858),
    ("Eccentricity", 0.9691831748987298),
    ("Asphericity", 0.43550226984521134),
    ("SpherocityIndex", 0.10009742302790876),
    ("PBF", 0.2224031689355117),
];

fn propanol() -> Mol {
    mol_from_smiles(concat!(
        "CCCO |(1.44534,-0.585581,0.158885;0.667797,0.646552,-0.278384;",
        "-0.741018,0.544094,0.296045;-1.37212,-0.605065,-0.176546)|"
    ))
    .unwrap()
}

#[test]
fn calc_mol_descriptors_3d_matches_python() {
    // Without coordinates.
    assert!(calc_mol_descriptors_3d(&mol_from_smiles("CCCO").unwrap()).is_err());

    let names: Vec<&str> = desc_list().iter().map(|(name, _)| *name).collect();
    assert_eq!(names, PROPANOL.map(|(name, _)| name));

    let descs = calc_mol_descriptors_3d(&propanol()).unwrap();
    // In descList's order, as Python's dictionary has them.
    let calculated: Vec<&str> = descs.iter().map(|(name, _)| *name).collect();
    assert_eq!(calculated, names);
    let descs: BTreeMap<&str, f64> = descs.into_iter().collect();
    assert!((descs["PMI1"] - 20.9583).abs() <= 1e-4);
    for (name, expected) in PROPANOL {
        assert_close(descs[name], expected, name);
    }

    let mol = multi_conf_mol();
    assert_eq!(mol.get_num_conformers(), 10);
    let descs: BTreeMap<&str, f64> = calc_mol_descriptors_3d(&mol).unwrap().into_iter().collect();
    assert!((descs["PMI1"] - 20.9595).abs() <= 1e-4);
    for (name, expected) in MULTICONF_DEFAULT {
        assert_close(descs[name], expected, name);
    }
}

#[test]
fn functions_match_python() {
    let m = propanol();
    assert_close(pmi1(&m).unwrap(), PROPANOL[0].1, "PMI1");
    assert_close(pmi2(&m).unwrap(), PROPANOL[1].1, "PMI2");
    assert_close(pmi3(&m).unwrap(), PROPANOL[2].1, "PMI3");
    assert_close(npr1(&m).unwrap(), PROPANOL[3].1, "NPR1");
    assert_close(npr2(&m).unwrap(), PROPANOL[4].1, "NPR2");
    assert_close(
        radius_of_gyration(&m).unwrap(),
        PROPANOL[5].1,
        "RadiusOfGyration",
    );
    assert_close(inertial_shape_factor(&m).unwrap(), PROPANOL[6].1, "ISF");
    assert_close(eccentricity(&m).unwrap(), PROPANOL[7].1, "Eccentricity");
    assert_close(asphericity(&m).unwrap(), PROPANOL[8].1, "Asphericity");
    assert_close(
        spherocity_index(&m).unwrap(),
        PROPANOL[9].1,
        "SpherocityIndex",
    );
    assert_close(pbf(&m).unwrap(), PROPANOL[10].1, "PBF");
    let unweighted = CalcPMI1Params {
        use_atomic_masses: false,
        ..CalcPMI1Params::default()
    };
    assert_close(
        pmi1_with_params(&m, &unweighted).unwrap(),
        1.6436518727564156,
        "PMI1",
    );
    let unweighted = CalcNPR1Params {
        use_atomic_masses: false,
        ..CalcNPR1Params::default()
    };
    assert_close(
        npr1_with_params(&m, &unweighted).unwrap(),
        0.2571896642281666,
        "NPR1",
    );
    assert_close(
        spherocity_index_with_conf_id(&m, -1).unwrap(),
        PROPANOL[9].1,
        "conf",
    );
    assert_close(pbf_with_conf_id(&m, -1).unwrap(), PROPANOL[10].1, "conf");
    assert!(pbf_with_conf_id(&m, 3).is_err());
}
