// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::path::PathBuf;

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
use rdkit_cxx::chem::{Atom, Mol, RWMol, mol_from_smiles, mol_to_smiles};

// Expected values come from Python's rdMolStandardize.

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn smiles(mol: &Mol) -> String {
    mol_to_smiles(mol).unwrap()
}

fn is_value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

type Copy = fn(&Mol) -> rdkit_cxx::Result<Mol>;

const FUNCTIONS: [(&str, Copy); 11] = [
    ("cleanup", cleanup),
    ("normalize", normalize),
    ("reionize", reionize),
    ("remove_fragments", remove_fragments),
    ("canonical_tautomer", canonical_tautomer),
    ("tautomer_parent", tautomer_parent),
    ("fragment_parent", fragment_parent),
    ("stereo_parent", stereo_parent),
    ("isotope_parent", isotope_parent),
    ("charge_parent", charge_parent),
    ("super_parent", super_parent),
];

// Each input with the SMILES of FUNCTIONS' results.
const CASES: [(&str, [&str; 11]); 9] = [
    (
        "CC(=O)[O-].[Na+]",
        [
            "CC(=O)[O-].[Na+]",
            "CC(=O)[O-].[Na+]",
            "CC(=O)[O-].[Na+]",
            "CC(=O)[O-]",
            "CC(=O)[O-].[Na+]",
            "CC(=O)[O-].[Na+]",
            "CC(=O)[O-]",
            "CC(=O)[O-].[Na+]",
            "CC(=O)[O-].[Na+]",
            "CC(=O)O",
            "CC(=O)O",
        ],
    ),
    (
        "[Na]OC(=O)c1ccccc1",
        [
            "O=C([O-])c1ccccc1.[Na+]",
            "O=C([O][Na])c1ccccc1",
            "O=C([O][Na])c1ccccc1",
            "O=C([O][Na])c1ccccc1",
            "O=C([O][Na])c1ccccc1",
            "O=C([O-])c1ccccc1.[Na+]",
            "O=C([O-])c1ccccc1",
            "O=C([O-])c1ccccc1.[Na+]",
            "O=C([O-])c1ccccc1.[Na+]",
            "O=C(O)c1ccccc1",
            "O=C(O)c1ccccc1",
        ],
    ),
    ("C[S+](C)[O-]", ["C[S+](C)[O-]"; 11]),
    (
        "Oc1ccccn1",
        [
            "Oc1ccccn1",
            "Oc1ccccn1",
            "Oc1ccccn1",
            "Oc1ccccn1",
            "O=c1cccc[nH]1",
            "O=c1cccc[nH]1",
            "Oc1ccccn1",
            "Oc1ccccn1",
            "Oc1ccccn1",
            "Oc1ccccn1",
            "O=c1cccc[nH]1",
        ],
    ),
    (
        "C[C@H](N)C(=O)O",
        [
            "C[C@H](N)C(=O)O",
            "C[C@H](N)C(=O)O",
            "C[C@H](N)C(=O)O",
            "C[C@H](N)C(=O)O",
            "CC(N)C(=O)O",
            "CC(N)C(=O)O",
            "C[C@H](N)C(=O)O",
            "CC(N)C(=O)O",
            "C[C@H](N)C(=O)O",
            "C[C@H](N)C(=O)O",
            "CC(N)C(=O)O",
        ],
    ),
    (
        "[13CH3]CO",
        [
            "[13CH3]CO",
            "[13CH3]CO",
            "[13CH3]CO",
            "[13CH3]CO",
            "[13CH3]CO",
            "[13CH3]CO",
            "[13CH3]CO",
            "[13CH3]CO",
            "CCO",
            "[13CH3]CO",
            "CCO",
        ],
    ),
    (
        "[Na+].[O-]C(=O)[C@@H](N)Cc1c[nH]cn1",
        [
            "N[C@@H](Cc1c[nH]cn1)C(=O)[O-].[Na+]",
            "N[C@@H](Cc1c[nH]cn1)C(=O)[O-].[Na+]",
            "N[C@@H](Cc1c[nH]cn1)C(=O)[O-].[Na+]",
            "N[C@@H](Cc1c[nH]cn1)C(=O)[O-]",
            "NC(Cc1c[nH]cn1)C(=O)[O-].[Na+]",
            "NC(Cc1c[nH]cn1)C(=O)[O-].[Na+]",
            "N[C@@H](Cc1c[nH]cn1)C(=O)[O-]",
            "NC(Cc1c[nH]cn1)C(=O)[O-].[Na+]",
            "N[C@@H](Cc1c[nH]cn1)C(=O)[O-].[Na+]",
            "N[C@@H](Cc1c[nH]cn1)C(=O)O",
            "NC(Cc1c[nH]cn1)C(=O)O",
        ],
    ),
    (
        "O=C([O-])c1ccccc1.C[NH3+]",
        [
            "C[NH3+].O=C([O-])c1ccccc1",
            "C[NH3+].O=C([O-])c1ccccc1",
            "C[NH3+].O=C([O-])c1ccccc1",
            "O=C([O-])c1ccccc1",
            "C[NH3+].O=C([O-])c1ccccc1",
            "C[NH3+].O=C([O-])c1ccccc1",
            "O=C([O-])c1ccccc1",
            "C[NH3+].O=C([O-])c1ccccc1",
            "C[NH3+].O=C([O-])c1ccccc1",
            "O=C(O)c1ccccc1",
            "O=C(O)c1ccccc1",
        ],
    ),
    (
        "CCO.Cl",
        [
            "CCO.Cl", "CCO.Cl", "CCO.Cl", "CCO", "CCO.Cl", "CCO.Cl", "CCO", "CCO.Cl", "CCO.Cl",
            "CCO", "CCO",
        ],
    ),
];

#[test]
fn module_functions() {
    for (input, expected) in CASES {
        let m = mol(input);
        for ((name, func), want) in FUNCTIONS.iter().zip(expected) {
            assert_eq!(smiles(&func(&m).unwrap()), want, "{name}({input})");
        }
        // The input is unchanged.
        assert_eq!(smiles(&m), smiles(&mol(input)));
    }
}

#[test]
fn in_place_functions() {
    type InPlace = fn(&mut Mol) -> rdkit_cxx::Result<()>;
    let functions: [(usize, InPlace); 10] = [
        (0, cleanup_in_place),
        (2, reionize_in_place),
        (3, remove_fragments_in_place),
        (5, tautomer_parent_in_place),
        (6, fragment_parent_in_place),
        (7, stereo_parent_in_place),
        (8, isotope_parent_in_place),
        (9, charge_parent_in_place),
        (10, super_parent_in_place),
        (1, normalize_in_place),
    ];
    for (input, expected) in CASES {
        for (i, func) in functions {
            let mut m = mol(input);
            func(&mut m).unwrap();
            // Normalizing in place transforms the whole molecule, not each
            // fragment, but these agree.
            assert_eq!(smiles(&m), expected[i], "{}({input})", FUNCTIONS[i].0);
        }
    }
}

#[test]
fn skipping_the_cleanup() {
    type Parent = fn(&Mol, bool) -> rdkit_cxx::Result<Mol>;
    let parents: [Parent; 6] = [
        |m, skip| {
            tautomer_parent_with_params(
                m,
                &TautomerParentParams {
                    skip_standardize: skip,
                    ..Default::default()
                },
            )
        },
        |m, skip| {
            fragment_parent_with_params(
                m,
                &FragmentParentParams {
                    skip_standardize: skip,
                    ..Default::default()
                },
            )
        },
        |m, skip| {
            stereo_parent_with_params(
                m,
                &StereoParentParams {
                    skip_standardize: skip,
                    ..Default::default()
                },
            )
        },
        |m, skip| {
            isotope_parent_with_params(
                m,
                &IsotopeParentParams {
                    skip_standardize: skip,
                    ..Default::default()
                },
            )
        },
        |m, skip| {
            charge_parent_with_params(
                m,
                &ChargeParentParams {
                    skip_standardize: skip,
                    ..Default::default()
                },
            )
        },
        |m, skip| {
            super_parent_with_params(
                m,
                &SuperParentParams {
                    skip_standardize: skip,
                    ..Default::default()
                },
            )
        },
    ];
    let m = mol("[Na]OC(=O)c1ccccc1");
    let skipped = [
        "O=C([O-])c1ccccc1.[Na+]",
        "O=C([O][Na])c1ccccc1",
        "O=C([O][Na])c1ccccc1",
        "O=C([O][Na])c1ccccc1",
        "O=C([O-])c1ccccc1.[Na+]",
        "O=C([O-])c1ccccc1.[Na+]",
    ];
    for (i, (parent, want)) in parents.iter().zip(skipped).enumerate() {
        assert_eq!(smiles(&parent(&m, true).unwrap()), want, "parent {i}");
        assert_eq!(
            smiles(&parent(&m, false).unwrap()),
            CASES[1].1[5 + i],
            "parent {i}"
        );
    }
    let mut m = mol("C[C@H](N)C(=O)O");
    tautomer_parent_in_place_with_params(
        &mut m,
        &TautomerParentInPlaceParams {
            skip_standardize: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(smiles(&m), "CC(N)C(=O)O");
}

#[test]
fn several_molecules() {
    // RDKit's own threads loop forever for 0 threads and race for negative
    // counts; these split the molecules among the threads they run.
    for num_threads in [1, 0, -1, 2, 3, 16] {
        let mut mols: Vec<Mol> = CASES.iter().map(|(input, _)| mol(input)).collect();
        super_parent_in_place_mols(&mut mols, num_threads).unwrap();
        for (m, (input, expected)) in mols.iter().zip(CASES) {
            assert_eq!(smiles(m), expected[10], "{input} on {num_threads} threads");
        }
        let mut mols: Vec<Mol> = CASES.iter().map(|(input, _)| mol(input)).collect();
        cleanup_in_place_mols(&mut mols, num_threads).unwrap();
        normalize_in_place_mols(&mut mols, num_threads).unwrap();
        for (m, (input, expected)) in mols.iter().zip(CASES) {
            assert_eq!(smiles(m), expected[0], "{input} on {num_threads} threads");
        }
    }
    let mut none: Vec<Mol> = Vec::new();
    cleanup_in_place_mols(&mut none, 4).unwrap();

    type Mols = fn(&mut [Mol], i32) -> rdkit_cxx::Result<()>;
    let functions: [(usize, Mols); 9] = [
        (2, reionize_in_place_mols),
        (3, remove_fragments_in_place_mols),
        (5, tautomer_parent_in_place_mols),
        (6, fragment_parent_in_place_mols),
        (7, stereo_parent_in_place_mols),
        (8, isotope_parent_in_place_mols),
        (9, charge_parent_in_place_mols),
        (0, cleanup_in_place_mols),
        (10, super_parent_in_place_mols),
    ];
    for (i, func) in functions {
        let mut mols: Vec<Mol> = CASES.iter().map(|(input, _)| mol(input)).collect();
        func(&mut mols, 4).unwrap();
        for (m, (input, expected)) in mols.iter().zip(CASES) {
            assert_eq!(smiles(m), expected[i], "{}({input})", FUNCTIONS[i].0);
        }
    }
}

#[test]
fn errors_in_threads_are_returned() {
    // RDKit's threads let the exception escape, which ends the process.
    let params = CleanupParameters {
        normalizations_file: PathBuf::from("/nonexistent/normalizations.txt"),
        ..Default::default()
    };
    for num_threads in [1, 2, 4] {
        let mut mols: Vec<Mol> = (0..8).map(|_| mol("CC")).collect();
        assert!(matches!(
            cleanup_in_place_mols_with_params(&mut mols, num_threads, &params),
            Err(Error::Os { .. })
        ));
    }
    assert!(matches!(
        cleanup_with_params(&mol("CC"), &params),
        Err(Error::Os { .. })
    ));
}

#[test]
fn standardizing_smiles() {
    assert_eq!(
        standardize_smiles("[Na]OC(=O)c1ccccc1").unwrap(),
        "O=C([O-])c1ccccc1.[Na+]"
    );
    assert_eq!(standardize_smiles("C[S+](C)[O-]").unwrap(), "C[S+](C)[O-]");
    assert!(is_value(standardize_smiles("C1CC")));
}

#[test]
fn cleanup_parameters_defaults() {
    let params = CleanupParameters::default();
    assert_eq!(params.normalizations_file, PathBuf::new());
    assert_eq!(params.acidbase_file, PathBuf::new());
    assert_eq!(params.fragment_file, PathBuf::new());
    assert_eq!(params.tautomer_transforms_file, PathBuf::new());
    assert_eq!(params.max_restarts, 200);
    assert!(!params.prefer_organic);
    assert!(params.do_canonical);
    assert_eq!(params.max_tautomers, 1000);
    assert_eq!(params.max_transforms, 1000);
    assert!(params.tautomer_remove_sp3_stereo);
    assert!(params.tautomer_remove_bond_stereo);
    assert!(params.tautomer_remove_isotopic_hs);
    assert!(params.tautomer_reassign_stereo);
    assert!(params.largest_fragment_chooser_use_atom_count);
    assert!(!params.largest_fragment_chooser_count_heavy_atoms_only);
    assert!(params.normalization_data.is_empty());
    assert!(params.fragment_data.is_empty());
    assert!(params.acidbase_data.is_empty());
    assert!(params.tautomer_transform_data.is_empty());
    for (params, skip) in [
        (
            TautomerParentParams::default().params,
            TautomerParentParams::default().skip_standardize,
        ),
        (
            SuperParentParams::default().params,
            SuperParentParams::default().skip_standardize,
        ),
        (
            ChargeParentInPlaceParams::default().params,
            ChargeParentInPlaceParams::default().skip_standardize,
        ),
    ] {
        assert_eq!(params, CleanupParameters::default());
        assert!(!skip);
    }
}

#[test]
fn parameters_from_json() {
    let mut params = CleanupParameters::default();
    update_params_from_json(
        &mut params,
        r#"{"maxRestarts": 5, "preferOrganic": true, "doCanonical": false,
            "maxTautomers": 3, "maxTransforms": 7, "tautomerRemoveSp3Stereo": false,
            "tautomerRemoveBondStereo": false, "tautomerRemoveIsotopicHs": false,
            "tautomerReassignStereo": false, "normalizations": "/x/n",
            "acidbaseFile": "/x/a", "fragmentFile": "/x/f", "tautomerTransforms": "/x/t"}"#,
    )
    .unwrap();
    assert_eq!(params.normalizations_file, PathBuf::from("/x/n"));
    assert_eq!(params.acidbase_file, PathBuf::from("/x/a"));
    assert_eq!(params.fragment_file, PathBuf::from("/x/f"));
    assert_eq!(params.tautomer_transforms_file, PathBuf::from("/x/t"));
    assert_eq!(params.max_restarts, 5);
    assert!(params.prefer_organic);
    assert!(!params.do_canonical);
    assert_eq!(params.max_tautomers, 3);
    assert_eq!(params.max_transforms, 7);
    assert!(!params.tautomer_remove_sp3_stereo);
    assert!(!params.tautomer_remove_bond_stereo);
    assert!(!params.tautomer_remove_isotopic_hs);
    assert!(!params.tautomer_reassign_stereo);

    // A value that does not convert keeps the current one; empty JSON
    // changes nothing; JSON that does not parse is RuntimeError in Python.
    let before = params.clone();
    update_params_from_json(&mut params, r#"{"maxRestarts": "x"}"#).unwrap();
    update_params_from_json(&mut params, "").unwrap();
    assert_eq!(params, before);
    assert!(matches!(
        update_params_from_json(&mut params, "{bad"),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(params, before);
}

#[test]
fn data_from_json() {
    let mut params = CleanupParameters::default();
    update_params_from_json(
        &mut params,
        r#"{"normalizationData": [{"name": "Sulfoxide to S=O",
             "smarts": "[S+:1][O-:2]>>[S+0:1]=[O-0:2]"}, {"name": "", "smarts": "C>>N"}]}"#,
    )
    .unwrap();
    // Entries without a name are skipped.
    assert_eq!(
        params.normalization_data,
        [(
            "Sulfoxide to S=O".to_string(),
            "[S+:1][O-:2]>>[S+0:1]=[O-0:2]".to_string()
        )]
    );
    assert_eq!(
        smiles(&normalize_with_params(&mol("C[S+](C)[O-]"), &params).unwrap()),
        "CS(C)=O"
    );
    assert_eq!(
        smiles(&normalize_with_params(&mol("C[N+](=O)[O-]"), &params).unwrap()),
        "C[N+](=O)[O-]"
    );
    // Entries are appended.
    update_params_from_json(
        &mut params,
        r#"{"normalizationData": [{"name": "x", "smarts": "[C:1]>>[C:1]"}]}"#,
    )
    .unwrap();
    assert_eq!(params.normalization_data.len(), 2);

    let mut params = CleanupParameters::default();
    update_params_from_json(
        &mut params,
        r#"{"fragmentData": [{"name": "hydrogen", "smarts": "[H]"},
             {"name": "methanol", "smarts": "[CH3][OH]"}]}"#,
    )
    .unwrap();
    let m = mol("CO.CCO.Cl");
    assert_eq!(
        smiles(&remove_fragments_with_params(&m, &params).unwrap()),
        "CCO.Cl"
    );
    assert_eq!(smiles(&remove_fragments(&m).unwrap()), "CO");

    let mut params = CleanupParameters::default();
    update_params_from_json(
        &mut params,
        r#"{"acidbaseData": [{"name": "-OSO3H", "acid": "OS(=O)(=O)[OH]",
             "base": "OS(=O)(=O)[O-]"}, {"name": "carbox", "acid": "C(=O)[OH]",
             "base": "C(=O)[O-]"}]}"#,
    )
    .unwrap();
    assert_eq!(params.acidbase_data.len(), 2);
    let m = mol("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O");
    assert_eq!(
        smiles(&reionize_with_params(&m, &params).unwrap()),
        "O=S([O-])c1ccc(S(=O)(=O)O)cc1"
    );
    assert_eq!(
        smiles(&reionize(&m).unwrap()),
        "O=S(O)c1ccc(S(=O)(=O)[O-])cc1"
    );

    let mut params = CleanupParameters::default();
    update_params_from_json(
        &mut params,
        r#"{"tautomerTransformData": [{"name": "1,3 (thio)keto/enol f",
             "smarts": "[CX4!H0]-[C]=[O,S,Se,Te;X1]", "bonds": "", "charges": ""}]}"#,
    )
    .unwrap();
    assert_eq!(
        smiles(&canonical_tautomer_with_params(&mol("Oc1ccccn1"), &params).unwrap()),
        "Oc1ccccn1"
    );
    assert_eq!(
        smiles(&canonical_tautomer_with_params(&mol("CC(C)=O"), &params).unwrap()),
        "CC(C)=O"
    );
}

#[test]
fn prefer_organic() {
    let m = mol("O=S(=O)([O-])[O-].C");
    let params = FragmentParentParams {
        params: CleanupParameters {
            prefer_organic: true,
            ..Default::default()
        },
        skip_standardize: true,
    };
    assert_eq!(
        smiles(&fragment_parent_with_params(&m, &params).unwrap()),
        "C"
    );
    assert_eq!(
        smiles(
            &fragment_parent_with_params(
                &m,
                &FragmentParentParams {
                    skip_standardize: true,
                    ..Default::default()
                }
            )
            .unwrap()
        ),
        "O=S(=O)([O-])[O-]"
    );
}

#[test]
fn short_tautomer_transforms_are_refused() {
    // RDKit indexes a transform's bonds and charges by its SMARTS's bonds
    // and atoms without checking their lengths, and reads past their ends.
    let short_bonds = (
        "x".to_string(),
        "[CX4!H0]-[C]=[O]".to_string(),
        "=".to_string(),
        String::new(),
    );
    let short_charges = (
        "x".to_string(),
        "[CX4!H0]-[C]=[O]".to_string(),
        "=-".to_string(),
        "+".to_string(),
    );
    let m = mol("CC(C)=O");
    for transform in [short_bonds, short_charges] {
        let params = CleanupParameters {
            tautomer_transform_data: vec![transform],
            ..Default::default()
        };
        assert!(is_value(canonical_tautomer_with_params(&m, &params)));
        let parent = TautomerParentParams {
            params: params.clone(),
            skip_standardize: false,
        };
        assert!(is_value(tautomer_parent_with_params(&m, &parent)));
        let parent = SuperParentInPlaceParams {
            params,
            skip_standardize: false,
        };
        let mut mols = vec![mol("CC(C)=O"), mol("CCO")];
        assert!(is_value(super_parent_in_place_mols_with_params(
            &mut mols, 2, &parent
        )));
        // Refused before any molecule is changed.
        assert_eq!(smiles(&mols[0]), "CC(C)=O");
    }
    // Lists as long as the SMARTS's bonds and atoms are used.
    let params = CleanupParameters {
        tautomer_transform_data: vec![(
            "x".to_string(),
            "[CX4!H0]-[C]=[O]".to_string(),
            "=-".to_string(),
            "000".to_string(),
        )],
        ..Default::default()
    };
    assert_eq!(
        smiles(&canonical_tautomer_with_params(&m, &params).unwrap()),
        "CC(C)=O"
    );
}

#[test]
fn tautomer_transforms_from_a_file() {
    let path = std::env::temp_dir().join(format!("rdkit-cxx-{}-tautomers.in", std::process::id()));
    std::fs::write(
        &path,
        "// Name\tSMARTS\tBonds\tCharges\n1,3 (thio)keto/enol f\t[CX4!H0]-[C]=[O,S,Se,Te;X1]\n",
    )
    .unwrap();
    let params = CleanupParameters {
        tautomer_transforms_file: path.clone(),
        ..Default::default()
    };
    assert_eq!(
        smiles(&canonical_tautomer_with_params(&mol("Oc1ccccn1"), &params).unwrap()),
        "Oc1ccccn1"
    );
    std::fs::write(&path, "x\t[CX4!H0]-[C]=[O]\t=\n").unwrap();
    assert!(is_value(canonical_tautomer_with_params(
        &mol("CC(C)=O"),
        &params
    )));
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(
        canonical_tautomer_with_params(&mol("CC(C)=O"), &params),
        Err(Error::Os { .. })
    ));
}

#[test]
fn file_names_with_nul_are_refused() {
    let params = CleanupParameters {
        acidbase_file: PathBuf::from("acid\0base"),
        ..Default::default()
    };
    assert!(is_value(reionize_with_params(&mol("CC"), &params)));
}

#[test]
fn large_molecules_are_refused() {
    // RDKit's reionizer sums charges, which its corrections can set to 128
    // in magnitude, in int, and subtracts the sums, so more than
    // i32::MAX / 255 atoms are refused (U9).
    let atom = Atom::from_num(6);
    let mut big = RWMol::default();
    for _ in 0..i32::MAX as u32 / 255 + 1 {
        big.add_atom(&atom);
    }
    assert!(is_value(cleanup(&big)));
    assert!(is_value(reionize_in_place(&mut big)));
    assert!(is_value(charge_parent(&big)));
    assert!(is_value(fragment_parent(&big)));
    let mut mols = vec![mol("CC"), Mol::from(big)];
    assert!(is_value(super_parent_in_place_mols(&mut mols, 2)));
    // Checked before any molecule is changed.
    assert_eq!(mols[0].get_num_atoms(), 2);
}

#[test]
fn large_tautomer_scores_are_refused() {
    // RDKit's tautomer score adds up to 382 per atom for hydrogens on P, S,
    // Se, and Te in int, besides its ring and substructure terms.
    let atom = Atom::from_num(6);
    let mut big = RWMol::default();
    for _ in 0..(i32::MAX as u32 - 47000) / 383 + 1 {
        big.add_atom(&atom);
    }
    assert!(is_value(canonical_tautomer(&big)));
}
