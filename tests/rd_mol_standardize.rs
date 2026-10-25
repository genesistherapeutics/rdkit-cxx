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
fn fragment_remover() {
    let default = FragmentRemover::default();
    let no_last = FragmentRemover::new_with_params(&FragmentRemoverParams {
        leave_last: false,
        ..Default::default()
    })
    .unwrap();
    let skip = FragmentRemover::new_with_params(&FragmentRemoverParams {
        skip_if_all_match: true,
        ..Default::default()
    })
    .unwrap();
    for (input, expected) in [
        ("CN(C)C.Cl.Cl.Br", ["CN(C)C", "CN(C)C", "CN(C)C"]),
        ("Cl.Br", ["Br", "", "Br.Cl"]),
        ("c1ccccc1.Cl", ["c1ccccc1", "", "Cl.c1ccccc1"]),
        ("[Na+].[Cl-]", ["[Na+]", "", "[Cl-].[Na+]"]),
    ] {
        let m = mol(input);
        for (remover, want) in [&default, &no_last, &skip].into_iter().zip(expected) {
            assert_eq!(smiles(&remover.remove(&m).unwrap()), want, "{input}");
        }
    }
    let mut m = mol("Cl.Br.CCO");
    default.remove_in_place(&mut m).unwrap();
    assert_eq!(smiles(&m), "CCO");

    let data = "// comment\nhydrogen\t[H]\nchloride\t[Cl]\n";
    let from_data = fragment_remover_from_data(data).unwrap();
    assert_eq!(
        smiles(&from_data.remove(&mol("CN(C)C.Cl.Br")).unwrap()),
        "Br.CN(C)C"
    );
    assert_eq!(smiles(&from_data.remove(&mol("Cl.Cl")).unwrap()), "Cl.Cl");
    let all = fragment_remover_from_data_with_params(
        data,
        &FragmentRemoverFromDataParams {
            leave_last: false,
            skip_if_all_match: false,
        },
    )
    .unwrap();
    assert_eq!(smiles(&all.remove(&mol("Cl.Cl")).unwrap()), "");

    assert!(is_value(fragment_remover_from_data("x\t[C")));
    assert!(is_value(fragment_remover_from_data("onlyname")));
    assert!(matches!(
        FragmentRemover::new_with_params(&FragmentRemoverParams {
            fragment_filename: PathBuf::from("/nonexistent"),
            ..Default::default()
        }),
        Err(Error::Os { .. })
    ));
}

#[test]
fn largest_fragment_chooser() {
    let heavy = CleanupParameters {
        largest_fragment_chooser_count_heavy_atoms_only: true,
        ..Default::default()
    };
    let no_count = CleanupParameters {
        largest_fragment_chooser_use_atom_count: false,
        ..Default::default()
    };
    let choosers = [
        LargestFragmentChooser::default(),
        LargestFragmentChooser::new_with_prefer_organic(true),
        LargestFragmentChooser::from_params(&heavy),
        LargestFragmentChooser::from_params(&no_count),
    ];
    for (input, expected) in [
        (
            "O=S(=O)([O-])[O-].C",
            [
                "O=S(=O)([O-])[O-]",
                "C",
                "O=S(=O)([O-])[O-]",
                "O=S(=O)([O-])[O-]",
            ],
        ),
        ("CCCCC.c1ccccc1", ["CCCCC", "CCCCC", "c1ccccc1", "c1ccccc1"]),
        (
            "[Na+].CC(=O)[O-]",
            ["CC(=O)[O-]", "CC(=O)[O-]", "CC(=O)[O-]", "CC(=O)[O-]"],
        ),
        ("O.O.CC", ["CC", "CC", "CC", "CC"]),
        (
            "[2H]C([2H])([2H])[2H].CC",
            ["CC", "CC", "[2H]C([2H])([2H])[2H]", "CC"],
        ),
        ("C.N", ["C", "C", "N", "N"]),
    ] {
        let m = mol(input);
        for (chooser, want) in choosers.iter().zip(expected) {
            assert_eq!(smiles(&chooser.choose(&m).unwrap()), want, "{input}");
        }
    }
    let mut m = mol("CCCCC.c1ccccc1");
    choosers[0].choose_in_place(&mut m).unwrap();
    assert_eq!(smiles(&m), "CCCCC");
    assert_eq!(
        LargestFragmentChooser::from_params(&CleanupParameters::default()),
        LargestFragmentChooser::default()
    );
}

#[test]
fn reionizer() {
    assert_eq!(
        charge_corrections(),
        [
            ChargeCorrection::new("[Li,Na,K]", "[Li,Na,K;X0+0]", 1),
            ChargeCorrection::new("[Mg,Ca]", "[Mg,Ca;X0+0]", 2),
            ChargeCorrection::new("[Cl]", "[Cl;X0+0]", -1),
        ]
    );
    let cc = ChargeCorrection::new("x", "[C]", 3);
    assert_eq!(
        (cc.name.as_str(), cc.smarts.as_str(), cc.charge),
        ("x", "[C]", 3)
    );

    let default = Reionizer::default();
    let no_corrections = Reionizer::new_with_ccs("", &[]).unwrap();
    let negative_sodium =
        Reionizer::new_with_ccs("", &[ChargeCorrection::new("x", "[Na]", -1)]).unwrap();
    for (input, expected) in [
        (
            "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O",
            [
                "O=S(O)c1ccc(S(=O)(=O)[O-])cc1",
                "O=S(O)c1ccc(S(=O)(=O)[O-])cc1",
                "O=S(O)c1ccc(S(=O)(=O)[O-])cc1",
            ],
        ),
        (
            "[Na].CC(=O)O",
            ["CC(=O)[O-].[Na+]", "CC(=O)O.[Na]", "CC(=O)O.[Na-]"],
        ),
        (
            "[Mg].OC(=O)CC(=O)O",
            [
                "O=C([O-])CC(=O)[O-].[Mg+2]",
                "O=C(O)CC(=O)O.[Mg]",
                "O=C(O)CC(=O)O.[Mg]",
            ],
        ),
        ("[Cl].CN", ["CN.[Cl-]", "CN.[Cl]", "CN.[Cl]"]),
    ] {
        let m = mol(input);
        for (reionizer, want) in [&default, &no_corrections, &negative_sodium]
            .into_iter()
            .zip(expected)
        {
            assert_eq!(smiles(&reionizer.reionize(&m).unwrap()), want, "{input}");
        }
    }
    let mut m = mol("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O");
    default.reionize_in_place(&mut m).unwrap();
    assert_eq!(smiles(&m), "O=S(O)c1ccc(S(=O)(=O)[O-])cc1");
    assert_eq!(
        smiles(
            &Reionizer::new("")
                .unwrap()
                .reionize(&mol("[Na].CC(=O)O"))
                .unwrap()
        ),
        "CC(=O)[O-].[Na+]"
    );

    let data = "// acid base\nsulfonic\tOS(=O)(=O)[OH]\tOS(=O)(=O)[O-]\n\
                carboxylic\tC(=O)[OH]\tC(=O)[O-]\n";
    let from_data = reionizer_from_data(data).unwrap();
    assert_eq!(
        smiles(
            &from_data
                .reionize(&mol("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O"))
                .unwrap()
        ),
        "O=S([O-])c1ccc(S(=O)(=O)O)cc1"
    );
    // Python's default is no charge corrections.
    assert_eq!(
        smiles(&from_data.reionize(&mol("[Na].CC(=O)O")).unwrap()),
        "CC(=O)O.[Na]"
    );
    let with_corrections =
        reionizer_from_data_with_charge_corrections(data, &charge_corrections()).unwrap();
    assert_eq!(
        smiles(&with_corrections.reionize(&mol("[Na].CC(=O)O")).unwrap()),
        "CC(=O)[O-].[Na+]"
    );

    assert!(matches!(
        Reionizer::new("/nonexistent"),
        Err(Error::Os { .. })
    ));
    assert!(is_value(reionizer_from_data("x\t[C\tC")));
}

#[test]
fn charge_corrections_that_do_not_parse_are_refused() {
    // RDKit parses each correction's SMARTS as it reionizes and dereferences
    // the null that SMARTS that do not parse give.
    let bad = [ChargeCorrection::new("x", "[C", 1)];
    assert!(is_value(Reionizer::new_with_ccs("", &bad)));
    assert!(is_value(reionizer_from_data_with_charge_corrections(
        "carboxylic\tC(=O)[OH]\tC(=O)[O-]\n",
        &bad
    )));
}

#[test]
fn uncharger() {
    let unchargers = [
        Uncharger::default(),
        Uncharger::new_with_params(&UnchargerParams {
            force: true,
            ..Default::default()
        }),
        Uncharger::new_with_params(&UnchargerParams {
            protonation_only: true,
            ..Default::default()
        }),
        Uncharger::new_with_params(&UnchargerParams {
            canonical_order: false,
            ..Default::default()
        }),
    ];
    for (input, expected) in [
        ("C[NH3+].O=C([O-])c1ccccc1", ["CN.O=C(O)c1ccccc1"; 4]),
        (
            "C[N+](C)(C)C.O=C([O-])c1ccccc1.O=C([O-])c1ccccc1",
            [
                "C[N+](C)(C)C.O=C(O)c1ccccc1.O=C([O-])c1ccccc1",
                "C[N+](C)(C)C.O=C(O)c1ccccc1.O=C(O)c1ccccc1",
                "C[N+](C)(C)C.O=C(O)c1ccccc1.O=C([O-])c1ccccc1",
                "C[N+](C)(C)C.O=C(O)c1ccccc1.O=C([O-])c1ccccc1",
            ],
        ),
        (
            "C[N+](C)(C)CC([O-])=O",
            [
                "C[N+](C)(C)CC(=O)[O-]",
                "C[N+](C)(C)CC(=O)O",
                "C[N+](C)(C)CC(=O)[O-]",
                "C[N+](C)(C)CC(=O)[O-]",
            ],
        ),
        ("CC[CH2+]", ["CCC", "CCC", "[CH2+]CC", "CCC"]),
    ] {
        let m = mol(input);
        for (uncharger, want) in unchargers.iter().zip(expected) {
            assert_eq!(smiles(&uncharger.uncharge(&m).unwrap()), want, "{input}");
        }
    }
    let mut m = mol("C[NH3+].O=C([O-])c1ccccc1");
    unchargers[0].uncharge_in_place(&mut m).unwrap();
    assert_eq!(smiles(&m), "CN.O=C(O)c1ccccc1");
}

#[test]
fn large_molecules_are_refused() {
    // One molecule for every check, to build it once. RDKit sums formal
    // charges in int (U9): the reionizer's corrections can set charges of
    // 128 in magnitude, and it subtracts its sums, so it refuses more than
    // i32::MAX / 255 atoms, as do the functions that clean up; the uncharger
    // only shrinks charges, so it refuses charges whose magnitudes sum past
    // i32::MAX, as 2^24 of -128 do. The tautomer score adds up to 382 per
    // atom, besides its ring and substructure terms, so the tautomer
    // functions refuse more than about 5.6 million atoms.
    let mut atom = Atom::from_num(6);
    atom.set_formal_charge(-128);
    let mut big = RWMol::default();
    for _ in 0..1 << 24 {
        big.add_atom(&atom);
    }
    assert!(is_value(Uncharger::default().uncharge(&big)));
    assert!(is_value(Uncharger::default().uncharge_in_place(&mut big)));
    assert!(is_value(Reionizer::default().reionize(&big)));
    assert!(is_value(Reionizer::default().reionize_in_place(&mut big)));
    assert!(is_value(cleanup(&big)));
    assert!(is_value(reionize_in_place(&mut big)));
    assert!(is_value(charge_parent(&big)));
    assert!(is_value(fragment_parent(&big)));
    assert!(is_value(canonical_tautomer(&big)));
    let mut mols = vec![mol("CC"), Mol::from(big)];
    assert!(is_value(super_parent_in_place_mols(&mut mols, 2)));
    // Checked before any molecule is changed.
    assert_eq!(mols[0].get_num_atoms(), 2);
}
