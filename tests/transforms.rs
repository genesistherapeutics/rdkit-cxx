// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn smarts(text: &str) -> Mol {
    mol_from_smarts(text).unwrap()
}

fn to_smiles(mol: &Mol) -> String {
    mol_to_smiles(mol).unwrap()
}

// Expected values below come from Python's rdmolops.

#[test]
fn delete_and_replace_substructs() {
    assert_eq!(
        to_smiles(&delete_substructs(&smiles("CC(=O)O.Cl"), &smarts("Cl")).unwrap()),
        "CC(=O)O"
    );
    let only_frags = DeleteSubstructsParams {
        only_frags: true,
        ..Default::default()
    };
    assert_eq!(
        to_smiles(
            &delete_substructs_with_params(&smiles("CCl.Cl"), &smarts("Cl"), &only_frags).unwrap()
        ),
        "CCl"
    );
    let results = replace_substructs(&smiles("OCCO"), &smarts("O"), &smiles("N")).unwrap();
    assert_eq!(
        results.iter().map(to_smiles).collect::<Vec<_>>(),
        vec!["NCCO", "NCCO"]
    );
    let all = ReplaceSubstructsParams {
        replace_all: true,
        ..Default::default()
    };
    let results =
        replace_substructs_with_params(&smiles("OCCO"), &smarts("O"), &smiles("N"), &all).unwrap();
    assert_eq!(
        results.iter().map(to_smiles).collect::<Vec<_>>(),
        vec!["NCCN"]
    );
    let bad_point = ReplaceSubstructsParams {
        replacement_connection_point: 5,
        ..Default::default()
    };
    assert!(
        replace_substructs_with_params(&smiles("OCCO"), &smarts("O"), &smiles("N"), &bad_point)
            .is_err()
    );
}

#[test]
fn replace_cores_and_sidechains() {
    let mol = smiles("c1ccccc1CCO");
    let core = smiles("c1ccccc1");
    assert_eq!(
        to_smiles(&replace_sidechains(&mol, &core).unwrap().unwrap()),
        "[1*]c1ccccc1"
    );
    assert!(replace_sidechains(&mol, &smiles("N")).unwrap().is_none());
    assert_eq!(
        to_smiles(&replace_core_query(&mol, &core).unwrap().unwrap()),
        "[1*]CCO"
    );
    assert!(replace_core_query(&mol, &smiles("N")).unwrap().is_none());
    let matches: Vec<(u32, u32)> = (0..6).map(|i| (i, i)).collect();
    assert_eq!(
        to_smiles(&replace_core(&mol, &core, &matches).unwrap().unwrap()),
        "[1*]CCO"
    );
    let by_index = ReplaceCoreParams {
        label_by_index: true,
        ..Default::default()
    };
    assert_eq!(
        to_smiles(
            &replace_core_with_params(&mol, &core, &matches, &by_index)
                .unwrap()
                .unwrap()
        ),
        "[5*]CCO"
    );
    for bad in [[(0, 99)], [(99, 0)], [(u32::MAX, 0)]] {
        assert!(matches!(
            replace_core(&mol, &core, &bad),
            Err(Error::Value { .. })
        ));
    }
    assert_eq!(
        ReplaceCoreQueryParams::default(),
        ReplaceCoreQueryParams {
            replace_dummies: true,
            label_by_index: false,
            require_dummy_match: false,
            use_chirality: false,
        }
    );
}

#[test]
fn murcko() {
    assert_eq!(
        to_smiles(&murcko_decompose(&smiles("c1ccccc1CCC1CC1CC")).unwrap()),
        "c1ccc(CCC2CC2)cc1"
    );
    // The scaffold computation leaves no stale matrix on its input.
    let mut mol = RWMol::from(smiles("C1CC1CC"));
    murcko_decompose(&mol).unwrap();
    let idx = mol.add_atom(&Atom::new("N").unwrap());
    mol.add_bond_with_order(4, idx, BondType::Single).unwrap();
    assert_eq!(get_distance_matrix(&mol).unwrap().len(), 6);
}

#[test]
fn combine() {
    let methane = mol_from_mol_block(&mol_to_mol_block(&smiles("C")).unwrap()).unwrap();
    let combined = combine_mols(&methane, &methane).unwrap();
    assert_eq!(to_smiles(&combined), "C.C");
    let moved = combine_mols_with_offset(&methane, &methane, Point3D::new(1.0, 2.0, 3.0)).unwrap();
    assert_eq!(
        moved.get_conformer().unwrap().get_atom_position(1).unwrap(),
        Point3D::new(1.0, 2.0, 3.0)
    );
}

#[test]
fn core_substitution() {
    let mol = smiles("CCC");
    let core = smiles("C*");
    let matches = vec![vec![0, 1], vec![1, 2]];
    assert_eq!(
        get_most_substituted_core_match(&mol, &core, &matches).unwrap(),
        vec![0, 1]
    );
    assert_eq!(
        sort_matches_by_degree_of_core_substitution(&mol, &core, &matches).unwrap(),
        matches
    );
    assert!(matches!(
        get_most_substituted_core_match(&mol, &core, &[]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        sort_matches_by_degree_of_core_substitution(&mol, &core, &[vec![]]),
        Err(Error::Value { .. })
    ));
    assert!(get_most_substituted_core_match(&mol, &core, &[vec![0, 9]]).is_err());
    assert!(get_most_substituted_core_match(&mol, &core, &[vec![0, 1, 2]]).is_err());
}

#[test]
fn recursive_queries() {
    let mut query = smarts("[C]");
    add_recursive_query(&mut query, &smiles("CO"), 0).unwrap();
    assert_eq!(mol_to_smarts(&query).unwrap(), "[C&$([#6]-[#8])]");
    assert_eq!(
        smiles("CCO").get_substruct_matches(&query).unwrap(),
        vec![vec![1]]
    );
    assert!(matches!(
        add_recursive_query(&mut query, &smiles("CO"), 1),
        Err(Error::Value { .. })
    ));
    let mut plain = smiles("CC");
    add_recursive_query_with_preserve_existing_query(&mut plain, &smiles("CO"), 0, false).unwrap();
    assert!(plain.has_query());

    let mut mol = smiles("CC");
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("query", "carbonyl");
    let queries = BTreeMap::from([("carbonyl".to_owned(), smarts("C=O"))]);
    mol_add_recursive_queries(&mut mol, &queries, "query").unwrap();
    assert_eq!(mol_to_smarts(&mol).unwrap(), "[#6&$(C=O)]-[#6]");
    mol.get_atom_with_idx_mut(1)
        .unwrap()
        .set_prop("query", "missing(");
    assert!(matches!(
        mol_add_recursive_queries(&mut mol, &queries, "query"),
        Err(Error::Key { .. })
    ));
}

#[test]
fn query_def_files() {
    let text = "// comment\nAcid\tC(=O)[OH]\nBad\tC(\nAmine\t[NX3]\n";
    let defs = parse_mol_query_def_file_text(text).unwrap();
    let as_smarts: BTreeMap<_, _> = defs
        .iter()
        .map(|(k, v)| (k.as_str(), mol_to_smarts(v).unwrap()))
        .collect();
    assert_eq!(
        as_smarts,
        BTreeMap::from([
            ("acid", "C(=O)[O&H1]".to_owned()),
            ("amine", "[N&X3]".to_owned())
        ])
    );
    let keep_case = ParseMolQueryDefFileParams {
        standardize: false,
        delimiter: ",".to_owned(),
        comment: "#".to_owned(),
        name_column: 1,
        smarts_column: 0,
    };
    let defs = parse_mol_query_def_file_text_with_params("#x\nC,Methyl\n", &keep_case).unwrap();
    assert_eq!(defs.keys().collect::<Vec<_>>(), vec!["Methyl"]);

    let dir = std::env::temp_dir().join(format!("rdkit-cxx-qdef-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("defs.txt");
    std::fs::write(&path, text).unwrap();
    assert_eq!(parse_mol_query_def_file(&path).unwrap().len(), 2);
    assert_eq!(
        parse_mol_query_def_file_with_params(&path, &ParseMolQueryDefFileParams::default())
            .unwrap()
            .len(),
        2
    );
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(matches!(
        parse_mol_query_def_file(dir.join("missing.txt")),
        Err(Error::Os { .. })
    ));
    assert!(matches!(
        parse_mol_query_def_file("defs.txt\0.bak"),
        Err(Error::Value { .. })
    ));
}

#[test]
fn fragment_bonds() {
    let mol = smiles("CCOC");
    let (frags, outputs) =
        fragment_on_bonds_with_outputs(&mol, &[1], &FragmentOnBondsParams::default()).unwrap();
    assert_eq!(to_smiles(&frags), "[1*]OC.[2*]CC");
    assert_eq!(outputs.cuts_per_atom, vec![0, 1, 1, 0]);
    let labeled = FragmentOnBondsParams {
        dummy_labels: Some(vec![(10, 20)]),
        bond_types: Some(vec![BondType::Double]),
        ..Default::default()
    };
    assert_eq!(
        to_smiles(&fragment_on_bonds_with_params(&mol, &[1], &labeled).unwrap()),
        "[10*]=OC.[20*]=CC"
    );
    let bare = FragmentOnBondsParams {
        add_dummies: false,
        ..Default::default()
    };
    assert_eq!(
        to_smiles(&fragment_on_bonds_with_params(&mol, &[1], &bare).unwrap()),
        "CC.CO"
    );
    assert!(matches!(
        fragment_on_bonds(&mol, &[]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        fragment_on_bonds(&mol, &[1, 1]),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        fragment_on_bonds(&mol, &[3]),
        Err(Error::Index { .. })
    ));
    for bad in [
        FragmentOnBondsParams {
            dummy_labels: Some(vec![]),
            ..Default::default()
        },
        FragmentOnBondsParams {
            bond_types: Some(vec![BondType::Single; 2]),
            ..Default::default()
        },
    ] {
        assert!(matches!(
            fragment_on_bonds_with_params(&mol, &[1], &bad),
            Err(Error::Value { .. })
        ));
    }
}

#[test]
fn fragment_some_bonds() {
    let mol = smiles("CCOC");
    let (frags, cuts) = fragment_on_some_bonds_return_cuts_per_atom(&mol, &[0, 2]).unwrap();
    assert_eq!(
        frags.iter().map(to_smiles).collect::<Vec<_>>(),
        vec!["*COC.[1*]C", "[2*]C.[3*]OCC"]
    );
    assert_eq!(cuts, vec![vec![1, 1, 0, 0], vec![0, 0, 1, 1]]);
    let two = FragmentOnSomeBondsParams {
        num_to_break: 2,
        ..Default::default()
    };
    assert_eq!(
        fragment_on_some_bonds_with_params(&mol, &[0, 1, 2], &two)
            .unwrap()
            .iter()
            .map(to_smiles)
            .collect::<Vec<_>>(),
        vec![
            "*C[2*].[1*]C.[1*]OC",
            "*CO[3*].[1*]C.[2*]C",
            "[1*]O[3*].[2*]C.[2*]CC"
        ]
    );
    assert_eq!(fragment_on_some_bonds(&mol, &[0]).unwrap().len(), 1);
    let (_, cuts) = fragment_on_some_bonds_return_cuts_per_atom_with_params(
        &mol,
        &[0],
        &FragmentOnSomeBondsParams::default(),
    )
    .unwrap();
    assert_eq!(cuts.len(), 1);
    assert!(matches!(
        fragment_on_some_bonds(&mol, &[9]),
        Err(Error::Index { .. })
    ));
    assert_eq!(
        FragmentOnSomeBondsParams::default(),
        FragmentOnSomeBondsParams {
            num_to_break: 1,
            add_dummies: true,
            dummy_labels: None,
            bond_types: None,
        }
    );
}

#[test]
fn brics() {
    assert_eq!(
        to_smiles(&fragment_on_brics_bonds(&smiles("CCCOCCC(=O)c1ccccc1")).unwrap()),
        "[16*]c1ccccc1.[3*]O[3*].[4*]CCC.[4*]CCC([6*])=O"
    );
}

#[test]
fn molzip_pairs() {
    assert_eq!(
        to_smiles(
            &molzip(&smiles("C[*:1]"), &smiles("N[*:1]"))
                .unwrap()
                .unwrap()
        ),
        "CN"
    );
    assert_eq!(
        to_smiles(&molzip_single(&smiles("C[*:1].N[*:1]")).unwrap().unwrap()),
        "CN"
    );
    let isotope = MolzipParams {
        label: MolzipLabel::Isotope,
        ..Default::default()
    };
    assert_eq!(
        to_smiles(
            &molzip_with_params(&smiles("C[1*]"), &smiles("N[1*]"), &isotope)
                .unwrap()
                .unwrap()
        ),
        "CN"
    );
    let mut atom_type = MolzipParams {
        label: MolzipLabel::AtomType,
        ..Default::default()
    };
    atom_type.set_atom_symbols(&["Xe", "Kr"]);
    assert_eq!(
        to_smiles(
            &molzip_with_params(&smiles("C[Xe]"), &smiles("N[Xe]"), &atom_type)
                .unwrap()
                .unwrap()
        ),
        "CN"
    );
    // RDKit 2024.09 cannot zip a linker of dummies alone: Python raises
    // RDKit's precondition ("no owner").
    assert!(matches!(
        molzip(&smiles("[*:1][*:2]"), &smiles("C[*:1].N[*:2]")),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(
        MolzipParams::default(),
        MolzipParams {
            label: MolzipLabel::AtomMapNumber,
            atom_symbols: vec![],
            enforce_valence_rules: true,
            generate_coordinates: false,
        }
    );
}

#[test]
fn molzip_fragment_on_bonds() {
    let fob = MolzipParams {
        label: MolzipLabel::FragmentOnBonds,
        ..Default::default()
    };
    let frags = fragment_on_bonds(&smiles("CCOC"), &[1]).unwrap();
    assert_eq!(
        to_smiles(&molzip_single_with_params(&frags, &fob).unwrap().unwrap()),
        "CCOC"
    );
    // The dummy's isotope names an atom with no matching dummy: RDKit gives
    // up and returns null.
    assert!(
        molzip_single_with_params(&smiles("C[1*]"), &fob)
            .unwrap()
            .is_none()
    );
    // get_other_atom returns null for a dummy without one neighbor, and
    // RDKit dereferences it.
    for bad in ["C.[0*]", "C[0*]C"] {
        assert!(matches!(
            molzip_single_with_params(&smiles(bad), &fob),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            molzip_with_params(&smiles("C"), &smiles(bad), &fob),
            Err(Error::Value { .. })
        ));
    }
}

#[test]
fn molzip_lists() {
    let mut mols = vec![smiles("c1ccccc1[*:1]"), smiles("O[*:1]")];
    assert_eq!(
        to_smiles(&molzip_fragments(&mut mols).unwrap().unwrap()),
        "Oc1ccccc1"
    );
    assert!(molzip_fragments(&mut []).unwrap().is_none());
    let coords = MolzipParams {
        generate_coordinates: true,
        ..Default::default()
    };
    let zipped = molzip_fragments_with_params(&mut mols, &coords)
        .unwrap()
        .unwrap();
    assert_eq!(zipped.get_num_conformers(), 1);
    assert!(mols.iter().all(|m| m.get_num_conformers() == 1));

    let mut row = BTreeMap::from([
        ("Core".to_owned(), smiles("c1ccccc1[*:1]")),
        ("R1".to_owned(), smiles("O[*:1]")),
    ]);
    assert_eq!(
        to_smiles(&molzip_row(&mut row).unwrap().unwrap()),
        "Oc1ccccc1"
    );
    let zipped = molzip_row_with_params(&mut row, &coords).unwrap().unwrap();
    assert_eq!(zipped.get_num_conformers(), 1);
    let mut no_core = BTreeMap::from([("R1".to_owned(), smiles("O[*:1]"))]);
    assert!(matches!(molzip_row(&mut no_core), Err(Error::Value { .. })));
}

#[test]
fn molzip_lists_that_rdkit_crashes_on() {
    // Zipping a list dereferences a null dummy for every linker, and a null
    // neighbor for a labeled atom without exactly one; Python segfaults on
    // both.
    let mut linker = vec![smiles("[*:1][*:2]"), smiles("C[*:1]"), smiles("N[*:2]")];
    assert!(matches!(
        molzip_fragments(&mut linker),
        Err(Error::Value { .. })
    ));
    let coords = MolzipParams {
        generate_coordinates: true,
        ..Default::default()
    };
    let mut lonely = vec![smiles("C[*:1]"), smiles("[*:1]")];
    assert!(matches!(
        molzip_fragments_with_params(&mut lonely, &coords),
        Err(Error::Value { .. })
    ));
    let mut row = BTreeMap::from([
        ("Core".to_owned(), smiles("C[*:1]")),
        ("R1".to_owned(), smiles("[*:1][*:2]")),
    ]);
    assert!(matches!(molzip_row(&mut row), Err(Error::Value { .. })));
}

#[test]
fn molzip_coordinates_refuse_what_the_depictor_would_mishandle() {
    // The zipped molecule keeps a square-planar atom without neighbors,
    // whose first neighbor RDKit's depictor reads unchecked; Python
    // segfaults generating its coordinates.
    let mut planar = smiles("[Pt].C[*:1]");
    planar
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiSquareplanar);
    let mut mols = vec![planar, smiles("[*:1]O")];
    let coords = MolzipParams {
        generate_coordinates: true,
        ..Default::default()
    };
    assert!(matches!(
        molzip_fragments_with_params(&mut mols, &coords),
        Err(Error::Value { message }) if message.contains("square-planar")
    ));
    assert!(mols.iter().all(|mol| mol.get_num_conformers() == 0));
    let zipped = molzip_fragments(&mut mols).unwrap().unwrap();
    assert_eq!(mol_to_smiles(&zipped).unwrap(), "CO.[Pt@SP]");
}
