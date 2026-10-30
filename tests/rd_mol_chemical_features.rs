// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_chemical_features::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::rd_geometry::Point3D;

// Expected values come from Python's rdMolChemicalFeatures.

const FDEF: &str = "
# a comment
AtomType NDonor [N&!H0]
AtomType ODonor [O&!H0]
DefineFeature HDonor1 [{NDonor},{ODonor}]
  Family Donor
  Weights 1.0
EndFeature
DefineFeature Arom6 a1aaaaa1
  Family Aromatic
  Weights 1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
DefineFeature Carbonyl C=O
  Family Acceptor
  Weights 1.0,3.0
EndFeature
";

// OCc1ccccc1C(=O)N, embedded.
const MOL_BLOCK: &str = "
     RDKit          3D

 11 11  0  0  0  0  0  0  0  0999 V2000
   -2.0686    2.3403   -0.4582 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7736    1.8627   -0.3598 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6680    0.4143   -0.2051 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8335   -0.3264   -0.1725 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8403   -1.7037   -0.0295 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6118   -2.3225    0.0817 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5593   -1.5787    0.0490 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5622   -0.2020   -0.0941 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8252    0.5278   -0.1218 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8779    1.7931   -0.2512 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0815   -0.1297   -0.0056 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
  9 11  1  0
  8  3  1  0
M  END
";

fn factory() -> MolChemicalFeatureFactory {
    build_feature_factory_from_string(FDEF).unwrap()
}

fn mol() -> Mol {
    mol_from_mol_block(MOL_BLOCK).unwrap()
}

fn assert_close(p: Point3D, expected: (f64, f64, f64)) {
    // The mol block keeps four decimals of Python's coordinates.
    for (a, b) in [(p.x, expected.0), (p.y, expected.1), (p.z, expected.2)] {
        assert!((a - b).abs() < 1e-3, "{p:?} is not {expected:?}");
    }
}

#[test]
fn factory_definitions() {
    let factory = factory();
    assert_eq!(factory.get_num_feature_defs(), 3);
    assert_eq!(
        factory.get_feature_families(),
        vec!["Donor", "Aromatic", "Acceptor"]
    );
    let expected: BTreeMap<String, String> = [
        ("Donor.HDonor1", "[$([N&!H0]),$([O&!H0])]"),
        ("Aromatic.Arom6", "a1aaaaa1"),
        ("Acceptor.Carbonyl", "C=O"),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_string(), v.to_string()))
    .collect();
    assert_eq!(factory.get_feature_defs(), expected);
    assert!(format!("{factory:?}").contains("Aromatic.Arom6"));
}

#[test]
fn features_of_a_molecule() {
    let factory = factory();
    let mol = mol();
    assert_eq!(factory.get_num_mol_features(&mol).unwrap(), 4);
    assert_eq!(
        factory
            .get_num_mol_features_with_include_only(&mol, "Donor")
            .unwrap(),
        2
    );
    let expected = [
        (
            1,
            "Donor",
            "HDonor1",
            vec![0],
            (-2.068623, 2.340276, -0.458233),
        ),
        (
            2,
            "Donor",
            "HDonor1",
            vec![10],
            (3.081476, -0.129748, -0.005636),
        ),
        (
            3,
            "Aromatic",
            "Arom6",
            vec![2, 3, 4, 5, 6, 7],
            (-0.638669, -0.953158, -0.061753),
        ),
        (
            4,
            "Acceptor",
            "Carbonyl",
            vec![8, 9],
            (1.864755, 1.476758, -0.218835),
        ),
    ];
    for (i, (id, family, type_, atoms, pos)) in expected.into_iter().enumerate() {
        let feat = factory.get_mol_feature(&mol, i as i32).unwrap();
        assert_eq!(feat.get_id(), id);
        assert_eq!(feat.get_family(), family);
        assert_eq!(feat.get_type(), type_);
        assert_eq!(feat.get_atom_ids(), atoms);
        assert_eq!(feat.get_active_conformer(), -1);
        assert_close(feat.get_pos().unwrap(), pos);
        assert_close(feat.get_pos_with_conf_id(-1).unwrap(), pos);
        assert!(std::ptr::eq(feat.get_mol(), &mol));
        assert!(std::ptr::eq(feat.get_factory(), &factory));
    }
    let params = GetMolFeatureParams {
        include_only: "Aromatic".to_string(),
        ..GetMolFeatureParams::default()
    };
    let feat = factory
        .get_mol_feature_with_params(&mol, 0, &params)
        .unwrap();
    assert_eq!(feat.get_family(), "Aromatic");
    assert_eq!(feat.get_id(), 1);
    assert!(matches!(
        factory.get_mol_feature(&mol, 4),
        Err(Error::Index { index: 4, len: 4 })
    ));
    assert!(matches!(
        factory.get_mol_feature(&mol, -1),
        Err(Error::Index { index: -1, len: 4 })
    ));
}

#[test]
fn defaults() {
    let params = GetMolFeatureParams::default();
    assert_eq!(params.include_only, "");
    assert_eq!(params.conf_id, -1);
}

#[test]
fn conformers() {
    let factory = factory();
    let mut mol = mol();
    let mut conf = Conformer::from(mol.get_conformer().unwrap());
    for i in 0..mol.get_num_atoms() {
        let p = conf.get_atom_position(i).unwrap();
        conf.set_atom_position(
            i,
            Point3D {
                x: 2.0 * p.x,
                y: 2.0 * p.y,
                z: 2.0 * p.z,
            },
        )
        .unwrap();
    }
    conf.set_id(5).unwrap();
    mol.add_conformer(&conf).unwrap();
    let params = GetMolFeatureParams {
        conf_id: 5,
        ..GetMolFeatureParams::default()
    };
    let mut feat = factory
        .get_mol_feature_with_params(&mol, 0, &params)
        .unwrap();
    assert_eq!(feat.get_active_conformer(), 5);
    assert_close(feat.get_pos().unwrap(), (-4.137246, 4.680552, -0.916467));
    assert_close(
        feat.get_pos_with_conf_id(0).unwrap(),
        (-2.068623, 2.340276, -0.458233),
    );
    feat.set_active_conformer(0);
    feat.clear_cache();
    assert_close(feat.get_pos().unwrap(), (-2.068623, 2.340276, -0.458233));
    assert!(matches!(
        feat.get_pos_with_conf_id(9),
        Err(Error::Value { .. })
    ));

    // Without conformers there is no position.
    let flat = mol_from_smiles("OCc1ccccc1C(=O)N").unwrap();
    let feat = factory.get_mol_feature(&flat, 0).unwrap();
    assert!(matches!(feat.get_pos(), Err(Error::Runtime { .. })));
    // A feature of several atoms needs weights.
    let factory =
        build_feature_factory_from_string("DefineFeature X [#6][#6]\n Family A\nEndFeature\n")
            .unwrap();
    let plain = mol_from_mol_block(MOL_BLOCK).unwrap();
    let feat = factory.get_mol_feature(&plain, 0).unwrap();
    assert_eq!(feat.get_atom_ids(), vec![1, 2]);
    assert!(matches!(feat.get_pos(), Err(Error::Runtime { .. })));
}

#[test]
fn atom_matches() {
    let factory = factory();
    let mol = mol();
    let feats: Vec<_> = (0..4)
        .map(|i| factory.get_mol_feature(&mol, i).unwrap())
        .collect();
    assert_eq!(
        get_atom_match(&feats[..2]).unwrap(),
        vec![vec![0], vec![10]]
    );
    let shared = [
        factory.get_mol_feature(&mol, 0).unwrap(),
        factory.get_mol_feature(&mol, 0).unwrap(),
    ];
    assert_eq!(get_atom_match(&shared).unwrap(), Vec::<Vec<u32>>::new());
    assert_eq!(
        get_atom_match(&feats).unwrap(),
        vec![vec![0], vec![10], vec![2, 3, 4, 5, 6, 7], vec![8, 9]]
    );
    // RDKit indexes its maxAts flags by atom without a check.
    assert!(matches!(
        get_atom_match_with_max_ats(&feats, 10),
        Err(Error::Index { index: 10, len: 10 })
    ));
    assert!(matches!(
        get_atom_match_with_max_ats(&feats, -1),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        get_atom_match_with_max_ats(&feats[..1], 1).unwrap(),
        vec![vec![0]]
    );
}

#[test]
fn parse_errors() {
    for (text, kind) in [
        ("Foo bar\n", "bad or missing keyword"),
        (
            "DefineFeature X [C\n Family A\nEndFeature\n",
            "invalid SMARTS",
        ),
        ("DefineFeature X C\nEndFeature\n", "did not find Family"),
    ] {
        match build_feature_factory_from_string(text) {
            Err(Error::Parse { message }) => assert!(message.contains(kind), "{message}"),
            other => panic!("{text:?}: {other:?}"),
        }
    }
    assert!(matches!(
        build_feature_factory_from_string(
            "DefineFeature X CC\n Family A\n Weights 1\nEndFeature\n"
        ),
        Err(Error::Value { .. })
    ));
    // Lines that hold only whitespace RDKit's tokenizer drops.
    assert!(matches!(
        build_feature_factory_from_string("\u{b}\n"),
        Err(Error::Parse { .. })
    ));
    // A definition that runs to the end of the input without EndFeature,
    // on which RDKit's parser rereads its last line forever, including
    // when a continuation line swallows the EndFeature.
    for text in [
        "DefineFeature HDonor1 [N,O;!H0]\n  Family HBondDonor3  Wei",
        "DefineFeature D [N,O;!H0]\n Family D",
        "DefineFeature D [N,O;!H0]\n Weights 1.0",
        "DefineFeature D [N,O;!H0]\n Family D \\\nEndFeature",
    ] {
        match build_feature_factory_from_string(text) {
            Err(Error::Parse { message }) => {
                assert!(message.contains("without EndFeature"), "{message}")
            }
            other => panic!("{text:?}: {other:?}"),
        }
    }
    // The same definitions with EndFeature, at the end without a newline.
    let factory =
        build_feature_factory_from_string("DefineFeature D [N,O;!H0]\n Family D\nEndFeature")
            .unwrap();
    assert_eq!(factory.get_num_feature_defs(), 1);
}

#[test]
fn factory_from_a_file() {
    let dir = std::env::temp_dir().join(format!("rdkit-cxx-fdef-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("features.fdef");
    std::fs::write(&path, FDEF).unwrap();
    let factory = build_feature_factory(&path).unwrap();
    assert_eq!(factory.get_num_feature_defs(), 3);
    assert!(matches!(
        build_feature_factory(dir.join("missing.fdef")),
        Err(Error::Os { .. })
    ));
    std::fs::write(&path, b"# \xff\n").unwrap();
    assert!(matches!(
        build_feature_factory(&path),
        Err(Error::Value { .. })
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn shared_factory_across_threads() {
    // The factory is Sync; its matches take the Ctrl-C lock one at a time.
    let factory = factory();
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                let mol = mol();
                for _ in 0..20 {
                    assert_eq!(factory.get_num_mol_features(&mol).unwrap(), 4);
                    let feat = factory.get_mol_feature(&mol, 2).unwrap();
                    assert_eq!(feat.get_atom_ids(), vec![2, 3, 4, 5, 6, 7]);
                }
            });
        }
    });
}
