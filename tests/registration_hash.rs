// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/UnitTestRegistrationHash.py, with exact layers
// and hashes from Python's RegistrationHash in RDKit 2024.03.6.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use rdkit_cxx::chem::registration_hash::{
    GetMolLayersParams, HashLayer, HashScheme, get_mol_hash, get_mol_hash_with_hash_scheme,
    get_mol_layers, get_mol_layers_with_params,
};
use rdkit_cxx::chem::{
    BondDir, BondStereo, CXSmilesFields, Mol, SmilesParserParams, mol_from_mol_block,
    mol_from_mol_file, mol_from_smiles, mol_from_smiles_with_params, mol_to_mol_block,
    remove_stereochemistry,
};

type Layers = BTreeMap<HashLayer, String>;

fn layers(expected: [&str; 7]) -> Layers {
    use HashLayer::*;
    [
        CanonicalSmiles,
        Escape,
        Formula,
        NoStereoSmiles,
        NoStereoTautomerHash,
        SgroupData,
        TautomerHash,
    ]
    .into_iter()
    .zip(expected.map(str::to_owned))
    .collect()
}

fn hash_sdf(molblock: &str, params: &GetMolLayersParams) -> Layers {
    get_mol_layers_with_params(&mol_from_mol_block(molblock).unwrap(), params).unwrap()
}

fn smiles_layers(smiles: &str) -> Layers {
    get_mol_layers(&mol_from_smiles(smiles).unwrap()).unwrap()
}

fn unsanitized_layers(smiles: &str) -> Layers {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    get_mol_layers(&mol_from_smiles_with_params(smiles, &params).unwrap()).unwrap()
}

fn data_field_names(names: &[&str]) -> GetMolLayersParams {
    GetMolLayersParams {
        data_field_names: Some(names.iter().map(|name| name.to_string()).collect()),
        ..Default::default()
    }
}

/// Checks that `new_layers` matches `old_layers` except for `expected_diffs`.
fn compare_layers_with_expected_diffs(
    new_layers: &Layers,
    old_layers: &Layers,
    expected_diffs: &[(HashLayer, &str)],
) {
    for (k, old) in old_layers {
        match expected_diffs.iter().find(|(layer, _)| layer == k) {
            Some((_, expected)) => assert_eq!(new_layers[k], *expected),
            None => assert_eq!(new_layers[k], *old),
        }
    }
}

#[test]
fn example_structure() {
    let structure = r#"
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 21 23 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.00838 1.55562 0 0
M  V30 2 C -2.00918 0.555619 0 0
M  V30 3 C -2.87558 0.056219 0 0
M  V30 4 C -2.87638 -0.943781 0 0
M  V30 5 C -3.74278 -1.44298 0 0
M  V30 6 N -2.01078 -1.44438 0 0
M  V30 7 C -1.14438 -0.945181 0 0
M  V30 8 C -0.278781 -1.44578 0 0
M  V30 9 C 0.587619 -0.946381 0 0
M  V30 10 N 1.45322 -1.44718 0 0
M  V30 11 C 2.31962 -0.947781 0 0
M  V30 12 O 3.18522 -1.44858 0 0
M  V30 13 C 2.32042 0.0522191 0 0
M  V30 14 C 1.45482 0.552819 0 0
M  V30 15 C 1.45562 1.55282 0 0
M  V30 16 F 2.45562 1.55202 0 0
M  V30 17 F 1.95642 2.41842 0 0
M  V30 18 F 0.590019 2.05362 0 0
M  V30 19 C 0.588419 0.0536191 0 0
M  V30 20 C -0.277181 0.554219 0 0
M  V30 21 C -1.14358 0.054819 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5 CFG=3
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 2 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 2 11 13
M  V30 13 1 13 14
M  V30 14 1 14 15
M  V30 15 1 15 16
M  V30 16 1 15 17
M  V30 17 1 15 18
M  V30 18 2 14 19
M  V30 19 1 19 20
M  V30 20 2 20 21
M  V30 21 1 21 2
M  V30 22 1 21 7
M  V30 23 1 19 9
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEREL1 ATOMS=(2 2 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END

$$$$
"#;

    let layers1 = hash_sdf(structure, &GetMolLayersParams::default());
    let expected_layers = layers([
        "C[C@H]1C[C@@H](C)c2cc3c(C(F)(F)F)cc(O)nc3cc2N1 |o1:1,3|",
        "",
        "C15H15F3N2O",
        "CC1CC(C)c2cc3c(C(F)(F)F)cc(O)nc3cc2N1",
        "CC1CC(C)[C]2[CH][C]3[C]([CH][C]2[N]1)[N][C]([O])[CH][C]3C(F)(F)F_2_0",
        "[]",
        "C[C@H]1C[C@@H](C)[C]2[CH][C]3[C]([CH][C]2[N]1)[N][C]([O])[CH][C]3C(F)(F)F_2_0 |o1:1,3|",
    ]);
    assert_eq!(layers1, expected_layers);

    // The structure with stereo removed has the same NoStereoSmiles.
    let mut mol = mol_from_mol_block(structure).unwrap();
    remove_stereochemistry(&mut mol).unwrap();
    let stereo_insensitive_layers = hash_sdf(
        &mol_to_mol_block(&mol).unwrap(),
        &GetMolLayersParams::default(),
    );
    assert_eq!(
        stereo_insensitive_layers[&HashLayer::NoStereoSmiles],
        expected_layers[&HashLayer::NoStereoSmiles]
    );
}

#[test]
fn xbhead_xbcorr() {
    // rdkit/Code/GraphMol/FileParsers/sgroup_test_data/repeat_groups_query1.mol
    let structure = r#"
  Mrv1824 06192020192D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.25 0.7484 0 0
M  V30 2 C -2.5837 -0.0216 0 0
M  V30 3 C -2.5837 -1.5617 0 0
M  V30 4 C -1.25 -2.3317 0 0
M  V30 5 C 0.0837 -1.5617 0 0
M  V30 6 C 0.0837 -0.0216 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 1 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(1 1) XBONDS=(2 1 6) XBHEAD=(2 6 1) XBCORR=(4 6 6 1 1) -
M  V30 BRKXYZ=(9 -2.02 -0.0216 0 -2.02 1.5184 0 0 0 0) BRKXYZ=(9 -0.48 1.5184 -
M  V30 0 -0.48 -0.0216 0 0 0 0) CONNECT=HT LABEL="1-3"
M  V30 END SGROUP
M  V30 END CTAB
M  END
"#;

    let expected_layers = layers([
        "C1CCCCC1",
        "",
        "C6H12",
        "C1CCCCC1",
        "C1CCCCC1_0_0",
        r#"[{"type": "SRU", "atoms": [2], "bonds": [[1, 2], [2, 3]], "index": 1, "connect": "HT", "label": "1-3", "XBHEAD": [[0, 5], [1, 0]], "XBCORR": [[[0, 5], [0, 5]], [[1, 0], [1, 0]]]}]"#,
        "C1CCCCC1_0_0",
    ]);
    assert_eq!(
        hash_sdf(structure, &GetMolLayersParams::default()),
        expected_layers
    );

    // The order of the XBHEAD bonds does not matter. (Python's test discards
    // the replaced structure.)
    let reordered = structure.replace("XBHEAD=(2 6 1)", "XBHEAD=(2 1 6)");
    assert_ne!(reordered, structure);
    assert_eq!(
        hash_sdf(&reordered, &GetMolLayersParams::default()),
        expected_layers
    );
}

#[test]
fn data_sgroups() {
    // see SHARED-7879
    let structure = r#"
   Mrv1808 05312108392D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.6458 2.7017 0 0
M  V30 2 C -4.3121 1.9317 0 0
M  V30 3 C -2.9785 2.7017 0 0
M  V30 4 C -1.6448 1.9317 0 0
M  V30 5 C -0.3111 2.7017 0 0
M  V30 6 C 1.0226 1.9317 0 0
M  V30 7 C 2.3563 2.7017 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 7) FIELDNAME="lambda max" FIELDINFO=nm -
M  V30 FIELDDISP="    1.4300   -0.7700    DRU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=250
M  V30 END SGROUP
M  V30 END CTAB
M  END
"#;

    let expected_layers = layers([
        "CCCCCCC",
        "",
        "C7H16",
        "CCCCCCC",
        "CCCCCCC_0_0",
        "[]",
        "CCCCCCC_0_0",
    ]);
    let expected_layers_v2 = layers([
        "CCCCCCC",
        "",
        "C7H16",
        "CCCCCCC",
        "[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]_0_0",
        "[]",
        "[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]_0_0",
    ]);
    assert_eq!(
        hash_sdf(structure, &GetMolLayersParams::default()),
        expected_layers
    );
    let v2 = GetMolLayersParams {
        enable_tautomer_hash_v2: true,
        ..Default::default()
    };
    assert_eq!(hash_sdf(structure, &v2), expected_layers_v2);
}

#[test]
fn stereo_imine_hash() {
    let stereo = unsanitized_layers(r"[H]\N=C\C");
    let no_stereo = unsanitized_layers("N=CC");
    let hydrogen_as_is = unsanitized_layers("N=CC[H]");

    use HashLayer::*;
    assert_ne!(stereo[&CanonicalSmiles], no_stereo[&CanonicalSmiles]);
    assert_eq!(
        no_stereo[&CanonicalSmiles],
        hydrogen_as_is[&CanonicalSmiles]
    );
    for layer in [NoStereoSmiles, NoStereoTautomerHash, TautomerHash] {
        assert_eq!(stereo[&layer], no_stereo[&layer]);
        assert_eq!(no_stereo[&layer], hydrogen_as_is[&layer]);
    }
}

#[test]
fn empty() {
    let layers1 = get_mol_layers(&Mol::default()).unwrap();
    assert_eq!(layers1, layers(["", "", "", "", "_0_0", "[]", "_0_0"]));
    assert_eq!(
        get_mol_hash(&layers1),
        "efaeaf876d44bef651a85875fb869046cdd612c9"
    );
}

#[test]
fn shared_7984() {
    let smi = "[2H]C([2H])=C([2H])C(=O)O";
    let no_stereo_mol = mol_from_smiles(smi).unwrap();
    assert_eq!(no_stereo_mol.get_num_atoms(), 8);
    let mut stereo_mol = no_stereo_mol.clone();

    // A double bond with these features is equal to a STEREONONE one.
    let mut bond = stereo_mol.get_bond_with_idx_mut(2).unwrap();
    bond.set_stereo(BondStereo::Stereoany);
    bond.set_bond_dir(BondDir::Eitherdouble);

    assert_eq!(
        get_mol_layers(&no_stereo_mol).unwrap(),
        get_mol_layers(&stereo_mol).unwrap()
    );
}

#[test]
fn non_matching_pseudoatom_label() {
    let pol_sdf = r#"
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Pol 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"#;
    let mod_sdf = r#"
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Mod 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"#;
    let carbon_sdf = r#"
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"#;

    let hashes: BTreeSet<String> = [pol_sdf, mod_sdf, carbon_sdf]
        .iter()
        .map(|sdf| get_mol_hash(&hash_sdf(sdf, &GetMolLayersParams::default())))
        .collect();
    assert_eq!(hashes.len(), 3);
}

#[test]
fn matching_pseudoatom_label() {
    let pol_sdf = r#"
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Pol 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"#;
    let pol_sdf_2 = r#"
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 Pol 0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"#;

    let hashes: BTreeSet<String> = [pol_sdf, pol_sdf_2]
        .iter()
        .map(|sdf| get_mol_hash(&hash_sdf(sdf, &GetMolLayersParams::default())))
        .collect();
    assert_eq!(hashes.len(), 1);
}

fn get_hash(m: &Mol, scheme: HashScheme) -> String {
    let hash = get_mol_hash_with_hash_scheme(&get_mol_layers(m).unwrap(), scheme);
    assert_eq!(hash.len(), 40); // SHA-1 hashes are 40 hex digits long
    hash
}

fn same_hash(m1: &Mol, m2: &Mol, scheme: HashScheme) -> bool {
    get_hash(m1, scheme) == get_hash(m2, scheme)
}

#[test]
fn hash_schemes() {
    use HashScheme::*;
    let mol1 =
        mol_from_smiles("C[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,11,o2:5,7|").unwrap();
    assert_eq!(mol1.get_num_atoms(), 14);

    let all_layers_hash = get_hash(&mol1, AllLayers);
    let stereo_insensitive_layers_hash = get_hash(&mol1, StereoInsensitiveLayers);
    let tautomer_insensitive_layers_hash = get_hash(&mol1, TautomerInsensitiveLayers);
    assert_ne!(all_layers_hash, tautomer_insensitive_layers_hash);
    assert_ne!(all_layers_hash, stereo_insensitive_layers_hash);
    assert_ne!(
        stereo_insensitive_layers_hash,
        tautomer_insensitive_layers_hash
    );

    // Compare stereoisomers
    let mol2 = mol_from_smiles("CC(O)C(C)CC1CC(C)CC(C)C1").unwrap();
    assert!(!same_hash(&mol1, &mol2, AllLayers));
    assert!(same_hash(&mol1, &mol2, StereoInsensitiveLayers));
    assert!(!same_hash(&mol1, &mol2, TautomerInsensitiveLayers));

    // Compare tautomers
    let mol1 = mol_from_smiles("N1C=NC2=CN=CN=C12").unwrap();
    let mol2 = mol_from_smiles("N1C=NC2=NC=NC=C12").unwrap();
    assert!(!same_hash(&mol1, &mol2, AllLayers));
    assert!(!same_hash(&mol1, &mol2, StereoInsensitiveLayers));
    assert!(same_hash(&mol1, &mol2, TautomerInsensitiveLayers));

    // Every scheme sees the S-group data (DAT/SRU/COP) the hash
    // canonicalizes.
    let smiles = "CNCC(=O)OC";
    let mol1 = mol_from_smiles(smiles).unwrap();
    let mol2 = mol_from_smiles(&format!("{smiles} |Sg:n:5,3,4::ht|")).unwrap();
    for scheme in [
        AllLayers,
        StereoInsensitiveLayers,
        TautomerInsensitiveLayers,
    ] {
        assert!(!same_hash(&mol1, &mol2, scheme));
    }
}

#[test]
fn no_stereo_hash() {
    // None of the StereoInsensitiveLayers layers has a SMILES stereo marker.
    let mol = mol_from_smiles(r"C\C=C/[C@@H](C)Cl").unwrap();
    assert_eq!(mol.get_num_atoms(), 6);
    let layers1 = get_mol_layers(&mol).unwrap();
    for layer in HashScheme::StereoInsensitiveLayers.value() {
        for marker in ['@', '/', '\\'] {
            assert!(!layers1[layer].contains(marker));
        }
    }
}

#[test]
fn atom_map_numbers() {
    // Atom map numbers are stripped (SHARED-8007).
    let atom_map_hash_layers = smiles_layers("Oc1cc([*:1])ccn1");
    let hash_layers = smiles_layers("Oc1cc([*])ccn1");
    assert_eq!(
        atom_map_hash_layers[&HashLayer::CanonicalSmiles],
        hash_layers[&HashLayer::CanonicalSmiles]
    );
}

#[test]
fn compare_stereoisomers() {
    let stereo_layers =
        smiles_layers("C[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,11,o2:5,7|");
    let stereo_insensitive_layers = smiles_layers("CC(O)C(C)CC1CC(C)CC(C)C1");
    assert_eq!(
        get_mol_hash_with_hash_scheme(&stereo_layers, HashScheme::StereoInsensitiveLayers),
        get_mol_hash_with_hash_scheme(
            &stereo_insensitive_layers,
            HashScheme::StereoInsensitiveLayers
        )
    );
    assert_eq!(
        stereo_layers[&HashLayer::NoStereoSmiles],
        "CC1CC(C)CC(CC(C)C(C)O)C1"
    );
    assert_eq!(
        stereo_insensitive_layers[&HashLayer::NoStereoSmiles],
        "CC1CC(C)CC(CC(C)C(C)O)C1"
    );
}

#[test]
fn sgroup_data_hashes() {
    let structure1 = r#"
  Mrv2108 11252113552D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 30 34 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.78 -9.2456 0 0
M  V30 2 C 7.2925 -8.8469 0 0
M  V30 3 C 6.2033 -9.9358 0 0
M  V30 4 C 6.6021 -11.4235 0 0
M  V30 5 C 8.0896 -11.8222 0 0
M  V30 6 C 9.1786 -10.733 0 0
M  V30 7 C 10.6662 -11.1317 0 0
M  V30 8 C 11.7552 -10.0426 0 0
M  V30 9 C 11.3565 -8.5552 0 0
M  V30 10 C 9.8688 -8.1565 0 0
M  V30 11 C 13.2428 -10.441 0 0
M  V30 12 C 14.3317 -9.3519 0 0 CFG=2
M  V30 13 C 15.8524 -9.5927 0 0
M  V30 14 C 16.5517 -8.2207 0 0
M  V30 15 C 15.463 -7.1316 0 0
M  V30 16 C 14.0905 -7.8307 0 0
M  V30 17 C 16.5515 -10.9649 0 0
M  V30 18 C 8.4882 -13.3096 0 0
M  V30 19 C 5.5131 -12.5124 0 0
M  V30 20 C 5.9117 -14.0001 0 0
M  V30 21 C 4.8227 -15.0888 0 0
M  V30 22 C 3.3352 -14.6902 0 0
M  V30 23 C 2.9365 -13.2027 0 0
M  V30 24 C 4.0255 -12.1139 0 0
M  V30 25 C 3.6269 -10.6262 0 0
M  V30 26 Cl 4.7157 -9.5372 0 0
M  V30 27 C 5.2472 -16.5691 0 0
M  V30 28 C 6.7414 -16.9417 0 0
M  V30 29 C 7.8112 -15.8339 0 0
M  V30 30 C 7.3867 -14.3535 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 2 7 8
M  V30 9 1 8 9
M  V30 10 2 9 10
M  V30 11 1 10 1
M  V30 12 1 12 11 CFG=1
M  V30 13 1 12 13
M  V30 14 1 13 14
M  V30 15 1 14 15
M  V30 16 1 15 16
M  V30 17 1 16 12
M  V30 18 1 13 17
M  V30 19 1 5 18
M  V30 20 1 4 19
M  V30 21 1 19 20
M  V30 22 2 20 21
M  V30 23 1 21 22
M  V30 24 2 22 23
M  V30 25 1 23 24
M  V30 26 2 24 19
M  V30 27 1 24 25
M  V30 28 1 3 26
M  V30 29 1 21 27
M  V30 30 1 20 30
M  V30 31 2 29 30
M  V30 32 1 28 29
M  V30 33 2 27 28
M  V30 34 1 8 11
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 12) FIELDNAME=stereolabel -
M  V30 FIELDDISP="   14.8391  -10.2748    DA    ALL  1       5" FIELDDATA=R
M  V30 2 DAT 0 ATOMS=(2 4 19) FIELDNAME=STEREOLABEL -
M  V30 FIELDDISP="    6.4455  -12.6102    DA    ALL  1       5" FIELDDATA=R
M  V30 END SGROUP
M  V30 END CTAB
M  END
"#;
    let structure2 = structure1.replace("FIELDDATA=R", "FIELDDATA=S");
    let structure3 = structure1.replacen("FIELDDATA=R", "FIELDDATA=S", 1);

    // By default these are all the same:
    let layers1 = hash_sdf(structure1, &GetMolLayersParams::default());
    let layers2 = hash_sdf(&structure2, &GetMolLayersParams::default());
    let expected_layers1 = layers([
        "Cc1ccc2ccccc2c1-c1c(Cl)cc2ccc(C[C@H]3CCCC3C)cc2c1C",
        "",
        "C29H29Cl",
        "Cc1ccc2ccccc2c1-c1c(Cl)cc2ccc(CC3CCCC3C)cc2c1C",
        "C[C]1[CH][CH][C]2[CH][CH][CH][CH][C]2[C]1[C]1[C](Cl)[CH][C]2[CH][CH][C](CC3CCCC3C)[CH][C]2[C]1C_0_0",
        "[]",
        "C[C]1[CH][CH][C]2[CH][CH][CH][CH][C]2[C]1[C]1[C](Cl)[CH][C]2[CH][CH][C](C[C@H]3CCCC3C)[CH][C]2[C]1C_0_0",
    ]);
    assert_eq!(layers1, expected_layers1);
    assert_eq!(layers2, layers1);

    // With both stereo labels they all differ:
    let both = data_field_names(&["STEREOLABEL", "stereolabel"]);
    compare_layers_with_expected_diffs(
        &hash_sdf(structure1, &both),
        &layers1,
        &[(
            HashLayer::SgroupData,
            r#"[{"fieldName": "STEREOLABEL", "atom": [10, 11], "bonds": [], "value": "R"}, {"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "R"}]"#,
        )],
    );
    compare_layers_with_expected_diffs(
        &hash_sdf(&structure2, &both),
        &layers1,
        &[(
            HashLayer::SgroupData,
            r#"[{"fieldName": "STEREOLABEL", "atom": [10, 11], "bonds": [], "value": "S"}, {"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "S"}]"#,
        )],
    );
    compare_layers_with_expected_diffs(
        &hash_sdf(&structure3, &both),
        &layers1,
        &[(
            HashLayer::SgroupData,
            r#"[{"fieldName": "STEREOLABEL", "atom": [10, 11], "bonds": [], "value": "R"}, {"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "S"}]"#,
        )],
    );

    // The data field to use can be picked:
    compare_layers_with_expected_diffs(
        &hash_sdf(structure1, &data_field_names(&["stereolabel"])),
        &layers1,
        &[(
            HashLayer::SgroupData,
            r#"[{"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "R"}]"#,
        )],
    );
}

#[test]
fn escape_layer() {
    let mol = mol_from_smiles("c1ccccc1").unwrap();
    let layers1 = get_mol_layers(&mol).unwrap();
    let params = GetMolLayersParams {
        escape: Some("make it unique".to_owned()),
        ..Default::default()
    };
    let layers_with_escape = get_mol_layers_with_params(&mol, &params).unwrap();
    assert_eq!(layers_with_escape[&HashLayer::Escape], "make it unique");
    for (k, v) in &layers_with_escape {
        if *k != HashLayer::Escape {
            assert_eq!(*v, layers1[k]);
        }
    }
}

#[test]
fn shared_8516() {
    assert_eq!(
        unsanitized_layers("C[C@H](N)C(=O)O"),
        unsanitized_layers("[H][C@@](C)(N)C(=O)O")
    );
}

#[test]
fn mol_layers_distinguish_carbon_isotopes() {
    assert_ne!(smiles_layers("c1ccccc1"), smiles_layers("[13c]1ccccc1"));
}

#[test]
fn mol_layers_distinguish_non_tautomeric_hydrogen_isotopes() {
    assert_eq!(
        smiles_layers("[2H]c1ccccc1"),
        layers([
            "[2H]c1ccccc1",
            "",
            "C6H6",
            "[2H]c1ccccc1",
            "[2H][C]1[CH][CH][CH][CH][CH]1_0_0",
            "[]",
            "[2H][C]1[CH][CH][CH][CH][CH]1_0_0",
        ])
    );
}

#[test]
fn mol_layers_distinguish_tautomeric_hydrogen_isotopes() {
    assert_eq!(
        smiles_layers("[2H]N1C=NC2=CN=CN=C12"),
        layers([
            "[2H]n1cnc2cncnc21",
            "",
            "C5H4N4",
            "[2H]n1cnc2cncnc21",
            "[2H]N1[CH][N][C]2[CH][N][CH][N][C]21_0_0",
            "[]",
            "[2H]N1[CH][N][C]2[CH][N][CH][N][C]21_0_0",
        ])
    );
    assert_eq!(
        smiles_layers("[2H]N1C=NC2=NC=NC=C12"),
        layers([
            "[2H]n1cnc2ncncc21",
            "",
            "C5H4N4",
            "[2H]n1cnc2ncncc21",
            "[2H]N1[CH][N][C]2[N][CH][N][CH][C]21_0_0",
            "[]",
            "[2H]N1[CH][N][C]2[N][CH][N][CH][C]21_0_0",
        ])
    );
}

#[test]
fn iodine() {
    // Stereo group canonicalization keeps isotopes.
    let layers1 = smiles_layers("CC[C@H](C)[999C@H](C)O  |o1:2,4|");
    assert_eq!(
        layers1[&HashLayer::CanonicalSmiles],
        "CC[C@H](C)[999C@H](C)O |o1:2,4|"
    );
}

#[test]
fn bad_bond_dir() {
    let mol_block = r#"
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 12 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 5.192672 -7.044349 0.000000 0
M  V30 2 C 5.272538 -5.675621 0.000000 0
M  V30 3 C 6.415595 -5.133214 0.000000 0
M  V30 4 C 7.487292 -5.854749 0.000000 0
M  V30 5 C 6.572156 -3.950482 0.000000 0
M  V30 6 C 5.935220 -2.698719 0.000000 0
M  V30 7 C 4.531603 -2.325867 0.000000 0
M  V30 8 C 4.159340 -0.962544 0.000000 0
M  V30 9 C 2.686034 -0.657857 0.000000 0
M  V30 10 C 5.166201 0.053779 0.000000 0
M  V30 11 C 6.621450 -0.316643 0.000000 0
M  V30 12 C 6.973858 -1.673599 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 2 3 5
M  V30 5 1 6 5 CFG=3
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 2 8 9
M  V30 9 1 8 10
M  V30 10 1 11 10
M  V30 11 1 6 12
M  V30 12 1 11 12
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC1 ATOMS=(1 6)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$
"#;
    let mut mol = mol_from_mol_block(mol_block).unwrap();
    // Simulate reapplying the wedging from the original mol block.
    let idx = mol.get_bond_between_atoms(4, 5).unwrap().unwrap().get_idx();
    mol.get_bond_with_idx_mut(idx)
        .unwrap()
        .set_bond_dir(BondDir::Beginwedge);
    get_mol_layers(&mol).unwrap();
}

#[test]
fn tautomer_v2_hash() {
    use HashLayer::{NoStereoTautomerHash, TautomerHash};
    let enol = mol_from_smiles("CC=CO").unwrap();
    let keto = mol_from_smiles("CCC=O").unwrap();

    // Version 1 of the tautomer hash, the default:
    let enol_layers = get_mol_layers(&enol).unwrap();
    let keto_layers = get_mol_layers(&keto).unwrap();
    for layer in [TautomerHash, NoStereoTautomerHash] {
        assert_ne!(enol_layers[&layer], keto_layers[&layer]);
    }
    assert_ne!(
        get_mol_hash_with_hash_scheme(&enol_layers, HashScheme::TautomerInsensitiveLayers),
        get_mol_hash_with_hash_scheme(&keto_layers, HashScheme::TautomerInsensitiveLayers)
    );

    // Version 2:
    let v2 = GetMolLayersParams {
        enable_tautomer_hash_v2: true,
        ..Default::default()
    };
    let enol_layers = get_mol_layers_with_params(&enol, &v2).unwrap();
    let keto_layers = get_mol_layers_with_params(&keto, &v2).unwrap();
    for layer in [TautomerHash, NoStereoTautomerHash] {
        assert_eq!(enol_layers[&layer], keto_layers[&layer]);
    }
    assert_eq!(
        get_mol_hash_with_hash_scheme(&enol_layers, HashScheme::TautomerInsensitiveLayers),
        get_mol_hash_with_hash_scheme(&keto_layers, HashScheme::TautomerInsensitiveLayers)
    );
}

#[test]
fn atropisomer() {
    // From RDKit's Code/GraphMol/FileParsers/test_data/atropisomers.
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/RP-6306_atrop1.sdf");
    let mut mol = mol_from_mol_file(path).unwrap();
    assert_eq!(
        mol.get_bond_with_idx(3).unwrap().get_stereo(),
        BondStereo::Stereoatropcw
    );

    let split = |text: &str| {
        let (smiles, extensions) = text.split_once(' ').unwrap();
        (smiles.to_owned(), extensions.to_owned())
    };
    let layers_cw = get_mol_layers(&mol).unwrap();
    let (smiles1, smi_ext1) = split(&layers_cw[&HashLayer::CanonicalSmiles]);
    let (taut1, taut_ext1) = split(&layers_cw[&HashLayer::TautomerHash]);
    assert_eq!(smiles1, "Cc1cc2c(C(N)=O)c(N)n(-c3c(C)ccc(O)c3C)c2nc1C");
    assert_eq!(
        taut1,
        "C[C]1[CH][C]2[C]([C]([N])[O])[C]([N])N([C]3[C](C)[CH][CH][C]([O])[C]3C)[C]2[N][C]1C_5_0"
    );
    assert_eq!(smi_ext1, "|wD:10.9|");
    assert_eq!(smi_ext1, taut_ext1);

    // The other atropisomer hashes differently:
    mol.get_bond_with_idx_mut(3)
        .unwrap()
        .set_stereo(BondStereo::Stereoatropccw);
    let layers_ccw = get_mol_layers(&mol).unwrap();
    assert_ne!(layers_cw, layers_ccw);

    // The SMILES and tautomer hash are the same, but not the extensions.
    let (smiles2, smi_ext2) = split(&layers_ccw[&HashLayer::CanonicalSmiles]);
    let (taut2, taut_ext2) = split(&layers_ccw[&HashLayer::TautomerHash]);
    assert_eq!(smiles2, smiles1);
    assert_eq!(taut2, taut1);
    assert_eq!(smi_ext2, "|wD:10.20|"); // the same atom, wedged down on another bond
    assert_eq!(smi_ext2, taut_ext2);
}

#[test]
fn github_8654() {
    let m1 = mol_from_smiles(
        "O=C(N[C@H]1C[C@@H](C(=O)O)[C@@H]2C[C@H]12)C1CC(=O)N(Cc2ccccn2)C1 |&1:3,5,9,11|",
    )
    .unwrap();
    let m2 = mol_from_smiles(
        "N(C(=O)C1CN(CC2=CC=CC=N2)C(=O)C1)[C@@H]3[C@@]4([C@@](C4)([C@H](C(O)=O)C3)[H])[H]",
    )
    .unwrap();
    let no_extensions = GetMolLayersParams {
        cxflag: CXSmilesFields::CX_NONE,
        ..Default::default()
    };
    // RDKit 2024.03.6 canonicalizes m1's stereo group even when cxflag does
    // not write it, which flips its stereo, so the molecules' layers differ
    // (RDKit 2026.03's fix for #8654 makes them equal).
    let layers1 = get_mol_layers_with_params(&m1, &no_extensions).unwrap();
    let layers2 = get_mol_layers_with_params(&m2, &no_extensions).unwrap();
    assert_eq!(
        layers1,
        layers([
            "O=C(N[C@H]1C[C@@H](C(=O)O)[C@@H]2C[C@H]12)C1CC(=O)N(Cc2ccccn2)C1",
            "",
            "C18H21N3O4",
            "O=C(NC1CC(C(=O)O)C2CC12)C1CC(=O)N(Cc2ccccn2)C1",
            "[O][C]([N]C1CC([C]([O])[O])C2CC12)C1C[C]([O])N(C[C]2[CH][CH][CH][CH][N]2)C1_2_0",
            "[]",
            "[O][C]([N][C@H]1C[C@@H]([C]([O])[O])[C@@H]2C[C@H]12)C1C[C]([O])N(C[C]2[CH][CH][CH][CH][N]2)C1_2_0",
        ])
    );
    assert_eq!(
        layers2[&HashLayer::CanonicalSmiles],
        "O=C(N[C@H]1C[C@@H](C(=O)O)[C@@H]2C[C@@H]21)C1CC(=O)N(Cc2ccccn2)C1"
    );
    assert_eq!(
        layers2[&HashLayer::TautomerHash],
        "[O][C]([N][C@H]1C[C@@H]([C]([O])[O])[C@@H]2C[C@@H]21)C1C[C]([O])N(C[C]2[CH][CH][CH][CH][N]2)C1_2_0"
    );
    assert_ne!(layers1, layers2);
    // With the default cxflag, m1's CXSmiles keeps its stereo group.
    let layers1 = get_mol_layers(&m1).unwrap();
    assert_eq!(
        layers1[&HashLayer::CanonicalSmiles],
        "O=C(N[C@H]1C[C@@H](C(=O)O)[C@@H]2C[C@H]12)C1CC(=O)N(Cc2ccccn2)C1 |&1:3,5,9,11|"
    );
    assert_eq!(
        layers1[&HashLayer::TautomerHash],
        "[O][C]([N][C@H]1C[C@@H]([C]([O])[O])[C@@H]2C[C@H]12)C1C[C]([O])N(C[C]2[CH][CH][CH][CH][N]2)C1_2_0 |&1:3,5,9,11|"
    );
}

#[test]
fn hashes_match_python() {
    use HashScheme::*;
    for (smiles, expected) in [
        (
            "C[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,11,o2:5,7|",
            [
                "cbdf21717a5eb862ddd8ad6265f646127215369f",
                "0bfd5ad6ea33688175caf9b446f27e5c6faaff95",
                "df3cda0c47d6b9e598e42521cc11f0e96bf7e4bf",
            ],
        ),
        (
            "CNCC(=O)OC |Sg:n:5,3,4::ht|",
            [
                "dcf6f8b9cf7dc5d20e8ad7597a77c6130f251660",
                "348e675ed6b69c48e88e055a2c0f09a018665567",
                "21a0fa47bec25d98cdb8864afdfc7b5d2455a2bb",
            ],
        ),
        (
            "Oc1cc([*:1])ccn1",
            [
                "917e9523255e9971be97249ce05247c500b66984",
                "cbd7aaa661d2f9c78d07ef459be53bb6741fa9c9",
                "ac6a7f66204a269213841b00551421b31d25af16",
            ],
        ),
        (
            "N1C=NC2=CN=CN=C12",
            [
                "627e9c086def7fccf3b01d9776782a75989499cf",
                "a20c4920df3ead309dfae9900c8f52cbd84a3048",
                "e315fad0fb9370ee31bb9be34b81a51129e91b72",
            ],
        ),
    ] {
        let layers1 = smiles_layers(smiles);
        for (scheme, hash) in [
            AllLayers,
            StereoInsensitiveLayers,
            TautomerInsensitiveLayers,
        ]
        .into_iter()
        .zip(expected)
        {
            assert_eq!(
                get_mol_hash_with_hash_scheme(&layers1, scheme),
                hash,
                "{smiles}"
            );
        }
    }
}

#[test]
fn sgroups_sort_like_python() {
    let polymer = smiles_layers("CNCC(=O)OC |Sg:n:5,3,4::ht|");
    assert_eq!(
        polymer[&HashLayer::SgroupData],
        r#"[{"type": "SRU", "atoms": [3, 4, 5], "bonds": [[5, 6]], "index": 1, "connect": "HT", "label": "", "XBHEAD": [[5, 6]], "XBCORR": []}]"#
    );

    let two_srus = smiles_layers(
        "*CC(*)C(*)C* |$star_e;;;star_e;;;star_e;star_e$,Sg:n:1,2::ht:::,Sg:n:4,5::ht:::|",
    );
    assert_eq!(
        two_srus[&HashLayer::SgroupData],
        r#"[{"type": "SRU", "atoms": [2, 3], "bonds": [[1, 2], [2, 4]], "index": 1, "connect": "HT", "label": "", "XBHEAD": [[4, 2], [6, 4]], "XBCORR": []}, {"type": "SRU", "atoms": [4, 6], "bonds": [[2, 4], [4, 5], [6, 7]], "index": 2, "connect": "HT", "label": "", "XBHEAD": [[1, 0]], "XBCORR": [[[1, 0], [2, 3]]]}]"#
    );
    assert_eq!(
        get_mol_hash(&two_srus),
        "3097685216cef3ae9b027dcc4c33c88ff4767d81"
    );

    // A data S-group sorts first and keeps no index; the others are
    // renumbered by position, and PARENT follows its group. The data value
    // is escaped as json.dumps does.
    let structure = r#"
  test

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 4 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 3 O 2 0 0 0
M  V30 4 C 3 0 0 0
M  V30 5 N 4 0 0 0
M  V30 6 C 5 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(1 2) XBONDS=(2 1 2) CONNECT=HT LABEL=n
M  V30 2 SRU 0 ATOMS=(1 4) XBONDS=(2 3 4) PARENT=1 CONNECT=HT LABEL=m
M  V30 3 COP 0 ATOMS=(4 2 3 4 5) XBONDS=(2 1 5) LABEL=co
M  V30 4 DAT 0 ATOMS=(1 6) FIELDNAME=Atrop FIELDDATA="x ""y"" \ é"
M  V30 END SGROUP
M  V30 END CTAB
M  END
"#;
    let layers1 = hash_sdf(structure, &GetMolLayersParams::default());
    assert_eq!(
        layers1[&HashLayer::SgroupData],
        r#"[{"fieldName": "Atrop", "atom": [5], "bonds": [], "value": "x \"y\" \\ \u00e9"}, {"type": "COP", "atoms": [1, 2, 3, 4], "index": 2}, {"type": "SRU", "atoms": [1], "bonds": [[0, 1], [1, 2]], "index": 3, "connect": "HT", "label": "n"}, {"type": "SRU", "atoms": [3], "bonds": [[2, 3], [3, 4]], "index": 4, "connect": "HT", "label": "m", "PARENT": 3}]"#
    );
    assert_eq!(
        get_mol_hash(&layers1),
        "c0b91acdff8c3439fb2e04079ef3cac7bb2a6011"
    );
}
