// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn smiles(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn unsanitized(smiles: &str) -> Mol {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    mol_from_smiles_with_params(smiles, &params).unwrap()
}

fn bond_dirs(mol: &Mol) -> Vec<BondDir> {
    mol.get_bonds().map(|b| b.get_bond_dir()).collect()
}

#[test]
fn chiral_centers() {
    let mut mol = smiles("C[C@H](F)Cl");
    let params = AssignStereochemistryParams {
        clean_it: true,
        force: true,
        ..Default::default()
    };
    assign_stereochemistry_with_params(&mut mol, &params).unwrap();
    assert_eq!(
        mol.get_atom_with_idx(1)
            .unwrap()
            .get_prop("_CIPCode")
            .unwrap(),
        "R"
    );

    assert!(
        find_meso_centers(&smiles("C[C@H](F)[C@@H](F)C"))
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        find_meso_centers(&smiles("C[C@H](F)[C@H](F)C")).unwrap(),
        [(1, 3)]
    );

    let mut mol = mol_from_mol_block(CHLOROFLUOROETHANE_3D).unwrap();
    remove_stereochemistry(&mut mol).unwrap();
    let center = |mol: &Mol| mol.get_atom_with_idx(1).unwrap().get_chiral_tag();
    assert_eq!(center(&mol), ChiralType::ChiUnspecified);
    assign_atom_chiral_tags_from_structure(&mut mol).unwrap();
    assert_ne!(center(&mol), ChiralType::ChiUnspecified);
    remove_stereochemistry(&mut mol).unwrap();
    assign_stereochemistry_from_3d(&mut mol).unwrap();
    assert_eq!(
        mol.get_atom_with_idx(1)
            .unwrap()
            .get_prop("_CIPCode")
            .unwrap(),
        "R"
    );
    remove_non_explicit_3d_chirality(&mut mol).unwrap();
    remove_stereochemistry(&mut mol).unwrap();
    assign_atom_chiral_tags_from_mol_parity(&mut mol).unwrap();
    assert_ne!(center(&mol), ChiralType::ChiUnspecified);
    // A missing conformer is an error.
    let missing = AssignStereochemistryFrom3DParams {
        conf_id: 7,
        ..Default::default()
    };
    assert!(matches!(
        assign_stereochemistry_from_3d_with_params(&mut mol, &missing),
        Err(Error::Value { .. })
    ));
}

#[test]
fn double_bonds() {
    let mut mol = smiles("C/C=C/C");
    assert_eq!(
        mol.get_bond_with_idx(1).unwrap().get_stereo(),
        BondStereo::Stereoe
    );
    remove_stereochemistry(&mut mol).unwrap();
    assert_eq!(
        mol.get_bond_with_idx(1).unwrap().get_stereo(),
        BondStereo::Stereonone
    );
    find_potential_stereo_bonds(&mut mol).unwrap();
    assert_eq!(
        mol.get_bond_with_idx(1).unwrap().get_stereo(),
        BondStereo::Stereoany
    );

    let mut mol = mol_from_mol_block(&mol_to_mol_block(&smiles("C/C=C/C")).unwrap()).unwrap();
    detect_bond_stereochemistry(&mut mol).unwrap();
    assert_eq!(
        bond_dirs(&mol),
        [BondDir::Endupright, BondDir::None, BondDir::Endupright]
    );
    set_bond_stereo_from_directions(&mut mol).unwrap();
    assert_eq!(
        mol.get_bond_with_idx(1).unwrap().get_stereo(),
        BondStereo::Stereotrans
    );
    set_double_bond_neighbor_directions(&mut mol).unwrap();
    let conf = Conformer::from(mol.get_conformer().unwrap());
    set_double_bond_neighbor_directions_with_conf(&mut mol, &conf).unwrap();
    // RDKit reads a position for each atom, through its own range check.
    assert!(set_double_bond_neighbor_directions_with_conf(&mut mol, &Conformer::new(2)).is_err());

    let mut any = smiles("CC=CC");
    any.get_bond_with_idx_mut(1)
        .unwrap()
        .set_stereo(BondStereo::Stereoany);
    add_wavy_bonds_for_stereo_any(&mut any).unwrap();
    assert_eq!(
        any.get_bond_with_idx(0).unwrap().get_bond_dir(),
        BondDir::Unknown
    );
    assert_eq!(
        any.get_bond_with_idx(1).unwrap().get_stereo(),
        BondStereo::Stereonone
    );
}

#[test]
fn wedging() {
    let block = mol_to_mol_block(&smiles("C[C@H](F)Cl")).unwrap();
    let unwedged = || {
        let mut mol = mol_from_mol_block(&block).unwrap();
        for i in 0..mol.get_num_bonds() {
            mol.get_bond_with_idx_mut(i)
                .unwrap()
                .set_bond_dir(BondDir::None);
        }
        mol
    };
    let wedged = [BondDir::Beginwedge, BondDir::None, BondDir::None];

    let mut mol = unwedged();
    let conf = Conformer::from(mol.get_conformer().unwrap());
    wedge_mol_bonds(&mut mol, &conf).unwrap();
    assert_eq!(bond_dirs(&mol), wedged);
    assert!(wedge_mol_bonds(&mut mol, &Conformer::new(1)).is_err());

    let mut mol = unwedged();
    wedge_bond(&mut mol.get_bond_with_idx_mut(0).unwrap(), 1, -1).unwrap();
    assert_eq!(bond_dirs(&mol), wedged);
    assert!(matches!(
        wedge_bond(&mut mol.get_bond_with_idx_mut(0).unwrap(), 1, 5),
        Err(Error::Value { .. })
    ));
    reapply_mol_block_wedging(&mut mol).unwrap();
    assert_eq!(bond_dirs(&mol), wedged);
    assert!(!BondWedgingParameters::default().wedge_two_bonds_if_possible);
}

#[test]
fn enhanced_stereo() {
    let and = smiles("C[C@H](F)[C@H](F)C |&1:1,3|");
    assert_eq!(and.get_stereo_groups().len(), 1);

    let mut or = smiles("C[C@H](F)[C@@H](F)C |o1:1,3|");
    simplify_enhanced_stereo(&mut or).unwrap();
    assert!(or.get_stereo_groups().is_empty());

    let mut annotated = smiles("C[C@H](F)Cl");
    add_stereo_annotations(&mut annotated).unwrap();
    assert_eq!(
        annotated
            .get_atom_with_idx(1)
            .unwrap()
            .get_prop("atomNote")
            .unwrap(),
        "(R)"
    );

    let mut flagged =
        mol_from_mol_block(&mol_to_mol_block(&smiles("C[C@H](F)Cl")).unwrap()).unwrap();
    flagged.set_int_prop("_MolFileChiralFlag", 0);
    translate_chiral_flag_to_stereo_groups(&mut flagged).unwrap();
    let groups = flagged.get_stereo_groups();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].get_group_type(), StereoGroupType::StereoAnd);
}

#[test]
fn param_defaults() {
    let assign = AssignStereochemistryParams::default();
    assert!(!assign.clean_it && !assign.force && !assign.flag_possible_stereo_centers);
    let from_3d = AssignStereochemistryFrom3DParams::default();
    assert_eq!((from_3d.conf_id, from_3d.replace_existing_tags), (-1, true));
    let meso = FindMesoCentersParams::default();
    assert!(meso.include_isotopes && !meso.include_atom_maps);
    let labels = AddStereoAnnotationsParams::default();
    assert_eq!(labels.abs_label, "abs ({cip})");
    assert_eq!(labels.or_label, "or{id}");
    let wavy = AddWavyBondsForStereoAnyParams::default();
    assert!(wavy.clear_double_bond_flags);
    assert_eq!(
        wavy.add_when_impossible,
        StereoBondThresholds::DBL_BOND_NO_STEREO
    );
    assert_eq!(
        [
            StereoBondThresholds::DBL_BOND_NO_STEREO,
            StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO,
            StereoBondThresholds::CHIRAL_ATOM,
            StereoBondThresholds::DIRECTION_SET
        ],
        [1000, 10000, 100000, 1000000]
    );
    let expand = ExpandAttachmentPointsParams::default();
    assert!(expand.add_as_queries && expand.add_coords);
}

#[test]
fn attachment_points() {
    let mut mol = mol_from_mol_block(ATTACHMENT_POINT).unwrap();
    expand_attachment_points(&mut mol).unwrap();
    assert_eq!(mol.get_num_atoms(), 3);
    collapse_attachment_points(&mut mol).unwrap();
    assert_eq!(mol.get_num_atoms(), 2);
}

#[test]
fn chemistry_problems() {
    let problems = detect_chemistry_problems(&unsanitized("CC(C)(C)(C)(C)C")).unwrap();
    assert!(
        matches!(&problems[..], [Error::AtomValence { atom_idx: 1, .. }]),
        "{problems:?}"
    );
    let problems = detect_chemistry_problems(&unsanitized("c1cccc1")).unwrap();
    assert!(
        matches!(&problems[..], [Error::Kekulize { atom_indices, .. }] if atom_indices == &[0, 1, 2, 3, 4]),
        "{problems:?}"
    );
    assert!(
        detect_chemistry_problems(&unsanitized("CCO"))
            .unwrap()
            .is_empty()
    );
    let only_cleanup = detect_chemistry_problems_with_sanitize_ops(
        &unsanitized("c1cccc1"),
        SanitizeFlags::SANITIZE_CLEANUP,
    )
    .unwrap();
    assert!(only_cleanup.is_empty());
}

const CHLOROFLUOROETHANE_3D: &str = "
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.7970    0.0813   -0.0509 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6796   -0.1669    0.1335 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.0371   -1.2693   -0.6039 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5439    1.3066   -0.3680 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.8902    0.8090   -0.8849 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2793    0.4953    0.8547 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2604   -0.8913   -0.2635 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9662   -0.3646    1.1829 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  2  8  1  1
M  END
";

const ATTACHMENT_POINT: &str = "
  Mrv2014

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0 ATTCHPT=1
M  V30 2 O 1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";
