// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdmolfiles' other formats, CXSMILES, fragments, and ranking, with the
//! inputs the wrappers refuse. Expected values come from RDKit's Python
//! wrappers.

mod common;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chem_draw::{
    CDXFormat, mol_to_chem_draw_block, mol_to_chem_draw_block_format,
};
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

fn smiles(s: &str) -> Mol {
    mol_from_smiles(s).unwrap()
}

/// Methanol with explicit hydrogens and made-up 3D coordinates.
fn methanol_3d() -> Mol {
    let mut mol = add_hs(&smiles("OC")).unwrap();
    let mut conf = Conformer::new(mol.get_num_atoms());
    for i in 0..mol.get_num_atoms() {
        let f = f64::from(i);
        conf.set_atom_position(i, Point3D::new(f, 0.5 * f, 0.1 * f * f))
            .unwrap();
    }
    mol.add_conformer(&conf).unwrap();
    mol
}

fn is_value<T: std::fmt::Debug>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

fn is_parse<T: std::fmt::Debug>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Parse { .. }))
}

#[test]
fn defaults_match_python() {
    let cx = MolToCXSmilesParams::default();
    assert_eq!(cx.flags, CXSmilesFields::CX_ALL);
    assert_eq!(cx.flags.bits(), 0x7fff_ffff);
    assert_eq!(
        cx.restore_bond_dirs,
        RestoreBondDirOption::RestoreBondDirOptionClear
    );
    let scsr = MolFromSCSRParams::default();
    assert!(scsr.include_leaving_groups && scsr.sanitize && scsr.remove_hs);
    assert_eq!(scsr.scsr_template_names, SCSRTemplateNames::AsEntered);
    assert_eq!(
        scsr.scsr_base_hbond_options,
        SCSRBaseHbondOptions::UseSapAll
    );
    let pdb = MolFromPDBBlockParams::default();
    assert!(pdb.sanitize && pdb.remove_hs && pdb.proximity_bonding);
    assert_eq!(pdb.flavor, 0);
    assert_eq!(MolToXYZBlockParams::default().precision, 6);
    assert_eq!(
        MolToTPLBlockParams::default().partial_charge_prop,
        "_GasteigerCharge"
    );
    assert_eq!(CreatePropertyListParams::default().line_size, 190);
    assert!(MolFragmentToSmartsParams::default().isomeric_smarts);
    assert!(!CanonicalRankAtomsParams::default().include_chiral_presence);
    assert_eq!(
        CXSmilesFields::CX_ALL_BUT_COORDS,
        CXSmilesFields::CX_ALL & !CXSmilesFields::CX_COORDS
    );
}

#[test]
fn cx_smiles_and_smarts() {
    let mol = smiles("C[C@H](F)Cl |&1:1|");
    assert_eq!(mol_to_cx_smiles(&mol).unwrap(), "C[C@H](F)Cl |&1:1|");
    let none = MolToCXSmilesParams {
        flags: CXSmilesFields::CX_NONE,
        ..Default::default()
    };
    assert_eq!(
        mol_to_cx_smiles_with_params(&mol, &none).unwrap(),
        "C[C@H](F)Cl"
    );
    let query = mol_from_smarts("[#6]-[#8] |$R;$|").unwrap();
    assert_eq!(mol_to_cx_smarts(&query).unwrap(), "[#6]-[#8] |$R;$|");
    assert_eq!(
        format!("{:?}", CXSmilesFields::CX_COORDS),
        "CXSmilesFields::CX_COORDS"
    );
}

#[test]
fn fragments() {
    let mol = smiles("OCCN");
    assert_eq!(mol_fragment_to_smiles(&mol, &[1, 2, 3]).unwrap(), "CCN");
    let symbols = MolFragmentToSmilesParams {
        atom_symbols: Some(["a", "b", "c", "d"].map(String::from).to_vec()),
        ..Default::default()
    };
    assert_eq!(
        mol_fragment_to_smiles_with_params(&mol, &[1, 2], &symbols).unwrap(),
        "bc"
    );
    assert_eq!(mol_fragment_to_cx_smiles(&mol, &[0, 1]).unwrap(), "CO");
    assert_eq!(mol_fragment_to_smarts(&mol, &[1, 2]).unwrap(), "[#6]-[#6]");
    assert_eq!(
        mol_fragment_to_cx_smarts(&mol, &[0, 1]).unwrap(),
        "[#8]-[#6]"
    );
    // RDKit refuses an empty bond list, as in Python.
    let no_bonds = MolFragmentToSmartsParams {
        bonds_to_use: Some(vec![]),
        ..Default::default()
    };
    assert!(matches!(
        mol_fragment_to_smarts_with_params(&mol, &[1, 2], &no_bonds),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn fragment_arguments_are_checked() {
    let mol = smiles("OCCN");
    assert!(is_value(mol_fragment_to_smiles(&mol, &[])));
    assert!(matches!(
        mol_fragment_to_smiles(&mol, &[4]),
        Err(Error::Index { index: 4, len: 4 })
    ));
    let bonds = MolFragmentToSmilesParams {
        bonds_to_use: Some(vec![3]),
        ..Default::default()
    };
    assert!(matches!(
        mol_fragment_to_cx_smiles_with_params(&mol, &[0], &bonds),
        Err(Error::Index { index: 3, len: 3 })
    ));
    let short = MolFragmentToSmilesParams {
        bond_symbols: Some(vec!["-".into()]),
        ..Default::default()
    };
    assert!(is_value(mol_fragment_to_smiles_with_params(
        &mol,
        &[0],
        &short
    )));
    assert!(matches!(
        mol_fragment_to_cx_smarts(&mol, &[9]),
        Err(Error::Index { .. })
    ));
    assert!(matches!(
        canonical_rank_atoms_in_fragment(&mol, &[0, 7]),
        Err(Error::Index { index: 7, len: 4 })
    ));
}

#[test]
fn random_smiles() {
    let mol = smiles("OCCN");
    let seeded = MolToRandomSmilesVectParams {
        random_seed: 42,
        ..Default::default()
    };
    let expected = ["C(CO)N", "C(N)CO", "C(CO)N"];
    assert_eq!(
        mol_to_random_smiles_vect_with_params(&mol, 3, &seeded).unwrap(),
        expected
    );
    // The generator is process-wide; the call holds it, so a seeded call
    // repeats even with other threads drawing from it.
    let threads: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(move || {
                let mol = mol_from_smiles("OCCN").unwrap();
                let random = SmilesWriteParams {
                    do_random: true,
                    canonical: false,
                    ..Default::default()
                };
                for _ in 0..50 {
                    mol_to_smiles_with_params(&mol, &random).unwrap();
                }
                let seeded = MolToRandomSmilesVectParams {
                    random_seed: 42,
                    ..Default::default()
                };
                mol_to_random_smiles_vect_with_params(&mol, 3, &seeded).unwrap()
            })
        })
        .collect();
    for thread in threads {
        assert_eq!(thread.join().unwrap(), expected);
    }
}

#[test]
fn canonical_ranks() {
    let ties = CanonicalRankAtomsParams {
        break_ties: false,
        ..Default::default()
    };
    assert_eq!(
        canonical_rank_atoms_with_params(&smiles("C1NCN1"), &ties).unwrap(),
        [0, 2, 0, 2]
    );
    assert_eq!(canonical_rank_atoms(&smiles("OCCN")).unwrap(), [1, 3, 2, 0]);
    let params = CanonicalRankAtomsInFragmentParams {
        break_ties: false,
        ..Default::default()
    };
    assert_eq!(
        canonical_rank_atoms_in_fragment_with_params(
            &smiles("C1NCN1.C1NCN1"),
            &[0, 1, 2, 3],
            &params
        )
        .unwrap(),
        [4, 6, 4, 6, -1, -1, -1, -1]
    );
    let unsanitized = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let raw = mol_from_smiles_with_params("C1CC1", &unsanitized).unwrap();
    assert_eq!(canonical_rank_atoms(&raw).unwrap(), [0, 1, 2]);
    let mut stereo = smiles("C[C@H](F)[C@@H](C)Cl |&1:1,3|");
    canonicalize_enhanced_stereo(&mut stereo).unwrap();
    assert_eq!(
        mol_to_cx_smiles(&stereo).unwrap(),
        "C[C@H](F)[C@@H](C)Cl |&1:1,3|"
    );
}

#[test]
fn property_lists() {
    let mut mol = smiles("OCCN");
    for i in 0..4 {
        mol.get_atom_with_idx_mut(i)
            .unwrap()
            .set_int_prop("n", i as i32);
    }
    create_atom_int_property_list(&mut mol, "n").unwrap();
    assert_eq!(mol.get_prop("atom.iprop.n").unwrap(), "0 1 2 3");
    mol.get_atom_with_idx_mut(0).unwrap().clear_prop("n");
    let wrapped = CreatePropertyListParams {
        missing_value_marker: "x".into(),
        line_size: 5,
    };
    create_atom_int_property_list_with_params(&mut mol, "n", &wrapped).unwrap();
    assert_eq!(mol.get_prop("atom.iprop.n").unwrap(), "[x]\nx 1\n2 3");
    for i in 0..3 {
        mol.get_bond_with_idx_mut(i)
            .unwrap()
            .set_double_prop("d", 1.5);
    }
    create_bond_double_property_list(&mut mol, "d").unwrap();
    assert_eq!(mol.get_prop("bond.dprop.d").unwrap(), "1.5 1.5 1.5");
    create_bond_string_property_list(&mut mol, "missing").unwrap();
    assert_eq!(mol.get_prop("bond.prop.missing").unwrap(), "n/a n/a n/a");
    create_atom_bool_property_list(&mut mol, "missing").unwrap();
    create_atom_double_property_list(&mut mol, "missing").unwrap();
    create_atom_string_property_list(&mut mol, "missing").unwrap();
    create_bond_int_property_list(&mut mol, "missing").unwrap();
    create_bond_bool_property_list(&mut mol, "missing").unwrap();
}

#[test]
fn property_list_line_size_is_checked() {
    // Without a marker RDKit drops a character from the empty first line when
    // the first value fills it.
    let mut mol = smiles("CC");
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_prop("s", "abcdefghij");
    let short = CreatePropertyListParams {
        line_size: 5,
        ..Default::default()
    };
    assert!(is_value(create_atom_string_property_list_with_params(
        &mut mol, "s", &short
    )));
    let marked = CreatePropertyListParams {
        missing_value_marker: "x".into(),
        line_size: 5,
    };
    create_atom_string_property_list_with_params(&mut mol, "s", &marked).unwrap();
}

#[test]
fn mol_block_versions() {
    let mol = smiles("CO");
    assert!(
        mol_to_v2k_mol_block(&mol)
            .unwrap()
            .ends_with("  1  2  1  0\nM  END\n")
    );
    assert!(mol_to_v3k_mol_block(&mol).unwrap().contains("V3000"));
    let path = common::temp_path("v3k.mol");
    mol_to_v3k_mol_file(&mol, &path).unwrap();
    assert!(std::fs::read_to_string(&path).unwrap().contains("V3000"));
    let read = mol_from_mol_file(&path).unwrap();
    assert_eq!(mol_to_smiles(&read).unwrap(), "CO");
    mol_to_v3k_mol_file_with_params(&mol, &path, &MolWriterParams::default()).unwrap();
    assert!(
        mol_to_v2k_mol_block_with_params(&mol, &MolWriterParams::default())
            .unwrap()
            .contains("V2000")
    );
    let precise = MolWriterParams {
        precision: 1001,
        ..Default::default()
    };
    assert!(is_value(mol_to_v3k_mol_block_with_params(&mol, &precise)));
}

/// A V3000 block for `C[C@H](F)Cl |&1:1|` whose stereo-group line has `pad`
/// extra spaces, written on one line or across continuation lines.
fn stereo_group_block(pad: usize, continued: bool) -> String {
    let block = mol_to_v3k_mol_block(&smiles("C[C@H](F)Cl |&1:1|")).unwrap();
    let line = format!("MDLV30/STERAC1 ATOMS=(1 2{})", " ".repeat(pad));
    let line = if continued {
        let chars: Vec<char> = line.chars().collect();
        chars
            .chunks(60)
            .map(|chunk| chunk.iter().collect::<String>())
            .collect::<Vec<_>>()
            .join("-\nM  V30 ")
    } else {
        line
    };
    block.replace("MDLV30/STERAC1 ATOMS=(1 2)", &line)
}

#[test]
fn long_stereo_group_lines_are_refused() {
    // RDKit's regex would recurse once per character of these lines.
    assert!(mol_from_mol_block(&stereo_group_block(3000, false)).is_ok());
    assert!(mol_from_mol_block(&stereo_group_block(3000, true)).is_ok());
    assert!(is_value(mol_from_mol_block(&stereo_group_block(
        5000, false
    ))));
    assert!(is_value(mol_from_mol_block(&stereo_group_block(
        100_000, true
    ))));
    let path = common::temp_path("stereo-group.mol");
    std::fs::write(&path, stereo_group_block(100_000, false)).unwrap();
    assert!(is_value(mol_from_mol_file(&path)));
    // Other long lines are fine.
    let title = format!("{}\n", "x".repeat(100_000));
    let block = mol_to_mol_block(&smiles("CO")).unwrap();
    assert!(mol_from_mol_block(&format!("{title}{}", block.split_once('\n').unwrap().1)).is_ok());
}

/// A V3000 block for `*CC` with an explicit hydrogen, whose `*` bond has the
/// properties `props`.
fn end_points_block(props: &str) -> String {
    format!(
        "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
         M  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\n\
         M  V30 1 C 0 0 0 0\nM  V30 2 C 1.5 0 0 0\nM  V30 3 H 2 1 0 0\n\
         M  V30 4 * 0 1 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\n\
         M  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 4 1 {props}\n\
         M  V30 END BOND\nM  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn bond_end_points_are_checked() {
    // Removing the hydrogen rewrites the ENDPTS value; Python's values.
    for (props, end_points) in [
        ("ENDPTS=(1 2) ATTACH=ANY", "(1 2)"),
        ("ENDPTS=( 2 2  3 ) ATTACH=ANY", "(1 2)"),
        ("endpts=(2 2 3)\tATTACH=ANY", "(1 2)"),
        ("ENDPTS=(1 4294967295) ATTACH=ANY", "(1 4294967294)"),
        ("ENDPTS=(2 2 -\nM  V30 3) ATTACH=ANY", "(1 2)"),
    ] {
        let mol = mol_from_mol_block(&end_points_block(props)).unwrap();
        assert_eq!(mol_to_smiles(&mol).unwrap(), "*CC");
        let bond = mol.get_bond_with_idx(1).unwrap();
        assert_eq!(bond.get_prop("_MolFileBondEndPts").unwrap(), end_points);
    }
    // RWMol converts these numbers to unsigned int unchecked, reads past an
    // empty list, or reads the first character of an empty value.
    for props in [
        "ENDPTS=(1 1e20) ATTACH=ANY",
        "ENDPTS=(1 -3) ATTACH=ANY",
        "ENDPTS=(5e9 1) ATTACH=ANY",
        "ENDPTS=(1 4294967296) ATTACH=ANY",
        "ENDPTS=() ATTACH=ANY",
        "ENDPTS=(  ) ATTACH=ANY",
        "ENDPTS= ATTACH=ANY",
        "ATTACH=ANY ENDPTS=",
        "ATTACH=ANY ENDPTS=\r",
        "ATTACH=ANY \"ENDPTS=\"",
        "endpts=(1 1e20)",
        "ENDPTS=(1 2)x) ATTACH=ANY",
        "ENDPTS=(1 2-\nM  V30 e9) ATTACH=ANY",
    ] {
        let block = end_points_block(props);
        assert!(is_value(mol_from_mol_block(&block)), "{props:?}");
        let unsanitized = MolFromMolBlockParams {
            sanitize: false,
            ..Default::default()
        };
        assert!(
            is_value(mol_from_mol_block_with_params(&block, &unsanitized)),
            "{props:?}"
        );
    }
}

#[test]
fn scsr_blocks() {
    // A hydrogen bond between atoms that are not template atoms, in a block
    // without templates: RDKit indexes its empty template list.
    let block = "\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
                 M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n\
                 M  V30 1 N 0.0 0.0 0.0 0\nM  V30 2 O 1.5 0.0 0.0 0\n\
                 M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 10 1 2\n\
                 M  V30 END BOND\nM  V30 END CTAB\nM  V30 BEGIN TEMPLATE\n\
                 M  V30 END TEMPLATE\nM  END\n";
    assert!(is_parse(mol_from_scsr_block(block)));
    let ignore = MolFromSCSRParams {
        scsr_base_hbond_options: SCSRBaseHbondOptions::Ignore,
        ..Default::default()
    };
    let mol = mol_from_scsr_block_with_params(block, &ignore).unwrap();
    assert_eq!(mol.get_num_atoms(), 2);
    let path = common::temp_path("hbond.scsr");
    std::fs::write(&path, block).unwrap();
    assert!(is_parse(mol_from_scsr_file(&path)));
    assert!(mol_from_scsr_file_with_params(&path, &ignore).is_ok());
    assert!(matches!(
        mol_from_scsr_file("does-not-exist.mol"),
        Err(Error::Os { .. })
    ));
}

#[test]
fn tpl_round_trip() {
    let mol = methanol_3d();
    let block = mol_to_tpl_block(&mol).unwrap();
    assert!(block.starts_with("BioCAD format, all rights reserved\n"));
    assert_eq!(
        mol_to_smiles(&mol_from_tpl_block(&block).unwrap()).unwrap(),
        "[H]OC([H])([H])[H]"
    );
    let path = common::temp_path("m.tpl");
    mol_to_tpl_file(&mol, &path).unwrap();
    let read = mol_from_tpl_file(&path).unwrap();
    assert_eq!(read.get_num_conformers(), 1);
    let twice = MolToTPLFileParams {
        write_first_conf_twice: true,
        ..Default::default()
    };
    mol_to_tpl_file_with_params(&mol, &path, &twice).unwrap();
    let skip = MolFromTPLFileParams {
        skip_first_conf: true,
        ..Default::default()
    };
    assert_eq!(
        mol_from_tpl_file_with_params(&path, &skip)
            .unwrap()
            .get_num_conformers(),
        1
    );
    assert!(matches!(
        mol_from_tpl_file("does-not-exist.tpl"),
        Err(Error::Os { .. })
    ));
    assert_eq!(mol_to_tpl_block(&smiles("CO")).unwrap(), "");
}

#[test]
fn tpl_short_lines_are_refused() {
    // RDKit reads the second field of the counts and CONFS lines unchecked.
    let block = mol_to_tpl_block(&methanol_3d()).unwrap();
    assert!(is_parse(mol_from_tpl_block(
        &block.replace("\n6 5\n", "\n6\n")
    )));
    assert!(is_parse(mol_from_tpl_block(
        &block.replace("CONFS 0", "CONFS")
    )));
    // A count the block cannot hold is refused before RDKit allocates for it.
    assert!(is_parse(mol_from_tpl_block(
        &block.replace("\n6 5\n", "\n4000000000 5\n")
    )));
}

#[test]
fn marvin_cml_xyz_pdb() {
    let phenol = smiles("c1ccccc1O");
    let mrv = mol_to_mrv_block(&phenol).unwrap();
    assert!(mrv.contains("chemaxon"));
    assert_eq!(
        mol_to_smiles(&mol_from_mrv_block(&mrv).unwrap()).unwrap(),
        "Oc1ccccc1"
    );
    let path = common::temp_path("m.mrv");
    mol_to_mrv_file(&phenol, &path).unwrap();
    assert_eq!(
        mol_to_smiles(&mol_from_mrv_file(&path).unwrap()).unwrap(),
        "Oc1ccccc1"
    );
    assert!(
        mol_to_cml_block(&smiles("CO"))
            .unwrap()
            .contains("<molecule id=\"m-1\"")
    );
    mol_to_cml_file(&smiles("CO"), common::temp_path("m.cml")).unwrap();
    // RDKit 2025.09 reads four neighbors for a tetrahedral atom's parity.
    let chiral = smiles("N[C@@H](C)C(=O)O");
    assert!(is_value(mol_to_cml_block(&chiral)));
    assert!(is_value(mol_to_cml_file(
        &chiral,
        common::temp_path("c.cml")
    )));
    assert!(
        mol_to_cml_block(&add_hs(&chiral).unwrap())
            .unwrap()
            .contains("atomParity")
    );

    let mol = methanol_3d();
    let xyz = mol_to_xyz_block_with_params(
        &mol,
        &MolToXYZBlockParams {
            precision: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(xyz.starts_with("6\n\nO      0.000    0.000    0.000\n"));
    assert_eq!(
        mol_to_smiles(&mol_from_xyz_block(&xyz).unwrap()).unwrap(),
        "C.O.[HH].[HH].[HH].[HH]"
    );
    let path = common::temp_path("m.xyz");
    mol_to_xyz_file(&mol, &path).unwrap();
    assert_eq!(mol_from_xyz_file(&path).unwrap().get_num_atoms(), 6);
    assert_eq!(mol_from_xyz_block("").unwrap().get_num_atoms(), 0);
    let precise = MolToXYZBlockParams {
        precision: 1001,
        ..Default::default()
    };
    assert!(is_value(mol_to_xyz_block_with_params(&mol, &precise)));

    assert_eq!(
        mol_to_pdb_block(&smiles("CO")).unwrap(),
        "HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \n\
         HETATM    2  O1  UNL     1       0.000   0.000   0.000  1.00  0.00           O  \n\
         CONECT    1    2\nEND\n"
    );
    let pdb = mol_to_pdb_block(&mol).unwrap();
    assert_eq!(
        mol_to_smiles(&mol_from_pdb_block(&pdb).unwrap()).unwrap(),
        "CO"
    );
    let path = common::temp_path("m.pdb");
    mol_to_pdb_file(&mol, &path).unwrap();
    // Proximity bonding gives the carbon a fifth bond, as in Python.
    let keep_hs = MolFromPDBFileParams {
        remove_hs: false,
        ..Default::default()
    };
    assert!(matches!(
        mol_from_pdb_file_with_params(&path, &keep_hs),
        Err(Error::AtomValence { atom_idx: 1, .. })
    ));
    let no_proximity = MolFromPDBFileParams {
        remove_hs: false,
        proximity_bonding: false,
        ..Default::default()
    };
    assert_eq!(
        mol_from_pdb_file_with_params(&path, &no_proximity)
            .unwrap()
            .get_num_atoms(),
        6
    );
}

#[test]
fn mrv_bonds_need_two_atoms() {
    // RDKit copies the second atom reference before it counts them.
    let mrv = mol_to_mrv_block(&smiles("CO")).unwrap();
    let bad = mrv.replace("atomRefs2=\"a1 a2\"", "atomRefs2=\"a1\"");
    assert!(is_parse(mol_from_mrv_block(&bad)));
    let path = common::temp_path("bad.mrv");
    std::fs::write(&path, bad).unwrap();
    assert!(is_parse(mol_from_mrv_file(&path)));
}

#[test]
fn pdb_proximity_coordinates_are_checked() {
    // RDKit hashes coordinates in int for proximity bonding; its PDB fields
    // hold digits, signs, and points, up to 1e8.
    let line = "HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  ";
    let at = |x: &str| format!("{}{x}{}\nEND\n", &line[..30], &line[38..]);
    assert_eq!(
        mol_from_pdb_block(&at("-9999999")).unwrap().get_num_atoms(),
        1
    );
    for x in ["99999999", "20000000", "10000001"] {
        let block = format!("{}{x}{}\nEND\n", &line[..30], &line[38..]);
        assert!(is_value(mol_from_pdb_block(&block)), "{x}");
        let no_proximity = MolFromPDBBlockParams {
            proximity_bonding: false,
            ..Default::default()
        };
        assert_eq!(
            mol_from_pdb_block_with_params(&block, &no_proximity)
                .unwrap()
                .get_num_atoms(),
            1
        );
    }
}

#[test]
fn residue_info_writers() {
    // RDKit's PDB and sequence writers cast monomer infos of type PDBRESIDUE
    // to AtomPDBResidueInfo, which only such infos can have (U10); others are
    // written without residues.
    assert!(matches!(
        AtomMonomerInfo::new(AtomMonomerType::Pdbresidue),
        Err(Error::Value { .. })
    ));
    let mut other = smiles("CO");
    other
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_monomer_info(&AtomMonomerInfo::new(AtomMonomerType::Other).unwrap());
    assert!(mol_to_pdb_block(&other).is_ok());
    assert!(mol_to_pdb_file(&other, common::temp_path("other.pdb")).is_ok());
}

#[test]
fn sequences() {
    let peptide = mol_from_sequence("GAV").unwrap();
    assert_eq!(
        mol_to_smiles(&peptide).unwrap(),
        "CC(C)[C@H](NC(=O)[C@H](C)NC(=O)CN)C(=O)O"
    );
    assert_eq!(mol_to_sequence(&peptide).unwrap(), "GAV");
    assert_eq!(mol_to_fasta(&peptide).unwrap(), ">\nGAV\n");
    assert_eq!(mol_to_helm(&peptide).unwrap(), "PEPTIDE1{G.A.V}$$$$");
    let fasta = mol_from_fasta(">pep\nGAV").unwrap();
    assert_eq!(fasta.get_prop("_Name").unwrap(), "pep");
    assert_eq!(
        mol_to_sequence(&mol_from_helm("PEPTIDE1{G.A.V}$$$$").unwrap()).unwrap(),
        "GAV"
    );
    let dna = MolFromSequenceParams {
        flavor: 6,
        ..Default::default()
    };
    assert_eq!(
        mol_from_sequence_with_params("ACG", &dna)
            .unwrap()
            .get_num_atoms(),
        59
    );
    assert!(is_parse(mol_from_sequence("AA!")));
    assert!(is_parse(mol_from_helm_with_sanitize("nonsense", false)));
    // RDKit's FASTAToMol names the null result of a titled sequence it cannot
    // parse.
    assert!(is_parse(mol_from_fasta(">t\nGA!")));
}

#[test]
fn chem_draw() {
    assert!(has_chem_draw_cdx_support());
    let cdxml = mol_to_chem_draw_block(&smiles("CO")).unwrap();
    assert!(cdxml.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>"));
    let cdx = mol_to_chem_draw_block_format(&smiles("CO"), CDXFormat::Cdx).unwrap();
    assert!(cdx.starts_with(b"VjCD0100"));
    // RDKit scales 2D coordinates by the mean bond length, which a molecule
    // without bonds lacks, and casts them to 32-bit fixed point.
    assert!(is_value(mol_to_chem_draw_block(&smiles("C"))));
    let mut far = methanol_3d();
    far.get_conformer_mut()
        .unwrap()
        .set_atom_position(0, Point3D::new(1e6, 0.0, 0.0))
        .unwrap();
    assert!(is_value(mol_to_chem_draw_block(&far)));
    assert!(mol_to_chem_draw_block(&methanol_3d()).is_ok());
}

#[test]
fn smiles_debug_parse_must_be_zero() {
    // RDKit keeps the parser's debug level in a global every parse reads.
    let debug = SmilesParserParams {
        debug_parse: 1,
        ..Default::default()
    };
    assert!(is_value(mol_from_smiles_with_params("C", &debug)));
    let debug = SmartsParserParams {
        debug_parse: 1,
        ..Default::default()
    };
    assert!(is_value(mol_from_smarts_with_params("C", &debug)));
}

#[test]
fn depicting_writers_run_on_many_threads() {
    // Each writer depicts a molecule without a conformer under the depiction
    // lock (U4); the results match a single thread's.
    let mol = smiles("C[C@H](F)C1=CC=CC=C1");
    let expected = (
        mol_to_v2k_mol_block(&mol).unwrap(),
        mol_to_v3k_mol_block(&mol).unwrap(),
        mol_to_mrv_block(&mol).unwrap(),
        mol_to_chem_draw_block(&mol).unwrap(),
    );
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for _ in 0..20 {
                    let mol = smiles("C[C@H](F)C1=CC=CC=C1");
                    assert_eq!(mol_to_v2k_mol_block(&mol).unwrap(), expected.0);
                    assert_eq!(mol_to_v3k_mol_block(&mol).unwrap(), expected.1);
                    assert_eq!(mol_to_mrv_block(&mol).unwrap(), expected.2);
                    assert_eq!(mol_to_chem_draw_block(&mol).unwrap(), expected.3);
                }
            });
        }
    });
}

const MOL2_ETHANOL: &str = "@<TRIPOS>MOLECULE
ethanol
 3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3     1  MOL       0.0000
      2 C2          1.5000    0.0000    0.0000 C.3     1  MOL       0.0000
      3 O1          2.0000    1.4000    0.0000 O.3     1  MOL       0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     2     END    1
";

#[test]
fn mol2_bond_atoms_in_range() {
    // RDKit 2025.09.6's Python gives [C][C][O-] for this ethanol.
    let mol = mol_from_mol2_block(&MOL2_ETHANOL.replace("END", "3")).unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "[C][C][O-]");
    // RDKit 2025.09 reads past its atom table for a bond with one atom out
    // of range (and 0 wraps to past it); the wrapper refuses both.
    for atom in ["4", "0", "50000000"] {
        let block = MOL2_ETHANOL.replace("END", atom);
        assert!(is_value(mol_from_mol2_block(&block)), "{atom}");
        let path = std::env::temp_dir().join(format!("rdkit-cxx-mol2-{atom}.mol2"));
        std::fs::write(&path, &block).unwrap();
        assert!(is_value(mol_from_mol2_file(&path)), "{atom}");
        std::fs::remove_file(&path).unwrap();
    }
    // Without a BOND header RDKit 2025.09 reads the bond lines from the start
    // of the data (valgrind: invalid reads, for this block's first line).
    let headless = format!(
        "# 1 99 1\n{}",
        MOL2_ETHANOL
            .replace("END", "3")
            .replace("@<TRIPOS>BOND\n", "")
    );
    assert!(is_value(mol_from_mol2_block(&headless)));
}

/// A V3000 block for `CO` with a superatom SGroup over the O whose line ends
/// with `fields`.
fn sgroup_block(fields: &str) -> String {
    format!(
        "\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
         M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 1 0 0\nM  V30 BEGIN ATOM\n\
         M  V30 1 C 0 0 0 0\nM  V30 2 O 1 0 0 0\nM  V30 END ATOM\n\
         M  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\n\
         M  V30 BEGIN SGROUP\nM  V30 1 SUP 0 ATOMS=(1 2) XBONDS=(1 1) LABEL=OH {fields}\n\
         M  V30 END SGROUP\nM  V30 END CTAB\nM  END\n"
    )
}

fn attach_points(mol: &Mol) -> Vec<(u32, i32, String)> {
    get_mol_substance_groups(mol)[0]
        .get_attach_points()
        .unwrap()
        .into_iter()
        .map(|a| (a.a_idx, a.lv_idx, a.id))
        .collect()
}

const UNSANITIZED: MolFromMolBlockParams = MolFromMolBlockParams {
    sanitize: false,
    remove_hs: true,
    strict_parsing: true,
};

#[test]
fn sgroup_fields_are_checked() {
    // Python's values.
    for (fields, expected) in [
        ("SAP=(3 2 1 1)", (1, 0, "1")),
        ("SAP=(3 2 AIDX 1) CSTATE=(4 1 0 0 0)", (1, 1, "1")),
        ("SAP=( 3  2 1   1x)", (1, 0, "1x")),
    ] {
        let mol = mol_from_mol_block_with_params(&sgroup_block(fields), &UNSANITIZED).unwrap();
        let (a, lv, id) = expected;
        assert_eq!(attach_points(&mol), [(a, lv, id.to_string())], "{fields}");
    }
    // RDKit matches labels case-sensitively, so this is a string property.
    let mol = mol_from_mol_block_with_params(&sgroup_block("sap=(3 2 1)"), &UNSANITIZED).unwrap();
    assert!(attach_points(&mol).is_empty());
    // RDKit 2025.09 pops the last character of an empty SAP id, and reads
    // CSTATE numbers and an array's first value that were never read; its
    // Python returns None for these, the string's length having been
    // corrupted first.
    for fields in [
        "SAP=(3 2 1)",
        "SAP=(3 2",
        "SAP=(",
        "SAP=",
        "SAP",
        "SAP=(3x 2 1 1)",
        "SAP=(3 4294967296 1 1)",
        "SAP=(3 2 -\nM  V30 1)",
        "CSTATE=(4",
        "CSTATE=(4 ",
        "CSTATE",
        "XBCORR=(1",
        "XBCORR=(1 ",
        "COMPNO=",
        "COMPNO= ",
    ] {
        let block = sgroup_block(fields);
        assert!(
            is_value(mol_from_mol_block_with_params(&block, &UNSANITIZED)),
            "{fields:?}"
        );
        assert!(is_value(mol_from_mol_block(&block)), "{fields:?}");
    }
    // Values RDKit reads all of, or refuses before it reads the unset part:
    // the SAP id is "1"; a CSTATE count that does not read is 0, which RDKit
    // refuses for a superatom; an array whose count exceeds the bonds is
    // skipped.
    for fields in ["SAP=(3 2x 1 1)", "CSTATE=(x 1)", "XBCORR=(2", "COMPNO=3"] {
        let block = sgroup_block(fields);
        assert!(
            !is_value(mol_from_mol_block_with_params(&block, &UNSANITIZED)),
            "{fields:?}"
        );
    }
    // The end of the input is the end of a line.
    let block = sgroup_block("SAP=(3 2 1 1)");
    let cut = &block[..block.find("1 1)").unwrap() + 1];
    assert!(is_value(mol_from_mol_block_with_params(cut, &UNSANITIZED)));
}

/// An SD file whose second record's `field` starts `before` bytes before
/// the first 64 KiB boundary (the chunk a checked stream reads in), after a
/// first record padded with a data field.
fn straddling_sd_file(field: &str, before: usize) -> String {
    let second = format!("{}$$$$\n", sgroup_block(field));
    let head = format!("{}>  <pad>\n", sgroup_block("SAP=(3 2 1 1)"));
    let tail = "\n\n$$$$\n";
    let at = 1 << 16;
    let key = &field[..field.find('=').unwrap()];
    let offset = second.find(&format!(" {key}=")).unwrap();
    let pad = at - before - offset - head.len() - tail.len();
    format!("{head}{}{tail}{second}", "x".repeat(pad))
}

#[test]
fn sgroup_fields_across_chunks() {
    let dir = std::env::temp_dir();
    for before in [1, 4, 8, 11] {
        let path = dir.join(format!("rdkit-cxx-sap-{before}.sdf"));
        std::fs::write(&path, straddling_sd_file("SAP=(3 2 1 1)", before)).unwrap();
        let mut suppl = ForwardSDMolSupplier::new(&path).unwrap();
        let mols: Vec<Mol> = suppl.mols().map(|m| m.unwrap()).collect();
        assert_eq!(mols.len(), 2);
        assert_eq!(attach_points(&mols[1]), [(1, 0, "1".to_string())]);
        let mut suppl = SDMolSupplier::new(&path).unwrap();
        assert_eq!(
            suppl
                .mols()
                .collect::<rdkit_cxx::Result<Vec<_>>>()
                .unwrap()
                .len(),
            2
        );

        // Refused where the field starts, so RDKit never parses the line
        // cut off at the boundary.
        std::fs::write(&path, straddling_sd_file("SAP=(3 2 1)", before)).unwrap();
        let mut suppl = ForwardSDMolSupplier::new(&path).unwrap();
        let results: Vec<_> = suppl.mols().collect();
        assert!(results[0].is_ok());
        assert!(
            results
                .iter()
                .any(|r| matches!(r, Err(Error::Value { .. }))),
            "{before}"
        );
        // A COMPNO with no value, which RDKit reads unset.
        std::fs::write(&path, straddling_sd_file("COMPNO=", before)).unwrap();
        let mut suppl = ForwardSDMolSupplier::new(&path).unwrap();
        let results: Vec<_> = suppl.mols().collect();
        assert!(results[0].is_ok());
        assert!(
            results
                .iter()
                .any(|r| matches!(r, Err(Error::Value { .. }))),
            "{before}"
        );
        std::fs::remove_file(&path).unwrap();
    }
    // A file that ends in the field.
    let path = dir.join("rdkit-cxx-sap-end.sdf");
    let block = sgroup_block("SAP=(3 2 1 1)");
    std::fs::write(&path, &block[..block.find("1 1)").unwrap() + 1]).unwrap();
    let mut suppl = ForwardSDMolSupplier::new(&path).unwrap();
    assert!(suppl.mols().any(|r| matches!(r, Err(Error::Value { .. }))));
    std::fs::remove_file(&path).unwrap();
}

/// A reader that cannot seek and hands out at most 4093 bytes a read, as a
/// decompressing reader may.
struct Trickle(std::io::Cursor<Vec<u8>>);

impl std::io::Read for Trickle {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = buf.len().min(4093);
        self.0.read(&mut buf[..n])
    }
}

// The scanner holds back a CTAB line until it ends, reading the next chunk
// ahead into its own store, so a source that cannot seek back gives every
// byte once: a field straddling the 64 KiB boundary reads as from a file.
#[test]
fn sgroup_fields_across_chunks_of_unseekable_readers() {
    let read = |text: String| {
        let reader = Trickle(std::io::Cursor::new(text.into_bytes()));
        let mut suppl = ForwardSDMolSupplier::from_fileobj(reader).unwrap();
        suppl.mols().collect::<Vec<_>>()
    };
    for before in [1, 2, 4, 8, 11, 13] {
        let mols = read(straddling_sd_file("SAP=(3 2 1 1)", before));
        assert_eq!(mols.len(), 2, "{before}");
        let second = mols[1].as_ref().unwrap();
        assert_eq!(attach_points(second), [(1, 0, "1".to_string())], "{before}");
        let results = read(straddling_sd_file("COMPNO=", before));
        assert!(results[0].is_ok(), "{before}");
        assert!(
            results
                .iter()
                .any(|r| matches!(r, Err(Error::Value { .. }))),
            "{before}"
        );
    }
}

#[test]
fn sgroup_headers_are_checked() {
    let lax = MolFromMolBlockParams {
        sanitize: false,
        strict_parsing: false,
        ..Default::default()
    };
    let block = sgroup_block("");
    let header = "M  V30 1 SUP 0 ATOMS=(1 2) XBONDS=(1 1) LABEL=OH \n";
    assert!(block.contains(header));
    // RDKit 2025.09 reads the external index, or the type and external
    // index, of a header missing them without checking that it read them
    // (rdkit#9417 initializes them in 2026.03.5).
    for line in [
        "M  V30 1\n",
        "M  V30 1 SUP\n",
        "M  V30 1 SUP \n",
        "M  V30 \n",
    ] {
        let short = block.replace(header, line);
        assert!(
            is_value(mol_from_mol_block_with_params(&short, &lax)),
            "{line:?}"
        );
        // The SCSR reader always parses without strict parsing.
        let scsr = short.replace(
            "M  V30 END CTAB\n",
            "M  V30 END CTAB\nM  V30 BEGIN TEMPLATE\nM  V30 END TEMPLATE\n",
        );
        assert!(is_value(mol_from_scsr_block(&scsr)), "{line:?}");
    }
    // Python's value: RDKit 2025.09 reads the group without strict parsing.
    let mol = mol_from_mol_block_with_params(&block, &lax).unwrap();
    assert_eq!(get_mol_substance_groups(&mol).len(), 1);
    // A COUNTS line without SGroups: RDKit reads SGroup lines until a header
    // does not read, which the END SGROUP line ends.
    let uncounted = block.replace("COUNTS 2 1 1 0 0", "COUNTS 2 1 0 0 0");
    assert!(mol_from_mol_block(&uncounted).is_err());
    let mol = mol_from_mol_block_with_params(&uncounted, &lax).unwrap();
    assert_eq!(get_mol_substance_groups(&mol).len(), 1);
}

/// A V3000 block for `CC` with the enhanced-stereo collection line `line`.
fn collection_block(line: &str) -> String {
    format!(
        "\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
         M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n\
         M  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\n\
         M  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\n\
         M  V30 BEGIN COLLECTION\nM  V30 {line}\nM  V30 END COLLECTION\n\
         M  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn stereo_collections_are_checked() {
    // Python's value: the collection reads, and its group goes, as neither
    // atom is a stereocenter.
    let mol = mol_from_mol_block(&collection_block("MDLV30/STERAC1 ATOMS=(2 1 2)")).unwrap();
    assert_eq!(mol.get_stereo_groups().len(), 0);
    // RDKit 2025.09 reads the first atom number unset for a collection
    // shorter than its count, and keeps an atom listed twice, which dangles
    // once removed (valgrind: an invalid read).
    for line in [
        "MDLV30/STERAC1 ATOMS=(1 )",
        "MDLV30/STERAC1 ATOMS=(2 )",
        "MDLV30/STERAC1 ATOMS=(3 1 1 2)",
        "MDLV30/STEABS ATOMS=(2 1 1)",
    ] {
        assert!(
            is_value(mol_from_mol_block(&collection_block(line))),
            "{line}"
        );
    }
}
