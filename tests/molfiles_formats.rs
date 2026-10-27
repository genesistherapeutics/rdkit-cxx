// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdmolfiles' other formats, CXSMILES, fragments, and ranking, with the
//! inputs the wrappers refuse. Expected values come from RDKit's Python
//! wrappers.

mod common;

use rdkit_cxx::Error;
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
    create_atom_bool_property_list(&mut mol, "missing").unwrap();
    create_atom_double_property_list(&mut mol, "missing").unwrap();
    create_atom_string_property_list(&mut mol, "missing").unwrap();
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
    );
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for _ in 0..20 {
                    let mol = smiles("C[C@H](F)C1=CC=CC=C1");
                    assert_eq!(mol_to_v2k_mol_block(&mol).unwrap(), expected.0);
                    assert_eq!(mol_to_v3k_mol_block(&mol).unwrap(), expected.1);
                    assert_eq!(mol_to_mrv_block(&mol).unwrap(), expected.2);
                }
            });
        }
    });
}
