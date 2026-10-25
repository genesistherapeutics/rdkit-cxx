// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

#[test]
fn parser_and_writer_defaults() {
    let p = SmilesParserParams::default();
    assert!(p.allow_cxsmiles && p.parse_name && p.remove_hs && p.sanitize && p.strict_cxsmiles);
    assert_eq!((p.debug_parse, p.replacements.len()), (0, 0));
    let s = SmartsParserParams::default();
    assert!(s.allow_cxsmiles && s.parse_name && s.strict_cxsmiles && !s.merge_hs);
    let w = SmilesWriteParams::default();
    assert!(w.canonical && w.clean_stereo && w.do_isomeric_smiles && w.include_dative_bonds);
    assert!(
        !(w.all_bonds_explicit
            || w.all_hs_explicit
            || w.do_kekule
            || w.do_random
            || w.ignore_atom_map_numbers)
    );
    assert_eq!(w.rooted_at_atom, -1);
    let m = MolWriterParams::default();
    assert!(m.include_stereo && m.kekulize && !m.force_v3000);
    assert_eq!((m.precision, m.conf_id), (6, -1));
    let r = MolFromMolBlockParams::default();
    assert!(r.sanitize && r.remove_hs && r.strict_parsing);
    let f = MolFromMolFileParams::default();
    assert!(f.sanitize && f.remove_hs && f.strict_parsing);
}

#[test]
fn smiles_params() {
    let unsanitized = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let mol = mol_from_smiles_with_params("C(C)(C)(C)(C)C", &unsanitized).unwrap();
    assert_eq!(mol.get_num_atoms(), 6);
    let mut params = SmilesParserParams::default();
    params.replacements.insert("{X}".into(), "O".into());
    let mol = mol_from_smiles_with_params("C{X}", &params).unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "CO");
    let kekule = SmilesWriteParams {
        do_kekule: true,
        ..Default::default()
    };
    let benzene = mol_from_smiles("c1ccccc1").unwrap();
    assert_eq!(
        mol_to_smiles_with_params(&benzene, &kekule).unwrap(),
        "C1=CC=CC=C1"
    );
    let rooted = SmilesWriteParams {
        rooted_at_atom: 1,
        ..Default::default()
    };
    assert_eq!(
        mol_to_smiles_with_params(&mol_from_smiles("OCC").unwrap(), &rooted).unwrap(),
        "C(C)O"
    );
}

#[test]
fn random_smiles_threads() {
    // Random SMILES draw from RDKit's global random number generator, which
    // the wrapper locks.
    let random = SmilesWriteParams {
        do_random: true,
        ..Default::default()
    };
    std::thread::scope(|s| {
        for _ in 0..4 {
            s.spawn(|| {
                let mol = mol_from_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
                let canonical = mol_to_smiles(&mol).unwrap();
                let mut written = std::collections::BTreeSet::new();
                for _ in 0..50 {
                    let smiles = mol_to_smiles_with_params(&mol, &random).unwrap();
                    let parsed = mol_from_smiles(&smiles).unwrap();
                    assert_eq!(mol_to_smiles(&parsed).unwrap(), canonical);
                    written.insert(smiles);
                }
                assert!(written.len() > 1);
            });
        }
    });
}

#[test]
fn smarts() {
    assert_eq!(
        mol_to_smarts(&mol_from_smarts("[OX2H]").unwrap()).unwrap(),
        "[O&X2&H1]"
    );
    assert!(matches!(mol_from_smarts("[C"), Err(Error::Parse { .. })));
    let mut params = SmartsParserParams::default();
    params.replacements.insert("{N}".into(), "[#7]".into());
    assert_eq!(
        mol_to_smarts(&mol_from_smarts_with_params("C{N}", &params).unwrap()).unwrap(),
        "C[#7]"
    );
    let non_isomeric = SmilesWriteParams {
        do_isomeric_smiles: false,
        ..Default::default()
    };
    let chiral = mol_from_smiles("F[C@H](Cl)Br").unwrap();
    assert!(
        !mol_to_smarts_with_params(&chiral, &non_isomeric)
            .unwrap()
            .contains('@')
    );
}

#[test]
fn mol_blocks_and_files() {
    let mol = mol_from_smiles("CO").unwrap();
    let block = mol_to_mol_block(&mol).unwrap();
    assert_eq!(
        mol_to_smiles(&mol_from_mol_block(&block).unwrap()).unwrap(),
        "CO"
    );
    assert!(matches!(
        mol_from_mol_block("garbage"),
        Err(Error::Parse { .. })
    ));
    let v3000 = MolWriterParams {
        force_v3000: true,
        ..Default::default()
    };
    assert!(
        mol_to_mol_block_with_params(&mol, &v3000)
            .unwrap()
            .contains("V3000")
    );
    let path = common::temp_path("co.mol");
    mol_to_mol_file(&mol, &path).unwrap();
    assert_eq!(
        mol_to_smiles(&mol_from_mol_file(&path).unwrap()).unwrap(),
        "CO"
    );
    mol_to_mol_file_with_params(&mol, &path, &v3000).unwrap();
    let lax = MolFromMolFileParams {
        sanitize: false,
        ..Default::default()
    };
    assert_eq!(
        mol_from_mol_file_with_params(&path, &lax)
            .unwrap()
            .get_num_atoms(),
        2
    );
    assert!(matches!(mol_from_mol_file("does-not-exist.mol"),
                     Err(Error::Os { message }) if message.contains("does-not-exist.mol")));
}

#[test]
fn mol_block_params() {
    let with_hs = mol_to_mol_block(&add_explicit_h(&mol_from_smiles("C").unwrap())).unwrap();
    assert_eq!(mol_from_mol_block(&with_hs).unwrap().get_num_atoms(), 1);
    let keep = MolFromMolBlockParams {
        remove_hs: false,
        ..Default::default()
    };
    assert_eq!(
        mol_from_mol_block_with_params(&with_hs, &keep)
            .unwrap()
            .get_num_atoms(),
        2
    );
}

fn add_explicit_h(mol: &Mol) -> Mol {
    let mut rw = RWMol::from(mol.clone());
    let h = rw.add_atom(&Atom::from_num(1));
    rw.add_bond_with_order(0, h, BondType::Single).unwrap();
    rw.update_property_cache().unwrap();
    rw.into()
}

#[test]
fn latin1_title_is_a_value_error_not_an_abort() {
    let mol = mol_from_mol_file(common::latin1_mol_file()).unwrap();
    assert!(matches!(mol_to_mol_block(&mol), Err(Error::Value { .. })));
}

#[test]
fn mol_blocks_from_threads() {
    // Writing a molecule without a conformer computes 2D coordinates, whose
    // settings are globals the wrappers lock (U4).
    let expected = mol_to_mol_block(&mol_from_smiles("C[C@H](F)Cl").unwrap()).unwrap();
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for _ in 0..20 {
                    let mol = mol_from_smiles("C[C@H](F)Cl").unwrap();
                    assert_eq!(mol_to_mol_block(&mol).unwrap(), expected);
                }
            });
        }
    });
}
