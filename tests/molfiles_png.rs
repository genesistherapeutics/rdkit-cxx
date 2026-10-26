// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Molecules and metadata in PNGs. Expected values come from RDKit's Python
//! wrappers.

mod common;

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

/// A 1x1 grayscale PNG, as Python's zlib and struct write it.
const PNG: &[u8] = &[
    137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73, 72, 68, 82, 0, 0, 0, 1, 0, 0, 0, 1, 8, 0, 0,
    0, 0, 58, 126, 155, 85, 0, 0, 0, 10, 73, 68, 65, 84, 120, 156, 99, 96, 0, 0, 0, 2, 0, 1, 72,
    175, 164, 113, 0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130,
];

fn smiles(s: &str) -> Mol {
    mol_from_smiles(s).unwrap()
}

fn keys(png: &[u8]) -> Vec<String> {
    metadata_from_png_string_as_list(png)
        .unwrap()
        .into_iter()
        .map(|(key, _)| key)
        .collect()
}

fn to_smiles(mols: Vec<Option<Mol>>) -> Vec<String> {
    mols.iter()
        .map(|mol| mol_to_smiles(mol.as_ref().unwrap()).unwrap())
        .collect()
}

/// `png` with a chunk of `kind` holding `data` before its IEND chunk.
fn with_chunk(png: &[u8], kind: &[u8; 4], data: &[u8]) -> Vec<u8> {
    let (body, end) = png.split_at(png.len() - 12);
    let mut out = body.to_vec();
    out.extend_from_slice(&(data.len() as u32).to_be_bytes());
    out.extend_from_slice(kind);
    out.extend_from_slice(data);
    out.extend_from_slice(&[0, 0, 0, 0]); // RDKit does not check the CRC
    out.extend_from_slice(end);
    out
}

#[test]
fn defaults_match_python() {
    let params = PNGMetadataParams::default();
    assert!(params.include_pkl && params.include_smiles && !params.include_mol);
    assert_eq!(params.property_flags, PropertyPickleOptions::NoProps);
    assert_eq!(params.cx_smiles_flags, CXSmilesFields::CX_ALL);
    assert_eq!(
        params.restore_bond_dirs,
        RestoreBondDirOption::RestoreBondDirOptionClear
    );
    assert_eq!(params.smiles_write_params, SmilesWriteParams::default());
    assert_eq!(MolsFromPNGStringParams::default().tag, "rdkitPKL");
    assert_eq!(MolsFromPNGFileParams::default().tag, "rdkitPKL");
}

#[test]
fn molecules() {
    assert!(metadata_from_png_string(PNG).unwrap().is_empty());
    let png = mol_metadata_to_png_string(&smiles("CO"), PNG).unwrap();
    assert_eq!(
        keys(&png),
        ["rdkitPKL rdkit 2026.09.1", "SMILES rdkit 2026.09.1"]
    );
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        mol_to_smiles(&unsafe { mol_from_png_string(&png) }.unwrap()).unwrap(),
        "CO"
    );
    let two = mol_metadata_to_png_string(&smiles("CCN"), &png).unwrap();
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        to_smiles(unsafe { mols_from_png_string(&two) }.unwrap()),
        ["CO", "CCN"]
    );
    let tagged = MolsFromPNGStringParams {
        tag: "SMILES".into(),
        ..Default::default()
    };
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        to_smiles(unsafe { mols_from_png_string_with_params(&two, &tagged) }.unwrap()),
        ["CO", "CCN"]
    );
    let mol_block = PNGMetadataParams {
        include_pkl: false,
        include_mol: true,
        ..Default::default()
    };
    let png = mol_metadata_to_png_string_with_params(&smiles("CO"), PNG, &mol_block).unwrap();
    assert_eq!(
        keys(&png),
        ["SMILES rdkit 2026.09.1", "MOL rdkit 2026.09.1"]
    );
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        mol_to_smiles(&unsafe { mol_from_png_string(&png) }.unwrap()).unwrap(),
        "CO"
    );

    let path = common::temp_path("two.png");
    std::fs::write(&path, &two).unwrap();
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        mol_to_smiles(&unsafe { mol_from_png_file(&path) }.unwrap()).unwrap(),
        "CO"
    );
    // SAFETY: a pickle RDKit wrote.
    assert_eq!(
        to_smiles(unsafe { mols_from_png_file(&path) }.unwrap()),
        ["CO", "CCN"]
    );
    let three = mol_metadata_to_png_file(&smiles("O"), &path).unwrap();
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        to_smiles(unsafe { mols_from_png_string(&three) }.unwrap()),
        ["CO", "CCN", "O"]
    );
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { mol_from_png_file("does-not-exist.png") },
        Err(Error::Os { .. })
    ));
}

#[test]
fn metadata() {
    let entries = BTreeMap::from([("a".into(), "b".into()), ("c".into(), "d".into())]);
    let png = add_metadata_to_png_string(&entries, PNG).unwrap();
    assert_eq!(
        metadata_from_png_string_as_list(&png).unwrap(),
        [("a".into(), b"b".to_vec()), ("c".into(), b"d".to_vec())]
    );
    assert_eq!(
        metadata_from_png_string(&png).unwrap(),
        BTreeMap::from([("a".into(), b"b".to_vec()), ("c".into(), b"d".to_vec())])
    );
    // No molecule: Python returns None and an empty tuple.
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { mol_from_png_string(&png) },
        Err(Error::Parse { .. })
    ));
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { mols_from_png_string(&png) }.unwrap().is_empty());
    let path = common::temp_path("metadata.png");
    std::fs::write(&path, &png).unwrap();
    assert_eq!(metadata_from_png_file(&path).unwrap().len(), 2);
    assert_eq!(metadata_from_png_file_as_list(&path).unwrap().len(), 2);
    let more =
        add_metadata_to_png_file(&BTreeMap::from([("e".into(), "f".into())]), &path).unwrap();
    assert_eq!(metadata_from_png_string(&more).unwrap().len(), 3);
}

#[test]
fn malformed_pngs_are_refused() {
    // RDKit walks the chunks without checking their lengths.
    let truncated = &PNG[..PNG.len() - 12];
    let no_keyword = with_chunk(PNG, b"tEXt", b"abc");
    // The IDAT chunk's length, at bytes 33 to 36, past the end of the data.
    let mut long = PNG.to_vec();
    long[33..37].copy_from_slice(&200u32.to_be_bytes());
    let mut header = PNG.to_vec();
    header[1] = b'Q';
    for png in [truncated, &no_keyword, &long, &header] {
        assert!(matches!(
            metadata_from_png_string(png),
            Err(Error::Parse { .. })
        ));
        // SAFETY: malformed PNG data, whose chunks the wrapper refuses before RDKit reads any
        // metadata.
        assert!(matches!(
            unsafe { mol_from_png_string(png) },
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            // SAFETY: as above.
            unsafe { mols_from_png_string(png) },
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            mol_metadata_to_png_string(&smiles("C"), png),
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            add_metadata_to_png_string(&BTreeMap::new(), png),
            Err(Error::Parse { .. })
        ));
    }
    let path = common::temp_path("truncated.png");
    std::fs::write(&path, truncated).unwrap();
    assert!(matches!(
        metadata_from_png_file(&path),
        Err(Error::Parse { .. })
    ));
    // A chunk with a keyword and an empty value is fine.
    let empty = with_chunk(PNG, b"tEXt", b"key\0");
    assert!(metadata_from_png_string(&empty).unwrap().is_empty());
}

#[test]
fn recorded_text_is_checked() {
    // A recorded mol block goes through the mol-block checks (U12).
    let block = mol_to_v3k_mol_block(&smiles("C[C@H](F)Cl |&1:1|")).unwrap();
    let long = block.replace(
        "MDLV30/STERAC1 ATOMS=(1 2)",
        &format!("MDLV30/STERAC1 ATOMS=(1 2{})", " ".repeat(10_000)),
    );
    let png = add_metadata_to_png_string(&BTreeMap::from([("MOL".into(), long)]), PNG).unwrap();
    assert!(matches!(
        // SAFETY: metadata with a mol block and no pickle.
        unsafe { mol_from_png_string(&png) },
        Err(Error::Value { .. })
    ));
    let png = add_metadata_to_png_string(&BTreeMap::from([("MOL".into(), block)]), PNG).unwrap();
    assert_eq!(
        // SAFETY: a pickle RDKit wrote.
        mol_to_smiles(&unsafe { mol_from_png_string(&png) }.unwrap()).unwrap(),
        "C[C@H](F)Cl"
    );
    let debug = SmilesParserParams {
        debug_parse: 1,
        ..Default::default()
    };
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { mol_from_png_string_with_params(&png, &debug) },
        Err(Error::Value { .. })
    ));
}

#[test]
fn depiction_refuses_what_the_depictor_would_mishandle() {
    // The mol block in the metadata depicts a molecule without a conformer,
    // and RDKit's depictor reads the first neighbor of a square-planar atom
    // unchecked (Python segfaults).
    let mut mol = smiles("[Pt]");
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiSquareplanar);
    let mol_block = PNGMetadataParams {
        include_mol: true,
        ..Default::default()
    };
    assert!(matches!(
        mol_metadata_to_png_string_with_params(&mol, PNG, &mol_block),
        Err(Error::Value { message }) if message.contains("square-planar")
    ));
    mol_metadata_to_png_string(&mol, PNG).unwrap();
}

#[test]
fn pickles_cannot_give_plain_info_the_pdb_type() {
    // A pickle recorded in the metadata goes through the monomer-type checks
    // that molecule pickles do (U10): RDKit casts any info of type PDBRESIDUE
    // to AtomPDBResidueInfo, reading past the end of a plain one.
    let mut mol = smiles("C");
    let params = AtomMonomerInfoParams {
        name: "MARK".into(),
        ..Default::default()
    };
    mol.get_atom_with_idx_mut(0).unwrap().set_monomer_info(
        &AtomMonomerInfo::new_with_params(AtomMonomerType::Other, &params).unwrap(),
    );
    let mut pkl = mol.to_binary();
    // The info's name (a u32 length, then its bytes), then its type as a byte.
    let at = pkl
        .windows(9)
        .position(|w| w == b"\x04\0\0\0MARK\x02")
        .unwrap()
        + 8;
    let chunk = |pkl: &[u8]| {
        let mut data = b"rdkitPKL\0".to_vec();
        data.extend_from_slice(pkl);
        with_chunk(PNG, b"tEXt", &data)
    };
    let png = chunk(&pkl);
    // SAFETY: a pickle RDKit wrote.
    let read = unsafe { mol_from_png_string(&png) }.unwrap();
    let info = read
        .get_atom_with_idx(0)
        .unwrap()
        .get_monomer_info()
        .unwrap();
    assert_eq!(info.get_monomer_type(), AtomMonomerType::Other);
    pkl[at] = 1;
    let png = chunk(&pkl);
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote but for a monomer type byte, which the wrapper checks (U10).
        unsafe { mol_from_png_string(&png) },
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        // SAFETY: as above.
        unsafe { mols_from_png_string(&png) },
        Err(Error::Value { .. })
    ));
}
