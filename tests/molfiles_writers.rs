// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Molecule writers. Expected values come from RDKit's Python wrappers,
//! writing to `io.StringIO`.

mod common;

use std::io::Write;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

/// Methanol and benzene, each named and tagged, without conformers.
fn two() -> Vec<Mol> {
    [("CO", "methanol"), ("c1ccccc1", "benzene")]
        .into_iter()
        .map(|(smi, name)| {
            let mut mol = mol_from_smiles(smi).unwrap();
            mol.set_prop("_Name", name);
            mol.set_prop("tag", &name.to_uppercase());
            mol.set_int_prop("n", smi.len() as i32);
            mol
        })
        .collect()
}

/// Methanol with 2D coordinates, named and tagged.
fn methanol_2d() -> Mol {
    let mut mol = mol_from_smiles("CO").unwrap();
    mol.set_prop("_Name", "methanol");
    let mut mol = mol_from_mol_block(&mol_to_mol_block(&mol).unwrap()).unwrap();
    mol.set_prop("tag", "METHANOL");
    mol
}

fn ethane() -> Mol {
    let mut mol = mol_from_smiles("CC").unwrap();
    mol.set_prop("_Name", "ethane");
    mol
}

const METHANOL_BLOCK: &str = "methanol\n     RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n";

const BENZENE_BLOCK: &str = "benzene\n     RDKit          2D\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  2  0\n  6  1  1  0\nM  END\n";

fn text(bytes: Vec<u8>) -> String {
    String::from_utf8(bytes).unwrap()
}

#[test]
fn writers_are_send() {
    fn send<T: Send>() {}
    send::<SDWriter<'_>>();
    send::<SmilesWriter<'_>>();
    send::<TDTWriter<'_>>();
    send::<PDBWriter<'_>>();
    send::<MaeWriter<'_>>();
}

#[test]
fn sd_writer() {
    let mut out = Vec::new();
    let mut writer = SDWriter::from_file_obj(&mut out).unwrap();
    assert!(writer.get_kekulize());
    assert!(!writer.get_force_v3000());
    for mol in two() {
        writer.write(&mol).unwrap();
    }
    assert_eq!(writer.num_mols(), 2);
    writer.close().unwrap();
    // As Python's: writing or flushing after close is a RuntimeError, and
    // closing again does nothing.
    assert!(matches!(
        writer.write(&ethane()),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(writer.flush(), Err(Error::Runtime { .. })));
    writer.close().unwrap();
    assert_eq!(writer.num_mols(), 2);
    drop(writer);
    assert_eq!(
        text(out),
        format!(
            "{METHANOL_BLOCK}>  <tag>  (1) \nMETHANOL\n\n>  <n>  (1) \n2\n\n$$$$\n\
             {BENZENE_BLOCK}>  <tag>  (2) \nBENZENE\n\n>  <n>  (2) \n8\n\n$$$$\n"
        )
    );

    let mut out = Vec::new();
    let mut writer = SDWriter::from_file_obj(&mut out).unwrap();
    writer.set_props(&["tag"]);
    writer.set_kekulize(false);
    writer.set_force_v3000(true);
    writer.write(&two()[1]).unwrap();
    writer.close().unwrap();
    drop(writer);
    let out = text(out);
    assert!(out.starts_with(
        "benzene\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n"
    ));
    assert!(out.contains("M  V30 1 4 1 2\n"));
    assert!(out.ends_with("M  END\n>  <tag>  (1) \nBENZENE\n\n$$$$\n"));

    // A failed write leaves the writer usable.
    let mut out = Vec::new();
    let mut writer = SDWriter::from_file_obj(&mut out).unwrap();
    let bad = mol_from_smiles_with_params(
        "c1cccc1",
        &SmilesParserParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(matches!(writer.write(&bad), Err(Error::Kekulize { .. })));
    assert!(writer.write_with_conf_id(&methanol_2d(), 3).is_err());
    writer.write_with_conf_id(&ethane(), -1).unwrap();
    assert_eq!(writer.num_mols(), 1);
    drop(writer);
    assert!(text(out).starts_with("ethane\n     RDKit          2D\n"));
}

#[test]
fn sd_writer_files() {
    let path = common::temp_path("written.sdf");
    let mut writer = SDWriter::new(&path).unwrap();
    writer.write(&two()[0]).unwrap();
    writer.flush().unwrap();
    writer.close().unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    assert!(written.starts_with(METHANOL_BLOCK));
    // Dropping a writer closes it.
    let mut writer = SDWriter::new(&path).unwrap();
    writer.write(&ethane()).unwrap();
    drop(writer);
    assert!(
        std::fs::read_to_string(&path)
            .unwrap()
            .starts_with("ethane\n")
    );
    assert!(matches!(
        SDWriter::new("/nonexistent/dir/x.sdf"),
        Err(Error::Os { message }) if message == "Bad output file /nonexistent/dir/x.sdf"
    ));
    // "-" is standard output.
    SDWriter::new("-").unwrap().close().unwrap();
}

#[test]
fn sd_text() {
    let methanol = &two()[0];
    assert_eq!(
        SDWriter::get_text(methanol).unwrap(),
        format!("{METHANOL_BLOCK}>  <tag>  \nMETHANOL\n\n>  <n>  \n2\n\n$$$$\n")
    );
    let params = SDWriterGetTextParams {
        kekulize: false,
        molid: 4,
        ..Default::default()
    };
    assert_eq!(
        SDWriter::get_text_with_params(methanol, &params).unwrap(),
        format!("{METHANOL_BLOCK}>  <tag>  (5) \nMETHANOL\n\n>  <n>  (5) \n2\n\n$$$$\n")
    );
    assert_eq!(
        SDWriterGetTextParams::default(),
        SDWriterGetTextParams {
            conf_id: -1,
            kekulize: true,
            force_v3000: false,
            molid: -1,
        }
    );
}

#[test]
fn sd_molids_are_checked() {
    // RDKit numbers properties as molid + 1 in int.
    let methanol = &two()[0];
    let at = |molid| SDWriterGetTextParams {
        molid,
        ..Default::default()
    };
    assert!(matches!(
        SDWriter::get_text_with_params(methanol, &at(i32::MAX)),
        Err(Error::Value { .. })
    ));
    assert!(
        SDWriter::get_text_with_params(methanol, &at(i32::MAX - 1))
            .unwrap()
            .contains(">  <tag>  (2147483647) \n")
    );
}

#[test]
fn smiles_writer() {
    let mut out = Vec::new();
    let mut writer = SmilesWriter::from_file_obj(&mut out).unwrap();
    writer.set_props(&["tag", "missing"]);
    for mol in two() {
        writer.write(&mol).unwrap();
    }
    writer.close().unwrap();
    drop(writer);
    assert_eq!(
        text(out),
        "SMILES Name tag missing\nCO methanol METHANOL \nc1ccccc1 benzene BENZENE \n"
    );

    let params = SmilesWriterParams {
        delimiter: "\t".into(),
        name_header: String::new(),
        include_header: false,
        isomeric_smiles: false,
        kekule_smiles: true,
    };
    let mut out = Vec::new();
    let mut writer = SmilesWriter::from_file_obj_with_params(&mut out, &params).unwrap();
    for mol in two() {
        writer.write(&mol).unwrap();
    }
    drop(writer);
    assert_eq!(text(out), "CO\nC1=CC=CC=C1\n");

    // Properties set after a molecule are ignored, as RDKit logs.
    let mut out = Vec::new();
    let mut writer = SmilesWriter::from_file_obj(&mut out).unwrap();
    writer.write(&ethane()).unwrap();
    writer.set_props(&["x"]);
    writer.write(&ethane()).unwrap();
    drop(writer);
    assert_eq!(text(out), "SMILES Name \nCC ethane\nCC ethane\n");
    assert_eq!(
        SmilesWriterParams::default(),
        SmilesWriterParams {
            delimiter: " ".into(),
            name_header: "Name".into(),
            include_header: true,
            isomeric_smiles: true,
            kekule_smiles: false,
        }
    );

    let path = common::temp_path("written.smi");
    let mut writer = SmilesWriter::new_with_params(&path, &params).unwrap();
    writer.write(&ethane()).unwrap();
    writer.close().unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "CC\n");
}

#[test]
fn tdt_writer() {
    let mut out = Vec::new();
    let mut writer = TDTWriter::from_file_obj(&mut out).unwrap();
    assert_eq!(writer.get_num_digits(), 4);
    assert!(!writer.get_write_2d());
    assert!(writer.get_write_names());
    writer.set_props(&["tag"]);
    writer.write(&methanol_2d()).unwrap();
    writer.set_write_2d();
    writer.set_num_digits(3).unwrap();
    writer.write(&methanol_2d()).unwrap();
    writer.set_write_names_with_state(false);
    writer.write(&ethane()).unwrap();
    assert!(writer.get_write_2d() && !writer.get_write_names());
    writer.close().unwrap();
    drop(writer);
    assert_eq!(
        text(out),
        "$SMI<CO>\nNAME<methanol>\n3D<0,0,0,1.299,0.75,0;>\ntag<METHANOL>\n|\n\
         $SMI<CO>\nNAME<methanol>\n2D<0,0,1.3,0.75;>\ntag<METHANOL>\n|\n$SMI<CC>\n|\n"
    );
    // Closing without molecules writes nothing.
    let mut out = Vec::new();
    TDTWriter::from_file_obj(&mut out).unwrap().close().unwrap();
    assert!(out.is_empty());
}

#[test]
fn tdt_digits_are_checked() {
    // libstdc++ formats each coordinate in a stack buffer as long as its
    // digits.
    let mut writer = TDTWriter::from_file_obj(Vec::new()).unwrap();
    assert!(matches!(
        writer.set_num_digits(1001),
        Err(Error::Value { .. })
    ));
    assert_eq!(writer.get_num_digits(), 4);
    writer.set_num_digits(1000).unwrap();
    writer.write(&methanol_2d()).unwrap();
}

#[test]
fn pdb_writer() {
    let record = "COMPND    methanol\n\
                  HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \n\
                  HETATM    2  O1  UNL     1       1.299   0.750   0.000  1.00  0.00           O  \n\
                  CONECT    1    2\nEND\n";
    let mut out = Vec::new();
    let mut writer = PDBWriter::from_file_obj(&mut out).unwrap();
    writer.write(&methanol_2d()).unwrap();
    writer.write(&methanol_2d()).unwrap();
    assert_eq!(writer.num_mols(), 2);
    drop(writer);
    assert_eq!(text(out), format!("{record}{record}"));

    let mut out = Vec::new();
    let mut writer = PDBWriter::from_file_obj_with_flavor(&mut out, 1).unwrap();
    writer.write(&methanol_2d()).unwrap();
    writer.write(&methanol_2d()).unwrap();
    drop(writer);
    assert_eq!(
        text(out),
        format!("MODEL        1\n{record}ENDMDL\nMODEL        2\n{record}ENDMDL\n")
    );

    let path = common::temp_path("written.pdb");
    let mut writer = PDBWriter::new_with_flavor(&path, 1).unwrap();
    writer.write(&methanol_2d()).unwrap();
    writer.close().unwrap();
    assert!(
        std::fs::read_to_string(&path)
            .unwrap()
            .starts_with("MODEL        1\n")
    );
}

#[test]
fn mae_writer() {
    let block = "f_m_ct {\n  i_m_ct_stereo_status\n  s_m_title\n  s_rdkit_tag\n  :::\n  1\n  methanol\n  METHANOL\n  m_atom[2] {\n    # First column is Index #\n    r_m_x_coord\n    r_m_y_coord\n    r_m_z_coord\n    i_m_atomic_number\n    i_m_formal_charge\n    :::\n    1 0.000000 0.000000 0.000000 6 0\n    2 1.299038 0.750000 0.000000 8 0\n    :::\n  }\n  m_bond[1] {\n    # First column is Index #\n    i_m_from\n    i_m_order\n    i_m_to\n    :::\n    1 1 1 2\n    :::\n  }\n}\n\n";
    let mut out = Vec::new();
    let mut writer = MaeWriter::from_fileobj(&mut out).unwrap();
    writer.set_props(&["tag"]);
    writer.write(&two()[0]).unwrap();
    assert_eq!(writer.num_mols(), 1);
    writer.close().unwrap();
    drop(writer);
    assert_eq!(
        text(out),
        format!(" {{\n  s_m_m2io_version\n  :::\n  2.0.0\n}}\n\n{block}")
    );
    let all = MaeWriter::get_text(&two()[0]).unwrap();
    assert!(all.starts_with("f_m_ct {\n  i_m_ct_stereo_status\n  i_rdkit_n\n  s_m_title\n  s_rdkit_tag\n  :::\n  1\n  2\n"));
    let params = MaeWriterGetTextParams {
        props_list: vec!["n".into()],
        ..Default::default()
    };
    assert!(
        MaeWriter::get_text_with_params(&two()[0], &params)
            .unwrap()
            .starts_with("f_m_ct {\n  i_m_ct_stereo_status\n  i_rdkit_n\n  s_m_title\n  :::\n")
    );
    // RDKit 2024.03 logs that a molecule without atoms cannot be exported,
    // and writes nothing, as Python's GetText returns ''.
    assert_eq!(MaeWriter::get_text(&Mol::default()).unwrap(), "");
    // Closing without molecules writes nothing.
    let mut out = Vec::new();
    MaeWriter::from_fileobj(&mut out).unwrap().close().unwrap();
    assert!(out.is_empty());

    let path = common::temp_path("written.mae");
    let mut writer = MaeWriter::new(&path).unwrap();
    writer.set_props(&["tag"]);
    writer.write(&two()[0]).unwrap();
    writer.close().unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    assert!(written.ends_with(block), "{written:?}");
}

#[test]
fn mae_inputs_are_checked() {
    // RDKit's Maestro writer casts any monomer info to a PDB residue info
    // (Python raises MemoryError, or worse).
    let mut mol = mol_from_smiles("CO").unwrap();
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_monomer_info(&AtomMonomerInfo::new(AtomMonomerType::Other).unwrap());
    assert!(matches!(
        MaeWriter::get_text(&mol),
        Err(Error::Value { .. })
    ));
    let mut writer = MaeWriter::from_fileobj(Vec::new()).unwrap();
    assert!(matches!(writer.write(&mol), Err(Error::Value { .. })));

    // It matches property names with a regex that recurses once per character
    // (Python crashes at 40,000).
    for (key, ok) in [("x".repeat(1000), true), ("x".repeat(1001), false)] {
        let mut mol = mol_from_smiles("CO").unwrap();
        mol.get_atom_with_idx_mut(1).unwrap().set_prop(&key, "v");
        assert_eq!(MaeWriter::get_text(&mol).is_ok(), ok);
        assert_eq!(writer.write(&mol).is_ok(), ok);
    }
    // Names it does not write are not checked.
    let mut mol = mol_from_smiles("CO").unwrap();
    mol.set_prop(&"x".repeat(5000), "v");
    let params = MaeWriterGetTextParams {
        props_list: vec!["other".into()],
        ..Default::default()
    };
    assert!(MaeWriter::get_text_with_params(&mol, &params).is_ok());
    mol.set_prop_with_computed(&"y".repeat(5000), "v", true);
    mol.clear_prop(&"x".repeat(5000));
    assert!(MaeWriter::get_text(&mol).is_ok());
}

/// A writer that fails after `limit` bytes.
struct Failing {
    limit: usize,
    written: usize,
}

impl Write for Failing {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if self.written + buf.len() > self.limit {
            return Err(std::io::Error::other("disk full"));
        }
        self.written += buf.len();
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// A writer that panics.
struct Panicking;

impl Write for Panicking {
    fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
        panic!("no writing here");
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn writer_errors_are_os_errors() {
    let mut writer = SDWriter::from_file_obj(Failing {
        limit: 10,
        written: 0,
    })
    .unwrap();
    // Writes are buffered, so flush finds the error, and every later call
    // reports it.
    writer.write(&ethane()).unwrap();
    let err = writer.flush().unwrap_err();
    assert!(
        matches!(&err, Error::Os { message } if message.contains("disk full")),
        "{err:?}"
    );
    assert!(matches!(writer.write(&ethane()), Err(Error::Os { .. })));
    assert!(matches!(writer.close(), Err(Error::Os { .. })));

    let mut writer = SmilesWriter::from_file_obj(Panicking).unwrap();
    writer.write(&ethane()).unwrap();
    assert!(matches!(
        writer.close(),
        Err(Error::Os { message }) if message.contains("panicked")
    ));
    // Dropping a writer that failed reports nothing.
    let mut writer = TDTWriter::from_file_obj(Panicking).unwrap();
    writer.write(&ethane()).unwrap();
    drop(writer);
}

#[test]
fn writers_borrow_or_own_their_output() {
    // An owned Write value lives as long as the writer.
    let mut writer = SDWriter::from_file_obj(Vec::new()).unwrap();
    writer.write(&ethane()).unwrap();
    // A writer may move to another thread.
    std::thread::spawn(move || {
        writer.write(&ethane()).unwrap();
        writer.close().unwrap();
    })
    .join()
    .unwrap();
    let mut out = Vec::new();
    {
        let mut writer = SmilesWriter::from_file_obj(&mut out).unwrap();
        writer.write(&ethane()).unwrap();
    }
    assert_eq!(text(out), "SMILES Name \nCC ethane\n");
}
