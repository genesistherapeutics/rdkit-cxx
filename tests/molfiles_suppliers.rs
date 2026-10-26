// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Molecule suppliers. Expected values come from RDKit's Python wrappers.

mod common;

use std::io::{Cursor, Read};

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

/// Methanol, cyclopropane, and ethane as SD records, each named and tagged.
fn three_sd() -> String {
    let mut sd = String::new();
    for (smi, name) in [("CO", "methanol"), ("C1CC1", "ring"), ("CC", "ethane")] {
        let mut mol = mol_from_smiles(smi).unwrap();
        mol.set_prop("_Name", name);
        sd += &mol_to_mol_block(&mol).unwrap();
        sd += &format!("> <tag>\n{}\n\n$$$$\n", name.to_uppercase());
    }
    sd
}

fn names(mols: impl Iterator<Item = Result<Mol, Error>>) -> Vec<String> {
    mols.map(|mol| mol.unwrap().get_prop("_Name").unwrap())
        .collect()
}

fn write(name: &str, text: &str) -> std::path::PathBuf {
    let path = common::temp_path(name);
    std::fs::write(&path, text).unwrap();
    path
}

const SMILES: &str = "SMILES Name\nCO methanol\nC1CC1 ring\nCC ethane\n";

#[test]
fn suppliers_are_send() {
    fn send<T: Send>() {}
    send::<SDMolSupplier>();
    send::<ForwardSDMolSupplier>();
    send::<SmilesMolSupplier>();
    send::<TDTMolSupplier>();
    send::<MultithreadedSDMolSupplier>();
    send::<MultithreadedSmilesMolSupplier>();
}

#[test]
fn sd_supplier() {
    let path = write("three.sdf", &three_sd());
    let mut suppl = SDMolSupplier::new(&path).unwrap();
    assert_eq!(suppl.len().unwrap(), 3);
    assert_eq!(names(suppl.mols()), ["methanol", "ring", "ethane"]);
    assert!(suppl.at_end());
    assert!(
        suppl
            .get_item_text(1)
            .unwrap()
            .starts_with("ring\n     RDKit     ")
    );
    let ethane = suppl.get(2).unwrap();
    assert_eq!(mol_to_smiles(&ethane).unwrap(), "CC");
    assert_eq!(ethane.get_prop("tag").unwrap(), "ETHANE");
    assert!(matches!(
        suppl.get(3),
        Err(Error::Index { index: 3, len: 3 })
    ));
    assert!(matches!(
        suppl.get_item_text(3),
        Err(Error::Index { index: 3, len: 3 })
    ));
    // Reading continues after a random access, as Python's next does.
    suppl.get(0).unwrap();
    assert_eq!(names(suppl.mols_from_here()), ["ring", "ethane"]);
    // Many random accesses: RDKit's line count restarts with each.
    for _ in 0..1000 {
        suppl.get(1).unwrap();
    }

    let mut text = SDMolSupplier::default();
    assert_eq!(text.len().unwrap(), 0);
    // As Python's SetData(""), which yields one None.
    let mols: Vec<_> = text.mols().collect();
    assert!(matches!(mols[..], [Err(Error::Parse { .. })]));
    text.set_data(&three_sd()).unwrap();
    assert_eq!(names(text.mols()), ["methanol", "ring", "ethane"]);
    let unsanitized = SDMolSupplierParams {
        sanitize: false,
        ..Default::default()
    };
    text.set_data_bytes_with_params(three_sd().as_bytes(), &unsanitized)
        .unwrap();
    assert_eq!(text.len().unwrap(), 3);
    text.set_data_bytes(b"").unwrap();
    assert!(text.is_empty().unwrap());

    assert!(matches!(
        SDMolSupplier::new("does-not-exist.sdf"),
        Err(Error::Os { .. })
    ));
    let empty = write("empty.sdf", "");
    assert!(matches!(SDMolSupplier::new(&empty), Err(Error::Os { .. })));
}

#[test]
fn sd_property_lists() {
    let mut mol = mol_from_smiles("CO").unwrap();
    mol.set_prop("_Name", "m");
    let sd = mol_to_mol_block(&mol).unwrap() + "> <atom.prop.x>\na b\n\n$$$$\n";
    let mut suppl = SDMolSupplier::default();
    suppl.set_data(&sd).unwrap();
    assert!(suppl.get_process_property_lists());
    let read = suppl.get(0).unwrap();
    assert_eq!(
        read.get_atom_with_idx(1).unwrap().get_prop("x").unwrap(),
        "b"
    );
    suppl.set_process_property_lists(false);
    assert!(
        !suppl
            .get(0)
            .unwrap()
            .get_atom_with_idx(1)
            .unwrap()
            .has_prop("x")
    );
}

#[test]
fn bad_records() {
    // A record RDKit cannot read is an Err, where Python yields None.
    let sd = three_sd().replacen("  1  2  1  0\n", "  1  9  1  0\n", 1);
    let mut suppl = SDMolSupplier::default();
    suppl.set_data(&sd).unwrap();
    let mols: Vec<_> = suppl.mols().collect();
    assert_eq!(mols.len(), 3);
    assert!(matches!(mols[0], Err(Error::Parse { .. })));
    assert!(mols[1].is_ok() && mols[2].is_ok());
}

#[test]
fn forward_supplier() {
    let path = write("forward.sdf", &three_sd());
    let mut suppl = ForwardSDMolSupplier::new(&path).unwrap();
    assert_eq!(names(suppl.mols()), ["methanol", "ring", "ethane"]);
    assert!(suppl.at_end());
    assert!(suppl.get_eof_hit_on_read());
    let mut cursor = ForwardSDMolSupplier::from_fileobj(Cursor::new(three_sd())).unwrap();
    assert!(cursor.get_process_property_lists());
    cursor.set_process_property_lists(false);
    assert_eq!(names(cursor.mols()), ["methanol", "ring", "ethane"]);
    let params = ForwardSDMolSupplierParams {
        sanitize: false,
        ..Default::default()
    };
    let mut file = ForwardSDMolSupplier::from_fileobj_with_params(
        std::fs::File::open(&path).unwrap(),
        &params,
    )
    .unwrap();
    assert_eq!(file.mols().count(), 3);
    assert!(matches!(
        ForwardSDMolSupplier::new("does-not-exist.sdf"),
        Err(Error::Os { .. })
    ));
}

/// A reader that fails, or panics, after its first record.
struct Failing {
    text: Cursor<String>,
    panic: bool,
}

impl Read for Failing {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.text.position() > 0 {
            if self.panic {
                panic!("failing reader");
            }
            return Err(std::io::Error::other("failing reader"));
        }
        let n = buf.len().min(400);
        self.text.read(&mut buf[..n])
    }
}

#[test]
fn reader_errors_are_os_errors() {
    for panic in [false, true] {
        let reader = Failing {
            text: Cursor::new(three_sd()),
            panic,
        };
        let mut suppl = ForwardSDMolSupplier::from_fileobj(reader).unwrap();
        let mols: Vec<_> = suppl.mols().collect();
        assert!(
            mols.iter().any(|mol| matches!(mol, Err(Error::Os { .. }))),
            "{mols:?}"
        );
    }
}

#[test]
fn smiles_supplier() {
    let path = write("three.smi", SMILES);
    let mut suppl = SmilesMolSupplier::new(&path).unwrap();
    assert_eq!(suppl.len().unwrap(), 3);
    assert_eq!(names(suppl.mols()), ["methanol", "ring", "ethane"]);
    assert_eq!(suppl.get_item_text(1).unwrap(), "C1CC1 ring");
    assert_eq!(mol_to_smiles(&suppl.get(2).unwrap()).unwrap(), "CC");
    assert!(matches!(suppl.get(3), Err(Error::Index { index: 3, .. })));
    let mut text = smiles_mol_supplier_from_text(SMILES).unwrap();
    assert_eq!(names(text.mols()), ["methanol", "ring", "ethane"]);
    let tabs = SmilesMolSupplierParams {
        delimiter: "\t".into(),
        title_line: false,
        ..Default::default()
    };
    let mut tabbed =
        smiles_mol_supplier_from_text_with_params("CO\tmethanol\nCC\tethane\n", &tabs).unwrap();
    assert_eq!(names(tabbed.mols()), ["methanol", "ethane"]);
    // Python's next raises FileParseException for empty data, which ends
    // the iteration.
    let mut empty = SmilesMolSupplier::default();
    let mols: Vec<_> = empty.mols().collect();
    assert!(matches!(mols[..], [Err(Error::Parse { .. })]));
    empty.set_data(SMILES).unwrap();
    assert_eq!(empty.len().unwrap(), 3);
}

#[test]
fn smiles_name_column_is_checked() {
    // RDKit would read before each line's fields.
    let bad = SmilesMolSupplierParams {
        name_column: -2,
        ..Default::default()
    };
    assert!(matches!(
        smiles_mol_supplier_from_text_with_params(SMILES, &bad),
        Err(Error::Value { .. })
    ));
    let path = write("bad-column.smi", SMILES);
    assert!(matches!(
        SmilesMolSupplier::new_with_params(&path, &bad),
        Err(Error::Value { .. })
    ));
    let threaded = MultithreadedSmilesMolSupplierParams {
        name_column: -2,
        ..Default::default()
    };
    assert!(matches!(
        MultithreadedSmilesMolSupplier::new_with_params(&path, &threaded),
        Err(Error::Value { .. })
    ));
    let numbered = SmilesMolSupplierParams {
        name_column: -1,
        ..Default::default()
    };
    let mut suppl = smiles_mol_supplier_from_text_with_params(SMILES, &numbered).unwrap();
    assert_eq!(suppl.mols().count(), 3);
}

#[test]
fn tdt_supplier() {
    let tdt = "$SMI<CO>\nNAME<methanol>\n|\n$SMI<C1CC1>\nNAME<ring>\n|\n";
    let named = TDTMolSupplierParams {
        name_record: "NAME".into(),
        ..Default::default()
    };
    let mut suppl = TDTMolSupplier::default();
    suppl.set_data_with_params(tdt, &named).unwrap();
    assert_eq!(suppl.len().unwrap(), 2);
    assert_eq!(names(suppl.mols()), ["methanol", "ring"]);
    assert_eq!(
        suppl.get_item_text(1).unwrap(),
        "$SMI<C1CC1>\nNAME<ring>\n|\n"
    );
    assert_eq!(mol_to_smiles(&suppl.get(0).unwrap()).unwrap(), "CO");
    // RDKit indexes records in int; an index past the end is an IndexError.
    for idx in [2, u32::MAX] {
        assert!(matches!(suppl.get(idx), Err(Error::Index { len: 2, .. })));
        assert!(matches!(
            suppl.get_item_text(idx),
            Err(Error::Index { len: 2, .. })
        ));
    }
    let path = write("two.tdt", tdt);
    let mut file = TDTMolSupplier::new_with_params(&path, &named).unwrap();
    assert_eq!(names(file.mols()), ["methanol", "ring"]);
    assert_eq!(TDTMolSupplier::new(&path).unwrap().mols().count(), 2);
    // A SMILES RDKit cannot parse is an Err, and reading goes on.
    suppl
        .set_data("$SMI<C(>\nNAME<bad>\n|\n$SMI<CC>\n|\n")
        .unwrap();
    let mols: Vec<_> = suppl.mols().collect();
    assert_eq!(mols.len(), 2);
    assert!(matches!(mols[0], Err(Error::Parse { .. })));
    assert!(mols[1].is_ok());
}

#[test]
fn tdt_records() {
    // Expected values from Python, which reads the same records.
    let smiles = |mols: Vec<Result<Mol, Error>>| -> Vec<Option<String>> {
        mols.into_iter()
            .map(|mol| mol.ok().map(|mol| mol_to_smiles(&mol).unwrap()))
            .collect()
    };
    let mut suppl = TDTMolSupplier::default();
    // A blank line after a record ends the iteration, though not the count,
    // and text runs to the next record.
    suppl
        .set_data("$SMI<CO>\nPN<methanol>\n|\n\n$SMI<C1CC1>\nPN<ring>\n|\n")
        .unwrap();
    assert_eq!(smiles(suppl.mols().collect()), [Some("CO".into())]);
    assert_eq!(suppl.len().unwrap(), 2);
    assert_eq!(
        suppl.get_item_text(0).unwrap(),
        "$SMI<CO>\nPN<methanol>\n|\n\n"
    );
    assert_eq!(mol_to_smiles(&suppl.get(1).unwrap()).unwrap(), "C1CC1");
    assert!(suppl.at_end());
    // Lines between records belong to the one before; lines after the last
    // read as one more record without a molecule.
    suppl
        .set_data("$SMI<CO>\n|\njunk\n$SMI<CC>\n|\njunk\n")
        .unwrap();
    assert_eq!(suppl.get_item_text(0).unwrap(), "$SMI<CO>\n|\njunk\n");
    assert_eq!(suppl.get_item_text(1).unwrap(), "$SMI<CC>\n|\n");
    assert_eq!(
        smiles(suppl.mols().collect()),
        [Some("CO".into()), Some("CC".into()), None]
    );
    // A molecule RDKit cannot sanitize skips to its "|" line, after which a
    // blank line ends the iteration.
    suppl
        .set_data("$SMI<c1cc1>\nPN<x>\n|\n\n$SMI<CC>\n|\n")
        .unwrap();
    assert_eq!(smiles(suppl.mols().collect()), [None]);
    // Input without records reads as one record without a molecule.
    suppl.set_data("").unwrap();
    assert_eq!(suppl.len().unwrap(), 0);
    assert_eq!(smiles(suppl.mols().collect()), [None]);
    // An exception ends the iteration; RDKit does not move past the record.
    suppl.set_data("$SMI<CO>\nPN<x\n|\n$SMI<CC>\n|\n").unwrap();
    let mols: Vec<_> = suppl.mols().collect();
    assert!(matches!(mols[..], [Err(Error::Parse { .. })]));
    assert!(suppl.mols_from_here().next().unwrap().is_err());
}

#[test]
fn tdt_record_positions_are_checked() {
    // RDKit's supplier records the position after a last line without a
    // newline as -1, then counts it as a record whose text it reads past its
    // positions (Python returns garbage or raises MemoryError); on input
    // without records its len never returns.
    let mut suppl = TDTMolSupplier::default();
    suppl.set_data("$SMI<CO>\n|").unwrap();
    assert_eq!(suppl.len().unwrap(), 1);
    assert_eq!(suppl.get_item_text(0).unwrap(), "$SMI<CO>\n|");
    assert_eq!(suppl.len().unwrap(), 1);
    assert!(matches!(
        suppl.get_item_text(1),
        Err(Error::Index { index: 1, len: 1 })
    ));
    suppl.set_data("no records\n").unwrap();
    assert_eq!(suppl.len().unwrap(), 0);
    // Many reads of the last text: RDKit's would count lines without bound.
    suppl.set_data("$SMI<CO>\n|\n$SMI<CC>\nPN<x>\n").unwrap();
    for _ in 0..1000 {
        assert_eq!(suppl.get_item_text(1).unwrap(), "$SMI<CC>\nPN<x>\n");
    }
}

#[test]
fn files_that_change_are_refused() {
    // Bytes RDKit reads again must be the bytes it was checked against (U12).
    let sd = three_sd().repeat(500);
    assert!(sd.len() > 3 << 16);
    let path = write("changing.sdf", &sd);
    let mut suppl = SDMolSupplier::new(&path).unwrap();
    assert_eq!(suppl.len().unwrap(), 1500);
    assert_eq!(suppl.get(0).unwrap().get_prop("_Name").unwrap(), "methanol");
    std::fs::write(&path, sd.replace("methanol", "METHANOL")).unwrap();
    assert!(matches!(suppl.get(1499), Err(Error::Value { .. })));
    // The supplier stays failed.
    assert!(matches!(suppl.get(0), Err(Error::Value { .. })));
    assert!(matches!(suppl.get_item_text(0), Err(Error::Value { .. })));

    let path = write("shrinking.sdf", &sd);
    let mut suppl = SDMolSupplier::new(&path).unwrap();
    assert_eq!(suppl.len().unwrap(), 1500);
    std::fs::write(&path, &sd[..sd.len() / 2]).unwrap();
    let text = suppl.get_item_text(1499);
    assert!(matches!(text, Err(Error::Value { .. })), "{text:?}");

    // A file read once may change behind the reader.
    let path = write("growing.sdf", &sd);
    let mut forward = ForwardSDMolSupplier::new(&path).unwrap();
    assert_eq!(
        forward.mols().take(10).filter(|mol| mol.is_ok()).count(),
        10
    );
}

#[test]
fn multithreaded_sd_supplier() {
    let path = write("threaded.sdf", &three_sd());
    let mut suppl = MultithreadedSDMolSupplier::new(&path).unwrap();
    assert_eq!(suppl.get_last_record_id(), 0);
    let mut got = names(suppl.mols());
    got.sort();
    assert_eq!(got, ["ethane", "methanol", "ring"]);
    assert!(suppl.at_end());
    assert!((1..=3).contains(&suppl.get_last_record_id()));
    assert!(suppl.get_last_item_text().unwrap().contains("$$$$"));
    // RDKit's threads read the setting without a lock.
    assert!(matches!(
        suppl.set_process_property_lists(false),
        Err(Error::Value { .. })
    ));
    let params = MultithreadedSDMolSupplierParams {
        num_writer_threads: 3,
        size_input_queue: 1,
        size_output_queue: 1,
        ..Default::default()
    };
    let mut threaded = MultithreadedSDMolSupplier::new_with_params(&path, &params).unwrap();
    threaded.set_process_property_lists(false).unwrap();
    assert!(!threaded.get_process_property_lists());
    assert_eq!(threaded.mols().count(), 3);
    // Dropping a supplier mid-read stops its threads.
    let mut early = MultithreadedSDMolSupplier::new(&path).unwrap();
    early.mols().next().unwrap().unwrap();
    drop(early);
    let mut empty = MultithreadedSDMolSupplier::default();
    assert_eq!(empty.mols().count(), 0);
}

#[test]
fn multithreaded_parameters_are_checked() {
    let path = write("params.sdf", &three_sd());
    for params in [
        // RDKit negates the thread count as an int, and starts every
        // thread at the first read.
        MultithreadedSDMolSupplierParams {
            num_writer_threads: 1 << 31,
            ..Default::default()
        },
        MultithreadedSDMolSupplierParams {
            num_writer_threads: 1025,
            ..Default::default()
        },
        // A queue of size 0 never lets RDKit's threads finish.
        MultithreadedSDMolSupplierParams {
            size_input_queue: 0,
            ..Default::default()
        },
        MultithreadedSDMolSupplierParams {
            size_output_queue: 1 << 32,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            MultithreadedSDMolSupplier::new_with_params(&path, &params),
            Err(Error::Value { .. })
        ));
    }
}

#[test]
fn multithreaded_smiles_supplier() {
    let path = write("threaded.smi", SMILES);
    let mut suppl = MultithreadedSmilesMolSupplier::new(&path).unwrap();
    let mols: Vec<_> = suppl.mols().collect();
    // RDKit reads the end of the file as an empty last record, as Python's
    // trailing None.
    assert_eq!(mols.len(), 4);
    let mut got: Vec<String> = mols
        .iter()
        .filter_map(|mol| mol.as_ref().ok())
        .map(|mol| mol.get_prop("_Name").unwrap())
        .collect();
    got.sort();
    assert_eq!(got, ["ethane", "methanol", "ring"]);
    assert_eq!(suppl.get_last_record_id(), 4);
    assert_eq!(suppl.get_last_item_text().unwrap(), "");
    // A title line alone reads as one empty record, as Python's [None].
    let title = write("title.smi", "SMILES Name\n");
    let mut title = MultithreadedSmilesMolSupplier::new(&title).unwrap();
    let mols: Vec<_> = title.mols().collect();
    assert!(matches!(mols[..], [Err(Error::Parse { .. })]));
    let mut empty = MultithreadedSmilesMolSupplier::default();
    assert!(empty.mols().all(|mol| mol.is_err()));
}

#[test]
fn long_stereo_group_lines_are_refused() {
    // RDKit's mol-block regex would recurse once per character (U12).
    let block = mol_to_v3k_mol_block(&mol_from_smiles("C[C@H](F)Cl |&1:1|").unwrap()).unwrap();
    let long = block.replace(
        "MDLV30/STERAC1 ATOMS=(1 2)",
        &format!("MDLV30/STERAC1 ATOMS=(1 2{})", " ".repeat(100_000)),
    );
    let sd = format!("{block}$$$$\n{long}$$$$\n");
    let path = write("long.sdf", &sd);

    let mut suppl = SDMolSupplier::new(&path).unwrap();
    assert!(suppl.get(0).is_ok());
    assert!(matches!(suppl.get(1), Err(Error::Value { .. })));
    let mut suppl = SDMolSupplier::new(&path).unwrap();
    assert!(matches!(suppl.len(), Err(Error::Value { .. })));

    let mut forward = ForwardSDMolSupplier::from_fileobj(Cursor::new(sd.clone())).unwrap();
    let mols: Vec<_> = forward.mols().collect();
    assert!(mols[0].is_ok());
    assert!(matches!(mols[1], Err(Error::Value { .. })));

    let mut threaded = MultithreadedSDMolSupplier::new(&path).unwrap();
    let mols: Vec<_> = threaded.mols().collect();
    assert!(
        mols.iter()
            .any(|mol| matches!(mol, Err(Error::Value { .. })))
    );

    let mut text = SDMolSupplier::default();
    text.set_data(&sd).unwrap();
    assert!(matches!(text.len(), Err(Error::Value { .. })));
}
