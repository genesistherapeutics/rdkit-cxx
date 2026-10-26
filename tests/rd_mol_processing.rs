// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! `rdMolProcessing`. Expected values come from RDKit's Python wrappers.

mod common;

use std::path::PathBuf;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_fingerprint_generator::{
    GetAtomPairGeneratorParams, GetMorganGeneratorParams, GetRDKitFPGeneratorParams,
    get_atom_pair_generator_with_params, get_morgan_feature_atom_inv_gen,
    get_morgan_generator_with_params, get_rdkit_fp_generator_with_params,
};
use rdkit_cxx::chem::rd_mol_processing::*;
use rdkit_cxx::chem::{mol_from_smiles, mol_to_mol_block, mol_to_v3k_mol_block};
use rdkit_cxx::data_structs::ExplicitBitVect;

const SMILES: [&str; 6] = [
    "CCO",
    "c1ccccc1",
    "C1CC1N",
    "CC(=O)O",
    "CCN(CC)CC",
    "O=C(O)c1ccccc1O",
];

/// Python's on bits for SMILES with the default generator (Morgan, radius 3).
const MORGAN3: [&[u32]; 6] = [
    &[80, 222, 294, 807, 1057, 1410],
    &[389, 1088, 1232, 1873],
    &[926, 1019, 1053, 1171, 1284, 1325, 1832],
    &[389, 508, 650, 807, 1017, 1057, 1917],
    &[80, 251, 294, 881, 932, 1057, 1561, 1981],
    &[
        266, 389, 456, 589, 650, 783, 807, 1047, 1088, 1199, 1211, 1380, 1447, 1517, 1519, 1602,
        1607, 1750, 1831, 1873, 1917, 1991,
    ],
];

fn write(name: &str, bytes: &[u8]) -> PathBuf {
    let path = common::temp_path(name);
    std::fs::write(&path, bytes).unwrap();
    path
}

fn smiles_file(smiles: &[&str]) -> String {
    let mut text = String::from("SMILES\tName\n");
    for (i, smi) in smiles.iter().enumerate() {
        text += &format!("{smi}\tm{i}\n");
    }
    text
}

fn sd_file(smiles: &[&str]) -> String {
    smiles
        .iter()
        .map(|smi| mol_to_mol_block(&mol_from_smiles(smi).unwrap()).unwrap() + "$$$$\n")
        .collect()
}

fn tdt_file(smiles: &[&str]) -> String {
    smiles
        .iter()
        .enumerate()
        .map(|(i, smi)| format!("$SMI<{smi}>\nNAME<m{i}>\n|\n"))
        .collect()
}

/// CRC-32 (IEEE), for gzip's trailer.
fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = if crc & 1 != 0 {
                (crc >> 1) ^ 0xEDB8_8320
            } else {
                crc >> 1
            };
        }
    }
    !crc
}

/// `data` as a gzip file of stored (uncompressed) deflate blocks.
fn gzip(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    let mut chunks = data.chunks(0xffff).peekable();
    if chunks.peek().is_none() {
        out.extend_from_slice(&[1, 0, 0, 0xff, 0xff]);
    }
    while let Some(chunk) = chunks.next() {
        out.push(u8::from(chunks.peek().is_none()));
        let len = u16::try_from(chunk.len()).unwrap();
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(&(!len).to_le_bytes());
        out.extend_from_slice(chunk);
    }
    out.extend_from_slice(&crc32(data).to_le_bytes());
    out.extend_from_slice(&u32::try_from(data.len()).unwrap().to_le_bytes());
    out
}

fn on_bits(fps: &[Option<ExplicitBitVect>]) -> Vec<Option<Vec<u32>>> {
    fps.iter()
        .map(|fp| fp.as_ref().map(|fp| fp.get_on_bits()))
        .collect()
}

fn counts(fps: &[Option<ExplicitBitVect>]) -> Vec<Option<u32>> {
    fps.iter()
        .map(|fp| fp.as_ref().map(|fp| fp.get_num_on_bits()))
        .collect()
}

fn morgan3() -> Vec<Option<Vec<u32>>> {
    MORGAN3.iter().map(|bits| Some(bits.to_vec())).collect()
}

fn with_threads(num_threads: i32) -> GetFingerprintsForMolsInFileParams<'static> {
    GetFingerprintsForMolsInFileParams {
        options: SupplierOptions {
            num_threads,
            ..Default::default()
        },
        ..Default::default()
    }
}

fn read(path: &PathBuf, num_threads: i32) -> Result<Vec<Option<ExplicitBitVect>>, Error> {
    get_fingerprints_for_mols_in_file_with_params(path, &with_threads(num_threads))
}

#[test]
fn defaults_match_python() {
    let options = SupplierOptions::default();
    assert_eq!(
        options,
        SupplierOptions {
            num_threads: 0,
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
            delimiter: "\t".to_string(),
            smiles_column: 0,
            name_column: 1,
            title_line: true,
            name_record: String::new(),
            conf_id_2d: -1,
            conf_id_3d: 0,
        }
    );
    let params = GetFingerprintsForMolsInFileParams::default();
    assert!(params.generator.is_none());
    assert_eq!(params.options, options);
}

#[test]
fn smiles_sd_and_tdt_files() {
    let smi = write("six.smi", smiles_file(&SMILES).as_bytes());
    assert_eq!(
        on_bits(&get_fingerprints_for_mols_in_file(&smi).unwrap()),
        morgan3()
    );
    let sdf = write("six.sdf", sd_file(&SMILES).as_bytes());
    let tdt = write("six.tdt", tdt_file(&SMILES).as_bytes());
    for num_threads in [1, 2, 4, 0, -1, i32::MIN, 8] {
        for path in [&smi, &sdf, &tdt] {
            assert_eq!(on_bits(&read(path, num_threads).unwrap()), morgan3());
        }
    }
    // A comma-separated file without a title line.
    let csv = write("two.csv", b"CCO,a\nCCN,b\n");
    for num_threads in [1, 2] {
        let params = GetFingerprintsForMolsInFileParams {
            options: SupplierOptions {
                num_threads,
                delimiter: ",".to_string(),
                title_line: false,
                ..Default::default()
            },
            ..Default::default()
        };
        let fps = get_fingerprints_for_mols_in_file_with_params(&csv, &params).unwrap();
        assert_eq!(counts(&fps), [Some(6), Some(6)]);
        // With the default tab delimiter the first line is the title, and
        // "CCN,b" is not SMILES.
        assert_eq!(counts(&read(&csv, num_threads).unwrap()), []);
    }
}

#[test]
fn gzipped_files() {
    let sdf = write("six.sdf.gz", &gzip(sd_file(&SMILES).as_bytes()));
    let sdfgz = write("six.sdfgz", &gzip(sd_file(&SMILES).as_bytes()));
    let smi = write("six.smi.gz", &gzip(smiles_file(&SMILES[..4]).as_bytes()));
    let tdt = write("six.tdt.gz", &gzip(tdt_file(&SMILES[..4]).as_bytes()));
    for num_threads in [1, 4] {
        assert_eq!(on_bits(&read(&sdf, num_threads).unwrap()), morgan3());
        assert_eq!(on_bits(&read(&sdfgz, num_threads).unwrap()), morgan3());
    }
    // RDKit's one-thread SMILES and TDT suppliers seek, which a gzip stream
    // cannot, and read nothing, as in Python.
    assert_eq!(counts(&read(&smi, 1).unwrap()), []);
    assert_eq!(
        counts(&read(&smi, 2).unwrap()),
        [Some(6), Some(4), Some(7), Some(7)]
    );
    assert_eq!(counts(&read(&tdt, 1).unwrap()), []);
    assert_eq!(counts(&read(&tdt, 2).unwrap()), []);
}

#[test]
fn unreadable_records_depend_on_the_threads() {
    // As Python: one thread keeps every record but a final unreadable one;
    // several keep records up to the last with a fingerprint.
    let trailing = write("trailing.smi", b"SMILES\tName\nCCO\ta\nfoo\tb\nbar\tc\n");
    assert_eq!(counts(&read(&trailing, 1).unwrap()), [Some(6), None]);
    assert_eq!(counts(&read(&trailing, 2).unwrap()), [Some(6)]);
    let leading = write("leading.smi", b"SMILES\tName\nfoo\ta\nCCO\tb\n");
    for num_threads in [1, 2] {
        assert_eq!(
            counts(&read(&leading, num_threads).unwrap()),
            [None, Some(6)]
        );
    }
    let middle = write(
        "middle.smi",
        b"SMILES\tName\nCCO\ta\nC1CC\tb\nCCN\tc\nfoo\td\n",
    );
    for num_threads in [1, 2, 0] {
        assert_eq!(
            counts(&read(&middle, num_threads).unwrap()),
            [Some(6), None, Some(6)]
        );
    }
    let bad = write("bad.smi", b"SMILES\tName\nfoo\ta\nbar\tb\n");
    assert_eq!(counts(&read(&bad, 1).unwrap()), [None]);
    assert_eq!(counts(&read(&bad, 2).unwrap()), []);
}

#[test]
fn generators() {
    let smi = write("gen.smi", smiles_file(&SMILES).as_bytes());
    let morgan = get_morgan_generator_with_params(&GetMorganGeneratorParams::default()).unwrap();
    let rdkit = get_rdkit_fp_generator_with_params(&GetRDKitFPGeneratorParams {
        max_path: 3,
        fp_size: 512,
        ..Default::default()
    })
    .unwrap();
    // Feature invariants match substructures, so each fingerprint holds the
    // crate's Ctrl-C lock (U4).
    let features = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 2,
        atom_invariants_generator: Some(get_morgan_feature_atom_inv_gen()),
        ..Default::default()
    })
    .unwrap();
    let feature_bits: [&[u32]; 6] = [
        &[0, 3, 805, 1085, 2041],
        &[4, 598, 2017],
        &[0, 19, 519, 792, 901, 1639],
        &[0, 1, 2, 32, 349, 1048, 1117, 1950],
        &[0, 18, 468, 546, 671, 822, 1085],
        &[
            1, 2, 3, 4, 32, 321, 380, 535, 598, 808, 809, 880, 1011, 1048, 1502, 1925, 1950, 2017,
        ],
    ];
    for num_threads in [1, 4] {
        let params = |generator| GetFingerprintsForMolsInFileParams {
            generator: Some(generator),
            options: SupplierOptions {
                num_threads,
                ..Default::default()
            },
        };
        let fps = get_fingerprints_for_mols_in_file_with_params(&smi, &params(&morgan)).unwrap();
        // The 32- and 64-bit Morgan generators set the same bits.
        assert_eq!(on_bits(&fps), morgan3());
        let fps = get_fingerprints_for_mols_in_file_with_params(&smi, &params(&rdkit)).unwrap();
        assert_eq!(
            counts(&fps),
            [Some(6), Some(5), Some(14), Some(14), Some(12), Some(37)]
        );
        let fps = get_fingerprints_for_mols_in_file_with_params(&smi, &params(&features)).unwrap();
        let expected: Vec<_> = feature_bits.iter().map(|b| Some(b.to_vec())).collect();
        assert_eq!(on_bits(&fps), expected);
    }
}

#[test]
fn fingerprint_errors_depend_on_the_threads() {
    // Atom pairs from 3D distances need a conformer, which SMILES lack. One
    // thread reports the error, as Python raises; several leave the records
    // without fingerprints, as Python's threads do.
    let smi = write("ap3d.smi", smiles_file(&SMILES[..2]).as_bytes());
    let generator = get_atom_pair_generator_with_params(&GetAtomPairGeneratorParams {
        use_2d: false,
        ..Default::default()
    })
    .unwrap();
    let params = |num_threads| GetFingerprintsForMolsInFileParams {
        generator: Some(&generator),
        options: SupplierOptions {
            num_threads,
            ..Default::default()
        },
    };
    assert!(matches!(
        get_fingerprints_for_mols_in_file_with_params(&smi, &params(1)),
        Err(Error::Value { .. })
    ));
    let fps = get_fingerprints_for_mols_in_file_with_params(&smi, &params(2)).unwrap();
    assert!(fps.is_empty());
}

#[test]
fn files_rdkit_cannot_read() {
    let is_os = |result: Result<_, Error>, text: &str| match result {
        Err(Error::Os { message }) => assert!(message.contains(text), "{message}"),
        other => panic!("{other:?}"),
    };
    is_os(
        get_fingerprints_for_mols_in_file(common::temp_path("missing.sdf")),
        "Bad input file",
    );
    is_os(
        get_fingerprints_for_mols_in_file(common::temp_path("missing.sdf.gz")),
        "Invalid input file",
    );
    is_os(
        get_fingerprints_for_mols_in_file(write("empty.sdf", b"")),
        "Invalid input file",
    );
    is_os(
        get_fingerprints_for_mols_in_file(write("notgz.sdf.gz", b"not gzip at all\n")),
        "Invalid input file",
    );
    is_os(
        get_fingerprints_for_mols_in_file(write("x.mol", b"x\n")),
        "Unsupported structure or compression extension",
    );
    is_os(
        get_fingerprints_for_mols_in_file(write("x.zst", b"x\n")),
        "Unsupported compression extension",
    );
}

#[test]
fn mae_files_are_refused() {
    // RDKit would read them with its MaeMolSupplier, which the crate does not
    // bind.
    for (name, bytes) in [
        ("x.mae", b"garbage\n".to_vec()),
        ("x.maegz", gzip(b"garbage\n")),
    ] {
        let path = write(name, &bytes);
        assert!(matches!(
            get_fingerprints_for_mols_in_file(&path),
            Err(Error::Value { .. })
        ));
    }
}

#[test]
fn name_columns_below_minus_one_are_refused() {
    // RDKit's SMILES suppliers would read before a line's fields.
    let smi = write("names.smi", smiles_file(&SMILES).as_bytes());
    for num_threads in [1, 2] {
        let params = GetFingerprintsForMolsInFileParams {
            options: SupplierOptions {
                num_threads,
                name_column: -2,
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(matches!(
            get_fingerprints_for_mols_in_file_with_params(&smi, &params),
            Err(Error::Value { .. })
        ));
    }
}

#[test]
fn long_stereo_group_lines_are_refused() {
    // RDKit's mol-block regex would recurse once per character (U12), on
    // RDKit's threads too.
    let block = mol_to_v3k_mol_block(&mol_from_smiles("C[C@H](F)Cl |&1:1|").unwrap()).unwrap();
    let long = block.replace(
        "MDLV30/STERAC1 ATOMS=(1 2)",
        &format!("MDLV30/STERAC1 ATOMS=(1 2{})", " ".repeat(100_000)),
    );
    let sd = format!("{block}$$$$\n{long}$$$$\n");
    let plain = write("long.sdf", sd.as_bytes());
    let gzipped = write("long.sdf.gz", &gzip(sd.as_bytes()));
    for path in [&plain, &gzipped] {
        for num_threads in [1, 2] {
            assert!(matches!(read(path, num_threads), Err(Error::Value { .. })));
        }
    }
}

/// A V3000 record for `*CC` with an explicit hydrogen, whose `*` bond has the
/// ENDPTS value `end_points`.
fn end_points_record(end_points: &str) -> String {
    format!(
        "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
         M  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\n\
         M  V30 1 C 0 0 0 0\nM  V30 2 C 1.5 0 0 0\nM  V30 3 H 2 1 0 0\n\
         M  V30 4 * 0 1 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\n\
         M  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 4 1 ENDPTS={end_points} ATTACH=ANY\n\
         M  V30 END BOND\nM  V30 END CTAB\nM  END\n$$$$\n"
    )
}

#[test]
fn bad_bond_end_points_are_refused() {
    // Removing the hydrogen rewrites each ENDPTS value, which RWMol converts
    // to unsigned int unchecked (U12), on RDKit's threads too. Python's bits.
    let good = format!(
        "{}{}",
        end_points_record("(1 2)"),
        end_points_record("(2 2 3)")
    );
    let bits = Some(vec![80, 226, 247, 294, 1057, 1665]);
    let path = write("end-points.sdf", good.as_bytes());
    for num_threads in [1, 2] {
        let fps = read(&path, num_threads).unwrap();
        assert_eq!(on_bits(&fps), [bits.clone(), bits.clone()]);
    }
    let bad = format!(
        "{}{}",
        end_points_record("(1 2)"),
        end_points_record("(1 1e20)")
    );
    let plain = write("bad-end-points.sdf", bad.as_bytes());
    let gzipped = write("bad-end-points.sdf.gz", &gzip(bad.as_bytes()));
    for path in [&plain, &gzipped] {
        for num_threads in [1, 2] {
            assert!(matches!(read(path, num_threads), Err(Error::Value { .. })));
        }
    }
}
