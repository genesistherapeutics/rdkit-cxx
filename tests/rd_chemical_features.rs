// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chemical_features::*;
use rdkit_cxx::geometry::rd_geometry::Point3D;

// Expected values come from Python's rdChemicalFeatures.

fn point(x: f64, y: f64, z: f64) -> Point3D {
    Point3D { x, y, z }
}

// FreeChemicalFeature("Donor", "HDonor1", Point3D(1, 2, 3.5), 7).__getstate__()
const DONOR: &str = "200000000700000006000000446f6e6f72000800000048446f6e6f72310000\
                     0000000000f03f00000000000000400000000000000c40";

fn hex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn le32(v: u32) -> [u8; 4] {
    v.to_le_bytes()
}

fn doubles(values: [f64; 3]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

#[test]
fn constructors_and_accessors() {
    let feat = FreeChemicalFeature::default();
    assert_eq!(feat.get_id(), -1);
    assert_eq!(feat.get_family(), "");
    assert_eq!(feat.get_type(), "");
    assert_eq!(feat.get_pos(), point(0.0, 0.0, 0.0));

    let mut feat = FreeChemicalFeature::new_with_id("Donor", "HDonor1", point(1.0, 2.0, 3.5), 7);
    assert_eq!(feat.get_id(), 7);
    assert_eq!(feat.get_family(), "Donor");
    assert_eq!(feat.get_type(), "HDonor1");
    assert_eq!(feat.get_pos(), point(1.0, 2.0, 3.5));
    assert_eq!(feat.to_pickle(), hex(DONOR));
    let copy = feat.clone();
    feat.set_id(3);
    feat.set_family("Acceptor");
    feat.set_type("A1");
    feat.set_pos(point(-1.0, 0.5, 0.0));
    assert_eq!(feat.get_id(), 3);
    assert_eq!(feat.get_family(), "Acceptor");
    assert_eq!(feat.get_type(), "A1");
    assert_eq!(feat.get_pos(), point(-1.0, 0.5, 0.0));
    assert_eq!(copy.get_family(), "Donor");

    let feat = FreeChemicalFeature::new("Aromatic", "Arom6", point(0.0, 1.0, 0.0));
    assert_eq!(feat.get_id(), -1);
    let feat = FreeChemicalFeature::from_family("Acceptor", point(-1.0, 0.0, 0.0));
    assert_eq!(
        (feat.get_id(), feat.get_family(), feat.get_type()),
        (-1, "Acceptor".to_string(), String::new())
    );
    assert!(format!("{feat:?}").contains("Acceptor"));
}

#[test]
fn pickles_round_trip() {
    let feat = FreeChemicalFeature::from_pickle(&hex(DONOR)).unwrap();
    assert_eq!(feat.get_id(), 7);
    assert_eq!(feat.get_family(), "Donor");
    assert_eq!(feat.get_type(), "HDonor1");
    assert_eq!(feat.get_pos(), point(1.0, 2.0, 3.5));

    // Version 1 has no id.
    let mut v1 = le32(0x10).to_vec();
    v1.extend(le32(4));
    v1.extend(b"Fam\0");
    v1.extend(le32(2));
    v1.extend(b"T\0");
    v1.extend(doubles([1.0, 2.0, 3.0]));
    let feat = FreeChemicalFeature::from_pickle(&v1).unwrap();
    assert_eq!(
        (feat.get_id(), feat.get_family(), feat.get_type()),
        (-1, "Fam".to_string(), "T".to_string())
    );

    // Strings end at their first NUL, and trailing bytes are ignored.
    let mut odd = le32(0x20).to_vec();
    odd.extend(3i32.to_le_bytes());
    odd.extend(le32(5));
    odd.extend(b"ab\0cd");
    odd.extend(le32(1));
    odd.extend(b"\0");
    odd.extend(doubles([1.0, 2.0, 3.0]));
    odd.extend(b"extra");
    let feat = FreeChemicalFeature::from_pickle(&odd).unwrap();
    assert_eq!(feat.get_id(), 3);
    assert_eq!(feat.get_family(), "ab");
    assert_eq!(feat.get_type(), "");
    assert_eq!(feat.get_pos(), point(1.0, 2.0, 3.0));
}

fn refused(pickle: &[u8]) -> bool {
    matches!(
        FreeChemicalFeature::from_pickle(pickle),
        Err(Error::Value { .. })
    )
}

#[test]
fn refuses_pickles_rdkit_reads_past() {
    // RDKit throws a const char * for an unknown version.
    assert!(refused(&le32(1)));
    assert!(refused(b""));
    assert!(refused(&le32(0x20)));
    // A family longer than the pickle: RDKit takes the string up to a NUL
    // in bytes it never read.
    let mut short = le32(0x20).to_vec();
    short.extend(3i32.to_le_bytes());
    short.extend(le32(100));
    short.extend(b"abc");
    assert!(refused(&short));
    // A family without a NUL: RDKit reads past its buffer.
    let mut no_nul = le32(0x20).to_vec();
    no_nul.extend(3i32.to_le_bytes());
    no_nul.extend(le32(3));
    no_nul.extend(b"abc");
    no_nul.extend(le32(1));
    no_nul.extend(b"\0");
    no_nul.extend(doubles([1.0, 2.0, 3.0]));
    assert!(refused(&no_nul));
    // A zero-length string reads past an empty buffer.
    let mut empty = le32(0x10).to_vec();
    empty.extend(le32(0));
    empty.extend(le32(1));
    empty.extend(b"\0");
    empty.extend(doubles([1.0, 2.0, 3.0]));
    assert!(refused(&empty));
    // The position is missing.
    let mut no_pos = le32(0x10).to_vec();
    no_pos.extend(le32(1));
    no_pos.extend(b"\0");
    no_pos.extend(le32(1));
    no_pos.extend(b"\0");
    no_pos.extend([0u8; 23]);
    assert!(refused(&no_pos));
    // Names that are not UTF-8.
    let mut latin1 = le32(0x10).to_vec();
    latin1.extend(le32(2));
    latin1.extend(b"\xe9\0");
    latin1.extend(le32(1));
    latin1.extend(b"\0");
    latin1.extend(doubles([1.0, 2.0, 3.0]));
    assert!(refused(&latin1));
}

#[test]
fn features_cross_threads() {
    let feat = FreeChemicalFeature::new("Donor", "D", point(1.0, 0.0, 0.0));
    let handles: Vec<_> = (0..4)
        .map(|i| {
            let mut feat = feat.clone();
            std::thread::spawn(move || {
                feat.set_id(i);
                FreeChemicalFeature::from_pickle(&feat.to_pickle())
                    .unwrap()
                    .get_id()
            })
        })
        .collect();
    let ids: Vec<i32> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    assert_eq!(ids, vec![0, 1, 2, 3]);
}
