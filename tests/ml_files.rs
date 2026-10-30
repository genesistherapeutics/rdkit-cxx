// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Checks the port of rdkit.ML.files against Python's in RDKit 2024.09.6,
// whose ReadDataFile fails (it assigns a map object to a NumPy row).

mod common;

use common::temp_path;
use rdkit_cxx::Error;
use rdkit_cxx::ml::files::*;

#[test]
fn re_file_matches_python() {
    let path = temp_path("refile.txt");
    std::fs::write(
        &path,
        "# comment\n1.0 2.0 3.0 # x\n\n  \n4 5 6\r\n7;8 9 10\n",
    )
    .unwrap();
    let mut f = ReFile::new(&path).unwrap();
    assert_eq!(f.readline().unwrap(), "1.0 2.0 3.0");
    assert_eq!(f.readlines().unwrap(), ["4 5 6", "7;8 9 10"]);
    assert_eq!(f.readline().unwrap(), "");
    f.rewind().unwrap();
    assert_eq!(f.readline().unwrap(), "1.0 2.0 3.0");
    let mut f = ReFile::new_with_params(
        &path,
        &ReFileParams {
            comment: ";".to_owned(),
            trailer: Some(String::new()),
        },
    )
    .unwrap();
    assert_eq!(
        f.readlines().unwrap(),
        ["# comment", "1.0 2.0 3.0 # x", "4 5 6", "7"]
    );
    assert!(matches!(
        ReFile::new(temp_path("missing.txt")),
        Err(Error::Os { .. })
    ));
}

#[test]
fn read_data_file() {
    let path = temp_path("data.txt");
    std::fs::write(&path, "# y x1 x2\n1.5 2 3 # first\n\n-2 4.25 6\n").unwrap();
    let (ind, dep) = rdkit_cxx::ml::files::read_data_file(&path).unwrap();
    assert_eq!(ind, [vec![2.0, 3.0], vec![4.25, 6.0]]);
    assert_eq!(dep, [1.5, -2.0]);
    let (ind, dep) = read_data_file_with_params(
        &path,
        &ReadDataFileParams {
            dep_var_col: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(ind, [vec![1.5, 2.0], vec![-2.0, 4.25]]);
    assert_eq!(dep, [3.0, 6.0]);
    let ints = ReadDataFileParams {
        data_type: DataType::Int,
        ..Default::default()
    };
    assert!(matches!(
        read_data_file_with_params(&path, &ints),
        Err(Error::Value { .. })
    ));
    std::fs::write(&path, "1 2\n3 4 5\n").unwrap();
    assert!(matches!(
        rdkit_cxx::ml::files::read_data_file(&path),
        Err(Error::Value { .. })
    ));
}
