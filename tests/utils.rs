// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/utils/UnitTestUtils.py, with the doctests of listutils, and
// checks rdkit_cxx::rd_paths, with values from Python's modules in RDKit
// 2024.09.6.

use std::io::Cursor;

use rdkit_cxx::Error;
use rdkit_cxx::rd_paths;
use rdkit_cxx::utils::chemutils::*;
use rdkit_cxx::utils::fileutils::*;
use rdkit_cxx::utils::listutils::compact_list_repr;

#[test]
fn fileutils() {
    let text = "class NoMatchFoundError(RuntimeError):\n  pass\n\n\ndef MoveToMatchingLine():\n";
    let mut in_file = Cursor::new(text);
    assert_eq!(
        move_to_matching_line(&mut in_file, "pass").unwrap(),
        "  pass\n"
    );
    assert!(matches!(
        move_to_matching_line(&mut in_file, "pass"),
        Err(Error::Runtime { .. })
    ));
    let mut in_file = Cursor::new(text);
    assert!(matches!(
        move_to_matching_line_with_full_match(&mut in_file, "pass", true),
        Err(Error::Runtime { .. })
    ));
    let mut in_file = Cursor::new("a\r\nmatch here\rnext\n");
    assert_eq!(
        move_to_matching_line_with_full_match(&mut in_file, "match", true).unwrap(),
        "match here\n"
    );
}

#[test]
fn chemutils() {
    let one = |e: &str| (e.to_owned(), 1.0);
    assert_eq!(split_composition("Fe").unwrap(), [one("Fe")]);
    assert_eq!(
        split_composition("Fe3Al").unwrap(),
        [("Fe".to_owned(), 3.0), one("Al")]
    );
    assert_eq!(
        split_composition("Fe99PdAl").unwrap(),
        [("Fe".to_owned(), 99.0), one("Pd"), one("Al")]
    );
    assert_eq!(
        split_composition("TiNiSiSO12P").unwrap(),
        [
            one("Ti"),
            one("Ni"),
            one("Si"),
            one("S"),
            ("O".to_owned(), 12.0),
            one("P")
        ]
    );
    assert_eq!(
        split_composition("Ca2.5Ox").unwrap(),
        [("Ca".to_owned(), 2.5), one("Ox")]
    );
    assert!(matches!(
        split_composition("C1.2.3"),
        Err(Error::Value { .. })
    ));

    let temp = [
        "[Xe] 4f^12 6s^2",
        "[Xe] 4f^14 5d^6 6s^2",
        "[Xe] 4f^14 5d^10 6s^2",
        "[Xe] 4f^14 5d^10 6s^2 6p^1",
        "[Xe] 5d^10",
    ];
    for (d, f, expected) in [
        (true, true, [14, 8, 2, 3, 10]),
        (false, true, [14, 8, 12, 13, 10]),
        (true, false, [14, 22, 16, 17, 10]),
        (false, false, [14, 22, 26, 27, 10]),
    ] {
        let params = ConfigToNumElectronsParams {
            ignore_full_d: d,
            ignore_full_f: f,
        };
        let got: Vec<i64> = temp
            .iter()
            .map(|t| config_to_num_electrons_with_params(t, &params).unwrap())
            .collect();
        assert_eq!(got, expected);
    }
    assert_eq!(config_to_num_electrons("[Ne] 3s^2 3p^4").unwrap(), 6);
    assert!(matches!(
        config_to_num_electrons("[Ne] 3s2"),
        Err(Error::Index { .. })
    ));
    assert!(matches!(
        config_to_num_electrons("[Ne] 3s^x"),
        Err(Error::Value { .. })
    ));
}

#[test]
fn listutils() {
    assert_eq!(compact_list_repr(&[0, 1, 1, 1, 1, 0]), "[0]+[1]*4+[0]");
    assert_eq!(
        compact_list_repr(&[0, 1, 1, 2, 1, 1]),
        "[0]+[1]*2+[2]+[1]*2"
    );
    assert_eq!(compact_list_repr::<i32>(&[]), "[]");
    assert_eq!(compact_list_repr(&[0, 1, 1, 1, 1]), "[0]+[1]*4");
    assert_eq!(compact_list_repr(&[1.0, 1.0, 2.5]), "[1.0]*2+[2.5]");
    // Python gives "['f']+['o']*2".
    assert_eq!(compact_list_repr(&['f', 'o', 'o']), "['f']+['o']*2");
    assert_eq!(compact_list_repr(&["f", "o", "o"]), "[\"f\"]+[\"o\"]*2");
}

#[test]
fn paths() {
    let data = rd_paths::rd_data_dir();
    assert!(data.ends_with("share/RDKit/Data"));
    assert!(data.join("Crippen.txt").is_file(), "{}", data.display());
    assert_eq!(data.parent(), rd_paths::rd_docs_dir().parent());
    assert!(rd_paths::rd_proj_dir().ends_with("share/RDKit/Projects"));
    assert!(rd_paths::rd_contrib_dir().ends_with("share/RDKit/Contrib"));
}
