// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::error::{__raise_for_tests as raise, __raise_unexpected_for_tests};

fn message(kind: &str) -> String {
    format!("test {kind}")
}

#[test]
fn each_rdkit_exception_maps_to_its_variant() {
    assert_eq!(
        raise("atom_valence").unwrap_err(),
        Error::AtomValence {
            atom_idx: 3,
            message: message("atom_valence")
        }
    );
    assert_eq!(
        raise("atom_kekulize").unwrap_err(),
        Error::AtomKekulize {
            atom_idx: 2,
            message: message("atom_kekulize")
        }
    );
    assert_eq!(
        raise("kekulize").unwrap_err(),
        Error::Kekulize {
            atom_indices: vec![0, 1, 2, 3, 4],
            message: message("kekulize")
        }
    );
    assert_eq!(
        raise("mol_sanitize").unwrap_err(),
        Error::MolSanitize {
            message: message("mol_sanitize")
        }
    );
    assert_eq!(
        raise("index").unwrap_err(),
        Error::Index { index: 9, len: 3 }
    );
    assert_eq!(raise("key").unwrap_err(), Error::Key { key: "nope".into() });
    assert_eq!(
        raise("value").unwrap_err(),
        Error::Value {
            message: message("value")
        }
    );
    assert_eq!(
        raise("conformer").unwrap_err(),
        Error::Value {
            message: message("conformer")
        }
    );
    assert_eq!(
        raise("bad_file").unwrap_err(),
        Error::Os {
            message: message("bad_file")
        }
    );
    assert_eq!(
        raise("smiles_parse").unwrap_err(),
        Error::Parse {
            message: message("smiles_parse")
        }
    );
    assert_eq!(
        raise("std").unwrap_err(),
        Error::Runtime {
            message: message("std")
        }
    );
    assert!(matches!(raise("invariant").unwrap_err(),
                     Error::Runtime { message } if message.contains("test invariant")));
}

#[test]
fn display_includes_the_message() {
    assert!(
        raise("value")
            .unwrap_err()
            .to_string()
            .contains("test value")
    );
}

#[test]
#[should_panic(expected = "unexpected RDKit exception: test std")]
fn unexpected_exceptions_panic() {
    __raise_unexpected_for_tests();
}
