// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! A synthon space whose fingerprint build RDKit 2025.09 cuts short. Its own
//! test binary, since it raises SIGINT.

use std::path::PathBuf;

use rdkit_cxx::chem::rd_fingerprint_generator::{
    FingerprintGenerator64, GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::chem::rd_synthon_space_search::*;
use rdkit_cxx::chem::*;

unsafe extern "C" {
    fn raise(sig: i32) -> i32;
}

const SIGINT: i32 = 2;

fn morgan(radius: u32) -> FingerprintGenerator64 {
    get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius,
        fp_size: 2048,
        ..Default::default()
    })
    .unwrap()
}

// RDKit 2025.09's searches install its SIGINT handler with
// ControlCHandler::reset(), and nothing clears the flag afterwards, so a
// SIGINT during a search leaves it set. buildSynthonFingerprints sets the
// space's fingerprint type and then returns at once for a set flag, and the
// next fingerprint search, which builds only for another type, dereferenced
// the synthons' missing fingerprints (a segfault). The wrapper clears the
// flag around each call, so the build completes.
#[test]
fn a_build_cut_short_is_not_searched() {
    let mut space = SynthonSpace::default();
    space
        .read_text_file(
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/amide_space.txt"),
        )
        .unwrap();
    let query = mol_from_smiles("c1ccccc1C(=O)N1CCCC1").unwrap();
    let pattern = mol_from_smarts("c1ccccc1C(=O)N").unwrap();
    let mut calls = 0;
    space
        .substructure_search_incremental(&pattern, |_| {
            calls += 1;
            // SAFETY: raise only sends this thread a signal, which RDKit's
            // handler, installed for the search, takes.
            unsafe { raise(SIGINT) };
            Ok(false)
        })
        .unwrap();
    assert!(calls > 0);
    let other = morgan(1);
    space.build_synthon_fingerprints(&other).unwrap();
    assert_eq!(
        space.get_synthon_fingerprint_type(),
        other.get_info_string().unwrap()
    );
    // Python's hits for the generator.
    let hits = space.fingerprint_search(&query, &other).unwrap();
    let mut names: Vec<String> = hits
        .get_hit_molecules()
        .iter()
        .map(|m| m.get_prop("_Name").unwrap())
        .collect();
    names.sort();
    assert_eq!(names, ["1-1;2-1;amide-1", "1-3;2-1;amide-1"]);
}
