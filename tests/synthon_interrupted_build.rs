// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// A SIGINT left in RDKit's Ctrl-C flag, and the synthon fingerprint build it
// would cut short. Raising a signal is process-wide, so this is its own test
// binary.

use std::path::PathBuf;

use rdkit_cxx::chem::rd_fingerprint_generator::{
    GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::chem::rd_synthon_space_search::*;
use rdkit_cxx::chem::*;

unsafe extern "C" {
    fn raise(sig: i32) -> i32;
}

const SIGINT: i32 = 2;

fn amide_space() -> SynthonSpace {
    let mut space = SynthonSpace::default();
    space
        .read_text_file(
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/amide_space.txt"),
        )
        .unwrap();
    space
}

fn morgan() -> rdkit_cxx::chem::rd_fingerprint_generator::FingerprintGenerator64 {
    get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 2,
        fp_size: 2048,
        ..Default::default()
    })
    .unwrap()
}

// RDKit 2026.03's fingerprint search calls ControlCHandler::reset(), which
// installs RDKit's handler and constructs none whose destructor would clear
// its flag, and calls the closure with it installed: the handler takes the
// SIGINT and sets the flag, which stays set after the search. Then
// buildSynthonFingerprints, which only reads the flag, returns at it after
// recording the fingerprint type, so the next fingerprint search takes the
// space's fingerprints as built and dereferences the synthons' missing ones.
#[test]
fn a_stale_ctrl_c_does_not_cut_the_fingerprint_build_short() {
    let generator = morgan();
    let query = mol_from_smiles("c1ccccc1C(=O)N1CCCC1").unwrap();
    let mut space = amide_space();
    let mut other = amide_space();
    let _ = other.fingerprint_search_incremental(&query, &generator, |_| {
        // SAFETY: raise has no preconditions; RDKit's handler takes it.
        unsafe { raise(SIGINT) };
        Ok(true)
    });
    space.build_synthon_fingerprints(&generator).unwrap();
    let results = space.fingerprint_search(&query, &generator).unwrap();
    assert_eq!(results.get_hit_molecules().len(), 1);
}
