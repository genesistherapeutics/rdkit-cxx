// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// A SIGINT after a synthon search, in its own test binary: RDKit's nested
// ControlCHandlers leave its signal handler installed once they are gone, so
// the signal sets the Ctrl-C flag with no handler to clear it. RDKit's
// buildSynthonFingerprints then returns before any fingerprint, with the
// space's fingerprint type set, and its fingerprint search, which clears the
// flag first, dereferences the missing fingerprints (RDKit 2026.09.1's
// Python segfaults). The wrapper clears the flag around the build.

use std::path::PathBuf;

use rdkit_cxx::chem::rd_fingerprint_generator::{
    GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::chem::rd_synthon_space_search::*;
use rdkit_cxx::chem::*;

const SIGINT: i32 = 2;
const SIG_DFL: usize = 0;
const SIG_IGN: usize = 1;

unsafe extern "C" {
    fn signal(signum: i32, handler: usize) -> usize;
    fn raise(sig: i32) -> i32;
}

/// SIGINT's disposition, read by setting and restoring it.
fn sigint_handler() -> usize {
    // SAFETY: SIG_IGN is a valid disposition, and the previous one is put
    // back as it was.
    unsafe {
        let current = signal(SIGINT, SIG_IGN);
        signal(SIGINT, current);
        current
    }
}

#[test]
fn a_sigint_after_a_search_does_not_cut_a_fingerprint_build_short() {
    // Should RDKit leave no handler of its own, the signal is ignored.
    // SAFETY: as in sigint_handler.
    unsafe { signal(SIGINT, SIG_IGN) };
    let mut space = SynthonSpace::default();
    space
        .read_text_file(
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/amide_space.txt"),
        )
        .unwrap();
    let query = mol_from_smarts("c1ccccc1C(=O)N").unwrap();
    assert!(
        !space
            .substructure_search(&query)
            .unwrap()
            .get_hit_molecules()
            .is_empty()
    );
    let handler = sigint_handler();
    assert!(
        handler != SIG_IGN && handler != SIG_DFL,
        "RDKit's handler is left installed"
    );
    // SAFETY: the handler installed is RDKit's, which only sets its flag.
    assert_eq!(unsafe { raise(SIGINT) }, 0);

    let generator = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 2,
        fp_size: 2048,
        ..Default::default()
    })
    .unwrap();
    space.build_synthon_fingerprints(&generator).unwrap();
    let query = mol_from_smiles("c1ccccc1C(=O)N1CCCC1").unwrap();
    let hits = space
        .fingerprint_search(&query, &generator)
        .unwrap()
        .get_hit_molecules();
    assert_eq!(hits.len(), 1);
}
