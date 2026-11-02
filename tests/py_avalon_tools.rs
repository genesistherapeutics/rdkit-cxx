// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::avalon::py_avalon_tools::*;

// Expected values come from Python's pyAvalonTools.

#[test]
fn constants() {
    assert_eq!(avalon_sss_bits(), 32767);
    assert_eq!(avalon_similarity_bits(), 15761407);
}

#[test]
fn struchk_flags() {
    let flags = [
        (StruChkFlag::bad_molecule, 1),
        (StruChkFlag::alias_conversion_failed, 2),
        (StruChkFlag::transformed, 4),
        (StruChkFlag::fragments_found, 8),
        (StruChkFlag::either_warning, 16),
        (StruChkFlag::stereo_error, 32),
        (StruChkFlag::dubious_stereo_removed, 64),
        (StruChkFlag::atom_clash, 128),
        (StruChkFlag::atom_check_failed, 256),
        (StruChkFlag::size_check_failed, 512),
        (StruChkFlag::recharged, 1024),
        (StruChkFlag::stereo_forced_bad, 2048),
        (StruChkFlag::stereo_transformed, 4096),
        (StruChkFlag::template_transformed, 8192),
    ];
    for (flag, bits) in flags {
        assert_eq!(flag.bits(), bits);
    }
    assert_eq!(
        (StruChkFlag::bad_molecule | StruChkFlag::transformed).bits(),
        5
    );
    assert_eq!(
        format!("{:?}", StruChkFlag::atom_clash),
        "StruChkFlag::atom_clash"
    );
    assert_eq!(StruChkResult::success.bits(), 0);
    assert_eq!(StruChkResult::bad_set.bits(), 2979);
    assert_eq!(StruChkResult::transformed_set.bits(), 29788);
    // The masks group the flags as struchk.h does.
    let bad = StruChkFlag::bad_molecule
        | StruChkFlag::alias_conversion_failed
        | StruChkFlag::stereo_error
        | StruChkFlag::stereo_forced_bad
        | StruChkFlag::atom_clash
        | StruChkFlag::atom_check_failed
        | StruChkFlag::size_check_failed;
    assert_eq!(bad.bits(), StruChkResult::bad_set.bits());
    assert_eq!(
        format!("{:?}", StruChkResult::from_bits_retain(7)),
        "StruChkResult(0x7)"
    );
}
