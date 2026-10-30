// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Checks the port of rdkit.RDLogger against Python's in RDKit 2025.09.6.
// RDKit's logs are process-wide, so this file holds one test.

use rdkit_cxx::Error;
use rdkit_cxx::rd_base::log_status;
use rdkit_cxx::rd_logger::*;

#[test]
fn logger() {
    assert_eq!(
        (debug(), info(), warning(), error(), critical()),
        (0, 1, 2, 3, 4)
    );
    let lg = Logger;
    lg.set_level(error()).unwrap();
    assert_eq!(
        log_status(),
        "rdApp.debug:disabled\nrdApp.info:disabled\nrdApp.warning:disabled\nrdApp.error:enabled"
    );
    lg.set_level(debug()).unwrap();
    assert_eq!(
        log_status(),
        "rdApp.debug:enabled\nrdApp.info:enabled\nrdApp.warning:enabled\nrdApp.error:enabled"
    );
    // RDKit 2025.09 has no CaptureErrorLog to read the messages back.
    lg.error("first");
    lg.critical("second");
    lg.log_it("rdApp.error", "third");
    // As Python indexes its list, -1 is the error log, and below -4 fails.
    assert!(matches!(lg.set_level(-5), Err(Error::Index { .. })));
    lg.set_level(warning()).unwrap();
    assert!(
        log_status()
            .starts_with("rdApp.debug:disabled\nrdApp.info:disabled\nrdApp.warning:enabled")
    );
    lg.set_level(info()).unwrap();
}
