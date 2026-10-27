// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::rd_base;

#[test]
fn rdkit_version_is_the_supported_series() {
    assert!(rd_base::rdkit_version().starts_with("2024.03."));
}

#[test]
fn boost_version_is_empty_or_boost_lib_version() {
    // RDKit fills this from CMake's Boost_LIB_VERSION ("1_90"), which a
    // CONFIG-mode find_package(Boost) leaves unset, as in conda-forge's build.
    let version = rd_base::boost_version();
    assert!(version.is_empty() || version.contains('_'), "{version}");
}

#[test]
fn build_is_described() {
    assert!(!rd_base::rdkit_build().is_empty());
}

// The logs are process-wide, so everything that touches them is in one test.
#[test]
fn logs() {
    rd_base::disable_log("rdApp.error");
    assert!(rd_base::log_status().contains("rdApp.error:disabled"));
    rd_base::enable_log("rdApp.error");
    assert!(rd_base::log_status().contains("rdApp.error:enabled"));

    // RDKit 2024.03 has no way to capture a log, so these only run.
    rd_base::log_error_msg("rd_base test error");
    rd_base::log_message("rdApp.error", "second message");
    rd_base::log_message("not.a.log", "ignored");
    // RDKit reports an unparsable SMILES to its error log.
    assert!(mol_from_smiles("C1CC").is_err());

    {
        let _blocked = rd_base::BlockLogs::default();
        assert!(rd_base::log_status().contains("rdApp.error:disabled"));
        rd_base::log_error_msg("not shown");
    }
    assert!(rd_base::log_status().contains("rdApp.error:enabled"));

    rd_base::log_to_cpp_streams();
    rd_base::log_debug_msg("debug is disabled by default");
    rd_base::log_info_msg("info is disabled by default");
    rd_base::log_warning_msg("rd_base test warning");
}

#[test]
fn seeding_is_accepted() {
    rd_base::seed_random_number_generator(42);
}
