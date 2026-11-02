// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.rdBase`.
//!
//! RDKit's logs write to the C++ streams (warnings and errors to stderr) from
//! the time the library loads, as they do once Python imports `rdBase`.

use cxx::UniquePtr;

use crate::error::expect_ok;
use crate::ffi::rd_base::ffi;

fn ascii(text: UniquePtr<cxx::CxxString>) -> String {
    text.to_str()
        .expect("RDKit's build strings are UTF-8")
        .to_owned()
}

/// The RDKit version. Mirrors `rdBase.rdkitVersion`.
pub fn rdkit_version() -> String {
    ascii(expect_ok(ffi::rdkit_version()))
}

/// The Boost version RDKit was built with. Mirrors `rdBase.boostVersion`.
pub fn boost_version() -> String {
    ascii(expect_ok(ffi::boost_version()))
}

/// The system and compiler RDKit was built with. Mirrors `rdBase.rdkitBuild`.
pub fn rdkit_build() -> String {
    ascii(expect_ok(ffi::rdkit_build()))
}

/// Disables the log named by `spec`, such as `"rdApp.error"`. Mirrors
/// `rdBase.DisableLog`.
pub fn disable_log(spec: &str) {
    expect_ok(ffi::disable_log(spec));
}

/// Enables the log named by `spec`, such as `"rdApp.error"`. Mirrors
/// `rdBase.EnableLog`.
pub fn enable_log(spec: &str) {
    expect_ok(ffi::enable_log(spec));
}

/// Which logs are enabled. Mirrors `rdBase.LogStatus`.
pub fn log_status() -> String {
    ascii(expect_ok(ffi::log_status()))
}

/// Points the logs back at the C++ streams, resetting them. Mirrors
/// `rdBase.LogToCppStreams`.
pub fn log_to_cpp_streams() {
    expect_ok(ffi::log_to_cpp_streams());
}

/// Mirrors `rdBase.LogDebugMsg`.
pub fn log_debug_msg(msg: &str) {
    expect_ok(ffi::log_debug_msg(msg));
}

/// Mirrors `rdBase.LogInfoMsg`.
pub fn log_info_msg(msg: &str) {
    expect_ok(ffi::log_info_msg(msg));
}

/// Mirrors `rdBase.LogWarningMsg`.
pub fn log_warning_msg(msg: &str) {
    expect_ok(ffi::log_warning_msg(msg));
}

/// Mirrors `rdBase.LogErrorMsg`.
pub fn log_error_msg(msg: &str) {
    expect_ok(ffi::log_error_msg(msg));
}

/// Logs to any `rdApp.*` log; other specs are ignored. Mirrors
/// `rdBase.LogMessage`.
pub fn log_message(spec: &str, msg: &str) {
    expect_ok(ffi::log_message(spec, msg));
}

/// Seeds the C random number generator some RDKit code uses. Mirrors
/// `rdBase.SeedRandomNumberGenerator`.
pub fn seed_random_number_generator(seed: u32) {
    expect_ok(ffi::seed_random_number_generator(seed));
}

/// Disables every log until dropped, then restores their state. Mirrors
/// `rdBase.BlockLogs`, a context manager in Python.
///
/// The logs are process-wide, so the guard is neither `Send` nor `Sync`.
pub struct BlockLogs {
    _guard: UniquePtr<ffi::BlockLogs>,
}

impl std::fmt::Debug for BlockLogs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BlockLogs").finish_non_exhaustive()
    }
}

impl Default for BlockLogs {
    fn default() -> BlockLogs {
        BlockLogs {
            _guard: expect_ok(ffi::block_logs_new()),
        }
    }
}

/// Captures what RDKit writes to its error log until dropped. Mirrors
/// `rdBase.CaptureErrorLog`, a context manager in Python.
pub struct CaptureErrorLog {
    capture: UniquePtr<ffi::CaptureErrorLog>,
}

impl std::fmt::Debug for CaptureErrorLog {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CaptureErrorLog").finish_non_exhaustive()
    }
}

impl Default for CaptureErrorLog {
    fn default() -> CaptureErrorLog {
        CaptureErrorLog {
            capture: expect_ok(ffi::capture_error_log_new()),
        }
    }
}

impl CaptureErrorLog {
    /// The messages captured so far, with RDKit's timestamps. Mirrors
    /// `CaptureErrorLog.messages`.
    pub fn messages(&self) -> String {
        let messages = expect_ok(ffi::capture_error_log_messages(
            self.capture.as_ref().unwrap(),
        ));
        String::from_utf8_lossy(messages.as_bytes()).into_owned()
    }
}
