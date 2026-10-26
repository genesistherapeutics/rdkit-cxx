// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.RDLogger`, a Python module: a logger over RDKit's own
//! logs, which `rdkit_cxx::rd_base` binds.
//!
//! Python formats a message with its extra arguments (`msg % args`); Rust
//! callers format the message first. Python's `exc_info` keyword appends
//! the Python exception being handled, which Rust has none of.
//!
//! ```no_run
//! use rdkit_cxx::rd_logger::{Logger, error};
//!
//! let lg = Logger;
//! lg.set_level(error()).unwrap();
//! lg.warning("not shown");
//! ```

use crate::rd_base::{disable_log, enable_log, log_message};
use crate::{Error, Result};

const LEVELS: [&str; 4] = ["rdApp.debug", "rdApp.info", "rdApp.warning", "rdApp.error"];

/// The level of debug messages. Mirrors `RDLogger.DEBUG`.
pub fn debug() -> i32 {
    0
}

/// The level of informational messages. Mirrors `RDLogger.INFO`.
pub fn info() -> i32 {
    1
}

/// The level of warnings. Mirrors `RDLogger.WARNING`.
pub fn warning() -> i32 {
    2
}

/// The level of errors. Mirrors `RDLogger.ERROR`.
pub fn error() -> i32 {
    3
}

/// The level of critical errors, which RDKit logs as errors. Mirrors
/// `RDLogger.CRITICAL`.
pub fn critical() -> i32 {
    4
}

/// A logger writing to RDKit's logs. It holds no state: its level is
/// RDKit's, shared by every logger. Mirrors `RDLogger.logger`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Logger;

impl Logger {
    /// Writes `msg` and a newline to the RDKit log `dest` (such as
    /// `rdApp.error`). Mirrors `logger.logIt`.
    pub fn log_it(&self, dest: &str, msg: &str) {
        log_message(dest, &format!("{msg}\n"));
    }

    /// Logs `DEBUG: msg` to `rdApp.debug`. Mirrors `logger.debug`.
    pub fn debug(&self, msg: &str) {
        self.log_it("rdApp.debug", &format!("DEBUG: {msg}"));
    }

    /// Logs `ERROR: msg` to `rdApp.error`. Mirrors `logger.error`.
    pub fn error(&self, msg: &str) {
        self.log_it("rdApp.error", &format!("ERROR: {msg}"));
    }

    /// Logs `INFO: msg` to `rdApp.info`. Mirrors `logger.info`.
    pub fn info(&self, msg: &str) {
        self.log_it("rdApp.info", &format!("INFO: {msg}"));
    }

    /// Logs `WARNING: msg` to `rdApp.warning`. Mirrors `logger.warning`.
    pub fn warning(&self, msg: &str) {
        self.log_it("rdApp.warning", &format!("WARNING: {msg}"));
    }

    /// Logs `CRITICAL: msg` to `rdApp.error`. Mirrors `logger.critical`.
    pub fn critical(&self, msg: &str) {
        self.log_it("rdApp.error", &format!("CRITICAL: {msg}"));
    }

    /// Enables RDKit's logs at level `val` (one of [`debug`] to
    /// [`critical`]) and above, and disables those below. Mirrors
    /// `logger.setLevel`.
    ///
    /// As Python indexes its list of logs, a negative level counts from the
    /// error log back; `Error::Index` below -4, as Python raises.
    pub fn set_level(&self, val: i32) -> Result<()> {
        let log = |i: i32| -> Result<&str> {
            let idx = if i < 0 { i + 4 } else { i };
            usize::try_from(idx)
                .ok()
                .and_then(|idx| LEVELS.get(idx).copied())
                .ok_or(Error::index(i, 4))
        };
        for i in val..4 {
            enable_log(log(i)?);
        }
        for i in 0..val {
            disable_log(log(i)?);
        }
        Ok(())
    }
}
