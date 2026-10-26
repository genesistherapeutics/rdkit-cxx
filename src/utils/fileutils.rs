// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.utils.fileutils`, a Python module: utilities for
//! reading files.

use std::io::BufRead;

use crate::utils::read_line_universal;
use crate::{Error, Result};

/// Reads lines until one contains `match_str`, and returns it (ending in
/// its newline, as Python reads lines). Mirrors
/// `fileutils.MoveToMatchingLine(inFile, matchStr)`.
///
/// `Error::Runtime` when no line matches (Python's `NoMatchFoundError`, a
/// `RuntimeError`); `Error::Os` when reading fails, and `Error::Value` for
/// a line that is not UTF-8.
pub fn move_to_matching_line<R: BufRead>(in_file: &mut R, match_str: &str) -> Result<String> {
    move_to_matching_line_with_full_match(in_file, match_str, false)
}

/// As [`move_to_matching_line`]; with `full_match`, the line must start
/// with `match_str`. Mirrors `fileutils.MoveToMatchingLine(inFile,
/// matchStr, fullMatch)`.
pub fn move_to_matching_line_with_full_match<R: BufRead>(
    in_file: &mut R,
    match_str: &str,
    full_match: bool,
) -> Result<String> {
    loop {
        let line = read_line(in_file)?;
        if line.is_empty() {
            return Err(Error::Runtime {
                message: match_str.to_owned(),
            });
        }
        let matched = match line.find(match_str) {
            Some(idx) => !full_match || idx == 0,
            None => false,
        };
        if matched {
            return Ok(line);
        }
    }
}

/// A line as Python's `readline` reads it, with I/O errors as `Error::Os`
/// and invalid UTF-8 as `Error::Value`.
pub(crate) fn read_line<R: BufRead + ?Sized>(reader: &mut R) -> Result<String> {
    read_line_universal(reader).map_err(|e| {
        if e.kind() == std::io::ErrorKind::InvalidData {
            Error::Value {
                message: e.to_string(),
            }
        } else {
            Error::Os {
                message: e.to_string(),
            }
        }
    })
}
