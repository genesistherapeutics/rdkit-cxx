// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.utils`.

pub mod chemutils;
pub mod fileutils;
pub mod listutils;

use std::io::{self, BufRead};

/// Reads a line as Python's text files do, with universal newlines: `\n`,
/// `\r\n`, and `\r` each end a line, which comes back ending in `\n`. An
/// empty string is the end of the input. `Error::Value`-worthy invalid
/// UTF-8 is an `io::Error` of kind `InvalidData`.
pub(crate) fn read_line_universal<R: BufRead + ?Sized>(reader: &mut R) -> io::Result<String> {
    let mut bytes = Vec::new();
    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        match buf.iter().position(|&b| b == b'\n' || b == b'\r') {
            Some(i) => {
                let ending = buf[i];
                bytes.extend_from_slice(&buf[..i]);
                reader.consume(i + 1);
                if ending == b'\r' {
                    let next = reader.fill_buf()?;
                    if next.first() == Some(&b'\n') {
                        reader.consume(1);
                    }
                }
                bytes.push(b'\n');
                break;
            }
            None => {
                let n = buf.len();
                bytes.extend_from_slice(buf);
                reader.consume(n);
            }
        }
    }
    String::from_utf8(bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
