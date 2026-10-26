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

/// Whether Python's `str.isspace` counts the character as whitespace.
pub(crate) fn is_py_space(c: char) -> bool {
    c.is_whitespace() || ('\x1c'..='\x1f').contains(&c)
}

/// Python's `str.strip()`.
pub(crate) fn py_strip(s: &str) -> &str {
    s.trim_matches(is_py_space)
}

/// Python's `str.split()`: the runs of non-whitespace.
pub(crate) fn py_split(s: &str) -> Vec<&str> {
    s.split(is_py_space).filter(|w| !w.is_empty()).collect()
}

/// Python's `int(s)` for a string: an optional sign and decimal digits,
/// single underscores allowed between digits, whitespace around.
pub(crate) fn py_int(s: &str) -> Option<i64> {
    let t = py_strip(s);
    let (sign, digits) = match t.strip_prefix('-') {
        Some(rest) => (-1, rest),
        None => (1, t.strip_prefix('+').unwrap_or(t)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
        || !digits.chars().all(|c| c.is_ascii_digit() || c == '_')
    {
        return None;
    }
    let value: i64 = digits.replace('_', "").parse().ok()?;
    Some(sign * value)
}

/// Python's `float(s)` for a string: Rust's parse, after stripping
/// whitespace, with Python's spellings of infinity and NaN, and single
/// underscores between digits.
pub(crate) fn py_float(s: &str) -> Option<f64> {
    let t = py_strip(s);
    let unsigned = t.trim_start_matches(['+', '-']);
    if t.len() - unsigned.len() > 1 {
        return None;
    }
    let lower = unsigned.to_ascii_lowercase();
    if ["inf", "infinity", "nan"].contains(&lower.as_str()) {
        let negative = t.starts_with('-');
        let v = if lower == "nan" {
            f64::NAN
        } else {
            f64::INFINITY
        };
        return Some(if negative { -v } else { v });
    }
    if !unsigned.starts_with(|c: char| c.is_ascii_digit() || c == '.') {
        return None;
    }
    let chars: Vec<char> = t.chars().collect();
    for (i, &c) in chars.iter().enumerate() {
        if c == '_' {
            let digit = |j: Option<&char>| j.is_some_and(char::is_ascii_digit);
            if !digit(i.checked_sub(1).and_then(|j| chars.get(j))) || !digit(chars.get(i + 1)) {
                return None;
            }
        }
    }
    t.replace('_', "").parse().ok()
}
