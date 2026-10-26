// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! RDKit's exceptions as a Rust error type.
//!
//! The C++ wrappers catch every exception and encode its type and payload in
//! the message cxx carries (`wrapper/include/rdkit-cxx/errors.h`), so typed
//! fields such as a sanitization error's atom indices survive the crossing.

use std::fmt;

use crate::ffi::errors::ffi;

/// The result of a fallible RDKit call. A Rust addition.
pub type Result<T> = std::result::Result<T, Error>;

/// An error RDKit raised, by the Python exception nanobind would raise.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Error {
    /// Unparsable input. Python returns `None`.
    Parse {
        /// What RDKit's parser reported.
        message: String,
    },
    /// An atom whose valence is not allowed (`AtomValenceException`).
    AtomValence {
        /// The atom's index.
        atom_idx: u32,
        /// RDKit's message.
        message: String,
    },
    /// An aromatic atom outside a ring (`AtomKekulizeException`).
    AtomKekulize {
        /// The atom's index.
        atom_idx: u32,
        /// RDKit's message.
        message: String,
    },
    /// Atoms that could not be kekulized (`KekulizeException`).
    Kekulize {
        /// The atoms' indices.
        atom_indices: Vec<u32>,
        /// RDKit's message.
        message: String,
    },
    /// Another sanitization failure (`MolSanitizeException`).
    MolSanitize {
        /// RDKit's message.
        message: String,
    },
    /// An index out of range (`IndexError`).
    Index {
        /// The index; negative for a Python-style index from the end that
        /// reaches past the start.
        index: i64,
        /// The length of what was indexed.
        len: usize,
    },
    /// A missing key (`KeyError`).
    Key {
        /// The key.
        key: String,
    },
    /// A bad value (`ValueError`), including invalid UTF-8.
    Value {
        /// What was wrong.
        message: String,
    },
    /// Any other failure (`RuntimeError`), including RDKit invariant
    /// violations.
    Runtime {
        /// RDKit's message.
        message: String,
    },
    /// A file that could not be opened (`OSError`).
    Os {
        /// What failed.
        message: String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Index { index, len } => {
                write!(f, "index {index} is out of range for length {len}")
            }
            Error::Key { key } => write!(f, "key not found: {key}"),
            Error::Parse { message }
            | Error::AtomValence { message, .. }
            | Error::AtomKekulize { message, .. }
            | Error::Kekulize { message, .. }
            | Error::MolSanitize { message }
            | Error::Value { message }
            | Error::Runtime { message }
            | Error::Os { message } => f.write_str(message),
        }
    }
}

impl std::error::Error for Error {}

const PREFIX: &str = "\u{1f}rdkit\u{1f}";

impl From<cxx::Exception> for Error {
    fn from(exception: cxx::Exception) -> Error {
        decode(exception.what())
    }
}

/// Decodes `\x1Frdkit\x1F<Variant>\x1F<payload>\x1F<what>`. Anything else is
/// `Runtime` whole.
pub(crate) fn decode(what: &str) -> Error {
    let runtime = || Error::Runtime {
        message: what.to_owned(),
    };
    let Some(rest) = what.strip_prefix(PREFIX) else {
        return runtime();
    };
    let mut parts = rest.splitn(3, '\u{1f}');
    let (Some(variant), Some(payload), Some(message)) = (parts.next(), parts.next(), parts.next())
    else {
        return runtime();
    };
    let message = message.to_owned();
    let number = |text: &str| text.parse::<u32>().ok();
    let error = match variant {
        "AtomValence" => number(payload).map(|atom_idx| Error::AtomValence { atom_idx, message }),
        "AtomKekulize" => number(payload).map(|atom_idx| Error::AtomKekulize { atom_idx, message }),
        "Kekulize" => payload
            .split_whitespace()
            .map(number)
            .collect::<Option<Vec<_>>>()
            .map(|atom_indices| Error::Kekulize {
                atom_indices,
                message,
            }),
        "MolSanitize" => Some(Error::MolSanitize { message }),
        "Index" => payload.split_once(' ').and_then(|(index, len)| {
            index
                .parse()
                .ok()
                .zip(len.parse().ok())
                .map(|(index, len)| Error::Index { index, len })
        }),
        "Key" => Some(Error::Key {
            key: payload.to_owned(),
        }),
        "Value" => Some(Error::Value { message }),
        "Os" => Some(Error::Os { message }),
        "Parse" => Some(Error::Parse { message }),
        "Runtime" => Some(Error::Runtime { message }),
        _ => None,
    };
    error.unwrap_or_else(runtime)
}

impl Error {
    /// An `Index` error for `index` into a sequence of `len` items; an index
    /// past `i64`'s range saturates.
    pub(crate) fn index<I: TryInto<i64>>(index: I, len: usize) -> Error {
        Error::Index {
            index: index.try_into().unwrap_or(i64::MAX),
            len,
        }
    }
}

/// Unwraps the result of a bridge call that cannot fail for any input the safe
/// API admits (CONVENTIONS.md E1). An exception here is a bug in the bindings.
pub(crate) fn expect_ok<T>(result: std::result::Result<T, cxx::Exception>) -> T {
    match result {
        Ok(value) => value,
        Err(exception) => panic!("unexpected RDKit exception: {}", Error::from(exception)),
    }
}

/// Unwraps an operator's result; operators panic where Python raises.
pub(crate) fn op<T>(result: std::result::Result<T, cxx::Exception>) -> T {
    result.unwrap_or_else(|err| panic!("{}", Error::from(err)))
}

/// Converts a C++ string, mapping invalid UTF-8 to `Error::Value`, as Python
/// raises `UnicodeDecodeError`, a `ValueError` (T2).
pub(crate) fn to_string(text: cxx::UniquePtr<cxx::CxxString>) -> Result<String> {
    let text = text.as_ref().expect("RDKit returned a null string");
    text.to_str()
        .map(str::to_owned)
        .map_err(|err| Error::Value {
            message: format!("invalid UTF-8 in a string from RDKit: {err}"),
        })
}

#[doc(hidden)]
pub fn __raise_for_tests(kind: &str) -> Result<()> {
    Ok(ffi::raise_for_tests(kind)?)
}

#[doc(hidden)]
pub fn __raise_unexpected_for_tests() {
    expect_ok(ffi::raise_for_tests("std"));
}
