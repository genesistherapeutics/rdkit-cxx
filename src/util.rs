// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Helpers the modules share.

use std::path::Path;

/// A file path as the bridge carries it (CONVENTIONS.md A8).
pub(crate) fn path_bytes(path: &Path) -> &[u8] {
    path.as_os_str().as_encoded_bytes()
}

/// A file path from the bytes the bridge carries (A8).
pub(crate) fn path_from_bytes(bytes: Vec<u8>) -> std::path::PathBuf {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStringExt;
        std::path::PathBuf::from(std::ffi::OsString::from_vec(bytes))
    }
    #[cfg(not(unix))]
    {
        std::path::PathBuf::from(String::from_utf8_lossy(&bytes).into_owned())
    }
}
