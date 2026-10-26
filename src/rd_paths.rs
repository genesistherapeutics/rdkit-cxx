// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.RDPaths`, which RDKit's build writes: the directories RDKit
//! installs its data in.
//!
//! As Python's are, they are fixed when the crate is built, under the
//! install prefix of the RDKit it builds against (`RDKIT_PREFIX`, else
//! `CONDA_PREFIX`), in its `share/RDKit` directory. Python's `RDConfig`
//! uses `$RDBASE` instead when it is set; these do not.

use std::path::PathBuf;

fn share() -> PathBuf {
    // Builds for the documentation alone (`DOCS_RS`) do not record it.
    let Some(prefix) = option_env!("RDKIT_CXX_RDKIT_PREFIX") else {
        panic!("the build script records RDKit's install prefix");
    };
    PathBuf::from(prefix).join("share").join("RDKit")
}

/// RDKit's data directory, `<prefix>/share/RDKit/Data`. Mirrors
/// `RDPaths.RDDataDir`.
pub fn rd_data_dir() -> PathBuf {
    share().join("Data")
}

/// RDKit's documentation directory, `<prefix>/share/RDKit/Docs`. Mirrors
/// `RDPaths.RDDocsDir`.
pub fn rd_docs_dir() -> PathBuf {
    share().join("Docs")
}

/// RDKit's projects directory, `<prefix>/share/RDKit/Projects`. Mirrors
/// `RDPaths.RDProjDir`.
pub fn rd_proj_dir() -> PathBuf {
    share().join("Projects")
}

/// RDKit's directory of contributed code, `<prefix>/share/RDKit/Contrib`.
/// Mirrors `RDPaths.RDContribDir`.
pub fn rd_contrib_dir() -> PathBuf {
    share().join("Contrib")
}
