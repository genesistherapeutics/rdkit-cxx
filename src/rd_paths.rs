// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.RDPaths`, which RDKit's build writes: the directories RDKit
//! installs its data in.
//!
//! As Python's `RDConfig`, through which RDKit's Python code reads them,
//! each is under `$RDBASE` when that is set. Otherwise it is under the
//! `share/RDKit` directory of the RDKit the crate was built against
//! (`RDKIT_PREFIX`, else `CONDA_PREFIX`, at build time), where conda
//! installs RDKit's data and where Python's `RDPaths` points.

use std::path::PathBuf;

fn share() -> PathBuf {
    if let Some(base) = std::env::var_os("RDBASE") {
        return PathBuf::from(base);
    }
    // Builds for the documentation alone (`DOCS_RS`) do not record it.
    let prefix = option_env!("RDKIT_CXX_RDKIT_PREFIX").unwrap_or("");
    PathBuf::from(prefix).join("share").join("RDKit")
}

/// RDKit's data directory, `$RDBASE/Data` or
/// `<prefix>/share/RDKit/Data`. Mirrors `RDPaths.RDDataDir`.
pub fn rd_data_dir() -> PathBuf {
    share().join("Data")
}

/// RDKit's documentation directory, `$RDBASE/Docs` or
/// `<prefix>/share/RDKit/Docs`. Mirrors `RDPaths.RDDocsDir`.
pub fn rd_docs_dir() -> PathBuf {
    share().join("Docs")
}

/// RDKit's projects directory, `$RDBASE/Projects` or
/// `<prefix>/share/RDKit/Projects`. Mirrors `RDPaths.RDProjDir`.
pub fn rd_proj_dir() -> PathBuf {
    share().join("Projects")
}

/// RDKit's directory of contributed code, `$RDBASE/Contrib` or
/// `<prefix>/share/RDKit/Contrib`. Mirrors `RDPaths.RDContribDir`.
pub fn rd_contrib_dir() -> PathBuf {
    share().join("Contrib")
}
