// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

fn main() {
    // Absent where nothing links, as when docs.rs builds the documentation.
    if let Ok(dir) = std::env::var("DEP_RDKIT_CXX_LIBDIR") {
        println!("cargo::rustc-link-arg=-Wl,-rpath,{dir}");
    }
}
