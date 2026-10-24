// Copyright (c) 2026 Steven Kearnes. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

fn main() {
    println!(
        "cargo::rustc-link-arg=-Wl,-rpath,{}",
        std::env::var("DEP_RDKIT_CXX_LIBDIR").unwrap()
    );
}
