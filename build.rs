// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Builds the C++ half of the crate: generates cxx's C++ for each bridge, builds
// it with the hand-written wrappers through CMake against an installed RDKit,
// and links the result.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// One name per `src/ffi/<name>.rs` (CONVENTIONS.md B4).
const BRIDGES: &[&str] = &[
    "rd_base",
    "errors",
    "mol",
    "molfiles",
    "png",
    "suppliers",
    "writers",
    "atom",
    "bond",
    "props",
    "geometry",
    "conformer",
    "ring_info",
    "molops",
    "substruct",
    "periodic_table",
    "monomer_info",
    "stereo",
    "substance_group",
    "mol_bundle",
    "query",
    "resonance",
    "data_structs",
    "fpb",
    "grid",
    "perception",
    "chirality",
    "fingerprints",
    "fingerprint_generator",
    "adjust_query",
    "graph",
    "transforms",
    "rdqueries",
    "partial_charges",
    "cip_labeler",
    "mol_hash",
    "mhfp",
    "reduced_graphs",
    "sln_parse",
    "conformer_parser",
    "mmpa",
    "deprotect",
    "descriptors",
    "descriptors_3d",
    "dclv",
    "properties",
    "standardize",
    "standardize_fragment",
    "standardize_charge",
    "standardize_normalize",
];

fn main() {
    for var in ["RDKIT_PREFIX", "CONDA_PREFIX", "CXX", "DOCS_RS"] {
        println!("cargo::rerun-if-env-changed={var}");
    }
    // docs.rs has no RDKit. Documentation needs no C++, so skip building it.
    if env::var_os("DOCS_RS").is_some() {
        return;
    }
    let prefix = rdkit_prefix();
    println!("cargo::rerun-if-changed=wrapper");
    println!("cargo::rerun-if-changed=src/ffi");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let (generated_include, generated_sources) = generate_bridges(&out_dir.join("cxxgen"));

    let sources: Vec<String> = generated_sources
        .iter()
        .map(|path| path.display().to_string())
        .collect();
    let dst = cmake::Config::new("wrapper")
        .define("CMAKE_PREFIX_PATH", &prefix)
        .define("RDKIT_CXX_SUPPORTED_SERIES", supported_series())
        .define("RDKIT_CXX_GENERATED_SOURCES", sources.join(";"))
        .define("RDKIT_CXX_GENERATED_INCLUDE_DIR", &generated_include)
        .profile("Release")
        .build_target("rdkitcxx")
        .build();
    link(&dst.join("build"));
}

/// The RDKit install prefix: `RDKIT_PREFIX`, else `CONDA_PREFIX`.
fn rdkit_prefix() -> PathBuf {
    match env::var_os("RDKIT_PREFIX").or_else(|| env::var_os("CONDA_PREFIX")) {
        Some(prefix) => PathBuf::from(prefix),
        None => {
            // Exit rather than panic, so the message is not buried under a
            // backtrace.
            eprintln!(
                "rdkit-cxx: no RDKit install found. Activate a conda environment with \
                 librdkit-dev, libboost-devel, cxx-compiler, and cmake, or set \
                 RDKIT_PREFIX to an RDKit install prefix."
            );
            std::process::exit(1);
        }
    }
}

/// The RDKit release series this crate version supports: crate 2026.9.x
/// supports RDKit 2026.09.
fn supported_series() -> String {
    let major = env::var("CARGO_PKG_VERSION_MAJOR").unwrap();
    let minor: u32 = env::var("CARGO_PKG_VERSION_MINOR")
        .unwrap()
        .parse()
        .unwrap();
    format!("{major}.{minor:02}")
}

/// Writes cxx's C++ for every bridge under `root`, returning the include
/// directory and the generated sources.
fn generate_bridges(root: &Path) -> (PathBuf, Vec<PathBuf>) {
    let include = root.join("include");
    let header_dir = include.join("rdkit-cxx").join("ffi");
    let source_dir = root.join("src");
    fs::create_dir_all(include.join("rust")).unwrap();
    fs::create_dir_all(&header_dir).unwrap();
    fs::create_dir_all(&source_dir).unwrap();
    write_if_changed(
        &include.join("rust").join("cxx.h"),
        cxx_gen::HEADER.as_bytes(),
    );

    let mut sources = Vec::new();
    for name in BRIDGES {
        let path = Path::new("src").join("ffi").join(format!("{name}.rs"));
        let text = fs::read_to_string(&path)
            .unwrap_or_else(|err| panic!("reading {}: {err}", path.display()));
        let tokens: proc_macro2::TokenStream = text
            .parse()
            .unwrap_or_else(|err| panic!("parsing {}: {err:?}", path.display()));
        let generated = cxx_gen::generate_header_and_cc(tokens, &cxx_gen::Opt::default())
            .unwrap_or_else(|err| panic!("generating C++ for {}: {err}", path.display()));
        write_if_changed(&header_dir.join(format!("{name}.rs.h")), &generated.header);
        let source = source_dir.join(format!("{name}.rs.cc"));
        write_if_changed(&source, &generated.implementation);
        sources.push(source);
    }
    (include, sources)
}

/// Writes `contents` unless the file already holds them, so CMake does not
/// rebuild unchanged sources.
fn write_if_changed(path: &Path, contents: &[u8]) {
    if fs::read(path).is_ok_and(|existing| existing == contents) {
        return;
    }
    fs::write(path, contents).unwrap_or_else(|err| panic!("writing {}: {err}", path.display()));
}

/// Links the wrapper library, the RDKit libraries CMake resolved, and the C++
/// standard library.
fn link(build_dir: &Path) {
    println!("cargo::rustc-link-search=native={}", build_dir.display());
    println!("cargo::rustc-link-lib=static=rdkitcxx");

    let listing = build_dir.join("link-libraries.txt");
    let listing = fs::read_to_string(&listing)
        .unwrap_or_else(|err| panic!("reading {}: {err}", listing.display()));
    let mut lib_dir = None;
    for line in listing.lines().filter(|line| !line.is_empty()) {
        let path = Path::new(line);
        let dir = path.parent().unwrap();
        let file = path.file_name().unwrap().to_str().unwrap();
        println!("cargo::rustc-link-search=native={}", dir.display());
        println!("cargo::rustc-link-lib=dylib={}", library_name(file));
        lib_dir.get_or_insert_with(|| dir.to_path_buf());
    }
    let lib_dir = lib_dir.expect("CMake listed no RDKit libraries");

    let (stdlib, stdlib_file) = if env::var("CARGO_CFG_TARGET_OS").unwrap() == "macos" {
        ("c++", "libc++.dylib")
    } else {
        ("stdc++", "libstdc++.so")
    };
    // Link the standard library the wrappers compiled against, not the
    // system's.
    let cxx = env::var("CXX").unwrap_or_else(|_| "c++".to_string());
    if let Ok(output) = Command::new(&cxx)
        .arg(format!("-print-file-name={stdlib_file}"))
        .output()
    {
        let printed = PathBuf::from(String::from_utf8_lossy(&output.stdout).trim());
        if printed.is_absolute() {
            if let Some(dir) = printed.parent() {
                println!("cargo::rustc-link-search=native={}", dir.display());
            }
        }
    }
    println!("cargo::rustc-link-lib=dylib={stdlib}");

    // The crate's own tests and examples find RDKit through an rpath. Other
    // crates add theirs from `DEP_RDKIT_CXX_LIBDIR` (see the README).
    println!("cargo::rustc-link-arg=-Wl,-rpath,{}", lib_dir.display());
    println!("cargo::metadata=libdir={}", lib_dir.display());
}

/// `libRDKitGraphMol.so.1.2026.09.1` → `RDKitGraphMol`.
fn library_name(file: &str) -> &str {
    let name = file.strip_prefix("lib").unwrap_or(file);
    let end = [".so", ".dylib"]
        .iter()
        .filter_map(|suffix| name.find(suffix))
        .min()
        .unwrap_or(name.len());
    &name[..end]
}
