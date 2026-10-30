// Copyright (c) 2026 Genesis Molecular AI, Inc.
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
    "depictor",
    "coord_gen",
    "force_field",
    "force_field_helpers",
    "embedder",
    "distance_geometry",
    "mol_transforms",
    "alignment",
    "determine_bonds",
    "shape_helpers",
    "mol_align",
    "chem_reactions",
    "enumerate",
    "standardize",
    "standardize_fragment",
    "standardize_charge",
    "standardize_normalize",
    "standardize_validate",
    "standardize_tautomer",
    "standardize_pipeline",
    "mol_catalog",
    "frag_catalog",
    "filter_catalog",
    "tautomer_query",
    "generalized_substruct",
    "rascal_mces",
    "fmcs",
    "scaffold_network",
    "substruct_library",
    "rgroup",
    "chemical_features",
    "mol_chemical_features",
    "free_sasa",
    "mif",
    "shape_align",
    "metric_matrix_calc",
    "cluster",
    "sim_div_pickers",
    "info_theory",
    "abbreviations",
    "mol_enumerator",
    "mol_processing",
    "mol_interchange",
    "draw",
    "inchi",
    "trajectory",
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
    // rdkit_cxx::rd_paths finds RDKit's data under the prefix, as Python's
    // RDPaths does under the one RDKit was installed to.
    println!("cargo::rustc-env=RDKIT_CXX_RDKIT_PREFIX={}", text(&prefix));
    println!("cargo::rerun-if-changed=wrapper");
    println!("cargo::rerun-if-changed=src/ffi");
    // conda installs packages' files with the times they were built at,
    // which can be older than the last build, so RDKit's files (see link) do
    // not show an update in place. Its installers (conda, mamba, pixi) all
    // write fresh records to conda-meta for every change to the environment,
    // and Cargo watches a directory's contents.
    let conda_meta = prefix.join("conda-meta");
    if conda_meta.is_dir() {
        println!("cargo::rerun-if-changed={}", text(&conda_meta));
    }
    // Clippy checks the Rust and links nothing, so it needs no C++. Cargo
    // runs this script for clippy apart from its run for `cargo build`
    // (clippy's rustc wrapper is part of the run's identity), so building the
    // C++ here would compile all of it a second time. Downstream build
    // scripts still read the library directory.
    if under_clippy() {
        println!("cargo::metadata=libdir={}", text(&prefix.join("lib")));
        return;
    }

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let (generated_include, generated_sources) = generate_bridges(&out_dir.join("cxxgen"));

    let sources: Vec<&str> = generated_sources
        .iter()
        .map(|path| cmake_list_item(path))
        .collect();
    let dst = cmake::Config::new("wrapper")
        .define("CMAKE_PREFIX_PATH", cmake_list_item(&prefix))
        .define("RDKIT_CXX_SUPPORTED_SERIES", supported_series())
        .define("RDKIT_CXX_GENERATED_SOURCES", sources.join(";"))
        .define("RDKIT_CXX_GENERATED_INCLUDE_DIR", text(&generated_include))
        .profile("Release")
        .build_target("rdkitcxx")
        .build();
    link(&dst.join("build"));
}

/// Whether Cargo runs this script for `cargo clippy`, which makes
/// clippy-driver the workspace's rustc wrapper.
fn under_clippy() -> bool {
    env::var_os("RUSTC_WORKSPACE_WRAPPER").is_some_and(|wrapper| {
        Path::new(&wrapper)
            .file_stem()
            .is_some_and(|stem| stem == "clippy-driver")
    })
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
                 librdkit-dev, libboost-devel, eigen, cxx-compiler, and cmake, or \
                 set RDKIT_PREFIX to an RDKit install prefix."
            );
            std::process::exit(1);
        }
    }
}

/// The RDKit release series this crate version supports: crate 2024.9.x
/// supports RDKit 2024.09.
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
    println!("cargo::rustc-link-search=native={}", text(build_dir));
    println!("cargo::rustc-link-lib=static=rdkitcxx");

    let listing = build_dir.join("link-libraries.txt");
    let listing = fs::read_to_string(&listing)
        .unwrap_or_else(|err| panic!("reading {}: {err}", listing.display()));
    let mut lib_dir = None;
    for line in listing.lines().filter(|line| !line.is_empty()) {
        let path = Path::new(line);
        let dir = path.parent().unwrap();
        let file = path.file_name().unwrap().to_str().unwrap();
        println!("cargo::rustc-link-search=native={}", text(dir));
        println!("cargo::rustc-link-lib=dylib={}", library_name(file));
        lib_dir.get_or_insert_with(|| dir.to_path_buf());
    }
    let lib_dir = lib_dir.expect("CMake listed no RDKit libraries");

    // RDKit's files are outside the package, so Cargo does not watch them.
    // An update in the same prefix replaces them, and the wrapper must then
    // be rebuilt, and its release check rerun, against the new headers. For
    // an install outside conda (whose files keep their packages' times; see
    // main), the files' times show it. Only existing paths: Cargo reruns
    // every build for a path that does not exist.
    let watched = build_dir.join("rdkit-files.txt");
    let watched = fs::read_to_string(&watched)
        .unwrap_or_else(|err| panic!("reading {}: {err}", watched.display()));
    for line in watched.lines().filter(|line| !line.is_empty()) {
        if Path::new(line).exists() {
            println!("cargo::rerun-if-changed={line}");
        }
    }

    let (stdlib, stdlib_file) = if env::var("CARGO_CFG_TARGET_OS").unwrap() == "macos" {
        ("c++", "libc++.dylib")
    } else {
        ("stdc++", "libstdc++.so")
    };
    // Link the standard library the wrappers compiled against, not the
    // system's: ask the compiler CMake used, which a cross build or an
    // environment without `CXX` may not name in `CXX`.
    let cxx = fs::read_to_string(build_dir.join("cxx-compiler.txt"))
        .ok()
        .map(|text| text.trim().to_string())
        .filter(|cxx| !cxx.is_empty())
        .or_else(|| env::var("CXX").ok())
        .unwrap_or_else(|| "c++".to_string());
    if let Ok(output) = Command::new(&cxx)
        .arg(format!("-print-file-name={stdlib_file}"))
        .output()
    {
        let printed = PathBuf::from(String::from_utf8_lossy(&output.stdout).trim());
        if printed.is_absolute()
            && let Some(dir) = printed.parent()
        {
            println!("cargo::rustc-link-search=native={}", text(dir));
        }
    }
    println!("cargo::rustc-link-lib=dylib={stdlib}");

    // The crate's own tests and examples find RDKit through an rpath. Other
    // crates add theirs from `DEP_RDKIT_CXX_LIBDIR` (see the README).
    println!("cargo::rustc-link-arg=-Wl,-rpath,{}", text(&lib_dir));
    println!("cargo::metadata=libdir={}", text(&lib_dir));
}

/// A path as text, for Cargo's directives and CMake's arguments, which
/// take paths as text: `display` would hand them a path that is not UTF-8
/// garbled.
fn text(path: &Path) -> &str {
    path.to_str()
        .unwrap_or_else(|| panic!("rdkit-cxx needs UTF-8 paths, but {} is not", path.display()))
}

/// A path as an element of a CMake list, which `;` separates.
fn cmake_list_item(path: &Path) -> &str {
    let item = text(path);
    assert!(
        !item.contains(';'),
        "CMake would split {item} at its ';'; build rdkit-cxx where its paths have none"
    );
    item
}

/// `libRDKitGraphMol.so.1.2024.09.6` → `RDKitGraphMol`.
fn library_name(file: &str) -> &str {
    let name = file.strip_prefix("lib").unwrap_or(file);
    let end = [".so", ".dylib"]
        .iter()
        .filter_map(|suffix| name.find(suffix))
        .min()
        .unwrap_or(name.len());
    &name[..end]
}
