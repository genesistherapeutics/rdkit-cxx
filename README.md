# rdkit-cxx

Rust bindings for [RDKit](https://github.com/rdkit/rdkit)'s C++ API. The API
mirrors RDKit's nanobind Python wrappers in Rust spelling, so Python code ports
line by line:

```rust,no_run
use rdkit_cxx::rd_base;

println!("RDKit {}", rd_base::rdkit_version());
```

The crate binds RDKit; it reimplements nothing. It compiles a small C++ wrapper
layer (with [cxx](https://cxx.rs) and CMake) against RDKit from conda-forge and
links RDKit's shared libraries.

Status: under development. Nothing is released before parity with the nanobind
wrappers. The design is in
[rdkit-logbook](https://github.com/skearnes/rdkit-logbook/blob/main/entries/2026-10-24-rdkit-rust-bindings/README.md).

## Environment

Each crate version supports one RDKit release series: 2026.9.x supports RDKit
2026.09. Install RDKit's development packages from conda-forge and activate the
environment before building:

```sh
mamba env create -f environment.yml
conda activate rdkit-cxx
cargo test
```

The build finds RDKit through `CONDA_PREFIX`. To build against another RDKit
install, such as a development build, set `RDKIT_PREFIX` to its prefix.

## Using the crate

RDKit's shared libraries live in the conda environment, and conda does not set
`LD_LIBRARY_PATH`. A binary that depends on this crate finds them through an
rpath, which its own build script adds with one line:

```rust,ignore
// build.rs
fn main() {
    println!("cargo::rustc-link-arg=-Wl,-rpath,{}", std::env::var("DEP_RDKIT_CXX_LIBDIR").unwrap());
}
```

This crate cannot add the rpath for you: a dependency's build script reaches
only its own package's binaries and tests. (conda-forge's Linux compilers add
an rpath to `$CONDA_PREFIX/lib` on their own, so without the line a binary may
run inside the environment and fail when linked any other way.) See
[`downstream/`](downstream) for a complete example.

## License

BSD-3-Clause, as RDKit.
