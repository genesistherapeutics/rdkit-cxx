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

## License

BSD-3-Clause, as RDKit.
