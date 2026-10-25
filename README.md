# rdkit-cxx

Rust bindings for [RDKit](https://github.com/rdkit/rdkit)'s C++ API. The API
mirrors RDKit's nanobind Python wrappers in Rust spelling, so Python code ports
line by line:

```rust,no_run
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles, mol_to_smiles};

let mol = mol_from_smiles("Oc1ccccc1O")?;
let hydroxyl = mol_from_smarts("[OX2H]")?;
println!("{} {:?}", mol_to_smiles(&mol)?, mol.get_substruct_matches(&hydroxyl)?);
# Ok::<(), rdkit_cxx::Error>(())
```

The crate binds RDKit; it reimplements nothing. It compiles a C++ wrapper layer
(with [cxx](https://cxx.rs) and CMake) against RDKit from conda-forge and links
RDKit's shared libraries.

Status: under development. The core slice (molecules, atoms, bonds,
conformers, ring info, properties, SMILES, SMARTS, and molfile I/O,
sanitization, and substructure matching) is bound, and its patterns are
recorded in [CONVENTIONS.md](CONVENTIONS.md). Nothing is released before
parity with the nanobind wrappers. The design is in
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

### From Python to Rust

- Module paths follow Python's: `rdkit.Chem.rdmolfiles.MolFromSmiles` is
  `rdkit_cxx::chem::rdmolfiles::mol_from_smiles`, also re-exported from
  `rdkit_cxx::chem`.
- Optional arguments become `_with_` functions: `GetNumAtoms(onlyExplicit=False)`
  is `get_num_atoms_with_only_explicit(false)`, and calls with several optional
  arguments take a params struct with `Default`:
  `mol_to_smiles_with_params(&mol, &SmilesWriteParams { do_kekule: true, ..Default::default() })`.
- Exceptions become `Result<T, rdkit_cxx::Error>`. A parse that returns `None`
  in Python returns `Error::Parse`, and sanitization errors keep their atom
  indices (`Error::AtomValence { atom_idx, .. }`).
- Atoms, bonds, and conformers reached through a molecule are views that
  borrow it (`AtomRef`, `AtomMut`). Editing goes through `&mut Mol` or
  `RWMol`, one atom at a time, by index.
- `Mol` is `Send` but not `Sync`: RDKit updates cached data inside read-only
  molecules, so parallel code gives each thread its own molecule.

## Parity

[`inventory/`](inventory) measures coverage of the nanobind modules.
`inventory-2026.09.json` lists every class, function, method, property, and
enum value they expose; `manifest.toml` maps each to Rust paths or to the
reason it has none. CI runs:

```sh
python3 -m pytest inventory
python3 inventory/check.py --write-paths tests/manifest_paths.rs
```

`check.py` fails on misnamed or malformed entries and on unmapped names in
modules marked complete, and writes a Rust test that names every mapped path,
so `cargo test` fails if one does not exist.

conda-forge's `rdkit` is built with Boost.Python, so the inventory comes from
an RDKit build with nanobind wrappers, regenerated once per RDKit release:

```sh
git clone --branch Release_2026_09_1 https://github.com/rdkit/rdkit rdkit-release
cmake -S rdkit-release -B rdkit-release/build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DRDK_BUILD_NANOBIND_WRAPPERS=ON -DRDK_BUILD_BOOST_PYTHON_WRAPPERS=OFF \
  -DRDK_INSTALL_INTREE=ON -DRDK_BUILD_CPP_TESTS=OFF \
  -DRDK_BUILD_AVALON_SUPPORT=ON -DRDK_BUILD_CAIRO_SUPPORT=ON -DRDK_BUILD_CONTRIB=ON \
  -DRDK_BUILD_FREESASA_SUPPORT=ON -DRDK_BUILD_INCHI_SUPPORT=ON \
  -DRDK_BUILD_XYZ2MOL_SUPPORT=ON -DRDK_BUILD_YAEHMOP_SUPPORT=ON
cmake --build rdkit-release/build --target install
python3 inventory/inventory.py --rdkit-path rdkit-release --output inventory/inventory-2026.09.json
```

The build needs RDKit's dependencies and `nanobind` 2.x (conda-forge:
`nanobind<3 numpy libboost-devel eigen cairo freetype cxx-compiler cmake
ninja`). The optional modules match conda-forge's build.

## Development

[pre-commit](https://pre-commit.com) hooks format the Rust, C++, and Python,
add license headers, and regenerate `tests/manifest_paths.rs`; CI runs them on
every file.

```sh
pip install pre-commit
pre-commit install
```

### Releases

Run the Release workflow from the Actions tab, on main. It publishes the
version in `Cargo.toml` if that version has no tag yet and otherwise bumps the
patch number first, so moving to a new RDKit series is a pull request that sets
the version to the series' `.0` (2027.3.0 for RDKit 2027.03). It pushes the
release commit to main before publishing, then pushes a `v` tag and creates a
GitHub release, so a failed publish leaves an untagged version that a rerun
ships.

The workflow publishes through crates.io trusted publishing for `release.yml`
and pushes with the `ACTIONS_PUSH` secret, a token whose owner can bypass the
main branch ruleset.

## License

BSD-3-Clause, as RDKit.
