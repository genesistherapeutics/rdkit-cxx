# rdkit-cxx

Rust bindings for [RDKit](https://github.com/rdkit/rdkit)'s C++ API, with Rust
ports for much of the functionality RDKit currently offers only in Python. The
API mirrors RDKit's nanobind wrappers and its Python modules in Rust spelling,
so Python code ports line by line:

```rust,no_run
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles, mol_to_smiles};

fn main() -> Result<(), rdkit_cxx::Error> {
    let mol = mol_from_smiles("Oc1ccccc1O")?;
    let hydroxyl = mol_from_smarts("[OX2H]")?;
    println!("{} {:?}", mol_to_smiles(&mol)?, mol.get_substruct_matches(&hydroxyl)?);
    Ok(())
}
```

The crate binds RDKit's C++ and reimplements none of it. It compiles a C++
wrapper layer (with [cxx](https://cxx.rs) and CMake) against RDKit from
conda-forge and links RDKit's shared libraries. Much of what RDKit implements
only in Python (descriptors such as QED, BRICS and Recap, pharmacophore
fingerprints, `Draw`, Butina clustering, and more) is ported to Rust on top of
those bindings, and computes what Python computes: `rdkit.Chem.QED` is
`rdkit_cxx::chem::qed`.

Coverage: every name in RDKit 2026.09's Python interface, in its nanobind
modules and its Python modules, has a Rust counterpart or a recorded reason it
has none ([Parity](#parity)). Of the 3,910 names the nanobind modules expose,
3,478 are bound. Of the 2,605 its Python modules define, 1,017 are ported; the
rest are Python-ecosystem glue (pandas, Jupyter, PyMOL, Windows COM), legacy
code whose replacement is bound (the sping canvases, the Python MCS search),
command-line plumbing, Python protocol hooks, or deprecated.

The patterns every wrapper follows are in [CONVENTIONS.md](CONVENTIONS.md), and
the limits that remain are under [Known limitations](#known-limitations).
[`examples/`](examples) ports the
[RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)'s recipes
(`cargo run --example cookbook_hybridization`). The
[design](https://github.com/skearnes/rdkit-logbook/blob/main/entries/2026-10-24-rdkit-rust-bindings/README.md)
and a [report on the work](https://github.com/skearnes/rdkit-logbook/blob/main/entries/2026-10-26-rdkit-rust-bindings-report/README.md),
including the RDKit bugs it found, are in rdkit-logbook.

## Environment

Each crate version supports one RDKit release series (2026.9.x supports RDKit
2026.09), and within it the releases it has checked: 2026.9.0 supports RDKit
2026.09.1. The wrapper compiles against copies of a few headers RDKit does not
install, whose layouts can change between patch releases, so the build refuses
any other release ([`wrapper/vendor/`](wrapper/vendor/README.md)), and a
program built against one stops when it loads with another, as after an
update in place. It builds with Rust 1.88 or newer.

Install RDKit's development packages from conda-forge and activate the
environment before building:

```sh
mamba create -n rdkit-cxx -c conda-forge "librdkit-dev 2026.09.1" libboost-devel cxx-compiler "cmake>=3.22"
conda activate rdkit-cxx
```

In a checkout of this repository, `mamba env create -f environment.yml` makes
the same environment, with Python and pytest for `inventory/`, and
`cargo test` runs the tests.

The build finds RDKit through `CONDA_PREFIX`. To build against another install
of a supported release, such as one built from RDKit's release tag, set
`RDKIT_PREFIX` to its prefix.

## Using the crate

Depend on the crate with a tilde requirement, which stays within one RDKit
series:

```toml
[dependencies]
rdkit-cxx = "~2026.9"
```

Cargo's default caret requirement can cross series: `"2027.3"` would also
accept 2027.9, which needs another RDKit.

RDKit's shared libraries live in the conda environment, and conda does not set
`LD_LIBRARY_PATH`. A binary that depends on this crate finds them through an
rpath, which its own build script adds:

```rust,ignore
// build.rs
fn main() {
    // Absent where nothing links, as when docs.rs builds the documentation.
    if let Ok(dir) = std::env::var("DEP_RDKIT_CXX_LIBDIR") {
        println!("cargo::rustc-link-arg=-Wl,-rpath,{dir}");
    }
}
```

This crate cannot add the rpath for you: a dependency's build script reaches
only its own package's binaries and tests. Cargo gives `DEP_RDKIT_CXX_LIBDIR`
only to packages that depend on rdkit-cxx directly, so a package whose
binaries or tests use it through another library also lists rdkit-cxx in its
own `[dependencies]`. (conda-forge's Linux compilers add an rpath to
`$CONDA_PREFIX/lib` on their own, so without the line a binary may run inside
the environment and fail when linked any other way.) See
[`downstream/`](downstream) for a complete example.

The crate pins cxx exactly (`cxx = "=1.0.205"`), because the C++ its build
generates must match the cxx runtime, so a project that also uses cxx needs
the same version.

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

## How the wrappers work

[CONVENTIONS.md](CONVENTIONS.md) holds the full rules. Its rule IDs are cited
below in parentheses.

### Layers

Each RDKit area, such as `substruct` or `chem_reactions`, is bound in three
layers:

| Layer | Path | Role |
| --- | --- | --- |
| Public API | `src/chem/…` | Safe Rust that mirrors the Python module: params structs, iterators, views, `Result`s. |
| Bridge | `src/ffi/<area>.rs` | One `#[cxx::bridge]` per area. It declares opaque C++ types and the wrapper functions, all returning `Result` (B2, B3). |
| C++ wrapper | `wrapper/include/rdkit-cxx/<area>.h`, `wrapper/src/<area>.cpp` | Free functions in namespace `rdkit_cxx` that check their inputs, take locks, call RDKit, and copy results out (B1). |

Bridges never bind RDKit's members directly. A C++ wrapper stands between
every bridge function and RDKit, so each RDKit call has a place for its
checks.

### Building

`build.rs` runs `cxx-gen` on every bridge listed in `BRIDGES`. Then CMake
(`wrapper/CMakeLists.txt`) compiles the generated C++ and the hand-written
wrappers as C++20 into one static library, `rdkitcxx`.

- CMake finds RDKit with `find_package(rdkit CONFIG)`, under `RDKIT_PREFIX`
  or `CONDA_PREFIX`.
- It refuses an RDKit from a release series other than the crate's.
- It links the RDKit libraries named in `RDKIT_CXX_RDKIT_TARGETS`.
- With `DOCS_RS` set, the build script skips the C++, so the docs build
  without RDKit.

### One call, end to end

Take `mol.get_substruct_matches(&query)`:

1. **Rust:** the public method calls its bridge function with the C++
   molecules the two `Mol`s own.
2. **C++:** the wrapper checks whatever RDKit would trust, takes the locks the
   call needs, calls `RDKit::SubstructMatch`, and copies the matches into a
   `rust::Vec`.
3. **Errors:** if anything throws, the `rust::behavior::trycatch` in
   `rdkit-cxx/errors.h` catches it. Every bridge includes that header first
   (E5). It encodes the exception's type and payload, such as an atom index
   or a bad index and length, into the message cxx carries across.
4. **Back in Rust:** `src/error.rs` decodes the message into an `Error`
   variant, most derived type first (E4):
   - `AtomValence { atom_idx, .. }` and the other sanitization variants;
   - `Index`, `Key`, `Value`, `Os`, `Parse`;
   - `Runtime` for RDKit invariants and anything else.

   A C++ exception never unwinds into Rust.

### Ownership

- **Molecules:** every `Mol` owns a C++ `RWMol` through a `UniquePtr` (W1).
  Functions that make molecules copy them out, so every molecule they return
  is a Rust-owned copy (W6).
- **Views:** atoms, bonds, and conformers are views that borrow their
  molecule. `XRef<'m>` is shared; `XMut<'m>` is exclusive and holds an index,
  so two mutable views can never alias (W2, W5).
- **Borrowed containers:** when RDKit wants a container of molecules, such as
  a `MolBundle`, a reaction's templates, or a list of fragments, the wrapper
  builds one per call from `shared_ptr`s with no-op deleters. These borrow
  the Rust molecules, and the container is dropped before the call returns
  (W8).
- **Lifetimes:** C++ objects that keep raw pointers into others carry the
  matching Rust lifetimes. For example, a `MolChemicalFeature<'a>` borrows
  both its molecule and its factory.
- **Shared objects:** objects RDKit shares through `shared_ptr`s between
  owners Rust can't see (the filter catalog's entries and matchers) are
  copied before they change (W9).

### Input RDKit trusts

The rule is that no input reachable from safe Rust reaches undefined
behavior in RDKit, even where Python's bindings allow it (U1). Where RDKit
trusts an argument, the C++ wrapper checks it first and returns an `Error`,
with a comment naming what RDKit would do without the check. Examples: an
index, a size it computes in `int`, a value it casts, a pickle's lengths, a
file's line count.

- **Binary formats:** RDKit's own binary formats (bit-vector, grid, and
  catalog pickles, FPB files, PNG chunks) are checked before RDKit reads
  them (U2). Boost.Serialization archives can't be checked, so the functions
  that read them are `unsafe fn`s whose contract is "bytes RDKit wrote".
- **Cached data:** RDKit caches distance matrices and descriptor results in
  a molecule's properties and reads stale ones past their end. Wrappers that
  compute them clear them before and after the call, so no molecule Rust
  holds carries one (U6, U8).
- **Settings in environment variables:** RDKit keeps some settings in
  environment variables, so their setters are `unsafe fn`s with
  `std::env::set_var`'s contract (U5).
- **Unbound names:** some RDKit functions can't be made safe from outside,
  such as a parser that corrupts memory on malformed input, or code that
  reads files after closing them. These are left unbound. Their entries in
  `inventory/manifest.toml` give the reason, and the crash it was confirmed
  by in the Python build.

### Threads and RDKit's global state

- **Molecules:** `Mol` is `Send` and not `Sync` (S1). The crate needs an
  RDKit built with `RDK_BUILD_THREADSAFE_SSS`, as conda-forge's is (S2).
- **Global locks:** RDKit writes some process-wide state without a lock.
  Every wrapper that reaches such state holds a mutex for the whole RDKit
  call (U4):
  - **Random number generator:** `random_mutex()`.
  - **2D depiction settings and CoordGen templates:** `lock_depiction()`.
  - **Ctrl-C handler:** every substructure match constructs a
    `ControlCHandler`, which swaps the `SIGINT` handler through static
    variables. Every wrapper that can match therefore holds
    `control_c_mutex()`.
  - **Single modules:** some modules' statics have their own locks, such as
    BRICS, R-group decomposition, the shape modules, and Murtagh's
    clustering.

  Wrappers that need several locks take them in one fixed order, so they
  can't deadlock.
- **RDKit's own threads:** a locked call can't order the threads RDKit would
  start for it. Such work runs on the calling thread, or on the crate's own
  threads, which take the lock in turn.
- **Race checks:** CI's `races` job runs `tests/threads.rs` under Valgrind's
  DRD. Each of its phases runs a wrapper that writes the guarded state and a
  plain wrapper holding the same lock, on threads at once. A missing lock
  leaves their accesses unordered, and the job fails.

### Rust closures called from C++

Some RDKit APIs take callbacks: filter matchers, property functors, FMCS
typers, tautomer scores, and validation methods. A Rust closure crosses as
an `extern "Rust"` type that the C++ subclass calls.

- **Its own molecule:** the closure gets its own copy of the molecule.
- **Panics:** a panic is caught at the boundary with `catch_unwind` and
  becomes an exception in C++. Where RDKit lets exceptions through, it
  reaches the caller as an `Error`. Where RDKit would not survive one, as in
  FMCS's search, the wrapper records it, stops the search, and returns it
  when RDKit is done.
- **Locks:** RDKit calls most callbacks between matches, so the call's
  `HoldControlC` releases the lock while the closure runs
  (`ReleaseControlC`). The closure can then call other wrappers that match.
  A closure RDKit calls inside a match runs with the lock held, and a wrapper
  it calls that needs the lock fails instead of deadlocking.

## Known limitations

Safe Rust can still reach undefined behavior in RDKit in two ways. Neither
corrupts memory on the supported platforms, and both need fixes in RDKit:

- **RDKit's log flags.** Each logger's enabled flag is a plain `bool` that
  every log statement reads. RDKit writes it inside some calls
  (`MolOps::hasQueryHs`, drawing's `prepareMolForDrawing`, deprotection,
  R-group decomposition, link-node enumeration, the embedder's torsion
  preferences), and `BlockLogs`, `enable_log`, and `disable_log` write it
  too, so two threads using RDKit at once race on it. The wrappers order
  RDKit's own writers with one another (under `control_c_mutex()`), but not
  with logging on other threads: that would mean locking every call.
- **A cast in molecule pickles.** `MolPickler` converts a float it reads to
  an `int` without a range check when it computes an atom's isotope, which
  is undefined for a value that does not fit. Pickles RDKit wrote never
  hold one; a crafted pickle can.

Other limits:

- Substructure matching runs one call at a time across the process, since
  every match swaps RDKit's `SIGINT` handler through static variables
  (`control_c_mutex()`), so matching does not scale across threads.
- Some inputs abort the process by overflowing the stack, which cannot
  corrupt memory: an SLN string with a `{` followed by about 100,000
  characters (RDKit's regular expression), and SMILES of very deep
  molecules on small thread stacks.
- Checks refuse some input Python accepts where RDKit would read it
  unsafely, such as substance-group parent chains that span molecules, or
  InChI options outside structure perception.
- Windows is not supported yet.

## Parity

[`inventory/`](inventory) measures coverage of the nanobind modules and of
RDKit's Python modules. `inventory-2026.09.json` lists every class,
function, method, property, and enum value the nanobind modules expose, and
`python-2026.09.json` what each Python module defines; `manifest.toml`
maps each to Rust paths or to the reason it has none, given for the name or,
for Python-ecosystem glue and legacy code, for its whole module. CI runs:

```sh
python3 -m pytest inventory
python3 inventory/check.py --strict --write-paths tests/manifest_paths.rs
```

`check.py` fails on misnamed or malformed entries, on unmapped names in
modules marked complete, and, with `--strict`, on any module neither complete
nor covered by a reason. It writes a Rust test that names every mapped path,
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
python3 inventory/inventory.py --rdkit-path rdkit-release --python --output inventory/python-2026.09.json
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

The workflows pin each action to a commit. Dependabot proposes updates to the
actions, the Cargo dependencies, and the pre-commit hooks every week, a week
after each release. A new RDKit series opens an issue instead, since moving to
one takes more than its pin ([Releases](#releases)).

## Releases

Releases are cut by hand: run the
[Release workflow](https://github.com/skearnes/rdkit-cxx/actions/workflows/release.yml)
on `main` from the Actions tab. Merging to `main` publishes nothing.

The workflow releases the version in `Cargo.toml` if no tag has it yet (the
first release of an RDKit series), and otherwise bumps the patch level with
[bump-my-version](https://github.com/callowayproject/bump-my-version). It
checks that CI passed on the release commit, tests it, builds the packaged
crate, publishes it to crates.io through trusted publishing, pushes the version
bump and a `v<version>` tag to `main`, and creates a GitHub release with
generated notes. Add the release's notable changes to
[CHANGELOG.md](CHANGELOG.md) before running it. A release that fails partway
finishes when its failed jobs are re-run: each step skips what an earlier
attempt finished, including the upload.

The patch level moves only through the workflow. A new RDKit series sets the
major and minor levels by hand, in the pull request that updates
`environment.yml`:

```sh
bump-my-version bump --new-version 2027.3.0
```

The [RDKit series workflow](.github/workflows/rdkit-series.yml) checks
conda-forge daily and opens an issue, assigned to the repository owner, with
the steps, when it has a `librdkit-dev` release newer than the one
`environment.yml` pins. A new patch release of the series needs a patch
release of the crate too: the build accepts it only once `wrapper/vendor/`
holds its headers ([`wrapper/vendor/`](wrapper/vendor/README.md)).

## License

BSD-3-Clause, as RDKit. The package also carries files from other projects,
each with its license: RDKit's headers (`wrapper/vendor/`), data files
(`wrapper/third_party/rdkit/Data/`), and test data (`tests/data/`), under
RDKit's BSD license (`wrapper/third_party/rdkit/LICENSE`), and coordgenlibs's
files (`wrapper/third_party/coordgen/`), under its BSD license.
