# rdkit-cxx

Rust bindings for [RDKit](https://github.com/rdkit/rdkit)'s C++ API, with Rust
ports for much of the functionality RDKit currently offers only in Python. The
API mirrors RDKit's Python interface, as conda-forge ships it, in Rust
spelling, so Python code ports line by line:

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

Coverage: every name in RDKit 2025.03's Python interface, in its compiled
modules and its Python modules, has a Rust counterpart or a recorded reason it
has none ([Parity](#parity)). Of the 3,702 names the compiled modules expose,
3,204 are bound. Of the 2,603 its Python modules define, 1,008 are ported; the
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

Each crate version supports one RDKit release series (2025.3.x supports RDKit
2025.03), and within it the releases it has checked: 2025.3.0 supports RDKit
2025.03.6, the series' last. The wrapper compiles against copies of a few
headers RDKit does not install, whose layouts can change between patch
releases, so the build refuses any other release
([`wrapper/vendor/`](wrapper/vendor/README.md)), and a program built against
one stops when it loads with another, as after an update in place. It builds
with Rust 1.88 or newer.

RDKit's version scheme and the crate's differ: RDKit numbers each series by
year and zero-padded month (2025.03), but crate versions are SemVer, which
doesn't allow leading zeros, so they drop the zero (2025.3). The crate's patch
number is its own, not RDKit's.

Install RDKit's development packages from conda-forge and activate the
environment before building:

```sh
mamba create -n rdkit-cxx -c conda-forge "librdkit-dev 2025.03.6" eigen libboost-devel cxx-compiler "cmake>=3.22"
conda activate rdkit-cxx
```

In a checkout of this repository, `mamba env create -f environment.yml` makes
the same environment, with Python and pytest for `inventory/`, and
`cargo test` runs the tests.

The build finds RDKit through `CONDA_PREFIX`. To build against another install
of a supported release, such as one built from RDKit's release tag, set
`RDKIT_PREFIX` to its prefix.

## Using the crate

Choose the crate version that matches the RDKit you build against, which is
not necessarily the newest crate. Each RDKit series has its own line of crate
versions, and each builds only against its own series
([Environment](#environment)). The newest version on crates.io supports only
the newest series, and the lines for earlier series get patch releases of
their own:

| RDKit release | Crate requirement | Branch |
| --- | --- | --- |
| 2026.09.1 | `~2026.9` | `main` |
| 2026.03.6 | `~2026.3` | `rdkit-2026.03` |
| 2025.09.6 | `~2025.9` | `rdkit-2025.09` |
| 2025.03.6 | `~2025.3` | `rdkit-2025.03` |
| 2024.09.6 | `~2024.9` | `rdkit-2024.09` |
| 2024.03.6 | `~2024.3` | `rdkit-2024.03` |

A line binds what its RDKit release has, and its branch's CHANGELOG.md lists
what it lacks compared with the newest line.

Depend on the crate with a tilde requirement, which stays within one RDKit
series:

```toml
[dependencies]
rdkit-cxx = "~2025.3"
```

`cargo add rdkit-cxx` picks the newest version, so name the series:
`cargo add rdkit-cxx@~2025.9`. Cargo's default caret requirement can cross
series: `"2027.3"` would also accept 2027.9, which needs another RDKit.

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
  Params structs gain fields only with a new RDKit series, never in a patch
  release, so that idiom keeps compiling: an option RDKit adds in a patch
  release waits for the next series.
- Enums, and the structs functions return, are `#[non_exhaustive]`: match
  them with a `_` arm, and read their fields rather than building them.
- Callbacks Python takes as any callable take closures: `impl Fn`
  arguments, or `&dyn Fn` fields in params structs.
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
  them (U2). Molecule pickles, and everything that carries them (reaction,
  bundle, catalog, and query pickles, and PNG metadata), can't be checked:
  RDKit's pickle reader trusts the query trees and values it reads. Neither
  can Boost.Serialization archives. So the functions that read them are
  `unsafe fn`s whose contract is "bytes RDKit wrote".
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
  - **Ctrl-C handler:** the embedder, `DetermineBondOrders`, the CIP
    labeler, and the synthon search reset RDKit's `ControlCHandler`, which
    swaps the `SIGINT` handler through static variables and leaves its own
    installed. Every wrapper that can reach one holds `control_c_mutex()`,
    which restores the process's handler when it is released; wrappers that
    match substructures hold it too, as later RDKit releases reset the
    handler from every match.
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

- **Its molecules:** most closures get their own copy of the molecule.
  FMCS's typers and checks get the molecules the caller passed in, borrowed
  as `&Mol`: the wrapper finds each one RDKit calls back with in the
  caller's list. A molecule without ring information, as one from SMARTS
  is, has RDKit's temporary ring information for the search, so ring
  queries in the closure see that.
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

Safe Rust can still reach undefined behavior in RDKit in a few ways. None
corrupts memory on the supported platforms, and each needs a fix in RDKit:

- **RDKit's log flags.** Each logger's enabled flag is a plain `bool` that
  every log statement reads. RDKit writes it inside some calls
  (`MolOps::hasQueryHs`, drawing's `prepareMolForDrawing`, deprotection,
  R-group decomposition, link-node enumeration, the embedder's torsion
  preferences), and `BlockLogs`, `enable_log`, and `disable_log` write it
  too, so two threads using RDKit at once race on it. The wrappers order
  RDKit's own writers with one another (under `control_c_mutex()`), but not
  with logging on other threads: that would mean locking every call.
- **Oversized shifts.** RDKit 2025.03's kekulization shifts a 32-bit one by
  the number of dummy "question" atoms in a ring system (`Kekulize.cpp`),
  undefined from 32; its CIP labeler shifts a 64-bit one by the number of
  descriptors a `PairList` holds, undefined past 64; and FMCS shifts a
  64-bit one by a seed's new external bonds, undefined at exactly 64 (a
  search that would first run 2^63 steps). `mol_from_smiles`,
  `sanitize_mol`, `kekulize`, and `kekulize_if_possible` refuse 32 or more
  such atoms in one ring system, and `assign_cip_labels` refuses more than
  64 stereo atoms and bonds. Not checked are RDKit's other parsers, the
  calls that sanitize inside RDKit, and fused ring systems in which one
  stereo atom is several nodes of a branch of the CIP digraph (a sound bound
  would refuse common polycycles). x86-64 and arm64 take the shift count
  modulo the width, so the result is wrong (or, for kekulization, very
  long), not a memory error.
- **Trajectory readers.** The AMBER and GROMOS trajectory readers add each
  snapshot through `Trajectory::insertSnapshot`, which subtracts an iterator
  its own `insert` may have invalidated (the two operands' order is
  unspecified); the result is an index RDKit ignores.
- **R-group decomposition** with `include_target_mol_in_results` inserts a
  range of a vector into itself when an R group has several fragments
  (`RGroupData::mergeIntoCombinedMol`), which libstdc++ copies from storage
  the insert leaves in place.
Other limits:

- Substructure matching runs one call at a time across the process, under
  `control_c_mutex()`, which later RDKit releases need since their matches
  swap the `SIGINT` handler; matching does not scale across threads.
- Some inputs abort the process by overflowing the stack, which cannot
  corrupt memory: an SLN string with a `{` followed by about 100,000
  characters (RDKit's regular expression), and SMILES of very deep
  molecules on small thread stacks.
- Checks refuse some input Python accepts where RDKit would read it
  unsafely, such as substance-group parent chains that span molecules, or
  InChI options outside structure perception.
- conda-forge's macOS build of RDKit 2025.03.6 has no SLN parser, so on
  macOS the `rd_sln_parse` functions are `Error::Value`.
- Windows is not supported yet.

## Parity

[`inventory/`](inventory) measures coverage of RDKit's compiled Python modules
and of its Python (`.py`) modules, as conda-forge builds them.
`inventory-2025.03.json` lists every class, function, method, property, and
enum value the compiled modules expose, and `python-2025.03.json` what each
Python module defines; `manifest.toml`
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

The reference is conda-forge's `rdkit`, whose compiled modules use
Boost.Python: it is what Python users install. (RDKit can also build nanobind
wrappers, which no package ships yet; they differ from Boost.Python's in
places, and the crate follows Boost.Python's.) The inventories are
regenerated once per RDKit release:

```sh
micromamba create -y -n rdkit-py -c conda-forge rdkit=2025.03.6 python=3.12
site=$(micromamba run -n rdkit-py python -c 'import rdkit, pathlib; print(pathlib.Path(rdkit.__file__).parents[1])')
micromamba run -n rdkit-py python -I inventory/inventory.py --rdkit-path "$site" --output inventory/inventory-2025.03.json
micromamba run -n rdkit-py python -I inventory/inventory.py --rdkit-path "$site" --python --output inventory/python-2025.03.json
```

## Development

[pre-commit](https://pre-commit.com) hooks format the Rust, C++, and Python,
add license headers, and regenerate `tests/manifest_paths.rs`; CI runs them on
every file.

```sh
pip install pre-commit
pre-commit install
```

The workflows pin each action to a commit. Dependabot and the RDKit series
workflow run on `main`, which follows the newest RDKit series.

## Releases

This is the crate's line for RDKit 2025.03, whose last release is 2025.03.6.
`main` follows the newest series; each earlier series is released from a
branch of its own, `rdkit-2025.03` for this one.

Each RDKit series has a long-running branch, `rdkit-<series>`
(`rdkit-2025.09` for RDKit 2025.09), from which its crate versions are
released, and to which fixes for that series go, through pull requests as on
`main`. `main` follows the newest series and releases it; the series gets its
branch when `main` moves on to the next one. Pushing to a branch publishes
nothing; CI runs on every push to `main` and the series' branches.

Releases are cut by hand: run the
[Release workflow](https://github.com/skearnes/rdkit-cxx/actions/workflows/release.yml)
from the Actions tab, from `main` (it always runs `main`'s copy), and choose
the branch to release.

The workflow releases the version in `Cargo.toml` if no tag has it yet (the
first release of an RDKit series), and otherwise bumps the patch level with
[bump-my-version](https://github.com/callowayproject/bump-my-version). It
checks that CI passed on the branch's head and that a series' branch has its
series' version, tests the release commit, builds the packaged crate,
publishes it to crates.io through trusted publishing, pushes the version bump
(if any) to the branch and a `v<version>` tag, and creates a GitHub release,
which it marks as the latest only for the newest version of all. The release's
notes are its section of [CHANGELOG.md](CHANGELOG.md): for a patch release,
the `## Unreleased` section, which the bump renames to the version, and for
the first release of a series, the section for its version
(`## 2027.3.0`). Write them before running it; the workflow stops before it
publishes anything when they are missing. A release that fails partway
finishes when its failed jobs are re-run: each step skips what an earlier
attempt finished, including the upload.

The patch level moves only through the workflow. A new RDKit series sets the
major and minor levels by hand, in the pull request that updates
`environment.yml`:

```sh
bump-my-version bump --new-version 2027.3.0
```

Before that pull request, create the old series' branch from `main`
(`rdkit-2026.09`), which keeps the series and its releases from then on, and
add it to the Release workflow's branch choices.

The [RDKit series workflow](.github/workflows/rdkit-series.yml) checks
conda-forge daily and opens an issue, assigned to the repository owner, with
the steps, when it has a `librdkit-dev` release newer than the one `main`'s
`environment.yml` pins: a new series, or a patch release of `main`'s. The
series' branches need no such check, since RDKit stops patching a series once
the next one is out. A new patch release of the series needs a patch release
of the crate too: the build accepts it only once `wrapper/vendor/`
holds its headers ([`wrapper/vendor/`](wrapper/vendor/README.md)).

## License

BSD-3-Clause, as RDKit. The package also carries files from other projects,
each with its license: RDKit's headers (`wrapper/vendor/`), data files
(`wrapper/third_party/rdkit/Data/`), and test data (`tests/data/`), under
RDKit's BSD license (`wrapper/third_party/rdkit/LICENSE`), and coordgenlibs's
files (`wrapper/third_party/coordgen/`), under its BSD license.
