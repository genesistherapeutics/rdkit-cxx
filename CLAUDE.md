# Working on rdkit-cxx

Rust bindings for RDKit's C++ API: cxx bridges (`src/ffi`), a C++ wrapper
built by CMake (`wrapper/src`, `wrapper/include`), and Rust ports of RDKit's
Python-only code (`src/chem`). Read [CONVENTIONS.md](CONVENTIONS.md) before
changing wrapper code; its rules bind every change, above all U1: no input
reachable from safe Rust may reach undefined behavior. The README covers the
environment, the design, parity, and releases.

## Branches

- Each crate line binds one RDKit release. `main` binds the newest series
  (the release `environment.yml` pins). Each earlier series has a
  long-running branch, `rdkit-<series>` (`rdkit-2025.09` binds RDKit
  2025.09.6), whose API, ports, tests, and CHANGELOG follow that release.
  The newest series gets its branch from `main` when `main` moves on to the
  next series (README, Releases).
- Every change, to `main` or a series' branch, goes through a pull request,
  squash-merged; the `no-commit-to-branch` hook refuses commits on them.
- Never merge `main` into a series' branch: that brings a newer RDKit's API
  and vendored headers. Port a fix by cherry-picking (`git cherry-pick -x`)
  and adapting it to the branch's release, or by writing it again.
- A fix belongs on every line whose RDKit release has the defect. Before
  porting one, check that release's code (below), and say in the commit why
  a line is left out.
- The owner runs the Release workflow, tags, and publishes. Don't, unless
  asked.

## Environment

Build in a conda environment made from `environment.yml`. In a cloud session,
micromamba needs no installer:

```sh
mkdir -p ~/mm && curl -sSL https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.10.0-0.tar.bz2 | tar -xj -C ~/mm
export MAMBA_ROOT_PREFIX=~/mmroot
~/mm/bin/micromamba env create -y -f environment.yml -n rdkit-cxx
~/mm/bin/micromamba run -n rdkit-cxx cargo test --test <name>
```

Expected values in tests come from RDKit's Python for the same release:
`micromamba create -y -n rdkit-py -c conda-forge rdkit=<release> python=3.12`.
Before trusting a script that generates expected values (drawing hashes,
descriptor tables, pickles), check that it reproduces the values already in
the tests.

RDKit's sources, for checking what a release does: a blobless clone
(`git clone --filter=blob:none --no-checkout https://github.com/rdkit/rdkit`),
read with `git show Release_2025_09_6:Code/GraphMol/MolPickler.cpp`; find a
fix's first release with `git tag --contains <commit>`.

## Before pushing

Run what CI runs:

```sh
cargo fmt --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
for e in examples/*.rs; do cargo run -q --example "$(basename "$e" .rs)" > /dev/null || echo "FAIL $e"; done
cmake -P tests/cmake/version_check.cmake
SKIP=no-commit-to-branch pre-commit run --all-files
python3 -m pytest inventory
python3 inventory/check.py --strict --write-paths tests/manifest_paths.rs
env -u CONDA_PREFIX -u RDKIT_PREFIX DOCS_RS=1 cargo doc --no-deps -p rdkit-cxx
```

- A test binary that crashes prints no `test result` line. Check cargo's
  exit status and grep for "process didn't exit successfully".
  `tests/threads.rs` has `harness = false` and prints no summary at all.
- CI also runs the multithreaded tests under Valgrind's DRD (the `races`
  job). Known false positives are suppressed in `tests/drd.supp`.
- Disk is limited: delete stale test binaries (older hashes of the same
  test) in `target/debug/deps` when it runs low.

## Pitfalls

- **Arguments RDKit's releases disagree on.** A call that passes arguments
  positionally still compiles against a release whose parameter list
  differs. On RDKit 2025.09, which has no `canonical` parameter, Kekulize's
  `canonical` bool became `maxBackTracks`. When binding or porting a call,
  read that release's declaration.
- **cxx shared enums** must match the C++ enum's underlying type
  (`#[repr(u32)]` for a plain enum).
- **`RDGeneral/RDConfig.h` macros** such as `RDK_USE_URF` change class
  layouts, and some releases' headers test them without including it. The
  series' branches from 2026.03 down force-include it
  (`wrapper/CMakeLists.txt`), so the wrapper sees the layout RDKit was built
  with.
- **conda-forge's packages differ by platform.** Its osx-arm64 builds of
  RDKit 2026.03 and earlier have no SLN parser, which those branches make
  optional (`RDKIT_CXX_HAVE_SLN_PARSE`), and RDKit 2024.03 has no
  linux-aarch64 build, so that branch's CI has no arm64 job. When moving a
  pin, diff the `lib/libRDKit*` and `include/rdkit/` file lists of the
  linux-64, linux-aarch64, and osx-arm64 `librdkit` and `librdkit-dev`
  packages.
- **Vendored headers** (`wrapper/vendor/<release>/`) are RDKit's own copies of
  headers conda-forge doesn't install. The build accepts only the releases
  that have them (`wrapper/vendor/README.md`).

## Changelog and releases

- Add user-facing changes to CHANGELOG.md under `## Unreleased`. A patch
  release's bump renames that section to the version, and its notes are
  that section.
- On a series' branch, the CHANGELOG has only that line's sections. The first
  one says what the line lacks compared with `main`'s release
  ("Compared with 2026.9.0").
