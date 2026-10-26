# RDKit data files

Data files from RDKit's sources (<https://github.com/rdkit/rdkit>, tag
`Release_2026_09_1`), unmodified, at their paths in that tree, under RDKit's
BSD 3-Clause license (`LICENSE`).

The ports of RDKit's pure-Python modules read these where Python reads
`RDConfig.RDDataDir`, so a Rust program needs no RDKit data directory at run
time. The crate embeds them with `include_str!`:

- `Data/FragmentDescriptors.csv`: `rdkit.Chem.Fragments`' patterns
  (`src/chem/fragments.rs`).

The pre-commit hooks skip this directory.
