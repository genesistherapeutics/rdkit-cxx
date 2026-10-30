# RDKit data files

Data files from RDKit's sources (<https://github.com/rdkit/rdkit>, tag
`Release_2025_09_6`), unmodified, at their paths in that tree, under RDKit's
BSD 3-Clause license (`LICENSE`).

The ports of RDKit's Python-only functionality read these where Python reads
`RDConfig.RDDataDir`, so a Rust program needs no RDKit data directory at run
time. The crate embeds them with `include_str!`:

- `Data/BaseFeatures.fdef`: the feature definitions of
  `rdkit.Chem.Pharm2D.DefaultSigFactory` (`src/chem/pharm_2d/init.rs`).
- `Data/FragmentDescriptors.csv`: `rdkit.Chem.Fragments`' patterns
  (`src/chem/fragments.rs`).
- `Data/FunctionalGroups.txt`: the default functional groups of
  `rdkit.Chem.BuildFragmentCatalog.BuildCatalog`
  (`src/chem/build_fragment_catalog.rs`).
- `Data/Functional_Group_Hierarchy.txt`: the default hierarchy of
  `rdkit.Chem.FunctionalGroups.BuildFuncGroupHierarchy`
  (`src/chem/functional_groups.rs`).
- `Data/Salts.txt`: the default salts of `rdkit.Chem.SaltRemover.SaltRemover`
  (`src/chem/salt_remover.rs`).
- `Data/SmartsLib/patty_rules.txt`: `rdkit.Chem.AtomPairs.Sheridan`'s atom
  types (`src/chem/atom_pairs/sheridan.rs`).

The pre-commit hooks skip this directory.
