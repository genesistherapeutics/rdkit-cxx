# Data

RDKit's data files that ported pure-Python modules read by default, from
RDKit's `Data` directory (`RDConfig.RDDataDir` in Python), under RDKit's BSD
license, with trailing whitespace trimmed from their comments. The crate
embeds them, so the ports need no RDKit data directory at run time.

- `Salts.txt`: `SaltRemover`'s default salts.
- `Functional_Group_Hierarchy.txt`: `FunctionalGroups`' default hierarchy,
  which `SimpleEnum.Enumerator` also reads.
- `FunctionalGroups.txt`: `BuildFragmentCatalog.BuildCatalog`'s default
  functional groups.
