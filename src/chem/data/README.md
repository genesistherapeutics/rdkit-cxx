# Data

From RDKit's source tree, release 2026.09.1, under RDKit's BSD license:

- `FragmentDescriptors.csv`: `Data/FragmentDescriptors.csv`, with trailing
  whitespace trimmed. `rdkit_cxx::chem::fragments` embeds it, as Python's
  `rdkit.Chem.Fragments` reads it from RDKit's data directory.
