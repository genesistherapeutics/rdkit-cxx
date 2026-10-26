# Changelog

## 2026.9.0

The first release. It supports RDKit 2026.09.1 from conda-forge on Linux
(x86-64 and arm64) and macOS (arm64), and builds with Rust 1.88 or newer. The
build refuses other RDKit releases: the wrapper compiles against copies of
headers RDKit does not install, whose layouts can change between patch
releases (`wrapper/vendor/`).

### Bindings

- Binds RDKit's C++ API with the coverage of its nanobind Python wrappers:
  3,478 of the 3,910 names they expose. Each of the other 432 has a reason
  in `inventory/manifest.toml`: Python protocol methods and nanobind's
  container types, 52 names left unbound because RDKit's code cannot be made
  safe from outside, and deprecated names.
- Module paths follow Python's in Rust spelling:
  `rdkit.Chem.rdmolfiles.MolFromSmiles` is
  `rdkit_cxx::chem::rdmolfiles::mol_from_smiles`.
- Safe Rust cannot reach undefined behavior in the wrappers. Inputs RDKit
  trusts are checked before RDKit sees them, process-wide state RDKit writes
  without a lock is locked, and readers of Boost.Serialization archives are
  `unsafe fn`s. Two gaps remain in RDKit itself; see the README's
  "Known limitations".

### Ports of Python-only functionality

Much of what RDKit offers only in Python is ported to Rust on top of the
bindings and computes what Python computes: 1,017 of the 2,605 names RDKit's
Python modules define. The ported modules' Python tests are ported with them,
with expected values from RDKit 2026.09.1.

- **Descriptors:** `Descriptors` (`desc_list` in Python's order,
  `calc_mol_descriptors`), `Descriptors3D`, `Crippen`, `MolSurf`,
  `GraphDescriptors`, `EState`, `Fragments`, `Lipinski`, `QED`,
  `ReducedGraphs`, `SpacialScore`, `ChemUtils.DescriptorUtilities`, and
  `ML.Descriptors.MoleculeDescriptors`. `calc_mol_descriptors` matches
  Python's 217 descriptors on a set of 40 molecules to 1e-12.
- **Fragments and enumeration:** `BRICS` (with a lazy `brics_build` in
  Python's seeded order), `Recap`, `Fraggle`, `Scaffolds.MurckoScaffold`,
  `SaltRemover`, `FunctionalGroups`, `FragmentMatcher`, `FragmentCatalog`,
  `BuildFragmentCatalog`, `EnumerateHeterocycles`, `EnumerateStereoisomers`,
  and `SimpleEnum`.
- **Fingerprints and pharmacophores:** `AtomPairs`, `MACCSkeys`,
  `TorsionFingerprints`, `Pharm2D`, `Pharm3D`, `FeatMaps`, `Features`,
  `FeatFinderCLI`, and `Subshape`.
- **Drawing:** `Draw` (images as PNG bytes, grids, reactions, files) and
  `Draw.SimilarityMaps`.
- **Molecule helpers:** the helpers of `rdkit.Chem` and `AllChem`
  (`canon_smiles`, `constrained_embed`, `enumerate_library_from_reaction`,
  and others), `RegistrationHash`, `MolKey` (except its Avalon-checker
  functions), `Randomize`, `TemplateAlign`, `SATIS`, `ChemicalFeatures`,
  `Graphs`, and `ChemUtils`' `AlignDepict`, `SDFToCSV`, and `TemplateExpand`.
- **Machine learning and data structures:** `ML.Cluster` (Butina
  clustering, cluster trees, resemblance metrics), `ML.Data`, `ML.InfoTheory`,
  `ML.Scoring`, `ML.SLT`, `DataStructs`' Python modules (`TopNContainer`,
  `VectCollection`, `BitEnsemble`, and others), `SimDivFilters`, `RDLogger`,
  `RDPaths`, and `RDRandom` (Python's Mersenne Twister, draw for draw).
- The data files the ports read (`FragmentDescriptors.csv`, `Salts.txt`,
  the functional-group hierarchies, `BaseFeatures.fdef`, and others) are
  vendored, so the defaults need no RDKit data directory.

Not ported, each with a recorded reason: Python-ecosystem glue (`PandasTools`,
`IPythonConsole`, `PyMol`, Windows COM), legacy code whose replacement is
bound (`sping`, `MolDrawing`, the MCS search in `fmcs`), deprecated code
(`Chem.Fingerprints`, `PropertyMol`, `FastSDMolSupplier`), command-line
plumbing, and Python protocol hooks.

### Where the ports differ from Python

- Arguments that only steer RDKit's descriptor caches (`force`) are not
  ported, and results are `f64`.
- Functions returning PIL images return PNG bytes, Python callables are Rust
  closures, and generators are iterators or lists.
- Some Python functions are broken in RDKit 2026.09.1; the ports implement
  what they intend:
  - Under NumPy 2 or Python 3.12: `ML.InfoTheory` `BitRank`,
    `ML.Data.Stats.StandardizeMatrix` on 2-D input,
    `ML.Cluster.Resemblance.EuclideanDistance` and `CalcMetricMatrix`,
    `ML.files.ReadDataFile`, and `RDRandom.shuffle` without `random`.
  - `ML.Data.Stats.PrincipalComponents` returns rows of LAPACK's matrix
    rather than eigenvectors; the port returns eigenvectors.
  - `SupplierFromFilename` fails on `.tdt` files.
- `AUTOCORR2D` refuses molecules with dummy atoms, for which RDKit's C++
  reads past the end of its property tables, so the descriptor calculators
  give NaN for it there.
- Where Python's results depend on the platform, as `InfoBitRanker`'s order
  of tied gains does, the tests accept any of them.
