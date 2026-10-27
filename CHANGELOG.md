# Changelog

## 2026.3.0

The first release for RDKit 2026.03: the bindings and ports of 2026.9.0, for
RDKit 2026.03.6 (see Compared with 2026.9.0, below). It
supports RDKit 2026.03.6 from conda-forge on Linux (x86-64 and arm64) and
macOS (arm64), and builds with Rust 1.88 or newer. The build refuses other
RDKit releases: the wrapper compiles against copies of headers RDKit does not
install, whose layouts can change between patch releases (`wrapper/vendor/`).
A program built against one release stops when it loads with another.

### Compared with 2026.9.0

- What RDKit 2026.03 lacks is not bound: the embedder's initial-embedding
  modes, embedding force field, `onlyInitialEmbedding`, `ETKDGv4`, and
  Gaussian torsion terms; `GetSymmSSSR`'s `algorithm` and `recalcSSSR`, the
  legacy ring-finding setting, and the attachment-point label helpers;
  `Compute2DCoords`' `usePathAngleExpansion`; `GetSurfacePoints`'
  `allPoints`; `GetAllConformerBestRMSToRef`; and
  `MolCatalogEntry.Serialize` (an entry still pickles through `serialize`).
- Results are RDKit 2026.03.6's: its 2D depictions (and the drawings and mol
  blocks made from them), conformers, bounds matrices, and ring atom orders
  differ from 2026.09's. Parsing does not find ring families.
- RDKit 2026.03 behaviors the bindings keep, as Python has them: the
  `Properties` registry computes `chi2v` and `chi2n` as `chi3v` and `chi3n`
  (RDKit #9643), and `DoubleCubicLatticeVolume.get_surface_points` returns
  the points the surface-area methods tested, once per call.
- The ports follow RDKit 2026.03.6's Python: `BertzCT` takes a non-aromatic
  bond's order as its `BondType`'s number (17 for a dative bond), and
  `fr_Ar_COO`'s description is 2026.03's.
- Checks for what RDKit 2026.09 fixed in RDKit:
  - Every wrapper source reads `RDConfig.h` first, since `RingInfo.h`
    declares its ring-family members only under a macro it does not
    include.
  - `RingInfo.num_relevant_cycles` is `Error::Runtime` before the ring
    families are found, where RDKit 2026.03 converts RingDecomposerLib's
    failure value out of range, and `Error::Value` for more relevant cycles
    than `u32` holds.
  - `MolCatalogEntry::from_pickle` builds the entry as the catalog's reader
    does: RDKit 2026.03's pickle constructor deletes an uninitialized
    pointer.
  - `assign_cip_labels` leaves unlabeled, as RDKit 2026.09 does, a double
    bond whose stereo atom is not a neighbor of its own end and an
    atropisomer bond with an end of no other bonds, where RDKit 2026.03
    reads past an empty vector.
  - A molecule of at most 255 atoms whose substance groups or stereo groups
    name a bond past 255, which RDKit 2026.03 pickles truncated, is refused
    wherever it would be pickled (`Error::Value`; `Mol::to_binary` and
    `ChemicalReaction::to_binary` panic). So are a substance group or stereo
    group list of more than 255 entries, more than 255 stereo groups, and a
    query of more than 255 children, which RDKit 2026.09 pickles truncated
    too: RDKit writes these pickles and cannot read them back.
  - A multithreaded supplier's queues hold at most 2^20 records: RDKit
    2026.03 dereferences null when it cannot allocate one.
  - Drawing refuses an atom of two bonds with a neighbor at its position,
    for which RDKit 2026.03 reads an uninitialized bond type.
  - RDKit 2026.03's substructure matches, embedder, CIP labeler,
    `DetermineBondOrders`, and synthon search leave RDKit's `SIGINT` handler
    installed; the wrappers restore the process's after each call.

### Bindings

- Binds RDKit's C++ API with the coverage of its Python wrappers, as
  conda-forge builds them with Boost.Python: 3,442 of the 4,011 names they
  expose. Each of the other 569 has a reason in `inventory/manifest.toml`:
  Python protocol methods and Boost.Python's container types, 50 names left
  unbound because RDKit's code cannot be made safe from outside, and
  deprecated names.
- Module paths follow Python's in Rust spelling:
  `rdkit.Chem.rdmolfiles.MolFromSmiles` is
  `rdkit_cxx::chem::rdmolfiles::mol_from_smiles`.
- `rdkit.Chem.rdtrajectory`, which only the Boost.Python wrappers build:
  `Trajectory` and `Snapshot` are Rust values (`ForceField`'s
  `minimize_trajectory` returns these snapshots), and
  `read_amber_trajectory` and `read_gromos_trajectory` bind RDKit's readers.
  Snapshots with fewer coordinates than their trajectory's points, which
  Python reads past, are refused, as are point counts that hang RDKit's
  AMBER reader or overflow its array sizes.
- Safe Rust cannot reach undefined behavior in the wrappers. Inputs RDKit
  trusts are checked before RDKit sees them, and process-wide state RDKit
  writes without a lock is locked. Readers of molecule pickles, of the
  formats that carry them (reaction, bundle, catalog, and query pickles, and
  PNG metadata), and of Boost.Serialization archives are `unsafe fn`s whose
  contract is bytes RDKit wrote. One gap remains in RDKit itself; see the
  README's "Known limitations".

### Where the bindings differ from Python

- The reference is conda-forge's build, with Boost.Python wrappers. RDKit's
  nanobind wrappers, which no package ships yet, differ from it in places;
  the crate follows Boost.Python's.
- Boost.Python takes an empty list as `None` for some arguments (`AddHs`'
  `onlyOnAtoms`, `MolFragmentToSmiles`' and `MolFragmentToSmarts`'
  `bondsToUse`, `FragmentOnBonds`' `bondTypes`, the fingerprint
  generators' `countBounds`, `AssignCIPLabels`' atoms and bonds,
  `AlignPoints`' weights). The crate keeps `None` and an empty list apart:
  an empty list of atoms or bonds is none of them.
- `CanonicalRankAtomsInFragment`'s Python binding passes
  `includeChiralPresence` where RDKit takes `includeAtomMaps`, so Python
  ranks by atom maps only with `includeChiralPresence` and never by chiral
  presence. The crate passes each option where RDKit takes it.
- `rdGaussianShape.ShapeInput` takes the C++ constructor's defaults; Python's
  binding requires all four arguments.
- `PointND.Distance`, which Python binds to `Point3D`'s function so that it
  always raises, is not bound.

### Ports of Python-only functionality

Much of what RDKit offers only in Python is ported to Rust on top of the
bindings and computes what Python computes: 1,016 of the 2,606 names RDKit's
Python modules define. The ported modules' Python tests are ported with them,
with expected values from RDKit 2026.03.6.

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
- Some Python functions are broken in RDKit 2026.03.6; the ports implement
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
- `ML.Cluster.Clusters.Cluster` caches its points when it first gains
  children and never updates the cache, so after `AddChild`, `AddChildren`,
  or `RemoveChild` on a tree built by hand, Python's `GetPoints` can return
  stale points; the port recomputes them. Trees Murtagh's clustering builds
  agree.
- `GraphDescriptors.Ipc` and `AvgIpc` (and `Graphs.CharacteristicPolynomial`)
  lose precision above about 80 heavy atoms, in Python and Rust alike, and
  there the two differ by up to about 1%, since NumPy's BLAS sums in another
  order.
- `Quantize.FindVarMultQuantBounds` refuses NaN values, whose places in
  Python's sort decide its result.
- `Randomize.RandomizeMol` returns the empty molecule as it is, where
  Python's binding refuses its empty atom order.
