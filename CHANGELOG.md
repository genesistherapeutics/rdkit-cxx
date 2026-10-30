# Changelog

## Unreleased

### Changes

- The copyright holder is now Genesis Molecular AI, in `LICENSE` and
  in every source file's header. The license is still BSD-3-Clause.

### Fixes

- `RingInfo.num_relevant_cycles` is `Error::Value` for more relevant cycles
  than `u32` holds. RDKit converts RingDecomposerLib's count, a `double`, to
  `unsigned int` without a check, which is undefined out of range; Python
  returns 33 for a ring of 33 diamonds, which has 2^33 + 33.
- Sanitization (`sanitize_mol`, `mol_from_smiles` and the other SMILES
  parsers, `sanitize_rxn`, `sanitize_rxn_as_mols`, `validate_smiles`),
  `get_symm_sssr`, 2D depiction with RDKit's depictor, and the wrappers that
  symmetrize a molecule's rings themselves (drawing, abbreviations, tautomer
  scores) are `Error::Value` for a molecule with more than `i32::MAX`
  relevant cycles. RDKit's symmetrized SSSR reserves room for them all,
  converting RingDecomposerLib's `double` count to `size_t` unchecked, which
  is undefined from 2^64 (a ring of 64 diamonds, 192 atoms); Python runs out
  of memory sanitizing a ring of 33. Molecules with at most 31 independent
  cycles skip the count, so sanitization costs about 1-4% more. Other paths
  to the conversion remain; see the README's "Known limitations".
- Molecules RDKit's pickler would write but could not read back are
  refused wherever a molecule is pickled: `Error::Value`, or a panic from
  `Mol::to_binary` and `ChemicalReaction::to_binary`, which stay
  infallible since Python's `ToBinary` raises for none of them. MolPickler
  writes some counts in one byte and truncates larger ones: each query
  node's number of children, and, in a molecule of at most 255 atoms and
  255 bonds, each substance group's numbers of atoms, parent atoms, bonds,
  brackets, CStates, and attachment points, and the number of stereo
  groups. The reader then takes the rest of the list (bracket coordinates,
  say) for the pickle's next fields and parses them unchecked, which safe
  Rust reached through `CachedMolHolder`, whose getters read back the
  pickles its `add_mol` writes. In Python, 256 brackets on a substance
  group make `Chem.Mol(m.ToBinary())` fail with "ENDMOL tag not found". The
  check covers molecule, reaction, bundle, tautomer-query, and extended-query
  pickles (and extended-query JSON), the substructure library's holders and
  archives, the molecule and filter catalogs, the fragment catalog's
  generator (whose fragments are pickled with one-byte counts whatever the
  molecule's size), `EnumerateLibrary` archives, synthon-space database
  files, and PNG metadata, from the metadata functions and from Cairo
  drawers with `include_metadata`.
- V3000 mol blocks, and the SD, mol-file, and rxn readers that read them,
  are `Error::Value` for an SGroup line (or its `DEFAULT` line) that ends
  after `COMPNO=` or a bare `COMPNO` label, and for an enhanced-stereo
  collection that lists no atoms for a nonzero count (`ATOMS=(1 )`). RDKit
  reads both values with extractions that write nothing there and then
  uses the unset variable; Valgrind reports a conditional jump on an
  uninitialized value in Python. Since RDKit parses a line cut off by a
  refusal or the end of the input as a whole one, a checked stream now
  holds back a COMPNO field or an enhanced-stereo line until it has read
  past it, reading ahead without seeking, and refuses one that spans more
  than 1 MiB.
- `json_to_mols` reads an attachment point without `lvIdx` as -1, the value
  RDKit's JSON writer leaves out. RDKit's reader left the point's
  `lvIdx` uninitialized, so a molecule written by `mol_to_json` read back
  an indeterminate leaving atom.
- `mol_from_smiles_with_params`, `mol_from_smarts_with_params`, and the PNG
  metadata readers are `Error::Value`, with `strict_cxsmiles` off, for a
  CXSMILES data substance group (`SgD:`) whose field has no `:` after it or
  whose coordinates have no `)`, and for a polymer group (`Sg:`) whose type
  has no `:` after it. RDKit's CXSMILES parser stops past the end of the
  extension there, and the non-strict parse then copies the text up to
  that point, up to 2 bytes past its end. Strict parsing, the default,
  rethrows RDKit's own parse error instead and is unchanged.
- Molecules with more than 65535 properties on the molecule, a conformer,
  an atom, a bond, or a substance group are refused wherever a molecule is
  pickled, as the one-byte counts above are (a panic from
  `Mol::to_binary` and `ChemicalReaction::to_binary`). MolPickler counts
  properties in a `u16`, which wraps, and its reader then takes the
  remaining property bytes for the pickle's next fields: in Python, 65536
  properties make every read fail ("Bond index out of range"), and 65535
  user properties beside the computed ones are silently dropped. Every
  property counts, whether or not the property flags write it; substance
  groups' properties are always pickled.
- `FragCatGenerator::add_frags_from_mol` counts every stereo group of the
  molecule against the one-byte limit, not only those a fragment can hold
  whole: from RDKit 2025.09 on, a fragment keeps each stereo group that
  shares an atom with it. With 256 three-atom groups and fragments of one
  bond, RDKit wrote a catalog it could not read ("ENDMOL tag not found").
- `expand_query` is `Error::Value` for a query atom or bond that has no
  query of its own: RDKit's `QueryBond::expandQuery` dereferences the null
  query, and `QueryAtom`'s fails an invariant. No safe Rust path makes such
  an atom or bond, so this only guards a path RDKit might add.
- `SynthonSpace::build_synthon_fingerprints` clears RDKit's Ctrl-C flag
  before and after the build, and the fingerprint searches and the
  database writers are `Error::Value` for a space whose fingerprints are
  incomplete. RDKit leaves its `SIGINT` handler installed after a
  substructure match, so a Ctrl-C afterwards sets its flag with nothing to
  clear it; the build then returned before fingerprinting anything, with
  the space's fingerprint type set, and the next fingerprint search
  dereferenced the missing fingerprints (a segfault, as in Python after a
  substructure search, a SIGINT, and `BuildSynthonFingerprints`).
- `DoubleCubicLatticeVolume::get_polar_surface_area_with_params` and
  `get_polar_volume_with_params` panic, with `include_hs`, for a molecule
  with a bond between two hydrogens, instead of overflowing the stack:
  RDKit 2026.09.1 recurses through the hydrogens forever (rdkit/rdkit#9687,
  fixed after that release; Python segfaults on H2). They stay infallible,
  as Python's raise for nothing else.

## 2026.9.0

The first release. It supports RDKit 2026.09.1 from conda-forge on Linux
(x86-64 and arm64) and macOS (arm64), and builds with Rust 1.88 or newer. The
build refuses other RDKit releases: the wrapper compiles against copies of
headers RDKit does not install, whose layouts can change between patch
releases (`wrapper/vendor/`). A program built against one release stops
when it loads with another.

### Bindings

- Binds RDKit's C++ API with the coverage of its Python wrappers, as
  conda-forge builds them with Boost.Python: 3,467 of the 4,039 names they
  expose. Each of the other 572 has a reason in `inventory/manifest.toml`:
  Python protocol methods and Boost.Python's container types, 52 names left
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
bindings and computes what Python computes: 1,016 of the 2,607 names RDKit's
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
