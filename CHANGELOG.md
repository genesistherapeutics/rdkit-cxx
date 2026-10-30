# Changelog

## Unreleased

### Changes

- The copyright holder is now Genesis Molecular AI, in `LICENSE` and
  in every source file's header. The license is still BSD-3-Clause.

## 2024.3.0

The first release for RDKit 2024.03: the bindings and ports of 2026.9.0, for
RDKit 2024.03.6, the series' last release (see Compared with 2026.9.0,
below). It supports RDKit 2024.03.6 from conda-forge on Linux (x86-64) and
macOS (arm64), and builds with Rust 1.88 or newer; conda-forge has no
linux-aarch64 build of RDKit 2024.03. The build refuses other RDKit
releases: the wrapper compiles against copies of headers RDKit does not
install, whose layouts can change between patch releases
(`wrapper/vendor/`). A program built against one release stops when it
loads with another. The environment needs conda-forge's `eigen` as well,
since RDKit 2024.03's `MolTransforms.h` includes Eigen.

### Compared with 2026.9.0

- What RDKit 2024.03 lacks is not bound, nor ported where Python had it:
  - Modules: `rdChemDraw`, `rdConformerParser`, `rdEnumerateStereoisomers`
    (the Python `EnumerateStereoisomers` port stays), `rdGaussianShape`,
    `rdMIF`, `rdMolProcessing`, `rdShapeAlign`, `rdSynthonSpaceSearch`, and
    the standardization `Pipeline`.
  - Molecules and atoms: `Atom.GetValence` and `ValenceType` (2024.03's
    `GetExplicitValence` and `GetImplicitValence` are bound instead;
    `get_implicit_valence` returns `i32`, as 2024.03 gives -1 before the
    valences are computed), `ClearPropertyCache`, `Mol.GetName` and
    `SetName`, `Atom.NOATOM` (`StereoInfo::NOATOM` stays), the monomer
    class and plain monomer info's residue fields (residue name, number, and
    chain are PDB residue info's), `CopyMolSubset`, the CIP ranks, stereo
    group canonicalization and cleanup, `CaptureErrorLog`, the periodic
    table's rows and `GetMaxAtomicNumber`, `RealValueVect`,
    `UniformRealValueGrid3D`, `AtomCoordsMatcher`, and newer substructure
    and RASCAL options.
  - File formats: SCSR, CDXML writing and ChemDraw support, CX reaction
    SMILES and SMARTS writers, `ReactionFromSmiles`, `clean_stereo`,
    `CX_ZERO_BONDS`, bond property lists, `MetadataFromPNGString`'s
    `asList`, and the PNG metadata parameters beyond 2024.03's
    `includePkl`, `includeSmiles`, and `includeMol`.
  - Descriptors: `NumAmideBonds`, `NumAtomStereoCenters`,
    `NumBridgeheadAtoms`, `NumHeterocycles`, `NumSpiroAtoms`,
    `NumUnspecifiedAtomStereoCenters`, and `Phi` (`desc_list` has 2024.03's
    210 descriptors); the double cubic lattice's per-atom, polar, and
    surface-point results; `CalcMolDescriptors3D`'s `confId`, which 2024.03
    ignores.
  - Embedding and force fields: the embedder's `timeout`, `basinThresh`,
    and legacy implementation, the newer failure causes, the `ETDGv2` and
    `DG` presets, the parameters' JSON, the bounds matrix of embedding
    parameters, and the MMFF property getters.
  - Drawing: `DrawElement`, the legend's position and vertical text, note
    colours, stereo-group labels, reagent padding, contour fill thresholds,
    the atom palette getter, the fingerprint-generator similarity maps, and
    `Draw`'s `MolToSVG` port, `DrawMolWithMatches`, and
    `MolsMatrixToGridImage`'s SVG output.
  - Reactions from CDXML (`ReactionsFromCDXMLBlock` and
    `ReactionsFromCDXMLFile`): RDKit 2024.03 reads them with its
    molecule-level CDXML parser, which reads atoms of fragments it has
    freed, as `MolsFromCDXML` and `MolsFromCDXMLFile` do.
  - Python ports: `Chem.MultiConfMolFromSDF`, `inchi.MolFromInchiAndAuxInfo`,
    and `DataStructs`' `getNForFlatMatrix` and `getElementFromFlatMatrix`.
  - Fingerprint generators' JSON and atoms-per-bit output, and the other
    options 2024.03's Python lacks (MolHash's tautomerism exclusion,
    R-group decomposition's multiple cores per molecule and target in
    results, scaffold networks' names, molzip's coordinate alignment,
    Kekulize's `canonical`, `DetermineBonds`' `maxIterations`, and the
    newer tautomer scores and validations).
- Where RDKit 2024.03 leaves a function one optional argument, that
  argument has a `_with_<param>` form rather than 2026.9.0's one-field
  params struct (CONVENTIONS.md A2): `kekulize_with_clear_aromatic_flags`,
  `kekulize_if_possible_with_clear_aromatic_flags`,
  `get_sssr_with_include_dative_bonds`,
  `get_symm_sssr_with_include_dative_bonds`,
  `reapply_mol_block_wedging_with_all_bond_types`, and
  `AtomMonomerInfo::new_with_name` replace `KekulizeParams`,
  `KekulizeIfPossibleParams`, `GetSSSRParams`, `GetSymmSSSRParams`,
  `ReapplyMolBlockWedgingParams`, and `AtomMonomerInfoParams` with their
  `_with_params` functions.
- Kept from 2024.03, which later releases correct: its Python's
  misspelled `FingeprintGenerator32` and `FingeprintGenerator64` classes
  (the Rust types keep their names), `Properties::register_property`'s
  return of the registry's size for a new property, NaN for normalizing a
  zero-length vector, and the Python ports' 2024.03 behavior
  (`FragmentDescriptors.csv`, whose `fr_prisulfonamd` pattern matches
  nothing; `TorsionFingerprints`' beta weights; `RegistrationHash`'s
  enhanced stereo; `FindMolChiralCenters`' legacy path, which keeps stale
  results; `GetUniqueCombinations`' recursive order, with
  `get_unique_combinations_new` for the sorted one; Butina clustering's set
  order, with `cluster_data_with_reordering` in place of
  `ClusterDataParams`; `EnumerateStereoisomers`, which keeps UNKNOWN bond
  directions; `MolToFile`'s messages).
- `RascalOptions`' default timeout is -1, not Python's 60: RDKit 2024.03's
  timeout path reads an emptied vector, so any other timeout is refused,
  and so is any `max_frag_separation` but -1: RDKit 2024.03 measures
  distances from a fragment it dropped, through a null pointer (#9424).
  `CalcGETAWAY` is not bound: RDKit 2024.03 writes through an empty
  vector's `begin()` on every call.
- conda-forge's macOS build of RDKit 2024.03.6 ships no SLN parser, so on
  macOS `rd_sln_parse`'s functions are `Error::Value`.
- Substructure matching still runs one call at a time under
  `control_c_mutex()`, though RDKit 2024.03 has no Ctrl-C handler; the
  mutex no longer saves and restores `SIGINT`.
- Checks for what RDKit 2024.03 gets wrong that later releases fixed, where
  safe Rust would reach it (`Error::Value` unless noted):
  - Commonchem JSON is checked for every value 2024.03's rapidjson parser
    reads unchecked, and an attachment point without `lvIdx` is read as -1.
  - Multithreaded suppliers, whose threads 2024.03 starts in the
    constructor and joins only after the derived supplier is destroyed, are
    read to their end before they are dropped, and their input is checked
    only once the threads are done; the property-list setting cannot change.
  - Removing atoms and bonds takes them out of stereo groups, which 2024.03
    leaves pointing at freed bonds; RDKit's own removals refuse groups that
    hold bonds or repeat atoms, and stereo groups refuse repeats.
  - Stereo groups that list an atom twice, which 2024.03 keeps and a
    removal leaves dangling (`[H][C@](C)(F)Cl |&1:0,0,1|` crashes Python),
    wherever groups come in: CXSMILES, V3000 enhanced-stereo collections
    (which must also list as many atoms as their count: 2024.03 reads a
    missing first atom uninitialized and repeats the last), Marvin atom
    arrays that give two atoms one id (`Error::Parse`), Commonchem JSON,
    molecule pickles, and `set_stereo_groups`.
  - Mol2 bond lines must name two atoms of the molecule (#9417): 2024.03
    refuses a bond only when both are out of range, and reads its atom
    table past the end for the other.
  - The PDB reader with proximity bonding, and `determine_connectivity` and
    `determine_bonds`, refuse a hydrogen RDKit would bond to several atoms,
    whose extra bonds RDKit 2024.03 removes while it iterates over them
    (#9444).
  - Molecules RDKit's pickler would write but could not read back are
    refused wherever a molecule is pickled: `Error::Value`, or a panic from
    `Mol::to_binary` and `ChemicalReaction::to_binary`. MolPickler writes
    each query node's number of children in one byte, and, in a molecule of
    at most 255 atoms, substance groups' list sizes, stereo groups' counts
    and sizes, and the bond indices both name; it truncates larger values.
    It writes property counts in two bytes, so more than 65,535 properties
    on a molecule, atom, bond, conformer, or substance group are refused
    too. The fragment catalog's generator checks its molecules' queries and
    property counts, but not their substance or stereo groups, which
    2024.03's fragments do not keep.
  - `RingInfo::num_relevant_cycles` is `Error::Value` for a count past
    `u32`, which RDKit converts unchecked.
  - `expand_query` on a query atom or bond without a query; RDKit 2024.03
    dereferences it.
  - Drawing: atom labels whose colon RDKit's colon adjustment indexes past
    the label's characters (`"A:"`), legends of only newlines, comic-mode
    arrowheads of no length, and a reaction layout checked as 2024.03 lays
    it out.
  - Molzip refuses linkers and handles the coordinates path's null and
    repeated fragments; `RenumberAtoms` needs a permutation;
    `FragmentOnBonds` refuses repeated bonds; `FragmentOnSomeBonds` of more
    than 63 bonds; `ReplaceCore` with an atom matched twice.
  - CXSMILES fields RDKit would read past the end of, wherever it parses
    SMILES or SMARTS with CXSMILES on (the readers, FMCS's initial seed,
    RASCAL's equivalent atoms, filter-catalog SMARTS, deprotections,
    scaffold-network reactions, the normalizations the standardizers read,
    query-definition files, recursive-query labels, `M  MRV SMA` lines and
    SMARTSQ data in mol blocks, and feature definitions); a polymer S-group
    only where strictCXSMILES is off. A normalizations file is read once and
    handed to RDKit as checked text or data, so a file of no transformations
    is refused where RDKit would read it from cleanup parameters. CML of tetrahedral atoms with fewer than four
    neighbors, InChI of atoms with more than 20 bonds or stereo it would
    overrun, and V3000 reaction agents; V3000 S-group values and S-group
    lines RDKit 2024.03 reads uninitialized when the line ends inside them
    or lacks their ids.
  - CoordGen of double bonds marked E or Z without stereo atoms, the
    Hückel method where YAeHMOP cannot open the output files RDKit 2024.03
    writes in the current directory (`Error::Runtime`), FMCS's `Verbose`
    with `StoreAll`, more than 32 CIP stereo
    configurations, MolHash's Extended Murcko scaffold of more than 4000
    atoms, and, where SMILES are parsed and molecules sanitized or kekulized,
    Kekulize's enumeration of 32 or more dummy atoms (README, Known
    limitations, for the other paths).
  - RDKit 2024.03's C-string exceptions (a bond query SMARTS cannot write)
    are `Error::Runtime`, as Python's RuntimeError, without reading the
    dangling string; `MaeWriter`
    writes through a stream of its own, since 2024.03's shared-stream
    constructor always throws; the Marvin reaction reader goes through
    2024.03's stream reader, whose block reader is declared but not
    exported.

### Bindings

- Binds RDKit's C++ API with the coverage of its Python wrappers, as
  conda-forge builds them with Boost.Python: 2,885 of the 3,322 names they
  expose. Each of the other 437 has a reason in `inventory/manifest.toml`:
  Python protocol methods and Boost.Python's container types, names left
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
  contract is bytes RDKit wrote. A few gaps remain in RDKit itself; see the
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
- `PointND.Distance`, which Python binds to `Point3D`'s function so that it
  always raises, is not bound.

### Ports of Python-only functionality

Much of what RDKit offers only in Python is ported to Rust on top of the
bindings and computes what Python computes: 1,000 of the 2,643 names RDKit's
Python modules define. The ported modules' Python tests are ported with them,
with expected values from RDKit 2024.03.6.

- **Descriptors:** `Descriptors` (`desc_list` in Python's order,
  `calc_mol_descriptors`), `Descriptors3D`, `Crippen`, `MolSurf`,
  `GraphDescriptors`, `EState`, `Fragments`, `Lipinski`, `QED`,
  `ReducedGraphs`, `SpacialScore`, `ChemUtils.DescriptorUtilities`, and
  `ML.Descriptors.MoleculeDescriptors`. `calc_mol_descriptors` matches
  Python's 210 descriptors on a set of 40 molecules to 1e-12.
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
- Some Python functions are broken in RDKit 2024.03.6; the ports implement
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
