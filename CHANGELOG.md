# Changelog

## 2024.9.0

The first release for RDKit 2024.09: the bindings and ports of 2026.9.0, for
RDKit 2024.09.6 (see Compared with 2026.9.0, below). It supports RDKit
2024.09.6 from conda-forge on Linux (x86-64 and arm64) and macOS (arm64),
and builds with Rust 1.88 or newer. The build refuses other RDKit releases:
the wrapper compiles against copies of headers RDKit does not install, whose
layouts can change between patch releases (`wrapper/vendor/`). A program
built against one release stops when it loads with another. RDKit 2024.09's
CMake package needs Eigen 3, so the environment installs `eigen 3.4.*`.

### Compared with 2026.9.0

- conda-forge's macOS build of RDKit 2024.09.6 ships no SLN parser (and no
  `rdSLNParse` for Python), so on macOS `rd_sln_parse`'s functions are
  `Error::Value`; its Linux builds have the parser.
- What RDKit 2024.09 lacks is not bound:
  - Modules: `rdGaussianShape`, `rdEnumerateStereoisomers` (the Python
    `EnumerateStereoisomers` port stays), and `rdChemDraw`, with the
    ChemDraw writer.
  - Functions: the SCSR readers, the bond property lists, `CopyMolSubset`,
    `CleanupStereoGroups`, `ComputeAtomCIPRanks`, `CaptureErrorLog`,
    `GetInchiVersion`, `inchi.MolFromInchiAndAuxInfo`, the
    `ClearPropertyCache` methods, `Mol.GetName` and `SetName`,
    `AtomCoordsMatcher`, `Draw.MolToSVG` and `DrawMolWithMatches`, the
    fingerprint generators' JSON and similarity maps,
    `OrderedEmbedFailureCauses`, `EmbedParametersToJSON`, the `ETDGv2` and
    `DG` presets, `GetAllConformerBestRMSToRef`, the attachment-point label
    helpers, and `MolCatalogEntry.Serialize` (an entry still pickles
    through `serialize`).
  - Values: `Atom.NOATOM`, `CX_ZERO_BONDS`, and the embedding failure
    causes after `EXCEEDED_TIMEOUT`.
  - Options: the drawing's legend position, drawing extents, stereo group
    labels, and single-colour bonds; contours' fill thresholds and the
    similarity maps' grid options; the embedder's initial-embedding modes,
    embedding force field, `onlyInitialEmbedding`, `ETKDGv4`, Gaussian
    torsion terms, and legacy implementation; MMFF's term, variant, and
    dielectric getters; RascalMCES' `completeSmallestRings` and
    `maxBestMCESs`; R-group decomposition's `allowMultipleCoresInSameMol`;
    `MolzipParams.alignCoordinates`; the depictor's `useRingTemplates` and
    `usePathAngleExpansion`; `DetermineBonds`' `maxIterations`;
    `Kekulize`'s and `KekulizeIfPossible`'s `canonical` (`KekulizeParams`
    has `clear_aromatic_flags` alone, and kekulizing backtracks up to 100
    times, as Python's does);
    `GetSSSR`'s and `GetSymmSSSR`'s `includeHydrogenBonds`, `algorithm`,
    and `recalcSSSR`, and the legacy ring-finding setting;
    `CanonicalizeStereoGroups`' `maxStereoGroups`;
    `ReapplyMolBlockWedging`'s `verify`;
    `AbbreviationDefinition.includesXBonds`; `SubstructMatchParameters`'
    extra checks and `aromaticMatchesSingleOrDouble`; MolHash's
    `excludeFromTautomerismProp`; and `GetSurfacePoints`.
- API in 2024.09's forms: `Atom.GetExplicitValence` and
  `GetImplicitValence` (no `GetValence`); `AddHs`, `GetBestRMS`,
  `GetBestAlignmentTransform`, `GetAllConformerBestRMS`, and
  `MolMetadataToPNG*` take keyword arguments, so their parameter structs
  (`AddHsParams` and others) replace `AddHsParameters`,
  `BestAlignmentParams`, and `PNGMetadataParams`; `AtomMonomerInfo` has a
  type and a name, and the residue fields are `AtomPDBResidueInfo`'s;
  `DoubleCubicLatticeVolume` computes its results when it is made;
  `rd_shape_align` has `prepare_conformer`, `align_mol`, and
  `align_mol_shape_input`; `GetMoleculeBoundsMatrix` takes `set15bounds`,
  `scaleVDW`, `doTriangleSmoothing`, and `useMacrocycle14config`.
- `rdSynthonSpaceSearch`, new in RDKit 2024.09, is not bound: its readers
  and searches are not audited (U1). Later releases of the crate bind it.
- `CalcGETAWAY` is not bound: RDKit 2024.09 writes through an empty
  vector's `begin()` on every call (fixed in 2025.03.4, RDKit #8601).
- `ReactionsFromCDXMLBlock` and `ReactionsFromCDXMLFile` are not bound:
  they run RDKit 2024.09's own CDXML parser, which reads atoms of fragments
  it has freed. Nor are `MolsFromCDXML` and `MolsFromCDXMLFile`.
- Results are RDKit 2024.09.6's: depictions, conformers, drawings, pickles,
  descriptor values, and the ports' outputs, where they differ from
  2026.09's. RDKit 2024.09 behaviors the bindings keep, as Python has them:
  the `Properties` registry computes `chi2v` and `chi2n` as `chi3v` and
  `chi3n` (RDKit #9643); the unpickler reads monomer type 1 as PDB residue
  info and refuses unknown types; the SD suppliers read the text after the
  last `$$$$` as a record; and `fr_prisulfonamd` matches nothing.
- The ports follow RDKit 2024.09.6's Python: `CalcMolDescriptors3D`
  ignores `confId`; legacy `FindMolChiralCenters` keeps
  `_ChiralityPossible`; `EnumerateStereoisomers` clears only
  `EITHERDOUBLE` bond directions and sets no chiral flag;
  `RegistrationHash` canonicalizes enhanced stereo whatever the flags;
  TorsionFingerprints' weights test only a bond's end atom for a terminal
  bond; `RandomizeActivities` without shuffling, which raises `NameError`,
  is `Error::Value`; and `BertzCT` takes a non-aromatic bond's order as its
  `BondType`'s number (17 for a dative bond).
- Checks for what later RDKit releases fixed:
  - Stereo groups (RDKit #8558): RDKit 2024.09's bond removals leave a
    group that holds the bond pointing at it, and its groups may list an
    atom twice, which a removal takes out once. `RWMol`'s removals take the
    bond out of the groups, as RDKit 2025.09 does; RDKit's own removals
    (the SMILES writer for several fragments, `get_mol_frags_as_mols`,
    `remove_hs`, the transforms, the standardizer, and the rest) refuse a
    molecule whose groups they would leave dangling (`Error::Value`), or
    that has an atom or bond in two ABS groups, which 2024.09's
    `insertMol` (molzip, `replace_substructs`) merges into one group that
    repeats it, and
    the pickle readers, `set_stereo_groups`, and the CXSMILES, V3000,
    Marvin, and JSON readers refuse groups that list an atom or bond twice
    (the Marvin reader, two atoms with one id), since RDKit's readers
    remove hydrogens before any later check could run.
  - The V3000 reader refuses an enhanced-stereo collection with fewer atoms
    than its count, for which RDKit 2024.09 reads an atom number no
    extraction wrote (still so in 2026.09).
  - `fragment_on_bonds` refuses a repeated bond (#9430);
    `fragment_on_some_bonds` makes no fragments for more bonds than it is
    given (#9451).
  - The CML writer refuses a tetrahedral atom of fewer than four neighbors
    (#9472), and the Mol2 reader a bond line that names an atom the
    molecule lacks (#9417).
  - The input scanner reads each V3000 substance-group line as RDKit does
    and refuses a field that would leave a value uninitialized or pop an
    empty string (#9417).
  - The PDB reader with proximity bonding, and `determine_connectivity`,
    refuse a hydrogen RDKit would bond to several atoms, whose extra bonds
    RDKit 2024.09 removes while it iterates over them (#9444).
  - FMCS refuses `verbose` (#9040), RascalMCES any `max_frag_separation`
    but -1 (#9424), and CoordGen double-bond stereo without two stereo
    atoms (#8415).
  - `molzip_fragments` with `generate_coordinates` refuses repeated R
    groups and `FragmentOnBonds` with several molecules, and returns `None`
    where RDKit would read a null molecule; `replace_core` refuses a match
    missing a shared core atom's neighbor (#9450).
  - Drawing refuses a legend of only newlines.
  - RDKit 2024.09 reads JSON with rapidjson, its assertions off:
    `json_to_mols` refuses every value of another type, missing key, and
    short array it would read unchecked, and gives a missing attachment
    point `lvIdx` -1, which 2024.09 leaves uninitialized.
  - Pickling refuses more than 65,535 properties on a molecule, atom,
    bond, conformer, or substance group, which RDKit writes in two bytes.
  - `renumber_atoms` refuses anything but a permutation (#9450), and
    `DoubleCubicLatticeVolume` what 2024.09's reads unchecked (a dot
    density, monomer info of another type, short PDB atom names, huge
    coordinates, and a protein without its ligand or a standard residue).
  - A multithreaded supplier stops its input and drains its queues before
    RDKit's destructor joins its threads (#8556); dropping one mid-read
    otherwise calls a pure virtual method or waits forever.
  - The Maestro writer is made through a stream it may own, since RDKit
    2024.09's shared-pointer constructor always fails.
  - The CXSMILES readers (SMILES, SMARTS, reaction SMARTS, and the patterns
    the standardizer, filters, MMPA, and depiction templates parse) refuse
    an `atomProp` block that runs to the end of the text (#9361), which
    RDKit 2024.09 reads past, and, without `strict_cxsmiles`, a
    substance-group field that does.
  - Drawing refuses a label whose colon RDKit 2024.09's `adjustColons`
    would read past or loop on.
  - The `ExtendedMurcko` hash refuses more than 4000 atoms: RDKit 2024.09
    takes a byte of stack per atom with `alloca`, which past a page can
    step over the stack's guard page.
  - RascalMCES refuses any timeout but -1, which is its default here
    (Python's is 60 seconds): when a search runs past it, RDKit 2024.09
    reads its best cliques after clearing them (#9450). The clusterings
    run their searches without one. The query expanders refuse a query
    atom or bond without a query.
- Checks the bindings for RDKit 2026.03 made, as they stand for 2024.09:
  - `MolCatalogEntry::from_pickle` builds the entry as the catalog's reader
    does, since 2024.09's pickle constructor deletes an uninitialized
    pointer too.
  - `assign_cip_labels` leaves out an atropisomer bond with a bare end,
    clearing any label it had, as RDKit 2026.09 does, since 2024.09 reads
    past its empty carriers for one, and refuses more than 64 stereo atoms
    and bonds at a time, past which 2024.09's descriptor pairing shifts by
    a negative count (`Error::Value`).
  - `RingInfo.num_relevant_cycles` is `Error::Runtime` before the ring
    families are found, and every wrapper source reads `RDConfig.h` first,
    as for 2026.03.
  - A multithreaded supplier's queues hold at most 2^20 records.
  - Drawing an atom of two bonds with a neighbor at its position is
    RDKit's `Error::Runtime`: 2024.09 normalizes the direction, and throws,
    where 2026.03 read an uninitialized bond type.
  - `draw_arrow` checks the head RDKit 2024.09 draws: without a polygon
    its sides are lines, drawn by hand in comic mode, which refuses sides
    of no length; 2024.09 has no precondition on the angle, so a flat
    arrow is drawn, as Python draws it. The reaction's arrow is drawn as
    lines, with no miter, and comic mode checks its shaft as such.
  - The pickle-width checks (`check_picklable`) hold for 2024.09's
    MolPickler, which writes the same widths.
  - Only the embedder leaves RDKit's `SIGINT` handler installed in 2024.09;
    the wrappers restore the process's after it.

### Bindings

- Binds RDKit's C++ API with the coverage of its Python wrappers, as
  conda-forge builds them with Boost.Python: 3,079 of the 3,591 names they
  expose. Each of the other 511 has a reason in `inventory/manifest.toml`:
  Python protocol methods and Boost.Python's container types, 67 names left
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
  contract is bytes RDKit wrote. Some gaps remain in RDKit itself, none of
  which corrupts memory; see the README's "Known limitations".

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
bindings and computes what Python computes: 1,008 of the 2,604 names RDKit's
Python modules define. The ported modules' Python tests are ported with them,
with expected values from RDKit 2024.09.6.

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
- Some Python functions are broken in RDKit 2024.09.6; the ports implement
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
