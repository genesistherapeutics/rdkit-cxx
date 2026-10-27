# Changelog

## 2025.3.0

The first release for RDKit 2025.03: the bindings and ports of 2026.9.0, for
RDKit 2025.03.6 (see Compared with 2026.9.0, below). It
supports RDKit 2025.03.6 from conda-forge on Linux (x86-64 and arm64) and
macOS (arm64), and builds with Rust 1.88 or newer. The build refuses other
RDKit releases: the wrapper compiles against copies of headers RDKit does not
install, whose layouts can change between patch releases (`wrapper/vendor/`).
A program built against one release stops when it loads with another.

### Compared with 2026.9.0

- conda-forge's macOS build of RDKit 2025.03.6 ships no SLN parser (and no
  `rdSLNParse` for Python), so on macOS `rd_sln_parse`'s functions are
  `Error::Value`; its Linux builds have the parser. RDKit 2025.03's CMake
  package also needs Eigen, which `environment.yml` installs.
- What RDKit 2025.03 lacks is not bound (and nor is anything its Python
  lacks):
  - `rdGaussianShape`, `rdBase.CaptureErrorLog`, and
    `rdMolHash.excludeFromTautomerismProp`.
  - Embedding: the initial-embedding modes, embedding force field,
    `onlyInitialEmbedding`, `ETKDGv4`, Gaussian torsion terms, the `DG`
    preset, `useLegacyImplementation`, the parameters' JSON form, the
    ordered failure causes, and the failure causes after `END_OF_ENUM` (12):
    minimization, k-term violations, clashes, and the timeout.
    `GetMoleculeBoundsMatrix` takes its four flags rather than embedding
    parameters.
  - Rings: `GetSSSR`'s and `GetSymmSSSR`'s options other than
    `includeDativeBonds`, `FindRingFamilies`' parameters, the legacy
    ring-finding setting, and `ReapplyMolBlockWedging`'s `verify`.
  - Molecules and atoms: `Mol.GetName`/`SetName`, `Atom.NOATOM` (the
    stereo-info constant remains), the monomer class (so
    `AtomMonomerInfo::new_with_name` takes the one optional argument left),
    the residue fields of plain `AtomMonomerInfo` (they are
    `AtomPDBResidueInfo`'s), subset copying
    (`CopyMolSubset`), the attachment-point label helpers, the coordinate
    atom matcher, and the substructure matcher's extra atom and bond checks.
  - Molecule files: bond property lists, and the CDXML writer of
    `rdmolfiles` (2025.03 writes ChemDraw from `rdChemDraw`:
    `rd_chem_draw::mol_to_chem_draw_block`).
  - Drawing: legend position and vertical legends, the drawing-extent
    elements, the stereo-group label options and annotation, single-colour
    bonds, and the atom palette getter; `Draw.MolToSVG` and
    `DrawMolWithMatches`, which 2025.03's `Draw` takes from `rdMolDraw2D`
    or lacks, and the similarity maps' fingerprint-generator functions.
  - `MMFFMolProperties`' term, dielectric and variant getters; MolAlign's
    `BestAlignmentParams` (the best-alignment functions take their
    arguments one by one, and keep hydrogens, as 2025.03 does); the
    fingerprint generators' atoms-per-bit output and JSON form; the
    abbreviations' `includesXBonds`; R-group decomposition's
    `allowMultipleCoresInSameMol`; molzip's `alignCoordinates`;
    `Compute2DCoords`' `usePathAngleExpansion`; `GetSurfacePoints` and
    the radii form of `DoubleCubicLatticeVolume`; `GetAllConformerBestRMSToRef`;
    `MolCatalogEntry.Serialize` (an entry still pickles through
    `serialize`); and in `rdShapeAlign`, custom features, normalization, the
    inertial rotation, the shape-options overloads, and the score functions.
  - The synthon search's shapes, incremental and line-range searches,
    possible-hits files, best hit, and synthon counts; its database files
    are format 3.0 or 3.1, and `SynthonSpaceSearchParams` gains
    `num_random_sweeps`.
  - Python helpers 2025.03's Python lacks: `Chem.MultiConfMolFromSDF` and
    `DataStructs.getNForFlatMatrix` and `getElementFromFlatMatrix`.
- Results are RDKit 2025.03.6's: its 2D depictions (and the drawings and mol
  blocks made from them), conformers, bounds matrices, kekulé forms (a ring's
  first bond is double), ring atom orders, pickles, and MolInterchange JSON
  (RapidJSON's number format, and H counts kept on query atoms) differ from
  2026.09's. Reaction CX SMILES and SMARTS carry atom maps as `atomProp`;
  bit vectors of no bits are fully similar; `path_to_submol` refuses a bond
  index past the end, as 2025.03 throws; a query bond SMARTS cannot write is
  `Error::Runtime`; a pickled monomer type of `PDBRESIDUE` reads as PDB
  residue info. Parsing does not find ring families.
- RDKit 2025.03 behaviors the bindings keep, as Python has them: the
  `Properties` registry computes `chi2v` and `chi2n` as `chi3v` and `chi3n`
  (RDKit #9643).
- The ports follow RDKit 2025.03.6's Python: `BertzCT` takes a non-aromatic
  bond's order as its `BondType`'s number (17 for a dative bond),
  `fr_Ar_COO`'s description is 2025.03's, `fr_prisulfonamd`'s SMARTS
  (`[OX1;+0;-1]`) never matches, `TorsionFingerprints` tests the end atom's
  neighbors twice where it means each end's, `CalcMolDescriptors3D` has no
  conformer argument (2025.03's ignores it), and `RandomizeActivities`
  shuffles (2025.03's other path raises `NameError`).
- Checks for what RDKit 2026.03 or 2026.09 fixed in RDKit:
  - Every wrapper source reads `RDConfig.h` first, since `RingInfo.h`
    declares its ring-family members only under a macro it does not
    include.
  - `RingInfo.num_relevant_cycles` is `Error::Runtime` before the ring
    families are found, where RDKit 2025.03 converts RingDecomposerLib's
    failure value out of range, and `Error::Value` for more relevant cycles
    than `u32` holds.
  - `MolCatalogEntry::from_pickle` builds the entry as the catalog's reader
    does: RDKit 2025.03's pickle constructor deletes an uninitialized
    pointer.
  - `assign_cip_labels` leaves unlabeled, as RDKit 2026.09 does, an
    atropisomer bond with an end of no other bonds, where RDKit 2025.03
    reads past an empty vector.
  - A molecule of at most 255 atoms whose substance groups or stereo groups
    name a bond past 255, which RDKit 2025.03 pickles truncated, is refused
    wherever it would be pickled (`Error::Value`; `Mol::to_binary` and
    `ChemicalReaction::to_binary` panic). So are a substance group or stereo
    group list of more than 255 entries, more than 255 stereo groups, and a
    query of more than 255 children, which RDKit 2026.09 pickles truncated
    too: RDKit writes these pickles and cannot read them back.
  - A multithreaded supplier's queues hold at most 2^20 records: RDKit
    2025.03 dereferences null when it cannot allocate one.
  - Drawing refuses an atom of two bonds with a neighbor at its position,
    for which RDKit 2025.03 throws from the middle of laying out the
    drawing, and a legend of only newlines, for which its `extractLegend`
    reads past an empty vector.
  - `renumber_atoms` takes only a permutation: RDKit 2025.03 reads past the
    end for an index equal to the atom count and adds an atom twice for a
    repeated one.
  - `fragment_on_bonds` and `fragment_on_some_bonds` refuse a repeated
    bond, which RDKit 2025.03 uses after removing it, and
    `fragment_on_some_bonds` returns no fragments when asked to break more
    bonds than it is given, where RDKit shifts past 64 bits.
  - `replace_core` with a match refuses two core atoms on one molecule
    atom whose core neighbor the match lacks (RDKit 2025.03 dereferences an
    end iterator), and `molzip_fragments` and `molzip_row` with
    `generate_coordinates` refuse a zip that fails or drops an input atom,
    which 2025.03 then dereferences.
  - `calc_getaway` with a custom atom property refuses a molecule with no
    two atoms 8 bonds apart, whose terms RDKit 2025.03 reads uninitialized.
  - RASCAL takes only -1 for `timeout` and `max_frag_separation`: RDKit
    2025.03 reads an emptied list when a search times out, and dereferences
    a dropped fragment when it applies a separation. The default timeout
    (and the clustering pairs') is -1, where Python's is 60 seconds.
  - `DoubleCubicLatticeVolume` refuses a nonzero dot density, whose dot
    counts RDKit 2025.03 leaves uninitialized, monomer info other than PDB
    residue info, which it casts to that, and probe radii and coordinates
    beyond 1e6, and, with `is_protein` and without the ligand, a molecule
    with no atom of a standard residue, whose grid RDKit 2025.03 leaves
    uninitialized.
  - Parsers refuse what RDKit 2025.03 misreads: Mol2 bonds naming atom 0
    or an atom past the end; V3000 SGroup lines it would read unset values
    from or pop an empty SAP id from, and `BEGIN SGROUP` with no SGroups in
    COUNTS (checked streams hand RDKit no part of a refused line); CXSMILES
    `atomProp:` lists with no `|` or `,` after them; stereo groups that list
    an atom or bond twice (from Rust, pickles, CXSMILES, V3000, Marvin, and
    JSON), which RDKit uses after freeing; any MolInterchange JSON value
    RDKit 2025.03's RapidJSON reader takes without a type check (a Parse
    error), and a ranged `AtomRingSize` query, which it casts to an
    equality query (#9172); and in CDXML reactions, a fragment's bond ordering longer than
    its connection points, binary CDX passed as text, and XML entity
    declarations.
  - `json_to_mols` gives an attachment point without `lvIdx` -1, which
    RDKit 2025.03 leaves uninitialized.
  - FMCS refuses `verbose` with `store_all`, for which RDKit 2025.03 matches
    against a null query molecule (#9040).
  - `mol_hash`'s `ExtendedMurcko` refuses more than 4,000 atoms, for which
    RDKit 2025.03 allocates its arrays on the stack.
  - `mol_from_smiles`, `sanitize_mol`, `kekulize`, and
    `kekulize_if_possible` refuse 32 or more dummy atoms in one aromatic
    ring system, and `assign_cip_labels` more than 64 stereo atoms and
    bonds, which RDKit 2025.03 counts with shifts past the width of their
    masks (README, Known limitations).
  - `mol_to_cml_block` and `mol_to_cml_file` refuse a tetrahedral atom with
    fewer than four neighbors, whose missing neighbor RDKit 2025.03 reads
    (add hydrogens first), and drawing refuses atom labels whose last drawn
    character is a colon or with a `<` and no `>` after it.
  - The `MaeWriter` takes its stream as Python's does: RDKit 2025.03's
    `shared_ptr` constructor always fails its own precondition.
  - The embedder, `DetermineBondOrders`, the CIP labeler (and the
    fingerprints that run it), and the synthon search leave RDKit's
    `SIGINT` handler installed; the wrappers restore the process's after
    each call.

### Bindings

- Binds RDKit's C++ API with the coverage of its Python wrappers, as
  conda-forge builds them with Boost.Python: 3,204 of the 3,702 names they
  expose. Each of the other 498 has a reason in `inventory/manifest.toml`:
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
  contract is bytes RDKit wrote. A few gaps remain in RDKit itself (its log flags, two
  oversized shifts, and SMARTS it parses from parameter files); see the
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
bindings and computes what Python computes: 1,008 of the 2,603 names RDKit's
Python modules define. The ported modules' Python tests are ported with them,
with expected values from RDKit 2025.03.6.

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
- Some Python functions are broken in RDKit 2025.03.6; the ports implement
  what they intend:
  - Under NumPy 2 or Python 3.12: `ML.InfoTheory` `BitRank`,
    `ML.Data.Stats.StandardizeMatrix` on 2-D input,
    `ML.Cluster.Resemblance.EuclideanDistance` and `CalcMetricMatrix`,
    `ML.files.ReadDataFile`, and `RDRandom.shuffle` without `random`.
  - `ML.Data.Stats.PrincipalComponents` returns rows of LAPACK's matrix
    rather than eigenvectors; the port returns eigenvectors.
  - `SupplierFromFilename` fails on `.tdt` files.
  - `ML.Cluster.ClusterUtils.FindClusterCentroidFromDists` and
    `ML.Cluster.Resemblance.ShowMetricMat` index lists with floats, so they
    raise `TypeError` for two or more points.
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
