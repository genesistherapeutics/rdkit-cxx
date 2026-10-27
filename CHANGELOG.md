# Changelog

## 2025.9.0

The first release for RDKit 2025.09: the bindings and ports of 2026.9.0, for
RDKit 2025.09.6 (see Compared with 2026.9.0, below). It
supports RDKit 2025.09.6 from conda-forge on Linux (x86-64 and arm64) and
macOS (arm64), and builds with Rust 1.88 or newer. The build refuses other
RDKit releases: the wrapper compiles against copies of headers RDKit does not
install, whose layouts can change between patch releases (`wrapper/vendor/`).
A program built against one release stops when it loads with another.

### Compared with 2026.9.0

- conda-forge's macOS build of RDKit 2025.09.6 ships no SLN parser (and no
  `rdSLNParse` for Python), so on macOS `rd_sln_parse`'s functions are
  `Error::Value`; its Linux builds have the parser.
- What RDKit 2025.09 lacks is not bound (API it lacks is removed, not
  stubbed):
  - `rdGaussianShape`, which RDKit 2025.09 does not have.
  - The embedder's initial-embedding modes, embedding force field,
    `onlyInitialEmbedding`, `ETKDGv4`, Gaussian torsion terms, and
    `useLegacyImplementation`; the DG preset; the `MINIMIZATION`,
    `KTERM_VIOLATION`, and `CLASH` failure causes; and
    `GetMoleculeBoundsMatrix`'s `EmbedParameters` form and its
    `forceTransAmides`, `set14bounds`, and `set13bounds`.
  - `GetSymmSSSR`'s `algorithm` and `recalcSSSR`, the legacy ring-finding
    setting, the attachment-point label helpers, `Compute2DCoords`'
    `usePathAngleExpansion`, `GetSurfacePoints`' `allPoints`,
    `GetAllConformerBestRMSToRef`, and `MolCatalogEntry.Serialize` (an
    entry still pickles through `serialize`).
  - `Kekulize`'s and `KekulizeIfPossible`'s `canonical`; `Mol.GetName` and
    `SetName`; `CaptureErrorLog`; `rdMolHash.excludeFromTautomerismProp`;
    the MMFF term, dielectric, and variant getters of `MMFFMolProperties`;
    and the residue fields and `monomerClass` of `AtomMonomerInfo` (they are
    `AtomPDBResidueInfo`'s alone).
  - Drawing's `legendPosition`, `legendVerticalText`, stereo-group label
    options, and `singleColourBonds`; stereo groups are labeled `abs`,
    `or<n>`, and `and<n>`.
  - `MolToCDXMLBlock`: `rd_chem_draw` writes CDXML and CDX with RDKit
    2025.09's `MolToChemDrawBlock`.
  - The synthon search's shape search, possible-hits files, line-range
    searches, best hit, and progress bar; `ConvertTextToDBFile`'s shape
    parameters. Its database files are format 3.1.
  - SMARTS' `k` ring-size primitive, and the Python-only `Draw.MolToSVG`,
    `Chem.MultiConfMolFromSDF`, and `inchi.MolFromInchiAndAuxInfo`.
- Results are RDKit 2025.09.6's: its Kekulé forms (2025.09's kekulization is
  not canonical), 2D depictions (and the drawings and mol blocks made from
  them), conformers, bounds matrices, ring atom orders, embedding failure
  counts, the synthon search's order of hits, pickles (version 16.2), and
  JSON's `toolkitVersion`. Parsing does not find ring families.
- RDKit 2025.09 behaviors the bindings keep, as Python has them: the
  `Properties` registry computes `chi2v` and `chi2n` as `chi3v` and `chi3n`
  (RDKit #9643); `DoubleCubicLatticeVolume.get_surface_points` returns the
  points the surface-area methods tested, once per call; the synthon
  search's `max_hit_mol_wt` defaults to 0 (no limit); and a monomer-info
  pickle of type `PDBRESIDUE` reads as an `AtomPDBResidueInfo`.
- The ports follow RDKit 2025.09.6's Python: `BertzCT` takes a non-aromatic
  bond's order as its `BondType`'s number (17 for a dative bond),
  `fr_Ar_COO`'s description is 2025.09's, and `CalculateTorsionLists` checks
  only a bond's end atom for heavy neighbors.
- RASCAL's `timeout` and `max_frag_separation` take only -1, as on the
  rdkit-2025.03 line, and other values are `Error::Value`: RDKit 2025.09
  reads an emptied list of cliques when a search times out, and measures
  fragment separation from a fragment it has just pruned. The timeout's
  default is -1, not Python's 60 seconds, and the clusterings run their
  searches without one.
- Checks for what newer RDKit releases fixed (each refused input is
  `Error::Value`):
  - `renumber_atoms` takes only a permutation: RDKit 2025.09 writes past its
    array for an index equal to the atom count and reads unplaced atoms for
    a repeated one.
  - `fragment_on_bonds` and `fragment_on_some_bonds` refuse repeated bond
    indices, which RDKit 2025.09 uses after freeing; `fragment_on_some_bonds`
    returns no fragments for more bonds to break than it has, where RDKit
    shifts past a word.
  - Stereo groups may not hold an atom or bond twice: RDKit 2025.09's
    constructor accepts duplicates, and removing such an atom leaves a
    dangling pointer. They are refused where they are made
    (`create_stereo_group`, `from_data`), before the wrappers' removals, and
    in the inputs whose readers remove hydrogens themselves before any check
    could run: CXSMILES enhanced stereo, V3000 collections, Marvin files
    (whose stereo groups take an atom by id, so ids must differ), and JSON.
  - `replace_core` with matches checks that each core atom's neighbors are
    matched; `molzip` with depiction refuses fragment lists RDKit 2025.09
    reads past.
  - The GETAWAY descriptors with custom properties refuse molecules of
    topological diameter below 8, for which RDKit 2025.09 reads
    uninitialized arrays.
  - The CML writer refuses a tetrahedral atom of fewer than four explicit
    neighbors, whose parity RDKit 2025.09 reads past them for.
  - The Mol2 readers refuse a bond whose atom number is outside the
    molecule's atoms.
  - Drawing refuses an atom label whose first colon is the last character
    drawn (`C:` segfaults in RDKit 2025.09's Python), and a legend of only
    newlines, which RDKit 2025.09 lays out with no lines.
  - Each V3000 SGroup line is replayed as RDKit 2025.09 consumes it, as on
    the rdkit-2025.03 line, and refused where RDKit would read a value no
    extraction wrote: a header without its type or external index, a
    `CSTATE` count or bond, a `COMPNO`, an array's first value, or an empty
    `SAP` identifier, whose last character RDKit pops. So is an
    enhanced-stereo collection shorter than its count. A CTAB line is
    scanned to its end before any of it is delivered from a stream, so a
    refused line is cut before its first byte.
  - CXSMILES extensions are checked as RDKit 2025.09 reads them, as on the
    rdkit-2024.09 line: an `atomProp` block that runs to the end of the
    text, which RDKit reads past, and, without `strict_cxsmiles`, a data or
    polymer substance group's field without its `:`, which RDKit copies
    past the end. They are checked in SMILES, SMARTS, reactions, PNG
    metadata, deprotections, normalizations (from data, files, and cleanup
    parameters), cleanup parameters' other SMARTS, query definitions,
    recursive-query properties, tautomer score terms, `SmartsMatcher`
    patterns (a refused one leaves the matcher without a pattern, as one
    that does not parse), SMILES mol holders, scaffold-network bond
    breakers, and the wrappers' other SMILES and SMARTS inputs. RDKit opens
    normalization and query-definition files itself; the wrappers read them
    once, check them, and hand RDKit the text.
  - `path_to_submol` and `copy_mol_subset` clear a double bond's CIS or
    TRANS (or an atropisomer's) whose stereo atoms they left out, which
    RDKit 2025.09 keeps and `MolToSmiles` then reads past.
  - Proximity bonding (PDB input and `determine_connectivity`) refuses a
    hydrogen that could be bonded to three or more atoms, whose bonds RDKit
    2025.09 prunes while it iterates over them.
  - JSON refuses a ranged `AtomRingSize` query, which RDKit 2025.09 reads as
    an `EqualityQuery`. A JSON attachment point without `lvIdx`, which RDKit
    2025.09 leaves uninitialized, gets -1, as RDKit's writer means by leaving
    it out.
  - The `ExtendedMurcko` hash takes fewer than 4,096 atoms: RDKit 2025.09
    takes a byte per atom from the stack.
  - Ring-system templates are checked as SMILES, as RDKit 2025.09 parses
    them.
  - Every wrapper source reads `RDConfig.h` first, since `RingInfo.h`
    declares its ring-family members only under a macro it does not
    include.
  - `RingInfo.num_relevant_cycles` is `Error::Runtime` before the ring
    families are found, where RDKit 2025.09 converts RingDecomposerLib's
    failure value out of range, and `Error::Value` for more relevant cycles
    than `u32` holds.
  - `MolCatalogEntry::from_pickle` builds the entry as the catalog's reader
    does: RDKit 2025.09's pickle constructor deletes an uninitialized
    pointer.
  - `assign_cip_labels` leaves unlabeled, as RDKit 2026.09 does, a double
    bond whose stereo atom is not a neighbor of its own end and an
    atropisomer bond with an end of no other bonds, where RDKit 2025.09
    reads past an empty vector.
  - A molecule of at most 255 atoms whose substance groups or stereo groups
    name a bond past 255, which RDKit 2025.09 pickles truncated, is refused
    wherever it would be pickled (`Error::Value`; `Mol::to_binary` and
    `ChemicalReaction::to_binary` panic). So are a substance group or stereo
    group list of more than 255 entries, more than 255 stereo groups, and a
    query of more than 255 children, which RDKit 2026.09 pickles truncated
    too, and a molecule, conformer, atom, bond, or substance group of more
    than 65535 properties, whose count RDKit writes in two bytes: RDKit
    writes these pickles and cannot read them back. The fragment catalog's
    generator checks the groups a fragment of the molecule can hold as a
    small molecule's pickle writes them, whatever the molecule's size,
    counting every stereo group, since a fragment keeps each one that shares
    an atom with it.
  - A multithreaded supplier's queues hold at most 2^20 records: RDKit
    2025.09 dereferences null when it cannot allocate one.
  - Drawing refuses an atom of two bonds with a neighbor at its position,
    for which RDKit 2025.09 throws partway through a layout.
  - RDKit 2025.09's embedder, CIP labeler, `DetermineBondOrders`, and
    synthon search leave RDKit's `SIGINT` handler installed; the wrappers
    restore the process's after each call.

### Bindings

- Binds RDKit's C++ API with the coverage of its Python wrappers, as
  conda-forge builds them with Boost.Python: 3,303 of the 3,851 names they
  expose. Each of the other 548 has a reason in `inventory/manifest.toml`:
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
bindings and computes what Python computes: 1,013 of the 2,598 names RDKit's
Python modules define. The ported modules' Python tests are ported with them,
with expected values from RDKit 2025.09.6.

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
- Some Python functions are broken in RDKit 2025.09.6; the ports implement
  what they intend:
  - Under NumPy 2 or Python 3.12: `ML.InfoTheory` `BitRank`,
    `ML.Data.Stats.StandardizeMatrix` on 2-D input,
    `ML.Cluster.Resemblance.EuclideanDistance` and `CalcMetricMatrix`,
    `ML.files.ReadDataFile`, and `RDRandom.shuffle` without `random`.
  - `ML.Data.Stats.PrincipalComponents` returns rows of LAPACK's matrix
    rather than eigenvectors; the port returns eigenvectors.
  - `SupplierFromFilename` fails on `.tdt` files.
  - `ML.Cluster.ClusterUtils.FindClusterCentroidFromDists` and
    `ML.Cluster.Resemblance.ShowMetricMat` index with `/`, a float, and
    raise `TypeError` (RDKit 2026.03 uses `//`).
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
- `ML.Data.DataUtils.RandomizeActivities` without `shuffle` raises
  `NameError` in RDKit 2025.09.6 (RDKit 2026.03 gave it another meaning), so
  `randomize_activities_with_shuffle` is `Error::Value` without it, and the
  shuffle-free `randomize_activities` is not ported.
