# Conventions

Every rule here binds all later wrapper code. A change that cannot follow a
rule changes the rule first.

"Python" in these rules is RDKit's Python interface as conda-forge ships it:
its Boost.Python modules and its `.py` modules.

## Files

- **L1.** Every source file starts with the BSD header in
  `.license-header.txt`, with "All rights reserved." on its own line. The
  pre-commit hook adds it (`addlicense -c "Genesis Molecular AI, Inc." -f
  .license-header.txt`), another hook rejects the one-line form, and CI runs
  the hooks. Generators write the same header themselves
  (`inventory/check.py`). Vendored third-party files
  (`wrapper/third_party/`, and RDKit's headers in `wrapper/vendor/<release>/`)
  stay as their upstream ships them, under their own licenses, and the hooks
  skip them.
- **L2.** rustfmt formats Rust, clang-format (`.clang-format`: Google's style
  with RDKit's pointer alignment) formats C++, and ruff (`ruff.toml`) lints
  and formats Python. The pre-commit hooks in `.pre-commit-config.yaml` run
  them all.

## Bridges

- **B1.** Every bridge function is a free function in namespace `rdkit_cxx`,
  declared in `wrapper/include/rdkit-cxx/<area>.h` and defined in
  `wrapper/src/<area>.cpp`. Bridges never bind RDKit members directly.
- **B2.** Every bridge function returns `Result`. A C++ exception that reaches
  cxx through a function without `Result` aborts the process.
- **B3.** One bridge, `src/ffi/<area>.rs`, per wrapper header. An opaque type
  is declared in one bridge and reused in others as
  `type X = crate::ffi::<area>::ffi::X;`.
- **B4.** A new area appends to `BRIDGES` in `build.rs`, to `RDKIT_CXX_SOURCES`
  in `wrapper/CMakeLists.txt`, and to `RDKIT_CXX_RDKIT_TARGETS` when it needs
  another RDKit library.
- **B5.** Wrapper functions for a class's methods are named
  `<class>_<method>` in snake_case (`mol_get_num_atoms`), so names stay unique
  across bridges.

## Errors

- **E1.** A public function returns `Result` when the Python call can raise
  for some input reachable through the safe Rust API: an index, a key,
  unparsable text, a failed sanitization, a missing property cache or ring
  info, a missing conformer, invalid UTF-8. Otherwise it returns the value,
  and an unexpected C++ exception panics through `expect_ok`.
- **E2.** Parsers that return `None` in Python return `Err(Error::Parse)`. A
  parse that fails sanitization returns the sanitization variant.
- **E3.** Index arguments are checked in the C++ wrapper, which throws
  `rdkit_cxx::IndexError`. (Python raises `RuntimeError` from RDKit's own
  range check there.)
- **E4.** C++ exceptions map to `Error` variants as follows, most derived
  first:

  | C++ exception | Variant |
  | --- | --- |
  | `RDKit::AtomValenceException` | `AtomValence` |
  | `RDKit::AtomKekulizeException` | `AtomKekulize` |
  | `RDKit::KekulizeException` | `Kekulize` |
  | `RDKit::MolSanitizeException` | `MolSanitize` |
  | `rdkit_cxx::IndexError` | `Index` |
  | `KeyErrorException` | `Key` |
  | `ValueErrorException`, `RDKit::ConformerException`, `RDKit::ChemicalReactionException`, `RDDepict::DepictException` | `Value` |
  | `RDKit::BadFileException` | `Os` |
  | `RDKit::SmilesParse::SmilesParseException`, `RDKit::FileParseException`, `RDKit::SLNParseException`, `RDKit::ChemicalReactionParserException` | `Parse` |
  | `Invar::Invariant` (with `toUserString()`), any other exception | `Runtime` |

- **E5.** Every bridge's `include!` list starts with `"rdkit-cxx/errors.h"`,
  which defines the `rust::behavior::trycatch` that encodes the exception.

## Ownership and threads

- **W1.** Every Rust molecule owns a C++ `RDKit::RWMol`, so in-place
  operations can take `RWMol&` without a cast.
- **W2.** Shared views are `XRef<'m>` (`Copy`, holding `&'m Mol` and
  `&'m ffi::X`). Exclusive views are `XMut<'m>` (holding `&'m mut Mol` and
  an index). An owned `X` exists only when Python has a constructor (for
  `Bond`, the `rdqueries` functions that build query bonds).
- **W3.** Getters are written once, in a `macro_rules!` that stamps them on
  every type that has them (owned, `Ref`, `Mut`), all calling wrappers that
  take `&ffi::X`. Setters work the same way, on owned and `Mut` types.
- **W4.** Methods that need the owning molecule (degree, hydrogen counts,
  ring membership, neighbors, bonds, `get_owning_mol`) exist only on views.
- **W5.** Python sequences become iterators named for the result (`Atoms`,
  `Bonds`) that implement `ExactSizeIterator` and `DoubleEndedIterator`.
  There are no mutable iterators. Code mutates by index, because two live
  `XMut`s into one molecule would alias. Filtered sequences (`GetAromaticAtoms`) return a
  `Vec` of views.
- **W6.** Functions that return new molecules copy into an `RWMol`; vectors
  of molecules cross as `MolList`, moved out one at a time by
  `mol_list_into_vec`.
- **W7.** `RWMol` wraps a `Mol` and implements `Deref<Target = Mol>` and
  `DerefMut`, so every `Mol` method works on it; `From<Mol> for RWMol` and
  `From<RWMol> for Mol` move without copying. Methods Python defines on
  `RWMol` live on `RWMol`. In-place operations (`sanitize_mol`, `kekulize`,
  property-cache updates) take `&mut Mol`, which an `RWMol` also provides;
  operations that return a new molecule take `&Mol`.
- **W8.** Classes that hold molecules (`MolBundle`) own them as `Vec<Mol>`,
  so they can lend `&Mol`. Each call into C++ builds a container whose
  `shared_ptr`s borrow those molecules (a no-op deleter) and drops it before
  the call returns. Functions that take Python lists or dicts of molecules
  lend them the same way, through `BorrowedMols`: `&[Mol]` or `&BTreeMap`
  when RDKit only reads them, `&mut` when it writes them (`molzipFragments`
  gives its inputs coordinates).
- **W9.** Objects RDKit shares through `shared_ptr`s, between owners Rust
  cannot see (`FilterCatalog`'s entries and matchers, shared by catalogs,
  entries, and match results), are shared by Rust handles too, and changed
  only when no other owner exists: a handle copies a shared object before
  it changes it (`FilterCatalogEntry`, by its `use_count`), or owns it alone
  because RDKit copies whatever another object holds (the concrete
  matchers, whose `Clone` is RDKit's `copy`). A borrowed child in such a
  tree (`FilterHierarchyMatcherMut`) borrows its root mutably, as W2's views
  do.
- **S1.** `Mol` and `RWMol` are `Send` and not `Sync`. RDKit writes computed
  properties through `const` molecules (`MolToSmiles` sets
  `_smilesAtomOutputOrder` in `SmilesWrite.cpp`), so two threads reading one
  molecule race. Parallel code gives each thread its own molecule. Owned
  `Atom`s and `Bond`s, outside any molecule, own all their data and are
  `Send` too.

- **S2.** The crate needs an RDKit built with `RDK_BUILD_THREADSAFE_SSS`, as
  conda-forge's is: process-wide settings such as the default pickle
  properties are safe to set from safe Rust only because RDKit guards them
  with a mutex. The wrappers refuse to compile otherwise.
- **S3.** Module-level functions that take any atom are generic over the
  sealed traits `AsAtom` (`Atom`, `AtomRef`, `AtomMut`) and `AsAtomMut`
  (`Atom`, `AtomMut`). Methods take `&dyn AsAtom` instead, so they are not
  generic and `tests/manifest_paths.rs` can name them. Methods that take any
  bond take `&dyn AsBond` (`Bond`, `BondRef`, `BondMut`).
- **S4.** A method that may change state RDKit keeps for later calls takes
  `&mut self`, even when Python's runs it on a const object.
  `ChemicalReaction`'s runs initialize the reaction first when it is not
  (validation writes its product templates' properties), so they, `initialize`,
  and `validate` take `&mut self`; a reaction is not initialized behind a
  shared reference. Classes that hold molecules RDKit may write through const
  (`ChemicalReaction`, `EnumerateLibrary`) are `Send` and not `Sync`, as S1
  says of molecules, and their bridges record why in a `SAFETY` comment.
- **S5.** Classes whose `const` methods RDKit itself runs from several
  threads (`FilterCatalog`, its entries, and its matchers, which
  `RunFilterCatalog` shares) are `Send` and `Sync`. Their query molecules
  are matched concurrently, which writes recursive queries' matches under
  each query's lock (`RDK_BUILD_THREADSAFE_SSS`), so a wrapper that copies
  such a molecule, which reads those matches, holds the same locks, taken in
  `SubstructMatch`'s order after `control_c_mutex()`, as matches take them
  (`copy_pattern` in `filter_catalog.cpp`). Matching holds
  `control_c_mutex()` (U4) and releases it while a Rust matcher's methods
  run, which RDKit calls between its matches.
- **S6.** An exception that leaves a `std::thread` calls `std::terminate`,
  so RDKit's own threads run wrapper work only where RDKit catches what it
  throws: `MultithreadedMolSupplier`'s parsing threads catch everything
  their write callback throws (`MultithreadedMolSupplier::writer`), so
  `rdMolProcessing` fingerprints each molecule there, as RDKit's
  `getFingerprintsForMolsInFile` does. Work that needs one of U4's locks
  does not run there: the calling thread, holding the lock, does it.

## Unchecked input

- **U1.** No input reachable from safe Rust may reach undefined behavior in
  RDKit, even where Python's bindings do. When RDKit trusts an argument (an
  index, a size, a value it casts or sums in `int`, a pickle, a file), the
  C++ wrapper checks it first and throws as E3 and E4 describe, usually
  `ValueErrorException`. Each such check carries a comment naming what RDKit
  would do without it.
- **U2.** Binary formats RDKit parses without validating them (bit-vector,
  value-vector, sparse-int-vector, and grid pickles) are checked by
  `rdkit-cxx/pickle_checks.h` before RDKit reads them, or, where parsing
  itself is safe, on the object right after. FPB files are read by the
  wrapper, checked, and handed to RDKit as those checked bytes, so RDKit
  never reads a file that changed after the check; PNG files are read and
  their chunks checked the same way. Parsers RDKit fuzzes
  (`Code/Fuzz`: SMILES and mol blocks) are trusted, except for the limits
  U12 checks. RDKit fuzzes molecule pickles too, but only their parsing:
  MolPickler takes a query's class from its tag and its matching function
  from its description independently, as the JSON reader below does, and
  converts values it reads without range checks, so a pickle RDKit did not
  write reaches undefined behavior while it is read or in code that reads
  the molecule later. Every function that hands caller bytes to MolPickler
  is therefore an `unsafe fn` whose contract is a pickle RDKit wrote,
  unmodified: molecule, bundle, reaction, tautomer-query, and
  extended-query pickles, extended-query JSON (which holds molecule
  pickles in base64), catalog pickles, PNG metadata, and
  `CachedMolHolder::add_binary`, whose bytes the holder's safe methods
  read later. MolPickler writes some counts, and in RDKit 2026.03 a small
  molecule's bond indices, in one byte, and reads a larger value's pickle
  out of step, and counts each object's properties in two bytes, so every
  wrapper that hands RDKit a molecule to pickle, alone or inside another
  format (`git grep -l 'pickleMol\|MolPickler::pickle\|pickleReaction'`
  in RDKit's sources), first calls `check_picklable`
  (`rdkit-cxx/mol_checks.h`), or one of the variants beside it, and
  refuses one it would truncate: RDKit could not read the pickle back, as
  the contract and the holders' safe reads assume. The wrappers still check what they checked before, for
  pickles that meet the contract but come from another RDKit: a molecule
  pickle's enum values, which MolPickler casts unchecked (bond types,
  directions, and stereo, atoms' chiral tags and hybridizations, stereo
  groups' types) and which the Rust enums cannot hold past RDKit's last
  value, `check_pickled_mol` refuses where U10 checks monomer types; and
  catalog pickles (`MolCatalog`, `FragCatalog`, and their entries), which
  frame molecule pickles with lengths and indices RDKit trusts,
  `rdkit-cxx/catalog_checks.h` reads as RDKit will, with RDKit's readers
  for the molecules and parameters, checking those and the molecules'
  monomer infos (U10), before RDKit reads the same bytes.
  Boost.Serialization archives of polymorphic objects (`EnumerateLibrary`'s
  `Serialize` and `GetState`, `FilterCatalog`'s pickles) are not checked:
  Boost trusts their class and
  object ids, so the functions that read them are `unsafe fn`s whose
  contract is bytes RDKit wrote, and the wrappers check only what such
  bytes can still get wrong (a state from another library; matchers that
  RDKit dereferences null in or recurses through forever, and patterns'
  monomer infos). Synthon-space database files
  (`SynthonSpace::read_db_file`), records of molecule and bit-vector
  pickles and Boost archives of shapes at offsets RDKit trusts, are read by
  an `unsafe fn` too, whose wrapper checks the header and offsets first.
  RDKit's JSON reader (`MolInterchange/Parser.cpp`), which
  no fuzzer covers, trusts the indices, counts, and query trees it reads,
  so `json_to_mols` parses the document first, with the same boost::json
  parser, checks what RDKit would read unchecked (`DocumentChecker` in
  `mol_interchange.cpp`, which follows RDKit's reads and leaves RDKit's own
  errors to RDKit), and hands RDKit the checked document serialized again.
  Its queries must be ones RDKit makes, since RDKit takes a query's class
  from its tag and its meaning from its description, while the code that
  reads queries takes the description to imply the class: each
  description with a class RDKit gives it, AND, OR, and XOR with two or
  more children, and no bond XOR or null bond in an AND or OR, which
  RDKit's own bond matcher and SMARTS writer cast. Text a wrapper checks
  before RDKit parses it
  is read once: a file is read through a checked stream (U12) or into
  memory, so RDKit parses the bytes that were checked; where RDKit reads a
  file only by name (a synthon search's possible hits), it reads a private
  copy of the checked lines.
- **U3.** A sized object whose `int` arithmetic would overflow (more than
  `INT_MAX` bits or grid points) is refused at construction: constructors
  that return values panic, documented under `# Panics`, and fallible ones
  return `Error::Value`. The exception is the unfolded fingerprints RDKit
  makes itself (`GetSparseFingerprint`'s `SparseBitVect` of 2^32 − 1 bits):
  the operations whose `int` arithmetic they would overflow (pickling,
  base64, conversion to `ExplicitBitVect`) refuse them instead.
- **U4.** Global state RDKit writes without a lock is guarded by a mutex in
  the wrapper (the grid-index cache in `GridUtils.cpp`, the BRICS patterns
  `fragmentOnBRICSBonds` fills on first use, the descriptor property
  registry in `Property.cpp`, `DCLV.cpp`'s `recordCache`, and R-group
  decomposition's random number singleton and fingerprint-variance scoring
  static, under `rgroup_mutex()`, and the shape modules' pharmacophore
  patterns and PubChem shape's overlap cutoffs, under `pubchem_shape_mutex()`
  and `gaussian_shape_mutex()`, and the static work variables of Murtagh's
  hierarchical clustering, under `murtagh_mutex()`, and the ChemDraw
  library's attribute-name maps, which it fills behind an unlocked check,
  under `chemdraw_mutex()`). State RDKit makes once,
  on first use, is
  made under `std::call_once` before any other call (the MACCS patterns in
  `MACCS.cpp`). 2D depiction's settings and
  CoordGen's templates are globals every coordinate generation reads, so
  every wrapper whose RDKit call may compute 2D coordinates holds
  `lock_depiction()` (`depiction_lock.h`): the depiction functions; the
  mol-block (V2000 and V3000), Marvin, ChemDraw, PNG, SD, and Maestro
  writers for a molecule without a conformer, and the V2000 rxn-block
  writer built on them (also in reaction PNG metadata) for a reaction with
  such a template; and `molzip` with
  `generateCoordinates`. They also check, under the lock, what the depictor
  that will run trusts about the molecule (`check_compute_2d_coords_mol` in
  `depiction_checks.h`). RDKit's random number generator
  (`RDGeneral/utils.cpp`) is shared by every module, so it has one lock,
  `random_mutex()` in `rdkit-cxx/global_state.h`, held for the whole call by
  every wrapper whose RDKit call can seed it or draw from it
  (`grep -rn 'getRandomGenerator\|getRandomVal\|getDoubleRandomSource' Code`
  in RDKit's sources): the SMILES and CXSMILES writers (for molecules and
  reactions) with `doRandom`, `MolToRandomSmilesVect`, PNG metadata with
  random SMILES, the depictor's sampling (`nFlipsPerSample` and random bond
  flips, in `Compute2DCoordsForReaction` too; the writers' default depiction
  never draws), the embedder without a seed, DistGeom's `EmbedBoundsMatrix`,
  and MolAlign's `RandomTransform`. RDKit's Ctrl-C handler (`ControlCHandler`, which swaps
  the SIGINT handler through static variables) has its own lock,
  `control_c_mutex()` in `rdkit-cxx/control_c_lock.h`, held for the whole call
  by every wrapper whose RDKit call constructs one. Every substructure match
  does (`Substruct/vf2.hpp`), so that is the match methods and every wrapper
  whose RDKit call matches patterns: the transforms that match, reactions,
  alignment, the force-field builders, Crippen and the other descriptors and
  fingerprints that type atoms by SMARTS, depiction with templates,
  standardization, the catalogs, MCS, SubstructLibrary, R-group
  decomposition, scaffold networks, chemical feature factories, the
  shape modules, and InChI conversion, besides the embedder,
  DetermineBonds, and the CIP labeler. Such a call runs on the calling
  thread the work RDKit would split across its own threads, which the lock
  cannot order (the wrapper's own threads take it in turn), and releases the
  lock while a Rust closure it calls between matches runs (a property
  functor, tautomer callback or score, validation method, filter matcher, or
  FMCS typer or callback; new wrappers use `HoldControlC` and
  `ReleaseControlC`). A closure RDKit calls inside a match (FMCS's final
  match checks) runs with the lock held, so a wrapper it calls that takes
  this lock, or one ordered before it, fails instead
  (`check_not_in_control_c_callback`), as does one a `ReleaseControlC`
  runs while the thread holds the lock through anything but `HoldControlC`,
  which it cannot release. MolAlign's
  `symmetrizeTerminalAtoms`, which the best-alignment functions and the
  embedder's pruning reach, assigns a static on every call, so callers hold
  `symmetrize_mutex()` (`rdkit-cxx/symmetrize_lock.h`). A wrapper that needs
  several of these locks takes them in this order, after any of the locks
  above for one module's state: `lock_depiction()`, `control_c_mutex()`,
  `symmetrize_mutex()`, then `random_mutex()`. `tests/threads.rs`, which CI
  runs under Valgrind's DRD (the `races` job), failing on any race in RDKit's
  global state or in RDKit or wrapper code (`tests/check_drd_log.py`), calls,
  from every module that writes the state one of these locks guards,
  wrappers that write it (drawing from the generator and constructing a
  `ControlCHandler` both write), and at least one that only reads it.
- **U5.** RDKit settings kept in environment variables
  (`SetUseLegacyRingFinding`, `SetAllowNontetrahedralChirality`,
  `SetUseLegacyStereoPerception`) are set with `setenv`, which races with
  every `getenv` in the process, as `std::env::set_var` does. Their setters
  are `unsafe fn` with `set_var`'s safety contract; the getters are safe.
- **U6.** RDKit caches distance, adjacency, and path matrices as unsized
  arrays in a molecule's computed properties (`Matrices.cpp`), and reads
  them back whenever they exist, while `RWMol::addAtom`, `addBond`,
  `removeBond`, and `insertMol` keep them, so a stale matrix is read past its
  end. No molecule Rust holds carries one between calls: a wrapper whose
  RDKit call computes one (`grep getDistanceMat\|get3DDistanceMat\|getAdjacencyMatrix`
  in the RDKit sources) holds a `MatrixCacheGuard` (`rdkit-cxx/mol_checks.h`)
  on the molecules it touches, which clears them before the call and after
  it, even when it throws; molecules copied from RDKit-owned ones, and
  RWMol's edits, clear them too. Python's `force` and `prefix` arguments,
  which only steer the cache, are not bound.
- **U7.** RDKit's queries compare an atom's or bond's value with theirs by
  subtracting in `int` (`Queries::queryCmp`). The `rdqueries` value queries
  read values below 2^30 in magnitude (atom fields are 8 or 16 bits; masses,
  scaled by 1000, stay below 6.6e7), so they refuse query values beyond
  2^30 (after the mass scaling) with `Error::Value`. Integer property
  queries read values the target sets, so the wrapper builds them from a
  subclass of RDKit's that compares in 64 bits; it pickles as RDKit's own.
- **U8.** RDKit's descriptors cache their results, and per-atom
  contributions, in a molecule's computed properties (`Crippen.cpp`,
  `MolSurf.cpp`, `ConnectivityDescriptors.cpp`, `PMI.cpp`, and the principal
  axes `PMI.cpp` and `PBF.cpp` read from `MolTransforms.cpp`), and read them
  back unless forced. The caches survive the edits U6 lists and ignore the
  conformer, and the Chi descriptors index their cached per-atom vectors by
  atom without checking the length, so a stale one is read past its end. As
  with matrices, no molecule Rust holds carries one between calls: a wrapper
  whose RDKit call reads or writes one holds a `DescriptorCacheGuard`
  (`rdkit-cxx/mol_checks.h`), which removes them by name before the call and
  after it. Python's `force` arguments, which only steer these caches, are
  not bound, and every result is computed afresh.
- **U9.** RDKit stores an atom's formal charge in 8 bits: `setFormalCharge`
  wraps any `int` into −128 to 127, as Python's does, so `set_formal_charge`
  takes any `i32`, and RDKit's per-atom arithmetic on charges (valences,
  effective atomic numbers, atom invariants) stays in range. Sums over atoms
  do not: `MolOps::getFormalCharge`, `getMolFormula` (which negates the
  sum), MolHash's formula, mesomer, redox-pair, net-charge, and tautomer
  hashes (the tautomer and protomer hashes subtract the charge from a
  hydrogen count also summed in `int`, at most 255 per atom),
  `ResonanceMolSupplier`'s conjugated groups, and the CML writer add charges
  in `int`, which only 2^24 or more atoms can overflow. `get_formal_charge`
  adds the charges in 64 bits, in RDKit's order, and refuses the molecule
  when the running sum leaves `int`'s range; the others call
  `check_formal_charge_sum` (`rdkit-cxx/mol_checks.h`), which refuses a
  molecule whose charges' magnitudes, plus that hydrogen bound where RDKit
  adds one, exceed `INT_MAX`. Both refuse with `Error::Value`. (EEM's sum
  is bounded by its 46339-atom limit.) A new binding whose RDKit call sums
  charges (`grep 'getFormalCharge()'` in the RDKit sources) checks the same
  way.
- **U10.** RDKit takes a monomer info's type for its class: code that finds
  `PDBRESIDUE` (`grep PDBRESIDUE` in the RDKit sources: `MolPickler.cpp`,
  `AddHs.cpp`, the PDB and sequence writers, `ProximityBonds.cpp`) casts the
  info to `AtomPDBResidueInfo`, which reads past the end of a plain
  `AtomMonomerInfo`. So no atom holds a plain info of that type, though
  Python allows one: `AtomMonomerInfo::new` refuses the type with
  `Error::Value`, as does `set_monomer_type` unless `dynamic_cast` finds PDB
  residue info, and molecule, bundle, and PNG-metadata pickles (and the
  reaction pickles, `EnumerateLibrary` archives, and catalog pickles and
  archives that hold them) whose atoms have such an info, or a type
  `AtomMonomerType` lacks, are refused right after RDKit reads them
  (`check_pickled_mol` in `rdkit-cxx/mol_checks.h`, which also checks the
  molecules of recursive queries, as MolPickler pickles them whole, and
  walks the queries with its own stack). RDKit code may then trust the
  type; wrapper code that casts uses `dynamic_cast` all the same. RDKit code
  that casts whatever the type (`DCLV.cpp`) still needs its own check.
- **U11.** RDKit keeps the SMILES and SMARTS parsers' debug levels in globals
  (`yysmiles_debug`, `yysmarts_debug`) that every parse reads. A parse writes
  them only when its `debugParse` differs, so the wrappers refuse a nonzero
  `debugParse` with `Error::Value`, and no parse writes them. RDKit's own
  `AtomFromSmiles`, `BondFromSmiles`, `AtomFromSmarts`, and `BondFromSmarts`
  (and `RDKitSVGToMol`, which calls them) write them on every call, racing
  with parses the wrapper cannot see, so their Rust functions are `unsafe fn`
  whose contract is that no other thread parses SMILES or SMARTS meanwhile.
- **U12.** Text RDKit reads with its line-oriented readers (mol blocks, SD,
  SMILES, and TDT files, and rxn blocks, whose templates are mol blocks) is
  checked as RDKit reads it, by the streams of
  `rdkit-cxx/streams.h` (or `check_input_text` for text in memory): RDKit
  counts lines in `int`, so input is refused past 2^31 - 2^20 lines, and it
  matches V3000 enhanced-stereo lines with `std::regex`, which recurses once
  per character (about 300 bytes of stack each), so such lines are refused
  past 4000 characters. RWMol's atom removals (`removeHs`, which the readers
  run) rewrite each V3000 bond's `ENDPTS` value without checking it, so a
  value that is empty, or in parentheses holds anything but numbers below
  2^32 between spaces, is refused too. A checked stream stops delivering bytes before the
  first one past a limit, and the wrapper throws `ValueErrorException` once
  RDKit asks for that byte, from that call and every later one: RDKit's
  readers trust the positions they record, and a cut-off read can record
  `tellg`'s −1. RDKit's random-access suppliers reread their files, so a
  checked stream also keeps a keyed hash of each 64 KiB chunk it scanned and
  fails when a chunk rereads differently or the file ends sooner than it did,
  so RDKit parses only the bytes that were checked. Python file objects are
  bound as `Read` values, read through the same streams, and `Write`
  values, which a supplier or writer may borrow for its lifetime. Gzipped files are
  decompressed as RDKit's `gzstream` decompresses them and checked after;
  the stream RDKit reads them through throws for every seek and tell, as
  `gzstream`'s does. RDKit's writers ignore
  a stream that fails, so the wrappers own the output (a file they open, or
  the Rust writer), hand RDKit the stream without ownership, and check it
  after every call: a failed write is `Error::Os` from then on.
- **U13.** `RWMol::replaceBond` deletes the bond it replaces but, unlike
  `replaceAtom`, leaves the stereo groups that held it (an atropisomer's)
  pointing at it, for copying, pickling, and the writers to read. So the
  wrappers replace bonds with `replace_bond` (`rdkit-cxx/mol_checks.h`),
  which points those groups at the new bond, and an RDKit call that replaces
  bonds of a molecule Rust will see runs under a `StereoGroupBondGuard`,
  which takes the bonds out of the groups and puts back those then at the
  same indices (or on a copy whose groups hold no bonds, if nothing returns
  it); `mol_checks.h` lists RDKit's calls and why the others need neither.

## Properties

- **P1.** Properties follow A2 and A4 like any method. Missing keys are
  `Error::Key`, type mismatches and invalid UTF-8 are `Error::Value`, and the
  getters and setters come from `prop_getters!` and `prop_setters!` rather
  than being written per class. Every class with properties derives from
  `RDKit::RDProps`, so one set of wrappers takes `const RDProps&`, and each
  class has an upcast, `<class>_as_props`. The `default=` overloads of the
  getters are not bound: `has_prop` covers them. A class whose Python
  `GetPropsAsDict` has other defaults (`SubstanceGroup`) passes them to
  `prop_getters!`; `GetPropsAsDictParams::default()` stays `Mol`'s.
- **P2.** Properties RDKit reads back as data, which any caller can set,
  are input under U1. When RDKit removes atoms or bonds (`RWMol`'s removals
  and batch edits, and the RDKit calls built on them) it parses every
  bond's `_MolFileBondEndPts` without checks and follows the substance
  groups' `PARENT` chains (by their `index`) with no end for a cycle, so a
  wrapper whose RDKit call removes atoms or bonds from a molecule, or from
  a copy of one, calls `check_removals` (`rdkit-cxx/mol_checks.h`) first,
  before any edit in place. Such calls include, besides the removals and
  the functions named for them, everything that splits a molecule into its
  fragments (`getMolFrags` with molecules, and so the SMILES writer for a
  molecule of several fragments, and what writes SMILES: hashes, tautomer
  enumeration, SMILES holders), the standardizer, and the embedder. Where
  RDKit combines molecules before it removes from them (`molzip`,
  `replaceSubstructs`, a synthon space's products), their chains are
  checked together, since RDKit does not renumber the groups' indices.
  The enumerator's link nodes (`_molLinkNodes`), position variation bonds,
  and SRU groups are checked in `mol_enumerator.cpp` as RDKit reads them.

## Documentation

- **D1.** Doc examples are `no_run` or `compile_fail`; behavior is tested
  under `tests/`. Each public item's doc comment names the Python object it
  mirrors, or says it is a Rust addition (an iterator, a view, a trait, or a
  type that names what Python leaves untyped). Hooks that only the crate's
  tests need (a C++ value to pin, an exception to raise) are `#[cfg(test)]`
  and crate-private, with their tests in the crate's unit tests, so they are
  never public API.

## Arguments and overloads

- **A1.** Required Python parameters become positional Rust parameters, in
  Python's order.
- **A2.** A Python signature with optional parameters is bound as `name`
  (using Python's defaults) plus, for exactly one optional parameter,
  `name_with_<param>`, or for two or more, `name_with_params(&<Name>Params)`.
  `<Name>Params` is named for the Python function or method, has the optional
  parameters as snake_case fields, and implements `Default` with Python's
  defaults, pinned by a test. A Python name may itself contain "With"
  (`SaltRemover.StripMolWithDeleted` is `strip_mol_with_deleted`); A2's
  suffix then follows it (`strip_mol_with_deleted_with_dont_remove_everything`,
  `lazy_pick_with_threshold_with_params`), and only the last `_with_` is
  A2's.
- **A3.** When Python also accepts a params object (an overload taking, say,
  `SmilesParserParams`), `name_with_params` takes the Rust mirror of that
  class instead. Keyword arguments the class lacks become extra fields on it.
  The mirror's `Default` reads a default-constructed C++ object, so it cannot
  drift.
- **A4.** An optional parameter that changes the return type gets its own
  function, `name_<param>` (`GetMolFrags(asMols=True)` →
  `get_mol_frags_as_mols`), and does not count toward A2. Python output
  arguments (lists the call fills in) are not bound as arguments; the return
  value carries that information, from the function itself or, when that
  would change its return type, from a `_with_outputs` form (A9).
- **A5.** Overloads that differ by an argument's type: the first that Python
  tries keeps the name, and the others append `_<type>` in snake_case
  (Boost.Python tries overloads last registered first, the reverse of the
  order its docstrings and the inventory list them). Constructors: a
  zero-argument `__init__` is `Default`; a copy `__init__` is `Clone`; the
  first other `__init__` is `new` (with `new_with_params` for its optional
  parameters, per A2); the rest are `from_<first parameter>`, or a `From` impl
  when the parameter is another bound class (`RWMol(mol)` → `From<Mol>`),
  which the manifest records in `overloads_note`. Overloads that differ in
  which parameters they take, not only in a type, get a suffix naming what
  sets them apart, recorded in the manifest's `overloads_note`:
  `ReplaceCore(mol, coreQuery)`, which matches the core itself, is
  `replace_core_query`; `CopyMolSubset(mol, atomIndices, bondIndices)` is
  `copy_mol_subset_atoms_bonds`; `molzip(a)` is `molzip_single`, and
  `molzip(row)` is `molzip_row`.
- **A6.** `_mut` comes last: `get_conformer_with_id_mut`.
- **A7.** Parameter structs are plain Rust structs with public fields. A
  mirror of a Python class crosses the bridge as a cxx shared struct
  `<Name>Data` with the same fields (maps as parallel key and value vectors),
  converted by a private `to_ffi`; its `Default` comes from a wrapper that
  returns the C++ defaults. Structs that A2 defines, with no Python class,
  write Python's defaults in Rust. Parameter structs stay exhaustive, so
  callers can write `S { x, ..Default::default() }`, and so a patch release
  never adds a field to one: an option RDKit adds in a patch release waits
  for the crate's next series. That is safe because A3's `<Name>Data`
  starts from a default-constructed C++ object, so a C++ field Rust does
  not know keeps RDKit's default. Structs only the crate builds (`*Outputs`,
  `*Result`, and the other values functions return) are
  `#[non_exhaustive]`, so they can gain fields; a struct callers also
  build, to pass back in (`ExcludedVolume`, `ChargeCorrection`,
  `CustomFeature`, or `BitInfo`, which a caller's metric returns), stays
  exhaustive.
- **A8.** File paths cross the bridge as bytes
  (`Path::as_os_str().as_encoded_bytes()`), never as `&str`.
- **A9.** Python output arguments that carry information nothing else
  returns (fingerprints' `atomBits`, `bitInfo`, and `atomCounts`) come back
  from a `<name>_with_outputs` function, which takes the same parameters and
  returns the result with an outputs struct. An in-out argument is an
  `Option` field of the parameters, and its updated value is in the outputs.
- **A10.** A function or closure the caller supplies is taken as `impl Fn`
  in a free function's arguments and as `&'a dyn Fn` in a parameter
  struct's field or a method's arguments (S3's reason: the path test names
  methods), so a closure that captures will do. Static tables of RDKit's or
  Python's own functions (`desc_list`, `similarity_functions`, the
  fragments) keep `fn` pointers. An alias that names a callback's
  signature is the trait object type (`MetricFunc = dyn Fn(...)`).
- **A11.** Streams a call reads or writes whole are taken by `&mut`: `&mut
  impl Read` and `&mut impl Write` in free functions, a named `R: Read +
  ?Sized` or `W: Write + ?Sized` in methods (with `rust_generics` in the
  manifest, T1). Readers and writers that keep a stream take it by value
  (T1).

## Types

- **N6.** C++ enums cross the bridge as cxx shared enums declared with the C++
  names, so cxx checks every value against RDKit's header. The public enum is
  a separate Rust enum (N4), converted by `rdkit_enum!`. Every public enum is
  `#[non_exhaustive]`, so a later series can add values without a breaking
  change, except namespace enums with no values (`AtomPairsParameters`).
  `from_ffi` panics on a value the Rust enum lacks: the build accepts only
  the RDKit releases whose headers were checked against these enums, so
  such a value is a bug in the crate, not a reachable input. `Error`'s
  variants stay exhaustive, since callbacks build them.
- **N9.** Traits callers are not meant to implement are sealed, by a
  supertrait in a private module (`AsAtom`, `ConvertToNumpyArray`). An open
  trait (`CustomFilterMatcher`, `FeatureLike`) says it is meant to be
  implemented, and a method added to it later has a default body.
- **T1.** Atom and bond indices are `u32`; conformer ids are `i32`, since
  they can be −1. A Python `int` otherwise follows the C++ parameter type.
  `float` → `f64`. `str` → `&str` in, `String` out. File names →
  `impl AsRef<Path>`, or a named `P: AsRef<Path>` in methods and associated
  functions, whose manifest entries set `rust_generics = "<&str>"` so the
  path test can name them (or a list, one per Rust path, with `""` for a
  path that is not generic). Python file objects → `R: Read + Send +
  'r` or `W: Write + Send + 'w` (U12), which the reader or writer borrows
  for its lifetime parameter (`ForwardSDMolSupplier<'r>`, `SDWriter<'w>`);
  one made from a file name has `'static`. Tuples of indices →
  `Vec<u32>`, nested → `Vec<Vec<u32>>`. `None`-able results → `Option`.
  NumPy coordinate arrays → `Vec<Point3D>`; other NumPy matrices → rows,
  `Vec<Vec<f64>>` (or `i32` for integer matrices). `dict` →
  `BTreeMap<String, _>`, unless its order is Python's, by insertion
  (`CalcMolDescriptors`), which a `Vec` of pairs keeps. A returned tuple of
  two or three values → a Rust tuple; of four or more → a named struct whose
  fields take the names RDKit's or Python's documentation gives them
  (A7 decides whether it is `#[non_exhaustive]`).

- **T3.** Nested index lists cross the bridge as `Vec<IndexList>`, and plain
  value types (`Point3D`) as cxx shared structs; both live in the `geometry`
  bridge, and bridges that use them `include!` its generated header and
  request `impl Vec<T> {}` there.
- **T2.** C++ strings cross the bridge as `UniquePtr<CxxString>`, never
  `rust::String` (whose constructor throws on invalid UTF-8). Rust converts
  them with `to_str()`.

## Names

- **N1.** Module paths mirror Python's, each component in snake_case:
  `rdkit.rdBase` → `rdkit_cxx::rd_base`, `rdkit.Chem.rdchem` →
  `rdkit_cxx::chem::rdchem`, `rdkit.Geometry.rdGeometry` →
  `rdkit_cxx::geometry::rd_geometry`. `rdkit_cxx::chem` re-exports its
  submodules the way `rdkit.Chem` does.
- **N2.** Functions, methods, properties, and module attributes take the
  snake_case of the Python name (`snake_case` in `inventory/check.py`) and
  keep prefixes such as `get_`. Module attributes become zero-argument
  functions.
- **N3.** Classes keep the Python class name: `Mol`, not `ROMol`.
- **N4.** Enum values become UpperCamelCase variants (`BondType.SINGLE` →
  `BondType::Single`). An acronym is one word, with only its first letter
  capitalized, in variants as in snake_case names (`FPType.MorganFP` →
  `MorganFp`, `SmallWorldIndexBR` → `SmallWorldIndexBr`); `RDKit` is one
  word too (`RDKitFP` → `RdkitFp`, `GetRDKitFPGenerator` →
  `get_rdkit_fp_generator`), under `name_exception`. Variants of enums the
  crate defines follow the same rule (`FileMolSupplier::Sd`). Flag values
  keep Python's spelling as associated constants
  (`SanitizeFlags::SANITIZE_ALL`).
- **N5.** Exceptions to N2–N4 are recorded in the manifest with
  `name_exception`.
- **N7.** Enum values Python also exports to the module
  (`rdchem.CHI_OTHER`) map to their variants (`ChiralType::ChiOther`) under
  `name_exception`; there are no module-level copies in Rust.
- **N8.** A name that is a Rust keyword takes a trailing underscore:
  `Atom.Match` → `match_`, `StereoInfo.type` → `type_`.
- **F1.** Python flag enums become a `u32` newtype with associated constants,
  bit operators, `bits()`, and `from_bits_retain`; `flag_ops!` writes all but
  the constants. Boost.Python binds RDKit's
  flag enums (`SanitizeFlags`) as plain enums, so the inventory lists their
  values as enum values; they keep Python's spelling under `name_exception`.

## Parity

- **I1.** A commit that binds a Python name adds its manifest entry and
  regenerates `tests/manifest_paths.rs`
  (`python3 inventory/check.py --write-paths tests/manifest_paths.rs`).
- **I2.** `reason` marks a name with no Rust counterpart by design (Python
  protocol methods, pickling hooks, APIs RDKit has deprecated), or one left
  unbound because no check could keep its RDKit code from undefined
  behavior (U1): the reason, which starts "not bound:", names what RDKit
  trusts (`MaeMolSupplier`, the ChemDraw readers). A name that is not bound
  yet has no entry. `overloads_note` explains overloads a binding
  covers with fewer Rust paths, or does not bind yet.
- **I3.** A module joins `complete` once every inventory entry in it is
  mapped; CI then holds it there. Ports of Python-only functionality follow
  R3.

## Ports of Python-only functionality

- **R1.** Functionality RDKit offers only in Python, in its `.py` modules,
  is ported to Rust on top of the bound API, at the paths N1 gives its
  modules: `rdkit.Chem.RegistrationHash` →
  `rdkit_cxx::chem::registration_hash`. Functions in a package's
  `__init__.py` go in a private `init` module that the package's Rust
  module re-exports (`Chem.FindMolChiralCenters` →
  `rdkit_cxx::chem::find_mol_chiral_centers`). Ports follow the naming and
  argument rules (N2–N4, A1–A2), and leave Python's private names (a
  leading underscore) private.
- **R2.** A port computes what Python computes, byte for byte where Python
  returns text or a hash (`RegistrationHash`'s `json.dumps` and SHA-1), and
  its tests pin values from Python's module in the supported RDKit release.
  Where Python changes process-wide state for a call
  (`SetUseLegacyStereoPerception`), the port requires the state instead
  and is `Error::Value` without it (U5).
- **R3.** Ports are measured as bindings are. `inventory/python-<series>.json`
  (`inventory.py --python`) lists what each of RDKit's Python modules defines,
  leaving out tests and what it imports. Each name gets a manifest entry,
  with Rust paths or a reason (I1, I2), and a module joins `complete` the
  same way (I3). A module or package with no Rust counterpart by design,
  such as Python-ecosystem glue (pandas, Jupyter, PyMOL, Windows COM) or a
  legacy subsystem (`sping`, `rdkit.Dbase`), takes one reason under
  `[modules.reasons]` instead, which covers every name in it.
  `inventory/manifest_add.py` merges entries, completions, and module
  reasons into the manifest. `check.py --strict`, which CI runs, fails on a
  module that is neither complete nor covered by a reason. Deprecated code
  is not ported: a module, class, or function RDKit deprecates (one that
  warns on import or use, or is documented as deprecated) takes a reason
  that names what replaces it (I2). A ported class generic over a
  Python value of any type (a container's payload, a callable) gives its
  type parameters defaults, so the manifest can name it and its methods:
  `tests/manifest_paths.rs` names associated items as `<Type>::item`, which
  applies them.
