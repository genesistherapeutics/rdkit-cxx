# Conventions

Every rule here binds all later wrapper code. A change that cannot follow a
rule changes the rule first.

## Files

- **L1.** Every source file starts with addlicense's BSD header. Add it with
  `addlicense -c "Steven Kearnes" -l bsd -ignore "target/**" .`
  (`go install github.com/google/addlicense@v1.2.0`); CI checks it.
  Generators write the same header themselves (`inventory/check.py`).

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

- **E1.** A public function returns `Result` when the nanobind call can raise
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
  | `ValueErrorException`, `RDKit::ConformerException` | `Value` |
  | `RDKit::BadFileException` | `Os` |
  | `RDKit::SmilesParse::SmilesParseException`, `RDKit::FileParseException` | `Parse` |
  | `Invar::Invariant` (with `toUserString()`), any other exception | `Runtime` |

- **E5.** Every bridge's `include!` list starts with `"rdkit-cxx/errors.h"`,
  which defines the `rust::behavior::trycatch` that encodes the exception.

## Ownership and threads

- **W1.** Every Rust molecule owns a C++ `RDKit::RWMol`, so in-place
  operations can take `RWMol&` without a cast.
- **W2.** Shared views are `XRef<'m>` (`Copy`, holding `&'m Mol` and
  `&'m ffi::X`). Exclusive views are `XMut<'m>` (holding `&'m mut Mol` and
  an index). An owned `X` exists only when Python has a constructor.
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
  the call returns.
- **S1.** `Mol` and `RWMol` are `Send` and not `Sync`. RDKit writes computed
  properties through `const` molecules (`MolToSmiles` sets
  `_smilesAtomOutputOrder` in `SmilesWrite.cpp`), so two threads reading one
  molecule race. Parallel code gives each thread its own molecule.

- **S2.** The crate needs an RDKit built with `RDK_BUILD_THREADSAFE_SSS`, as
  conda-forge's is: process-wide settings such as the default pickle
  properties are safe to set from safe Rust only because RDKit guards them
  with a mutex. The wrappers refuse to compile otherwise.
- **S3.** Module-level functions that take any atom are generic over the
  sealed traits `AsAtom` (`Atom`, `AtomRef`, `AtomMut`) and `AsAtomMut`
  (`Atom`, `AtomMut`). Methods take `&dyn AsAtom` instead, so they are not
  generic and `tests/manifest_paths.rs` can name them.

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

## Documentation

- **D1.** Doc examples are `no_run` or `compile_fail`; behavior is tested
  under `tests/`. Each public item's doc comment names the Python object it
  mirrors.

## Arguments and overloads

- **A1.** Required Python parameters become positional Rust parameters, in
  Python's order.
- **A2.** A Python signature with optional parameters is bound as `name`
  (using Python's defaults) plus, for exactly one optional parameter,
  `name_with_<param>`, or for two or more, `name_with_params(&<Name>Params)`.
  `<Name>Params` is named for the Python function or method, has the optional
  parameters as snake_case fields, and implements `Default` with Python's
  defaults, pinned by a test.
- **A3.** When Python also accepts a params object (an overload taking, say,
  `SmilesParserParams`), `name_with_params` takes the Rust mirror of that
  class instead. Keyword arguments the class lacks become extra fields on it.
  The mirror's `Default` reads a default-constructed C++ object, so it cannot
  drift.
- **A4.** An optional parameter that changes the return type gets its own
  function, `name_<param>` (`GetMolFrags(asMols=True)` →
  `get_mol_frags_as_mols`), and does not count toward A2. Python output
  arguments (lists the call fills in) are not bound; the return value carries
  that information.
- **A5.** Overloads that differ by an argument's type: the first in
  nanobind's order keeps the name, and the others append `_<type>` in
  snake_case. Constructors: a zero-argument `__init__` is `Default`; a copy
  `__init__` is `Clone`; the first other `__init__` is `new` (with
  `new_with_params` for its optional parameters, per A2); the rest are
  `from_<first parameter>`, or a `From` impl when the parameter is another
  bound class (`RWMol(mol)` → `From<Mol>`), which the manifest records in
  `overloads_note`. nanobind's `__new__` (from `nb::new_`) repeats
  `__init__`'s overloads and gets a `reason`.
- **A6.** `_mut` comes last: `get_conformer_with_id_mut`.
- **A7.** Parameter structs are plain Rust structs with public fields. A
  mirror of a Python class crosses the bridge as a cxx shared struct
  `<Name>Data` with the same fields (maps as parallel key and value vectors),
  converted by a private `to_ffi`; its `Default` comes from a wrapper that
  returns the C++ defaults. Structs that A2 defines, with no Python class,
  write Python's defaults in Rust.
- **A8.** File paths cross the bridge as bytes
  (`Path::as_os_str().as_encoded_bytes()`), never as `&str`.

## Types

- **N6.** C++ enums cross the bridge as cxx shared enums declared with the C++
  names, so cxx checks every value against RDKit's header. The public enum is
  a separate Rust enum (N4), converted by `rdkit_enum!`.
- **T1.** Atom and bond indices are `u32`; conformer ids are `i32`, since
  they can be −1. A Python `int` otherwise follows the C++ parameter type.
  `float` → `f64`. `str` → `&str` in, `String` out. File names →
  `impl AsRef<Path>`. Tuples of indices → `Vec<u32>`, nested →
  `Vec<Vec<u32>>`. `None`-able results → `Option`. NumPy coordinate arrays →
  `Vec<Point3D>`. `dict` → `BTreeMap<String, _>`.

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
  `BondType::Single`). Flag values keep Python's spelling as associated
  constants (`SanitizeFlags::SANITIZE_ALL`).
- **N5.** Exceptions to N2–N4 are recorded in the manifest with
  `name_exception`.
- **N7.** Enum values nanobind also exports to the module
  (`rdchem.CHI_OTHER`) map to their variants (`ChiralType::ChiOther`) under
  `name_exception`; there are no module-level copies in Rust.
- **N8.** A name that is a Rust keyword takes a trailing underscore:
  `Atom.Match` → `match_`, `StereoInfo.type` → `type_`.
- **F1.** Python flag enums become a `u32` newtype with associated constants,
  bit operators, `bits()`, and `from_bits_retain`. nanobind binds RDKit's
  flag enums (`SanitizeFlags`) as arithmetic `IntEnum`s, not `enum.Flag`, so
  the inventory lists their values as enum values; they keep Python's
  spelling under `name_exception`.

## Parity

- **I1.** A commit that binds a Python name adds its manifest entry and
  regenerates `tests/manifest_paths.rs`
  (`python3 inventory/check.py --write-paths tests/manifest_paths.rs`).
- **I2.** `reason` marks a name with no Rust counterpart by design (Python
  protocol methods, pickling hooks, APIs RDKit has deprecated). A name that is
  not bound yet has no entry. `overloads_note` explains overloads a binding
  covers with fewer Rust paths, or does not bind yet.
- **I3.** A module joins `complete` once every inventory entry in it is
  mapped; CI then holds it there.
