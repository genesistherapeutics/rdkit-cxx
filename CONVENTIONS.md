# Conventions

Every rule here binds all later wrapper code. A change that cannot follow a
rule changes the rule first.

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
  `XMut`s into one molecule would alias.
- **S1.** `Mol` and `RWMol` are `Send` and not `Sync`. RDKit writes computed
  properties through `const` molecules (`MolToSmiles` sets
  `_smilesAtomOutputOrder` in `SmilesWrite.cpp`), so two threads reading one
  molecule race. Parallel code gives each thread its own molecule.

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
  `__init__` is `Clone`; the first other `__init__` is `new`; the rest are
  `from_<first parameter>`.
- **A6.** `_mut` comes last: `get_conformer_with_id_mut`.

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
