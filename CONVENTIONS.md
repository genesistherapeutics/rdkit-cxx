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

## Types

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
