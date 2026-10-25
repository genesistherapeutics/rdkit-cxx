# coordgen headers

The headers of Schrödinger's coordgenlibs v3.0.2
(<https://github.com/schrodinger/coordgenlibs>, tag `v3.0.2`), unmodified,
under its BSD 3-Clause license (`LICENSE`).

RDKit builds this version into `libRDKitcoordgen` (`RDKit::coordgen`; see
`External/CoordGen/CMakeLists.txt` in RDKit's sources), and its
`GraphMol/CoordGen.h` includes these headers, but RDKit's packages do not
install them. `wrapper/src/coord_gen.cpp` includes `GraphMol/CoordGen.h`
with this directory on the include path. The pre-commit hooks skip it.
