# Vendored RDKit headers

The wrapper needs a few RDKit headers that conda-forge's `librdkit-dev` does
not install: the drawing code's internal classes (`GraphMol/MolDraw2D`, such
as `DrawMol.h` and `DrawText.h`) and `ML/InfoTheory`. Each directory here
holds RDKit's own copies of them, unmodified, from the release it is named for:

| Directory | Copied from |
| --- | --- |
| `2025.03.6/` | [`Release_2025_03_6`](https://github.com/rdkit/rdkit/tree/Release_2025_03_6) |

They are RDKit's code, under RDKit's BSD 3-clause license
([`../third_party/rdkit/LICENSE`](../third_party/rdkit/LICENSE)). The
pre-commit hooks skip these directories.

## Why one directory per release

The wrapper builds objects of these classes, such as the `DrawMol` that
checks a drawing's layout before RDKit draws it, so their layouts must match
the library the crate links. RDKit changes them within a series:
`DrawShape.h`'s `DrawShapeArrow` gained a member between `Release_2025_03_1`
and `Release_2025_03_6`. Nothing RDKit installs would catch a mismatch, so the
build (`wrapper/cmake/RDKitCxxVersion.cmake`) accepts only an RDKit release
with a directory here, and compiles against that directory.

## Supporting a new RDKit release

RDKit 2025.03.6 is the last release of the 2025.03 series. To support an
earlier one, copy the same files from its tag into a directory named for it
(`2025.03.5/`), check what changed against this one (`git diff --no-index`),
update the wrappers that use any changed class, add the row above, and
release a patch (README, Releases).
