# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Checks the parity manifest against the inventories.

The inventories, written by inventory.py, list every class, function, method,
property, and enum value RDKit's compiled modules expose (Boost.Python's, as
conda-forge builds them; inventory-<series>.json) and that its Python (`.py`)
modules define
(python-<series>.json). The manifest (manifest.toml) maps each to Rust paths,
or to a reason it has none, given for the name or for its whole module
(`[modules.reasons]`, CONVENTIONS.md R3). This script reports manifest entries
that are malformed, misnamed (N2-N5), or missing from the inventories, and
inventory entries left unmapped in modules marked complete (I3) or in no
module at all. With --write-paths it also writes a Rust test that names every
mapped path, so `cargo test` fails when a manifest path does not exist.
"""

import argparse
import json
import pathlib
import re
import sys
import tomllib

HERE = pathlib.Path(__file__).parent
DEFAULT_INVENTORIES = [HERE / "inventory-2024.03.json", HERE / "python-2024.03.json"]
DEFAULT_MANIFEST = HERE / "manifest.toml"

# addlicense's BSD header, so the generated file matches what it writes.
LICENSE_HEADER = [
    "// Copyright (c) 2026 Steven Kearnes",
    "// All rights reserved.",
    "// Use of this source code is governed by a BSD-style",
    "// license that can be found in the LICENSE file.",
]

CALLABLE_KINDS = {"module_function", "method", "static_method"}
SNAKE_CASE_KINDS = CALLABLE_KINDS | {"module_attribute", "property"}
INIT_NAMES = {"new", "default", "clone"}


def snake_case(name: str) -> str:
    name = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name)
    name = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", name)
    return name.lower()


def upper_camel_case(name: str) -> str:
    # ALL_CAPS names split only at underscores, so SP3D2 becomes Sp3d2.
    words = name.split("_") if name.isupper() else snake_case(name).split("_")
    return "".join(word.capitalize() for word in words)


def module_of(python: str, classes: set[str]) -> str:
    """The module an inventory name belongs to: its longest prefix that is not a class."""
    parent = python.rsplit(".", 1)[0]
    while parent in classes:
        parent = parent.rsplit(".", 1)[0]
    return parent


def expected_name(name: str, kind: str, rust_kind: str | None) -> str:
    if rust_kind == "field" or kind in SNAKE_CASE_KINDS:
        return snake_case(name)
    if kind == "enum_value":
        return upper_camel_case(name)
    return name  # classes and flag values keep Python's spelling


def name_problem(python: str, kind: str, item: dict, rust: str) -> str | None:
    name = python.rsplit(".", 1)[1]
    segment = rust.rsplit("::", 1)[-1]
    if name == "__init__":
        if segment in INIT_NAMES or segment.startswith(("from_", "new_with_")):
            return None
        return f"{python}: {rust} is not new, new_with_*, default, clone, or from_* (A5)"
    if name.startswith("__") or "name_exception" in item:
        return None
    expected = expected_name(name, kind, item.get("rust_kind"))
    if segment == expected or segment.startswith(expected + "_"):
        return None
    if kind == "class" and segment in (expected + "Ref", expected + "Mut"):
        return None  # views of the class (W2)
    return f"{python}: Rust name {segment!r} is neither {expected!r} nor {expected + '_'!r}... (N2-N4)"


def module_reason(module: str, reasons: dict[str, str]) -> str | None:
    """The reason given for a module or for a package holding it."""
    while module:
        if module in reasons:
            return reasons[module]
        module = module.rpartition(".")[0]
    return None


def merge(inventories: list[dict]) -> dict:
    versions = {inventory["rdkit_version"] for inventory in inventories}
    if len(versions) != 1:
        raise SystemExit(f"the inventories are for different RDKit versions: {sorted(versions)}")
    entries = [entry for inventory in inventories for entry in inventory["entries"]]
    return {"rdkit_version": versions.pop(), "entries": entries}


def check(inventory: dict, manifest: dict) -> list[str]:
    entries = {e["python"]: e for e in inventory["entries"]}
    classes = {name for name, e in entries.items() if e["kind"] == "class"}
    problems = []
    mapped = set()
    for item in manifest.get("entry", []):
        python = item["python"]
        if python in mapped:
            problems.append(f"{python}: listed twice")
            continue
        mapped.add(python)
        entry = entries.get(python)
        if entry is None:
            problems.append(f"{python}: not in the inventory")
            continue
        if ("rust" in item) == ("reason" in item):
            problems.append(f"{python}: needs exactly one of rust and reason")
            continue
        if "reason" in item:
            continue
        rust = item["rust"]
        if not rust:
            problems.append(f"{python}: rust is empty")
            continue
        for path in rust:
            if problem := name_problem(python, entry["kind"], item, path):
                problems.append(problem)
        generics = item.get("rust_generics")
        if isinstance(generics, list) and len(generics) != len(rust):
            problems.append(f"{python}: rust_generics lists {len(generics)} paths' generics for {len(rust)} paths")
        overloads = len(entry["overloads"])
        if len(rust) < overloads and "overloads_note" not in item:
            problems.append(
                f"{python}: {overloads} overloads but {len(rust)} Rust paths; bind each or add overloads_note"
            )
    modules = manifest.get("modules", {})
    complete = set(modules.get("complete", []))
    reasons = modules.get("reasons", {})
    known = {module_of(name, classes) for name in entries}
    for module in reasons:
        if not any(m == module or m.startswith(module + ".") for m in known):
            problems.append(f"{module}: has a reason but is not in the inventory")
        if module in complete:
            problems.append(f"{module}: has a reason and is also marked complete")
    for name in entries:
        module = module_of(name, classes)
        if name not in mapped and module in complete:
            problems.append(f"{name}: unmapped in complete module {module}")
    return problems


def path_statement(path: str, rust_kind: str | None, rust_generics: str | None = None) -> str:
    *parents, last = path.split("::")
    owner = "::".join(parents)
    if rust_kind == "field":
        return f"let _ = |value: &{owner}| {{ let _ = &value.{last}; }};"
    if parents[-1][0].isupper():
        # An associated function, method, constant, or enum variant, named
        # through a qualified path, which applies a generic type's default
        # parameters (R3). A generic item is named with the manifest's
        # instantiation (CONVENTIONS.md T1).
        return f"let _ = <{owner}>::{last}{'::' + rust_generics if rust_generics else ''};"
    if last[0].isupper():
        return f"let _: Option<{path}> = None;"
    # A free function. `use` also covers `impl Trait` parameters, whose types
    # `let _ = f;` cannot infer.
    return f"use {path} as _;"


def generics_of(item: dict) -> list[str | None]:
    """Each Rust path's generics: rust_generics names one instantiation for
    every path, or a list of them, one per path, "" for a path that is not
    generic."""
    generics = item.get("rust_generics")
    paths = item.get("rust", [])
    if isinstance(generics, list):
        return [g or None for g in generics]
    return [generics] * len(paths)


def uncovered_modules(inventory: dict, manifest: dict) -> list[str]:
    """Inventory modules neither complete nor covered by a reason (R3)."""
    entries = {e["python"]: e for e in inventory["entries"]}
    classes = {name for name, e in entries.items() if e["kind"] == "class"}
    modules = manifest.get("modules", {})
    complete = set(modules.get("complete", []))
    reasons = modules.get("reasons", {})
    found = {module_of(name, classes) for name in entries}
    return sorted(m for m in found if m not in complete and module_reason(m, reasons) is None)


def paths_source(manifest: dict) -> str:
    lines = [
        *LICENSE_HEADER,
        "",
        "// @generated by inventory/check.py from inventory/manifest.toml. Do not edit.",
        "//",
        "// Names every Rust path the manifest maps a Python name to, so the build",
        "// fails if one does not exist.",
        "",
        "#[test]",
        "#[rustfmt::skip]",
        "#[allow(unused_imports)]",
        "fn manifest_paths_exist() {",
    ]
    for item in sorted(manifest.get("entry", []), key=lambda item: item["python"]):
        for path, generics in zip(item.get("rust", []), generics_of(item), strict=True):
            lines.append(f"    {path_statement(path, item.get('rust_kind'), generics)}")
    lines.append("}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--inventory", type=pathlib.Path, action="append", help="an inventory to check against (default: both)"
    )
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--write-paths", type=pathlib.Path, help="write the Rust path test here")
    parser.add_argument(
        "--strict", action="store_true", help="also fail on modules neither complete nor covered by a reason"
    )
    args = parser.parse_args()
    inventory = merge([json.loads(path.read_text()) for path in args.inventory or DEFAULT_INVENTORIES])
    manifest = tomllib.loads(args.manifest.read_text())
    problems = check(inventory, manifest)
    if args.strict:
        problems += [
            f"{module}: neither complete nor covered by a reason" for module in uncovered_modules(inventory, manifest)
        ]
    for problem in problems:
        print(problem, file=sys.stderr)
    if args.write_paths:
        args.write_paths.write_text(paths_source(manifest))
    entries = {e["python"]: e for e in inventory["entries"]}
    classes = {name for name, e in entries.items() if e["kind"] == "class"}
    items = [item for item in manifest.get("entry", []) if item["python"] in entries]
    with_reason = sum("reason" in item for item in items)
    listed = {item["python"] for item in items}
    reasons = manifest.get("modules", {}).get("reasons", {})
    by_module = sum(
        name not in listed and module_reason(module_of(name, classes), reasons) is not None for name in entries
    )
    complete = len(manifest.get("modules", {}).get("complete", []))
    print(
        f"mapped {len(items) + by_module} of {len(entries)} entries ({with_reason} with a reason, "
        f"{by_module} by their module's); complete modules: {complete}, modules with a reason: {len(reasons)}"
    )
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
