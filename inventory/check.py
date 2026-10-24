# Copyright (c) 2026 Steven Kearnes All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Checks the parity manifest against the nanobind inventory.

The inventory (inventory-<series>.json, written by inventory.py) lists every
class, function, method, property, and enum value RDKit's nanobind modules
expose. The manifest (manifest.toml) maps each to Rust paths, or to a reason it
has none. This script reports manifest entries that are malformed, misnamed
(CONVENTIONS.md N2-N5), or missing from the inventory, and inventory entries
left unmapped in modules marked complete (I3). With --write-paths it also
writes a Rust test that names every mapped path, so `cargo test` fails when a
manifest path does not exist.
"""

import argparse
import json
import pathlib
import re
import sys
import tomllib

HERE = pathlib.Path(__file__).parent
DEFAULT_INVENTORY = HERE / "inventory-2026.09.json"
DEFAULT_MANIFEST = HERE / "manifest.toml"

# addlicense's BSD header, so the generated file matches what it writes.
LICENSE_HEADER = [
    "// Copyright (c) 2026 Steven Kearnes All rights reserved.",
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
        overloads = len(entry["overloads"])
        if len(rust) < overloads and "overloads_note" not in item:
            problems.append(
                f"{python}: {overloads} overloads but {len(rust)} Rust paths; bind each or add overloads_note"
            )
    complete = set(manifest.get("modules", {}).get("complete", []))
    for name, entry in entries.items():
        if name not in mapped and module_of(name, classes) in complete:
            problems.append(f"{name}: unmapped in complete module {module_of(name, classes)}")
    return problems


def path_statement(path: str, rust_kind: str | None) -> str:
    *parents, last = path.split("::")
    owner = "::".join(parents)
    if rust_kind == "field":
        return f"let _ = |value: &{owner}| {{ let _ = &value.{last}; }};"
    if parents[-1][0].isupper():
        # An associated function, method, constant, or enum variant.
        return f"let _ = {path};"
    if last[0].isupper():
        return f"let _: Option<{path}> = None;"
    # A free function. `use` also covers `impl Trait` parameters, whose types
    # `let _ = f;` cannot infer.
    return f"use {path} as _;"


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
        for path in item.get("rust", []):
            lines.append(f"    {path_statement(path, item.get('rust_kind'))}")
    lines.append("}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--inventory", type=pathlib.Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--manifest", type=pathlib.Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--write-paths", type=pathlib.Path, help="write the Rust path test here")
    args = parser.parse_args()
    inventory = json.loads(args.inventory.read_text())
    manifest = tomllib.loads(args.manifest.read_text())
    problems = check(inventory, manifest)
    for problem in problems:
        print(problem, file=sys.stderr)
    if args.write_paths:
        args.write_paths.write_text(paths_source(manifest))
    names = {e["python"] for e in inventory["entries"]}
    items = [item for item in manifest.get("entry", []) if item["python"] in names]
    with_reason = sum("reason" in item for item in items)
    complete = len(manifest.get("modules", {}).get("complete", []))
    print(
        f"mapped {len(items)} of {len(names)} entries ({with_reason} with a reason); complete modules: {complete}"
    )
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
