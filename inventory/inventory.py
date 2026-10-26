# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Writes the inventory of RDKit's nanobind Python modules.

The inventory lists every class, function, method, static method, property,
module attribute, and enum value the compiled modules expose, with nanobind's
signature for each overload. check.py measures the manifest against it.

With --python it writes the inventory of RDKit's pure-Python modules instead:
the functions, classes, methods, properties, enum values, and module
attributes each defines, leaving out tests and what it imports. Modules that
do not import (they need IPython, Qt, or Windows) are read from their source.

It needs an RDKit build with RDK_BUILD_NANOBIND_WRAPPERS (conda-forge's rdkit
uses Boost.Python); see the README.
"""

import argparse
import ast
import contextlib
import enum
import importlib
import importlib.machinery
import inspect
import io
import json
import logging
import pathlib
import re
import sys
import warnings

CALLABLE_TYPES = {"nb_func", "nb_method"}


def is_nb_callable(obj) -> bool:
    return type(obj).__name__ in CALLABLE_TYPES


def signatures(obj) -> list[str]:
    return [overload[0] for overload in getattr(obj, "__nb_signature__", ())]


def extension_modules(root: pathlib.Path) -> list[str]:
    """The dotted names of the compiled modules under root/rdkit."""
    # Longest first, so foo.cpython-312-x86_64-linux-gnu.so is not taken as foo.cpython-312-x86_64-linux-gnu.
    suffixes = sorted(importlib.machinery.EXTENSION_SUFFIXES, key=len, reverse=True)
    names = set()
    for path in (root / "rdkit").rglob("*"):
        suffix = next((s for s in suffixes if path.name.endswith(s)), None)
        if suffix and path.is_file():
            relative = path.relative_to(root).with_name(path.name[: -len(suffix)])
            names.add(".".join(relative.parts))
    return sorted(names)


def is_public(name: str, obj) -> bool:
    """Public names, plus the protocol methods nanobind defines."""
    if name.startswith("__") and name.endswith("__"):
        return is_nb_callable(obj)
    return not name.startswith("_")


def class_entries(cls, qualified: str) -> list[dict]:
    entries = [{"python": qualified, "kind": "class", "overloads": []}]
    if issubclass(cls, enum.Enum):
        kind = "flag_value" if issubclass(cls, enum.Flag) else "enum_value"
        for name in cls.__members__:
            entries.append({"python": f"{qualified}.{name}", "kind": kind, "overloads": []})
        return entries
    for name, obj in vars(cls).items():
        if not is_public(name, obj):
            continue
        python = f"{qualified}.{name}"
        if isinstance(obj, staticmethod):
            entries.append({"python": python, "kind": "static_method", "overloads": signatures(obj.__func__)})
        elif isinstance(obj, property) or type(obj).__name__ == "nb_static_property":
            entries.append({"python": python, "kind": "property", "overloads": []})
        elif is_nb_callable(obj):
            entries.append({"python": python, "kind": "method", "overloads": signatures(obj)})
        elif inspect.isclass(obj) and obj.__module__ == cls.__module__:
            entries.extend(class_entries(obj, python))
    return entries


def module_entries(module) -> list[dict]:
    entries = []
    for name, obj in vars(module).items():
        if not is_public(name, obj) or inspect.ismodule(obj):
            continue
        python = f"{module.__name__}.{name}"
        if inspect.isclass(obj):
            # Classes defined here; others are imports.
            if obj.__module__ == module.__name__ and obj.__qualname__ == name:
                entries.extend(class_entries(obj, python))
        elif is_nb_callable(obj):
            if getattr(obj, "__module__", module.__name__) == module.__name__:
                entries.append({"python": python, "kind": "module_function", "overloads": signatures(obj)})
        elif isinstance(obj, enum.Enum):
            # Enum values exported to the module (nanobind's export_values).
            if type(obj).__module__ == module.__name__:
                entries.append({"python": python, "kind": "module_attribute", "overloads": []})
        elif not callable(obj):
            entries.append({"python": python, "kind": "module_attribute", "overloads": []})
    return entries


def is_test(part: str) -> bool:
    return (
        part.startswith(("UnitTest", "test", "Test"))
        or part.endswith("_test")
        or part in {"conftest", "nbtests", "test_data"}
    )


def python_modules(root: pathlib.Path) -> list[str]:
    """The dotted names of the pure-Python modules under root/rdkit, without tests."""
    names = set()
    for path in (root / "rdkit").rglob("*.py"):
        parts = path.relative_to(root).with_suffix("").parts
        if any(is_test(part) for part in parts):
            continue
        names.add(".".join(parts[:-1] if parts[-1] == "__init__" else parts))
    return sorted(names)


def module_path(root: pathlib.Path, name: str) -> pathlib.Path:
    path = root.joinpath(*name.split("."))
    return path / "__init__.py" if path.is_dir() else path.with_suffix(".py")


def assigned_names(tree: ast.Module) -> list[str]:
    """Names a module's top level assigns, in order, and those its functions
    declare global (Descriptors.descList)."""
    names = []
    for node in tree.body:
        targets = (
            node.targets if isinstance(node, ast.Assign) else [node.target] if isinstance(node, ast.AnnAssign) else []
        )
        for target in targets:
            names.extend(n.id for n in ast.walk(target) if isinstance(n, ast.Name))
        if isinstance(node, ast.FunctionDef):
            for sub in ast.walk(node):
                if isinstance(sub, ast.Global):
                    names.extend(sub.names)
    return names


def is_logger(obj) -> bool:
    """A module's own logger, which is how it logs, not what it offers."""
    return isinstance(obj, logging.Logger) or type(obj).__module__ == "rdkit.RDLogger"


def python_signature(name: str, obj) -> list[str]:
    """Python's signature, named for the attribute (a lambda's is <lambda>),
    without the addresses of default values' reprs."""
    try:
        signature = str(inspect.signature(obj))
    except (TypeError, ValueError):
        signature = "(...)"
    return [f"def {name}{re.sub(r' at 0x[0-9a-f]+', '', signature)}"]


def python_class_entries(cls, qualified: str) -> list[dict]:
    entries = [{"python": qualified, "kind": "class", "overloads": []}]
    if issubclass(cls, enum.Enum):
        kind = "flag_value" if issubclass(cls, enum.Flag) else "enum_value"
        for name in cls.__members__:
            entries.append({"python": f"{qualified}.{name}", "kind": kind, "overloads": []})
        return entries
    for name, obj in vars(cls).items():
        dunder = name.startswith("__") and name.endswith("__")
        if name.startswith("_") and not (dunder and inspect.isfunction(obj)):
            continue
        python = f"{qualified}.{name}"
        if isinstance(obj, (staticmethod, classmethod)):
            entries.append(
                {"python": python, "kind": "static_method", "overloads": python_signature(name, obj.__func__)}
            )
        elif isinstance(obj, property):
            entries.append({"python": python, "kind": "property", "overloads": []})
        elif inspect.isfunction(obj):
            entries.append({"python": python, "kind": "method", "overloads": python_signature(name, obj)})
        elif inspect.isclass(obj):
            if obj.__module__ == cls.__module__ and obj.__qualname__ == f"{cls.__qualname__}.{name}":
                entries.extend(python_class_entries(obj, python))
        elif not callable(obj):
            entries.append({"python": python, "kind": "property", "overloads": []})
    return entries


def python_module_entries(module, tree: ast.Module) -> list[dict]:
    """What a pure-Python module defines: functions and classes whose
    __module__ is the module (so also those it makes at import), and the
    names its top level assigns that hold neither."""
    entries = []
    assigned = set(assigned_names(tree))
    for name, obj in vars(module).items():
        if name.startswith("_") or inspect.ismodule(obj):
            continue
        python = f"{module.__name__}.{name}"
        if inspect.isclass(obj):
            if obj.__module__ == module.__name__ and obj.__qualname__ == name:
                entries.extend(python_class_entries(obj, python))
        elif inspect.isfunction(obj):
            if obj.__module__ == module.__name__:
                entries.append({"python": python, "kind": "module_function", "overloads": python_signature(name, obj)})
        elif name in assigned and not callable(obj) and not is_logger(obj):
            entries.append({"python": python, "kind": "module_attribute", "overloads": []})
    return entries


def source_entries(name: str, tree: ast.Module) -> list[dict]:
    """As python_module_entries, from the source alone, for a module that
    does not import."""
    entries = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and not node.name.startswith("_"):
            entries.append(
                {"python": f"{name}.{node.name}", "kind": "module_function", "overloads": [f"def {node.name}(...)"]}
            )
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            entries.append({"python": f"{name}.{node.name}", "kind": "class", "overloads": []})
            for sub in node.body:
                dunder = sub.name.startswith("__") if isinstance(sub, ast.FunctionDef) else False
                if isinstance(sub, ast.FunctionDef) and (not sub.name.startswith("_") or dunder):
                    entries.append(
                        {
                            "python": f"{name}.{node.name}.{sub.name}",
                            "kind": "method",
                            "overloads": [f"def {sub.name}(...)"],
                        }
                    )
    for target in assigned_names(tree):
        if not target.startswith("_"):
            entries.append({"python": f"{name}.{target}", "kind": "module_attribute", "overloads": []})
    unique = {entry["python"]: entry for entry in entries}
    return list(unique.values())


def python_entries(root: pathlib.Path) -> tuple[list[dict], list[str]]:
    entries, unimported = [], []
    warnings.simplefilter("ignore")
    for name in python_modules(root):
        tree = ast.parse(module_path(root, name).read_text(encoding="utf-8", errors="replace"))
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                module = importlib.import_module(name)
        except Exception:
            unimported.append(name)
            entries.extend(source_entries(name, tree))
            continue
        entries.extend(python_module_entries(module, tree))
    return entries, unimported


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--rdkit-path", type=pathlib.Path, required=True, help="directory holding the rdkit package")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--python", action="store_true", help="inventory the pure-Python modules instead")
    args = parser.parse_args()
    sys.path.insert(0, str(args.rdkit_path))
    from rdkit import rdBase

    entries = []
    if args.python:
        entries, unimported = python_entries(args.rdkit_path)
        for name in unimported:
            print(f"{name} does not import; read from its source")
    else:
        for name in extension_modules(args.rdkit_path):
            entries.extend(module_entries(importlib.import_module(name)))
    entries.sort(key=lambda entry: entry["python"])
    version = rdBase.rdkitVersion
    args.output.write_text(json.dumps({"rdkit_version": version, "entries": entries}, indent=1) + "\n")
    print(f"wrote {len(entries)} entries for RDKit {version} to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
