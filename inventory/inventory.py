# Copyright (c) 2026 Genesis Molecular AI, Inc.
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Writes the inventory of RDKit's compiled Python modules.

The inventory lists every class, function, method, static method, property,
module attribute, and enum value the compiled modules expose, with the
binding's signature for each overload: Boost.Python's, as conda-forge
builds RDKit, or nanobind's. check.py measures the manifest against it.

With --python it writes the inventory of RDKit's Python (.py) modules instead:
the functions, classes, methods, properties, enum values, and module
attributes each defines, leaving out tests and what it imports. Modules that
do not import (they need IPython, Qt, or Windows) are read from their source.

It needs an RDKit build with Python wrappers: conda-forge's rdkit, whose
Boost.Python modules are the crate's reference, or a build with
RDK_BUILD_NANOBIND_WRAPPERS; see the README.
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

NB_CALLABLE_TYPES = {"nb_func", "nb_method"}

# A Boost.Python overload's Python signature, the line before its "C++
# signature" line in the docstring.
BOOST_SIGNATURE = re.compile(r"^\s*(\w+)\((.*)\)\s*->\s*(.*?)\s*:\s*$")
# The address in a default value's repr, which changes from run to run.
ADDRESS = re.compile(r" at 0x[0-9a-f]+>")


def is_boost_callable(obj) -> bool:
    return type(obj).__name__ == "function" and type(obj).__module__ == "Boost.Python"


def is_compiled_callable(obj) -> bool:
    """A compiled function or method, of either binding."""
    return type(obj).__name__ in NB_CALLABLE_TYPES or is_boost_callable(obj)


def boost_signatures(obj) -> list[str]:
    lines = (getattr(obj, "__doc__", None) or "").splitlines()
    sigs = []
    for i, line in enumerate(lines):
        if "C++ signature" in line:
            sig = next((lines[j].strip() for j in range(i - 1, -1, -1) if BOOST_SIGNATURE.match(lines[j])), "?")
            sigs.append(ADDRESS.sub(">", sig))
    return sigs


def signatures(obj) -> list[str]:
    if is_boost_callable(obj):
        return boost_signatures(obj)
    return [overload[0] for overload in getattr(obj, "__nb_signature__", ())]


def is_boost_enum(cls) -> bool:
    return any(base.__name__ == "enum" and base.__module__ == "Boost.Python" for base in cls.__mro__)


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


def is_boost_reduce(name: str, obj) -> bool:
    """Boost.Python's generic __reduce__, one undocumented function every
    class holds, which pickles through the class's __getinitargs__ and
    __getstate__ (listed when it has them) or raises that pickling is not
    enabled."""
    return name == "__reduce__" and is_boost_callable(obj) and obj.__doc__ is None


def is_public(name: str, obj) -> bool:
    """Public names, plus the protocol methods the binding defines."""
    if name.startswith("__") and name.endswith("__"):
        return is_compiled_callable(obj) and not is_boost_reduce(name, obj)
    return not name.startswith("_")


def class_entries(cls, qualified: str) -> list[dict]:
    entries = [{"python": qualified, "kind": "class", "overloads": []}]
    if issubclass(cls, enum.Enum):
        kind = "flag_value" if issubclass(cls, enum.Flag) else "enum_value"
        for name in cls.__members__:
            entries.append({"python": f"{qualified}.{name}", "kind": kind, "overloads": []})
        return entries
    if is_boost_enum(cls):
        for name in cls.names:
            entries.append({"python": f"{qualified}.{name}", "kind": "enum_value", "overloads": []})
        return entries
    for name, obj in vars(cls).items():
        if not is_public(name, obj):
            continue
        python = f"{qualified}.{name}"
        if isinstance(obj, staticmethod):
            entries.append({"python": python, "kind": "static_method", "overloads": signatures(obj.__func__)})
        elif isinstance(obj, property) or type(obj).__name__ in {"nb_static_property", "property"}:
            entries.append({"python": python, "kind": "property", "overloads": []})
        elif is_compiled_callable(obj):
            entries.append({"python": python, "kind": "method", "overloads": signatures(obj)})
        elif inspect.isfunction(obj):
            # A method RDKit's Python adds to the class (Mol.GetAtoms under
            # Boost.Python).
            entries.append({"python": python, "kind": "method", "overloads": python_signature(name, obj)})
        elif inspect.isclass(obj):
            if obj.__module__ == cls.__module__:
                entries.extend(class_entries(obj, python))
        elif not callable(obj):
            # A plain class attribute (AtomPairsParameters' ints).
            entries.append({"python": python, "kind": "property", "overloads": []})
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
        elif is_compiled_callable(obj):
            if getattr(obj, "__module__", module.__name__) == module.__name__:
                entries.append({"python": python, "kind": "module_function", "overloads": signatures(obj)})
        elif isinstance(obj, enum.Enum) or is_boost_enum(type(obj)):
            # Enum values exported to the module (export_values).
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
    """The dotted names of the Python (.py) modules under root/rdkit, without tests."""
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
    """What a Python module defines: functions and classes whose
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
    parser.add_argument("--python", action="store_true", help="inventory the Python (.py) modules instead")
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
