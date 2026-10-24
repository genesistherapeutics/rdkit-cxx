"""Writes the inventory of RDKit's nanobind Python modules.

The inventory lists every class, function, method, static method, property,
module attribute, and enum value the compiled modules expose, with nanobind's
signature for each overload. check.py measures the manifest against it.

It needs an RDKit build with RDK_BUILD_NANOBIND_WRAPPERS (conda-forge's rdkit
uses Boost.Python); see the README.
"""

import argparse
import enum
import importlib
import importlib.machinery
import inspect
import json
import pathlib
import sys

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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--rdkit-path", type=pathlib.Path, required=True, help="directory holding the rdkit package")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.rdkit_path))
    from rdkit import rdBase

    entries = []
    for name in extension_modules(args.rdkit_path):
        entries.extend(module_entries(importlib.import_module(name)))
    entries.sort(key=lambda entry: entry["python"])
    version = rdBase.rdkitVersion
    args.output.write_text(json.dumps({"rdkit_version": version, "entries": entries}, indent=1) + "\n")
    print(f"wrote {len(entries)} entries for RDKit {version} to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
