# Copyright (c) 2026 Genesis Molecular AI
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Merges entries into manifest.toml, keeping it sorted.

ENTRIES is a JSON list of entry tables, as manifest.toml spells them; an entry
replaces the one with the same `python` name. --complete marks modules complete
(CONVENTIONS.md I3), and --reason gives a module's reason (R3) as MODULE=REASON.

    python3 inventory/manifest_add.py entries.json --complete rdkit.Chem.Lipinski
"""

import argparse
import json
import pathlib
import sys
import tomllib

HERE = pathlib.Path(__file__).parent
KEYS = ["python", "rust", "reason", "rust_kind", "rust_generics", "name_exception", "overloads_note"]


def value(v) -> str:
    if isinstance(v, list):
        return "[" + ", ".join(value(x) for x in v) + "]"
    return json.dumps(v, ensure_ascii=False)


def render(header: str, modules: dict, entries: dict[str, dict]) -> str:
    lines = [header.rstrip("\n"), "", "[modules]", "complete = " + value(sorted(modules.get("complete", [])))]
    reasons = modules.get("reasons", {})
    if reasons:
        lines += ["", "[modules.reasons]"]
        lines += [f"{value(name)} = {value(reasons[name])}" for name in sorted(reasons)]
    for name in sorted(entries):
        unknown = set(entries[name]) - set(KEYS)
        if unknown:
            raise SystemExit(f"{name}: unknown keys {sorted(unknown)}")
        lines += ["", "[[entry]]"] + [f"{k} = {value(entries[name][k])}" for k in KEYS if k in entries[name]]
    return "\n".join(lines).lstrip("\n") + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("entries", type=pathlib.Path, nargs="?", help="JSON list of entries to add or replace")
    parser.add_argument("--manifest", type=pathlib.Path, default=HERE / "manifest.toml")
    parser.add_argument("--complete", action="append", default=[], metavar="MODULE")
    parser.add_argument("--reason", action="append", default=[], metavar="MODULE=REASON")
    args = parser.parse_args()
    text = args.manifest.read_text()
    header = text[: text.index("[modules]")]
    data = tomllib.loads(text)
    modules = data.get("modules", {})
    entries = {e["python"]: e for e in data.get("entry", [])}
    for entry in json.loads(args.entries.read_text()) if args.entries else []:
        entries[entry["python"]] = entry
    modules["complete"] = sorted(set(modules.get("complete", [])) | set(args.complete))
    for item in args.reason:
        name, _, reason = item.partition("=")
        modules.setdefault("reasons", {})[name] = reason
    args.manifest.write_text(render(header, modules, entries))
    reasons = len(modules.get("reasons", {}))
    print(f"{len(entries)} entries, {len(modules['complete'])} complete modules, {reasons} with a reason")
    return 0


if __name__ == "__main__":
    sys.exit(main())
