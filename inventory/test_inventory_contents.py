# Copyright (c) 2026 Steven Kearnes. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Tests of the checked-in inventory."""

import json
import pathlib

import pytest

INVENTORY = pathlib.Path(__file__).parent / "inventory-2026.09.json"


@pytest.fixture(scope="module")
def inventory():
    return json.loads(INVENTORY.read_text())


@pytest.fixture(scope="module")
def by_name(inventory):
    return {entry["python"]: entry for entry in inventory["entries"]}


def test_core_entries(by_name, inventory):
    assert inventory["rdkit_version"] == "2026.09.1"
    assert by_name["rdkit.rdBase.rdkitVersion"]["kind"] == "module_attribute"
    assert by_name["rdkit.rdBase.DisableLog"]["kind"] == "module_function"
    assert by_name["rdkit.Chem.rdchem.Mol"]["kind"] == "class"
    assert by_name["rdkit.Chem.rdchem.Mol.GetAtomWithIdx"]["kind"] == "method"
    assert len(by_name["rdkit.Chem.rdmolfiles.MolFromSmiles"]["overloads"]) == 2
    assert by_name["rdkit.Chem.rdchem.BondType.SINGLE"]["kind"] == "enum_value"
    assert by_name["rdkit.Chem.rdmolfiles.SmilesParserParams.sanitize"]["kind"] == "property"


def test_sanitize_flags_are_an_int_enum(by_name):
    # nanobind binds SanitizeFlags with is_arithmetic, an IntEnum rather than
    # an enum.Flag, so its values are enum values (CONVENTIONS.md F1).
    assert by_name["rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL"]["kind"] == "enum_value"
    assert by_name["rdkit.Chem.rdmolops.SANITIZE_ALL"]["kind"] == "module_attribute"


def test_entries_are_sorted_and_unique(inventory):
    names = [entry["python"] for entry in inventory["entries"]]
    assert names == sorted(set(names))
