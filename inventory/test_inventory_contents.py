# Copyright (c) 2026 Genesis Molecular AI
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Tests of the checked-in inventory."""

import json
import pathlib

import pytest

INVENTORY = pathlib.Path(__file__).parent / "inventory-2025.03.json"


@pytest.fixture(scope="module")
def inventory():
    return json.loads(INVENTORY.read_text())


@pytest.fixture(scope="module")
def by_name(inventory):
    return {entry["python"]: entry for entry in inventory["entries"]}


def test_core_entries(by_name, inventory):
    assert inventory["rdkit_version"] == "2025.03.6"
    assert by_name["rdkit.rdBase.rdkitVersion"]["kind"] == "module_attribute"
    assert by_name["rdkit.rdBase.DisableLog"]["kind"] == "module_function"
    assert by_name["rdkit.Chem.rdchem.Mol"]["kind"] == "class"
    assert by_name["rdkit.Chem.rdchem.Mol.GetAtomWithIdx"]["kind"] == "method"
    assert len(by_name["rdkit.Chem.rdmolfiles.MolFromSmiles"]["overloads"]) == 2
    assert by_name["rdkit.Chem.rdchem.BondType.SINGLE"]["kind"] == "enum_value"
    assert by_name["rdkit.Chem.rdmolfiles.SmilesParserParams.sanitize"]["kind"] == "property"


def test_sanitize_flags_are_an_int_enum(by_name):
    # Boost.Python binds SanitizeFlags as a plain enum, so its values are enum
    # values (CONVENTIONS.md F1).
    assert by_name["rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL"]["kind"] == "enum_value"
    assert by_name["rdkit.Chem.rdmolops.SANITIZE_ALL"]["kind"] == "module_attribute"


def test_boost_python_entries(by_name):
    # Methods and class attributes RDKit's Python adds to compiled classes.
    assert by_name["rdkit.Chem.rdchem.Mol.GetAtoms"]["kind"] == "method"
    assert by_name["rdkit.DataStructs.cDataStructs.ExplicitBitVect.ToBitString"]["kind"] == "method"
    assert by_name["rdkit.Chem.rdMolDescriptors.AtomPairsParameters.codeSize"]["kind"] == "property"
    # Boost.Python's generic __reduce__ is left out; pickling hooks are not.
    assert "rdkit.Chem.rdchem.Atom.__reduce__" not in by_name
    assert by_name["rdkit.Chem.rdchem.Mol.__getinitargs__"]["kind"] == "method"
    # Overloads come from the docstring, with their defaults.
    assert by_name["rdkit.Geometry.rdGeometry.UniformGrid3D"]["overloads"][0].startswith(
        "UniformGrid3D( (float)dimX, (float)dimY, (float)dimZ [, (float)spacing=0.5"
    )


def test_entries_are_sorted_and_unique(inventory):
    names = [entry["python"] for entry in inventory["entries"]]
    assert names == sorted(set(names))


PYTHON_INVENTORY = pathlib.Path(__file__).parent / "python-2025.03.json"


@pytest.fixture(scope="module")
def python_inventory():
    return json.loads(PYTHON_INVENTORY.read_text())


@pytest.fixture(scope="module")
def python_by_name(python_inventory):
    return {entry["python"]: entry for entry in python_inventory["entries"]}


def test_python_entries(python_by_name, python_inventory):
    assert python_inventory["rdkit_version"] == "2025.03.6"
    assert python_by_name["rdkit.Chem.FindMolChiralCenters"]["kind"] == "module_function"
    assert python_by_name["rdkit.Chem.Lipinski.HDonorSmarts"]["kind"] == "module_attribute"
    assert python_by_name["rdkit.Chem.RegistrationHash.HashLayer.CANONICAL_SMILES"]["kind"] == "enum_value"
    assert python_by_name["rdkit.ML.Cluster.Butina.ClusterData"]["kind"] == "module_function"


def test_python_entries_made_at_import(python_by_name):
    # Lambdas and functions modules make as they import are named for the
    # attribute that holds them; descList is assigned through `global`.
    assert python_by_name["rdkit.Chem.Lipinski.NumHDonors"]["overloads"] == ["def NumHDonors(x)"]
    assert python_by_name["rdkit.Chem.Fragments.fr_Al_COO"]["kind"] == "module_function"
    assert python_by_name["rdkit.Chem.Descriptors.descList"]["kind"] == "module_attribute"


def test_python_entries_leave_out_imports_tests_and_loggers(python_by_name):
    assert "rdkit.Chem.AllChem.EmbedMolecule" not in python_by_name  # from rdDistGeom
    assert "rdkit.Chem.AllChem.logger" not in python_by_name
    assert not any(".UnitTest" in name or ".conftest" in name for name in python_by_name)


def test_python_signatures_hold_no_addresses(python_inventory):
    assert not any(" at 0x" in overload for entry in python_inventory["entries"] for overload in entry["overloads"])


def test_python_entries_are_sorted_and_unique(python_inventory):
    names = [entry["python"] for entry in python_inventory["entries"]]
    assert names == sorted(set(names))
