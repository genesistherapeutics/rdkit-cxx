# Copyright (c) 2026 Genesis Molecular AI
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Tests for check.py."""

import pytest

from check import check, merge, paths_source, snake_case, uncovered_modules, upper_camel_case


@pytest.mark.parametrize(
    "name,expected",
    [
        ("GetAtomWithIdx", "get_atom_with_idx"),
        ("allowCXSMILES", "allow_cxsmiles"),
        ("MolToSmiles", "mol_to_smiles"),
        ("rdkitVersion", "rdkit_version"),
        ("Is3D", "is3_d"),  # the case name_exception exists for
    ],
)
def test_snake_case(name, expected):
    assert snake_case(name) == expected


@pytest.mark.parametrize(
    "name,expected",
    [
        ("SINGLE", "Single"),
        ("CHI_TETRAHEDRAL_CW", "ChiTetrahedralCw"),
        ("SP3D2", "Sp3d2"),
    ],
)
def test_upper_camel_case(name, expected):
    assert upper_camel_case(name) == expected


def inventory(*entries):
    return {"rdkit_version": "2024.03.6", "entries": list(entries)}


def entry(python, kind, overloads=1):
    callable_kinds = {"module_function", "method", "static_method"}
    count = overloads if kind in callable_kinds else 0
    return {"python": python, "kind": kind, "overloads": [f"sig{i}" for i in range(count)]}


INVENTORY = inventory(
    entry("rdkit.rdBase.rdkitVersion", "module_attribute"),
    entry("rdkit.rdBase.DisableLog", "module_function"),
    entry("rdkit.Chem.rdchem.Mol", "class"),
    entry("rdkit.Chem.rdchem.Mol.__init__", "method", overloads=3),
    entry("rdkit.Chem.rdchem.Mol.__repr__", "method"),
    entry("rdkit.Chem.rdchem.Mol.GetConformer", "method"),
    entry("rdkit.Chem.rdchem.Mol.GetNumAtoms", "method"),
    entry("rdkit.Chem.rdchem.Mol.HasSubstructMatch", "method", overloads=2),
    entry("rdkit.Chem.rdchem.Conformer", "class"),
    entry("rdkit.Chem.rdchem.Conformer.Is3D", "method"),
    entry("rdkit.Chem.rdchem.BondType", "class"),
    entry("rdkit.Chem.rdchem.BondType.SINGLE", "enum_value"),
    entry("rdkit.Chem.rdmolops.SanitizeFlags", "class"),
    entry("rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL", "flag_value"),
    entry("rdkit.Geometry.rdGeometry.Point3D", "class"),
    entry("rdkit.Geometry.rdGeometry.Point3D.x", "property"),
)


def manifest(*entries, complete=()):
    return {"modules": {"complete": list(complete)}, "entry": list(entries)}


def problems(*entries, complete=()):
    return check(INVENTORY, manifest(*entries, complete=complete))


def test_a_name_missing_from_the_inventory():
    assert len(problems({"python": "rdkit.rdBase.Nope", "rust": ["rdkit_cxx::rd_base::nope"]})) == 1


def test_both_rust_and_reason():
    found = problems(
        {
            "python": "rdkit.rdBase.DisableLog",
            "rust": ["rdkit_cxx::rd_base::disable_log"],
            "reason": "no",
        }
    )
    assert len(found) == 1


def test_an_unmapped_entry_in_a_complete_module():
    found = problems(
        {"python": "rdkit.rdBase.DisableLog", "rust": ["rdkit_cxx::rd_base::disable_log"]},
        complete=["rdkit.rdBase"],
    )
    assert len(found) == 1
    assert "rdkit.rdBase.rdkitVersion" in found[0]


def test_a_wrong_rust_name():
    found = problems({"python": "rdkit.rdBase.DisableLog", "rust": ["rdkit_cxx::rd_base::disablelog"]})
    assert len(found) == 1


def test_two_overloads_and_one_rust_path():
    found = problems(
        {
            "python": "rdkit.Chem.rdchem.Mol.HasSubstructMatch",
            "rust": ["rdkit_cxx::chem::Mol::has_substruct_match"],
        }
    )
    assert len(found) == 1


@pytest.mark.parametrize(
    "item",
    [
        {"python": "rdkit.Chem.rdchem.Mol.GetConformer", "rust": ["rdkit_cxx::chem::Mol::get_conformer_with_id"]},
        {"python": "rdkit.Chem.rdchem.Mol.GetNumAtoms", "rust": ["rdkit_cxx::chem::Mol::get_num_atoms"]},
        {"python": "rdkit.rdBase.rdkitVersion", "rust": ["rdkit_cxx::rd_base::rdkit_version"]},
        {
            "python": "rdkit.Chem.rdchem.Conformer.Is3D",
            "rust": ["rdkit_cxx::chem::Conformer::is_3d"],
            "name_exception": "digits stay with their letters",
        },
        {
            "python": "rdkit.Chem.rdchem.Mol.HasSubstructMatch",
            "rust": ["rdkit_cxx::chem::Mol::has_substruct_match"],
            "overloads_note": "MolBundle is not bound yet",
        },
        {
            "python": "rdkit.Chem.rdchem.Mol.__init__",
            "rust": [
                "rdkit_cxx::chem::Mol::default",
                "rdkit_cxx::chem::Mol::clone",
                "rdkit_cxx::chem::Mol::from_mol_block",
            ],
        },
        {"python": "rdkit.Chem.rdchem.Mol.__init__", "rust": ["rdkit_cxx::chem::Mol::new"], "overloads_note": "x"},
        {
            "python": "rdkit.Chem.rdchem.Mol.__init__",
            "rust": ["rdkit_cxx::chem::Mol::new", "rdkit_cxx::chem::Mol::new_with_params"],
            "overloads_note": "x",
        },
        {"python": "rdkit.Chem.rdchem.Mol.__repr__", "reason": "Python protocol method"},
        {"python": "rdkit.Chem.rdchem.BondType.SINGLE", "rust": ["rdkit_cxx::chem::BondType::Single"]},
        {
            "python": "rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL",
            "rust": ["rdkit_cxx::chem::SanitizeFlags::SANITIZE_ALL"],
        },
        {"python": "rdkit.Chem.rdchem.Mol", "rust": ["rdkit_cxx::chem::Mol"]},
        {
            "python": "rdkit.Chem.rdchem.Conformer",
            "rust": ["rdkit_cxx::chem::ConformerRef", "rdkit_cxx::chem::ConformerMut"],
        },
        {
            "python": "rdkit.Geometry.rdGeometry.Point3D.x",
            "rust": ["rdkit_cxx::geometry::Point3D::x"],
            "rust_kind": "field",
        },
    ],
)
def test_no_problems(item):
    # rdBase is incomplete here, so its other unmapped names are fine.
    assert problems(item) == []


def test_classes_map_only_to_the_class_and_its_views():
    found = problems({"python": "rdkit.Chem.rdchem.Mol", "rust": ["rdkit_cxx::chem::MolHolder"]})
    assert len(found) == 1


def test_init_maps_only_to_constructor_names():
    found = problems(
        {"python": "rdkit.Chem.rdchem.Mol.__init__", "rust": ["rdkit_cxx::chem::Mol::make"], "overloads_note": "x"}
    )
    assert len(found) == 1


def test_generics_for_each_path():
    item = {
        "python": "rdkit.Chem.rdchem.Mol.__init__",
        "rust": ["rdkit_cxx::chem::Mol::default"],
        "rust_generics": ["", "<&str>"],
    }
    assert any("rust_generics lists 2" in p for p in check(INVENTORY, manifest(item)))


def test_a_duplicate_entry():
    item = {"python": "rdkit.rdBase.DisableLog", "rust": ["rdkit_cxx::rd_base::disable_log"]}
    assert len(problems(item, item)) == 1


def test_paths_source():
    source = paths_source(
        manifest(
            {"python": "rdkit.rdBase.DisableLog", "rust": ["rdkit_cxx::rd_base::disable_log"]},
            {"python": "rdkit.Chem.rdchem.Mol", "rust": ["rdkit_cxx::chem::Mol"]},
            {"python": "rdkit.Chem.rdchem.Mol.GetNumAtoms", "rust": ["rdkit_cxx::chem::Mol::get_num_atoms"]},
            {"python": "rdkit.Chem.rdchem.BondType.SINGLE", "rust": ["rdkit_cxx::chem::BondType::Single"]},
            {
                "python": "rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL",
                "rust": ["rdkit_cxx::chem::SanitizeFlags::SANITIZE_ALL"],
            },
            {
                "python": "rdkit.Geometry.rdGeometry.Point3D.x",
                "rust": ["rdkit_cxx::geometry::Point3D::x"],
                "rust_kind": "field",
            },
            {"python": "rdkit.Chem.rdchem.Mol.__repr__", "reason": "Python protocol method"},
            {
                "python": "rdkit.Chem.rdchem.Mol.GetConformer",
                "rust": ["rdkit_cxx::chem::Mol::from_file"],
                "rust_generics": "<&str>",
            },
            {
                "python": "rdkit.Chem.rdchem.Mol.__init__",
                "rust": ["rdkit_cxx::chem::Mol::default", "rdkit_cxx::chem::Mol::from_reader"],
                "rust_generics": ["", "<std::fs::File>"],
            },
        )
    )
    assert source.startswith("// Copyright (c) 2026 Genesis Molecular AI\n// All rights reserved.\n")
    assert source.count("#[test]") == 1
    assert "fn manifest_paths_exist()" in source
    # Module functions are imported, which also works for `impl Trait`
    # parameters, where `let _ = f;` cannot infer the type.
    assert "use rdkit_cxx::rd_base::disable_log as _;" in source
    assert "let _: Option<rdkit_cxx::chem::Mol> = None;" in source
    assert "let _ = <rdkit_cxx::chem::Mol>::get_num_atoms;" in source
    assert "let _ = <rdkit_cxx::chem::BondType>::Single;" in source
    assert "let _ = <rdkit_cxx::chem::SanitizeFlags>::SANITIZE_ALL;" in source
    assert "let _ = |value: &rdkit_cxx::geometry::Point3D| { let _ = &value.x; };" in source
    assert "let _ = <rdkit_cxx::chem::Mol>::from_file::<&str>;" in source
    # A list gives each path its own generics.
    assert "let _ = <rdkit_cxx::chem::Mol>::default;" in source
    assert "let _ = <rdkit_cxx::chem::Mol>::from_reader::<std::fs::File>;" in source
    assert "__repr__" not in source
    # A qualified path applies a generic type's default parameters.
    assert "<rdkit_cxx::chem::Mol>::get_num_atoms" in source


def test_a_module_reason_covers_its_names():
    found = problems(
        {"python": "rdkit.rdBase.DisableLog", "rust": ["rdkit_cxx::rd_base::disable_log"]},
        complete=["rdkit.rdBase"],
    )
    assert len(found) == 1
    found = check(
        INVENTORY,
        {"modules": {"complete": [], "reasons": {"rdkit.rdBase": "why"}}, "entry": []},
    )
    assert found == []


def test_a_module_reason_for_no_module():
    found = check(INVENTORY, {"modules": {"reasons": {"rdkit.Nope": "why"}}, "entry": []})
    assert len(found) == 1


def test_a_module_reason_and_complete():
    found = check(
        INVENTORY, {"modules": {"complete": ["rdkit.rdBase"], "reasons": {"rdkit.rdBase": "why"}}, "entry": []}
    )
    assert any("also marked complete" in problem for problem in found)


def test_uncovered_modules():
    manifest_ = {"modules": {"complete": ["rdkit.rdBase"], "reasons": {"rdkit.Chem": "why"}}, "entry": []}
    assert uncovered_modules(INVENTORY, manifest_) == ["rdkit.Geometry.rdGeometry"]


def test_merge():
    merged = merge(
        [
            inventory(entry("rdkit.rdBase.DisableLog", "module_function")),
            inventory(entry("rdkit.Chem.CanonSmiles", "module_function")),
        ]
    )
    assert [e["python"] for e in merged["entries"]] == ["rdkit.rdBase.DisableLog", "rdkit.Chem.CanonSmiles"]
    with pytest.raises(SystemExit):
        merge([inventory(), {"rdkit_version": "2027.03.1", "entries": []}])
