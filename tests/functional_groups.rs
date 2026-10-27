// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/Chem/UnitTestFunctionalGroups.py; values from Python's
// FunctionalGroups in RDKit 2025.03.6.

use rdkit_cxx::chem::functional_groups::*;
use rdkit_cxx::chem::{add_hs, mol_from_smiles};
use rdkit_cxx::{Error, Result};

const TXT: &str = "
AcidChloride\tC(=O)Cl\tAcid Chloride
  AcidChloride.Benzoyl\tC(=O)(Cl)c1ccccc1\tBenzoyl
Amine\tN\tAmine
  Amine.Primary\t[N;H2]\tPrimary
    Amine.Primary.Aromatic\t[N;H2][a]\tPrimary Aromatic
  Amine.Aromatic\tN[a]\tAromatic
";

fn from_data(data: &str) -> Result<Vec<FGHierarchyNode>> {
    build_func_group_hierarchy_with_params(&BuildFuncGroupHierarchyParams {
        data: Some(data.into()),
        ..Default::default()
    })
}

#[test]
fn basics() -> Result<()> {
    let hierarchy = from_data(TXT)?;
    assert_eq!(hierarchy.len(), 2);
    assert_eq!((hierarchy[0].len(), hierarchy[1].len()), (2, 4));
    assert_eq!(hierarchy[0].name, "Acid Chloride");
    assert_eq!(hierarchy[0].children[0].name, "Benzoyl");
    assert_eq!(hierarchy[0].label, "AcidChloride");
    assert_eq!(hierarchy[0].rxn_smarts, "");
    assert_eq!(hierarchy[0].parent, None);
    assert_eq!(
        hierarchy[0].children[0].parent.as_deref(),
        Some("AcidChloride")
    );
    for (smiles, fp) in [
        ("ClC(=O)CCCNc1ccccc1", [1, 0, 1, 0, 0, 1]),
        ("OC(=O)CCC", [0, 0, 0, 0, 0, 0]),
        ("Nc1ccccc1N", [0, 0, 2, 2, 2, 2]),
        ("NCCNC(N)Cl", [0, 0, 3, 2, 0, 0]),
    ] {
        assert_eq!(
            create_mol_fingerprint(&mol_from_smiles(smiles)?, &hierarchy)?,
            fp,
            "{smiles}"
        );
    }
    Ok(())
}

#[test]
fn comments() -> Result<()> {
    let txt = "
AcidChloride\tC(=O)Cl\tAcid Chloride
  AcidChloride.Benzoyl\tC(=O)(Cl)c1ccccc1\tBenzoyl
Amine\tN\tAmine
  Amine.Primary\t[N;H2]\tPrimary
    //Amine.Primary.Aromatic\t[N;H2][a]\tPrimary Aromatic
  Amine.Aromatic\tN[a]\tAromatic
";
    let hierarchy = from_data(txt)?;
    assert_eq!(hierarchy.len(), 2);
    assert_eq!((hierarchy[0].len(), hierarchy[1].len()), (2, 3));
    Ok(())
}

#[test]
fn reactions() -> Result<()> {
    let txt = "BoronicAcid\t[$(B-!@[#6])](O)(O)\tBoronic Acid\t[#6:1][B:2]([O:3])[O:4]>>[#6:1].[B:2]([O:3])[O:4]
  BoronicAcid.Aromatic\t[$(B-!@c)](O)(O)\tAromatic\t[c:1][B:2]([O:3])[O:4]>>[c:1].[B:2]([O:3])[O:4]
  BoronicAcid.Aliphatic\t[$(B-!@C)](O)(O)\tAliphatic\t[C:1][B:2]([O:3])[O:4]>>[C:1].[B:2]([O:3])[O:4]
  ";
    let hierarchy = from_data(txt)?;
    assert_eq!(hierarchy.len(), 1);
    assert_eq!(hierarchy[0].children.len(), 2);
    assert_ne!(hierarchy[0].rxn_smarts, "");
    assert_ne!(hierarchy[0].children[0].rxn_smarts, "");
    // Empty fields after the name, and comments after tabs.
    let hierarchy = from_data("A\tC\tX\t[C:1]>>[N:1]\t\nA.B\tCC\tY\t\t// c\n")?;
    assert_eq!(hierarchy[0].rxn_smarts, "[C:1]>>[N:1]");
    assert_eq!(hierarchy[0].children[0].rxn_smarts, "");
    Ok(())
}

#[test]
fn rdkit_hierarchy() -> Result<()> {
    let hierarchy = build_func_group_hierarchy()?;
    assert_eq!(hierarchy.len(), 12);
    assert_eq!(
        hierarchy.iter().map(FGHierarchyNode::len).sum::<usize>(),
        52
    );
    let labels: Vec<&str> = hierarchy[..6].iter().map(|n| n.label.as_str()).collect();
    assert_eq!(
        labels,
        [
            "AcidChloride",
            "CarboxylicAcid",
            "SulfonylChloride",
            "Amine",
            "BoronicAcid",
            "Isocyanate"
        ]
    );
    assert_eq!(hierarchy[0].children[0].label, "AcidChloride.Aromatic");
    let boronic = hierarchy.iter().find(|n| n.label == "BoronicAcid").unwrap();
    assert_eq!(
        boronic.rxn_smarts,
        "[#6:1]-!@[B:2]([O:3])[O:4]>>[#6:1][X].[B:2]([O:3])[O:4]"
    );
    assert_eq!(boronic.smarts, "[$(B-!@[#6])](O)(O)");
    for (smiles, on, counts) in [
        ("CC(=O)Cl", &[0, 2][..], &[1, 1][..]),
        ("c1ccccc1C(=O)O", &[3, 4], &[1, 1]),
        (
            "NCC(=O)O",
            &[3, 5, 6, 10, 11, 13, 21],
            &[1, 1, 1, 1, 1, 1, 1],
        ),
        ("OB(O)c1ccccc1", &[23, 24], &[1, 1]),
        ("CC(C)(C)OC(=O)N1CCCC1", &[], &[]),
    ] {
        let fp = create_mol_fingerprint(&mol_from_smiles(smiles)?, &hierarchy)?;
        let bits: Vec<usize> = (0..fp.len()).filter(|&i| fp[i] > 0).collect();
        assert_eq!(bits, on, "{smiles}");
        let nonzero: Vec<u32> = fp.iter().copied().filter(|&x| x > 0).collect();
        assert_eq!(nonzero, counts, "{smiles}");
    }
    Ok(())
}

#[test]
fn fingerprints_with_hydrogens() -> Result<()> {
    // test4Hs asserts that hydrogens never change a fingerprint, but in
    // RDKit 2025.03.6 they change 725 of these 4978 in Python, the first
    // for this molecule; the port computes what Python computes.
    let hierarchy = build_func_group_hierarchy()?;
    let csv = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/NCI_5K_TPSA.csv"
    ))
    .unwrap();
    let mut differ = Vec::new();
    let mut n = 0;
    for line in csv.lines().filter(|l| !l.starts_with('#')) {
        let smiles = line.split(',').next().unwrap();
        let m = mol_from_smiles(smiles)?;
        n += 1;
        if create_mol_fingerprint(&m, &hierarchy)?
            != create_mol_fingerprint(&add_hs(&m)?, &hierarchy)?
        {
            differ.push(smiles);
        }
    }
    assert_eq!((n, differ.len()), (4978, 725));
    assert_eq!(differ[0], "NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O");
    Ok(())
}

#[test]
fn malformed_files() {
    for (data, message) in [
        ("A\tC\n", r"Input line 1 ('A\tC') is not long enough."),
        ("A\tC\tX\nA\tN\tY\n", "Duplicate label on line 2."),
        ("A.B\tC\tX\n", "Hierarchy member A (line 1) not found."),
        ("A\tC(\tX\n", r#"Smarts "C(" (line 1) could not be parsed."#),
        (
            "\n\nA\tC\tX // ok\n B\tO\n",
            r"Input line 4 ('B\tO') is not long enough.",
        ),
    ] {
        match from_data(data) {
            Err(Error::Value { message: m }) => assert_eq!(m, message),
            other => panic!("{other:?}"),
        }
    }
    let missing = BuildFuncGroupHierarchyParams {
        file_nm: Some("/nonexistent/groups.txt".into()),
        ..Default::default()
    };
    assert!(matches!(
        build_func_group_hierarchy_with_params(&missing),
        Err(Error::Os { .. })
    ));
}

#[test]
fn nodes() -> Result<()> {
    let node = FGHierarchyNode::new("Amine", mol_from_smiles("N")?);
    assert_eq!(
        (node.len(), node.label.as_str(), node.parent.as_deref()),
        (1, "", None)
    );
    let params = FGHierarchyNodeParams {
        smarts: "[NH2]".into(),
        label: "Amine.Primary".into(),
        parent: Some("Amine".into()),
        ..Default::default()
    };
    let node = FGHierarchyNode::new_with_params("Primary", mol_from_smiles("N")?, &params);
    assert_eq!(node.parent.as_deref(), Some("Amine"));
    assert!(node.removal_reaction.is_none());
    Ok(())
}
