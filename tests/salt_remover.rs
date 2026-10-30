// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/Chem/UnitTestSaltRemover.py and the module's doctests; values
// from Python's SaltRemover in RDKit 2024.09.6.

use rdkit_cxx::chem::salt_remover::*;
use rdkit_cxx::chem::{
    Mol, SmilesParserParams, mol_from_smarts, mol_from_smiles, mol_from_smiles_with_params,
    mol_to_smarts, mol_to_smiles,
};
use rdkit_cxx::{Error, Result};

fn data(text: &str, format: InputFormat) -> Result<SaltRemover> {
    SaltRemover::new_with_params(&SaltRemoverParams {
        defn_data: Some(text.into()),
        defn_format: format,
        ..Default::default()
    })
}

fn smarts(mols: &[Mol]) -> Vec<String> {
    mols.iter().map(|m| mol_to_smarts(m).unwrap()).collect()
}

#[test]
fn definitions() -> Result<()> {
    let remover = SaltRemover::new();
    assert_eq!(remover.salts().len(), 15);
    assert_eq!(
        smarts(&remover.salts()[..5]),
        [
            "[Cl,Br,I]",
            "[#3,#11,#19,#20,#12]",
            "[O,N]",
            "N(=O)(O)O",
            "P(=O)(O)(O)O"
        ]
    );
    assert_eq!(remover.defn_filename(), None);
    assert_eq!(data("[Cl,Br]", InputFormat::Smarts)?.salts().len(), 1);
    assert_eq!(
        data("[Na+]\nCC(=O)O", InputFormat::Smiles)?.salts().len(),
        2
    );
    match data("[Cl,fail]", InputFormat::Smarts) {
        Err(Error::Value { message }) => assert_eq!(message, "[Cl,fail]"),
        other => panic!("{other:?}"),
    }
    // A blank line in data is an error, as in Python.
    match data("[Cl,Br]\n\n[Na]", InputFormat::Smarts) {
        Err(Error::Value { message }) => assert_eq!(message, "\n"),
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        data("C", InputFormat::Mol),
        Err(Error::Value { .. })
    ));
    Ok(())
}

#[test]
fn definition_files() -> Result<()> {
    let sdf = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/witch-salts.sdf");
    let remover = SaltRemover::new_with_params(&SaltRemoverParams {
        defn_filename: Some(sdf.into()),
        defn_format: InputFormat::Mol,
        ..Default::default()
    })?;
    assert_eq!(remover.salts().len(), 240);
    assert_eq!(remover.defn_filename(), Some(std::path::Path::new(sdf)));
    let m = mol_from_smiles("Cc1onc(-c2ccccc2)c1C([O-])=NC1C(=O)N2C1SC(C)(C)C2C(=O)O.O.[Na+]")?;
    let stripped = remover.strip_mol_with_deleted(&m)?;
    assert_eq!(
        mol_to_smiles(&stripped.mol)?,
        "Cc1onc(-c2ccccc2)c1C([O-])=NC1C(=O)N2C1SC(C)(C)C2C(=O)O.O"
    );
    assert_eq!(stripped.deleted.len(), 1);
    assert_eq!(mol_to_smiles(&stripped.deleted[0])?, "[Na+]");

    // test_withDontRemoveEverything
    let m = mol_from_smiles("Cc1ccccc1")?;
    let stripped = remover.strip_mol_with_deleted_with_dont_remove_everything(&m, true)?;
    assert!(stripped.deleted.is_empty());
    assert_eq!(mol_to_smiles(&stripped.mol)?, "Cc1ccccc1");

    let smi = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/c6h6-cdk.smi");
    let remover = SaltRemover::new_with_params(&SaltRemoverParams {
        defn_filename: Some(smi.into()),
        defn_format: InputFormat::Smiles,
        ..Default::default()
    })?;
    assert_eq!(remover.salts().len(), 216);
    Ok(())
}

#[test]
fn stripping() -> Result<()> {
    let remover = data("[Cl,Br]", InputFormat::Smarts)?;
    for (smiles, atoms) in [
        ("CN(C)C.Cl", 4),
        ("CN(C)C.Cl.Cl.Br", 4),
        ("CN(Br)Cl", 4),
        ("CN(Br)Cl.Cl", 4),
        ("C[NH+](C)(C).[Cl-]", 4),
    ] {
        let res = remover.strip_mol(&mol_from_smiles(smiles)?)?;
        assert_eq!(res.get_num_atoms(), atoms, "{smiles}");
    }
    let res = remover.strip_mol(&mol_from_smiles("CN(Br)Cl.Cl")?)?;
    assert_eq!(mol_to_smiles(&res)?, "CN(Cl)Br");
    let stripped = remover.strip_mol_with_deleted(&mol_from_smiles("CN(C)C.Cl.Br")?)?;
    assert_eq!(mol_to_smiles(&stripped.mol)?, "CN(C)C");
    assert_eq!(smarts(&stripped.deleted), ["[Cl,Br]"]);
    let stripped = remover.strip_mol_with_deleted(&mol_from_smiles("CN(C)C.Cl.Cl")?)?;
    assert_eq!(stripped.mol.get_num_atoms(), 4);
    assert_eq!(smarts(&stripped.deleted), ["[Cl,Br]"]);

    let remover = data("Cl", InputFormat::Smiles)?;
    let res = remover.strip_mol(&mol_from_smiles("CN(Br)Cl.Cl")?)?;
    assert_eq!(mol_to_smiles(&res)?, "CN(Cl)Br");

    let remover = data("[Na+]\nCC(=O)O", InputFormat::Smiles)?;
    assert_eq!(
        remover
            .strip_mol(&mol_from_smiles("CC(=O)O.[Na+]")?)?
            .get_num_atoms(),
        0
    );

    // RDKit's own salts.
    let remover = SaltRemover::default();
    let dont = StripMolParams {
        dont_remove_everything: true,
        ..Default::default()
    };
    for (smiles, stripped, kept, deleted) in [
        ("CN(C)C.Cl", "CN(C)C", "CN(C)C", &["[Cl,Br,I]"][..]),
        (
            "CC(=O)O.[Na]",
            "",
            "CC(=O)O",
            &["[#3,#11,#19,#20,#12]", "[C&H3]C(=O)O"],
        ),
        ("Cl.Cl", "", "Cl.Cl", &["[Cl,Br,I]"]),
        ("CCO", "CCO", "CCO", &[]),
        (
            "C[NH+](C)(C).[Cl-]",
            "C[NH+](C)C",
            "C[NH+](C)C",
            &["[Cl,Br,I]"],
        ),
        ("OC(=O)C(=O)O.CCN", "CCN", "CCN", &["OC(=O)C(=O)O"]),
    ] {
        let m = mol_from_smiles(smiles)?;
        assert_eq!(
            mol_to_smiles(&remover.strip_mol(&m)?)?,
            stripped,
            "{smiles}"
        );
        assert_eq!(
            mol_to_smiles(&remover.strip_mol_with_params(&m, &dont)?)?,
            kept,
            "{smiles}"
        );
        assert_eq!(
            smarts(&remover.strip_mol_with_deleted(&m)?.deleted),
            deleted,
            "{smiles}"
        );
    }
    Ok(())
}

#[test]
fn chirality() -> Result<()> {
    let defn = "OC(=O)\\C=C/C(O)=O\nOC(=O)/C=C/C(O)=O";
    let chiral = SaltRemover::new_with_params(&SaltRemoverParams {
        defn_data: Some(defn.into()),
        defn_format: InputFormat::Smiles,
        use_chirality: true,
        ..Default::default()
    })?;
    let remover = data(defn, InputFormat::Smiles)?;
    let maleic = mol_to_smiles(&mol_from_smiles("OC(=O)\\C=C/C(O)=O")?)?;
    let fumaric = mol_to_smiles(&mol_from_smiles("OC(=O)/C=C/C(O)=O")?)?;
    let deleted = |r: &SaltRemover, smiles: &str| -> Result<Vec<String>> {
        r.strip_mol_with_deleted(&mol_from_smiles(smiles)?)?
            .deleted
            .iter()
            .map(mol_to_smiles)
            .collect()
    };
    assert_eq!(deleted(&remover, "OC(=O)C=CC(O)=O")?, [maleic.as_str()]);
    assert!(deleted(&chiral, "OC(=O)C=CC(O)=O")?.is_empty());
    assert_eq!(deleted(&remover, &maleic)?, [maleic.as_str()]);
    assert_eq!(deleted(&chiral, &maleic)?, [maleic.as_str()]);
    assert_eq!(deleted(&remover, &fumaric)?, [maleic.as_str()]);
    assert_eq!(deleted(&chiral, &fumaric)?, [fumaric.as_str()]);
    Ok(())
}

#[test]
fn github_issues() -> Result<()> {
    // GitHub 4550
    let m = mol_from_smiles_with_params(
        "Cl.C[N]1=CC=CC=C1",
        &SmilesParserParams {
            sanitize: false,
            ..Default::default()
        },
    )?;
    assert_eq!(m.get_num_atoms(), 8);
    let params = StripMolParams {
        sanitize: false,
        ..Default::default()
    };
    let res = SaltRemover::new().strip_mol_with_params(&m, &params)?;
    assert_eq!(mol_to_smiles(&res)?, "CN1=CC=CC=C1");
    // GitHub 7327
    let m = SaltRemover::new().strip_mol(&mol_from_smiles("C=CC=O")?)?;
    assert_eq!(m.get_num_atoms(), 4);
    let q = mol_from_smarts("[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]")?;
    assert_eq!(m.get_substruct_matches(&q)?, [[1, 2]]);
    Ok(())
}
