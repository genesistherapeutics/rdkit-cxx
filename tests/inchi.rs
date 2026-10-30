// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/Chem/UnitTestInchi.py's tests of rdkit.Chem.inchi but those
// of MolBlockToInchi (not ported) and the PubChem sets (rdinchi's own
// tests); values from Python's module in RDKit 2025.09.6.

use rdkit_cxx::chem::{
    MolToInchiParams, inchi_available, inchi_to_inchi_key, mol_from_inchi, mol_from_smiles,
    mol_to_inchi, mol_to_inchi_and_aux_info, mol_to_inchi_key, mol_to_inchi_with_params,
};
use rdkit_cxx::{Error, Result};

#[test]
fn perchloric_acid() -> Result<()> {
    assert!(inchi_available());
    for (smiles, expected) in [
        ("OCl(=O)(=O)=O", "InChI=1S/ClHO4/c2-1(3,4)5/h(H,2,3,4,5)"),
        (
            "CC1=CC2=NCC(CN2C=C1)C(=O)c3ccc4cc(C)ccc4c3.OCl(=O)(=O)=O",
            "InChI=1S/C21H20N2O.ClHO4/c1-14-3-4-17-11-18(6-5-16(17)9-14)21(24)19-12-22-20-10-15(2)7-8-23(20)13-19;2-1(3,4)5/h3-11,19H,12-13H2,1-2H3;(H,2,3,4,5)",
        ),
        (
            "CNc1ccc2nc3ccccc3[n+](C)c2c1.[O-]Cl(=O)(=O)=O",
            "InChI=1S/C14H13N3.ClHO4/c1-15-10-7-8-12-14(9-10)17(2)13-6-4-3-5-11(13)16-12;2-1(3,4)5/h3-9H,1-2H3;(H,2,3,4,5)",
        ),
    ] {
        assert_eq!(mol_to_inchi(&mol_from_smiles(smiles)?)?, expected);
    }
    Ok(())
}

#[test]
fn options_and_keys() -> Result<()> {
    let m = mol_from_smiles("CC=C(N)C")?;
    let inchi1 = mol_to_inchi(&m)?;
    let suu = MolToInchiParams {
        options: "/SUU".into(),
        ..Default::default()
    };
    let inchi2 = mol_to_inchi_with_params(&m, &suu)?;
    assert_eq!(
        format!("{}/b4-3?", inchi1.split_once('/').unwrap().1),
        inchi2.split_once('/').unwrap().1
    );
    assert_eq!(
        inchi_to_inchi_key("InChI=1S/C9H12/c1-2-6-9-7-4-3-5-8-9/h3-5,7-8H,2,6H2,1H3").as_deref(),
        Some("ODLMAHJVESYWTB-UHFFFAOYSA-N")
    );
    assert_eq!(inchi_to_inchi_key(&inchi1), Some(mol_to_inchi_key(&m)?));
    assert_eq!(mol_to_inchi_key(&m)?, "KZFUWRMVPLBQKD-UHFFFAOYSA-N");
    assert_eq!(inchi_to_inchi_key("xx"), None);
    Ok(())
}

#[test]
fn warnings_as_errors() -> Result<()> {
    let m = mol_from_smiles("[CH3]")?;
    assert_eq!(mol_to_inchi(&m)?, "InChI=1S/CH3/h1H3");
    let strict = MolToInchiParams {
        treat_warning_as_error: true,
        ..Default::default()
    };
    match mol_to_inchi_with_params(&m, &strict) {
        Err(Error::Value { message }) => assert_eq!(message, "Accepted unusual valence(s): C(3)"),
        other => panic!("{other:?}"),
    }
    assert!(mol_from_inchi("InChI=1S/CH5/h1H5")?.is_none());
    Ok(())
}

#[test]
fn inchi_and_aux_info() -> Result<()> {
    // RDKit 2025.09 has no MolFromInchiAndAuxInfo, which reads the AuxInfo
    // back; the AuxInfo numbers the input's atoms.
    let (inchi, aux) = mol_to_inchi_and_aux_info(&mol_from_smiles("OCC")?)?;
    assert_eq!(inchi, "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
    assert!(aux.starts_with("AuxInfo=1/0/N:3,2,1/"), "{aux}");
    assert_eq!(mol_to_inchi(&mol_from_inchi(&inchi)?.unwrap())?, inchi);
    Ok(())
}
