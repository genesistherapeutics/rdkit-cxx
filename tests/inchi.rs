// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/Chem/UnitTestInchi.py's tests of rdkit.Chem.inchi but those
// of MolBlockToInchi (not ported) and the PubChem sets (rdinchi's own
// tests); values from Python's module in RDKit 2026.09.1.

use rdkit_cxx::chem::{
    MolFromInchiParams, MolToInchiParams, inchi_available, inchi_to_inchi_key, mol_from_inchi,
    mol_from_inchi_and_aux_info, mol_from_inchi_and_aux_info_with_params, mol_from_mol_block,
    mol_from_smiles, mol_to_inchi, mol_to_inchi_and_aux_info, mol_to_inchi_key,
    mol_to_inchi_with_params,
};
use rdkit_cxx::{Error, Result};

const MB: &str = "
  Mrv1824 02111920092D

  6  6  0  0  0  0            999 V2000
   -5.5134    3.5259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    3.1134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5134    1.8759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    3.1134    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  2  0  0  0  0
M  END";

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
fn round_trip_atom_order() -> Result<()> {
    for smi in ["c1ccccc1O", "CC(=O)O", "C(=O)(N)C", "c1cc(O)ccc1N"] {
        let mol = mol_from_smiles(smi)?;
        let (inchi, aux) = mol_to_inchi_and_aux_info(&mol)?;
        let mol2 = mol_from_inchi_and_aux_info(&inchi, &aux)?.unwrap();
        let nums = |m: &rdkit_cxx::chem::Mol| {
            m.get_atoms()
                .map(|a| a.get_atomic_num())
                .collect::<Vec<_>>()
        };
        assert_eq!(nums(&mol), nums(&mol2), "{smi}");
    }
    // Stereo survives the round trip.
    let mol = mol_from_smiles("[C@@H](O)(F)Cl")?;
    let (inchi, aux) = mol_to_inchi_and_aux_info(&mol)?;
    let mol2 = mol_from_inchi_and_aux_info(&inchi, &aux)?.unwrap();
    assert_eq!(mol_to_inchi_and_aux_info(&mol2)?.0, inchi);
    // Several fragments.
    let mol = mol_from_smiles("CC.OO")?;
    let (inchi, aux) = mol_to_inchi_and_aux_info(&mol)?;
    assert_eq!(
        mol_from_inchi_and_aux_info(&inchi, &aux)?
            .unwrap()
            .get_num_atoms(),
        4
    );
    // Without removing hydrogens.
    let (inchi, aux) =
        mol_to_inchi_and_aux_info(&rdkit_cxx::chem::add_hs(&mol_from_smiles("OCC")?)?)?;
    assert_eq!(
        aux,
        "AuxInfo=1/0/N:3,2,1/rA:9OCCHHHHHH/rB:s1;s2;s1;s2;s2;s3;s3;s3;/rC:;;;;;;;;;"
    );
    let keep_hs = MolFromInchiParams {
        remove_hs: false,
        ..Default::default()
    };
    let mol = mol_from_inchi_and_aux_info_with_params(&inchi, &aux, &keep_hs)?.unwrap();
    let symbols: Vec<String> = mol.get_atoms().map(|a| a.get_symbol().unwrap()).collect();
    assert_eq!(symbols, ["O", "C", "C"]);
    Ok(())
}

#[test]
fn missing_aux_info() -> Result<()> {
    let mol = mol_from_inchi_and_aux_info("InChI=1S/CH4/h1H4", "")?.unwrap();
    assert_eq!(mol.get_num_atoms(), 1);
    assert!(mol_from_inchi_and_aux_info("not_an_inchi", "AuxInfo=1/0/N:1/")?.is_none());
    Ok(())
}

#[test]
fn coordinates() -> Result<()> {
    let mol = mol_from_mol_block(MB)?;
    let (inchi, aux) = mol_to_inchi_and_aux_info(&mol)?;
    let rc = aux.split("/rC:").nth(1).unwrap().split('/').next().unwrap();
    assert!(rc.ends_with(';'));
    let mol2 = mol_from_inchi_and_aux_info(&inchi, &aux)?.unwrap();
    assert_eq!(mol2.get_num_conformers(), 1);
    let (_, aux2) = mol_to_inchi_and_aux_info(&mol2)?;
    assert!(aux2.contains("/rC:") && !aux2.contains("/rC:;"));

    // Empty and all-zero coordinates give no conformer.
    let (inchi, aux) = mol_to_inchi_and_aux_info(&mol_from_smiles("C")?)?;
    let replace_rc = |aux: &str, with: &str| {
        let (head, tail) = aux.split_once("/rC:").unwrap();
        let rest = tail.find('/').map_or("", |i| &tail[i..]);
        format!("{head}/rC:{with}{rest}")
    };
    let mol2 = mol_from_inchi_and_aux_info(&inchi, &replace_rc(&aux, ";"))?.unwrap();
    assert_eq!(mol2.get_num_conformers(), 0);
    let (inchi, aux) = mol_to_inchi_and_aux_info(&mol_from_smiles("CCO")?)?;
    let mol2 = mol_from_inchi_and_aux_info(&inchi, &replace_rc(&aux, ",,;,,;,,;"))?.unwrap();
    assert_eq!(mol2.get_num_conformers(), 0);

    // Coordinates follow their atoms through the reordering.
    let mb = "
  Mrv1824 02111920092D

  3  2  0  0  0  0            999 V2000
    1.0000    2.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    4.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    6.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END";
    let mol = mol_from_mol_block(mb)?;
    let (inchi, aux) = mol_to_inchi_and_aux_info(&mol)?;
    let mol2 = mol_from_inchi_and_aux_info(&inchi, &aux)?.unwrap();
    assert_eq!(mol2.get_num_conformers(), 1);
    for (i, (sym, pos)) in [("O", [1.0, 2.0]), ("C", [3.0, 4.0]), ("N", [5.0, 6.0])]
        .iter()
        .enumerate()
    {
        assert_eq!(mol2.get_atom_with_idx(i as u32)?.get_symbol()?, *sym);
        let p = mol2.get_conformer()?.get_atom_position(i as u32)?;
        assert!((p.x - pos[0]).abs() < 1e-4 && (p.y - pos[1]).abs() < 1e-4 && p.z.abs() < 1e-4);
    }
    Ok(())
}
