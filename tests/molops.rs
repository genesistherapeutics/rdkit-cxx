// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn unsanitized(smiles: &str) -> Mol {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    mol_from_smiles_with_params(smiles, &params).unwrap()
}

#[test]
fn sanitize() {
    assert_eq!(
        sanitize_mol(&mut unsanitized("CCO")).unwrap(),
        SanitizeFlags::SANITIZE_NONE
    );
    assert!(matches!(
        sanitize_mol(&mut unsanitized("C(C)(C)(C)(C)C")),
        Err(Error::AtomValence { atom_idx: 0, .. })
    ));
    let catch = SanitizeMolParams {
        catch_errors: true,
        ..Default::default()
    };
    assert_eq!(
        sanitize_mol_with_params(&mut unsanitized("C(C)(C)(C)(C)C"), &catch).unwrap(),
        SanitizeFlags::SANITIZE_PROPERTIES
    );
    let skip = SanitizeMolParams {
        sanitize_ops: SanitizeFlags::SANITIZE_ALL & !SanitizeFlags::SANITIZE_PROPERTIES,
        ..Default::default()
    };
    assert_eq!(
        sanitize_mol_with_params(&mut unsanitized("C(C)(C)(C)(C)C"), &skip).unwrap(),
        SanitizeFlags::SANITIZE_NONE
    );
}

#[test]
fn sanitize_flags() {
    assert_eq!(SanitizeFlags::SANITIZE_ALL.bits(), 0xFFFFFFF);
    assert_eq!(
        (SanitizeFlags::SANITIZE_KEKULIZE | SanitizeFlags::SANITIZE_CLEANUP).bits(),
        0x9
    );
    let mut flags = SanitizeFlags::SANITIZE_NONE;
    flags |= SanitizeFlags::SANITIZE_ADJUSTHS;
    assert_eq!(flags, SanitizeFlags::from_bits_retain(0x200));
    flags &= SanitizeFlags::SANITIZE_CLEANUP;
    assert_eq!(flags, SanitizeFlags::SANITIZE_NONE);
    assert_eq!(
        format!("{:?}", SanitizeFlags::SANITIZE_ALL),
        "SanitizeFlags::SANITIZE_ALL"
    );
    assert_eq!(
        format!("{:?}", SanitizeFlags::from_bits_retain(0x9)),
        "SanitizeFlags(0x9)"
    );
    assert_eq!(
        SanitizeMolParams::default(),
        SanitizeMolParams {
            sanitize_ops: SanitizeFlags::SANITIZE_ALL,
            catch_errors: false,
        }
    );
}

#[test]
fn kekulize_forms() {
    let mut benzene = mol_from_smiles("c1ccccc1").unwrap();
    kekulize(&mut benzene).unwrap();
    assert!(benzene.get_atoms().all(|a| a.get_is_aromatic()));
    let clear = KekulizeParams {
        clear_aromatic_flags: true,
        ..Default::default()
    };
    kekulize_with_params(&mut benzene, &clear).unwrap();
    assert!(benzene.get_atoms().all(|a| !a.get_is_aromatic()));
    let kekule = SmilesWriteParams {
        do_kekule: true,
        ..Default::default()
    };
    assert_eq!(
        mol_to_smiles_with_params(&benzene, &kekule).unwrap(),
        "C1=CC=CC=C1"
    );
    assert!(matches!(kekulize(&mut unsanitized("c1cccc1")),
                     Err(Error::Kekulize { atom_indices, .. }) if atom_indices == [0, 1, 2, 3, 4]));
    assert_eq!(
        KekulizeParams::default(),
        KekulizeParams {
            clear_aromatic_flags: false,
            canonical: true
        }
    );
}

#[test]
fn add_hs_returns_a_new_molecule() {
    let methane = mol_from_smiles("C").unwrap();
    assert_eq!(add_hs(&methane).unwrap().get_num_atoms(), 5);
    assert_eq!(methane.get_num_atoms(), 1);
    let only_o = AddHsParameters {
        only_on_atoms: Some(vec![1]),
        ..Default::default()
    };
    assert_eq!(
        add_hs_with_params(&mol_from_smiles("CO").unwrap(), &only_o)
            .unwrap()
            .get_num_atoms(),
        3
    );
    let bad = AddHsParameters {
        only_on_atoms: Some(vec![7]),
        ..Default::default()
    };
    assert_eq!(
        add_hs_with_params(&methane, &bad).unwrap_err(),
        Error::Index { index: 7, len: 1 }
    );
    let defaults = AddHsParameters::default();
    assert!(
        !(defaults.explicit_only
            || defaults.add_coords
            || defaults.add_residue_info
            || defaults.skip_queries)
    );
    assert_eq!(defaults.only_on_atoms, None);
}

#[test]
fn fragments() {
    let mol = mol_from_smiles("CC.O").unwrap();
    assert_eq!(get_mol_frags(&mol), [vec![0, 1], vec![2]]);
    let frags: Vec<_> = get_mol_frags_as_mols(&mol)
        .unwrap()
        .iter()
        .map(|m| mol_to_smiles(m).unwrap())
        .collect();
    assert_eq!(frags, ["CC", "O"]);
    let raw = get_mol_frags_as_mols_with_sanitize_frags(&mol, false).unwrap();
    assert_eq!(raw.len(), 2);
    let empty = mol_from_smiles("").unwrap();
    assert!(get_mol_frags(&empty).is_empty());
    assert!(get_mol_frags_as_mols(&empty).unwrap().is_empty());
}

#[test]
fn removing_hs() {
    let mol = add_hs(&mol_from_smiles("[2H]C([H])O").unwrap()).unwrap();
    assert_eq!(mol.get_num_atoms(), 6);
    assert_eq!(remove_hs(&mol).unwrap().get_num_atoms(), 3);
    assert_eq!(remove_all_hs(&mol).unwrap().get_num_atoms(), 2);
    let isotopes = RemoveHsParameters {
        remove_isotopes: true,
        ..Default::default()
    };
    assert_eq!(
        remove_hs_with_params(&mol, &isotopes)
            .unwrap()
            .get_num_atoms(),
        2
    );
    // Python's implicitOnly=True.
    let implicit_only = RemoveHsParameters {
        remove_nonimplicit: false,
        ..Default::default()
    };
    assert_eq!(
        remove_hs_with_params(&mol, &implicit_only)
            .unwrap()
            .get_num_atoms(),
        4
    );
    let defaults = RemoveHsParameters::default();
    assert!(
        defaults.remove_in_s_groups && defaults.remove_mapped && defaults.remove_with_wedged_bond
    );
    assert!(defaults.remove_nonimplicit && defaults.show_warnings && defaults.sanitize);
    assert!(
        !defaults.remove_isotopes
            && !defaults.remove_degree_zero
            && !defaults.update_explicit_count
    );
}

#[test]
fn query_hs() {
    let query = mol_from_smarts("[#6][H]").unwrap();
    assert_eq!(has_query_hs(&query), (true, false));
    assert_eq!(
        has_query_hs(&mol_from_smarts("CC").unwrap()),
        (false, false)
    );
    assert_eq!(
        has_query_hs(&mol_from_smarts("[#6][#1,#6]").unwrap()),
        (true, true)
    );
    let merged = merge_query_hs(&query).unwrap();
    assert_eq!(mol_to_smarts(&merged).unwrap(), "[#6&!H0]");
    let mapped = mol_from_smarts("[#6][H:1]").unwrap();
    let unmapped_only = MergeQueryHsParams {
        merge_unmapped_only: true,
        ..Default::default()
    };
    assert_eq!(
        merge_query_hs_with_params(&mapped, &unmapped_only)
            .unwrap()
            .get_num_atoms(),
        2
    );
    assert_eq!(merge_query_hs(&mapped).unwrap().get_num_atoms(), 1);
}
