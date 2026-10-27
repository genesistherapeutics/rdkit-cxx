// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// RDKit 2024.09's bond removals leave the stereo groups that hold the bond
// pointing at it, deleted (rdkit/rdkit#8558). RWMol's own removals take the
// bond out of the groups, as RDKit 2025.09 does, and RDKit's other removals
// refuse a molecule whose groups they would leave so.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

// The molecule with an AND group of `atoms` and `bonds`.
fn with_group(mol: Mol, atoms: &[u32], bonds: &[u32]) -> RWMol {
    let mut rw = RWMol::from(mol);
    let params = CreateStereoGroupParams {
        atom_ids: atoms.to_vec(),
        bond_ids: bonds.to_vec(),
        ..Default::default()
    };
    let group = create_stereo_group_with_params(StereoGroupType::StereoAnd, &rw, &params).unwrap();
    rw.set_stereo_groups(&[group]).unwrap();
    rw
}

// Each group's atoms and bonds.
fn groups(mol: &Mol) -> Vec<(Vec<u32>, Vec<u32>)> {
    mol.get_stereo_groups()
        .iter()
        .map(|g| (g.get_atoms().to_vec(), g.get_bonds().to_vec()))
        .collect()
}

// What reads the groups' bonds: a copy, a pickle, and the writers.
fn assert_groups_read(mol: &Mol) {
    let copy = mol.clone();
    assert_eq!(groups(&copy), groups(mol));
    // SAFETY: a pickle RDKit wrote.
    let read = unsafe { Mol::from_pkl_string(&mol.to_binary()) }.unwrap();
    assert_eq!(groups(&read), groups(mol));
    mol_to_cx_smiles(mol).unwrap();
    mol_to_v3k_mol_block(mol).unwrap();
}

fn assert_value_error<T: std::fmt::Debug>(result: Result<T, Error>) {
    match result {
        Err(Error::Value { message }) => assert!(message.contains("stereo group"), "{message}"),
        other => panic!("expected Error::Value, got {other:?}"),
    }
}

#[test]
fn removing_a_bond_takes_it_out_of_groups() {
    let mut rw = with_group(smiles("CCCC"), &[0], &[1, 2]);
    rw.remove_bond(2, 3).unwrap();
    assert_eq!(groups(&rw), [(vec![0], vec![1])]);
    assert_groups_read(&rw);

    // A group of only the bond goes with it.
    let mut rw = with_group(smiles("CCCC"), &[], &[1]);
    rw.remove_bond(1, 2).unwrap();
    assert!(groups(&rw).is_empty());
    assert_groups_read(&rw);
}

#[test]
fn removing_an_atom_renumbers_groups() {
    let mut rw = with_group(smiles("CCCCC"), &[4], &[3]);
    rw.remove_atom(0).unwrap();
    assert_eq!(groups(&rw), [(vec![3], vec![2])]);
    assert_groups_read(&rw);

    let mut rw = with_group(smiles("CCCCC"), &[0], &[3]);
    rw.remove_atom(4).unwrap();
    assert_eq!(groups(&rw), [(vec![0], vec![])]);
    assert_groups_read(&rw);
}

#[test]
fn batch_removals_update_groups_at_the_commit() {
    let mut rw = with_group(smiles("CCCCC"), &[0, 4], &[1, 3]);
    rw.begin_batch_edit().unwrap();
    rw.remove_bond(1, 2).unwrap();
    rw.remove_atom(0).unwrap();
    assert_eq!(groups(&rw), [(vec![0, 4], vec![1, 3])]);
    rw.commit_batch_edit().unwrap();
    assert_eq!(groups(&rw), [(vec![3], vec![1])]);
    assert_groups_read(&rw);

    let mut rw = with_group(smiles("CCCC"), &[0], &[1]);
    rw.begin_batch_edit().unwrap();
    rw.remove_bond(1, 2).unwrap();
    rw.rollback_batch_edit();
    assert_eq!(groups(&rw), [(vec![0], vec![1])]);
}

#[test]
fn splitting_fragments_refuses_groups_across_them() {
    // The group's bond is in the first fragment and its atom in the second,
    // whose copy would keep the group and lose the bond.
    let across = with_group(smiles("CC.CC"), &[2], &[0]);
    assert_value_error(mol_to_smiles(&across));
    assert_value_error(get_mol_frags_as_mols(&across));

    for mol in [
        with_group(smiles("CC.CC"), &[0], &[0]),
        with_group(smiles("CC.CC"), &[], &[0, 1]),
        with_group(smiles("CCCC"), &[0], &[1]),
    ] {
        mol_to_smiles(&mol).unwrap();
        for frag in get_mol_frags_as_mols(&mol).unwrap() {
            assert_groups_read(&frag);
        }
    }
}

#[test]
fn removing_hydrogens_refuses_groups_of_their_bonds() {
    // Bond 0 joins the carbons; bond 1 is a hydrogen's.
    let hs = add_hs(&smiles("CC")).unwrap();
    assert_value_error(remove_hs(&with_group(hs.clone(), &[0], &[1])));
    for mol in [
        with_group(hs.clone(), &[0], &[0]),
        with_group(hs.clone(), &[], &[1]),
    ] {
        assert_groups_read(&remove_hs(&mol).unwrap());
    }
}

#[test]
fn other_removals_refuse_groups_of_bonds() {
    let mol = with_group(smiles("CCCO"), &[0], &[1]);
    assert_value_error(delete_substructs(&mol, &smiles("O")));
    let atoms_only = with_group(smiles("CCCO"), &[0], &[]);
    let deleted = delete_substructs(&atoms_only, &smiles("O")).unwrap();
    assert_eq!(groups(&deleted), [(vec![0], vec![])]);
}
