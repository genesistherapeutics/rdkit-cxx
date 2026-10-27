// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Inputs for which RDKit 2024.03's editing calls would reach undefined
// behavior: stereo groups left pointing at removed bonds or at a repeated
// atom's second copy, molzip's linkers and coordinate lookups, replaceCore's
// match lookup, and fragmentOnSomeBonds' shift. Each is refused with
// Error::Value or handled as RDKit 2026.03 does, and the ordinary uses next
// to it still work.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn is_value<T: std::fmt::Debug>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

// An atropisomer whose AND group holds its chiral atom 8 and, as RDKit moves
// an atropisomer's atoms into its bond, its atropisomeric bond 2 (atoms 2
// and 3). `GROUP` is the group's CXSMILES field.
fn atropisomer(group: &str) -> RWMol {
    RWMol::from(smiles(&format!(
        "CC1=C(N2C=CC=C2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;\
         0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;\
         -0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;\
         -3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;\
         0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,{group}|"
    )))
}

fn groups(mol: &Mol) -> Vec<(Vec<u32>, Vec<u32>)> {
    mol.get_stereo_groups()
        .iter()
        .map(|g| (g.get_atoms().to_vec(), g.get_bonds().to_vec()))
        .collect()
}

// Reads every group's atoms and bonds back through RDKit: copying,
// pickling, and writing.
fn read_groups(mol: &Mol) {
    let copy = mol.clone();
    assert_eq!(groups(&copy), groups(mol));
    // SAFETY: a pickle RDKit wrote.
    let unpickled = unsafe { Mol::from_pkl_string(&mol.to_binary()) }.unwrap();
    assert_eq!(groups(&unpickled), groups(mol));
    // Breaking the pyrrole's bond leaves it unkekulizable.
    let params = MolWriterParams {
        kekulize: false,
        ..Default::default()
    };
    mol_to_v3k_mol_block_with_params(mol, &params).unwrap();
}

#[test]
fn atropisomer_groups_hold_their_bond() {
    let mol = atropisomer("&1:2,8");
    assert_eq!(groups(&mol), [(vec![8], vec![2])]);
    assert_eq!(groups(&atropisomer("&1:2")), [(vec![], vec![2])]);
    read_groups(&mol);
}

#[test]
fn removing_a_bond_takes_it_out_of_groups() {
    let mut mol = atropisomer("&1:2,8");
    mol.remove_bond(2, 3).unwrap();
    assert_eq!(groups(&mol), [(vec![8], vec![])]);
    read_groups(&mol);

    // A group left with neither atoms nor bonds goes, as in RDKit 2026.03.
    let mut mol = atropisomer("&1:2");
    mol.remove_bond(2, 3).unwrap();
    assert!(groups(&mol).is_empty());

    // Removing another bond leaves the group's.
    let mut mol = atropisomer("&1:2,8");
    mol.remove_bond(0, 1).unwrap();
    assert_eq!(groups(&mol), [(vec![8], vec![1])]);
    read_groups(&mol);
}

#[test]
fn removing_an_atom_takes_its_bonds_out_of_groups() {
    let mut mol = atropisomer("&1:2,8");
    mol.remove_atom(3).unwrap();
    assert_eq!(groups(&mol), [(vec![7], vec![])]);
    read_groups(&mol);

    // RDKit 2024.03 drops a group of bonds alone whenever it removes an
    // atom; it stays, as in RDKit 2026.03, while its bonds do.
    let mut mol = atropisomer("&1:2");
    mol.remove_atom(0).unwrap();
    assert_eq!(groups(&mol), [(vec![], vec![1])]);
    read_groups(&mol);

    let mut mol = atropisomer("&1:2,8");
    mol.remove_atom(8).unwrap();
    assert_eq!(groups(&mol), [(vec![], vec![2])]);
    read_groups(&mol);
}

#[test]
fn committing_a_batch_takes_its_removals_out_of_groups() {
    let mut mol = atropisomer("&1:2,8");
    mol.begin_batch_edit().unwrap();
    mol.remove_bond(2, 3).unwrap();
    mol.remove_atom(15).unwrap();
    // Nothing is removed until the commit.
    assert_eq!(groups(&mol), [(vec![8], vec![2])]);
    mol.commit_batch_edit().unwrap();
    assert_eq!(groups(&mol), [(vec![8], vec![])]);
    read_groups(&mol);

    // The atoms' bonds go too.
    let mut mol = atropisomer("&1:2,8");
    mol.begin_batch_edit().unwrap();
    mol.remove_atom(2).unwrap();
    mol.commit_batch_edit().unwrap();
    assert_eq!(groups(&mol), [(vec![7], vec![])]);
    read_groups(&mol);

    let mut mol = atropisomer("&1:2,8");
    mol.begin_batch_edit().unwrap();
    mol.remove_bond(2, 3).unwrap();
    mol.rollback_batch_edit();
    assert_eq!(groups(&mol), [(vec![8], vec![2])]);
}

#[test]
fn rdkit_removals_refuse_groups_with_bonds() {
    let mol = atropisomer("&1:2,8");
    assert!(is_value(fragment_on_bonds(&mol, &[0])));
    assert!(is_value(delete_substructs(
        &mol,
        &mol_from_smarts("Cl").unwrap()
    )));
    assert!(is_value(fragment_on_some_bonds(&mol, &[0, 1])));
    assert!(is_value(molzip_single(&mol)));

    // A group of atoms alone is kept, and fixed by RDKit's removal.
    let plain = atropisomer("&1:8");
    assert_eq!(groups(&plain), [(vec![8], vec![])]);
    let fragmented = fragment_on_bonds(&plain, &[0]).unwrap();
    assert_eq!(groups(&fragmented), [(vec![8], vec![])]);
    delete_substructs(&plain, &mol_from_smarts("Cl").unwrap()).unwrap();
}

#[test]
fn stereo_groups_refuse_repeated_atoms() {
    let mut mol = RWMol::from(smiles("C[C@H](F)Cl"));
    assert!(is_value(create_stereo_group(
        StereoGroupType::StereoOr,
        &mol,
        &[1, 1]
    )));
    let group = create_stereo_group(StereoGroupType::StereoOr, &mol, &[1]).unwrap();
    mol.set_stereo_groups(&[group]).unwrap();
    assert_eq!(groups(&mol), [(vec![1], vec![])]);
}

#[test]
fn merged_abs_groups_are_refused_or_scrubbed() {
    // Two ABS groups sharing an atom are allowed, but insertMol merges them
    // into one that repeats it.
    let mut mol = RWMol::from(smiles("C[C@H](F)Cl"));
    let abs = create_stereo_group(StereoGroupType::StereoAbsolute, &mol, &[1]).unwrap();
    mol.set_stereo_groups(&[abs.clone(), abs]).unwrap();
    // RDKit's own removals, which may merge them first, refuse them.
    assert!(is_value(fragment_on_bonds(&mol, &[0])));
    assert!(is_value(molzip(&mol, &smiles("C"))));

    mol.insert_mol(&smiles("CC"));
    let merged = mol.get_stereo_groups();
    assert_eq!(groups(&mol), [(vec![1, 1], vec![])]);
    assert!(is_value(fragment_on_bonds(&mol, &[0])));
    // Setting the merged group back is refused, as RDKit 2026.03 does.
    let mut other = mol.clone();
    assert!(is_value(other.set_stereo_groups(&merged)));

    // Removing the atom takes both copies out.
    mol.remove_atom(1).unwrap();
    assert!(groups(&mol).is_empty());
    read_groups(&mol);

    let mut mol = RWMol::from(smiles("C[C@H](F)Cl"));
    let abs = create_stereo_group(StereoGroupType::StereoAbsolute, &mol, &[1]).unwrap();
    mol.set_stereo_groups(&[abs.clone(), abs]).unwrap();
    mol.insert_mol(&smiles("CC"));
    mol.begin_batch_edit().unwrap();
    mol.remove_atom(1).unwrap();
    mol.commit_batch_edit().unwrap();
    assert!(groups(&mol).is_empty());
}

fn to_smiles(mol: &Mol) -> String {
    mol_to_smiles(mol).unwrap()
}

fn fob() -> MolzipParams {
    MolzipParams {
        label: MolzipLabel::FragmentOnBonds,
        ..Default::default()
    }
}

fn with_coordinates(label: MolzipLabel) -> MolzipParams {
    MolzipParams {
        label,
        generate_coordinates: true,
        ..Default::default()
    }
}

#[test]
fn molzip_refuses_linkers() {
    assert!(is_value(molzip_single(&smiles("C[*:1].[*:1][*:2].[*:2]N"))));
    assert!(is_value(molzip(
        &smiles("C[*:1].[*:2]N"),
        &smiles("[*:1][*:2]")
    )));
    let isotope = MolzipParams {
        label: MolzipLabel::Isotope,
        ..Default::default()
    };
    assert!(is_value(molzip_single_with_params(
        &smiles("C[1*].[1*][2*].[2*]N"),
        &isotope
    )));
    let mut atom_type = MolzipParams {
        label: MolzipLabel::AtomType,
        ..Default::default()
    };
    atom_type.set_atom_symbols(&["Xe", "Kr"]);
    assert!(is_value(molzip_single_with_params(
        &smiles("C[Xe].[Xe][Kr].[Kr]N"),
        &atom_type
    )));
    // With FragmentOnBonds, two dummies bonded to each other.
    assert!(is_value(molzip_single_with_params(
        &smiles("[1*][1*]"),
        &fob()
    )));
    // A dummy whose isotope names a dummy finds no dummy beside it naming
    // its own atom back, so RDKit gives up before it deletes anything.
    assert!(
        molzip_with_params(&smiles("C[2*]"), &smiles("[0*]N"), &fob())
            .unwrap()
            .is_none()
    );
    let mut list = vec![smiles("C[*:1]"), smiles("[*:1][*:2]"), smiles("N[*:2]")];
    assert!(is_value(molzip_fragments(&mut list)));

    // Ordinary zips.
    assert_eq!(
        to_smiles(&molzip_single(&smiles("C[*:1].N[*:1]")).unwrap().unwrap()),
        "CN"
    );
    assert_eq!(
        to_smiles(
            &molzip(&smiles("C[*:1]"), &smiles("[*:1]C[*:2].N[*:2]"))
                .unwrap()
                .unwrap()
        ),
        "CCN"
    );
    let fragmented = fragment_on_bonds(&smiles("CCN"), &[1]).unwrap();
    assert_eq!(
        to_smiles(
            &molzip_single_with_params(&fragmented, &fob())
                .unwrap()
                .unwrap()
        ),
        "CCN"
    );
}

#[test]
fn molzip_lists_with_coordinates() {
    // RDKit cannot zip this one (the dummy names atom 0, which has no dummy
    // naming it back): no result, rather than a null dereference.
    let mut list = vec![smiles("CC*")];
    assert!(
        molzip_fragments_with_params(&mut list, &fob())
            .unwrap()
            .is_none()
    );
    let mut list = vec![smiles("CC*")];
    assert!(
        molzip_fragments_with_params(&mut list, &with_coordinates(MolzipLabel::FragmentOnBonds))
            .unwrap()
            .is_none()
    );

    // A fragment after the first whose SMILES repeats one's.
    let repeated = || vec![smiles("C[*:1]"), smiles("N[*:1]"), smiles("N[*:1]")];
    assert!(is_value(molzip_fragments_with_params(
        &mut repeated(),
        &with_coordinates(MolzipLabel::AtomMapNumber)
    )));
    let mut list = repeated();
    assert_eq!(
        to_smiles(&molzip_fragments(&mut list).unwrap().unwrap()),
        "CN"
    );
    let mut row = BTreeMap::from([
        ("Core".to_owned(), smiles("C[*:1]")),
        ("R1".to_owned(), smiles("N[*:1]")),
        ("R2".to_owned(), smiles("N[*:1]")),
    ]);
    assert!(is_value(molzip_row_with_params(
        &mut row,
        &with_coordinates(MolzipLabel::AtomMapNumber)
    )));

    // More than one molecule with FragmentOnBonds.
    let fragmented = fragment_on_bonds(&smiles("CCN"), &[1]).unwrap();
    let mut list = vec![fragmented.clone(), smiles("CC")];
    assert!(is_value(molzip_fragments_with_params(
        &mut list,
        &with_coordinates(MolzipLabel::FragmentOnBonds)
    )));
    let mut list = vec![fragmented];
    let zipped =
        molzip_fragments_with_params(&mut list, &with_coordinates(MolzipLabel::FragmentOnBonds))
            .unwrap()
            .unwrap();
    assert_eq!(to_smiles(&zipped), "CCN");
    assert_eq!(zipped.get_num_conformers(), 1);
    assert_eq!(list[0].get_num_conformers(), 1);

    // Ordinary lists get coordinates.
    let mut list = vec![smiles("c1ccccc1[*:1]"), smiles("[*:1]CC")];
    let zipped =
        molzip_fragments_with_params(&mut list, &with_coordinates(MolzipLabel::AtomMapNumber))
            .unwrap()
            .unwrap();
    assert_eq!(to_smiles(&zipped), "CCc1ccccc1");
    assert!(list.iter().all(|mol| mol.get_num_conformers() == 1));
    let mut row = BTreeMap::from([
        ("Core".to_owned(), smiles("c1ccccc1[*:1]")),
        ("R1".to_owned(), smiles("[*:1]CC")),
    ]);
    let zipped = molzip_row_with_params(&mut row, &with_coordinates(MolzipLabel::AtomMapNumber))
        .unwrap()
        .unwrap();
    assert_eq!(zipped.get_num_conformers(), 1);
}

#[test]
fn replace_core_refuses_unmatched_neighbors() {
    let mol = smiles("CCCC");
    let core = smiles("CC");
    assert!(is_value(replace_core(&mol, &core, &[(0, 1), (0, 1)])));
    // An index out of range is RDKit's own error.
    assert!(is_value(replace_core(&mol, &core, &[(0, 1), (5, 1)])));
    assert_eq!(
        to_smiles(
            &replace_core(&mol, &core, &[(0, 1), (1, 2)])
                .unwrap()
                .unwrap()
        ),
        "[1*]C.[2*]C"
    );
}

#[test]
fn fragment_on_some_bonds_without_enough_bonds() {
    let mol = smiles("CCCC");
    for num_to_break in [0, 3, 10, 63, 64, 65, 1000, u32::MAX] {
        let params = FragmentOnSomeBondsParams {
            num_to_break,
            ..Default::default()
        };
        assert!(
            fragment_on_some_bonds_with_params(&mol, &[0, 1], &params)
                .unwrap()
                .is_empty()
        );
        let (frags, cuts) =
            fragment_on_some_bonds_return_cuts_per_atom_with_params(&mol, &[0, 1], &params)
                .unwrap();
        assert!(frags.is_empty() && cuts.is_empty());
    }
    let params = FragmentOnSomeBondsParams {
        num_to_break: 2,
        ..Default::default()
    };
    assert_eq!(
        fragment_on_some_bonds_with_params(&mol, &[0, 1], &params)
            .unwrap()
            .len(),
        1
    );
    assert_eq!(fragment_on_some_bonds(&mol, &[0, 1]).unwrap().len(), 2);
    // The other checks come first.
    let params = FragmentOnSomeBondsParams {
        num_to_break: 0,
        ..Default::default()
    };
    assert!(is_value(fragment_on_some_bonds_with_params(
        &mol,
        &[],
        &params
    )));
    assert!(is_value(fragment_on_some_bonds_with_params(
        &mol,
        &[0, 0],
        &params
    )));
}

#[test]
fn fragmenting_refuses_repeated_bonds() {
    let mol = smiles("CCCC");
    assert!(is_value(fragment_on_bonds(&mol, &[1, 1])));
    assert!(is_value(fragment_on_some_bonds(&mol, &[0, 1, 0])));
    assert_eq!(
        to_smiles(&fragment_on_bonds(&mol, &[1]).unwrap()),
        "[1*]CC.[2*]CC"
    );
}
