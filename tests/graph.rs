// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

#[test]
fn distance_matrix() {
    // Values from Python's rdmolops.GetDistanceMatrix.
    let mol = smiles("CC(=O)O");
    assert_eq!(
        get_distance_matrix(&mol).unwrap(),
        vec![
            vec![0.0, 1.0, 2.0, 2.0],
            vec![1.0, 0.0, 1.0, 1.0],
            vec![2.0, 1.0, 0.0, 2.0],
            vec![2.0, 1.0, 2.0, 0.0],
        ]
    );
    let bo = GetDistanceMatrixParams {
        use_bo: true,
        ..Default::default()
    };
    assert_eq!(
        get_distance_matrix_with_params(&mol, &bo).unwrap()[0],
        vec![0.0, 1.0, 1.5, 2.0]
    );
    let weights = GetDistanceMatrixParams {
        use_atom_wts: true,
        ..Default::default()
    };
    assert_eq!(
        get_distance_matrix_with_params(&mol, &weights).unwrap()[2][2],
        6.0 / 8.0
    );
    assert!(get_distance_matrix(&smiles("")).unwrap().is_empty());
}

#[test]
fn matrix_params_defaults() {
    assert_eq!(
        GetDistanceMatrixParams::default(),
        GetDistanceMatrixParams {
            use_bo: false,
            use_atom_wts: false,
        }
    );
    assert_eq!(
        Get3DDistanceMatrixParams::default(),
        Get3DDistanceMatrixParams {
            conf_id: -1,
            use_atom_wts: false,
        }
    );
}

#[test]
fn distance_matrix_after_graph_edits() {
    // RDKit caches the matrix in a computed property that RWMol's edits
    // keep, so a stale 3x3 matrix would be read for the bigger molecule.
    let mut mol = RWMol::from(smiles("CCO"));
    assert_eq!(get_distance_matrix(&mol).unwrap().len(), 3);
    let idx = mol.add_atom(&Atom::new("N").unwrap());
    mol.add_bond_with_order(2, idx, BondType::Single).unwrap();
    let matrix = get_distance_matrix(&mol).unwrap();
    assert_eq!(matrix.len(), 4);
    assert_eq!(matrix[0][3], 3.0);
    mol.remove_bond(1, 2).unwrap();
    assert_eq!(get_distance_matrix(&mol).unwrap()[0][3], 1e8);
    // findAllPathsOfLengthN(onlyShortestPaths) builds paths from the cached
    // matrix and dereferences the bond between each pair.
    let shortest = FindAllPathsOfLengthNParams {
        only_shortest_paths: true,
        ..Default::default()
    };
    assert_eq!(
        find_all_paths_of_length_n_with_params(&mol, 2, &shortest).unwrap(),
        Vec::<Vec<u32>>::new()
    );
    let mut combined = RWMol::from(smiles("CC"));
    assert_eq!(get_distance_matrix(&combined).unwrap().len(), 2);
    combined.insert_mol(&smiles("OO"));
    assert_eq!(get_distance_matrix(&combined).unwrap().len(), 4);
}

#[test]
fn distance_matrix_3d() {
    let mol = mol_from_mol_block(
        "
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    4.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
",
    )
    .unwrap();
    assert_eq!(
        get_3d_distance_matrix(&mol).unwrap(),
        vec![vec![0.0, 5.0], vec![5.0, 0.0]]
    );
    let params = Get3DDistanceMatrixParams {
        use_atom_wts: true,
        ..Default::default()
    };
    assert_eq!(
        get_3d_distance_matrix_with_params(&mol, &params).unwrap()[0][0],
        1.0
    );
    let missing = Get3DDistanceMatrixParams {
        conf_id: 7,
        ..Default::default()
    };
    assert!(get_3d_distance_matrix_with_params(&mol, &missing).is_err());
    assert!(get_3d_distance_matrix(&smiles("CO")).is_err());
}

#[test]
fn adjacency_matrix() {
    let mol = smiles("CC(=O)O");
    assert_eq!(
        get_adjacency_matrix(&mol).unwrap(),
        vec![
            vec![0, 1, 0, 0],
            vec![1, 0, 1, 1],
            vec![0, 1, 0, 0],
            vec![0, 1, 0, 0],
        ]
    );
    assert_eq!(
        get_adjacency_matrix_with_empty_val(&mol, -1).unwrap()[0],
        vec![-1, 1, -1, -1]
    );
    assert_eq!(
        get_adjacency_matrix_use_bo(&mol).unwrap()[1],
        vec![1.0, 0.0, 2.0, 1.0]
    );
    assert_eq!(
        get_adjacency_matrix_use_bo_with_empty_val(&mol, 9).unwrap()[2],
        vec![9.0, 2.0, 9.0, 9.0]
    );
}

#[test]
fn paths() {
    let mol = smiles("CCCC");
    assert_eq!(
        find_all_paths_of_length_n(&mol, 2).unwrap(),
        vec![vec![0, 1], vec![1, 2]]
    );
    let atoms = FindAllPathsOfLengthNParams {
        use_bonds: false,
        ..Default::default()
    };
    assert_eq!(
        find_all_paths_of_length_n_with_params(&mol, 2, &atoms).unwrap(),
        vec![vec![0, 1], vec![1, 2], vec![2, 3]]
    );
    // RDKit would loop to the requested length, and add one to it.
    assert!(
        find_all_paths_of_length_n(&mol, u32::MAX)
            .unwrap()
            .is_empty()
    );
    assert!(
        find_all_paths_of_length_n_with_params(&mol, u32::MAX, &atoms)
            .unwrap()
            .is_empty()
    );
    // A ring closure makes a path of numAtoms bonds.
    assert_eq!(
        find_all_paths_of_length_n(&smiles("C1CC1"), 3)
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        FindAllPathsOfLengthNParams::default(),
        FindAllPathsOfLengthNParams {
            use_bonds: true,
            use_hs: false,
            rooted_at_atom: -1,
            only_shortest_paths: false,
        }
    );
}

#[test]
fn subgraphs() {
    // Values from Python's rdmolops.
    let mol = smiles("C1CC1C");
    assert_eq!(
        find_all_subgraphs_of_length_n(&mol, 2).unwrap(),
        vec![vec![0, 1], vec![0, 3], vec![1, 3], vec![1, 2], vec![2, 3]]
    );
    assert!(
        find_all_subgraphs_of_length_n(&mol, u32::MAX)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        find_all_subgraphs_of_length_m_to_n(&mol, 1, 6).unwrap(),
        vec![
            vec![vec![0], vec![1], vec![2], vec![3]],
            vec![vec![0, 1], vec![0, 3], vec![1, 3], vec![1, 2], vec![2, 3]],
            vec![vec![0, 1, 3], vec![0, 1, 2], vec![0, 3, 2], vec![1, 3, 2]],
            vec![vec![0, 1, 3, 2]],
            vec![],
            vec![],
        ]
    );
    assert_eq!(
        find_all_subgraphs_of_length_m_to_n(&mol, 9, 10).unwrap(),
        vec![Vec::<Vec<u32>>::new(); 2]
    );
    assert!(matches!(
        find_all_subgraphs_of_length_m_to_n(&mol, 2, 1),
        Err(Error::Value { .. })
    ));
    assert_eq!(
        find_unique_subgraphs_of_length_n(&mol, 2).unwrap(),
        vec![vec![0, 1]]
    );
    assert_eq!(
        FindUniqueSubgraphsOfLengthNParams::default(),
        FindUniqueSubgraphsOfLengthNParams {
            use_hs: false,
            use_bo: true,
            rooted_at_atom: -1,
        }
    );
    assert_eq!(
        FindAllSubgraphsOfLengthNParams::default(),
        FindAllSubgraphsOfLengthNParams {
            use_hs: false,
            rooted_at_atom: -1,
        }
    );
    assert_eq!(
        FindAllSubgraphsOfLengthMToNParams::default(),
        FindAllSubgraphsOfLengthMToNParams {
            use_hs: false,
            rooted_at_atom: -1,
        }
    );
}

#[test]
fn atom_environment() {
    let mol = smiles("CCCC");
    assert_eq!(
        find_atom_environment_of_radius_n(&mol, 2, 1).unwrap(),
        vec![0, 1, 2]
    );
    let (path, outputs) = find_atom_environment_of_radius_n_with_outputs(
        &mol,
        2,
        1,
        &FindAtomEnvironmentOfRadiusNParams::default(),
    )
    .unwrap();
    assert_eq!(path, vec![0, 1, 2]);
    assert_eq!(
        outputs.atom_map,
        BTreeMap::from([(0, 1), (1, 0), (2, 1), (3, 2)])
    );
    assert!(
        find_atom_environment_of_radius_n(&mol, 9, 0)
            .unwrap()
            .is_empty()
    );
    let loose = FindAtomEnvironmentOfRadiusNParams {
        enforce_size: false,
        ..Default::default()
    };
    assert_eq!(
        find_atom_environment_of_radius_n_with_params(&mol, u32::MAX, 0, &loose).unwrap(),
        vec![0, 1, 2]
    );
    assert!(matches!(
        find_atom_environment_of_radius_n(&mol, 1, 4),
        Err(Error::Value { .. })
    ));
}

#[test]
fn path_to_submol() {
    let mol = smiles("CCCC");
    let (submol, outputs) = path_to_submol_with_outputs(&mol, &[1, 2], false).unwrap();
    assert_eq!(mol_to_smiles(&submol).unwrap(), "CCC");
    assert_eq!(outputs.atom_map, BTreeMap::from([(1, 0), (2, 1), (3, 2)]));
    // Bond indices past the end are refused, as RDKit 2025.03 throws.
    assert!(matches!(
        rdkit_cxx::chem::path_to_submol(&mol, &[0, 99]),
        Err(rdkit_cxx::Error::Index { index: 99, len: 3 })
    ));
    let query = path_to_submol_with_use_query(&mol, &[0], true).unwrap();
    assert!(query.has_query());
}

#[test]
fn split_by_pdb() {
    let mut mol = smiles("CCO.N");
    let mut info = AtomPDBResidueInfo::new(" C1 ");
    info.set_residue_name("ALA");
    info.set_chain_id("A");
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_pdb_residue_info(&info);
    mol.get_atom_with_idx_mut(1)
        .unwrap()
        .set_pdb_residue_info(&info);
    let residues = split_mol_by_pdb_residues(&mol).unwrap();
    assert_eq!(residues.keys().collect::<Vec<_>>(), vec!["", "ALA"]);
    assert_eq!(mol_to_smiles(&residues["ALA"]).unwrap(), "CC");
    assert_eq!(mol_to_smiles(&residues[""]).unwrap(), "N.O");
    let chains = split_mol_by_pdb_chain_id(&mol).unwrap();
    assert_eq!(chains.keys().collect::<Vec<_>>(), vec!["", "A"]);
    let only = SplitMolByPDBResiduesParams {
        white_list: Some(vec!["ALA".to_owned()]),
        negate_list: false,
    };
    let residues = split_mol_by_pdb_residues_with_params(&mol, &only).unwrap();
    assert_eq!(residues.keys().collect::<Vec<_>>(), vec!["ALA"]);
    let others = SplitMolByPDBChainIdParams {
        white_list: Some(vec!["A".to_owned()]),
        negate_list: true,
    };
    let chains = split_mol_by_pdb_chain_id_with_params(&mol, &others).unwrap();
    assert_eq!(chains.keys().collect::<Vec<_>>(), vec![""]);
}

#[test]
fn split_by_pdb_with_other_monomer_types() {
    // As in Python, only info of type PDBRESIDUE has a residue, even PDB
    // residue info whose type was changed.
    let mut mol = smiles("CCC");
    let mut info = AtomPDBResidueInfo::new(" C1 ");
    info.set_residue_name("ALA");
    info.set_chain_id("A");
    info.set_monomer_type(AtomMonomerType::Other);
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_pdb_residue_info(&info);
    let plain = AtomMonomerInfo::new(AtomMonomerType::Other).unwrap();
    mol.get_atom_with_idx_mut(1)
        .unwrap()
        .set_monomer_info(&plain);
    let residues = split_mol_by_pdb_residues(&mol).unwrap();
    assert_eq!(residues.keys().collect::<Vec<_>>(), vec![""]);
    let chains = split_mol_by_pdb_chain_id(&mol).unwrap();
    assert_eq!(chains.keys().collect::<Vec<_>>(), vec![""]);
}

#[test]
fn renumber() {
    let mol = smiles("CC(=O)O");
    let params = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    assert_eq!(
        mol_to_smiles_with_params(&renumber_atoms(&mol, &[3, 2, 1, 0]).unwrap(), &params).unwrap(),
        "OC(=O)C"
    );
    for bad in [
        &[0, 1, 2][..],
        &[0, 1, 2, 3, 4],
        &[0, 0, 1, 2],
        &[0, 1, 2, 9],
    ] {
        assert!(matches!(
            renumber_atoms(&mol, bad),
            Err(Error::Value { .. })
        ));
    }
}
