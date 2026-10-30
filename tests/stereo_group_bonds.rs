// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// RWMol::replaceBond deletes the bond it replaces and leaves stereo groups
// pointing at it. Each test runs one path that replaces a bond a stereo
// group holds, and checks that the group holds the new bond.

use rdkit_cxx::chem::rd_chem_reactions::*;
use rdkit_cxx::chem::rd_depictor::*;
use rdkit_cxx::chem::rd_generalized_substruct::*;
use rdkit_cxx::chem::rd_mol_align::*;
use rdkit_cxx::chem::rd_tautomer_query::TautomerQuery;
use rdkit_cxx::chem::rdqueries::atom_num_equals_query_atom;
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

// An atropisomer whose AND group holds its chiral atom 8 and, as RDKit
// moves an atropisomer's atoms into its bond, its atropisomeric bond 2
// (atoms 2 and 3).
const ATROPISOMER: &str = "CC1=C(N2C=CC=C2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;\
    0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;\
    -1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;\
    -2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,&1:2,8|";

// An AND group of atom 0 and `bond`. RDKit 2024.03's CreateStereoGroup
// takes no bonds, so the group is the atropisomer's, moved to these
// indices: its atom first (renumber_atoms keeps the bonds' order), and the
// bonds before its bond removed, or a chain's bonds inserted before them.
fn group_of_atom_0_and(bond: u32) -> StereoGroup {
    let atropisomer = smiles(ATROPISOMER);
    let mut donor = if bond > 2 {
        let mut chain = RWMol::from(smiles(&"C".repeat(bond as usize - 1)));
        chain.insert_mol(&atropisomer);
        chain
    } else {
        RWMol::from(atropisomer)
    };
    let atom = donor.get_stereo_groups()[0].get_atoms()[0];
    let order: Vec<u32> = std::iter::once(atom)
        .chain((0..donor.get_num_atoms()).filter(|&i| i != atom))
        .collect();
    donor = RWMol::from(renumber_atoms(&donor, &order).unwrap());
    for _ in bond..2 {
        let first = donor.get_bond_with_idx(0).unwrap();
        let (begin, end) = (first.get_begin_atom_idx(), first.get_end_atom_idx());
        donor.remove_bond(begin, end).unwrap();
    }
    let group = donor.get_stereo_groups().remove(0);
    assert_eq!(group.get_atoms(), [0]);
    assert_eq!(group.get_bonds(), [bond]);
    group
}

// The molecule with an AND group of atom 0 and `bond`.
fn with_group(mol: Mol, bond: u32) -> RWMol {
    let mut rw = RWMol::from(mol);
    rw.set_stereo_groups(&[group_of_atom_0_and(bond)]).unwrap();
    rw
}

fn group_bonds(mol: &Mol) -> Vec<Vec<u32>> {
    mol.get_stereo_groups()
        .iter()
        .map(|g| g.get_bonds().to_vec())
        .collect()
}

// Removing the bond takes it out of a group that holds it, which a group
// still pointing at the deleted bond would keep.
fn assert_removal_updates_group(mol: Mol, begin: u32, end: u32) {
    let mut rw = RWMol::from(mol);
    rw.remove_bond(begin, end).unwrap();
    assert_eq!(group_bonds(&rw), [Vec::<u32>::new()]);
}

fn generic_bonds() -> AdjustQueryParameters {
    AdjustQueryParameters {
        make_bonds_generic: true,
        ..AdjustQueryParameters::no_adjustments()
    }
}

#[test]
fn replace_bond_moves_groups_to_the_new_bond() {
    let mut rw = with_group(smiles("CCCC"), 1);
    let donor = smiles("C=C");
    rw.replace_bond(1, &donor.get_bond_with_idx(0).unwrap())
        .unwrap();
    assert_eq!(group_bonds(&rw), [[1]]);
    assert_eq!(
        rw.get_bond_with_idx(1).unwrap().get_bond_type(),
        BondType::Double
    );
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { Mol::from_pkl_string(&rw.to_binary()) }.unwrap();
    assert_eq!(group_bonds(&copy), [[1]]);
    assert!(mol_to_v3k_mol_block(&rw).unwrap().contains("MDLV30/STERAC"));
    assert_removal_updates_group(rw.into(), 1, 2);
}

#[test]
fn adjusting_queries_keeps_groups() {
    let mol = with_group(smiles("CCCC"), 1);
    let adjusted = adjust_query_properties_with_params(&mol, &generic_bonds()).unwrap();
    assert!(adjusted.get_bond_with_idx(1).unwrap().has_query());
    assert_eq!(group_bonds(&adjusted), [[1]]);
    assert_eq!(group_bonds(&mol), [[1]]);
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { Mol::from_pkl_string(&adjusted.to_binary()) }.unwrap();
    assert_eq!(group_bonds(&copy), [[1]]);
    assert_removal_updates_group(adjusted, 1, 2);

    let generic =
        adjust_query_properties_with_generic_groups_with_params(&mol, &generic_bonds()).unwrap();
    assert!(generic.get_bond_with_idx(1).unwrap().has_query());
    assert_eq!(group_bonds(&generic), [[1]]);
    assert_removal_updates_group(generic, 1, 2);
}

#[test]
fn extended_queries_keep_groups() {
    // Acetone's tautomers change bonds 0, 1, and 2, which the template
    // replaces with queries.
    let mol = with_group(smiles("CC(=O)C"), 1);
    let params = CreateExtendedQueryMolParams {
        do_enumeration: false,
        do_tautomers: true,
        ..Default::default()
    };
    let query = create_extended_query_mol_with_params(&mol, &params).unwrap();
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { ExtendedQueryMol::from_data(&query.to_binary().unwrap()) }.unwrap();
    assert!(!query.to_json().unwrap().is_empty());
    for q in [&query, &copy] {
        assert!(mol_has_substruct_match(&smiles("CC(O)=C"), q).unwrap());
    }

    let adjusted = CreateExtendedQueryMolParams {
        do_enumeration: false,
        do_tautomers: false,
        adjust_query_properties: true,
        adjust_query_parameters: generic_bonds(),
    };
    let query = create_extended_query_mol_with_params(&mol, &adjusted).unwrap();
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { ExtendedQueryMol::from_data(&query.to_binary().unwrap()) }.unwrap();
    for q in [&query, &copy, &query.clone()] {
        assert!(mol_has_substruct_match(&smiles("C=C(C)O"), q).unwrap());
    }
}

#[test]
fn tautomer_query_templates_keep_groups() {
    let mol = with_group(smiles("CC(=O)C"), 1);
    let query = TautomerQuery::new(&mol).unwrap();
    assert!(query.get_modified_bonds().contains(&1));
    let template = query.get_template_molecule();
    assert!(template.get_bond_with_idx(1).unwrap().has_query());
    assert_eq!(group_bonds(&template), [[1]]);
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { TautomerQuery::from_pickle(&query.to_binary().unwrap()) }.unwrap();
    assert_eq!(group_bonds(&copy.get_template_molecule()), [[1]]);
    assert_removal_updates_group(template, 1, 2);
}

#[test]
fn sanitizing_reactions_keeps_template_groups() {
    let mut rxn = ChemicalReaction::default();
    let template = with_group(mol_from_smarts("[C:1]-[C:2]-[C:3]-[C:4]").unwrap(), 0);
    rxn.add_reactant_template(&template);
    rxn.add_product_template(&mol_from_smarts("[C:1]-[C:2].[C:3]-[C:4]").unwrap());
    let params = SanitizeRxnParams {
        params: generic_bonds(),
        ..Default::default()
    };
    sanitize_rxn_with_params(&mut rxn, &params).unwrap();
    let reactant = rxn.get_reactant_template(0).unwrap();
    assert_eq!(group_bonds(&reactant), [[0]]);
    // SAFETY: a pickle RDKit wrote.
    let copy = unsafe { ChemicalReaction::from_bin_str(&rxn.to_binary()) }.unwrap();
    assert_eq!(group_bonds(&copy.get_reactant_template(0).unwrap()), [[0]]);
}

#[test]
fn depicting_with_r_groups_leaves_groups() {
    // The single bond from the aromatic query atom to the R group is one
    // RDKit's R-group handling replaces with a query.
    let mut reference = RWMol::from(smiles("c1ccccc1*"));
    compute_2d_coords(&mut reference).unwrap();
    let mut atom = atom_num_equals_query_atom(6).unwrap();
    atom.set_is_aromatic(true);
    reference.replace_atom(5, &atom).unwrap();
    let reference = with_group(reference.into(), 5);
    let params = ConstrainedDepictionParams {
        allow_r_groups: true,
        ..Default::default()
    };
    let mut mol = smiles("c1ccccc1C");
    let pairs = generate_depiction_matching_2d_structure_with_params(&mut mol, &reference, &params)
        .unwrap();
    assert_eq!(pairs.len(), 7);
    assert_eq!(group_bonds(&reference), [[5]]);

    // As a pattern, mapping the molecule to a plain reference.
    let mut target = smiles("c1ccccc1C");
    compute_2d_coords(&mut target).unwrap();
    let pattern = ConstrainedDepictionParams {
        ref_patt: Some(reference.get_mol()),
        ..params
    };
    let mut mol = smiles("c1ccccc1CC");
    assert_eq!(
        generate_depiction_matching_2d_structure_with_params(&mut mol, &target, &pattern)
            .unwrap()
            .len(),
        7
    );
    assert_eq!(group_bonds(pattern.ref_patt.as_ref().unwrap()), [[5]]);
}

#[test]
fn symmetrized_alignment_leaves_groups() {
    // Symmetrizing the carboxylate replaces its C-O bonds with queries.
    let mut reference = smiles("CC(=O)[O-]");
    compute_2d_coords(&mut reference).unwrap();
    let mut probe = with_group(reference.clone(), 2);
    let params = GetBestRMSParams {
        symmetrize_conjugated_terminal_groups: true,
        ..Default::default()
    };
    let rms = get_best_rms_with_params(&mut probe, &reference, &params).unwrap();
    assert!(rms < 1e-6);
    assert_eq!(group_bonds(&probe), [[2]]);
}
