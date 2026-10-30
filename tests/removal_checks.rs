// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// RWMol's atom and bond removals parse each bond's _MolFileBondEndPts
// without checks and follow the substance groups' PARENT chains with no end
// for a cycle (CONVENTIONS.md P2). Every wrapper whose RDKit call removes
// atoms or bonds, from the molecule or a copy, refuses such a molecule
// before it starts, and leaves an in-place molecule as it was.

use rdkit_cxx::Error;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::*;
use rdkit_cxx::chem::rd_chem_reactions::*;
use rdkit_cxx::chem::rd_depictor::*;
use rdkit_cxx::chem::rd_deprotect::deprotect_in_place;
use rdkit_cxx::chem::rd_determine_bonds::determine_connectivity;
use rdkit_cxx::chem::rd_dist_geom::embed_molecule;
use rdkit_cxx::chem::rd_enumerate_stereoisomers::StereoisomerEnumerator;
use rdkit_cxx::chem::rd_force_field_helpers::uff_optimize_molecule;
use rdkit_cxx::chem::rd_generalized_substruct::create_extended_query_mol;
use rdkit_cxx::chem::rd_mmpa::fragment_mol;
use rdkit_cxx::chem::rd_mol_descriptors::bcut_2d;
use rdkit_cxx::chem::rd_mol_hash::{HashFunction, mol_hash};
use rdkit_cxx::chem::rd_r_group_decomposition::r_group_decompose;
use rdkit_cxx::chem::rd_rascal_mces::find_mces;
use rdkit_cxx::chem::rd_reduced_graphs::generate_mol_extended_reduced_graph;
use rdkit_cxx::chem::rd_substruct_library::CachedSmilesMolHolder;
use rdkit_cxx::chem::rd_tautomer_query::TautomerQuery;
use rdkit_cxx::chem::scaffolds::rd_scaffold_network::{
    ScaffoldNetworkParams, create_scaffold_network,
};
use rdkit_cxx::chem::*;

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

// Two substance groups whose PARENTs name each other, and a third that holds
// atom 0, which a removal of that atom drops, and RDKit then follows the
// other two's chains.
fn with_parent_cycle(mut mol: Mol) -> Mol {
    for (i, parent) in [(0, None), (1, Some(3)), (2, Some(2))] {
        let mut group = create_mol_substance_group(&mut mol, "DAT");
        group.add_atom_with_idx(0).unwrap();
        group.set_unsigned_prop("index", i + 1);
        if let Some(parent) = parent {
            group.set_unsigned_prop("PARENT", parent);
        }
    }
    mol
}

// A bond whose end points RDKit would read past their end.
fn with_bad_end_points(mut mol: Mol) -> Mol {
    mol.get_bond_with_idx_mut(0)
        .unwrap()
        .set_prop("_MolFileBondEndPts", "()");
    mol
}

fn assert_refused<T: std::fmt::Debug>(what: &str, res: Result<T, Error>) {
    match res {
        Err(Error::Value { .. }) => {}
        other => panic!("{what}: expected Error::Value, got {other:?}"),
    }
}

// Each malformed form of a molecule.
fn malformed(text: &str) -> [Mol; 2] {
    [
        with_parent_cycle(smiles(text)),
        with_bad_end_points(smiles(text)),
    ]
}

#[test]
fn rwmol_removals_refuse_before_editing() {
    for mol in malformed("CCO.N") {
        let before = mol_to_cx_smarts(&mol).unwrap();
        let mut rw = RWMol::from(mol);
        assert_refused("remove_atom", rw.remove_atom(0));
        assert_refused("remove_bond", rw.remove_bond(0, 1));
        let donor = smiles("C=C");
        let bond = donor.get_bond_with_idx(0).unwrap();
        let drop_groups = ReplaceBondParams {
            keep_s_groups: false,
            ..Default::default()
        };
        assert_refused(
            "replace_bond",
            rw.replace_bond_with_params(0, &bond, &drop_groups),
        );
        assert_eq!(mol_to_cx_smarts(&rw).unwrap(), before);
        rw.begin_batch_edit().unwrap();
        rw.remove_atom(0).unwrap_err();
        assert_refused("commit_batch_edit", rw.commit_batch_edit());
        // The batch is still open, and nothing was removed.
        assert!(rw.begin_batch_edit().is_err());
        rw.rollback_batch_edit();
        assert_eq!(rw.get_num_atoms(), 4);
        assert_eq!(mol_to_cx_smarts(&rw).unwrap(), before);
    }
    // A well-formed molecule still edits.
    let mut rw = RWMol::from(smiles("CCO.N"));
    rw.begin_batch_edit().unwrap();
    rw.remove_atom(3).unwrap();
    rw.commit_batch_edit().unwrap();
    assert_eq!(rw.get_num_atoms(), 3);
}

// The molecule with a substance group of atom 0, with this index and PARENT.
fn with_group(mut mol: Mol, index: u32, parent: u32) -> Mol {
    let mut group = create_mol_substance_group(&mut mol, "DAT");
    group.add_atom_with_idx(0).unwrap();
    group.set_unsigned_prop("index", index);
    group.set_unsigned_prop("PARENT", parent);
    mol
}

#[test]
fn combined_molecules_are_checked_together() {
    // Each molecule's chain ends, but combined, each PARENT finds the other
    // molecule's group.
    let mol = with_group(smiles("CC[*:1]"), 1, 2);
    let other = with_group(smiles("[*:1]O"), 2, 1);
    assert!(mol_to_smiles(&mol).is_ok() && mol_to_smiles(&other).is_ok());
    assert_refused("molzip", molzip(&mol, &other));
    assert!(molzip(&mol, &smiles("[*:1]O")).unwrap().is_some());
    let target = with_group(smiles("CCO"), 1, 2);
    let replacement = with_group(smiles("N"), 2, 1);
    assert_refused(
        "replace_substructs",
        replace_substructs(&target, &smiles("O"), &replacement),
    );
    assert_eq!(
        replace_substructs(&target, &smiles("O"), &smiles("N"))
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn molops_refuse() {
    for mol in malformed("CC[H].O") {
        assert_refused("remove_hs", remove_hs(&mol));
        assert_refused("remove_all_hs", remove_all_hs(&mol));
        assert_refused("merge_query_hs", merge_query_hs(&mol));
        assert_refused("get_mol_frags_as_mols", get_mol_frags_as_mols(&mol));
        assert_refused("haptic_bonds_to_dative", haptic_bonds_to_dative(&mol));
        assert_refused("dative_bonds_to_haptic", dative_bonds_to_haptic(&mol));
        let mut copy = mol.clone();
        assert_refused(
            "collapse_attachment_points",
            collapse_attachment_points(&mut copy),
        );
    }
    assert_eq!(get_mol_frags_as_mols(&smiles("CC.O")).unwrap().len(), 2);
}

#[test]
fn writers_refuse() {
    // A molecule of several fragments is written a fragment at a time.
    for mol in malformed("CC.O") {
        assert_refused("mol_to_smiles", mol_to_smiles(&mol));
        assert_refused("mol_to_cx_smiles", mol_to_cx_smiles(&mol));
        let mut writer = SmilesWriter::from_file_obj(Vec::new()).unwrap();
        assert_refused("SmilesWriter::write", writer.write(&mol));
    }
}

#[test]
fn transforms_refuse() {
    for mol in malformed("CCOc1ccccc1.N") {
        assert_refused("delete_substructs", delete_substructs(&mol, &smiles("O")));
        assert_refused(
            "replace_substructs",
            replace_substructs(&mol, &smiles("O"), &smiles("S")),
        );
        assert_refused(
            "replace_sidechains",
            replace_sidechains(&mol, &smiles("c1ccccc1")),
        );
        assert_refused("murcko_decompose", murcko_decompose(&mol));
        assert_refused("fragment_on_bonds", fragment_on_bonds(&mol, &[1]));
        assert_refused("fragment_on_brics_bonds", fragment_on_brics_bonds(&mol));
        assert_refused("molzip_single", molzip_single(&mol));
    }
}

#[test]
fn reactions_refuse() {
    let mut rxn = reaction_from_smarts("[C:1]O>>[C:1]").unwrap();
    for mol in malformed("CCO") {
        let before = mol_to_cx_smarts(&mol).unwrap();
        let mut reactant = mol.clone();
        assert_refused(
            "run_reactant_in_place",
            rxn.run_reactant_in_place(&mut reactant),
        );
        assert_eq!(mol_to_cx_smarts(&reactant).unwrap(), before);
        assert_refused(
            "reduce_product_to_side_chains",
            reduce_product_to_side_chains(&mol),
        );
        assert_refused("reaction_from_molecule", reaction_from_molecule(&mol));
        let mut templated = ChemicalReaction::default();
        templated.add_reactant_template(&mol);
        templated.add_product_template(&smiles("C"));
        assert_refused("reaction_to_smiles", reaction_to_smiles(&templated));
        let merge = SanitizeRxnParams {
            sanitize_ops: rdkit_cxx::chem::rd_chem_reactions::SanitizeFlags::SANITIZE_MERGEHS,
            ..Default::default()
        };
        assert_refused(
            "sanitize_rxn",
            sanitize_rxn_with_params(&mut templated, &merge),
        );
    }
    let mut reactant = smiles("CCO");
    assert!(rxn.run_reactant_in_place(&mut reactant).unwrap());
}

#[test]
fn standardization_refuses() {
    for mol in malformed("CC[O-].[Na+]") {
        let before = mol_to_cx_smarts(&mol).unwrap();
        assert_refused(
            "cleanup",
            rdkit_cxx::chem::mol_standardize::rd_mol_standardize::cleanup(&mol),
        );
        assert_refused("fragment_parent", fragment_parent(&mol));
        let mut copy = mol.clone();
        assert_refused("cleanup_in_place", cleanup_in_place(&mut copy));
        assert_refused(
            "FragmentRemover::remove_in_place",
            FragmentRemover::default().remove_in_place(&mut copy),
        );
        assert_refused(
            "LargestFragmentChooser::choose_in_place",
            LargestFragmentChooser::default().choose_in_place(&mut copy),
        );
        assert_refused(
            "MetalDisconnector::disconnect_in_place",
            MetalDisconnector::default().disconnect_in_place(&mut copy),
        );
        assert_refused(
            "Normalizer::normalize",
            Normalizer::default().normalize(&mol),
        );
        assert_eq!(mol_to_cx_smarts(&copy).unwrap(), before);
        assert_refused(
            "TautomerEnumerator::enumerate",
            TautomerEnumerator::default().enumerate(&mol),
        );
        assert_refused(
            "TautomerEnumerator::canonicalize",
            TautomerEnumerator::default().canonicalize(&mol),
        );
        assert_refused(
            "FragmentValidation::validate",
            FragmentValidation.validate(&mol),
        );
    }
    // CXSMILES can write the cycle itself.
    let cyclic = "CC.O |Sg:n:2::ht,Sg:n:0::ht,Sg:n:0::ht,SgH:1:2,2:1|";
    assert_refused("standardize_smiles", standardize_smiles(cyclic));
    assert_refused("validate_smiles", validate_smiles(cyclic));
    assert_refused("mol_to_smiles", mol_to_smiles(&smiles(cyclic)));
}

#[test]
fn other_modules_refuse() {
    for mol in malformed("CCOc1ccccc1.CCN") {
        assert_refused("mol_hash", mol_hash(&mol, HashFunction::MurckoScaffold));
        assert_refused("fragment_mol", fragment_mol(&mol));
        assert_refused("find_mces", find_mces(&mol, &smiles("CCO")));
        assert_refused(
            "generate_mol_extended_reduced_graph",
            generate_mol_extended_reduced_graph(&mol),
        );
        assert_refused("bcut_2d", bcut_2d(&mol));
        assert_refused(
            "create_scaffold_network",
            create_scaffold_network(
                std::slice::from_ref(&mol),
                &ScaffoldNetworkParams::default(),
            ),
        );
        assert_refused("create_extended_query_mol", create_extended_query_mol(&mol));
        assert_refused("TautomerQuery::new", TautomerQuery::new(&mol));
        assert_refused(
            "StereoisomerEnumerator::new",
            StereoisomerEnumerator::new(&mol),
        );
        assert_refused(
            "r_group_decompose",
            r_group_decompose(&[smiles("c1ccccc1[*:1]")], std::slice::from_ref(&mol)),
        );
        assert_refused(
            "r_group_decompose core",
            r_group_decompose(std::slice::from_ref(&mol), &[smiles("CCO")]),
        );
        let mut holder = CachedSmilesMolHolder::default();
        assert_refused("CachedSmilesMolHolder::add_mol", holder.add_mol(&mol));
        let mut copy = mol.clone();
        assert_refused("embed_molecule", embed_molecule(&mut copy));
        assert_refused("uff_optimize_molecule", uff_optimize_molecule(&mut copy));
        assert_refused("deprotect_in_place", deprotect_in_place(&mut copy));
        assert_refused("determine_connectivity", determine_connectivity(&mut copy));
    }
}

#[test]
fn depiction_with_r_groups_refuses() {
    let mut reference = smiles("c1ccccc1*");
    compute_2d_coords(&mut reference).unwrap();
    for bad in [
        with_parent_cycle(reference.clone()),
        with_bad_end_points(reference.clone()),
    ] {
        let params = ConstrainedDepictionParams {
            allow_r_groups: true,
            ..Default::default()
        };
        let mut mol = smiles("c1ccccc1C");
        assert_refused(
            "generate_depiction_matching_2d_structure",
            generate_depiction_matching_2d_structure_with_params(&mut mol, &bad, &params),
        );
    }
}
