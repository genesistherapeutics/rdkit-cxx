// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Runs the wrappers that lock RDKit's global state (CONVENTIONS.md U4)
//! from several threads at once: every one that writes the state a lock
//! guards, and at least one that only reads it. CI also runs this binary
//! under Valgrind's DRD, which sees RDKit's accesses to its globals even in
//! conda's uninstrumented build, and fails on any the locks leave unordered.
//!
//! DRD reports two accesses only when no chain of locks orders them, and
//! Valgrind runs one thread at a time: a thread that took a lock after an
//! unlocked access, and another that took the same lock before its own,
//! would hide the race. So the test runs in phases. Each starts four
//! threads at once that make one kind of call each: the wrapper under test
//! on two, and on the other two a plain wrapper that holds the same lock
//! (the random SMILES writer, a substructure match, a settings reader or
//! writer), whose accesses a missing lock leaves unordered with the first's.
//! The inputs are built on the main thread before the threads start.
//!
//! It is a plain `main` (`harness = false` in Cargo.toml), so libtest's
//! threads add nothing to DRD's report. To run the check as CI does, with
//! Valgrind installed:
//!
//! ```sh
//! CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUNNER="valgrind --tool=drd --suppressions=tests/drd.supp \
//!   --log-file=target/drd.log" \
//!   cargo test --test threads
//! python3 tests/check_drd_log.py target/drd.log
//! ```

use std::collections::BTreeSet;
use std::path::PathBuf;
use std::thread;

use rdkit_cxx::chem::draw::rd_mol_draw_2d::{MolDraw2DSVG, MolDrawOptions};
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::{
    FragmentRemover, MetalDisconnector, MolVSValidation, Normalizer, Pipeline, TautomerEnumerator,
    Uncharger, cleanup_in_place_mols,
};
use rdkit_cxx::chem::rd_abbreviations::{
    AbbreviationDefinition, condense_mol_abbreviations, get_default_abbreviations,
};
use rdkit_cxx::chem::rd_chem_reactions::{
    ChemicalReaction, Compute2DCoordsForReactionParams, EnumerateLibrary,
    compute_2d_coords_for_reaction_with_params, reaction_from_smarts, reaction_to_rxn_block,
    reaction_to_smiles_with_params,
};
use rdkit_cxx::chem::rd_coord_gen::{
    CoordGenParams, add_coords, add_coords_with_params, set_default_template_file_dir,
};
use rdkit_cxx::chem::rd_depictor::{
    Compute2DCoordsParams, add_ring_system_templates, compute_2d_coords,
    compute_2d_coords_mimic_distmat, compute_2d_coords_with_params,
    generate_depiction_matching_2d_structure, generate_depiction_matching_3d_structure,
    is_coord_gen_support_available, load_default_ring_system_templates, set_prefer_coord_gen,
    set_ring_system_templates,
};
use rdkit_cxx::chem::rd_deprotect::deprotect;
use rdkit_cxx::chem::rd_determine_bonds::determine_bonds;
use rdkit_cxx::chem::rd_dist_geom::{
    EmbedParameters, GetMoleculeBoundsMatrixParams, embed_molecule,
    embed_multiple_confs_with_params, get_experimental_torsions,
    get_molecule_bounds_matrix_with_params,
};
use rdkit_cxx::chem::rd_enumerate_stereoisomers::{
    StereoEnumerationOptions, StereoisomerEnumerator,
};
use rdkit_cxx::chem::rd_fingerprint_generator::{
    FingerprintGenerator64, GetMorganGeneratorParams, get_morgan_feature_atom_inv_gen,
    get_morgan_generator_with_params,
};
use rdkit_cxx::chem::rd_fmcs::find_mcs;
use rdkit_cxx::chem::rd_force_field_helpers::uff_optimize_molecule;
use rdkit_cxx::chem::rd_generalized_substruct::{
    self, ExtendedQueryMol, create_extended_query_mol,
};
use rdkit_cxx::chem::rd_mif::HBond;
use rdkit_cxx::chem::rd_mmpa::fragment_mol;
use rdkit_cxx::chem::rd_mol_align::{
    RandomTransformParams, align_mol, get_best_rms, random_transform, random_transform_with_params,
};
use rdkit_cxx::chem::rd_mol_chemical_features::{
    MolChemicalFeatureFactory, build_feature_factory_from_string,
};
use rdkit_cxx::chem::rd_mol_descriptors::{
    DoubleCubicLatticeVolume, Properties, PropertyFunctor, calc_crippen_descriptors,
    get_maccs_keys_fingerprint, get_usrcat,
};
use rdkit_cxx::chem::rd_mol_enumerator::enumerate;
use rdkit_cxx::chem::rd_mol_hash::{HashFunction, mol_hash};
use rdkit_cxx::chem::rd_mol_processing::{
    GetFingerprintsForMolsInFileParams, SupplierOptions,
    get_fingerprints_for_mols_in_file_with_params,
};
use rdkit_cxx::chem::rd_r_group_decomposition::RGroupDecomposition;
use rdkit_cxx::chem::rd_rascal_mces::{find_mces, rascal_cluster};
use rdkit_cxx::chem::rd_reduced_graphs::get_erg_fingerprint;
use rdkit_cxx::chem::rd_shape_align;
use rdkit_cxx::chem::rd_substruct_library::{MolHolder, SubstructLibrary};
use rdkit_cxx::chem::rd_synthon_space_search::SynthonSpace;
use rdkit_cxx::chem::rd_tautomer_query::TautomerQuery;
use rdkit_cxx::chem::rdfiltercatalog::{FilterCatalog, FilterCatalogEntry, SmartsMatcher};
use rdkit_cxx::chem::rdfragcatalog::{FragCatGenerator, FragCatParams, FragCatalog};
use rdkit_cxx::chem::rdinchi::{inchi_to_mol, mol_to_inchi};
use rdkit_cxx::chem::scaffolds::rd_scaffold_network::{
    ScaffoldNetworkParams, create_scaffold_network,
};
use rdkit_cxx::chem::{
    Conformer, Mol, MolFragmentToSmilesParams, MolFromSCSRParams, MolToCXSmilesParams,
    MolToRandomSmilesVectParams, MolzipParams, PNGMetadataParams, ResonanceMolSupplier,
    SCSRBaseHbondOptions, SDWriter, SmilesWriteParams, assign_cip_labels, delete_substructs,
    fragment_on_brics_bonds, mol_fragment_to_cx_smiles_with_params,
    mol_fragment_to_smiles_with_params, mol_from_scsr_block_with_params, mol_from_smarts,
    mol_from_smiles, mol_from_xyz_block, mol_metadata_to_png_string_with_params,
    mol_to_cx_smiles_with_params, mol_to_mol_block,
    mol_to_random_smiles_vect_with_params, mol_to_smiles, mol_to_smiles_with_params,
    molzip_fragments_with_params, pattern_fingerprint,
};
use rdkit_cxx::distance_geometry::embed_bounds_matrix;
use rdkit_cxx::geometry::{
    Point3D, UniformGrid3D, compute_grid_centroid, find_grid_terminal_points,
};
use rdkit_cxx::ml::cluster::clustering::murtagh_cluster;
use rdkit_cxx::sim_div_filters::rd_sim_div_pickers::{ClusterMethod, HierarchicalClusterPicker};

const ASPIRIN: &str = "CC(=O)Oc1ccccc1C(=O)O";

/// Propionic acid and ammonia to propanamide.
const AMIDE: &str = "[CH3:4][CH2:5][C:1](=[O:2])O.[N:3]>>[CH3:4][CH2:5][C:1](=[O:2])[N:3]";

/// A 1x1 grayscale PNG for the metadata writer.
const PNG: &[u8] = &[
    137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73, 72, 68, 82, 0, 0, 0, 1, 0, 0, 0, 1, 8, 0, 0,
    0, 0, 58, 126, 155, 85, 0, 0, 0, 10, 73, 68, 65, 84, 120, 156, 99, 96, 0, 0, 0, 2, 0, 1, 72,
    175, 164, 113, 0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130,
];

/// Ethanol's atoms and 3D coordinates, without bonds, for DetermineBonds.
const ETHANOL_XYZ: &str = "9
ethanol
C -0.9640 0.0151 0.0410
C 0.5041 -0.3085 -0.2036
O 1.2237 0.3416 0.7792
H -1.3563 0.7407 -0.7043
H -1.5833 -0.9142 -0.0029
H -1.0287 0.5174 1.0296
H 0.7725 -0.0064 -1.2300
H 0.5980 -1.4177 -0.1686
H 1.8339 1.0320 0.4595
";

/// An RNA dimer whose base templates SCSR's Auto option matches to find
/// their hydrogen-bond sites.
const RIBOSE_SCSR: &str = include_str!("data/RiboseFullname.mol");

fn main() {
    // RDKit fills some static data on first use (the periodic table, the
    // SMILES parser's and writer's tables), ordered by C++'s thread-safe
    // initialization, whose fast path DRD cannot see. Filling it before any
    // thread starts keeps those reads out of the report. Globals that only
    // the locks order (the BRICS patterns) are left for the threads to fill.
    mol_to_smiles(&mol(ASPIRIN)).unwrap();
    // And the drawers' (boost::format's, and the text drawers').
    MolDraw2DSVG::new(100, 100)
        .unwrap()
        .draw_molecule(&depicted(ASPIRIN))
        .unwrap();
    std::fs::write(
        templates(),
        "C1CC2CCC1C2 |(1.55254,0.75,;1.55254,-0.75,;0.125954,-1.21353,;-1.30063,-0.75,;\
         -1.30063,0.75,;0.125954,1.21353,;-0.755724,0,)|\n",
    )
    .unwrap();
    std::fs::write(
        groups(),
        "-C(=O)O\t*-C(=O)[O;H1,-]\tcarboxylic acids\n-O\t*-[O;H1,-]\talcohols\n\n=O\t*=O\taldehydes\n",
    )
    .unwrap();
    std::fs::write(synthons(), SYNTHONS).unwrap();

    // RDKit's random number generator: each module's writers against the
    // SMILES writer.
    phase(varied_random_smiles(), varied_random_smiles());
    phase(random_cx_smiles(), random_smiles());
    phase(random_fragment_smiles(), random_smiles());
    phase(seeded_random_smiles(), random_smiles());
    phase(random_png_metadata(), random_smiles());
    phase(flipped_coords(), random_smiles());
    phase(mimicked_distances(), random_smiles());
    phase(matched_3d_depiction(), random_smiles());
    phase(random_reaction_smiles(), random_smiles());
    phase(flipped_reaction_coords(), random_smiles());
    phase(unseeded_embedding(), random_smiles());
    phase(bounds_matrix_embedding(), random_smiles());
    phase(random_transforms(), random_smiles());
    phase(embedded_stereoisomers(), random_smiles());
    // Comic-mode drawing draws its lines with std::rand, under the same
    // lock.
    phase(comic_drawing(), comic_drawing());

    // 2D depiction's settings: each writer against a reader, from every
    // module that reads them.
    phase(prefer_coord_gen(), mol_block());
    phase(ring_system_templates(), ring_template_coords());
    phase(added_ring_system_templates(), zipped_coords());
    phase(coord_gen_template_dir(), coord_gen_coords());
    phase(bond_length_coords(), sd_text());
    phase(default_ring_system_templates(), ring_template_coords());
    phase(prefer_coord_gen(), png_mol_block());
    phase(prefer_coord_gen(), rxn_block());
    phase(prefer_coord_gen(), drawn_without_coords());

    // ControlCHandler's statics: a substructure match from each module
    // against a plain one.
    phase(has_substruct_match(), substruct_matches());
    phase(resonance_matches(), substruct_matches());
    phase(deleted_substructs(), substruct_matches());
    phase(brics_fragments(), substruct_matches());
    phase(best_rms(), substruct_matches());
    phase(aligned(), substruct_matches());
    phase(crippen(), substruct_matches());
    phase(maccs_keys(), substruct_matches());
    phase(computed_properties(), substruct_matches());
    phase(reaction_products(), substruct_matches());
    phase(enumerated_library(), substruct_matches());
    phase(uff_optimized(), substruct_matches());
    phase(deprotected(), substruct_matches());
    phase(erg_fingerprint(), substruct_matches());
    phase(pattern_fp(), substruct_matches());
    phase(query_atoms(), substruct_matches());
    phase(hashed(), substruct_matches());
    phase(mmpa_fragments(), substruct_matches());
    phase(usrcat(), substruct_matches());
    phase(matched_2d_depiction(), substruct_matches());
    phase(templated_coord_gen(), substruct_matches());
    phase(seeded_embedding(), substruct_matches());
    phase(experimental_torsions(), substruct_matches());
    phase(feature_morgan_fp(), substruct_matches());
    phase(feature_fps_in_file(), substruct_matches());
    phase(fps_in_file(), substruct_matches());
    phase(scsr_hbonds(), substruct_matches());
    phase(cip_labels(), substruct_matches());
    phase(bonds_from_xyz(), substruct_matches());
    phase(cleaned_up_mols(), substruct_matches());
    phase(uncharged(), substruct_matches());
    phase(fragments_removed(), substruct_matches());
    phase(normalized(), substruct_matches());
    phase(metals_disconnected(), substruct_matches());
    phase(canonical_tautomer(), substruct_matches());
    phase(validated(), substruct_matches());
    phase(pipeline_run(), substruct_matches());
    phase(filter_catalog_match(), substruct_matches());
    phase(fragments_cataloged(), substruct_matches());
    phase(maximum_common_substructure(), substruct_matches());
    phase(maximum_common_edge_substructure(), substruct_matches());
    phase(rascal_clusters(), substruct_matches());
    phase(tautomer_query_matches(), substruct_matches());
    phase(extended_query_matches(), substruct_matches());
    phase(substruct_library_matches(), substruct_matches());
    phase(rgroups_decomposed(), substruct_matches());
    phase(scaffold_network(), substruct_matches());
    phase(pubchem_shape_alignment(), substruct_matches());
    phase(mif_h_bonds(), substruct_matches());
    phase(mol_features(), substruct_matches());
    phase(condensed_abbreviations(), substruct_matches());
    phase(enumerated_link_nodes(), substruct_matches());
    phase(embedded_stereoisomers(), substruct_matches());
    phase(synthon_substructure_search(), substruct_matches());
    phase(synthon_fingerprint_search(), substruct_matches());
    phase(synthon_rascal_search(), substruct_matches());
    phase(inchi_written(), substruct_matches());
    phase(inchi_read(), substruct_matches());

    // The rest: the grid centroids' cache, the property registry, DCLV's
    // surface cache, PubChem shape's globals, and Murtagh's statics.
    phase(grid_centroids(), grid_centroids());
    phase(registered_property(), available_properties());
    phase(dclv_surfaces(), dclv_surfaces());
    phase(pubchem_shape_alignment(), pubchem_shape_inputs());
    phase(murtagh_clusters(), hierarchical_picks());
    phase(cdxml_blocks(), cdxml_blocks());

    std::fs::remove_file(templates()).unwrap();
    std::fs::remove_file(groups()).unwrap();
    std::fs::remove_file(synthons()).unwrap();
}

/// One thread's call in a phase: `run` on the input `make` builds, on the
/// main thread, for thread `i`.
struct Call<T> {
    make: fn(u32) -> T,
    run: fn(T),
}

/// Runs `even` on threads 0 and 2 and `odd` on threads 1 and 3, at once.
fn phase<A: Send + 'static, B: Send + 'static>(even: Call<A>, odd: Call<B>) {
    let (a0, b1, a2, b3) = ((even.make)(0), (odd.make)(1), (even.make)(2), (odd.make)(3));
    let handles = [
        thread::spawn(move || (even.run)(a0)),
        thread::spawn(move || (odd.run)(b1)),
        thread::spawn(move || (even.run)(a2)),
        thread::spawn(move || (odd.run)(b3)),
    ];
    // Joining each thread orders its accesses before the next phase's, and
    // before the destructors that run at exit, for DRD.
    for handle in handles {
        handle.join().unwrap();
    }
}

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn random() -> SmilesWriteParams {
    SmilesWriteParams {
        do_random: true,
        ..Default::default()
    }
}

fn amide() -> ChemicalReaction {
    reaction_from_smarts(AMIDE).unwrap()
}

fn amide_reactants() -> Vec<Mol> {
    vec![mol("CCC(=O)O"), mol("N")]
}

/// Functional groups for the fragment catalog, written by `main`.
fn groups() -> PathBuf {
    std::env::temp_dir().join(format!(
        "rdkit-cxx-threads-{}-groups.txt",
        std::process::id()
    ))
}

/// Five of RDKit's CDK2 ligands, which Rascal clusters together.
const CDK2: [&str; 5] = [
    "CC(C)C(=O)COc1nc(N)nc2[nH]cnc12",
    "Nc1nc(OCC2CCCO2)c2nc[nH]c2n1",
    "Nc1nc(OCC2CCC(=O)N2)c2nc[nH]c2n1",
    "Nc1nc(OCC2CCCCC2)c2nc[nH]c2n1",
    "Nc1nc(OCC2CC=CCC2)c2nc[nH]c2n1",
];

/// Ring-system templates for the depictor, written by `main`.
fn templates() -> PathBuf {
    std::env::temp_dir().join(format!(
        "rdkit-cxx-threads-{}-templates.cxsmi",
        std::process::id()
    ))
}

/// A molecule with an embedded conformer.
fn embedded(smiles: &str) -> Mol {
    let mut mol = mol(smiles);
    let mut seeded = EmbedParameters {
        random_seed: 42,
        ..Default::default()
    };
    embed_multiple_confs_with_params(&mut mol, 1, &mut seeded).unwrap();
    mol
}

/// A molecule with 2D coordinates.
fn depicted(smiles: &str) -> Mol {
    let mut mol = mol(smiles);
    compute_2d_coords(&mut mol).unwrap();
    mol
}

/// Propanol in 3D, from Python's EmbedMolecule(AddHs(CCCO), randomSeed=42)
/// without the hydrogens.
fn propanol_3d() -> Mol {
    let mut mol = mol("CCCO");
    let mut conf = Conformer::new(4);
    let positions = [
        (-1.2784, -0.0795, -0.4157),
        (-0.1858, 0.7166, 0.2689),
        (0.9073, -0.314, 0.5667),
        (1.2802, -0.8242, -0.6921),
    ];
    for (atom, (x, y, z)) in (0..).zip(positions) {
        conf.set_atom_position(atom, Point3D::new(x, y, z)).unwrap();
    }
    mol.add_conformer(&conf).unwrap();
    mol
}

// RDKit's random number generator. Drawing from it, and seeding it, write
// it.

fn random_smiles() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            mol_to_smiles_with_params(&mol, &random()).unwrap();
        },
    }
}

fn varied_random_smiles() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            let canonical = mol_to_smiles(&mol).unwrap();
            let mut written = BTreeSet::new();
            for _ in 0..50 {
                let smiles = mol_to_smiles_with_params(&mol, &random()).unwrap();
                let parsed = mol_from_smiles(&smiles).unwrap();
                assert_eq!(mol_to_smiles(&parsed).unwrap(), canonical);
                written.insert(smiles);
            }
            assert!(written.len() > 1);
        },
    }
}

fn random_cx_smiles() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            let cx = MolToCXSmilesParams {
                params: random(),
                ..Default::default()
            };
            mol_to_cx_smiles_with_params(&mol, &cx).unwrap();
        },
    }
}

fn random_fragment_smiles() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            let fragment = MolFragmentToSmilesParams {
                params: random(),
                ..Default::default()
            };
            mol_fragment_to_smiles_with_params(&mol, &[0, 1, 2], &fragment).unwrap();
            mol_fragment_to_cx_smiles_with_params(&mol, &[0, 1, 2], &fragment).unwrap();
        },
    }
}

fn seeded_random_smiles() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            let seeded = MolToRandomSmilesVectParams {
                random_seed: 7,
                ..Default::default()
            };
            let smiles = mol_to_random_smiles_vect_with_params(&mol, 5, &seeded).unwrap();
            assert_eq!(smiles.len(), 5);
        },
    }
}

fn random_png_metadata() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            let png = PNGMetadataParams {
                smiles_write_params: random(),
                ..Default::default()
            };
            let written = mol_metadata_to_png_string_with_params(&mol, PNG, &png).unwrap();
            assert!(written.len() > PNG.len());
        },
    }
}

// Random bond flips, with rotatable bonds to flip.
fn flipped_coords() -> Call<Mol> {
    Call {
        make: |_| mol("c1ccccc1CCC(=O)O"),
        run: |mut mol| {
            let sampled = Compute2DCoordsParams {
                n_flips_per_sample: 2,
                n_sample: 5,
                ..Default::default()
            };
            compute_2d_coords_with_params(&mut mol, &sampled).unwrap();
        },
    }
}

fn mimicked_distances() -> Call<Mol> {
    Call {
        make: |_| mol("c1ccccc1CCC(=O)O"),
        run: |mut mol| {
            let n = mol.get_num_atoms() as usize;
            compute_2d_coords_mimic_distmat(&mut mol, &vec![1.5; n * (n - 1) / 2]).unwrap();
        },
    }
}

// Matching a 3D structure mimics its distances.
fn matched_3d_depiction() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol("CCCO"), propanol_3d()),
        run: |(mut mol, reference)| {
            generate_depiction_matching_3d_structure(&mut mol, &reference).unwrap();
            assert_eq!(mol.get_num_conformers(), 1);
        },
    }
}

fn random_reaction_smiles() -> Call<ChemicalReaction> {
    Call {
        make: |_| amide(),
        run: |reaction| {
            reaction_to_smiles_with_params(&reaction, &random()).unwrap();
        },
    }
}

fn flipped_reaction_coords() -> Call<ChemicalReaction> {
    Call {
        make: |_| amide(),
        run: |mut reaction| {
            let flips = Compute2DCoordsForReactionParams {
                n_flips_per_sample: 2,
                n_sample: 3,
                ..Default::default()
            };
            compute_2d_coords_for_reaction_with_params(&mut reaction, &flips).unwrap();
        },
    }
}

fn unseeded_embedding() -> Call<Mol> {
    Call {
        make: |_| mol("CCCO"),
        run: |mut mol| assert_eq!(embed_molecule(&mut mol).unwrap(), 0),
    }
}

fn bounds_matrix_embedding() -> Call<Vec<Vec<f64>>> {
    Call {
        make: |_| {
            let unsmoothed = GetMoleculeBoundsMatrixParams {
                do_triangle_smoothing: false,
                ..Default::default()
            };
            get_molecule_bounds_matrix_with_params(&mol("CCCC"), &unsmoothed).unwrap()
        },
        run: |bounds| assert_eq!(embed_bounds_matrix(&bounds).unwrap().len(), 4),
    }
}

fn random_transforms() -> Call<Mol> {
    Call {
        make: |_| propanol_3d(),
        run: |mut mol| {
            random_transform(&mut mol).unwrap();
            let seeded = RandomTransformParams {
                seed: 5,
                ..Default::default()
            };
            random_transform_with_params(&mut mol, &seeded).unwrap();
        },
    }
}

// 2D depiction's settings. Computing 2D coordinates reads them, as do the
// writers that compute them for a molecule without a conformer.

fn prefer_coord_gen() -> Call<()> {
    Call {
        make: |_| (),
        run: |()| {
            set_prefer_coord_gen(true);
            set_prefer_coord_gen(false);
        },
    }
}

fn mol_block() -> Call<Mol> {
    Call {
        make: |_| mol("F[C@H](Cl)Br"),
        run: |mol| assert!(mol_to_mol_block(&mol).unwrap().contains("V2000")),
    }
}

fn ring_system_templates() -> Call<PathBuf> {
    Call {
        make: |_| templates(),
        run: |path| set_ring_system_templates(&path).unwrap(),
    }
}

// Norbornane, which the templates match.
fn ring_template_coords() -> Call<Mol> {
    Call {
        make: |_| mol("C1CC2CCC1C2"),
        run: |mut mol| {
            let templated = Compute2DCoordsParams {
                use_ring_templates: true,
                ..Default::default()
            };
            compute_2d_coords_with_params(&mut mol, &templated).unwrap();
        },
    }
}

fn added_ring_system_templates() -> Call<PathBuf> {
    Call {
        make: |_| templates(),
        run: |path| add_ring_system_templates(&path).unwrap(),
    }
}

fn zipped_coords() -> Call<Vec<Mol>> {
    Call {
        make: |_| vec![mol("c1ccccc1[*:1]"), mol("O[*:1]")],
        run: |mut mols| {
            let coords = MolzipParams {
                generate_coordinates: true,
                ..Default::default()
            };
            let zipped = molzip_fragments_with_params(&mut mols, &coords)
                .unwrap()
                .unwrap();
            assert_eq!(zipped.get_num_conformers(), 1);
        },
    }
}

fn coord_gen_template_dir() -> Call<()> {
    Call {
        make: |_| (),
        run: |()| set_default_template_file_dir("").unwrap(),
    }
}

fn coord_gen_coords() -> Call<Mol> {
    Call {
        make: |_| mol("c1ccccc1CCC(=O)O"),
        run: |mut mol| {
            if is_coord_gen_support_available() {
                add_coords(&mut mol).unwrap();
            }
        },
    }
}

fn bond_length_coords() -> Call<Mol> {
    Call {
        make: |_| mol("c1ccccc1CCC(=O)O"),
        run: |mut mol| {
            let longer = Compute2DCoordsParams {
                bond_length: 1.2,
                ..Default::default()
            };
            compute_2d_coords_with_params(&mut mol, &longer).unwrap();
        },
    }
}

fn sd_text() -> Call<Mol> {
    Call {
        make: |_| mol("F[C@H](Cl)Br"),
        run: |mol| assert!(SDWriter::get_text(&mol).unwrap().contains("$$$$")),
    }
}

fn default_ring_system_templates() -> Call<()> {
    Call {
        make: |_| (),
        run: |()| load_default_ring_system_templates(),
    }
}

fn png_mol_block() -> Call<Mol> {
    Call {
        make: |_| mol("F[C@H](Cl)Br"),
        run: |mol| {
            let with_mol = PNGMetadataParams {
                include_mol: true,
                ..Default::default()
            };
            mol_metadata_to_png_string_with_params(&mol, PNG, &with_mol).unwrap();
        },
    }
}

fn rxn_block() -> Call<ChemicalReaction> {
    Call {
        make: |_| amide(),
        run: |reaction| {
            assert!(
                reaction_to_rxn_block(&reaction)
                    .unwrap()
                    .starts_with("$RXN")
            )
        },
    }
}

// Drawing a molecule without conformers depicts it first.
fn drawn_without_coords() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            let mut d = MolDraw2DSVG::new(200, 200).unwrap();
            d.draw_molecule(&mol).unwrap();
            assert!(d.get_drawing_text().unwrap().contains("<path"));
        },
    }
}

// Comic mode wiggles each line with std::rand (handdrawnLine in
// MolDraw2DDetails.cpp).
fn comic_drawing() -> Call<Mol> {
    Call {
        make: |_| depicted(ASPIRIN),
        run: |mol| {
            let mut d = MolDraw2DSVG::new(200, 200).unwrap();
            let comic = MolDrawOptions {
                comic_mode: true,
                ..Default::default()
            };
            d.set_draw_options(&comic).unwrap();
            d.draw_molecule(&mol).unwrap();
        },
    }
}

// ControlCHandler's statics, which every substructure match writes.

fn substruct_matches() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol(ASPIRIN), mol_from_smarts("C(=O)[OH]").unwrap()),
        run: |(mol, query)| assert_eq!(mol.get_substruct_matches(&query).unwrap().len(), 1),
    }
}

fn has_substruct_match() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol(ASPIRIN), mol_from_smarts("c1ccccc1").unwrap()),
        run: |(mol, query)| assert!(mol.has_substruct_match(&query).unwrap()),
    }
}

fn resonance_matches() -> Call<(ResonanceMolSupplier, Mol)> {
    Call {
        make: |_| {
            let acetate = mol("CC(=O)[O-]");
            (
                ResonanceMolSupplier::new(&acetate).unwrap(),
                mol_from_smarts("C=O").unwrap(),
            )
        },
        run: |(mut suppl, query)| assert!(!suppl.get_substruct_matches(&query).unwrap().is_empty()),
    }
}

fn deleted_substructs() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol(ASPIRIN), mol_from_smarts("C(=O)[OH]").unwrap()),
        run: |(mol, query)| {
            let deleted = delete_substructs(&mol, &query).unwrap();
            assert_eq!(deleted.get_num_atoms(), 10);
        },
    }
}

// The first BRICS call also fills RDKit's pattern globals.
fn brics_fragments() -> Call<Mol> {
    Call {
        make: |_| mol("CCCOCCC(=O)c1ccccc1"),
        run: |mol| assert_eq!(fragment_on_brics_bonds(&mol).unwrap().get_num_atoms(), 20),
    }
}

// Best-RMS alignment also symmetrizes terminal atoms.
fn best_rms() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (embedded("CCCO"), propanol_3d()),
        run: |(mut probe, reference)| assert!(get_best_rms(&mut probe, &reference).unwrap() >= 0.0),
    }
}

fn aligned() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (embedded("CCCO"), propanol_3d()),
        run: |(mut probe, reference)| assert!(align_mol(&mut probe, &reference).unwrap() >= 0.0),
    }
}

fn crippen() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| assert!(calc_crippen_descriptors(&mol).unwrap().0 > 0.0),
    }
}

fn maccs_keys() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            get_maccs_keys_fingerprint(&mol).unwrap();
        },
    }
}

// RDKit's own property functors match substructures.
fn computed_properties() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            let properties = Properties::new(&["NumHBD"]).unwrap();
            assert_eq!(properties.compute_properties(&mol).unwrap(), vec![1.0]);
        },
    }
}

fn reaction_products() -> Call<(ChemicalReaction, Vec<Mol>)> {
    Call {
        make: |_| (amide(), amide_reactants()),
        run: |(mut reaction, reactants)| {
            assert_eq!(reaction.run_reactants(&reactants).unwrap().len(), 1)
        },
    }
}

fn enumerated_library() -> Call<(ChemicalReaction, Vec<Vec<Mol>>)> {
    Call {
        make: |_| (amide(), vec![vec![mol("CCC(=O)O")], vec![mol("N")]]),
        run: |(reaction, blocks)| {
            let mut library = EnumerateLibrary::new(&reaction, &blocks).unwrap();
            assert_eq!(library.next_smiles().unwrap().len(), 1);
        },
    }
}

// Typing atoms for UFF matches substructures.
fn uff_optimized() -> Call<Mol> {
    Call {
        make: |_| embedded("CCO"),
        run: |mut mol| {
            uff_optimize_molecule(&mut mol).unwrap();
        },
    }
}

fn deprotected() -> Call<Mol> {
    Call {
        make: |_| mol("CC(C)(C)OC(=O)NCC(=O)O"),
        run: |boc| {
            assert_eq!(
                mol_to_smiles(&deprotect(&boc).unwrap()).unwrap(),
                "NCC(=O)O"
            )
        },
    }
}

fn erg_fingerprint() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            get_erg_fingerprint(&mol).unwrap();
        },
    }
}

fn pattern_fp() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            pattern_fingerprint(&mol).unwrap();
        },
    }
}

fn query_atoms() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol(ASPIRIN), mol_from_smarts("[OX2H]").unwrap()),
        run: |(mol, hydroxyl)| {
            let query = hydroxyl.get_atom_with_idx(0).unwrap();
            assert_eq!(mol.get_atoms_matching_query(&query).unwrap().len(), 1);
        },
    }
}

fn hashed() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| {
            mol_hash(&mol, HashFunction::HetAtomTautomerv2).unwrap();
        },
    }
}

fn mmpa_fragments() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| assert!(!fragment_mol(&mol).unwrap().is_empty()),
    }
}

fn usrcat() -> Call<Mol> {
    Call {
        make: |_| embedded("CCCO"),
        run: |mol| {
            get_usrcat(&mol).unwrap();
        },
    }
}

fn matched_2d_depiction() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol("c1ccccc1CCO"), depicted("c1ccccc1")),
        run: |(mut mol, reference)| {
            assert_eq!(
                generate_depiction_matching_2d_structure(&mut mol, &reference)
                    .unwrap()
                    .len(),
                6
            );
        },
    }
}

fn templated_coord_gen() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol("c1ccccc1CCO"), depicted("c1ccccc1")),
        run: |(mut mol, template)| {
            if is_coord_gen_support_available() {
                let templated = CoordGenParams {
                    template_mol: Some(template),
                    ..Default::default()
                };
                add_coords_with_params(&mut mol, &templated).unwrap();
            }
        },
    }
}

// Seeded, the embedder leaves the generator alone; RMS pruning also
// symmetrizes terminal atoms.
fn seeded_embedding() -> Call<Mol> {
    Call {
        make: |_| mol("CCCO"),
        run: |mut mol| {
            let mut pruned = EmbedParameters {
                prune_rms_thresh: 0.1,
                random_seed: 42,
                ..Default::default()
            };
            let ids = embed_multiple_confs_with_params(&mut mol, 3, &mut pruned).unwrap();
            assert!(!ids.is_empty());
        },
    }
}

fn experimental_torsions() -> Call<Mol> {
    Call {
        make: |_| mol("CCCO"),
        run: |mol| {
            get_experimental_torsions(&mol).unwrap();
        },
    }
}

// Feature invariants match their feature patterns.
fn feature_morgan_fp() -> Call<(FingerprintGenerator64, Mol)> {
    Call {
        make: |_| {
            let features = GetMorganGeneratorParams {
                atom_invariants_generator: Some(get_morgan_feature_atom_inv_gen()),
                ..Default::default()
            };
            (
                get_morgan_generator_with_params(&features).unwrap(),
                mol(ASPIRIN),
            )
        },
        run: |(generator, mol)| {
            generator.get_fingerprint(&mol).unwrap();
        },
    }
}

/// A SMILES file for thread `i`, which the call removes.
fn smiles_file(i: u32) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("rdkit-cxx-threads-{}-{i}.smi", std::process::id()));
    std::fs::write(
        &path,
        format!("SMILES\tName\n{ASPIRIN}\ta\nc1ccccc1O\tb\nCCN\tc\n"),
    )
    .unwrap();
    path
}

fn fps_in_file_with(generator: &FingerprintGenerator64, path: PathBuf) {
    let params = GetFingerprintsForMolsInFileParams {
        generator: Some(generator),
        options: SupplierOptions {
            num_threads: 2,
            ..Default::default()
        },
    };
    let fps = get_fingerprints_for_mols_in_file_with_params(&path, &params).unwrap();
    assert_eq!(fps.len(), 3);
    std::fs::remove_file(&path).unwrap();
}

/// rdMolProcessing with feature invariants: RDKit's threads parse, and the
/// calling thread fingerprints under the Ctrl-C lock.
fn feature_fps_in_file() -> Call<(FingerprintGenerator64, PathBuf)> {
    Call {
        make: |i| {
            let features = GetMorganGeneratorParams {
                atom_invariants_generator: Some(get_morgan_feature_atom_inv_gen()),
                ..Default::default()
            };
            (
                get_morgan_generator_with_params(&features).unwrap(),
                smiles_file(i),
            )
        },
        run: |(generator, path)| fps_in_file_with(&generator, path),
    }
}

/// rdMolProcessing without a lock: RDKit's threads fingerprint the molecules
/// they parse, in the wrapper's write callback.
fn fps_in_file() -> Call<(FingerprintGenerator64, PathBuf)> {
    Call {
        make: |i| {
            (
                get_morgan_generator_with_params(&GetMorganGeneratorParams::default()).unwrap(),
                smiles_file(i),
            )
        },
        run: |(generator, path)| fps_in_file_with(&generator, path),
    }
}

fn scsr_hbonds() -> Call<()> {
    Call {
        make: |_| (),
        run: |()| {
            let auto = MolFromSCSRParams {
                scsr_base_hbond_options: SCSRBaseHbondOptions::Auto,
                ..Default::default()
            };
            mol_from_scsr_block_with_params(RIBOSE_SCSR, &auto).unwrap();
        },
    }
}

fn cip_labels() -> Call<Mol> {
    Call {
        make: |_| mol("F[C@H](Cl)Br"),
        run: |mut mol| assign_cip_labels(&mut mol).unwrap(),
    }
}

fn bonds_from_xyz() -> Call<Mol> {
    Call {
        make: |_| mol_from_xyz_block(ETHANOL_XYZ).unwrap(),
        run: |mut mol| {
            determine_bonds(&mut mol).unwrap();
            assert_eq!(mol.get_num_bonds(), 8);
        },
    }
}

// Standardization. The multi-molecule forms run on the wrapper's own
// threads, which take the lock in turn.
fn cleaned_up_mols() -> Call<Vec<Mol>> {
    Call {
        make: |_| vec![mol("CC(=O)[O-].[Na+]"), mol("CN(=O)=O"), mol("Oc1ccccn1")],
        run: |mut mols| cleanup_in_place_mols(&mut mols, 2).unwrap(),
    }
}

fn uncharged() -> Call<(Uncharger, Mol)> {
    Call {
        make: |_| (Uncharger::default(), mol("CC(=O)[O-]")),
        run: |(uncharger, mol)| {
            let neutral = uncharger.uncharge(&mol).unwrap();
            assert_eq!(mol_to_smiles(&neutral).unwrap(), "CC(=O)O");
        },
    }
}

fn fragments_removed() -> Call<(FragmentRemover, Mol)> {
    Call {
        make: |_| (FragmentRemover::default(), mol("CC(=O)O.Cl")),
        run: |(remover, mol)| {
            let kept = remover.remove(&mol).unwrap();
            assert_eq!(mol_to_smiles(&kept).unwrap(), "CC(=O)O");
        },
    }
}

fn normalized() -> Call<(Normalizer, Mol)> {
    Call {
        make: |_| (Normalizer::default(), mol("CN(=O)=O")),
        run: |(normalizer, mol)| {
            let nitro = normalizer.normalize(&mol).unwrap();
            assert_eq!(mol_to_smiles(&nitro).unwrap(), "C[N+](=O)[O-]");
        },
    }
}

fn metals_disconnected() -> Call<(MetalDisconnector, Mol)> {
    Call {
        make: |_| (MetalDisconnector::new(), mol("CCC(=O)O[Na]")),
        run: |(disconnector, mol)| {
            let ions = disconnector.disconnect(&mol).unwrap();
            assert_eq!(mol_to_smiles(&ions).unwrap(), "CCC(=O)[O-].[Na+]");
        },
    }
}

fn canonical_tautomer() -> Call<(TautomerEnumerator, Mol)> {
    Call {
        make: |_| (TautomerEnumerator::default(), mol("Oc1ccccn1")),
        run: |(enumerator, mol)| {
            let pyridone = enumerator.canonicalize(&mol).unwrap();
            assert_eq!(mol_to_smiles(&pyridone).unwrap(), "O=c1cccc[nH]1");
        },
    }
}

fn validated() -> Call<(MolVSValidation, Mol)> {
    Call {
        make: |_| (MolVSValidation::default(), mol("CC(=O)O.Cl")),
        run: |(validation, mol)| assert!(!validation.validate(&mol).unwrap().is_empty()),
    }
}

fn pipeline_run() -> Call<(Pipeline, String)> {
    Call {
        make: |_| {
            (
                Pipeline::default(),
                mol_to_mol_block(&depicted("CC(=O)[O-].[Na+]")).unwrap(),
            )
        },
        run: |(pipeline, block)| {
            pipeline.run(&block).unwrap();
        },
    }
}

// The catalogs and MCS.
fn filter_catalog_match() -> Call<(FilterCatalog, Mol)> {
    Call {
        make: |_| {
            let mut catalog = FilterCatalog::default();
            let carboxyl = SmartsMatcher::from_smarts("carboxyl", "C(=O)[OH]");
            catalog.add_entry(&FilterCatalogEntry::new("carboxyl", &carboxyl).unwrap());
            (catalog, mol(ASPIRIN))
        },
        run: |(catalog, mol)| assert!(catalog.has_match(&mol).unwrap()),
    }
}

fn fragments_cataloged() -> Call<(FragCatalog, Mol)> {
    Call {
        make: |_| {
            let params = FragCatParams::new(1, 3, groups()).unwrap();
            (FragCatalog::new(&params), mol("OCC=CC(=O)O"))
        },
        run: |(mut catalog, mol)| {
            let added = FragCatGenerator
                .add_frags_from_mol(&mol, &mut catalog)
                .unwrap();
            assert_eq!(added, 3);
        },
    }
}

fn maximum_common_substructure() -> Call<Vec<Mol>> {
    Call {
        make: |_| vec![mol("c1ccccc1CCO"), mol("c1ccccc1CCN")],
        run: |mols| assert_eq!(find_mcs(&mols).unwrap().num_atoms(), 8),
    }
}

fn maximum_common_edge_substructure() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (mol(CDK2[1]), mol(CDK2[3])),
        run: |(mol1, mol2)| assert!(!find_mces(&mol1, &mol2).unwrap().is_empty()),
    }
}

// Splitting a cluster computes more MCESs.
fn rascal_clusters() -> Call<Vec<Mol>> {
    Call {
        make: |_| CDK2.iter().map(|smiles| mol(smiles)).collect(),
        run: |mols| assert_eq!(rascal_cluster(&mols).unwrap()[0], [3, 4, 1, 2, 0]),
    }
}

fn tautomer_query_matches() -> Call<(TautomerQuery, Mol)> {
    Call {
        make: |_| {
            let query = TautomerQuery::new(&mol("O=C1CCCCC1")).unwrap();
            (query, mol("OC1=CCCCC1"))
        },
        run: |(query, target)| assert_eq!(query.get_substruct_matches(&target).unwrap().len(), 1),
    }
}

fn extended_query_matches() -> Call<(ExtendedQueryMol, Mol)> {
    Call {
        make: |_| {
            let query = mol("COCc1n[nH]c(F)c1 |LN:1:1.3|");
            (
                create_extended_query_mol(&query).unwrap(),
                mol("COOCc1n[nH]c(F)c1"),
            )
        },
        run: |(query, target)| {
            assert_eq!(
                rd_generalized_substruct::mol_get_substruct_matches(&target, &query)
                    .unwrap()
                    .len(),
                1
            );
        },
    }
}

fn substruct_library_matches() -> Call<(SubstructLibrary, Mol)> {
    Call {
        make: |_| {
            let mut library = SubstructLibrary::new(MolHolder::default().into(), None, None);
            for smiles in [ASPIRIN, "CCO", "c1ccccc1"] {
                library.add_mol(&mol(smiles)).unwrap();
            }
            (library, mol_from_smarts("[OX2H]").unwrap())
        },
        run: |(library, query)| assert_eq!(library.get_matches(&query).unwrap(), [0, 1]),
    }
}

// Decomposing also scores, under rgroup_mutex.
fn rgroups_decomposed() -> Call<(RGroupDecomposition, Vec<Mol>)> {
    Call {
        make: |_| {
            let decomposition = RGroupDecomposition::new(&[mol("c1ccccc1")]).unwrap();
            (decomposition, vec![mol("c1ccccc1Cl"), mol("c1ccccc1CCO")])
        },
        run: |(mut decomposition, mols)| {
            for mol in &mols {
                decomposition.add(mol).unwrap();
            }
            assert!(decomposition.process());
        },
    }
}

fn scaffold_network() -> Call<Vec<Mol>> {
    Call {
        make: |_| vec![mol("c1ccccc1CC1NC(=O)CCC1"), mol("c1cccnc1CC1NC(=O)CCC1")],
        run: |mols| {
            let network =
                create_scaffold_network(&mols, &ScaffoldNetworkParams::default()).unwrap();
            assert!(!network.nodes.is_empty());
        },
    }
}

// PubChem's shape code sets its overlap cutoff in globals and fills a static
// list of color patterns (pubchem_shape_mutex), whose matches construct
// ControlCHandlers.
fn pubchem_shape_alignment() -> Call<(Mol, Mol)> {
    Call {
        make: |_| (embedded("CCO"), embedded("CCN")),
        run: |(reference, mut probe)| {
            let (shape, _) = rd_shape_align::align_mol(&reference, &mut probe).unwrap();
            assert!(shape > 0.0);
        },
    }
}

fn pubchem_shape_inputs() -> Call<Mol> {
    Call {
        make: |_| embedded("Oc1ccccc1"),
        run: |mol| {
            let shape = rd_shape_align::prepare_conformer(&mol).unwrap();
            assert!(shape.sov() > 0.0);
        },
    }
}

// An H-bond field finds its special atoms with a serine pattern.
fn mif_h_bonds() -> Call<Mol> {
    Call {
        make: |_| embedded("NC(CO)C(=O)O"),
        run: |mol| {
            HBond::new(&mol).unwrap();
        },
    }
}

fn mol_features() -> Call<(MolChemicalFeatureFactory, Mol)> {
    Call {
        make: |_| {
            let fdef =
                "DefineFeature Carbonyl C=O\n  Family Acceptor\n  Weights 1.0,3.0\nEndFeature\n";
            (
                build_feature_factory_from_string(fdef).unwrap(),
                mol(ASPIRIN),
            )
        },
        run: |(factory, mol)| assert_eq!(factory.get_num_mol_features(&mol).unwrap(), 2),
    }
}

// RDKit parses the default abbreviations into a function-local static,
// whose thread-safe initialization DRD cannot see, so `make` reads them.
fn condensed_abbreviations() -> Call<(Vec<AbbreviationDefinition>, Mol)> {
    Call {
        make: |_| (get_default_abbreviations(), mol("CCCCCCCCc1ccccc1C(F)(F)F")),
        run: |(abbrevs, mol)| {
            let condensed = condense_mol_abbreviations(&mol, &abbrevs).unwrap();
            assert!(condensed.get_num_atoms() < mol.get_num_atoms());
        },
    }
}

// Link nodes are enumerated with reactions.
fn enumerated_link_nodes() -> Call<Mol> {
    Call {
        make: |_| mol("COC |LN:1:1.3|"),
        run: |mol| assert_eq!(enumerate(&mol).unwrap().size(), 3),
    }
}

// Embedding each isomer draws from RDKit's generator, with a
// ControlCHandler alive.
fn embedded_stereoisomers() -> Call<Mol> {
    Call {
        make: |_| mol("CC1CC2CCC1C2"),
        run: |mol| {
            let options = StereoEnumerationOptions {
                try_embedding: true,
                ..Default::default()
            };
            let isomers = StereoisomerEnumerator::new_with_params(&mol, &options).unwrap();
            assert_eq!(isomers.count(), 4);
        },
    }
}

// MolToInchi's reverse cleanup matches a pattern on every call.
fn inchi_written() -> Call<Mol> {
    Call {
        make: |_| mol(ASPIRIN),
        run: |mol| assert!(mol_to_inchi(&mol).unwrap().inchi.starts_with("InChI=1S/")),
    }
}

// InchiToMol's cleanup matches patterns around neutral five-valent
// nitrogens, which InChI gives for nitro groups.
fn inchi_read() -> Call<String> {
    Call {
        make: |_| "InChI=1S/CH3NO2/c1-2(3)4/h1H3".to_string(),
        run: |inchi| assert!(inchi_to_mol(&inchi).unwrap().mol.is_some()),
    }
}

// The ChemDraw library fills its attribute-name maps on the first CDXML
// write (chemdraw_mutex), which this phase must make.
fn cdxml_blocks() -> Call<Mol> {
    Call {
        make: |_| depicted("CCO"),
        run: |mol| {
            let block = rdkit_cxx::chem::rd_chem_draw::mol_to_chem_draw_block(&mol).unwrap();
            assert!(block.contains("BondLength"));
        },
    }
}

/// A synthon space of benzamides, written by `main`.
const SYNTHONS: &str = "SMILES\tsynton_id\tsynton#\treaction_id
c1ccccc1C(=O)[1*]\t1-1\t1\tamide
[1*]N1CCCC1\t2-1\t2\tamide
[1*]NC(C)CC\t2-2\t2\tamide
";

fn synthons() -> PathBuf {
    std::env::temp_dir().join(format!(
        "rdkit-cxx-threads-{}-synthons.txt",
        std::process::id()
    ))
}

// Reading the space and searching it construct ControlCHandlers.
fn synthon_substructure_search() -> Call<Mol> {
    Call {
        make: |_| mol_from_smarts("c1ccccc1C(=O)N").unwrap(),
        run: |query| {
            let mut space = SynthonSpace::default();
            space.read_text_file(synthons()).unwrap();
            let hits = space.substructure_search(&query).unwrap();
            assert_eq!(hits.get_hit_molecules().len(), 2);
        },
    }
}

// The synthons' fingerprints, built before the search, use feature
// invariants, which match their feature patterns.
fn synthon_fingerprint_search() -> Call<(SynthonSpace, FingerprintGenerator64, Mol)> {
    Call {
        make: |_| {
            let mut space = SynthonSpace::default();
            space.read_text_file(synthons()).unwrap();
            let features = GetMorganGeneratorParams {
                atom_invariants_generator: Some(get_morgan_feature_atom_inv_gen()),
                ..Default::default()
            };
            (
                space,
                get_morgan_generator_with_params(&features).unwrap(),
                mol("c1ccccc1C(=O)N1CCCC1"),
            )
        },
        run: |(mut space, generator, query)| {
            let hits = space.fingerprint_search(&query, &generator).unwrap();
            assert!(!hits.get_hit_molecules().is_empty());
        },
    }
}

// RASCAL matches substructures.
fn synthon_rascal_search() -> Call<(SynthonSpace, Mol)> {
    Call {
        make: |_| {
            let mut space = SynthonSpace::default();
            space.read_text_file(synthons()).unwrap();
            (space, mol("c1ccccc1C(=O)NC(C)CC"))
        },
        run: |(mut space, query)| {
            let hits = space.rascal_search(&query, &Default::default()).unwrap();
            assert!(!hits.get_hit_molecules().is_empty());
        },
    }
}

// Murtagh's clustering keeps its state in statics (murtagh_mutex), which
// HierarchicalClusterPicker's clustering shares.
fn murtagh_clusters() -> Call<Vec<Vec<f64>>> {
    Call {
        make: |i| {
            (0..6)
                .map(|j| vec![f64::from(j), f64::from(i + j * j)])
                .collect()
        },
        run: |data| {
            murtagh_cluster(&data, 6, 2, 1).unwrap();
        },
    }
}

fn hierarchical_picks() -> Call<Vec<f64>> {
    Call {
        make: |i| (0..15).map(|j| f64::from(1 + (i + j) % 7)).collect(),
        run: |dist_mat| {
            let picker = HierarchicalClusterPicker::new(ClusterMethod::Ward);
            assert_eq!(picker.pick(&dist_mat, 6, 3).unwrap().len(), 3);
        },
    }
}

// The rest.

// The centroid functions share a cache in RDKit; each thread's radius adds
// its own entry.
fn grid_centroids() -> Call<(UniformGrid3D, f64)> {
    Call {
        make: |i| {
            let mut grid = UniformGrid3D::new(10.0, 10.0, 10.0).unwrap();
            grid.set_sphere_occupancy(&Point3D::default(), 1.5, 0.25)
                .unwrap();
            (grid, 1.0 + 0.5 * f64::from(i))
        },
        run: |(grid, radius)| {
            for _ in 0..5 {
                find_grid_terminal_points(&grid, radius, 0.5).unwrap();
            }
            let (weight_sum, _) =
                compute_grid_centroid(&grid, &Point3D::default(), radius).unwrap();
            assert!(weight_sum > 0.0);
        },
    }
}

// Registering a property writes RDKit's global registry.
fn registered_property() -> Call<(u32, Mol)> {
    Call {
        make: |i| (i, mol("CCO")),
        run: |(i, mol)| {
            let name = format!("rustThread{i}");
            let functor = PropertyFunctor::new(&name, "1.0.0", move |_| f64::from(i));
            Properties::register_property(&functor);
            let properties = Properties::new(&[name.as_str()]).unwrap();
            assert_eq!(
                properties.compute_properties(&mol).unwrap(),
                vec![f64::from(i)]
            );
        },
    }
}

fn available_properties() -> Call<()> {
    Call {
        make: |_| (),
        run: |()| assert!(Properties::get_available_properties().contains(&"NumHBD".to_string())),
    }
}

// DCLV keeps its state in each call in RDKit 2025.03, so calls run at once.
fn dclv_surfaces() -> Call<Mol> {
    Call {
        make: |_| {
            let mut mol = mol("CCO");
            let mut conf = Conformer::new(3);
            for atom in 0..3 {
                let x = 1.5 * f64::from(atom);
                conf.set_atom_position(atom, Point3D::new(x, 0.3 * x, 0.0))
                    .unwrap();
            }
            mol.add_conformer(&conf).unwrap();
            mol
        },
        run: |mol| {
            let dclv = DoubleCubicLatticeVolume::new(&mol).unwrap();
            assert!(dclv.get_surface_area() > 0.0);
            assert!(dclv.get_volume() > 0.0);
        },
    }
}
