// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Runs the wrappers that lock RDKit's global state (CONVENTIONS.md U4)
//! from several threads at once: every one that writes the state a lock
//! guards, and at least one that only reads it. CI also runs this binary
//! under Valgrind's DRD, which sees RDKit's accesses to its globals even in
//! conda's uninstrumented build, and fails on any the locks leave unordered.
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
use std::path::Path;
use std::thread;

use rdkit_cxx::chem::rd_coord_gen::{add_coords, set_default_template_file_dir};
use rdkit_cxx::chem::rd_depictor::{
    Compute2DCoordsParams, add_ring_system_templates, compute_2d_coords_mimic_distmat,
    compute_2d_coords_with_params, generate_depiction_matching_3d_structure,
    is_coord_gen_support_available, load_default_ring_system_templates, set_prefer_coord_gen,
    set_ring_system_templates,
};
use rdkit_cxx::chem::rd_mol_descriptors::{
    DoubleCubicLatticeVolume, Properties, PythonPropertyFunctor,
};
use rdkit_cxx::chem::{
    Conformer, Mol, MolFragmentToSmilesParams, MolToCXSmilesParams, MolToRandomSmilesVectParams,
    MolzipParams, PNGMetadataParams, SmilesWriteParams, fragment_on_brics_bonds,
    mol_fragment_to_cx_smiles_with_params, mol_fragment_to_smiles_with_params, mol_from_smiles,
    mol_metadata_to_png_string_with_params, mol_to_cx_smiles_with_params, mol_to_mol_block,
    mol_to_random_smiles_vect_with_params, mol_to_smiles, mol_to_smiles_with_params,
    molzip_fragments_with_params,
};
use rdkit_cxx::geometry::{
    Point3D, UniformGrid3D, compute_grid_centroid, find_grid_terminal_points,
};

/// A 1x1 grayscale PNG for the metadata writer.
const PNG: &[u8] = &[
    137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73, 72, 68, 82, 0, 0, 0, 1, 0, 0, 0, 1, 8, 0, 0,
    0, 0, 58, 126, 155, 85, 0, 0, 0, 10, 73, 68, 65, 84, 120, 156, 99, 96, 0, 0, 0, 2, 0, 1, 72,
    175, 164, 113, 0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130,
];

fn main() {
    // RDKit fills some static data on first use (the periodic table, the
    // SMILES parser's tables), ordered by C++'s thread-safe initialization,
    // whose fast path DRD cannot see. Filling it before any thread starts
    // keeps those reads out of the report. Globals that only the locks order
    // (the BRICS patterns) are left for the threads to fill.
    mol_to_smiles(&mol_from_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap()).unwrap();
    let mut grid = UniformGrid3D::new(10.0, 10.0, 10.0).unwrap();
    grid.set_sphere_occupancy(&Point3D::default(), 1.5, 0.25)
        .unwrap();
    let templates = std::env::temp_dir().join(format!(
        "rdkit-cxx-threads-{}-templates.cxsmi",
        std::process::id()
    ));
    std::fs::write(
        &templates,
        "C1CC2CCC1C2 |(1.55254,0.75,;1.55254,-0.75,;0.125954,-1.21353,;-1.30063,-0.75,;\
         -1.30063,0.75,;0.125954,1.21353,;-0.755724,0,)|\n",
    )
    .unwrap();
    thread::scope(|s| {
        let handles = [0, 1, 2, 3].map(|i| {
            let grid = &grid;
            let templates = &templates;
            s.spawn(move || {
                random_smiles();
                brics();
                grid_centroids(grid, 1.0 + 0.5 * f64::from(i));
                depiction();
                depiction_settings(i, templates);
                descriptors(i);
            })
        });
        // Joining each thread, rather than leaving it to the scope, orders its
        // accesses before the destructors that run at exit, for DRD.
        for handle in handles {
            handle.join().unwrap();
        }
    });
    std::fs::remove_file(&templates).unwrap();
}

// Every SMILES writer that draws from, or seeds, RDKit's global random
// number generator.
fn random_smiles() {
    let random = SmilesWriteParams {
        do_random: true,
        ..Default::default()
    };
    let mol = mol_from_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let canonical = mol_to_smiles(&mol).unwrap();
    let mut written = BTreeSet::new();
    for _ in 0..50 {
        let smiles = mol_to_smiles_with_params(&mol, &random).unwrap();
        let parsed = mol_from_smiles(&smiles).unwrap();
        assert_eq!(mol_to_smiles(&parsed).unwrap(), canonical);
        written.insert(smiles);
    }
    assert!(written.len() > 1);
    let cx = MolToCXSmilesParams {
        params: random.clone(),
        ..Default::default()
    };
    mol_to_cx_smiles_with_params(&mol, &cx).unwrap();
    let fragment = MolFragmentToSmilesParams {
        params: random.clone(),
        ..Default::default()
    };
    mol_fragment_to_smiles_with_params(&mol, &[0, 1, 2], &fragment).unwrap();
    mol_fragment_to_cx_smiles_with_params(&mol, &[0, 1, 2], &fragment).unwrap();
    let seeded = MolToRandomSmilesVectParams {
        random_seed: 7,
        ..Default::default()
    };
    let smiles = mol_to_random_smiles_vect_with_params(&mol, 5, &seeded).unwrap();
    assert_eq!(smiles.len(), 5);
    let png = PNGMetadataParams {
        smiles_write_params: random,
        ..Default::default()
    };
    assert!(
        mol_metadata_to_png_string_with_params(&mol, PNG, &png)
            .unwrap()
            .len()
            > PNG.len()
    );
}

// The first BRICS call fills RDKit's pattern globals.
fn brics() {
    let mol = mol_from_smiles("CCCOCCC(=O)c1ccccc1").unwrap();
    assert_eq!(fragment_on_brics_bonds(&mol).unwrap().get_num_atoms(), 20);
}

// The centroid functions share a cache in RDKit; each thread's radius adds
// its own entry.
fn grid_centroids(grid: &UniformGrid3D, radius: f64) {
    for _ in 0..5 {
        find_grid_terminal_points(grid, radius, 0.5).unwrap();
    }
    let (weight_sum, _) = compute_grid_centroid(grid, &Point3D::default(), radius).unwrap();
    assert!(weight_sum > 0.0);
}

// A mol block for a molecule without a conformer, and molzip's
// generate_coordinates, compute 2D coordinates, which read 2D depiction's
// global settings.
fn depiction() {
    let chiral = mol_from_smiles("F[C@H](Cl)Br").unwrap();
    assert!(mol_to_mol_block(&chiral).unwrap().contains("V2000"));
    let coords = MolzipParams {
        generate_coordinates: true,
        ..Default::default()
    };
    let mut mols = vec![
        mol_from_smiles("c1ccccc1[*:1]").unwrap(),
        mol_from_smiles("O[*:1]").unwrap(),
    ];
    let zipped = molzip_fragments_with_params(&mut mols, &coords)
        .unwrap()
        .unwrap();
    assert_eq!(zipped.get_num_conformers(), 1);
}

// The calls that write 2D depiction's settings (the CoordGen preference, the
// bond length, the ring-system and CoordGen templates), and the depictor's
// sampling, which draws from RDKit's generator: random bond flips, and the
// distance mimicking that matching a 3D structure runs.
fn depiction_settings(i: u32, templates: &Path) {
    set_prefer_coord_gen(i % 2 == 0);
    set_ring_system_templates(templates).unwrap();
    add_ring_system_templates(templates).unwrap();
    load_default_ring_system_templates();
    set_default_template_file_dir("").unwrap();
    // Rotatable bonds, so the samples have bonds to flip.
    let mut mol = mol_from_smiles("c1ccccc1CCC(=O)O").unwrap();
    let sampled = Compute2DCoordsParams {
        n_flips_per_sample: 2,
        n_sample: 5,
        bond_length: 1.2,
        force_rdkit: true,
        use_ring_templates: true,
        ..Default::default()
    };
    compute_2d_coords_with_params(&mut mol, &sampled).unwrap();
    let n = mol.get_num_atoms() as usize;
    compute_2d_coords_mimic_distmat(&mut mol, &vec![1.5; n * (n - 1) / 2]).unwrap();
    if is_coord_gen_support_available() {
        add_coords(&mut mol).unwrap();
    }
    let mut propanol = mol_from_smiles("CCCO").unwrap();
    generate_depiction_matching_3d_structure(&mut propanol, &propanol_3d()).unwrap();
    assert_eq!(propanol.get_num_conformers(), 1);
    set_prefer_coord_gen(false);
}

/// Propanol in 3D, from Python's EmbedMolecule(AddHs(CCCO), randomSeed=42)
/// without the hydrogens.
fn propanol_3d() -> Mol {
    let mut mol = mol_from_smiles("CCCO").unwrap();
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

// Registering a property writes RDKit's global property registry, which
// the other calls read; DCLV's surface tests share a global cache.
fn descriptors(i: u32) {
    let name = format!("rustThread{i}");
    let functor = PythonPropertyFunctor::new(&name, "1.0.0", move |_| f64::from(i));
    Properties::register_property(&functor);
    let mut mol = mol_from_smiles("CCO").unwrap();
    let properties = Properties::new(&[name.as_str()]).unwrap();
    assert_eq!(
        properties.compute_properties(&mol).unwrap(),
        vec![f64::from(i)]
    );
    assert!(Properties::get_available_properties().contains(&name));
    let mut conf = Conformer::new(3);
    for atom in 0..3 {
        let x = 1.5 * f64::from(atom);
        conf.set_atom_position(atom, Point3D::new(x, 0.3 * x, 0.0))
            .unwrap();
    }
    mol.add_conformer(&conf).unwrap();
    let mut dclv = DoubleCubicLatticeVolume::new(&mol).unwrap();
    assert!(dclv.get_surface_area() > 0.0);
    assert!(dclv.get_polar_surface_area() > 0.0);
}
