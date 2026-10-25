// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Runs every wrapper that locks RDKit's global state (CONVENTIONS.md U4)
//! from several threads at once. CI also runs this binary under Valgrind's
//! DRD, which sees RDKit's accesses to its globals even in conda's
//! uninstrumented build, and fails on any the locks leave unordered. It is a
//! plain `main` (`harness = false` in Cargo.toml), so libtest's threads add
//! nothing to DRD's report. To run the check as CI does, with Valgrind
//! installed:
//!
//! ```sh
//! CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUNNER="valgrind --tool=drd --log-file=target/drd.log" \
//!   cargo test --test threads
//! python3 tests/check_drd_log.py target/drd.log
//! ```

use std::collections::BTreeSet;
use std::thread;

use rdkit_cxx::chem::rd_mol_descriptors::{
    DoubleCubicLatticeVolume, Properties, PythonPropertyFunctor,
};
use rdkit_cxx::chem::{
    Conformer, MolzipParams, SmilesWriteParams, fragment_on_brics_bonds, mol_from_smiles,
    mol_to_mol_block, mol_to_smiles, mol_to_smiles_with_params, molzip_fragments_with_params,
};
use rdkit_cxx::geometry::{
    Point3D, UniformGrid3D, compute_grid_centroid, find_grid_terminal_points,
};

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
    thread::scope(|s| {
        let handles = [0, 1, 2, 3].map(|i| {
            let grid = &grid;
            s.spawn(move || {
                random_smiles();
                brics();
                grid_centroids(grid, 1.0 + 0.5 * f64::from(i));
                depiction();
                descriptors(i);
            })
        });
        // Joining each thread, rather than leaving it to the scope, orders its
        // accesses before the destructors that run at exit, for DRD.
        for handle in handles {
            handle.join().unwrap();
        }
    });
}

// Random SMILES draw from RDKit's global random number generator.
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
