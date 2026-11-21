// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Rings, Aromaticity, and Kekulization" chapter.

use std::collections::BTreeSet;

use rdkit_cxx::Result;
use rdkit_cxx::chem::rdqueries::{
    hybridization_equals_query_atom, is_aliphatic_query_atom, is_aromatic_query_atom,
};
use rdkit_cxx::chem::{Mol, mol_from_smarts, mol_from_smiles};

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[
    ("RDKitCB_3", "Count Ring Systems", count_ring_systems),
    ("RDKitCB_8", "Identify Aromatic Rings", identify_aromatic_rings),
    ("RDKitCB_9", "Identify Aromatic Atoms", identify_aromatic_atoms),
];

/// Groups rings that share atoms into ring systems. Rings that share one
/// atom, spiro rings, join a system only with `include_spiro`.
fn get_ring_systems(mol: &Mol, include_spiro: bool) -> Vec<BTreeSet<u32>> {
    let mut systems: Vec<BTreeSet<u32>> = Vec::new();
    for ring in mol.get_ring_info().atom_rings() {
        let mut ring_atoms: BTreeSet<u32> = ring.into_iter().collect();
        let mut new_systems = Vec::new();
        for system in systems {
            let num_in_common = ring_atoms.intersection(&system).count();
            if num_in_common > 0 && (include_spiro || num_in_common > 1) {
                ring_atoms.extend(system);
            } else {
                new_systems.push(system);
            }
        }
        new_systems.push(ring_atoms);
        systems = new_systems;
    }
    systems
}

/// Counts the ring systems in a molecule. By Greg Landrum.
fn count_ring_systems() -> Result<()> {
    let mol = mol_from_smiles("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")?;
    let systems = get_ring_systems(&mol, false);
    println!("{systems:?}");
    assert_eq!(
        systems,
        [
            BTreeSet::from([14, 15, 16, 17, 18, 19]),
            BTreeSet::from([1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
        ]
    );
    Ok(())
}

/// Whether every bond in a ring is aromatic.
fn is_ring_aromatic(mol: &Mol, bond_ring: &[u32]) -> Result<bool> {
    for &idx in bond_ring {
        if !mol.get_bond_with_idx(idx)?.get_is_aromatic() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Finds which rings are aromatic from their bonds. By Benjamin Datko and Greg
/// Landrum.
fn identify_aromatic_rings() -> Result<()> {
    let m = mol_from_smiles("c1cccc2c1CCCC2")?;
    let ri = m.get_ring_info();
    // The atoms that make up each ring, and the bonds.
    println!("{:?}", ri.atom_rings());
    assert_eq!(ri.atom_rings(), [[0, 1, 2, 3, 4, 5], [4, 5, 6, 7, 8, 9]]);
    let bond_rings = ri.bond_rings();
    println!("{bond_rings:?}");
    assert_eq!(bond_rings, [[0, 1, 2, 3, 4, 9], [4, 5, 6, 7, 8, 10]]);

    let aromatic = [
        is_ring_aromatic(&m, &bond_rings[0])?,
        is_ring_aromatic(&m, &bond_rings[1])?,
    ];
    println!("{aromatic:?}");
    assert_eq!(aromatic, [true, false]);
    Ok(())
}

/// Tells aromatic carbons from olefinic ones, with SMARTS and with query
/// atoms. By Paolo Tosco.
fn identify_aromatic_atoms() -> Result<()> {
    let mol = mol_from_smiles("c1ccccc1C=CCC")?;
    let aromatic_carbon = mol_from_smarts("c")?;
    let matches = mol.get_substruct_matches(&aromatic_carbon)?;
    println!("{matches:?}");
    assert_eq!(matches, [[0], [1], [2], [3], [4], [5]]);

    // RDKit's SMARTS extends atom primitives with hybridization: sp2
    // aliphatic carbons.
    let olefinic_carbon = mol_from_smarts("[C^2]")?;
    let matches = mol.get_substruct_matches(&olefinic_carbon)?;
    println!("{matches:?}");
    assert_eq!(matches, [[6], [7]]);

    // Query atoms skip parsing SMARTS.
    let indices = |query| -> Result<Vec<u32>> {
        let atoms = mol.get_atoms_matching_query(query)?;
        Ok(atoms.iter().map(|atom| atom.get_idx()).collect())
    };
    let aromatic = indices(&is_aromatic_query_atom())?;
    println!("{aromatic:?}");
    assert_eq!(aromatic, [0, 1, 2, 3, 4, 5]);

    // The query takes `HybridizationType`'s value in RDKit: 3 is SP2.
    let sp2 = hybridization_equals_query_atom(3)?;
    let sp2_indices = indices(&sp2)?;
    println!("{sp2_indices:?}");
    assert_eq!(sp2_indices, [0, 1, 2, 3, 4, 5, 6, 7]);

    let mut combined = is_aliphatic_query_atom();
    combined.expand_query(&sp2)?;
    let combined_indices = indices(&combined)?;
    println!("{combined_indices:?}");
    assert_eq!(combined_indices, [6, 7]);
    Ok(())
}
