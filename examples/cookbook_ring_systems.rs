// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Count Ring Systems (RDKitCB_3): group rings that share atoms into ring
//! systems.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#count-ring-systems),
//! by Greg Landrum.

use std::collections::BTreeSet;

use rdkit_cxx::Result;
use rdkit_cxx::chem::{Mol, mol_from_smiles};

/// Rings that share two or more atoms (or one, with `include_spiro`) belong
/// to the same system.
fn get_ring_systems(mol: &Mol, include_spiro: bool) -> Vec<BTreeSet<u32>> {
    let ri = mol.get_ring_info();
    let mut systems: Vec<BTreeSet<u32>> = Vec::new();
    for ring in ri.atom_rings() {
        let mut ring_ats: BTreeSet<u32> = ring.into_iter().collect();
        let mut n_systems = Vec::new();
        for system in systems {
            let n_in_common = ring_ats.intersection(&system).count();
            if n_in_common > 0 && (include_spiro || n_in_common > 1) {
                ring_ats.extend(system);
            } else {
                n_systems.push(system);
            }
        }
        n_systems.push(ring_ats);
        systems = n_systems;
    }
    systems
}

fn main() -> Result<()> {
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
