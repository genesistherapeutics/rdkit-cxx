// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Contiguous Rotatable Bonds (RDKitCB_22): find the largest number of
//! contiguous rotatable bonds.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#contiguous-rotatable-bonds),
//! by Paolo Tosco.

use std::collections::BTreeSet;

use rdkit_cxx::Result;
use rdkit_cxx::chem::lipinski::rotatable_bond_smarts;
use rdkit_cxx::chem::{Mol, mol_from_smiles};

/// Groups of contiguous rotatable bonds, largest first.
fn find_bond_groups(mol: &Mol) -> Result<Vec<Vec<u32>>> {
    let rot_atom_pairs = mol.get_substruct_matches(&rotatable_bond_smarts())?;
    let mut rot_bond_set = BTreeSet::new();
    for ap in &rot_atom_pairs {
        let bond = mol
            .get_bond_between_atoms(ap[0], ap[1])?
            .expect("a rotatable bond joins each matched pair");
        rot_bond_set.insert(bond.get_idx());
    }
    let mut rot_bond_groups = Vec::new();
    while let Some(i) = rot_bond_set.pop_first() {
        let mut connected_bond_set = BTreeSet::from([i]);
        let mut stack = vec![i];
        while let Some(i) = stack.pop() {
            let b = mol.get_bond_with_idx(i)?;
            let mut bonds = Vec::new();
            for a in [b.get_begin_atom(), b.get_end_atom()] {
                bonds.extend(
                    a.get_bonds().iter().map(|b| b.get_idx()).filter(|idx| {
                        rot_bond_set.contains(idx) && !connected_bond_set.contains(idx)
                    }),
                );
            }
            connected_bond_set.extend(&bonds);
            stack.extend(bonds);
        }
        rot_bond_set.retain(|idx| !connected_bond_set.contains(idx));
        rot_bond_groups.push(connected_bond_set.into_iter().collect::<Vec<_>>());
    }
    rot_bond_groups.sort_by_key(|group| std::cmp::Reverse(group.len()));
    Ok(rot_bond_groups)
}

fn main() -> Result<()> {
    let mol = mol_from_smiles("CCC(CC(C)CC1CCC1)C(CC(=O)O)N")?;

    // Bond groups are sorted by decreasing size, so the first group's size
    // (if any) is the largest number of contiguous rotatable bonds.
    let bond_groups = find_bond_groups(&mol)?;
    let largest_n_cont_rot_bonds = bond_groups.first().map_or(0, Vec::len);
    println!("{largest_n_cont_rot_bonds}");
    println!("{bond_groups:?}");
    assert_eq!(largest_n_cont_rot_bonds, 8);
    assert_eq!(bond_groups, [[1, 2, 3, 5, 6, 10, 11, 12]]);
    Ok(())
}
