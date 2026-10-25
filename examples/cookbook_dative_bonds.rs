// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Organometallics with Dative Bonds (RDKitCB_19): replace single bonds to
//! transition metals with dative bonds where the ligand atom has too many
//! bonds.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#organometallics-with-dative-bonds),
//! by Greg Landrum.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    AtomRef, BondType, Mol, RWMol, SmilesParserParams, ValenceType, get_periodic_table,
    mol_from_smiles_with_params, mol_to_smiles,
};

fn is_transition_metal(at: &AtomRef<'_>) -> bool {
    let n = at.get_atomic_num();
    (22..=29).contains(&n) || (40..=47).contains(&n) || (72..=79).contains(&n)
}

/// Replaces single bonds between metals and atoms whose atomic numbers are
/// in `from_atoms` with dative bonds, when the atom has "too many" bonds.
fn set_dative_bonds(mol: &Mol, from_atoms: &[i32]) -> Result<RWMol> {
    let pt = get_periodic_table();
    let mut rwmol = RWMol::from(mol.clone());
    rwmol.update_property_cache_with_strict(false)?;
    // Find the bonds first: editing them would invalidate the views.
    let mut to_replace = Vec::new();
    for metal in rwmol.get_atoms().filter(is_transition_metal) {
        for nbr in metal.get_neighbors() {
            let bond = rwmol
                .get_bond_between_atoms(nbr.get_idx(), metal.get_idx())?
                .expect("neighbors are bonded");
            if from_atoms.contains(&nbr.get_atomic_num())
                && i64::from(nbr.get_valence(ValenceType::Explicit)?)
                    > i64::from(pt.get_default_valence(nbr.get_atomic_num() as u32)?)
                && bond.get_bond_type() == BondType::Single
            {
                to_replace.push((nbr.get_idx(), metal.get_idx()));
            }
        }
    }
    for (nbr, metal) in to_replace {
        rwmol.remove_bond(nbr, metal)?;
        rwmol.add_bond_with_order(nbr, metal, BondType::Dative)?;
    }
    Ok(rwmol)
}

fn main() -> Result<()> {
    let unsanitized = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let m = mol_from_smiles_with_params("CN(C)(C)[Pt]", &unsanitized)?;
    let m2 = set_dative_bonds(&m, &[7, 8])?;

    // The bond between nitrogen and platinum:
    let bond_type = m2
        .get_bond_between_atoms(1, 4)?
        .expect("N and Pt are bonded")
        .get_bond_type();
    println!("{bond_type:?}");
    assert_eq!(bond_type, BondType::Dative);

    // It also shows up in the SMILES, an RDKit extension:
    println!("{}", mol_to_smiles(&m2)?);
    assert_eq!(mol_to_smiles(&m2)?, "C[N](C)(C)->[Pt]");
    Ok(())
}
