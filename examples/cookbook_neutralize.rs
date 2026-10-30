// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Neutralizing Molecules (RDKitCB_33): neutralize charged atoms by adding or
//! removing hydrogens, and compare with `rdMolStandardize.Uncharger`.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules),
//! by Noel O'Boyle, adapted for RDKit by Vincent Scalfani.
//!
//! `neutralize_atoms` neutralizes atoms with a +1 or -1 charge where it can,
//! ignoring the molecule's overall charge, and leaves charge-separated groups
//! (nitro groups, for example) alone. The `Uncharger` instead keeps a
//! negative charge that balances a positive charge it cannot neutralize.
//!
//! The recipe's comparison calls `un.uncharge(mol2)` and prints `mol2`, but
//! `uncharge` returns a new molecule and leaves its argument alone, so the
//! cookbook prints the inputs. This prints what the `Uncharger` returns,
//! which differs from `neutralize_atoms` only where a positive charge cannot
//! be neutralized.

use rdkit_cxx::Result;
use rdkit_cxx::chem::mol_standardize::rd_mol_standardize::Uncharger;
use rdkit_cxx::chem::{Mol, mol_from_smarts, mol_from_smiles, mol_to_smiles};

const SMI_LIST: [&str; 6] = [
    "CC(CNC[O-])[N+]([O-])=O",
    "C[N+](C)(C)CCC([O-])=O",
    "[O-]C1=CC=[N+]([O-])C=C1",
    "[O-]CCCN=[N+]=[N-]",
    "C[NH+](C)CC[S-]",
    "CP([O-])(=O)OC[NH3+]",
];

fn neutralize_atoms(mol: &mut Mol) -> Result<()> {
    let pattern = mol_from_smarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")?;
    let at_matches = mol.get_substruct_matches(&pattern)?;
    for at_idx in at_matches.iter().map(|m| m[0]) {
        let mut atom = mol.get_atom_with_idx_mut(at_idx)?;
        let chg = atom.get_formal_charge();
        let hcount = atom.get_total_num_hs()? as i32;
        atom.set_formal_charge(0);
        atom.set_num_explicit_hs((hcount - chg) as u32);
        atom.update_property_cache()?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let mut mols = SMI_LIST
        .iter()
        .map(|smi| mol_from_smiles(smi))
        .collect::<Result<Vec<_>>>()?;

    // Neutralize molecules by atom
    let mut neutralized = Vec::new();
    for mol in &mut mols {
        neutralize_atoms(mol)?;
        println!("{}", mol_to_smiles(mol)?);
        neutralized.push(mol_to_smiles(mol)?);
    }
    assert_eq!(
        neutralized,
        [
            "CC(CNCO)[N+](=O)[O-]",
            "C[N+](C)(C)CCC(=O)O",
            "[O-][n+]1ccc(O)cc1",
            "[N-]=[N+]=NCCCO",
            "CN(C)CCS",
            "CP(=O)(O)OCN",
        ]
    );

    // Compare to rdMolStandardize.Uncharger results:
    let un = Uncharger::default();
    let mut uncharged = Vec::new();
    for smi in SMI_LIST {
        let mol2 = un.uncharge(&mol_from_smiles(smi)?)?;
        println!("{}", mol_to_smiles(&mol2)?);
        uncharged.push(mol_to_smiles(&mol2)?);
    }
    assert_eq!(
        uncharged,
        [
            "CC(CNCO)[N+](=O)[O-]",
            "C[N+](C)(C)CCC(=O)[O-]",
            "[O-][n+]1ccc(O)cc1",
            "[N-]=[N+]=NCCCO",
            "CN(C)CCS",
            "CP(=O)(O)OCN",
        ]
    );
    Ok(())
}
