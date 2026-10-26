// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Enumerate SMILES (RDKitCB_24): write random SMILES for one molecule,
//! reproducibly with a seed.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#enumerate-smiles),
//! by Guillaume Godin and Greg Landrum.

use std::collections::BTreeSet;

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    MolToRandomSmilesVectParams, SmilesWriteParams, mol_from_smiles, mol_to_random_smiles_vect,
    mol_to_random_smiles_vect_with_params, mol_to_smiles, mol_to_smiles_with_params,
};

fn main() -> Result<()> {
    let mol = mol_from_smiles("CC(N)C1CC1")?;

    // 100 random SMILES, without duplicates. Their order is random.
    let random = SmilesWriteParams {
        do_random: true,
        canonical: false,
        ..Default::default()
    };
    let smis = (0..100)
        .map(|_| mol_to_smiles_with_params(&mol, &random))
        .collect::<Result<BTreeSet<_>>>()?;
    println!("{smis:?}");
    // Each is the same molecule.
    let canonical = mol_to_smiles(&mol)?;
    for smi in &smis {
        assert_eq!(mol_to_smiles(&mol_from_smiles(smi)?)?, canonical);
    }
    assert!(smis.len() > 1);

    // By default the results are not reproducible:
    let m = mol_from_smiles("Oc1ncc(OC(CC)C)cc1")?;
    println!("{:?}", mol_to_random_smiles_vect(&m, 5)?);

    // But a random number seed makes them so:
    let seeded = MolToRandomSmilesVectParams {
        random_seed: 0xf00d,
        ..Default::default()
    };
    let s1 = mol_to_random_smiles_vect_with_params(&m, 5, &seeded)?;
    println!("{s1:?}");
    assert_eq!(
        s1,
        [
            "Oc1ccc(OC(CC)C)cn1",
            "CC(CC)Oc1cnc(O)cc1",
            "c1(O)ncc(cc1)OC(C)CC",
            "c1cc(cnc1O)OC(CC)C",
            "c1c(OC(CC)C)cnc(c1)O",
        ]
    );
    let s2 = mol_to_random_smiles_vect_with_params(&m, 5, &seeded)?;
    println!("{}", s2 == s1);
    assert_eq!(s2, s1);
    Ok(())
}
