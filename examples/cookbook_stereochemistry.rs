// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Identifying Stereochemistry (RDKitCB_30): find chiral centers and double
//! bond stereochemistry.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#identifying-stereochemistry),
//! by Vincent Scalfani.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    BondStereo, BondType, FindMolChiralCentersParams, StereoDescriptor, StereoSpecified,
    StereoType, find_mol_chiral_centers_with_params, find_potential_stereo, mol_from_smiles,
};

fn main() -> Result<()> {
    let mut m = mol_from_smiles("C[C@H]1CCC[C@@H](C)[C@@H]1Cl")?;

    // legacy FindMolChiralCenters()
    let legacy_params = FindMolChiralCentersParams {
        force: true,
        include_unassigned: true,
        use_legacy_implementation: Some(true),
        ..Default::default()
    };
    let legacy = find_mol_chiral_centers_with_params(&mut m, &legacy_params)?;
    println!("{legacy:?}");
    let expected = |codes: [(u32, &str); 3]| codes.map(|(idx, code)| (idx, code.to_owned()));
    assert_eq!(legacy, expected([(1, "S"), (5, "R"), (7, "R")]));

    // new stereochemistry code
    let params = FindMolChiralCentersParams {
        use_legacy_implementation: Some(false),
        ..legacy_params
    };
    let centers = find_mol_chiral_centers_with_params(&mut m, &params)?;
    println!("{centers:?}");
    assert_eq!(centers, expected([(1, "S"), (5, "R"), (7, "r")]));

    // Identifying Double Bond Stereochemistry
    let mut mol = mol_from_smiles(r"C\C=C(/F)\C(=C\F)\C=C")?;

    // Using GetStereo()
    let mut stereo = Vec::new();
    for b in mol.get_bonds() {
        println!(
            "{} {} {:?} {:?}",
            b.get_begin_atom_idx(),
            b.get_end_atom_idx(),
            b.get_bond_type(),
            b.get_stereo()
        );
        stereo.push((b.get_bond_type(), b.get_stereo()));
    }
    use BondStereo::{Stereoe, Stereonone, Stereoz};
    use BondType::{Double, Single};
    assert_eq!(
        stereo,
        [
            (Single, Stereonone),
            (Double, Stereoz),
            (Single, Stereonone),
            (Single, Stereonone),
            (Double, Stereoe),
            (Single, Stereonone),
            (Single, Stereonone),
            (Double, Stereonone),
        ]
    );

    // Double bond configuration can also be identified with new
    // stereochemistry code using Chem.FindPotentialStereo()
    let si = find_potential_stereo(&mut mol)?;
    for element in &si {
        println!(
            "  Type: {:?}, Which: {}, Specified: {:?}, Descriptor: {:?}",
            element.type_, element.centered_on, element.specified, element.descriptor
        );
    }
    let found: Vec<_> = si
        .iter()
        .map(|e| (e.type_, e.centered_on, e.specified, e.descriptor))
        .collect();
    assert_eq!(
        found,
        [
            (
                StereoType::BondDouble,
                1,
                StereoSpecified::Specified,
                StereoDescriptor::BondCis
            ),
            (
                StereoType::BondDouble,
                4,
                StereoSpecified::Specified,
                StereoDescriptor::BondTrans
            ),
        ]
    );
    Ok(())
}
