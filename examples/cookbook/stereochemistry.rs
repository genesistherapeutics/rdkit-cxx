// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Stereochemistry" chapter.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{
    AssignCIPLabelsParams, BondStereo, BondType, Mol, StereoDescriptor, StereoSpecified,
    StereoType, assign_cip_labels_with_params, find_potential_stereo, mol_from_smiles,
};

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[(
    "RDKitCB_30",
    "Identifying Stereochemistry",
    identifying_stereochemistry,
)];

/// The tetrahedral stereocenters and their CIP labels, with `?` for
/// unspecified centers. Python's `Chem.FindMolChiralCenters(mol,
/// includeUnassigned=True, useLegacyImplementation=False)`, which lives in
/// Python rather than in RDKit's C++.
fn find_mol_chiral_centers(mol: &mut Mol) -> Result<Vec<(u32, String)>> {
    let stereo = find_potential_stereo(mol)?;
    let label = |type_| -> Vec<u32> {
        stereo
            .iter()
            .filter(|info| info.type_ == type_)
            .map(|info| info.centered_on)
            .collect()
    };
    assign_cip_labels_with_params(
        mol,
        &AssignCIPLabelsParams {
            atoms_to_label: Some(label(StereoType::AtomTetrahedral)),
            bonds_to_label: Some(label(StereoType::BondDouble)),
            ..Default::default()
        },
    )?;
    let mut centers = Vec::new();
    for info in stereo.iter().filter(|info| info.type_ == StereoType::AtomTetrahedral) {
        let atom = mol.get_atom_with_idx(info.centered_on)?;
        let code = if atom.has_prop("_CIPCode") {
            atom.get_prop("_CIPCode")?
        } else {
            "?".to_owned()
        };
        centers.push((info.centered_on, code));
    }
    Ok(centers)
}

/// Finds chiral centers and double-bond stereochemistry. By Vincent Scalfani.
fn identifying_stereochemistry() -> Result<()> {
    let mut m = mol_from_smiles("C[C@H]1CCC[C@@H](C)[C@@H]1Cl")?;
    let centers = find_mol_chiral_centers(&mut m)?;
    println!("{centers:?}");
    // Atom 7 is pseudoasymmetric, so its label is lowercase. (Python's legacy
    // implementation labels it `R`.)
    let expected = [(1, "S"), (5, "R"), (7, "r")].map(|(idx, code)| (idx, code.to_owned()));
    assert_eq!(centers, expected);

    // The stereo of each bond.
    let mut mol = mol_from_smiles(r"C\C=C(/F)\C(=C\F)\C=C")?;
    let bonds: Vec<(u32, u32, BondType, BondStereo)> = mol
        .get_bonds()
        .map(|bond| {
            (
                bond.get_begin_atom_idx(),
                bond.get_end_atom_idx(),
                bond.get_bond_type(),
                bond.get_stereo(),
            )
        })
        .collect();
    for (begin, end, bond_type, stereo) in &bonds {
        println!("{begin} {end} {bond_type:?} {stereo:?}");
    }
    use BondStereo::{Stereoe, Stereonone, Stereoz};
    use BondType::{Double, Single};
    assert_eq!(
        bonds,
        [
            (0, 1, Single, Stereonone),
            (1, 2, Double, Stereoz),
            (2, 3, Single, Stereonone),
            (2, 4, Single, Stereonone),
            (4, 5, Double, Stereoe),
            (5, 6, Single, Stereonone),
            (4, 7, Single, Stereonone),
            (7, 8, Double, Stereonone),
        ]
    );

    // RDKit's newer stereo perception also finds double-bond configurations.
    let stereo = find_potential_stereo(&mut mol)?;
    for info in &stereo {
        println!(
            "  Type: {:?}, Which: {}, Specified: {:?}, Descriptor: {:?}",
            info.type_, info.centered_on, info.specified, info.descriptor
        );
    }
    let summary: Vec<_> = stereo
        .iter()
        .map(|info| (info.type_, info.centered_on, info.specified, info.descriptor))
        .collect();
    use StereoType::BondDouble;
    assert_eq!(
        summary,
        [
            (BondDouble, 1, StereoSpecified::Specified, StereoDescriptor::BondCis),
            (BondDouble, 4, StereoSpecified::Specified, StereoDescriptor::BondTrans),
        ]
    );
    Ok(())
}
