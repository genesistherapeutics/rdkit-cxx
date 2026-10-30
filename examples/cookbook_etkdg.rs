// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Conformer Generation with ETKDG (RDKitCB_25): tricks for generating
//! conformers with ETKDG.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#conformer-generation-with-etkdg),
//! by Shuzhe Wang.
//!
//! ETKDG (Riniker and Landrum, J. Chem. Inf. Model. 2015, 55, 2562) corrects
//! distance-geometry conformers with torsion preferences from the Cambridge
//! Structural Database. Its options live in `EmbedParameters`.

use std::collections::BTreeMap;

use rdkit_cxx::Result;
use rdkit_cxx::chem::rd_dist_geom::{
    embed_multiple_confs_with_params, etkdg, etkdg_v2, etkdg_v3, get_molecule_bounds_matrix,
    sr_etkdg_v3,
};
use rdkit_cxx::chem::{add_hs, mol_from_smiles};
use rdkit_cxx::distance_geometry::dist_geom::do_triangle_smoothing;

fn main() -> Result<()> {
    // ETKDG is the default conformer generation routine:
    let _params = etkdg();

    // Newer torsion potentials (Guba et al., J. Chem. Inf. Model. 2016, 56,
    // 1):
    let _params = etkdg_v2();

    // Improvements for small rings and macrocycles (Wang et al., J. Chem.
    // Inf. Model. 2020, 60, 2044). This one adds small ring torsion
    // potentials:
    let _params = sr_etkdg_v3();
    // and this one macrocycle torsion potentials and macrocycle handles:
    let mut params = etkdg_v3();
    // To use the two together:
    params.use_small_ring_torsions = true;

    // a macrocycle attached to a small ring
    let mut mol = add_hs(&mol_from_smiles(
        "C(OC(CCCCCCC(OCCSC(CCCCCC1)=O)=O)OCCSC1=O)N1CCOCC1",
    )?)?;
    let ids = embed_multiple_confs_with_params(&mut mol, 3, &mut params)?;
    println!("{ids:?}");
    assert_eq!(ids.len(), 3);

    // Changing a molecule's bounds matrix confines its conformational space:
    let mut mol = add_hs(&mol_from_smiles("C1CCC1C")?)?;
    let mut bm = get_molecule_bounds_matrix(&mol)?;
    bm[0][3] = 1.21;
    bm[3][0] = 1.20;
    bm[2][3] = 1.21;
    bm[3][2] = 1.20;
    bm[4][3] = 1.21;
    bm[3][4] = 1.20;
    // These bounds are inconsistent, so smoothing reports failure (false),
    // in Python too; the recipe goes on regardless.
    let smoothed = do_triangle_smoothing(&mut bm)?;
    println!("{smoothed}");
    assert!(!smoothed);
    params.set_bounds_mat(&bm)?;

    // Custom pairwise Coulombic interactions (CPCIs) refine the embedded
    // conformers. This sets a repulsive interaction of 0.9 e^2 between atoms
    // 0 and 3, to keep them further apart:
    params.set_cpci(&BTreeMap::from([((0, 3), 0.9)]));

    params.use_random_coords = true;
    let ids = embed_multiple_confs_with_params(&mut mol, 3, &mut params)?;
    println!("{ids:?}");
    assert_eq!(ids.len(), 3);
    Ok(())
}
