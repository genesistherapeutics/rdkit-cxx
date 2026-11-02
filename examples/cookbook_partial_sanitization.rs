// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Explicit Valence Error - Partial Sanitization (RDKitCB_15): read a
//! molecule without the valence check, then sanitize it partially.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#explicit-valence-error-partial-sanitization),
//! by Greg Landrum. Use caution, and make sure your molecules actually make
//! sense before doing this!

use rdkit_cxx::Result;
use rdkit_cxx::chem::rdqueries::atom_num_equals_query_atom;
use rdkit_cxx::chem::{
    HybridizationType, SanitizeFlags, SanitizeMolParams, SmilesParserParams,
    mol_from_smiles_with_params, sanitize_mol_with_params,
};

fn main() -> Result<()> {
    // RDKit rejects hypervalent P by default, so skip sanitization:
    let unsanitized = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let mut m = mol_from_smiles_with_params("F[P-](F)(F)(F)(F)F.CN(C)C(F)=[N+](C)C", &unsanitized)?;

    // RDKit has not assigned a hybridization to the P (or any other atom).
    // Select the first and only P with a query:
    let q = atom_num_equals_query_atom(15)?;
    let phosphorus = m.get_atoms_matching_query(&q)?[0].get_idx();
    let hybridization = m.get_atom_with_idx(phosphorus)?.get_hybridization();
    println!("{hybridization:?}");
    assert_eq!(hybridization, HybridizationType::Unspecified);

    // Regenerate computed properties like implicit valence and ring
    // information
    m.update_property_cache_with_strict(false)?;

    // Apply several sanitization rules; setting the hybridization is key.
    let params = SanitizeMolParams {
        sanitize_ops: SanitizeFlags::SANITIZE_FINDRADICALS
            | SanitizeFlags::SANITIZE_KEKULIZE
            | SanitizeFlags::SANITIZE_SETAROMATICITY
            | SanitizeFlags::SANITIZE_SETCONJUGATION
            | SanitizeFlags::SANITIZE_SETHYBRIDIZATION
            | SanitizeFlags::SANITIZE_SYMMRINGS,
        catch_errors: true,
    };
    sanitize_mol_with_params(&mut m, &params)?;

    // The P is now SP3D2, so the six F would be octahedral in a depiction.
    let hybridization = m.get_atom_with_idx(phosphorus)?.get_hybridization();
    println!("{hybridization:?}");
    assert_eq!(hybridization, HybridizationType::Sp3d2);
    Ok(())
}
