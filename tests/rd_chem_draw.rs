// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_chem_draw::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdChemDraw.

#[test]
fn writes_cdxml_as_rdmolfiles_does() {
    let mol = mol_from_smiles("CCO").unwrap();
    let block = mol_to_chem_draw_block(&mol).unwrap();
    assert!(block.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>"));
    assert!(block.contains("<CDXML"));
    assert_eq!(block, mol_to_cdxml_block(&mol).unwrap());
    let bytes = mol_to_chem_draw_block_format(&mol, CDXFormat::Cdxml).unwrap();
    assert_eq!(bytes, block.into_bytes());
}

#[test]
fn writes_binary_cdx() {
    let mol = mol_from_smiles("CCO").unwrap();
    let bytes = mol_to_chem_draw_block_format(&mol, CDXFormat::Cdx).unwrap();
    // Python's MolToChemDrawBlock cannot decode it; rdmolfiles' bytes start
    // this way.
    assert!(bytes.starts_with(b"VjCD0100\x04\x03"));
    assert_eq!(
        bytes,
        mol_to_cdxml_block_format(&mol, CDXMLFormat::Cdx).unwrap()
    );
}

#[test]
fn refuses_coordinates_chem_draw_cannot_hold() {
    // A 2D molecule without bonds is scaled by 28.8 over a zero mean bond
    // length (rdmolfiles' check).
    let mut mol = mol_from_mol_block(
        "\n     RDKit          2D\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n    \
         0.0000    0.0000    0.0000 C   0  0\n    1.0000    0.0000    0.0000 O   0  0\n\
         M  END\n",
    )
    .unwrap();
    assert!(matches!(
        mol_to_chem_draw_block(&mol),
        Err(Error::Value { .. })
    ));
    // Computed coordinates are checked the same way.
    mol.remove_all_conformers();
    assert!(matches!(
        mol_to_chem_draw_block_format(&mol, CDXFormat::Cdx),
        Err(Error::Value { .. })
    ));
}

#[test]
fn writes_from_several_threads() {
    // Molecules without coordinates are depicted under the depiction lock.
    let handles: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                let mol = mol_from_smiles("c1ccccc1CCN").unwrap();
                for _ in 0..10 {
                    assert!(mol_to_chem_draw_block(&mol).unwrap().contains("<CDXML"));
                }
            })
        })
        .collect();
    for handle in handles {
        handle.join().unwrap();
    }
}
