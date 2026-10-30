// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports rdkit/Chem/UnitTestEnumerateHeterocycles.py's example tests and the
// module's doctests; values from Python's module in RDKit 2024.03.6.

use rdkit_cxx::Result;
use rdkit_cxx::chem::enumerate_heterocycles::*;
use rdkit_cxx::chem::{mol_from_smarts, mol_from_smiles, mol_to_smiles};

#[test]
fn smarts_match_examples() -> Result<()> {
    for row in get_heterocycle_reaction_smarts() {
        if row.smarts.is_empty() {
            continue;
        }
        let substruct = mol_from_smarts(row.smarts)?;
        for smiles in row.example.split(',') {
            assert!(
                mol_from_smiles(smiles)?.has_substruct_match(&substruct)?,
                "{} not in {smiles}",
                row.smarts
            );
        }
        for smiles in row.negative_example.split(',').filter(|s| !s.is_empty()) {
            assert!(
                !mol_from_smiles(smiles)?.has_substruct_match(&substruct)?,
                "{} in {smiles}",
                row.smarts
            );
        }
    }
    Ok(())
}

#[test]
fn something_matches_every_negative_example() -> Result<()> {
    let rows = get_heterocycle_reaction_smarts();
    let substructs: Vec<_> = rows
        .iter()
        .filter(|r| !r.smarts.is_empty())
        .map(|r| mol_from_smarts(r.smarts).unwrap())
        .collect();
    for row in &rows {
        for smiles in row.negative_example.split(',').filter(|s| !s.is_empty()) {
            let mol = mol_from_smiles(smiles)?;
            assert!(
                substructs
                    .iter()
                    .any(|s| mol.has_substruct_match(s).unwrap()),
                "{smiles}"
            );
        }
    }
    Ok(())
}

/// The input, how deep, how many, the first few in Python's order, and the
/// last few sorted.
type Case<'a> = (&'a str, Option<u32>, usize, &'a [&'a str], [&'a str; 3]);

#[test]
fn enumerations() -> Result<()> {
    assert_eq!(get_heterocycle_reactions().len(), 41);
    let cases: &[Case] = &[
        (
            "c1ccccc1",
            None,
            8,
            &[
                "c1ccncc1", "c1ccnnc1", "c1cncnc1", "c1cnccn1", "c1ccccc1", "c1cnnnc1",
            ],
            ["c1cnncn1", "c1cnnnc1", "c1ncncn1"],
        ),
        (
            "n1ccccc1",
            Some(1),
            4,
            &["c1ccnnc1", "c1cncnc1", "c1cnccn1", "c1ccccc1"],
            ["c1ccnnc1", "c1cnccn1", "c1cncnc1"],
        ),
        (
            "c1ccoc1",
            None,
            30,
            &[
                "c1cnoc1",
                "c1cocn1",
                "c1cc[nH]c1",
                "c1ccsc1",
                "c1cnon1",
                "c1conn1",
            ],
            ["n1nn[nH]n1", "n1nnon1", "n1nnsn1"],
        ),
        (
            "c1ccc2[nH]ccc2c1",
            Some(1),
            10,
            &[
                "c1ccc2[nH]ncc2c1",
                "c1ccc2[nH]cnc2c1",
                "c1ccc2ocnc2c1",
                "c1ccc2scnc2c1",
                "c1cc2[nH]ccc2cn1",
                "c1cc2cc[nH]c2cn1",
            ],
            ["c1ccc2scnc2c1", "c1cnc2[nH]ccc2c1", "c1cnc2cc[nH]c2c1"],
        ),
        (
            "Cc1ccc(-c2ccncc2)cc1",
            Some(2),
            22,
            &[
                "Cc1ccc(-c2ccncc2)cn1",
                "Cc1ccc(-c2ccncc2)nc1",
                "Cc1ccc(-c2ccncn2)cc1",
                "Cc1ccc(-c2ccnnc2)cc1",
                "Cc1ccc(-c2ccccc2)cc1",
                "Cc1ccc(-c2ccncc2)nn1",
            ],
            [
                "Cc1cnc(-c2ccncc2)cn1",
                "Cc1cnc(-c2ccncc2)nc1",
                "Cc1ncc(-c2ccncc2)cn1",
            ],
        ),
        (
            "c1cc[nH]c1",
            None,
            30,
            &[
                "c1cn[nH]c1",
                "c1cnoc1",
                "c1cnsc1",
                "c1c[nH]cn1",
                "c1cocn1",
                "c1cscn1",
            ],
            ["n1nn[nH]n1", "n1nnon1", "n1nnsn1"],
        ),
    ];
    for &(smiles, depth, n, first, last) in cases {
        let mut smis: Vec<String> =
            enumerate_heterocycles_with_depth(&mol_from_smiles(smiles)?, depth)
                .map(|m| mol_to_smiles(&m.unwrap()).unwrap())
                .collect();
        assert_eq!(smis.len(), n, "{smiles}");
        assert_eq!(&smis[..first.len()], first, "{smiles}");
        smis.sort();
        assert_eq!(&smis[n - 3..], last, "{smiles}");
    }
    let all: Vec<String> = enumerate_heterocycles(&mol_from_smiles("c1ccccc1")?)
        .map(|m| mol_to_smiles(&m.unwrap()).unwrap())
        .collect();
    assert_eq!(all.len(), 8);
    Ok(())
}
