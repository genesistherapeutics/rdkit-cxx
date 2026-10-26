// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports SimpleEnum.Enumerator's doctests; values from Python's module in
// RDKit 2026.09.1. The reactions with atom lists are inlined, their
// fixed-width lines' trailing spaces escaped.

mod common;

use rdkit_cxx::chem::rd_chem_reactions::{
    ChemicalReaction, reaction_from_rxn_block, reaction_from_rxn_file,
};
use rdkit_cxx::chem::simple_enum::enumerator::*;
use rdkit_cxx::chem::{Mol, mol_from_smiles, mol_to_smiles};
use rdkit_cxx::{Error, Result};

const AZIDE_REACTION: &str = "$RXN

  Marvin       071401090920

  2  1
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

  4  3  0  0  0  0            999 V2000
   -2.3792   -0.6875    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   -1.6647   -0.2750    0.0000 N   0  0  0  0  0  0  0  0  0  2  0  0
   -0.9542    0.1458    0.0000 N   0  3  0  0  0  0  0  0  0  3  0  0
   -0.2417    0.5625    0.0000 N   0  5  0  0  0  0  0  0  0  4  0  0
  2  3  2  0  0  0  0
  1  2  1  0  0  0  0
  3  4  2  0  0  0  0
V    2 azide
M  CHG  2   3   1   4  -1
M  RGP  1   1   1
M  END
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

  4  3  1  0  0  0            999 V2000
    1.8375   -0.4375    0.0000 R#  0  0  0  0  0  0  0  0  0  5  0  0
    2.5520   -0.0250    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    3.2664   -0.4375    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    2.5520    0.8000    0.0000 O   0  0  0  0  0  0  0  0  0  7  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  2  0  0  0  0
  3 F    2   8  17
V    2 carboxylicacid,acidchloride
M  ALS   3  2 F O   Cl\x20\x20
M  RGP  1   1   2
M  END
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

 10 10  0  0  0  0            999 V2000
   10.5791   -0.8125    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0
   11.4041   -0.8125    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   11.7068    0.0258    0.0000 N   0  0  0  0  0  0  0  0  0  2  0  0
   10.9916    0.4583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3266   -0.0284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5419    0.2262    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.3700    1.0331    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    8.5852    1.2877    0.0000 O   0  0  0  0  0  0  0  0  0  7  0  0
    9.9828    1.5854    0.0000 R#  0  0  0  0  0  0  0  0  0  5  0  0
   12.4459    0.2254    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  2  0  0  0  0
  5  6  1  0  0  0  0
  2  3  1  0  0  0  0
  6  7  1  0  0  0  0
  3  4  1  0  0  0  0
  7  8  2  0  0  0  0
  4  5  2  0  0  0  0
  7  9  1  0  0  0  0
  5  1  1  0  0  0  0
  3 10  1  0  0  0  0
M  RGP  2   9   2  10   1
M  END
";

const BAD_VALUE1: &str = "$RXN

      ISIS     052820091627

  2  1
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -3.2730   -7.0542    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -3.9875   -7.4667    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
V    1 halogen.bromine.aromatic
M  RGP  1   2   1
M  END
$MOL

  -ISIS-  05280916272D

  4  3  0  0  0  0  0  0  0  0999 V2000
    3.4375   -7.7917    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
    4.1520   -7.3792    0.0000 B   0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -6.5542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8664   -7.7917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  1  0  0  0  0
V    2 boromicacid
M  RGP  1   1   2
M  END
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   11.2667   -7.3417    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   11.9811   -6.9292    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
  1  2  1  0  0  0  0
M  RGP  2   1   1   2   2
M  END
";

const BAD_VALUE2: &str = "$RXN

  Marvin       071401090920

  2  1
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

  4  3  0  0  0  0            999 V2000
   -2.3792   -0.6875    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   -1.6647   -0.2750    0.0000 N   0  0  0  0  0  0  0  0  0  2  0  0
   -0.9542    0.1458    0.0000 N   0  3  0  0  0  0  0  0  0  3  0  0
   -0.2417    0.5625    0.0000 N   0  5  0  0  0  0  0  0  0  4  0  0
  2  3  2  0  0  0  0
  1  2  1  0  0  0  0
  3  4  2  0  0  0  0
V    2 azide
M  CHG  2   3   1   4  -1
M  RGP  1   1   1
M  END
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

  4  3  1  0  0  0            999 V2000
    1.8375   -0.4375    0.0000 R#  0  0  0  0  0  0  0  0  0  5  0  0
    2.5520   -0.0250    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    3.2664   -0.4375    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    2.5520    0.8000    0.0000 O   0  0  0  0  0  0  0  0  0  7  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  2  0  0  0  0
  3 F    2   8  17
V    2 carboxylicacid,acidchlroide
M  ALS   3  2 F O   Cl\x20\x20
M  RGP  1   1   2
M  END
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

 10 10  0  0  0  0            999 V2000
   10.5791   -0.8125    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0
   11.4041   -0.8125    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   11.7068    0.0258    0.0000 N   0  0  0  0  0  0  0  0  0  2  0  0
   10.9916    0.4583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3266   -0.0284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5419    0.2262    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.3700    1.0331    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    8.5852    1.2877    0.0000 O   0  0  0  0  0  0  0  0  0  7  0  0
    9.9828    1.5854    0.0000 R#  0  0  0  0  0  0  0  0  0  5  0  0
   12.4459    0.2254    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  2  0  0  0  0
  5  6  1  0  0  0  0
  2  3  1  0  0  0  0
  6  7  1  0  0  0  0
  3  4  1  0  0  0  0
  7  8  2  0  0  0  0
  4  5  2  0  0  0  0
  7  9  1  0  0  0  0
  5  1  1  0  0  0  0
  3 10  1  0  0  0  0
M  RGP  2   9   2  10   1
M  END
";

const BAD_VALUE3: &str = "$RXN

  Marvin       071401090920

  2  1
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

  4  3  0  0  0  0            999 V2000
   -2.3792   -0.6875    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   -1.6647   -0.2750    0.0000 N   0  0  0  0  0  0  0  0  0  2  0  0
   -0.9542    0.1458    0.0000 N   0  3  0  0  0  0  0  0  0  3  0  0
   -0.2417    0.5625    0.0000 N   0  5  0  0  0  0  0  0  0  4  0  0
  2  3  2  0  0  0  0
  1  2  1  0  0  0  0
  3  4  2  0  0  0  0
V    2 azide
M  CHG  2   3   1   4  -1
M  RGP  1   1   1
M  END
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

  4  3  1  0  0  0            999 V2000
    1.8375   -0.4375    0.0000 R#  0  0  0  0  0  0  0  0  0  5  0  0
    2.5520   -0.0250    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    3.2664   -0.4375    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    2.5520    0.8000    0.0000 O   0  0  0  0  0  0  0  0  0  7  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  2  0  0  0  0
  3 F    2   8  17
V    2 carboxyliccaid,acidchloride
M  ALS   3  2 F O   Cl\x20\x20
M  RGP  1   1   2
M  END
$MOL

  Marvin  07140909202D\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20

 10 10  0  0  0  0            999 V2000
   10.5791   -0.8125    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0
   11.4041   -0.8125    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   11.7068    0.0258    0.0000 N   0  0  0  0  0  0  0  0  0  2  0  0
   10.9916    0.4583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3266   -0.0284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5419    0.2262    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.3700    1.0331    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    8.5852    1.2877    0.0000 O   0  0  0  0  0  0  0  0  0  7  0  0
    9.9828    1.5854    0.0000 R#  0  0  0  0  0  0  0  0  0  5  0  0
   12.4459    0.2254    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  2  0  0  0  0
  5  6  1  0  0  0  0
  2  3  1  0  0  0  0
  6  7  1  0  0  0  0
  3  4  1  0  0  0  0
  7  8  2  0  0  0  0
  4  5  2  0  0  0  0
  7  9  1  0  0  0  0
  5  1  1  0  0  0  0
  3 10  1  0  0  0  0
M  RGP  2   9   2  10   1
M  END
";

fn boronic() -> ChemicalReaction {
    let mut rxn = reaction_from_rxn_file(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/boronic1.rxn"
    ))
    .unwrap();
    rxn.initialize().unwrap();
    rxn
}

fn block(text: &str) -> ChemicalReaction {
    let mut rxn = reaction_from_rxn_block(text).unwrap();
    rxn.initialize().unwrap();
    rxn
}

fn labels(
    result: &rdkit_cxx::chem::rd_chem_reactions::PreprocessReactionResult,
) -> Vec<Vec<(u32, String)>> {
    result.4.clone()
}

#[test]
fn preprocessing() -> Result<()> {
    let mut rxn = boronic();
    let result = preprocess_reaction(&mut rxn)?;
    assert_eq!((result.0, result.1, result.2, result.3), (0, 0, 2, 1));
    assert_eq!(
        labels(&result),
        [
            vec![(0, "halogen.bromine.aromatic".to_owned())],
            vec![(1, "boronicacid".to_owned())]
        ]
    );
    // Labeled atoms match only their groups afterward.
    let mut rxn = boronic();
    let (m1, m2) = (mol_from_smiles("CCBr")?, mol_from_smiles("c1ccccc1Br")?);
    {
        let r1 = rxn.get_reactant_template(0)?;
        assert!(m1.has_substruct_match(&r1)?);
        assert!(m2.has_substruct_match(&r1)?);
    }
    preprocess_reaction(&mut rxn)?;
    let r1 = (*rxn.get_reactant_template(0)?).clone();
    assert!(!m1.has_substruct_match(&r1)?);
    assert!(m2.has_substruct_match(&r1)?);

    // Labels separated by commas are or queries.
    let mut rxn = block(AZIDE_REACTION);
    let result = preprocess_reaction(&mut rxn)?;
    assert_eq!(
        labels(&result),
        [
            vec![(1, "azide".to_owned())],
            vec![(1, "carboxylicacid,acidchloride".to_owned())]
        ]
    );
    let r2 = (*rxn.get_reactant_template(1)?).clone();
    let matches: Vec<bool> = ["CC(=O)O", "CC(=O)Cl", "CC(=O)N"]
        .iter()
        .map(|s| {
            mol_from_smiles(s)
                .unwrap()
                .has_substruct_match(&r2)
                .unwrap()
        })
        .collect();
    assert_eq!(matches, [true, true, false]);

    for (text, key) in [
        (BAD_VALUE1, "boromicacid"),
        (BAD_VALUE2, "carboxylicacid,acidchlroide"),
        (BAD_VALUE3, "carboxyliccaid,acidchloride"),
    ] {
        match preprocess_reaction(&mut block(text)) {
            Err(Error::Key { key: k }) => assert_eq!(k, key),
            other => panic!("{other:?}"),
        }
    }
    let mut empty = ChemicalReaction::default();
    empty.initialize()?;
    let result = preprocess_reaction(&mut empty)?;
    assert_eq!((result.0, result.1, result.2, result.3), (0, 2, 0, 0));
    assert!(result.4.is_empty());
    Ok(())
}

#[test]
fn preprocessing_with_a_file() -> Result<()> {
    let file = common::temp_path("groups.txt");
    std::fs::write(
        &file,
        "halogen.bromine.aromatic\t[Br;$(Br-c)]\nboronicacid\t[$(B-!@[#6])](O)(O)\n",
    )
    .unwrap();
    let mut rxn = boronic();
    let params = PreprocessReactionParams {
        func_group_filename: Some(file.clone()),
        ..Default::default()
    };
    let result = preprocess_reaction_with_params(&mut rxn, &params)?;
    assert_eq!((result.0, result.1, result.2, result.3), (0, 0, 2, 1));
    let r1 = (*rxn.get_reactant_template(0)?).clone();
    assert!(!mol_from_smiles("CCBr")?.has_substruct_match(&r1)?);
    assert!(mol_from_smiles("c1ccccc1Br")?.has_substruct_match(&r1)?);
    std::fs::remove_file(&file).unwrap();
    let params = PreprocessReactionParams {
        func_group_filename: Some(file),
        ..Default::default()
    };
    assert!(matches!(
        preprocess_reaction_with_params(&mut boronic(), &params),
        Err(Error::Os { .. })
    ));
    Ok(())
}

#[test]
fn enumeration() -> Result<()> {
    let reacts1: Vec<Mol> = ["Brc1ccccc1", "Brc1ncccc1", "Brc1cnccc1"]
        .iter()
        .map(|s| mol_from_smiles(s).unwrap())
        .collect();
    let reacts2: Vec<Mol> = ["CCB(O)O", "CCCB(O)O"]
        .iter()
        .map(|s| mol_from_smiles(s).unwrap())
        .collect();
    let bbs = [reacts1.clone(), reacts2];
    let first = |prods: Vec<Result<Vec<Mol>>>| -> Vec<String> {
        prods
            .into_iter()
            .map(|p| mol_to_smiles(&p.unwrap()[0]).unwrap())
            .collect()
    };
    let mut rxn = boronic();
    let prods = first(enumerate_reaction(&mut rxn, &bbs)?.collect());
    assert_eq!(
        prods,
        [
            "CCc1ccccc1",
            "CCc1ccccc1",
            "CCCc1ccccc1",
            "CCCc1ccccc1",
            "CCc1ccccn1",
            "CCc1ccccn1",
            "CCCc1ccccn1",
            "CCCc1ccccn1",
            "CCc1cccnc1",
            "CCc1cccnc1",
            "CCCc1cccnc1",
            "CCCc1cccnc1"
        ]
    );
    let mut rxn = boronic();
    let prods =
        first(enumerate_reaction_with_unique_products_only(&mut rxn, &bbs, true)?.collect());
    assert_eq!(
        prods,
        [
            "CCc1ccccc1",
            "CCCc1ccccc1",
            "CCc1ccccn1",
            "CCCc1ccccn1",
            "CCc1cccnc1",
            "CCCc1cccnc1"
        ]
    );
    let mut rxn = boronic();
    match enumerate_reaction(&mut rxn, &bbs[..1]) {
        Err(Error::Value { message }) => {
            assert_eq!(message, "2 reactants in reaction, 1 bb lists supplied")
        }
        Err(e) => panic!("{e:?}"),
        Ok(_) => panic!("expected an error"),
    }
    let mut empty = ChemicalReaction::default();
    match enumerate_reaction(&mut empty, &[]) {
        Err(Error::Value { message }) => assert_eq!(message, "bad reaction"),
        Err(e) => panic!("{e:?}"),
        Ok(_) => panic!("expected an error"),
    }
    Ok(())
}
