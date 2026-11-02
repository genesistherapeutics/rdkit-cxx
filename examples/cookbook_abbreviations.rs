// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! With Abbreviations (RDKitCB_34): draw molecules with functional groups
//! abbreviated.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#with-abbreviations),
//! by Greg Landrum.

mod common;

use common::{MOL_SIZE, show, write_image};
use rdkit_cxx::chem::draw::rd_mol_draw_2d::MolDrawOptions;
use rdkit_cxx::chem::draw::{MolsToGridImageParams, mols_to_grid_image_use_svg_with_params};
use rdkit_cxx::chem::rd_abbreviations::{
    CondenseMolAbbreviationsParams, condense_mol_abbreviations,
    condense_mol_abbreviations_with_params, get_default_abbreviations,
};
use rdkit_cxx::chem::{mol_from_smiles, mol_to_smiles};

/// The table the recipe prints.
const ABBREVIATIONS: &str = "\
Abbrev     SMILES
--------   --------
CO2Et      *C(=O)OCC
COOEt      *C(=O)OCC
OiBu       *OCC(C)C
nDec       *CCCCCCCCCC
nNon       *CCCCCCCCC
nOct       *CCCCCCCC
nHept      *CCCCCCC
nHex       *CCCCCC
nPent      *CCCCC
iPent      *C(C)CCC
tBu        *C(C)(C)C
iBu        *C(C)CC
nBu        *CCCC
iPr        *C(C)C
nPr        *CCC
Et         *CC
NCF3       *NC(F)(F)F
CF3        *C(F)(F)F
CCl3       *C(Cl)(Cl)Cl
CN         *C#N
NC         *[N+]#[C-]
N(OH)CH3   *N(C)[OH]
NO2        *[N+](=O)[O-]
NO         *N=O
SO3H       *S(=O)(=O)[OH]
CO2H       *C(=O)[OH]
COOH       *C(=O)[OH]
OEt        *OCC
OAc        *OC(C)=O
NHAc       *NC(C)=O
Ac         *C(C)=O
CHO        *C=O
NMe        *NC
SMe        *SC
OMe        *OC
CO2-       *C(=O)[O-]
COO-       *C(=O)[O-]
";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let console = MolDrawOptions::default();

    let m = mol_from_smiles("COc1ccc(C(=O)[O-])cc1")?;
    write_image("abbreviations_plain.svg", show(&m, MOL_SIZE, &console)?)?;

    let abbrevs = get_default_abbreviations();
    let nm = condense_mol_abbreviations(&m, &abbrevs)?;
    write_image("abbreviations.svg", show(&nm, MOL_SIZE, &console)?)?;
    let labels: Vec<String> = nm
        .get_atoms()
        .filter(|atom| atom.has_prop("atomLabel"))
        .map(|atom| atom.get_prop("atomLabel"))
        .collect::<rdkit_cxx::Result<_>>()?;
    assert_eq!(labels, ["OMe", "CO2-"]);
    assert_eq!((m.get_num_atoms(), nm.get_num_atoms()), (11, 8));

    // Abbreviations that cover more than 40% of the molecule are not
    // applied by default.
    let m = mol_from_smiles("c1c[nH]cc1C(F)(F)F")?;
    let nm1 = condense_mol_abbreviations(&m, &abbrevs)?;
    let nm2 = condense_mol_abbreviations_with_params(
        &m,
        &abbrevs,
        &CondenseMolAbbreviationsParams {
            max_coverage: 0.8,
            ..Default::default()
        },
    )?;
    assert_eq!(
        [m.get_num_atoms(), nm1.get_num_atoms(), nm2.get_num_atoms()],
        [9, 9, 6]
    );
    let grid = MolsToGridImageParams {
        legends: Some(
            ["", "default", "maxCoverage=0.8"]
                .map(str::to_owned)
                .to_vec(),
        ),
        draw_options: Some(console),
        ..Default::default()
    };
    let svg = mols_to_grid_image_use_svg_with_params(&[m, nm1, nm2], &grid)?;
    write_image("abbreviations_coverage.svg", &svg)?;

    // The available abbreviations and their SMILES, where * is the dummy
    // atom the group attaches to.
    let mut table = format!("{:<10} {}\n", "Abbrev", "SMILES");
    table += &format!("{:<10} {}\n", "--------", "--------");
    for a in &abbrevs {
        let mol = a.mol.as_ref().expect("each abbreviation has a molecule");
        table += &format!("{:<10} {}\n", a.label, mol_to_smiles(mol)?);
    }
    print!("{table}");
    assert_eq!(table, ABBREVIATIONS);
    Ok(())
}
