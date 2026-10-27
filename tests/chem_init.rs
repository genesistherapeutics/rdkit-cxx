// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The Python functions of rdkit/Chem/__init__.py (FindMolChiralCenters has
// its own tests); values from Python's rdkit.Chem in RDKit 2026.03.6.

mod common;

use rdkit_cxx::chem::{
    FileMolSupplier, SDMolSupplierParams, canon_smiles, canon_smiles_with_use_chiral,
    mol_to_smiles, multi_conf_mol_from_sdf, multi_conf_mol_from_sdf_with_params,
    quick_smarts_match, quick_smarts_match_with_unique, supplier_from_filename,
};
use rdkit_cxx::{Error, Result};

/// Three conformers of ethanol with hydrogens, from RDKit's embedder with
/// seed 42, as `SDWriter` wrote them.
const ETHANOL_CONFS: &str = "
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.9640    0.0151    0.0410 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5041   -0.3085   -0.2036 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2237    0.3416    0.7792 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3563    0.7407   -0.7043 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5833   -0.9142   -0.0029 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0287    0.5174    1.0296 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7725   -0.0064   -1.2300 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5980   -1.4177   -0.1686 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8339    1.0320    0.4595 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
$$$$

     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.9469    0.2047    0.0073 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4345   -0.3899   -0.1941 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3430    0.0493    0.7831 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5855   -0.1732   -0.7986 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2972   -0.2011    0.9838 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9359    1.3064    0.0392 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2905   -1.4966   -0.0905 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8644   -0.1580   -1.1864 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8331    0.8584    0.4561 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
$$$$

     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8949    0.2106   -0.0171 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4408   -0.4375    0.1828 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3958    0.1053   -0.6585 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8903    1.3078    0.1428 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1585    0.0656   -1.0949 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7004   -0.2953    0.5556 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6964   -0.3666    1.2514 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3134   -1.5297   -0.0275 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7978    0.9398   -0.3346 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
$$$$
";

#[test]
fn quick_smarts_matches() -> Result<()> {
    assert_eq!(quick_smarts_match("c1ccccc1O", "[#6][#8]")?, [[5, 6]]);
    assert_eq!(
        quick_smarts_match_with_unique("OCCO", "[#8][#6]", false)?,
        [[0, 1], [3, 2]]
    );
    assert!(matches!(
        quick_smarts_match("C1CC", "C"),
        Err(Error::Parse { .. })
    ));
    Ok(())
}

#[test]
fn canonical_smiles() -> Result<()> {
    assert_eq!(canon_smiles("OC[C@@H](F)Cl")?, "OC[C@@H](F)Cl");
    assert_eq!(
        canon_smiles_with_use_chiral("OC[C@@H](F)Cl", false)?,
        "OCC(F)Cl"
    );
    assert_eq!(canon_smiles_with_use_chiral("[13CH3]O", false)?, "CO");
    Ok(())
}

#[test]
fn suppliers_from_filenames() -> Result<()> {
    let csv = common::temp_path("supplier.csv");
    std::fs::write(&csv, "smiles,name\nCCO,ethanol\nc1ccccc1,benzene\n").unwrap();
    let mut suppl = supplier_from_filename(&csv)?;
    assert!(matches!(suppl, FileMolSupplier::Smiles(_)));
    let mols: Vec<_> = suppl.mols().collect::<Result<_>>()?;
    let smiles: Vec<String> = mols.iter().map(mol_to_smiles).collect::<Result<_>>()?;
    assert_eq!(smiles, ["CCO", "c1ccccc1"]);
    let names: Vec<String> = mols
        .iter()
        .map(|m| m.get_prop("_Name"))
        .collect::<Result<_>>()?;
    assert_eq!(names, ["ethanol", "benzene"]);

    let sdf = common::temp_path("supplier.SDF");
    std::fs::write(&sdf, ETHANOL_CONFS).unwrap();
    let mut suppl = supplier_from_filename(&sdf)?;
    assert!(matches!(suppl, FileMolSupplier::Sd(_)));
    assert_eq!(suppl.mols().count(), 3);

    // Python's fails for every .tdt file; this one reads it.
    let tdt = common::temp_path("supplier.tdt");
    std::fs::write(&tdt, "$SMI<CCO>\n|\n$SMI<c1ccccc1>\n|\n").unwrap();
    let mut suppl = supplier_from_filename(&tdt)?;
    assert!(matches!(suppl, FileMolSupplier::Tdt(_)));
    assert_eq!(suppl.mols().count(), 2);

    match supplier_from_filename(common::temp_path("supplier.foo")) {
        Err(Error::Value { message }) => assert_eq!(message, "unrecognized extension: foo"),
        other => panic!("{other:?}"),
    }
    for path in [csv, sdf, tdt] {
        std::fs::remove_file(path).unwrap();
    }
    Ok(())
}

#[test]
fn multi_conformer_molecules() -> Result<()> {
    let sdf = common::temp_path("multi_conf.sdf");
    std::fs::write(&sdf, ETHANOL_CONFS).unwrap();
    let mol = multi_conf_mol_from_sdf(&sdf)?.unwrap();
    assert_eq!(mol.get_num_atoms(), 3);
    let ids: Vec<u32> = mol.get_conformers().map(|c| c.get_id()).collect();
    assert_eq!(ids, [0, 1, 2]);
    let p = mol.get_conformer_with_id(2)?.get_atom_position(0)?;
    assert_eq!((p.x, p.y, p.z), (-0.8949, 0.2106, -0.0171));
    let params = SDMolSupplierParams {
        remove_hs: false,
        ..Default::default()
    };
    let mol = multi_conf_mol_from_sdf_with_params(&sdf, &params)?.unwrap();
    assert_eq!((mol.get_num_atoms(), mol.get_num_conformers()), (9, 3));
    std::fs::remove_file(sdf).unwrap();
    Ok(())
}
