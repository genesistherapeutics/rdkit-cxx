// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// AllChem's helpers, with its doctests; values from Python's AllChem in
// RDKit 2026.03.6, run on these same mol blocks. Alignment and
// optimization are compared with tolerances, since their floating point can
// differ across platforms.

use rdkit_cxx::chem::all_chem::*;
use rdkit_cxx::chem::{
    Mol, MolFromMolBlockParams, mol_from_mol_block_with_params, mol_from_mol_file,
    mol_from_pdb_file, mol_from_smiles, mol_to_smiles,
};
use rdkit_cxx::{Error, Result};

const METHANE: &str = "
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0215    0.0033    0.0161 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6655    0.8882   -0.1076 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3797   -0.8611   -0.5824 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0951   -0.3071    1.0660 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9716    0.2767   -0.3921 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
";
const WATER: &str = "
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -0.0057    0.3844   -0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7970   -0.1943   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8026   -0.1901    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
";
const ETHANOL: &str = "
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
";
const PENTANOL_0: &str = "
     RDKit          3D

 15 14  0  0  0  0  0  0  0  0999 V2000
    1.5320    0.4531   -0.3775 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9666   -0.6991    0.3710 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4711   -0.8808   -0.0770 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1712    0.3911    0.2513 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5135    0.4178   -0.0881 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2243    0.5179   -1.4473 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.3032    1.3956    0.1602 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6555    0.3690   -0.3257 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5426   -1.6368    0.2846 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9407   -0.4198    1.4488 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5228   -1.1819   -1.1319 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9466   -1.7009    0.5195 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1208    0.5346    1.3518 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7050    1.3065   -0.1880 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7140    1.1336   -0.7517 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
  2  9  1  0
  2 10  1  0
  3 11  1  0
  3 12  1  0
  4 13  1  0
  4 14  1  0
  5 15  1  0
M  END
";
const PENTANOL_1: &str = "
     RDKit          3D

 15 14  0  0  0  0  0  0  0  0999 V2000
   -1.6585    0.5996   -0.0068 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7432   -0.3581    0.7402 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3672   -0.8688   -0.1278 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2571    0.2575   -0.6619 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8347    0.9523    0.4050 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8617    0.1122   -0.9902 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5961    0.6398    0.6114 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100    1.5911   -0.1534 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3807   -1.2167    1.0425 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2930    0.1656    1.6062 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1003   -1.4212   -0.9826 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9948   -1.5446    0.4544 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0580   -0.2590   -1.2513 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6795    0.9263   -1.3424 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6525    0.4241    0.6566 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
  2  9  1  0
  2 10  1  0
  3 11  1  0
  3 12  1  0
  4 13  1  0
  4 14  1  0
  5 15  1  0
M  END
";
const PENTANOL_2: &str = "
     RDKit          3D

 15 14  0  0  0  0  0  0  0  0999 V2000
    1.5003    0.4795    0.2130 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8391   -0.8590   -0.0063 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5088   -0.7019   -0.6547 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4838    0.1132    0.1353 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9284    1.3816    0.2966 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2160    0.9903    1.1412 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6071    0.3270    0.1695 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1994    1.1398   -0.6417 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5129   -1.3896   -0.7393 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7939   -1.4897    0.8993 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3139   -0.1488   -1.6144 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8915   -1.7067   -0.9148 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7633   -0.3322    1.1144 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4450    0.2862   -0.4208 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3340    1.9103    1.0227 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
  2  9  1  0
  2 10  1  0
  3 11  1  0
  3 12  1  0
  4 13  1  0
  4 14  1  0
  5 15  1  0
M  END
";
const PENTANOL_3: &str = "
     RDKit          3D

 15 14  0  0  0  0  0  0  0  0999 V2000
   -2.1275    0.0498   -0.1788 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7452    0.3327   -0.6601 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3287   -0.5719   -0.1046 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6976   -0.1860   -0.6685 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0574    1.1142   -0.3608 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1831   -0.2704    0.8697 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6996   -0.6419   -0.8161 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6788    1.0319   -0.2084 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7541    0.1332   -1.7711 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4488    1.3823   -0.5165 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3328   -0.6778    0.9795 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.1270   -1.5998   -0.5200 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4465   -0.9086   -0.2684 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6587   -0.2872   -1.7675 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9884    1.0996   -0.0251 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
  2  9  1  0
  2 10  1  0
  3 11  1  0
  3 12  1  0
  4 13  1  0
  4 14  1  0
  5 15  1  0
M  END
";
const TEMPLATE: &str = "
     RDKit          3D

 12 13  0  0  0  0  0  0  0  0999 V2000
   -2.9785    0.0780    1.2005 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9897   -0.7630    0.8460 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5520   -0.3866   -0.3400 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4555   -1.0327   -1.0430 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8492   -0.4756   -0.5854 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5105   -1.0478    0.5070 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7223   -0.5100    0.9515 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2712    0.6035    0.3078 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6075    1.1806   -0.7793 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3957    0.6431   -1.2242 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2235    0.6920   -0.7718 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1571    1.0185    0.1959 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
  3 11  1  0
 11 12  2  0
 12  1  1  0
 10  5  1  0
M  END
";

fn read(block: &str) -> Mol {
    mol_from_mol_block_with_params(
        block,
        &MolFromMolBlockParams {
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap()
}

/// Pentanol with its four conformers, ids 0 to 3.
fn pentanol() -> Mol {
    let mut m = read(PENTANOL_0);
    for block in [PENTANOL_1, PENTANOL_2, PENTANOL_3] {
        let conf = Conformer::from(read(block).get_conformer().unwrap());
        m.add_conformer_with_assign_id(&conf, true).unwrap();
    }
    m
}

fn assert_close(actual: &[f64], expected: &[f64], tol: f64) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() <= tol, "{actual:?} != {expected:?}");
    }
}

#[test]
fn volumes() -> Result<()> {
    let wide = ComputeMolVolumeParams {
        grid_spacing: 0.3,
        box_margin: 1.0,
        ..Default::default()
    };
    for (block, volume, wide_volume) in [
        (METHANE, 28.60000000000001, 28.376999999999995),
        (WATER, 19.928000000000004, 20.060999999999996),
        (ETHANOL, 54.064000000000014, 54.13499999999999),
    ] {
        let m = read(block);
        assert_close(&[compute_mol_volume(&m)?], &[volume], 1e-9);
        assert_close(
            &[compute_mol_volume_with_params(&m, &wide)?],
            &[wide_volume],
            1e-9,
        );
    }
    Ok(())
}

#[test]
fn shapes() -> Result<()> {
    let m = read(ETHANOL);
    let occupancy = |grid: &rdkit_cxx::geometry::rd_geometry::UniformGrid3D| {
        let v = grid.get_occupancy_vect();
        (0..v.len())
            .map(|i| v.get_val(i).unwrap())
            .collect::<Vec<u32>>()
    };
    let grid = compute_mol_shape(&m)?;
    assert_eq!(grid.get_num_x(), 40);
    let occ = occupancy(&grid);
    assert_eq!(occ.iter().sum::<u32>(), 917);
    assert_eq!(occ.iter().filter(|&&x| x == 3).count(), 189);
    let params = ComputeMolShapeParams {
        box_dim: (10.0, 12.0, 14.0),
        spacing: 0.4,
        encode_shape: EncodeShapeParams {
            vdw_scale: 1.0,
            ..Default::default()
        },
        ..Default::default()
    };
    let grid = compute_mol_shape_with_params(&m, &params)?;
    assert_eq!(
        (grid.get_num_x(), grid.get_num_y(), grid.get_num_z()),
        (25, 30, 35)
    );
    assert_eq!(occupancy(&grid).iter().sum::<u32>(), 2870);
    Ok(())
}

#[test]
fn conformer_rms() -> Result<()> {
    let mut m = pentanol();
    assert_close(
        &[get_conformer_rms(&mut m, 0, 1)?],
        &[1.5335098389207877],
        1e-9,
    );
    // The second conformer is left aligned to the first.
    let p = m.get_conformer_with_id(1)?.get_atom_position(0)?;
    assert_close(
        &[p.x, p.y, p.z],
        &[1.6818476626902419, 0.44606895996556667, -0.2874533967509153],
        1e-9,
    );
    let heavy = GetConformerRMSParams {
        atom_ids: Some(vec![0, 1, 2, 3, 4]),
        ..Default::default()
    };
    assert_close(
        &[get_conformer_rms_with_params(
            &mut pentanol(),
            0,
            2,
            &heavy,
        )?],
        &[1.2683239757921558],
        1e-9,
    );
    let prealigned = GetConformerRMSParams {
        prealigned: true,
        ..Default::default()
    };
    assert_close(
        &[get_conformer_rms_with_params(
            &mut pentanol(),
            0,
            1,
            &prealigned,
        )?],
        &[3.164152573965632],
        1e-9,
    );
    assert_close(
        &get_conformer_rms_matrix(&mut pentanol())?,
        &[
            1.5335098389207868,
            1.18516123453815,
            1.5691895558175155,
            1.5020222097350482,
            1.7443642365849272,
            1.5629532691344719,
        ],
        1e-9,
    );
    assert_close(
        &get_conformer_rms_matrix_with_params(&mut pentanol(), &heavy)?,
        &[
            0.626681981621143,
            0.6879966192060447,
            1.5582671927131386,
            0.7425089977040414,
            1.6944292540855466,
            1.6367338381581584,
        ],
        1e-9,
    );
    assert_close(
        &get_conformer_rms_matrix_with_params(&mut pentanol(), &prealigned)?,
        &[
            3.164152573965632,
            1.510884064380851,
            3.181272367151232,
            3.7269144905672307,
            1.6659053696613944,
            3.521553038835754,
        ],
        1e-9,
    );
    Ok(())
}

#[test]
fn transforms() -> Result<()> {
    let tform = vec![
        vec![1.0, 0.0, 0.0, 1.5],
        vec![0.0, 1.0, 0.0, 0.0],
        vec![0.0, 0.0, 1.0, -2.0],
        vec![0.0, 0.0, 0.0, 1.0],
    ];
    let position = |m: &Mol, id: i32| {
        let p = m
            .get_conformer_with_id(id)
            .unwrap()
            .get_atom_position(0)
            .unwrap();
        vec![p.x, p.y, p.z]
    };
    let ids = |m: &Mol| m.get_conformers().map(|c| c.get_id()).collect::<Vec<u32>>();
    let mut m = pentanol();
    transform_mol(&mut m, &tform)?;
    assert_eq!(ids(&m), [0]);
    assert_close(&position(&m, 0), &[3.032, 0.4531, -2.3775], 1e-12);
    let mut m = pentanol();
    let params = TransformMolParams {
        conf_id: 2,
        keep_confs: true,
    };
    transform_mol_with_params(&mut m, &tform, &params)?;
    assert_eq!(m.get_num_conformers(), 4);
    assert_close(&position(&m, 2), &[3.0003, 0.4795, -1.787], 1e-12);
    let mut m = pentanol();
    let params = TransformMolParams {
        conf_id: 2,
        keep_confs: false,
    };
    transform_mol_with_params(&mut m, &tform, &params)?;
    assert_eq!(ids(&m), [0]);
    assert_close(&position(&m, 0), &[3.0003, 0.4795, -1.787], 1e-12);
    Ok(())
}

#[test]
fn libraries_from_reactions() -> Result<()> {
    let s1 = vec![mol_from_smiles("NC")?, mol_from_smiles("NCC")?];
    let s2 = vec![mol_from_smiles("OC=O")?, mol_from_smiles("OC(=O)C")?];
    let mut rxn = reaction_from_smarts("[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]")?;
    let sets = [s2, s1];
    let products: Vec<String> = enumerate_library_from_reaction(&mut rxn, &sets)?
        .map(|p| mol_to_smiles(&p.unwrap()[0]).unwrap())
        .collect();
    assert_eq!(products, ["CNC=O", "CCNC=O", "CNC(C)=O", "CCNC(C)=O"]);
    let reactants: Vec<Vec<String>> =
        enumerate_library_from_reaction_return_reactants(&mut rxn, &sets)?
            .map(|pr| {
                let pr = pr.unwrap();
                assert_eq!(pr.products.len(), 1);
                pr.reactants
                    .iter()
                    .map(|m| mol_to_smiles(m).unwrap())
                    .collect()
            })
            .collect();
    assert_eq!(
        reactants,
        [
            ["O=CO", "CN"],
            ["O=CO", "CCN"],
            ["CC(=O)O", "CN"],
            ["CC(=O)O", "CCN"]
        ]
    );
    match enumerate_library_from_reaction(&mut rxn, &sets[..1]) {
        Err(Error::Value { message }) => assert_eq!(message, "1 sidechains provided, 2 required"),
        _ => panic!("expected an error"),
    }
    // The doctest's lazy library: the first products of a large one.
    let amines: Vec<Mol> = (0..100)
        .map(|x| mol_from_smiles(&format!("N{}", "C".repeat(x))).unwrap())
        .collect();
    let acids: Vec<Mol> = (0..100)
        .map(|x| mol_from_smiles(&format!("OC(=O){}", "C".repeat(x))).unwrap())
        .collect();
    let sets = [acids, amines];
    let first: Vec<String> = enumerate_library_from_reaction(&mut rxn, &sets)?
        .take(4)
        .map(|p| mol_to_smiles(&p.unwrap()[0]).unwrap())
        .collect();
    assert_eq!(first, ["NC=O", "CNC=O", "CCNC=O", "CCCNC=O"]);
    Ok(())
}

#[test]
fn constrained_embedding() -> Result<()> {
    let template = read(TEMPLATE);
    for (use_tethers, rms, deltas) in [
        (
            true,
            0.0063113184357433145,
            [
                [
                    -0.004940880371335332,
                    -0.0008436319453429819,
                    0.009285559242237618,
                ],
                [
                    -0.003063817118147183,
                    0.0027166997921834746,
                    0.0006283794415891109,
                ],
            ],
        ),
        (
            false,
            1.3337906193165248,
            [
                [
                    -0.056376171061598424,
                    -1.195030905901694,
                    -1.0579735813140485,
                ],
                [
                    0.0028233995731397066,
                    0.048180647262614285,
                    0.04040647901225547,
                ],
            ],
        ),
    ] {
        let mut mol = mol_from_smiles("c1nn(Cc2ccccc2)cc1-c3ccccc3")?;
        let params = ConstrainedEmbedParams {
            use_tethers,
            ..Default::default()
        };
        constrained_embed_with_params(&mut mol, &template, &params)?;
        let embed_rms: String = mol.get_prop("EmbedRMS")?;
        assert_close(&[embed_rms.parse().unwrap()], &[rms], 1e-3);
        for (i, delta) in deltas.iter().enumerate() {
            let p = mol.get_conformer()?.get_atom_position(i as u32)?;
            let t = template.get_conformer()?.get_atom_position(i as u32)?;
            assert_close(&[p.x - t.x, p.y - t.y, p.z - t.z], delta, 1e-3);
        }
    }
    // The doctest's default form.
    let mut mol = mol_from_smiles("c1nn(Cc2ccccc2)cc1-c3ccccc3")?;
    constrained_embed(&mut mol, &template)?;
    let p = mol.get_conformer()?.get_atom_position(0)?;
    let t = template.get_conformer()?.get_atom_position(0)?;
    assert!(p.distance(&t) < 0.02);
    match constrained_embed(&mut mol_from_smiles("CCCC")?, &template) {
        Err(Error::Value { message }) => assert_eq!(message, "molecule doesn't match the core"),
        other => panic!("{other:?}"),
    }
    Ok(())
}

#[test]
fn bond_orders_from_templates() -> Result<()> {
    let singles = |m: &Mol| {
        m.get_bonds()
            .filter(|b| b.get_bond_type_as_double().unwrap() == 1.0)
            .count()
    };
    let template = mol_from_smiles("CN1C(=NC(C1=O)(c2ccccc2)c3ccccc3)N")?;
    let mol = mol_from_pdb_file(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/4DJU_lig.pdb"
    ))?;
    assert_eq!((singles(&template), singles(&mol)), (8, 22));
    let new_mol = assign_bond_orders_from_template(&template, &mol)?;
    assert_eq!(singles(&new_mol), 8);
    assert_eq!(
        mol_to_smiles(&new_mol)?,
        "CN1C(=O)C(c2ccccc2)(c2ccccc2)N=C1N"
    );
    // GitHub 235: different formal charges.
    let template =
        mol_from_smiles("CN(C)C(=O)Cc1ccc2c(c1)NC(=O)c3ccc(cc3N2)c4ccc(c(c4)OC)[N+](=O)[O-]")?;
    let mol = mol_from_mol_file(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/4FTR_lig.mol"
    ))?;
    assert_eq!(
        mol_to_smiles(&mol)?,
        "COC1CC(C2CCC3C(O)NC4CC(CC(O)N(C)C)CCC4NC3C2)CCC1N(O)O"
    );
    assert_eq!(
        mol_to_smiles(&assign_bond_orders_from_template(&template, &mol)?)?,
        "COc1cc(-c2ccc3c(c2)Nc2ccc(CC(=O)N(C)C)cc2NC3=O)ccc1[N+](=O)[O-]"
    );
    match assign_bond_orders_from_template(&mol_from_smiles("c1ccccc1")?, &mol_from_smiles("CCO")?)
    {
        Err(Error::Value { message }) => assert_eq!(message, "No matching found"),
        other => panic!("{other:?}"),
    }
    let same = assign_bond_orders_from_template(&mol_from_smiles("CC")?, &mol_from_smiles("CCO")?)?;
    assert_eq!(mol_to_smiles(&same)?, "CCO");
    Ok(())
}
