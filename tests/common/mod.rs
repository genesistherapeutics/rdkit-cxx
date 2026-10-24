// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Helpers shared by the integration tests.

#![allow(dead_code)] // each test file uses some of these

use std::path::PathBuf;

use rdkit_cxx::chem::{mol_from_smiles, mol_to_mol_block};

/// A path in the temporary directory that no other test process uses.
pub fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rdkit-cxx-{}-{name}", std::process::id()))
}

/// Writes a mol file for methanol whose title line is Latin-1, not UTF-8.
pub fn latin1_mol_file() -> PathBuf {
    let block = mol_to_mol_block(&mol_from_smiles("CO").unwrap()).unwrap();
    let rest = block.split_once('\n').unwrap().1;
    let path = temp_path("latin1.mol");
    std::fs::write(&path, [b"\xe9t\xe9\n".as_slice(), rest.as_bytes()].concat()).unwrap();
    path
}

/// Ethanol with hydrogens and one conformer, as RDKit's embedder placed it
/// with seed 42, rounded to the mol block's four decimals. The Python
/// scripts that give the tests their expected values read the same block.
pub const ETHANOL: &str = "
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

/// A second ethanol conformer, the embedder's next with seed 42.
pub const ETHANOL_CONF2: [(f64, f64, f64); 9] = [
    (-0.9469, 0.2047, 0.0073),
    (0.4345, -0.3899, -0.1941),
    (1.343, 0.0493, 0.7831),
    (-1.5855, -0.1732, -0.7986),
    (-1.2972, -0.2011, 0.9838),
    (-0.9359, 1.3064, 0.0392),
    (0.2905, -1.4966, -0.0905),
    (0.8644, -0.158, -1.1864),
    (1.8331, 0.8584, 0.4561),
];

/// `ETHANOL`, keeping its hydrogens.
pub fn ethanol() -> rdkit_cxx::chem::Mol {
    use rdkit_cxx::chem::{MolFromMolBlockParams, mol_from_mol_block_with_params};
    let params = MolFromMolBlockParams {
        remove_hs: false,
        ..Default::default()
    };
    mol_from_mol_block_with_params(ETHANOL, &params).unwrap()
}

/// `ethanol()` with `ETHANOL_CONF2` as a second conformer.
pub fn ethanol_two_confs() -> rdkit_cxx::chem::Mol {
    use rdkit_cxx::chem::Conformer;
    use rdkit_cxx::geometry::Point3D;
    let mut mol = ethanol();
    let mut conf = Conformer::new(9);
    let positions: Vec<Point3D> = ETHANOL_CONF2
        .iter()
        .map(|&(x, y, z)| Point3D::new(x, y, z))
        .collect();
    conf.set_positions(&positions).unwrap();
    conf.set_3d(true);
    mol.add_conformer_with_assign_id(&conf, true).unwrap();
    mol
}

/// Asserts that two floats agree to 1e-9, relative to the larger.
pub fn assert_close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= 1e-9 * scale,
        "{actual} != {expected}"
    );
}
