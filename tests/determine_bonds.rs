// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use common::ethanol;
use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_determine_bonds::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::Point3D;

// Expected values below come from Python's rdDetermineBonds on the same
// molecules.

/// Acrylate's atoms, with coordinates from ETKDG and no bonds.
const ACRYLATE: &str = "
     RDKit          3D

  8  0  0  0  0  0  0  0  0  0999 V2000
   -1.2358    0.0575    0.1045 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0440   -0.5343    0.0873 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1617    0.2599   -0.1471 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0960    1.5116   -0.3360 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4330   -0.3171   -0.1739 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3611    1.1228   -0.0500 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1434   -0.5068    0.2756 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0936   -1.5936    0.2396 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";

/// A methyl radical's atoms.
const METHYL: &str = "
     RDKit          3D

  4  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0900    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5450    0.9440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5450   -0.9440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";

/// A hydrogen between two carbons, close enough to bond to both.
const BRIDGE: &str = "
     RDKit          3D

  3  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";

fn atoms(mol_block: &str) -> Mol {
    let params = MolFromMolBlockParams {
        remove_hs: false,
        ..Default::default()
    };
    mol_from_mol_block_with_params(mol_block, &params).unwrap()
}

fn bonds(mol: &Mol) -> Vec<(u32, u32)> {
    let mut pairs: Vec<(u32, u32)> = mol
        .get_bonds()
        .map(|b| {
            let (i, j) = (b.get_begin_atom_idx(), b.get_end_atom_idx());
            (i.min(j), i.max(j))
        })
        .collect();
    pairs.sort();
    pairs
}

fn value<T: std::fmt::Debug>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

#[test]
fn params_defaults() {
    assert_eq!(
        DetermineConnectivityParams::default(),
        DetermineConnectivityParams {
            use_hueckel: false,
            charge: 0,
            cov_factor: 1.3,
            use_vdw: false,
        }
    );
    assert_eq!(
        DetermineBondOrdersParams::default(),
        DetermineBondOrdersParams {
            charge: 0,
            allow_charged_fragments: true,
            embed_chiral: true,
            use_atom_map: false,
            max_iterations: 0,
        }
    );
    assert_eq!(
        DetermineBondsParams::default(),
        DetermineBondsParams {
            use_hueckel: false,
            charge: 0,
            cov_factor: 1.3,
            allow_charged_fragments: true,
            embed_chiral: true,
            use_atom_map: false,
            use_vdw: false,
            max_iterations: 0,
        }
    );
}

#[test]
fn connectivity() {
    let ethanol_bonds = vec![
        (0, 1),
        (0, 3),
        (0, 4),
        (0, 5),
        (1, 2),
        (1, 6),
        (1, 7),
        (2, 8),
    ];
    let mut mol = ethanol();
    determine_connectivity(&mut mol).unwrap();
    assert_eq!(bonds(&mol), ethanol_bonds);
    let vdw = DetermineConnectivityParams {
        use_vdw: true,
        ..Default::default()
    };
    let mut mol = ethanol();
    determine_connectivity_with_params(&mut mol, &vdw).unwrap();
    assert_eq!(bonds(&mol), ethanol_bonds);
    let mut mol = ethanol();
    let tight = DetermineConnectivityParams {
        cov_factor: 1.0,
        ..vdw.clone()
    };
    determine_connectivity_with_params(&mut mol, &tight).unwrap();
    assert_eq!(bonds(&mol), vec![(1, 2)]);

    // The atoms are left without implicit hydrogens.
    let mut mol = atoms(ACRYLATE);
    determine_connectivity(&mut mol).unwrap();
    assert_eq!(
        bonds(&mol),
        vec![(0, 1), (0, 5), (0, 6), (1, 2), (1, 7), (2, 3), (2, 4)]
    );
    assert!(mol.get_atoms().all(|atom| atom.get_no_implicit()));

    // Without a conformer, connect-the-dots adds no bonds, and the van der
    // Waals method has no distances.
    let mut mol = mol_from_smiles("CC").unwrap();
    determine_connectivity(&mut mol).unwrap();
    assert_eq!(mol.get_num_bonds(), 0);
    let mut mol = mol_from_smiles("CC").unwrap();
    assert!(value(determine_connectivity_with_params(&mut mol, &vdw)));
}

#[test]
fn bond_orders() {
    let connected = || {
        let mut mol = atoms(ACRYLATE);
        determine_connectivity(&mut mol).unwrap();
        mol
    };
    let anion = DetermineBondOrdersParams {
        charge: -1,
        ..Default::default()
    };
    let mut mol = connected();
    determine_bond_orders_with_params(&mut mol, &anion).unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "[H]C([H])=C([H])C(=O)[O-]");
    let mapped = DetermineBondOrdersParams {
        use_atom_map: true,
        embed_chiral: false,
        ..anion.clone()
    };
    let mut mol = connected();
    determine_bond_orders_with_params(&mut mol, &mapped).unwrap();
    assert_eq!(
        mol_to_smiles(&mol).unwrap(),
        "[C:1](=[C:2]([C:3](=[O:4])[O-:5])[H:8])([H:6])[H:7]"
    );
    // No ordering gives a neutral molecule; Python raises RDKit's
    // ValueError.
    assert!(value(determine_bond_orders(&mut connected())));
    let radicals = DetermineBondOrdersParams {
        allow_charged_fragments: false,
        ..anion.clone()
    };
    assert!(value(determine_bond_orders_with_params(
        &mut connected(),
        &radicals
    )));

    let methyl = |embed_chiral: bool| {
        let mut mol = atoms(METHYL);
        let params = DetermineBondsParams {
            allow_charged_fragments: false,
            embed_chiral,
            ..Default::default()
        };
        determine_bonds_with_params(&mut mol, &params).unwrap();
        mol
    };
    let mol = methyl(true);
    assert_eq!(mol_to_smiles(&mol).unwrap(), "[H][C]([H])[H]");
    assert_eq!(
        mol.get_atom_with_idx(0)
            .unwrap()
            .get_num_radical_electrons(),
        1
    );
    assert_eq!(mol_to_smiles(&methyl(false)).unwrap(), "[H]C([H])[H]");
}

#[test]
fn bonds_from_coordinates() {
    let anion = DetermineBondsParams {
        charge: -1,
        ..Default::default()
    };
    for params in [
        anion.clone(),
        DetermineBondsParams {
            use_vdw: true,
            ..anion.clone()
        },
    ] {
        let mut mol = atoms(ACRYLATE);
        determine_bonds_with_params(&mut mol, &params).unwrap();
        assert_eq!(mol_to_smiles(&mol).unwrap(), "[H]C([H])=C([H])C(=O)[O-]");
    }
    // Python raises MaxFindBondOrdersItersExceeded, a RuntimeError.
    for max_iterations in [1, 2] {
        let params = DetermineBondsParams {
            max_iterations,
            ..anion.clone()
        };
        assert!(matches!(
            determine_bonds_with_params(&mut atoms(ACRYLATE), &params),
            Err(Error::Runtime { .. })
        ));
    }
    // A molecule of one atom or none is left as it is.
    let mut one = mol_from_smiles("[Rf]").unwrap();
    determine_bonds(&mut one).unwrap();
    determine_bonds(&mut Mol::default()).unwrap();
}

/// Lines of YAeHMOP's built-in parameters for H, C, and O, then N's.
const EHT_PARAMETERS: &str = "\
;EHParameters for the tests
H        1      1      1     1     s  -13.600   1.3000   0.0000   1.0000   0.0000
C        6      4      1     2     s  -21.400   1.6250   0.0000   1.0000   0.0000
C        6      4      1     2     p  -11.400   1.6250   0.0000   1.0000   0.0000
O        8      6      1     2     s  -32.300   2.2750   0.0000   1.0000   0.0000
O        8      6      1     2     p  -14.800   2.2750   0.0000   1.0000   0.0000
N        7      5      1     2     s  -26.000   1.9500   0.0000   1.0000   0.0000
N        7      5      1     2     p  -13.400   1.9500   0.0000   1.0000   0.0000
";

fn hueckel_params(charge: i32) -> DetermineConnectivityParams {
    DetermineConnectivityParams {
        use_hueckel: true,
        charge,
        ..Default::default()
    }
}

#[test]
fn hueckel() {
    if !hueckel_enabled() {
        assert!(value(determine_connectivity_with_params(
            &mut ethanol(),
            &hueckel_params(0)
        )));
        return;
    }
    // Without RDBASE or BIND_PARM_FILE, YAeHMOP reads the parameter file it
    // was built with, and frees memory twice when that is missing, as
    // conda-forge's is: the process aborts.
    if std::env::var_os("RDBASE").is_none() && std::env::var_os("BIND_PARM_FILE").is_none() {
        assert!(value(determine_connectivity_with_params(
            &mut ethanol(),
            &hueckel_params(0)
        )));
    }
    // The rest run in a child process, which has the variables set.
    let good = common::temp_path("eht_parms.dat");
    std::fs::write(&good, EHT_PARAMETERS).unwrap();
    let long = common::temp_path("eht_parms_long.dat");
    std::fs::write(&long, format!(";{}\n{EHT_PARAMETERS}", "x".repeat(300))).unwrap();
    for (mode, file) in [
        ("file", good.as_path()),
        ("long", long.as_path()),
        (
            "missing",
            std::path::Path::new("/nonexistent/eht_parms.dat"),
        ),
    ] {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["hueckel_child", "--exact", "--ignored", "--test-threads=1"])
            .env("RDKIT_CXX_EHT_TEST", mode)
            .env("BIND_PARM_FILE", file)
            .env_remove("RDBASE")
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(
            output.status.success() && stdout.contains("1 passed"),
            "{mode}: {stdout}{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    std::fs::remove_file(good).unwrap();
    std::fs::remove_file(long).unwrap();
}

#[test]
#[ignore = "run by `hueckel` with its environment"]
fn hueckel_child() {
    let Some(mode) = std::env::var_os("RDKIT_CXX_EHT_TEST") else {
        return;
    };
    let hueckel = hueckel_params(0);
    if mode == "long" {
        // YAeHMOP reads each line into 240 bytes.
        assert!(value(determine_connectivity_with_params(
            &mut ethanol(),
            &hueckel
        )));
        return;
    }
    let mut mol = ethanol();
    determine_connectivity_with_params(&mut mol, &hueckel).unwrap();
    assert_eq!(mol.get_num_bonds(), 8);
    // YAeHMOP reads a dummy atom's parameters from a null file, and Python
    // crashes; it has none for elements past lawrencium, and Python exits.
    for atomic_num in [0, 104] {
        let mut mol = atoms(METHYL);
        mol.get_atom_with_idx_mut(1)
            .unwrap()
            .set_atomic_num(atomic_num);
        assert!(value(determine_connectivity_with_params(
            &mut mol, &hueckel
        )));
    }
    if mode == "missing" {
        return; // The built-in parameters, checked as far as they need.
    }

    // The Hückel method sets the first atom's charge.
    let mut mol = atoms(ACRYLATE);
    determine_connectivity_with_params(&mut mol, &hueckel_params(-1)).unwrap();
    assert_eq!(mol.get_num_bonds(), 7);
    assert_eq!(mol.get_atom_with_idx(0).unwrap().get_formal_charge(), -1);
    // Then the bond orders do not add up; Python raises RDKit's ValueError.
    let bonds = DetermineBondsParams {
        use_hueckel: true,
        charge: -1,
        ..Default::default()
    };
    assert!(value(determine_bonds_with_params(
        &mut atoms(ACRYLATE),
        &bonds
    )));

    // Methyl has 7 valence electrons in 7 orbitals, so its charge can be
    // from -7 to 6. Python's YAeHMOP exits the process past them.
    // The first atom keeps the charge in 8 bits: 262 is 6, and -200 is 56.
    for (charge, ok) in [
        (-8, false),
        (-7, true),
        (6, true),
        (7, false),
        (262, true),
        (-200, false),
    ] {
        let result =
            determine_connectivity_with_params(&mut atoms(METHYL), &hueckel_params(charge));
        assert_eq!(result.is_ok(), ok, "{charge}: {result:?}");
        assert!(ok || value(result));
    }
    // N's lines end the file, after which YAeHMOP rereads the last line
    // forever.
    let mut mol = atoms(METHYL);
    mol.get_atom_with_idx_mut(0).unwrap().set_atomic_num(7);
    assert!(value(determine_connectivity_with_params(
        &mut mol, &hueckel
    )));
}

#[test]
fn inputs_rdkit_trusts() {
    // Connect-the-dots casts each coordinate over 5.45 to int.
    for x in [f64::NAN, f64::INFINITY, 1e8] {
        let mut mol = ethanol();
        mol.get_conformer_mut()
            .unwrap()
            .set_atom_position(8, Point3D::new(x, 0.0, 0.0))
            .unwrap();
        assert!(value(determine_connectivity(&mut mol.clone())));
        assert!(value(determine_bonds(&mut mol.clone())));
        // The van der Waals method compares distances alone.
        let vdw = DetermineConnectivityParams {
            use_vdw: true,
            ..Default::default()
        };
        determine_connectivity_with_params(&mut mol, &vdw).unwrap();
    }

    // RDKit takes std::abs of the charge, undefined for i32::MIN.
    let radicals = DetermineBondsParams {
        allow_charged_fragments: false,
        charge: i32::MIN,
        ..Default::default()
    };
    assert!(value(determine_bonds_with_params(
        &mut atoms(METHYL),
        &radicals
    )));
}

#[test]
fn monomer_info_connect_the_dots_would_misread() {
    let pdb = |residue_number: i32| {
        AtomPDBResidueInfo::new_with_params(
            " X  ",
            &AtomPDBResidueInfoParams {
                residue_name: "ALA".to_owned(),
                residue_number,
                ..Default::default()
            },
        )
    };
    // Python's results: with no residues the hydrogen keeps its first
    // bond, and with PDB residues the bond within its residue.
    let mut mol = atoms(BRIDGE);
    determine_connectivity(&mut mol).unwrap();
    assert_eq!(bonds(&mol), vec![(0, 1)]);
    let mut mol = atoms(BRIDGE);
    for (i, residue) in [1, 2, 2].into_iter().enumerate() {
        mol.get_atom_with_idx_mut(i as u32)
            .unwrap()
            .set_pdb_residue_info(&pdb(residue));
    }
    determine_connectivity(&mut mol).unwrap();
    assert_eq!(bonds(&mol), vec![(1, 2)]);

    // RDKit would read other monomer info as PDB residue info, past its
    // end, whatever its type (U10); the call hides it, and puts it back.
    let mut mol = atoms(BRIDGE);
    for (i, name) in ["A", "B", "C"].into_iter().enumerate() {
        let info = AtomMonomerInfo::new_with_params(
            AtomMonomerType::Other,
            &AtomMonomerInfoParams {
                name: name.to_owned(),
            },
        )
        .unwrap();
        mol.get_atom_with_idx_mut(i as u32)
            .unwrap()
            .set_monomer_info(&info);
    }
    determine_connectivity(&mut mol).unwrap();
    assert_eq!(bonds(&mol), vec![(0, 1)]);
    for (i, name) in ["A", "B", "C"].into_iter().enumerate() {
        let info = mol
            .get_atom_with_idx(i as u32)
            .unwrap()
            .get_monomer_info()
            .unwrap();
        assert_eq!(info.get_name().unwrap(), name);
        assert_eq!(info.get_monomer_type(), AtomMonomerType::Other);
    }
}

#[test]
fn concurrent_bond_orders() {
    // Each call swaps the process's SIGINT handler through RDKit's static
    // variables; the wrapper runs one at a time.
    let handles: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                for _ in 0..5 {
                    let mut mol = atoms(ACRYLATE);
                    let params = DetermineBondsParams {
                        charge: -1,
                        ..Default::default()
                    };
                    determine_bonds_with_params(&mut mol, &params).unwrap();
                }
            })
        })
        .collect();
    for handle in handles {
        handle.join().unwrap();
    }
}
