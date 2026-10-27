// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's 2D descriptors. Expected values are from Python's
//! rdMolDescriptors, with `force=True` where Python would read a cache.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_descriptors::*;
use rdkit_cxx::chem::{Atom, BondType, Mol, RWMol, mol_from_smiles};

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

/// Paracetamol.
fn apap() -> Mol {
    smiles("CC(=O)Nc1ccc(O)cc1")
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn assert_all_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert_close(*a, *e);
    }
}

#[test]
fn crippen_labute_tpsa() {
    let mol = apap();
    let (logp, mr) = calc_crippen_descriptors(&mol).unwrap();
    assert_close(logp, 1.3505999999999998);
    assert_close(mr, 42.41050000000002);
    let (logp, mr) = calc_crippen_descriptors_with_include_hs(&mol, false).unwrap();
    assert_close(logp, 0.5430999999999998);
    assert_close(mr, 32.653800000000004);
    assert_close(calc_labute_asa(&mol).unwrap(), 64.66690538413799);
    assert_close(
        calc_labute_asa_with_include_hs(&mol, false).unwrap(),
        63.94513703834635,
    );
    assert_close(calc_tpsa(&mol).unwrap(), 49.33);
    let sulfonamide = smiles("CS(=O)(=O)NP");
    assert_close(calc_tpsa(&sulfonamide).unwrap(), 46.17);
    assert_close(
        calc_tpsa_with_include_s_and_p(&sulfonamide, true).unwrap(),
        54.550000000000004,
    );
}

#[test]
fn weights_and_formulas() {
    let mol = apap();
    assert_close(calc_exact_mol_wt(&mol).unwrap(), 151.063328528);
    assert_close(
        calc_exact_mol_wt_with_only_heavy(&mol, true).unwrap(),
        141.99290324,
    );
    assert_eq!(calc_mol_formula(&mol).unwrap(), "C8H9NO2");
    let labeled = smiles("[2H]C([2H])([2H])[13CH2]O");
    assert_eq!(calc_mol_formula(&labeled).unwrap(), "C2H6O");
    let separate = CalcMolFormulaParams {
        separate_isotopes: true,
        ..Default::default()
    };
    assert_eq!(
        calc_mol_formula_with_params(&labeled, &separate).unwrap(),
        "C[13C]H3D3O"
    );
    let unabbreviated = CalcMolFormulaParams {
        separate_isotopes: true,
        abbreviate_h_isotopes: false,
    };
    assert_eq!(
        calc_mol_formula_with_params(&labeled, &unabbreviated).unwrap(),
        "C[13C]H3[2H]3O"
    );
}

#[test]
fn params_defaults() {
    assert_eq!(
        CalcMolFormulaParams::default(),
        CalcMolFormulaParams {
            separate_isotopes: false,
            abbreviate_h_isotopes: true,
        }
    );
}

#[test]
fn counts() {
    let mol = apap();
    assert_eq!(calc_num_lipinski_hbd(&mol).unwrap(), 2);
    assert_eq!(calc_num_lipinski_hba(&mol), 3);
    assert_eq!(calc_num_hbd(&mol).unwrap(), 2);
    assert_eq!(calc_num_hba(&mol).unwrap(), 2);
    assert_eq!(calc_num_heavy_atoms(&mol), 11);
    assert_eq!(calc_num_atoms(&mol).unwrap(), 20);
    assert_eq!(calc_num_heteroatoms(&mol).unwrap(), 3);
    assert_eq!(calc_num_amide_bonds(&mol).unwrap(), 1);

    let rotors = smiles("CC(=O)NC1CCC(c2ccccc2)CC1C(=O)OC");
    assert_close(calc_fraction_csp3(&rotors).unwrap(), 0.5);
    assert_eq!(calc_num_rotatable_bonds(&rotors).unwrap(), 3);
    assert_eq!(
        calc_num_rotatable_bonds_with_strict(&rotors, true).unwrap(),
        3
    );
    assert_eq!(
        calc_num_rotatable_bonds_with_strict(&rotors, false).unwrap(),
        5
    );
    for (option, expected) in [
        (NumRotatableBondsOptions::NonStrict, 5),
        (NumRotatableBondsOptions::Strict, 3),
        (NumRotatableBondsOptions::StrictLinkages, 4),
        (NumRotatableBondsOptions::Default, 3),
    ] {
        assert_eq!(
            calc_num_rotatable_bonds_with_strict(&rotors, option).unwrap(),
            expected
        );
    }
    assert_eq!(
        NumRotatableBondsOptions::from(true),
        NumRotatableBondsOptions::Strict
    );
}

#[test]
fn ring_counts() {
    let mol = smiles("c1ccccc1C1CCNCC1c1ccncc1C1CCCC1C1=CCOC1");
    assert_eq!(calc_num_rings(&mol).unwrap(), 5);
    assert_eq!(calc_num_aromatic_rings(&mol), 2);
    assert_eq!(calc_num_saturated_rings(&mol), 2);
    assert_eq!(calc_num_heterocycles(&mol), 3);
    assert_eq!(calc_num_aromatic_heterocycles(&mol), 1);
    assert_eq!(calc_num_aromatic_carbocycles(&mol), 1);
    assert_eq!(calc_num_saturated_heterocycles(&mol), 1);
    assert_eq!(calc_num_saturated_carbocycles(&mol), 1);
    assert_eq!(calc_num_aliphatic_rings(&mol), 3);
    assert_eq!(calc_num_aliphatic_heterocycles(&mol), 2);
    assert_eq!(calc_num_aliphatic_carbocycles(&mol), 1);

    let polycycle = smiles("C1CCC12CCCC2C1CC2CCC1C2");
    let (count, outputs) = calc_num_spiro_atoms_with_outputs(&polycycle).unwrap();
    assert_eq!((count, outputs.atoms), (1, vec![3]));
    assert_eq!(calc_num_spiro_atoms(&polycycle).unwrap(), 1);
    let (count, outputs) = calc_num_bridgehead_atoms_with_outputs(&polycycle).unwrap();
    assert_eq!((count, outputs.atoms), (2, vec![10, 13]));
    assert_eq!(calc_num_bridgehead_atoms(&polycycle).unwrap(), 2);

    let stereo = smiles("C[C@H](O)C(F)C(Cl)Br");
    assert_eq!(calc_num_atom_stereo_centers(&stereo).unwrap(), 3);
    assert_eq!(
        calc_num_unspecified_atom_stereo_centers(&stereo).unwrap(),
        2
    );
}

#[test]
fn rings_not_found() {
    // A molecule built atom by atom has no ring info.
    let mut mol = RWMol::default();
    mol.add_atom(&Atom::from_num(6));
    assert!(matches!(calc_num_rings(&mol), Err(Error::Runtime { .. })));
    assert_eq!(calc_num_aromatic_rings(&mol), 0);
}

#[test]
fn vsa() {
    let mol = apap();
    assert_all_close(
        &slogp_vsa(&mol).unwrap(),
        &[
            5.316788604006331,
            11.013707124192212,
            4.794537184071822,
            0.0,
            6.923737199690624,
            24.26546827384644,
            0.0,
            0.0,
            0.0,
            5.687386274683562,
            5.749511833283905,
            0.0,
        ],
    );
    assert_all_close(
        &slogp_vsa_with_bins(&mol, &[0.0, 0.2]).unwrap(),
        &[21.125032912270363, 31.18920547353706, 11.436898107967467],
    );
    assert_eq!(
        slogp_vsa_with_bins(&mol, &[]).unwrap(),
        slogp_vsa(&mol).unwrap()
    );
    assert_all_close(
        &smr_vsa(&mol).unwrap(),
        &[
            9.901064578912528,
            0.0,
            0.0,
            0.0,
            6.923737199690624,
            5.316788604006331,
            24.26546827384644,
            0.0,
            5.749511833283905,
            11.594566004035068,
        ],
    );
    assert_eq!(peoe_vsa(&mol).unwrap().len(), 14);
    assert_close(peoe_vsa(&mol).unwrap()[0], 10.423315998847038);
    assert_all_close(
        &peoe_vsa_with_bins(&mol, &[-0.1, 0.1]).unwrap(),
        &[15.217853182918859, 36.876591748220626, 11.65669156263541],
    );
    let mut custom = apap();
    for idx in (1..custom.get_num_atoms()).step_by(2) {
        custom
            .get_atom_with_idx_mut(idx)
            .unwrap()
            .set_double_prop("x", 2.5);
    }
    assert_all_close(
        &custom_prop_vsa(&custom, "x", &[1.0, 2.0]).unwrap(),
        &[0.0, 34.64492219020993, 29.10621430356496],
    );
}

#[test]
fn vsa_bins_must_be_sorted() {
    // std::upper_bound is undefined for unsorted bins and for NaN.
    let mol = apap();
    for bins in [[0.2, 0.0], [0.0, f64::NAN]] {
        assert!(matches!(
            slogp_vsa_with_bins(&mol, &bins),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            smr_vsa_with_bins(&mol, &bins),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            peoe_vsa_with_bins(&mol, &bins),
            Err(Error::Value { .. })
        ));
        assert!(matches!(
            custom_prop_vsa(&mol, "x", &bins),
            Err(Error::Value { .. })
        ));
    }
    // Equal boundaries are sorted.
    assert_eq!(slogp_vsa_with_bins(&mol, &[0.0, 0.0]).unwrap().len(), 3);
}

#[test]
fn mqns_values() {
    assert_eq!(
        mqns(&apap()).unwrap(),
        vec![
            8, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 11, 4, 1, 0, 3, 3, 0, 1, 5, 3, 2, 2, 0, 0, 3, 1, 1, 0,
            4, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
        ]
    );
}

#[test]
fn connectivity() {
    let mol = apap();
    let chi = [
        6.164862962722325,
        3.2490981480278296,
        2.228060242189872,
        1.1876733148984595,
        0.7480265991539318,
    ];
    let v = [calc_chi0v, calc_chi1v, calc_chi2v, calc_chi3v, calc_chi4v];
    let n = [calc_chi0n, calc_chi1n, calc_chi2n, calc_chi3n, calc_chi4n];
    for i in 0..5 {
        assert_close(v[i](&mol).unwrap(), chi[i]);
        assert_close(n[i](&mol).unwrap(), chi[i]);
    }
    assert_close(calc_chi_nv(&mol, 5).unwrap(), 0.45410082788731426);
    assert_close(calc_chi_nn(&mol, 5).unwrap(), 0.45410082788731426);
    assert_close(calc_chi_nv(&mol, 0).unwrap(), 11.0);
    assert_eq!(calc_chi_nv(&mol, 11).unwrap(), 0.0);
    // Longer than any path, where n + 1 would wrap.
    assert_eq!(calc_chi_nv(&mol, 12).unwrap(), 0.0);
    assert_eq!(calc_chi_nn(&mol, u32::MAX).unwrap(), 0.0);

    let (alpha, outputs) = calc_hall_kier_alpha_with_outputs(&mol).unwrap();
    assert_close(alpha, -1.5099999999999998);
    assert_all_close(
        &outputs.atom_contribs,
        &[
            0.0, -0.13, -0.2, -0.2, -0.13, -0.13, -0.13, -0.13, -0.2, -0.13, -0.13,
        ],
    );
    assert_close(calc_hall_kier_alpha(&mol).unwrap(), alpha);
    assert_close(calc_kappa1(&mol).unwrap(), 7.595374077976819);
    assert_close(calc_kappa2(&mol).unwrap(), 3.0531381005525);
    assert_close(calc_kappa3(&mol).unwrap(), 2.292303972882069);
    assert_close(calc_phi(&mol).unwrap(), 2.1081569077654403);
}

#[test]
fn stale_chi_cache() {
    // Python caches each atom's delta and, after an edit, reads the cached
    // vector past its end (CalcChi1v gives the old molecule's 1.0233...).
    // Here every call computes afresh (U8).
    let mut mol = RWMol::from(smiles("CCO"));
    assert_close(calc_chi1v(&mol).unwrap(), 1.0233345472033855);
    let idx = mol.add_atom(&Atom::from_num(6));
    mol.add_bond_with_order(2, idx, BondType::Single).unwrap();
    mol.update_property_cache().unwrap();
    assert_close(calc_chi1v(&mol).unwrap(), 1.4040302062452235);
    assert_close(calc_chi0v(&mol).unwrap(), 3.1153550716504106);
    assert!(!mol.has_prop("_connectivityHKDeltas"));
}

#[test]
fn atom_features() {
    let mol = apap();
    let features = get_atom_features(&mol, 3).unwrap();
    assert_eq!(features.len(), 49);
    let on: Vec<usize> = (0..49).filter(|&i| features[i] == 1.0).collect();
    assert_eq!(on, vec![2, 16, 22, 27, 34, 44]);
    assert_close(features[48], 1.0 / 11.0);
    let chiral = get_atom_features_with_addchiral(&smiles("C[C@H](N)C(F)Cl"), 1, true).unwrap();
    assert_eq!(chiral.len(), 52);
    assert_eq!(chiral[49..], [0.0, 1.0, 0.0]);
    assert!(matches!(
        get_atom_features(&mol, 11),
        Err(Error::Index { index: 11, len: 11 })
    ));
}

#[test]
fn oxidation_numbers() {
    let mut mol = smiles("CC(=O)O.[Fe+3]");
    calc_oxidation_numbers(&mut mol).unwrap();
    let numbers: Vec<i32> = (0..mol.get_num_atoms())
        .map(|i| {
            mol.get_atom_with_idx(i)
                .unwrap()
                .get_int_prop("OxidationNumber")
                .unwrap()
        })
        .collect();
    assert_eq!(numbers, vec![-3, 3, -2, -2, 3]);
}

#[test]
fn autocorr_2d() {
    let mol = apap();
    let values = calc_autocorr_2d(&mol).unwrap();
    assert_eq!(values.len(), 192);
    assert_all_close(&values[..5], &[2.565, 2.836, 2.59, 2.319, 2.106]);
    assert_all_close(&values[100..103], &[0.017, -0.653, 1.377]);
    let mut custom = apap();
    for idx in 0..custom.get_num_atoms() {
        custom
            .get_atom_with_idx_mut(idx)
            .unwrap()
            .set_double_prop("x", f64::from(idx) * 0.5);
    }
    let values = calc_autocorr_2d_with_custom_atom_property(&custom, "x").unwrap();
    assert_eq!(values.len(), 32);
    assert_all_close(
        &values[..10],
        &[
            4.505, 4.718, 4.398, 3.239, 2.788, 1.872, 1.609, 0.0, 20.75, 19.75,
        ],
    );
    assert!(!mol.has_prop("DistanceMatrix"));
}

#[test]
fn property_tables_cover_1_to_110() {
    // RDKit reads table[atomic number - 1] from 110-entry tables: before the
    // start for a dummy atom, past the end for oganesson.
    for text in ["*O", "[Og]O"] {
        assert!(matches!(
            calc_autocorr_2d(&smiles(text)),
            Err(Error::Value { .. })
        ));
    }
    assert_eq!(calc_autocorr_2d(&smiles("[Ds]O")).unwrap().len(), 192);
    // A custom property reads no table.
    assert_eq!(
        calc_autocorr_2d_with_custom_atom_property(&smiles("*O"), "x")
            .unwrap()
            .len(),
        32
    );
}

#[test]
fn bcut() {
    let mol = apap();
    assert_all_close(
        &bcut_2d(&mol).unwrap(),
        &[
            16.254688714379643,
            10.280190859964746,
            1.9280360590359849,
            -1.9945540738765895,
            2.125247579160181,
            -2.044011312610291,
            5.88095478459792,
            -0.11373578338265988,
        ],
    );
    let props: Vec<f64> = (0..11).map(f64::from).collect();
    let (high, low) = bcut_2d_atom_props(&mol, &props).unwrap();
    assert_close(high, 10.592249293030436);
    assert_close(low, -0.7662125043492891);
    assert!(matches!(
        bcut_2d_atom_props(&mol, &props[1..]),
        Err(Error::Runtime { .. })
    ));
    let mut custom = apap();
    for idx in 0..custom.get_num_atoms() {
        custom
            .get_atom_with_idx_mut(idx)
            .unwrap()
            .set_double_prop("x", f64::from(idx) * 0.5);
    }
    let (high, low) = bcut_2d_atom_propname(&custom, "x").unwrap();
    assert_close(high, 5.891312409196216);
    assert_close(low, -1.0960723664111496);
    assert!(matches!(
        bcut_2d_atom_propname(&mol, "x"),
        Err(Error::Key { .. })
    ));
}

#[test]
fn maccs_keys() {
    let fp = get_maccs_keys_fingerprint(&apap()).unwrap();
    assert_eq!(fp.get_num_bits(), 167);
    assert_eq!(
        fp.get_on_bits(),
        vec![
            92, 110, 113, 117, 127, 131, 133, 135, 139, 143, 151, 152, 154, 156, 157, 158, 159,
            160, 161, 162, 163, 164, 165
        ]
    );
    assert!(
        get_maccs_keys_fingerprint(&Mol::default())
            .unwrap()
            .get_on_bits()
            .is_empty()
    );
}

#[test]
fn maccs_keys_in_threads() {
    // RDKit makes the MACCS patterns on first use, without a lock; the
    // wrapper makes them once before any call.
    let threads: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                get_maccs_keys_fingerprint(&apap())
                    .unwrap()
                    .get_num_on_bits()
            })
        })
        .collect();
    for thread in threads {
        assert_eq!(thread.join().unwrap(), 23);
    }
}

#[test]
fn atom_pair_codes() {
    let mol = apap();
    let codes: Vec<u32> = mol
        .get_atoms()
        .map(|atom| get_atom_pair_atom_code(&atom).unwrap())
        .collect();
    assert_eq!(codes, vec![33, 43, 105, 66, 43, 42, 42, 43, 97, 42, 42]);
    let branch = GetAtomPairAtomCodeParams {
        branch_subtract: 1,
        ..Default::default()
    };
    let codes: Vec<u32> = mol
        .get_atoms()
        .map(|atom| get_atom_pair_atom_code_with_params(&atom, &branch).unwrap())
        .collect();
    assert_eq!(codes, vec![32, 42, 104, 65, 42, 41, 41, 42, 96, 41, 41]);
    let chiral = smiles("C[C@H](N)O");
    let chirality = GetAtomPairAtomCodeParams {
        include_chirality: true,
        ..Default::default()
    };
    let codes: Vec<u32> = chiral
        .get_atoms()
        .map(|atom| get_atom_pair_atom_code_with_params(&atom, &chirality).unwrap())
        .collect();
    assert_eq!(codes, vec![33, 547, 65, 97]);
    assert!(matches!(
        get_atom_pair_atom_code(&Atom::from_num(6)),
        Err(Error::Runtime { .. })
    ));
    assert_eq!(GetAtomPairAtomCodeParams::default().branch_subtract, 0);

    assert_eq!(get_atom_pair_code(65, 98, 3).unwrap(), 1607715);
    assert_eq!(get_atom_pair_code(98, 65, 3).unwrap(), 1607715);
    assert_eq!(
        get_atom_pair_code_with_include_chirality(65, 98, 3, true).unwrap(),
        6424611
    );
    assert_eq!(
        get_atom_pair_code(u32::MAX, u32::MAX, 30).unwrap(),
        4294967294
    );
    assert!(matches!(
        get_atom_pair_code(1, 2, 31),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn morgan_invariants() {
    let mol = apap();
    assert_eq!(
        get_connectivity_invariants(&mol).unwrap(),
        vec![
            2246728737, 2246699815, 864942730, 847961216, 3217380708, 3218693969, 3218693969,
            3217380708, 864662311, 3218693969, 3218693969
        ]
    );
    assert_eq!(
        get_connectivity_invariants_with_include_ring_membership(&mol, false).unwrap(),
        vec![
            2246728737, 2246699815, 864942730, 847961216, 2246699815, 2246703798, 2246703798,
            2246699815, 864662311, 2246703798, 2246703798
        ]
    );
    assert_eq!(
        get_feature_invariants(&mol).unwrap(),
        vec![0, 0, 2, 1, 4, 4, 4, 4, 3, 4, 4]
    );
}

#[test]
fn atom_pairs_parameters() {
    assert_eq!(AtomPairsParameters::version(), "1.1.0");
    assert_eq!(AtomPairsParameters::num_type_bits(), 4);
    assert_eq!(AtomPairsParameters::num_pi_bits(), 2);
    assert_eq!(AtomPairsParameters::num_branch_bits(), 3);
    assert_eq!(AtomPairsParameters::num_chiral_bits(), 2);
    assert_eq!(AtomPairsParameters::code_size(), 9);
    assert_eq!(AtomPairsParameters::num_path_bits(), 5);
    assert_eq!(AtomPairsParameters::num_atom_pair_fingerprint_bits(), 23);
    // RDKit 2024.03 lists technetium (43) where later releases have
    // iodine (53).
    assert_eq!(
        AtomPairsParameters::atom_types(),
        vec![5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 43, 0]
    );
}
