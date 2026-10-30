// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rdinchi::*;
use rdkit_cxx::chem::{
    Atom, BondStereo, ChiralType, Conformer, Mol, RWMol, SmilesParserParams, mol_from_smarts,
    mol_from_smiles, mol_from_smiles_with_params, mol_to_smiles,
};
use rdkit_cxx::geometry::Point3D;

// Expected values come from Python's rdinchi.

const STANDARD_LOG: &str = "Generating standard InChI\nInput format: MOLfile\n\
                            Output format: Plain text\nFull Aux. info\nNo timeout\n\
                            Up to 1024 atoms per structure";

fn inchi(smiles: &str) -> MolToInchiResult {
    mol_to_inchi(&mol_from_smiles(smiles).unwrap()).unwrap()
}

fn inchi_with(smiles: &str, options: &str) -> Result<MolToInchiResult, Error> {
    mol_to_inchi_with_options(&mol_from_smiles(smiles).unwrap(), options)
}

fn is_value_error<T>(result: Result<T, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

#[test]
fn version() {
    assert_eq!(get_inchi_version(), "1.07.3");
}

#[test]
fn inchis_and_keys() {
    let cases = [
        (
            "CC(=O)O",
            "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)",
            "QTBSBXVTEAMEQO-UHFFFAOYSA-N",
        ),
        (
            "C[C@H](N)C(=O)O",
            "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
            "QNAYBMKLOCPYGJ-REOHCLBHSA-N",
        ),
        (
            "C/C=C/C",
            "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+",
            "IAQRGUVFOMOMEM-ONEGZZNKSA-N",
        ),
        (
            r"C[C@@H](Cl)/C=C\Br",
            "InChI=1S/C4H6BrCl/c1-4(6)2-3-5/h2-4H,1H3/b3-2-/t4-/m1/s1",
            "UHHVHSUMQAGKGR-BLDZNJJUSA-N",
        ),
        (
            "[13CH4]",
            "InChI=1S/CH4/h1H4/i1+1",
            "VNWKTOKETHGBQD-OUBTZVSYSA-N",
        ),
        (
            "[2H]C([2H])([2H])O",
            "InChI=1S/CH4O/c1-2/h2H,1H3/i1D3",
            "OKKJLVBELUTLKV-FIBGUPNXSA-N",
        ),
        (
            "c1ccccc1",
            "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
            "UHOVQNZJYSORNB-UHFFFAOYSA-N",
        ),
        ("[CH3]", "InChI=1S/CH3/h1H3", "WCYWZMWISLQXQU-UHFFFAOYSA-N"),
    ];
    for (smiles, expected, key) in cases {
        let result = inchi(smiles);
        assert_eq!(result.inchi, expected, "{smiles}");
        assert_eq!(
            mol_to_inchi_key(&mol_from_smiles(smiles).unwrap()).unwrap(),
            key,
            "{smiles}"
        );
        assert_eq!(inchi_to_inchi_key(expected), key, "{smiles}");
    }

    let acetic = inchi("CC(=O)O");
    assert_eq!(acetic.inchi, "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)");
    assert_eq!(acetic.retcode, 0);
    assert_eq!(acetic.message, "");
    assert_eq!(acetic.log, STANDARD_LOG);
    assert_eq!(
        acetic.aux_info,
        "AuxInfo=1/1/N:1,2,3,4/E:(3,4)/rA:4CCOO/rB:s1;d2;s2;/rC:;;;;"
    );
}

#[test]
fn charges_and_metals() {
    let methylammonium = inchi("C[NH3+]");
    assert_eq!(methylammonium.inchi, "InChI=1S/CH5N/c1-2/h2H2,1H3/p+1");
    assert_eq!(methylammonium.retcode, 1);
    assert_eq!(methylammonium.message, "Proton(s) added/removed");
    assert_eq!(
        methylammonium.log,
        format!("{STANDARD_LOG}\n\n\nWarning (Proton(s) added/removed)")
    );

    assert_eq!(inchi("[Na+].[Cl-]").inchi, "InChI=1S/ClH.Na/h1H;/q;+1/p-1");
    assert_eq!(
        inchi("[Cu+2].[O-]S(=O)(=O)[O-]").inchi,
        "InChI=1S/Cu.H2O4S/c;1-5(2,3)4/h;(H2,1,2,3,4)/q+2;/p-2"
    );
    let iron = inchi("CC(=O)O[Fe]");
    assert_eq!(
        iron.inchi,
        "InChI=1S/C2H4O2.Fe/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
    );
    assert_eq!(
        iron.message,
        "Accepted unusual valence(s): Fe(1); Metal was disconnected; Proton(s) added/removed"
    );
    let reconnected = inchi_with("CC(=O)O[Fe]", "/RecMet").unwrap();
    assert_eq!(
        reconnected.inchi,
        "InChI=1/C2H4O2.Fe/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1/rC2H3FeO2/c1-2(4)5-3/h1H3"
    );
}

#[test]
fn options() {
    let alanine = "C[C@H](N)C(=O)O";
    let cases = [
        (
            "/SNon",
            "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)",
            "QNAYBMKLOCPYGJ-UHFFFAOYSA-N",
        ),
        (
            "-snon",
            "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)",
            "QNAYBMKLOCPYGJ-UHFFFAOYSA-N",
        ),
        (
            "-SRel",
            "InChI=1/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/s2",
            "QNAYBMKLOCPYGJ-RUOZJJLTNA-N",
        ),
        (
            "/FixedH",
            "InChI=1/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1/f/h5H",
            "QNAYBMKLOCPYGJ-SNQCPAJUNA-N",
        ),
        (
            "\t/FixedH  /SUU ",
            "InChI=1/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1/f/h5H",
            "QNAYBMKLOCPYGJ-SNQCPAJUNA-N",
        ),
        (
            "/KET /15T",
            "InChI=1/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
            "QNAYBMKLOCPYGJ-REOHCLBHNA-N",
        ),
        (
            "   ",
            "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
            "QNAYBMKLOCPYGJ-REOHCLBHSA-N",
        ),
    ];
    let mol = mol_from_smiles(alanine).unwrap();
    for (options, expected, key) in cases {
        assert_eq!(
            mol_to_inchi_with_options(&mol, options).unwrap().inchi,
            expected,
            "{options:?}"
        );
        assert_eq!(
            mol_to_inchi_key_with_options(&mol, options).unwrap(),
            key,
            "{options:?}"
        );
    }

    let timed = mol_to_inchi_with_options(&mol, "/W5 /WM100 -AuxNone").unwrap();
    assert_eq!(timed.aux_info, "");
    assert!(timed.log.contains("Timeout per structure: 100 msec"));
    assert_eq!(
        inchi_with("OC1=NC=CC=C1", "/FixedH").unwrap().inchi,
        "InChI=1/C5H5NO/c7-5-3-1-2-4-6-5/h1-4H,(H,6,7)/f/h7H"
    );
    let no_hs = inchi_with("CCO", "/DoNotAddH").unwrap();
    assert_eq!(no_hs.inchi, "InChI=1S/C2O/c1-2-3");
    assert_eq!(
        no_hs.message,
        "Accepted unusual valence(s): C(1); C(2); O(1)"
    );
    for options in [
        "/SRac",
        "/SUCF",
        "/ChiralFlagON",
        "/ChiralFlagOFF",
        "/NEWPSOFF",
        "/SLUUD",
        "/LargeMolecules",
        "/Polymers",
        "/Polymers105",
        "/NPZz",
        "/SAtZZ",
        "/WarnOnEmptyStructure",
        "/NoWarnings",
        "/W0",
        "/WM0",
        "/w2.5",
        "/W1000000000",
        "/W1000000000.999",
        "/WM1000000000000",
    ] {
        let result = inchi_with("CCO", options).unwrap();
        assert!(
            result.inchi.ends_with("/C2H6O/c1-2-3/h3H,2H2,1H3"),
            "{options}"
        );
    }
}

#[test]
fn unsupported_options_are_refused() {
    // InChI ignores some of these and corrupts memory on others: "/OutputSDF
    // /Key" frees a pointer it never allocated.
    for options in [
        "/OutputSDF /Key",
        "/OutputSDF",
        "/Key",
        "/SNon /Key",
        "SNon",
        "/SNon\n/FixedH",
        "/",
        "-",
        "/SNonX",
        "/Tabbed",
        "/W",
        "/WM",
        "/Wnan",
        "/Winf",
        "/W1e3",
        "/W-5",
        "/W+5",
        "/W.5",
        "/W5.",
        "/W0x10",
        "/WM-5",
        "/WM1.5",
        "/W1000000001",
        "/W9223372036854775.808",
        "/WM1000000000001",
        "/WM99999999999999999999999999",
        "/SNon\0/Key",
        "\"/SNon\"",
    ] {
        assert!(is_value_error(inchi_with("CCO", options)), "{options:?}");
        assert!(
            is_value_error(mol_to_inchi_key_with_options(
                &mol_from_smiles("CCO").unwrap(),
                options
            )),
            "{options:?}"
        );
    }
}

#[test]
fn too_many_atoms_are_refused() {
    // InChI's atom indices are signed shorts.
    let mut mol = RWMol::from(Mol::default());
    for _ in 0..32767 {
        mol.add_atom(&Atom::new("He").unwrap());
    }
    assert!(is_value_error(mol_to_inchi(&mol)));
    assert!(is_value_error(mol_to_inchi_key(&mol)));
    mol.remove_atom(0).unwrap();
    assert_eq!(mol.get_num_atoms(), 32766);
    // As many atoms as InChI takes at all: it refuses them itself without
    // LargeMolecules.
    let result = mol_to_inchi(&mol).unwrap();
    assert_eq!(result.inchi, "");
    assert!(result.retcode >= 2);
}

#[test]
fn too_many_stereo_elements_are_refused() {
    // num_stereo0D is a signed short too. Each center has three neighbor
    // atoms, which InChI needs.
    let mut mol: Mol = mol_from_smiles(&"C(F)".repeat(16000)).unwrap();
    for idx in (2..31998).step_by(2) {
        mol.get_atom_with_idx_mut(idx)
            .unwrap()
            .set_chiral_tag(ChiralType::ChiTetrahedralCw);
    }
    for idx in 0..16770 {
        mol.get_bond_with_idx_mut(idx)
            .unwrap()
            .set_stereo(BondStereo::Stereoany);
    }
    assert!(is_value_error(mol_to_inchi(&mol)));
    assert!(is_value_error(mol_to_inchi_key(&mol)));
}

#[test]
fn stereo_centers_with_too_few_neighbor_atoms_are_refused() {
    // MolToInchi would pass InChI two neighbor indices it never set.
    let smarts = mol_from_smarts("[C@H](F)Cl").unwrap();
    let mut tagged = mol_from_smiles("C(F)Cl").unwrap();
    tagged
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiTetrahedralCw);
    for mol in [&smarts, &tagged] {
        assert!(matches!(
            mol_to_inchi(mol),
            Err(Error::Value { message }) if message.contains("atom 0")
        ));
        assert!(is_value_error(mol_to_inchi_key(mol)));
    }
    // With a hydrogen an atom, the center has three.
    let params = SmilesParserParams {
        remove_hs: false,
        ..Default::default()
    };
    let mut with_h = mol_from_smiles_with_params("C([H])(F)Cl", &params).unwrap();
    with_h
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiTetrahedralCw);
    assert!(mol_to_inchi(&with_h).is_ok());
}

#[test]
fn hydrogen_counts_inchi_cannot_hold_are_refused_on_reading() {
    // 200 deuteriums wrap to -56 in InChI's signed char, which RDKit reads
    // as about 4.3 billion atoms to add.
    for count in ["128", "200", "255", "1000", "4294967295"] {
        let text = format!("InChI=1S/C3H6/c1-2-3-1/h1-3H2/i1D{count}");
        assert!(is_value_error(inchi_to_mol(&text)), "{text}");
    }
    assert!(
        inchi_to_mol("InChI=1S/C3H6/c1-2-3-1/h1-3H2/i1D2")
            .unwrap()
            .mol
            .is_some()
    );
    // The formula's total count is not an atom's.
    let alkane = mol_from_smiles(&"C".repeat(70)).unwrap();
    let inchi = mol_to_inchi(&alkane).unwrap().inchi;
    assert!(inchi.contains("C70H142"));
    assert!(inchi_to_mol(&inchi).unwrap().mol.is_some());
}

#[test]
fn isotopes_inchi_cannot_represent_are_refused() {
    assert_eq!(inchi("[112CH4]").inchi, "InChI=1S/CH4/h1H4/i1+100");
    assert_eq!(inchi("[138U]").inchi, "InChI=1S/U/i1-100");
    assert_eq!(inchi("[2H]").inchi, "InChI=1S/H/i1+1");
    assert_eq!(inchi("[3H][H]").inchi, "InChI=1S/H2/h1H/i1+2");
    // Python gives /i1+105, /i1+104, and /i1-66 for these.
    for smiles in ["[113CH4]", "[70000CH4]", "[137U]"] {
        let mol = mol_from_smiles(smiles).unwrap();
        assert!(is_value_error(mol_to_inchi(&mol)), "{smiles}");
        assert!(is_value_error(mol_to_inchi_key(&mol)), "{smiles}");
    }
    let mut mol = mol_from_smiles("C").unwrap();
    mol.get_atom_with_idx_mut(0).unwrap().set_isotope(u32::MAX);
    assert!(is_value_error(mol_to_inchi(&mol)));
}

fn iron_with_hs(num_hs: u32) -> Mol {
    let mut mol = mol_from_smiles("[Fe]").unwrap();
    let mut atom = mol.get_atom_with_idx_mut(0).unwrap();
    atom.set_num_explicit_hs(num_hs);
    atom.set_no_implicit(true);
    mol
}

#[test]
fn hydrogen_counts_inchi_cannot_hold_are_refused() {
    let result = mol_to_inchi(&iron_with_hs(127)).unwrap();
    assert_eq!(result.inchi, "");
    assert_eq!(result.retcode, 3);
    assert_eq!(
        result.message,
        "Accepted unusual valence(s): Fe(127); Cannot disconnect metal error"
    );
    // Python raises for these. (RDKit keeps explicit counts in 8 bits, so
    // 256 is 0, in Python too.)
    for num_hs in [128, 200, 255] {
        assert!(
            is_value_error(mol_to_inchi(&iron_with_hs(num_hs))),
            "{num_hs}"
        );
        assert!(
            is_value_error(mol_to_inchi_key(&iron_with_hs(num_hs))),
            "{num_hs}"
        );
    }
    // InChI counts carbon's hydrogens itself, unless it is a radical.
    let mut carbon = mol_from_smiles("C").unwrap();
    {
        let mut atom = carbon.get_atom_with_idx_mut(0).unwrap();
        atom.set_num_explicit_hs(200);
        atom.set_no_implicit(true);
    }
    assert_eq!(mol_to_inchi(&carbon).unwrap().inchi, "InChI=1S/CH4/h1H4");
    carbon
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_num_radical_electrons(1);
    assert!(is_value_error(mol_to_inchi(&carbon)));
}

#[test]
fn coordinates_inchi_cannot_use_are_refused() {
    let mut mol = mol_from_smiles("C[C@H](O)/C=C/C").unwrap();
    let mut conf = Conformer::new(mol.get_num_atoms());
    for idx in 0..mol.get_num_atoms() {
        conf.set_atom_position(
            idx,
            Point3D::new(f64::from(idx), 0.5 * f64::from(idx % 2), 0.0),
        )
        .unwrap();
    }
    mol.add_conformer(&conf).unwrap();
    assert!(mol_to_inchi(&mol).is_ok());
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1.1e18, -1e300] {
        let mut copy = mol.clone();
        copy.get_conformer_mut()
            .unwrap()
            .set_atom_position(2, Point3D::new(0.0, bad, 0.0))
            .unwrap();
        assert!(is_value_error(mol_to_inchi(&copy)), "{bad}");
        assert!(is_value_error(mol_to_inchi_key(&copy)), "{bad}");
    }
}

#[test]
fn kekulization_errors() {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    let mol = mol_from_smiles_with_params("c1cccc1", &params).unwrap();
    // Python raises KekulizeException.
    assert!(matches!(mol_to_inchi(&mol), Err(Error::Kekulize { .. })));
    assert!(matches!(
        mol_to_inchi_key(&mol),
        Err(Error::Kekulize { .. })
    ));
}

#[test]
fn inchi_to_mol_reads_inchis() {
    let cases = [
        ("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)", "CC(=O)O"),
        (
            "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
            "C[C@H](N)C(=O)O",
        ),
        ("InChI=1S/CH4/h1H4/i1+1", "[13CH4]"),
        ("InChI=1S/C2H6/c1-2/h1-2H3/i1D", "[2H]CC"),
        ("InChI=1/C2H6/c1-2/h1-2H3", "CC"),
        ("InChI=1S/C/q-5", "[C-5]"),
    ];
    for (inchi, smiles) in cases {
        let result = inchi_to_mol(inchi).unwrap();
        assert_eq!(result.retcode, 0, "{inchi}");
        assert_eq!(result.message, "", "{inchi}");
        assert_eq!(
            result.log,
            "Input format: InChI (plain identifier)\nOutput format: Plain text\n\
             Aux. info suppressed\nNo timeout\nUp to 1024 atoms per structure",
            "{inchi}"
        );
        assert_eq!(mol_to_smiles(&result.mol.unwrap()).unwrap(), smiles);
    }

    let deuterated = "InChI=1S/C2H6/c1-2/h1-2H3/i1D";
    let kept = inchi_to_mol_with_params(
        deuterated,
        &InchiToMolParams {
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(mol_to_smiles(&kept.mol.unwrap()).unwrap(), "[2H]CC");
    let raw = inchi_to_mol_with_params(
        deuterated,
        &InchiToMolParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap()
    .mol
    .unwrap();
    assert_eq!(raw.get_num_atoms(), 3);
    assert_eq!(
        InchiToMolParams::default(),
        InchiToMolParams {
            sanitize: true,
            remove_hs: true
        }
    );

    // A round trip.
    let alanine = inchi("C[C@H](N)C(=O)O").inchi;
    let mol = inchi_to_mol(&alanine).unwrap().mol.unwrap();
    assert_eq!(mol_to_smiles(&mol).unwrap(), "C[C@H](N)C(=O)O");
}

#[test]
fn unreadable_inchis_have_no_molecule() {
    // Python returns (None, retcode, message, log).
    for (inchi, retcode, log) in [
        ("garbage", -1, ""),
        ("", -1, ""),
        (
            "InChI=1S/C2H6/c1-3",
            2,
            "Input format: InChI (plain identifier)\nOutput format: Plain text\n\
             Aux. info suppressed\nNo timeout\nUp to 1024 atoms per structure\n\n\n\n\
             Structure: 1 Syntax error (-2) in MOBILE_H_CONNECTIONS (1)",
        ),
    ] {
        let result = inchi_to_mol(inchi).unwrap();
        assert!(result.mol.is_none(), "{inchi:?}");
        assert_eq!(result.retcode, retcode, "{inchi:?}");
        assert_eq!(result.message, "", "{inchi:?}");
        assert_eq!(result.log, log, "{inchi:?}");
    }
    // A warning still gives the molecule.
    let mismatched = inchi_to_mol("InChI=1S/C2H6/c1-2/h1-2H3/q+9").unwrap();
    assert!(mismatched.mol.is_some());
    assert_eq!(mismatched.retcode, 1);
    assert_eq!(
        mismatched.message,
        " Problems/mismatches: Mobile-H( Charge(s): Do not match)"
    );
    // Sanitization errors are errors, as Python raises them.
    assert!(matches!(
        inchi_to_mol("InChI=1S/CH5/h1H5"),
        Err(Error::AtomValence { atom_idx: 0, .. })
    ));
    let unsanitized = inchi_to_mol_with_params(
        "InChI=1S/CH5/h1H5",
        &InchiToMolParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(unsanitized.mol.is_some());
}

#[test]
fn inchi_keys_of_strings() {
    assert_eq!(
        inchi_to_inchi_key("InChI=1/C2H6/c1-2/h1-2H3"),
        "OTMSDBZUPAUEDD-UHFFFAOYNA-N"
    );
    // InChI rejects these; RDKit logs why and returns an empty key.
    assert_eq!(inchi_to_inchi_key("garbage"), "");
    assert_eq!(inchi_to_inchi_key(""), "");
    // InChI checks only the layout of the string it hashes.
    assert_eq!(
        inchi_to_inchi_key("InChI=1S/C2H6/c1-3"),
        "HRLUEIONRMXDLA-UHFFFAOYSA-N"
    );
}

#[test]
fn from_threads() {
    let handles: Vec<_> = (0..4)
        .map(|i| {
            std::thread::spawn(move || {
                for _ in 0..20 {
                    let smiles = if i % 2 == 0 {
                        "C[C@H](N)C(=O)O"
                    } else {
                        "OC(=O)c1ccccc1"
                    };
                    let result = inchi(smiles);
                    let mol = inchi_to_mol(&result.inchi).unwrap().mol.unwrap();
                    assert_eq!(
                        mol_to_inchi_key(&mol).unwrap(),
                        inchi_to_inchi_key(&result.inchi)
                    );
                }
            })
        })
        .collect();
    for handle in handles {
        handle.join().unwrap();
    }
}
