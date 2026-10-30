// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports RDKit's rdkit/Chem/UnitTestFragmentDescriptors.py (RDKit
// 2024.09.6), and pins values from Python's rdkit.Chem.Fragments.

use rdkit_cxx::Result;
use rdkit_cxx::chem::{Mol, fragments, mol_from_smiles};

fn run(f: fn(&Mol) -> Result<u32>, data: &[(&str, u32)]) {
    for &(smi, tgt_val) in data {
        let m = mol_from_smiles(smi).unwrap();
        assert_eq!(f(&m).unwrap(), tgt_val, "bad value for smiles {smi}");
    }
}

#[test]
fn test1_fr_c_o() {
    run(
        fragments::fr_c_o,
        &[
            ("C=O", 1),
            ("O=CC=O", 2),
            ("O=CC(=O)O", 2),
            ("O=C(O)C(=O)O", 2),
        ],
    );
}

#[test]
fn test2_fr_c_o_no_coo() {
    run(
        fragments::fr_c_o_no_coo,
        &[
            ("C=O", 1),
            ("O=CC=O", 2),
            ("O=CC(=O)O", 1),
            ("O=C(O)C(=O)O", 0),
        ],
    );
}

#[test]
fn test3_fr_nitro_arom() {
    run(
        fragments::fr_nitro_arom,
        &[
            ("[O-][N+](C1=C([N+]([O-])=O)C=C([N+]([O-])=O)C=C1)=O", 3),
            ("[O-][N+](C1=CC=CC=C1)=O", 1),
            ("[N+]([O-])(=O)C1OC(=CC=1)/C=N/NC(N)=O", 0),
            ("O=C(O)C(=O)O", 0),
        ],
    );
}

#[test]
fn test4_fr_nitro_arom_nonortho() {
    run(
        fragments::fr_nitro_arom_nonortho,
        &[
            ("[O-][N+](C1=C([N+]([O-])=O)C=C([N+]([O-])=O)C=C1)=O", 1),
            ("[O-][N+](C1=CC=CC=C1)=O", 1),
            ("[N+]([O-])(=O)C1OC(=CC=1)/C=N/NC(N)=O", 0),
            ("O=C(O)C(=O)O", 0),
        ],
    );
}

#[test]
fn test5_fr_azide() {
    run(
        fragments::fr_azide,
        &[("C1(N=[N+]=[N-])=C(C(=C(C(=C1F)F)F)F)F", 1)],
    );
}

#[test]
fn test6_fr_dihydropyridine() {
    run(
        fragments::fr_dihydropyridine,
        &[
            ("C1C=CNC=C1", 1),
            ("C1=CCNC=C1", 1),
            ("C1=CCN=CC1", 1),
            ("C1=CC=NCC1", 1),
            (
                "C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O",
                1,
            ),
            ("O=C(O)C(=O)O", 0),
        ],
    );
}

#[test]
fn test7_fr_phenol_no_ortho_hbond() {
    run(
        fragments::fr_phenol_no_ortho_hbond,
        &[
            ("OC1=C(C(O)=O)C=CC=C1", 0),
            ("OC1=CC=CC=C1", 1),
            ("OC1=C(O)C=CC=C1", 2),
            ("O=C(O)C1=CC(O)=CC=C1", 1),
            ("OC1=C(O)C(C(O)=O)=CC=C1", 1),
            ("OC1=C(C(N)=O)C=CC=C1", 0),
            ("OC1=C(CO)C=CC=C1", 0),
            ("OC1=C(C(OC)=O)C=CC=C1", 1),
            ("C1=CC=NCC1", 0),
            (
                "C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O",
                0,
            ),
            ("O=C(O)C(=O)O", 0),
        ],
    );
}

#[test]
fn test8_fr_al_oh_no_tert() {
    run(
        fragments::fr_al_oh_no_tert,
        &[
            ("OC1=C(C(O)=O)C=CC=C1", 0),
            ("CC(C)(O)C", 0),
            ("CC(O)C", 1),
            ("C=C(O)C", 1),
            ("OC1=C(O)C=CC=C1", 0),
            ("OC1=C(CO)C=CC=C1", 1),
            (
                "C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O",
                0,
            ),
        ],
    );
}

#[test]
fn test9_fr_benzodiazepine() {
    run(
        fragments::fr_benzodiazepine,
        &[
            ("OC1=C(C(O)=O)C=CC=C1", 0),
            ("C12=CC=CC=C1NCCN=C2", 1),
            ("C1=C(C=C2C(=C1)NC(CN=C2C3=C(C=CC=C3)Cl)=O)[N+]([O-])=O", 1),
            ("C12=C(C(=NCC(N1C)=O)C3=CC=CC=C3)C=C(C=C2)Cl", 1),
            ("CCC1=CC2C(=NCC(=O)N(C)C=2S1)C3=CC=CC=C3Cl", 0),
            ("CN(C)CCN1C(=O)C2=CC=CC=C2NC3C=C(Cl)C=CC1=3", 0),
            (
                "C1C=CC(=CC=1[N+]([O-])=O)C2C(=C(NC(=C2C(OC)=O)C)C)C(OCCNCC3=CC=CC=C3)=O",
                0,
            ),
        ],
    );
}

#[test]
fn test10_fr_para_hydroxylation() {
    run(
        fragments::fr_para_hydroxylation,
        &[
            ("OC1=C(C(O)=O)C=CC=C1", 1),
            ("c1ccccc1Oc2ccccc2", 2),
            ("COC1=CC=CC=C1", 1),
            ("CN(C)C2=CC=CC=C2", 1),
            ("O=C(C4=CC=CC=C4)C3=CC=CC=C3", 0),
            ("COC1=CC=C(C)C(OC)=C1C", 0),
            ("CN(C)C1=CC=CC=C1OC", 2),
            ("O=C(C2=CC=CC=C2)CC1=CC=CC=C1", 0),
            ("O=C(NC2=CC=CC=C2)CN(C)C1=CC=CC=C1", 2),
            ("O=C(N(C)C2=CC=CC=C2)CC1=CC=CC=C1", 1),
            ("C12=CC=CC=C1OCCC2", 1),
            ("C12=CC=CC=C1OC=C2", 1),
            ("C12=CC=CC=C1NC=N2", 2),
            ("CC1=CC=C(C=CO2)C2=C1", 0),
        ],
    );
}

#[test]
fn test11_fr_allylic_oxid() {
    run(
        fragments::fr_allylic_oxid,
        &[
            (
                "C1C(C=C2[C@](C1)([C@]3([C@@](CC2)([C@]4([C@](CC3)([C@@H](CC4)O)C)[H])[H])[H])C)=O",
                1,
            ),
            ("CCC=CCCC", 2),
            ("C=CC", 1),
            ("C=CCO", 0),
        ],
    );
}

#[test]
fn test12_fr_aryl_methyl() {
    run(
        fragments::fr_aryl_methyl,
        &[
            ("c1ccccc1C", 1),
            ("c1ccccc1CC", 1),
            ("c1(C)cccc(C)c1CC", 2),
            ("c1(N)cccc(O)c1CC", 0),
            ("c1cccc(C)c1CC", 2),
            ("c1ccccc1CN", 0),
            ("c1ccccc1CCN", 0),
            ("c1(C)ccccc1CCN", 1),
        ],
    );
}

#[test]
fn test13_fr_ndealkylation1() {
    run(
        fragments::fr_ndealkylation1,
        &[
            ("CNC1(CCCCC1=O)C2(=CC=CC=C2Cl)", 1),
            ("C1=C(C(=C(C=C1)C)NC(CN(CC)CC)=O)C", 1),
            ("c1(C)ccccc1CCN", 0),
        ],
    );
}

#[test]
fn test14_fr_ndealkylation2() {
    run(
        fragments::fr_ndealkylation2,
        &[
            (
                "C12(=C(N=CC=C1[C@@H](O)C3(N4CC(C(C3)CC4)([H])C=C)[H])C=CC(=C2)OC)",
                0,
            ),
            ("N12CCC(CC2)CC1", 0),
            ("N1(CCC)CCCCC1", 1),
            ("CCC1(CCCCN(C)C1)C2=CC=CC(O)=C2", 1),
            (
                "C1(=C2C3=C(C=C1)C[C@@H]4[C@]5([C@@]3([C@H]([C@H](CC5)O)O2)CCN4CC6CCC6)O)O",
                1,
            ),
            ("CC1N=C2CCCCN2C(=O)C=1CCN3CCC(CC3)C4=NOC5C=C(F)C=CC4=5", 1),
            ("N1CCOCC1", 0),
            ("n1ccccc1", 0),
        ],
    );
}

#[test]
fn test15_fr_alkyl_carbamate() {
    run(
        fragments::fr_alkyl_carbamate,
        &[
            ("C(COC(NC(C)C)=O)(COC(N)=O)(CCC)C", 1),
            (
                "CN(CC3=CSC(C(C)C)=N3)C(N[C@@H]([C@H](C)C)C(N[C@@H](CC4=CC=CC=C4)C[C@H](O)[C@H](CC2=CC=CC=C2)NC(OCC1=CN=CS1)=O)=O)=O",
                1,
            ),
            ("c1(C)ccccc1CCN", 0),
        ],
    );
}

#[test]
fn test16_fr_lactone() {
    run(
        fragments::fr_lactone,
        &[
            (r"C1(CCC(O1)=O)(C)CC/C=C\CC", 1),
            (
                "CN(CC3=CSC(C(C)C)=N3)C(N[C@@H]([C@H](C)C)C(N[C@@H](CC4=CC=CC=C4)C[C@H](O)[C@H](CC2=CC=CC=C2)NC(OCC1=CN=CS1)=O)=O)=O",
                0,
            ),
            ("c1(C)ccccc1CCN", 0),
        ],
    );
}

#[test]
fn test17_fr_ketone_topliss() {
    run(
        fragments::fr_ketone_topliss,
        &[
            ("O=C(C=CC1=CC=CC=C1)C=CC2=CC=CC=C2", 0),
            ("c1ccccc1-C(=O)-c2ccccc2", 0),
            ("c1ccccc1-C(=O)-CCO", 1),
            ("CC(=O)NCCC", 0),
        ],
    );
}

#[test]
fn test18_fr_arn() {
    run(
        fragments::fr_arn,
        &[
            (
                "C1=CC(=CC=C1C(N[C@@H](CCC(O)=O)C(O)=O)=O)N(CC2=NC3C(=NC(=NC=3N=C2)N)N)C",
                2,
            ),
            ("S(NC1=NC=CC=N1)(=O)(=O)C2=CC=C(C=C2)N", 1),
            ("S(NC1=NC=CC=N1)(=O)(=O)C2=CC=C(C=C2)", 0),
            ("c1ccccc1-C(=O)-CCO", 0),
            ("NNC(=O)C1C2CCCCC=2N=C3C=CC=CC=13", 1),
            ("NC(=N)C1C2CCCCC=2N=C3C=CC=CC=13", 1),
            ("NNC1C2CCCCC=2N=C3C=CC=CC=13", 1),
            ("NC1C2CCCCC=2N=C3C=CC=CC=13", 1),
        ],
    );
}

#[test]
fn test19_fr_hoccn() {
    run(
        fragments::fr_hoccn,
        &[
            ("CC(C)(C)NCC(O)C1=CC(O)=CC(O)=C1", 1),
            (r"OCCN1CCN(CC/C=C2\C3=CC=CC=C3SC4C=CC(Cl)=CC2=4)CC1", 1),
            ("c1ccccc1-C(=O)-CCO", 0),
            ("CC(=O)NCCC", 0),
        ],
    );
}

#[test]
fn test20_fr_methoxy() {
    run(
        fragments::fr_methoxy,
        &[
            ("c1ccccc1OC", 1),
            ("c1ccccc1Oc2ccccc2", 0),
            ("c1ccccc1OCC", 0),
        ],
    );
}

#[test]
fn test21_fr_oxime() {
    run(
        fragments::fr_oxime,
        &[
            (r"C/C(C(C)=O)=N\O", 1),
            (
                "C(=N/OC(C(=O)O)(C)C)(/C1=CS[N+](=N1)C)C(N[C@@H]2C(N([C@@H]2C)S(O)(=O)=O)=O)=O",
                1,
            ),
            ("c1ccccc1OCC", 0),
        ],
    );
}

#[test]
fn test22_fr_alkyl_halide() {
    run(
        fragments::fr_alkyl_halide,
        &[
            ("CCCF", 1),
            ("c1ccccc1C(C)(C)Br", 1),
            ("c1ccccc1CC(C)(C)Cl", 1),
        ],
    );
}

/// Every fragment count, in Python's order, unique and not (countUnique),
/// from Python's Fragments.fns.
const FNS_CASES: [(&str, &[u32], &[u32]); 4] = [
    (
        "O=C(O)c1ccccc1C(=O)OCC",
        &[
            2, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
        &[
            2, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "CC(=O)Nc1ccc(O)cc1",
        &[
            1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
        &[
            1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "c1ccc2c(c1)[nH]c1ccccc12",
        &[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
        &[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 6, 24, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
    (
        "NC(=O)CC(N)C(=O)O",
        &[
            2, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
        &[
            2, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ],
    ),
];

#[test]
fn fns_match_python() {
    let fns = fragments::fns();
    let names: Vec<&str> = fns.iter().map(|(name, _)| *name).collect();
    assert_eq!(
        names,
        [
            "fr_C_O",
            "fr_C_O_noCOO",
            "fr_Al_OH",
            "fr_Ar_OH",
            "fr_methoxy",
            "fr_oxime",
            "fr_ester",
            "fr_Al_COO",
            "fr_Ar_COO",
            "fr_COO",
            "fr_COO2",
            "fr_ketone",
            "fr_ether",
            "fr_phenol",
            "fr_aldehyde",
            "fr_quatN",
            "fr_NH2",
            "fr_NH1",
            "fr_NH0",
            "fr_Ar_N",
            "fr_Ar_NH",
            "fr_aniline",
            "fr_Imine",
            "fr_nitrile",
            "fr_hdrzine",
            "fr_hdrzone",
            "fr_nitroso",
            "fr_N_O",
            "fr_nitro",
            "fr_azo",
            "fr_diazo",
            "fr_azide",
            "fr_amide",
            "fr_priamide",
            "fr_amidine",
            "fr_guanido",
            "fr_Nhpyrrole",
            "fr_imide",
            "fr_isocyan",
            "fr_isothiocyan",
            "fr_thiocyan",
            "fr_halogen",
            "fr_alkyl_halide",
            "fr_sulfide",
            "fr_SH",
            "fr_C_S",
            "fr_sulfone",
            "fr_sulfonamd",
            "fr_prisulfonamd",
            "fr_barbitur",
            "fr_urea",
            "fr_term_acetylene",
            "fr_imidazole",
            "fr_furan",
            "fr_thiophene",
            "fr_thiazole",
            "fr_oxazole",
            "fr_pyridine",
            "fr_piperdine",
            "fr_piperzine",
            "fr_morpholine",
            "fr_lactam",
            "fr_lactone",
            "fr_tetrazole",
            "fr_epoxide",
            "fr_unbrch_alkane",
            "fr_bicyclic",
            "fr_benzene",
            "fr_phos_acid",
            "fr_phos_ester",
            "fr_nitro_arom",
            "fr_nitro_arom_nonortho",
            "fr_dihydropyridine",
            "fr_phenol_noOrthoHbond",
            "fr_Al_OH_noTert",
            "fr_benzodiazepine",
            "fr_para_hydroxylation",
            "fr_allylic_oxid",
            "fr_aryl_methyl",
            "fr_Ndealkylation1",
            "fr_Ndealkylation2",
            "fr_alkyl_carbamate",
            "fr_ketone_Topliss",
            "fr_ArN",
            "fr_HOCCN"
        ]
    );
    for (smiles, unique, all) in FNS_CASES {
        let m = mol_from_smiles(smiles).unwrap();
        let counts: Vec<u32> = fns.iter().map(|(_, f)| f(&m).unwrap()).collect();
        assert_eq!(counts, unique, "{smiles}");
        let index = |name: &str| names.iter().position(|n| *n == name).unwrap();
        for (name, f) in [
            (
                "fr_benzene",
                fragments::fr_benzene_with_count_unique as fn(&Mol, bool) -> Result<u32>,
            ),
            ("fr_ester", fragments::fr_ester_with_count_unique),
            ("fr_amide", fragments::fr_amide_with_count_unique),
            ("fr_ArN", fragments::fr_arn_with_count_unique),
        ] {
            assert_eq!(f(&m, false).unwrap(), all[index(name)], "{smiles} {name}");
            assert_eq!(f(&m, true).unwrap(), unique[index(name)], "{smiles} {name}");
        }
    }
}
