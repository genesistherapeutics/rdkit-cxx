// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.AtomPairs in RDKit 2026.09.1, and the
// doctests of its modules (UnitTestDescriptors.py's regressions test
// rdMolDescriptors' deprecated fingerprints, which are not ported).

use rdkit_cxx::Error;
use rdkit_cxx::chem::atom_pairs::{pairs, sheridan, torsions, utils};
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_mol_descriptors::get_atom_pair_atom_code;

struct Expected {
    smiles: &'static str,
    bit_vect: Vec<u32>,
    patty: Vec<&'static str>,
    bp: Vec<(u64, i32)>,
    bt: Vec<(u64, i32)>,
    ids: Vec<u64>,
    ids3: Vec<u64>,
}

fn expected() -> Vec<Expected> {
    vec![
        Expected {
            smiles: "CCC",
            bit_vect: vec![541730, 558113],
            patty: vec!["HYD", "HYD", "HYD"],
            bp: vec![(98497, 2), (98498, 1)],
            bt: vec![],
            ids: vec![],
            ids3: vec![8405024],
        },
        Expected {
            smiles: "OCC(=O)O",
            bit_vect: vec![
                705601, 1590337, 1590338, 1590625, 1590626, 1592355, 1721410, 1721697, 1723426,
                1723427,
            ],
            patty: vec!["POL", "HYD", "OTH", "ANI", "ANI"],
            bp: vec![
                (32834, 1),
                (49219, 2),
                (98370, 2),
                (98401, 1),
                (114753, 2),
                (114786, 1),
                (114881, 1),
            ],
            bt: vec![(538447363, 2)],
            ids: vec![12895666272, 13969408096],
            ids3: vec![25182250, 25186849, 27284001, 27284064],
        },
        Expected {
            smiles: "OCC(N)O",
            bit_vect: vec![
                574529, 1066050, 1066081, 1590337, 1590338, 1590369, 1590370, 1591330, 1591331,
                1592355,
            ],
            patty: vec!["POL", "HYD", "HYD", "CAT", "POL"],
            bp: vec![
                (49186, 1),
                (49187, 1),
                (49251, 1),
                (98337, 1),
                (98338, 1),
                (98401, 2),
                (98402, 2),
                (98497, 1),
            ],
            bt: vec![(538446850, 1), (538446852, 1)],
            ids: vec![12893307456, 12893569120],
            ids3: vec![16794145, 25182242, 25182753, 25182784],
        },
        Expected {
            smiles: "c1ccccc1CC(=O)NC",
            bit_vect: vec![
                558115, 689189, 689190, 689191, 689218, 689219, 689220, 689473, 689474, 689475,
                705570, 705572, 705601, 705857, 705858, 705859, 705860, 705861, 705890, 1082401,
                1082434, 1082692, 1082693, 1082694, 1082721, 1082723, 1721379, 1721410, 1721668,
                1721669, 1721670, 1721697, 1721699, 1722434,
            ],
            patty: vec![
                "HYD", "HYD", "HYD", "HYD", "HYD", "HYD", "HYD", "OTH", "ACC", "DON", "HYD",
            ],
            bp: vec![
                (82050, 1),
                (98433, 1),
                (98434, 1),
                (98435, 1),
                (98436, 2),
                (98437, 2),
                (98438, 1),
                (98466, 1),
                (98467, 2),
                (98468, 2),
                (98469, 2),
                (98470, 1),
                (98497, 7),
                (98498, 8),
                (98499, 6),
                (98500, 2),
                (98501, 2),
                (98502, 2),
                (98503, 1),
                (114817, 1),
                (114849, 1),
                (114881, 1),
                (114882, 2),
                (114883, 2),
                (114884, 2),
                (114885, 1),
            ],
            bt: vec![
                (940576262, 1),
                (941100039, 8),
                (941100549, 1),
                (941100550, 1),
                (941361159, 1),
                (1075317767, 2),
            ],
            ids: vec![
                4439965728,
                5513433129,
                5513433129,
                5513433633,
                5513433633,
                5513695273,
                5513695273,
                5645554217,
                5645554217,
                5647650857,
                5647650857,
                8734916650,
                13969408042,
                13969424416,
            ],
            ids3: vec![
                10768425, 10768425, 10768425, 10768929, 10768929, 10768937, 11026474, 11030569,
                11030569, 11042848, 17060385, 27284001, 27284033,
            ],
        },
        Expected {
            smiles: "C[C@H](F)Cl",
            bit_vect: vec![574497, 2114594, 2114657, 4211746, 4211809, 4214818],
            patty: vec!["HYD", "HYD", "HYD", "HYD"],
            bp: vec![(98497, 3), (98498, 3)],
            bt: vec![],
            ids: vec![],
            ids3: vec![33571360, 67125792, 67125888],
        },
        Expected {
            smiles: "CC(=O)Nc1ccc(O)cc1",
            bit_vect: vec![
                689188, 689189, 689473, 689474, 689475, 705569, 705571, 705574, 705857, 705858,
                705859, 705860, 705890, 705891, 705893, 1082402, 1082690, 1082691, 1082721,
                1082724, 1590311, 1590594, 1590595, 1590625, 1590628, 1590630, 1591365, 1721378,
                1721668, 1721669, 1721697, 1721699, 1721702, 1722434, 1723431,
            ],
            patty: vec![
                "HYD", "OTH", "ACC", "DON", "OTH", "HYD", "HYD", "OTH", "POL", "HYD", "HYD",
            ],
            bp: vec![
                (65637, 1),
                (82023, 1),
                (82050, 1),
                (98402, 2),
                (98403, 2),
                (98407, 1),
                (98434, 3),
                (98435, 2),
                (98466, 1),
                (98468, 2),
                (98469, 2),
                (98497, 2),
                (98498, 2),
                (98499, 2),
                (98500, 2),
                (98501, 2),
                (114785, 1),
                (114788, 1),
                (114790, 1),
                (114817, 2),
                (114820, 1),
                (114849, 1),
                (114851, 1),
                (114854, 1),
                (114881, 5),
                (114882, 4),
                (114883, 3),
                (114884, 2),
                (114886, 1),
                (114914, 1),
                (114915, 1),
                (114917, 1),
            ],
            bt: vec![
                (941100548, 2),
                (941100549, 2),
                (941362183, 4),
                (1074793990, 1),
                (1074793991, 3),
                (1075317768, 2),
            ],
            ids: vec![
                5513695273,
                5513695273,
                5513695273,
                5513695273,
                5647650858,
                5647650858,
                5653942816,
                5653942825,
                5653942825,
                8734920745,
                8734920745,
                12895670313,
                12895670313,
                13969424426,
            ],
            ids3: vec![
                10768937, 10768937, 11030569, 11030569, 11030569, 11030569, 11042858, 17060384,
                17060393, 17060393, 25186857, 25186857, 27284000, 27284033,
            ],
        },
        Expected {
            smiles: "C1CCCCN1",
            bit_vect: vec![558145, 558146, 558147, 1082433, 1082434, 1082435],
            patty: vec!["HYD", "HYD", "HYD", "HYD", "HYD", "CAT"],
            bp: vec![
                (98337, 2),
                (98338, 2),
                (98339, 1),
                (98497, 4),
                (98498, 4),
                (98499, 2),
            ],
            bt: vec![(941097479, 2), (941100034, 2), (941100039, 2)],
            ids: vec![
                4437590049, 4437590049, 4445978657, 4445978657, 8732557345, 8732557345,
            ],
            ids3: vec![8667169, 8667169, 8667169, 8683553, 17055777, 17055777],
        },
        Expected {
            smiles: "CN1CCC[C@H]1c2cccnc2",
            bit_vect: vec![
                558114, 558115, 558145, 558146, 574498, 574529, 574530, 689188, 689189, 689190,
                689219, 689220, 689221, 689222, 689250, 689251, 689252, 689473, 689474, 689475,
                705571, 705602, 705603, 705633, 705857, 705858, 705859, 1098785, 1098817, 1098818,
                1098849, 1099075, 1099076, 1099077, 1099106, 1213477, 1213508, 1213509, 1213539,
                1213761, 1213762, 1213763, 1213794, 1214564,
            ],
            patty: vec![
                "HYD", "CAT", "HYD", "HYD", "HYD", "HYD", "HYD", "HYD", "HYD", "HYD", "ACC", "HYD",
            ],
            bp: vec![
                (81956, 1),
                (98337, 3),
                (98338, 3),
                (98339, 2),
                (98340, 1),
                (98341, 1),
                (98465, 2),
                (98466, 2),
                (98467, 2),
                (98468, 1),
                (98469, 3),
                (98497, 8),
                (98498, 12),
                (98499, 10),
                (98500, 8),
                (98501, 4),
                (98502, 3),
            ],
            bt: vec![
                (941097479, 6),
                (941099527, 2),
                (941100034, 4),
                (941100038, 3),
                (941100039, 7),
            ],
            ids: vec![
                4437606944, 4437869088, 4446241313, 4571807777, 4580458529, 5513433634, 5513691681,
                5513691681, 5513695273, 5521821737, 5645811745, 5645828640, 5645828641, 5647650857,
                5647667241, 8866775073, 8867037217, 8867041833, 8867041833, 9808400425, 9808400930,
                9808400937,
            ],
            ids3: vec![
                8667169, 8684064, 8929313, 8946208, 8946209, 10768425, 10768930, 10768930,
                10768937, 10784809, 11026977, 11030569, 17317921, 17318433, 17318442, 19157033,
                19157034,
            ],
        },
        Expected {
            smiles: "O=C([O-])CC[NH3+]",
            bit_vect: vec![
                558145, 705601, 705602, 1066049, 1066050, 1066339, 1590338, 1590339, 1590625,
                1591332, 1721410, 1721411, 1721697, 1722404, 1723426,
            ],
            patty: vec!["ANI", "OTH", "ANI", "HYD", "HYD", "CAT"],
            bp: vec![
                (32804, 2),
                (32834, 1),
                (98337, 1),
                (98338, 1),
                (98370, 2),
                (98371, 2),
                (98497, 1),
                (114723, 1),
                (114753, 2),
                (114881, 1),
                (114882, 1),
            ],
            bt: vec![(941100547, 2), (1075317762, 1)],
            ids: vec![8598339626, 12895666209, 13969408033],
            ids3: vec![11026465, 16793633, 25186849, 27284001, 27284064],
        },
        Expected {
            smiles: "c1ccc2[nH]ccc2c1",
            bit_vect: vec![
                689473, 689474, 689475, 689476, 705857, 705858, 705859, 705889, 1213761, 1213762,
                1213763, 1213764, 1213793, 1213794,
            ],
            patty: vec![
                "HYD", "HYD", "HYD", "OTH", "DON", "HYD", "HYD", "HYD", "HYD",
            ],
            bp: vec![
                (98433, 1),
                (98434, 3),
                (98435, 2),
                (98436, 1),
                (98497, 6),
                (98498, 6),
                (98499, 6),
                (98500, 3),
                (114817, 1),
                (114881, 2),
                (114882, 4),
                (114883, 1),
            ],
            bt: vec![
                (941100037, 1),
                (941100039, 4),
                (941100549, 3),
                (941361159, 2),
                (941362183, 3),
                (1074793479, 1),
                (1075317767, 3),
            ],
            ids: vec![
                5513433129, 5513695273, 5513695273, 5513695785, 5513695785, 5521822249, 5647650857,
                5647650857, 5647913001, 5647913001, 5647913001, 5647929385, 5656039465, 9808400426,
                9808662569, 9808663081, 9808663081,
            ],
            ids3: vec![
                10768425, 10768425, 10768937, 11030569, 11030569, 11030569, 11031081, 11031081,
                11031081, 11046953, 19157033, 19157545, 19157546,
            ],
        },
    ]
}

#[test]
fn fingerprints_match_python() {
    for e in expected() {
        let m = mol_from_smiles(e.smiles).unwrap();
        let s = e.smiles;
        let bv = pairs::get_atom_pair_fingerprint_as_bit_vect(&m).unwrap();
        assert_eq!(bv.get_num_bits(), pairs::fp_len(), "{s}");
        assert_eq!(bv.get_on_bits(), e.bit_vect, "{s}");
        assert_eq!(sheridan::assign_patty_types(&m).unwrap(), e.patty, "{s}");
        let bp: Vec<_> = sheridan::get_bp_fingerprint(&m)
            .unwrap()
            .get_nonzero_elements()
            .into_iter()
            .collect();
        assert_eq!(bp, e.bp, "{s}");
        let bt: Vec<_> = sheridan::get_bt_fingerprint(&m)
            .unwrap()
            .get_nonzero_elements()
            .into_iter()
            .collect();
        assert_eq!(bt, e.bt, "{s}");
        assert_eq!(
            torsions::get_topological_torsion_fingerprint_as_ids(&m).unwrap(),
            e.ids,
            "{s}"
        );
        assert_eq!(
            torsions::get_topological_torsion_fingerprint_as_ids_with_target_size(&m, 3).unwrap(),
            e.ids3,
            "{s}"
        );
    }
}

fn explained(symbol: &str, n_branch: u32, n_pi: u32) -> utils::AtomCodeExplanation {
    (symbol.to_owned(), n_branch, n_pi)
}

type Chiral = (String, u32, u32, String);

fn chiral(symbol: &str, n_branch: u32, n_pi: u32, label: &str) -> Chiral {
    (symbol.to_owned(), n_branch, n_pi, label.to_owned())
}

fn chiral_tuple(e: &utils::ChiralAtomCodeExplanation) -> Chiral {
    (e.symbol.clone(), e.n_branch, e.n_pi, e.chirality.clone())
}

fn chiral_pair(
    pair: (
        utils::ChiralAtomCodeExplanation,
        u32,
        utils::ChiralAtomCodeExplanation,
    ),
) -> (Chiral, u32, Chiral) {
    (chiral_tuple(&pair.0), pair.1, chiral_tuple(&pair.2))
}

#[test]
fn module_constants_match_python() {
    assert_eq!(pairs::num_path_bits(), 5);
    assert_eq!(pairs::num_fp_bits(), 23);
    assert_eq!(pairs::fp_len(), 8388608);
    assert_eq!(sheridan::num_path_bits(), 5);
    assert_eq!(sheridan::num_fp_bits(), 23);
    assert_eq!(sheridan::fp_len(), 8388608);
    let typ_map = sheridan::typ_map();
    assert_eq!(typ_map.len(), 7);
    assert_eq!((typ_map["CAT"], typ_map["HYD"], typ_map["OTH"]), (1, 6, 7));
}

#[test]
fn utils_doctests() {
    let m = mol_from_smiles("C=CC(=O)O").unwrap();
    let codes: Vec<u64> = m
        .get_atoms()
        .map(|a| get_atom_pair_atom_code(&a).unwrap().into())
        .collect();
    assert_eq!(codes, [41, 42, 43, 105, 97]);
    let expected = [
        explained("C", 1, 1),
        explained("C", 2, 1),
        explained("C", 3, 1),
        explained("O", 1, 1),
        explained("O", 1, 0),
    ];
    for (code, e) in codes.iter().zip(expected) {
        assert_eq!(utils::explain_atom_code(*code), e);
    }

    let code_with = |smiles: &str, idx: u32, include_chirality: bool| -> u64 {
        let m = mol_from_smiles(smiles).unwrap();
        let params = rdkit_cxx::chem::rd_mol_descriptors::GetAtomPairAtomCodeParams {
            branch_subtract: 0,
            include_chirality,
        };
        rdkit_cxx::chem::rd_mol_descriptors::get_atom_pair_atom_code_with_params(
            &m.get_atom_with_idx(idx).unwrap(),
            &params,
        )
        .unwrap()
        .into()
    };
    assert_eq!(
        utils::explain_atom_code(code_with("C[C@H](F)Cl", 1, false)),
        explained("C", 3, 0)
    );
    let code = code_with("C[C@H](F)Cl", 1, true);
    assert_eq!(
        chiral_tuple(&utils::explain_atom_code_include_chirality(code).unwrap()),
        chiral("C", 3, 0, "R")
    );
    assert_eq!(utils::explain_atom_code(code), explained("C", 3, 0));
    let code = code_with("C[C@H](F)Cl", 0, true);
    assert_eq!(
        chiral_tuple(&utils::explain_atom_code_include_chirality(code).unwrap()),
        chiral("C", 1, 0, "")
    );
    let code = code_with("C[C@@H](F)Cl", 1, true);
    assert_eq!(
        chiral_tuple(&utils::explain_atom_code_include_chirality(code).unwrap()),
        chiral("C", 3, 0, "S")
    );
    // Python's chirality dictionary has no 3, and the last type is a dummy.
    assert!(matches!(
        utils::explain_atom_code_include_chirality(3 << 9),
        Err(Error::Key { key }) if key == "3"
    ));
    assert_eq!(
        utils::explain_atom_code((1 << 20) - 1),
        explained("*", 7, 3)
    );

    assert_eq!(utils::bits_in_common(&[1, 2, 3, 4, 10], &[2, 4, 6]), 2);
    assert_eq!(utils::bits_in_common(&[1, 2, 2, 3, 4], &[2, 2, 4, 5, 6]), 3);

    assert_eq!(utils::dice_similarity(&[1, 2, 3], &[1, 2, 3]), 1.0);
    assert_eq!(utils::dice_similarity(&[1, 2, 3], &[5, 6]), 0.0);
    assert_eq!(utils::dice_similarity(&[1, 2, 3, 4], &[1, 3, 5, 7]), 0.5);
    assert_eq!(utils::dice_similarity(&[1, 2, 3, 4, 5, 6], &[1, 3]), 0.5);
    assert_eq!(utils::dice_similarity(&[1, 1, 3, 4, 5, 6], &[1, 1]), 0.5);
    assert_eq!(utils::dice_similarity(&[1, 1, 3, 4, 5, 6], &[1]), 2.0 / 7.0);
    assert_eq!(utils::dice_similarity::<u32>(&[], &[]), 0.0);
    let two_thirds = utils::dice_similarity(&[1, 1, 3, 4], &[1, 1]);
    assert_eq!(two_thirds, 4.0 / 6.0);
    assert_eq!(
        utils::dice_similarity_with_bounds(&[1, 1, 3, 4], &[1, 1], Some(0.3)),
        two_thirds
    );
    assert_eq!(
        utils::dice_similarity_with_bounds(&[1, 1, 3, 4], &[1, 1], Some(0.33)),
        two_thirds
    );
    assert_eq!(
        utils::dice_similarity_with_bounds(&[1, 1, 3, 4, 5, 6], &[1, 1], Some(0.34)),
        0.0
    );

    assert_eq!(utils::dot(&[1, 2, 3, 4, 10], &[2, 4, 6]), 2);
    assert_eq!(utils::dot(&[1, 2, 2, 3, 4], &[2, 2, 4, 5, 6]), 5);
    assert_eq!(utils::dot(&[1, 2, 2, 3, 4], &[2, 4, 5, 6]), 2);
    assert_eq!(utils::dot(&[1, 2, 2, 3, 4], &[5, 6]), 0);
    assert_eq!(utils::dot(&[], &[5, 6]), 0);

    let cosine = |a: &[u32], b: &[u32]| format!("{:.3}", utils::cosine_similarity(a, b));
    assert_eq!(cosine(&[1, 2, 3, 4, 10], &[2, 4, 6]), "0.516");
    assert_eq!(cosine(&[1, 2, 2, 3, 4], &[2, 2, 4, 5, 6]), "0.714");
    assert_eq!(cosine(&[1, 2, 2, 3, 4], &[1, 2, 2, 3, 4]), "1.000");
    assert_eq!(cosine(&[1, 2, 2, 3, 4], &[5, 6, 7]), "0.000");
    assert_eq!(cosine(&[1, 2, 2, 3, 4], &[]), "0.000");
}

#[test]
fn pairs_doctests() {
    let m = mol_from_smiles("CCCCC").unwrap();
    let atom = |i| m.get_atom_with_idx(i).unwrap();
    let code = |i| u64::from(get_atom_pair_atom_code(&atom(i)).unwrap());
    let (c1, c2, c3) = (code(0), code(1), code(2));
    let t = 1 | c1.min(c2) << 5 | c1.max(c2) << (9 + 5);
    assert_eq!(pairs::py_score_pair(&atom(0), &atom(1), 1).unwrap(), t);
    assert_eq!(pairs::py_score_pair(&atom(1), &atom(0), 1).unwrap(), t);
    let t = 2 | c1.min(c3) << 5 | c1.max(c3) << (9 + 5);
    assert_eq!(pairs::py_score_pair(&atom(0), &atom(2), 2).unwrap(), t);
    let params = pairs::PyScorePairParams {
        atom_codes: Some((c1 as u32, c3 as u32)),
        include_chirality: false,
    };
    assert_eq!(
        pairs::py_score_pair_with_params(&atom(0), &atom(2), 2, &params).unwrap(),
        t
    );
    // Distances wrap at 31, as Python's `% _maxPathLen`.
    let m2 = mol_from_smiles("F[C@@H](Cl)[C@H](F)Cl").unwrap();
    let a1 = m2.get_atom_with_idx(1).unwrap();
    let a3 = m2.get_atom_with_idx(3).unwrap();
    assert_eq!(pairs::py_score_pair(&a1, &a3, 40).unwrap(), 574569);

    let m = mol_from_smiles("C=CC").unwrap();
    let atom = |i| m.get_atom_with_idx(i).unwrap();
    let score =
        |i, j, d| pairs::explain_pair_score(pairs::py_score_pair(&atom(i), &atom(j), d).unwrap());
    assert_eq!(
        score(0, 1, 1),
        (explained("C", 1, 1), 1, explained("C", 2, 1))
    );
    assert_eq!(
        score(0, 2, 2),
        (explained("C", 1, 0), 2, explained("C", 1, 1))
    );
    assert_eq!(
        score(1, 2, 1),
        (explained("C", 1, 0), 1, explained("C", 2, 1))
    );
    assert_eq!(
        score(2, 1, 1),
        (explained("C", 1, 0), 1, explained("C", 2, 1))
    );

    let chiral_params = pairs::PyScorePairParams {
        atom_codes: None,
        include_chirality: true,
    };
    let m = mol_from_smiles("C[C@H](F)Cl").unwrap();
    let atom = |i| m.get_atom_with_idx(i).unwrap();
    let s = pairs::py_score_pair(&atom(0), &atom(1), 1).unwrap();
    assert_eq!(
        pairs::explain_pair_score(s),
        (explained("C", 1, 0), 1, explained("C", 3, 0))
    );
    let s = pairs::py_score_pair_with_params(&atom(0), &atom(1), 1, &chiral_params).unwrap();
    assert_eq!(
        chiral_pair(pairs::explain_pair_score_include_chirality(s).unwrap()),
        (chiral("C", 1, 0, ""), 1, chiral("C", 3, 0, "R"))
    );
    let s = pairs::py_score_pair_with_params(&a1, &a3, 1, &chiral_params).unwrap();
    assert_eq!(s, 69420129);
    assert_eq!(
        chiral_pair(pairs::explain_pair_score_include_chirality(s).unwrap()),
        (chiral("C", 3, 0, "R"), 1, chiral("C", 3, 0, "S"))
    );

    let m = mol_from_smiles("CCC").unwrap();
    let atom = |i| m.get_atom_with_idx(i).unwrap();
    let mut v = vec![
        pairs::py_score_pair(&atom(0), &atom(1), 1).unwrap() as u32,
        pairs::py_score_pair(&atom(0), &atom(2), 2).unwrap() as u32,
    ];
    v.sort();
    assert_eq!(
        pairs::get_atom_pair_fingerprint_as_bit_vect(&m)
            .unwrap()
            .get_on_bits(),
        v
    );
}

#[test]
fn sheridan_doctests() {
    let m = mol_from_smiles("OCC(=O)O").unwrap();
    assert_eq!(
        sheridan::assign_patty_types(&m).unwrap(),
        ["POL", "HYD", "OTH", "ANI", "ANI"]
    );
    let fp = sheridan::get_bp_fingerprint(&m).unwrap();
    assert_eq!(fp.get_total_val().unwrap(), 10);
    let m = mol_from_smiles("OCC(N)O").unwrap();
    assert_eq!(
        sheridan::assign_patty_types(&m).unwrap(),
        ["POL", "HYD", "HYD", "CAT", "POL"]
    );
    assert_eq!(
        sheridan::get_bt_fingerprint(&m)
            .unwrap()
            .get_total_val()
            .unwrap(),
        2
    );

    // Other rules and fingerprint functions, as Python's defns and fpfn.
    let defns = [
        (
            rdkit_cxx::chem::mol_from_smarts("[#6]").unwrap(),
            "HYD".to_owned(),
        ),
        (
            rdkit_cxx::chem::mol_from_smarts("[#8]").unwrap(),
            "ACC".to_owned(),
        ),
    ];
    assert_eq!(
        sheridan::assign_patty_types_with_defns(&m, &defns).unwrap(),
        ["ACC", "HYD", "HYD", "", "ACC"]
    );
    let typs = sheridan::get_bp_fingerprint_with_fpfn(&m, |_, typs| Ok(typs.to_vec())).unwrap();
    assert_eq!(typs, [3, 6, 6, 1, 3]);
    let typs = sheridan::get_bt_fingerprint_with_fpfn(&m, |_, typs| Ok(typs.len())).unwrap();
    assert_eq!(typs, 5);
}

#[test]
fn torsions_doctests() {
    let m = mol_from_smiles("CCCCC").unwrap();
    let code = |i, sub| {
        let params = rdkit_cxx::chem::rd_mol_descriptors::GetAtomPairAtomCodeParams {
            branch_subtract: sub,
            include_chirality: false,
        };
        let atom = m.get_atom_with_idx(i).unwrap();
        u64::from(
            rdkit_cxx::chem::rd_mol_descriptors::get_atom_pair_atom_code_with_params(
                &atom, &params,
            )
            .unwrap(),
        )
    };
    let t = code(0, 1) | code(1, 2) << 9 | code(2, 2) << 18 | code(3, 1) << 27;
    assert_eq!(torsions::py_score_path(&m, &[0, 1, 2, 3], 4).unwrap(), t);
    assert_eq!(torsions::py_score_path(&m, &[3, 2, 1, 0], 4).unwrap(), t);
    let score = torsions::py_score_path(&m, &[0, 1, 2, 3, 4], 5).unwrap();
    assert_eq!(
        torsions::explain_path_score_with_size(score, 5),
        [
            explained("C", 1, 0),
            explained("C", 2, 0),
            explained("C", 2, 0),
            explained("C", 2, 0),
            explained("C", 1, 0)
        ]
    );
    assert!(matches!(
        torsions::py_score_path(&m, &[0, 1, 2, 3, 4], 8),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        torsions::py_score_path(&m, &[0, 1, 2], 4),
        Err(Error::Value { .. })
    ));

    let m = mol_from_smiles("C=CC(=O)O").unwrap();
    assert_eq!(
        torsions::py_score_path(&m, &[0, 1, 2, 4], 4).unwrap(),
        12895670312
    );
    assert_eq!(
        torsions::py_score_path(&m, &[0, 1, 2], 3).unwrap(),
        11030568
    );
    let codes: Vec<u32> = m
        .get_atoms()
        .map(|a| get_atom_pair_atom_code(&a).unwrap())
        .collect();
    assert_eq!(
        torsions::py_score_path_with_atom_codes(&m, &[0, 1, 2, 4], 4, Some(&codes)).unwrap(),
        12895670312
    );

    let explain = |smiles: &str, path: &[u32]| {
        let m = mol_from_smiles(smiles).unwrap();
        let score = torsions::py_score_path(&m, path, path.len()).unwrap();
        torsions::explain_path_score_with_size(score, path.len())
    };
    let c10 = explained("C", 1, 0);
    let c11 = explained("C", 1, 1);
    let c21 = explained("C", 2, 1);
    let o10 = explained("O", 1, 0);
    assert_eq!(
        explain("C=CC", &[0, 1, 2]),
        [c10.clone(), c21.clone(), c11.clone()]
    );
    assert_eq!(
        explain("C=CC", &[2, 1, 0]),
        [c10.clone(), c21.clone(), c11.clone()]
    );
    assert_eq!(
        explain("C=CO", &[0, 1, 2]),
        [c11.clone(), c21.clone(), o10.clone()]
    );
    assert_eq!(
        explain("OC=CO", &[0, 1, 2, 3]),
        [o10.clone(), c21.clone(), c21.clone(), o10.clone()]
    );
    assert_eq!(
        explain("CC=CO", &[0, 1, 2, 3]),
        [c10.clone(), c21.clone(), c21.clone(), o10.clone()]
    );
    let c31 = explained("C", 3, 1);
    assert_eq!(
        explain("C=CC(=O)O", &[0, 1, 2, 3]),
        [c11.clone(), c21.clone(), c31.clone(), explained("O", 1, 1)]
    );
    assert_eq!(
        explain("C=CC(=O)O", &[0, 1, 2, 4]),
        [c11.clone(), c21.clone(), c31, o10.clone()]
    );
    let o20 = explained("O", 2, 0);
    assert_eq!(
        explain("OOOO", &[0, 1, 2]),
        [o10.clone(), o20.clone(), o20.clone()]
    );
    let score =
        torsions::py_score_path(&mol_from_smiles("OOOO").unwrap(), &[0, 1, 2, 3], 4).unwrap();
    assert_eq!(
        torsions::explain_path_score(score),
        [o10.clone(), o20.clone(), o20, o10]
    );

    // UnitTestDescriptors.testGetTopologicalTorsionFingerprintAsIds.
    let m = mol_from_smiles("C1CCCCN1").unwrap();
    assert_eq!(
        torsions::get_topological_torsion_fingerprint_as_ids(&m).unwrap(),
        [
            4437590049, 4437590049, 4445978657, 4445978657, 8732557345, 8732557345
        ]
    );
}
