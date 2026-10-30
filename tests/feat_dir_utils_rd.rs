// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.Features.FeatDirUtilsRD in RDKit
// 2024.03.6, for every heavy atom of a molecule in tests/data, and the tests
// of UnitTestFeatDirUtilsRD.py.

// Coordinates such as 0.3927 are not approximations of π/8.
#![allow(clippy::approx_constant)]

use rdkit_cxx::chem::features::feat_dir_utils_rd::*;
use rdkit_cxx::chem::rdmolfiles::{MolFromMolFileParams, mol_from_mol_file_with_params};
use rdkit_cxx::chem::{Conformer, Mol, RWMol, add_hs, mol_from_smiles};
use rdkit_cxx::geometry::rd_geometry::Point3D;

fn p(x: f64, y: f64, z: f64) -> Point3D {
    Point3D::new(x, y, z)
}

/// NumPy's matrix product (in the rotations) may round differently.
fn assert_close(a: Point3D, b: Point3D) {
    let d = ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs());
    assert!(d <= 1e-12, "{a:?} != {b:?}");
}

fn assert_vects(actual: FeatVects, kind: &str, expected: Vec<(Point3D, Point3D)>) {
    assert_eq!(actual.1, kind);
    assert_eq!(actual.0.len(), expected.len(), "{actual:?}");
    for ((a0, a1), (e0, e1)) in actual.0.into_iter().zip(expected) {
        assert_close(a0, e0);
        assert_close(a1, e1);
    }
}

fn molecule() -> Mol {
    let path = format!("{}/tests/data/feat_dirs.mol", env!("CARGO_MANIFEST_DIR"));
    let params = MolFromMolFileParams {
        remove_hs: false,
        ..Default::default()
    };
    mol_from_mol_file_with_params(path, &params).unwrap()
}

#[test]
fn feat_vects_match_python() {
    let mol = molecule();
    let conf = mol.get_conformer().unwrap();

    assert_vects(
        get_acceptor1_feat_vects(&conf, &[0]).unwrap(),
        "cone",
        vec![(
            p(7.6856, -0.885, 1.3033),
            p(8.065798854319553, -1.556883945145039, 2.5893873980528843),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[1]).unwrap(),
        "cone",
        vec![(
            p(7.3613, -0.3119, 0.2063),
            p(6.981101145680446, 0.3599839451450389, -1.0797873980528843),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[2]).unwrap(),
        "linear",
        vec![
            (
                p(8.4168, 0.2289, -0.5967),
                p(8.097903458202474, 0.9123305618212719, -1.8933216344015855),
            ),
            (
                p(8.4168, 0.2289, -0.5967),
                p(9.8411213376837, 0.11184909674361979, -0.14105999619680204),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[3]).unwrap(),
        "linear",
        vec![
            (
                p(5.9647, -0.2134, -0.2106),
                p(5.6134826221076, 0.4641046772373866, -1.5019689502957134),
            ),
            (
                p(5.9647, -0.2134, -0.2106),
                p(4.881861254600301, -0.7897629558108712, 0.6526878978823693),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[4]).unwrap(),
        "cone",
        vec![(
            p(5.5768, 0.3927, -1.3977),
            p(5.157652579897908, 1.0476245973804559, -2.6804272554864528),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[5]).unwrap(),
        "cone",
        vec![(
            p(4.236, 0.4482, -1.7047),
            p(2.775027011951434, 0.5086743442994446, -2.039215742341072),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[6]).unwrap(),
        "cone",
        vec![(
            p(3.2507, -0.0785, -0.8748),
            p(
                2.1887621985131407,
                -0.6461673500894439,
                0.019650605352628747,
            ),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[7]).unwrap(),
        "cone",
        vec![(
            p(3.6783, -0.6832, 0.3126),
            p(4.13662732307469, -1.3313537237213877, 1.5853265281077815),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[8]).unwrap(),
        "cone",
        vec![(
            p(5.0263, -0.7522, 0.6464),
            p(6.480529453706842, -0.8266375610576944, 1.0065051866820058),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[9]).unwrap(),
        "linear",
        vec![
            (
                p(1.8736, 0.0175, -1.2139),
                p(0.8340295840223051, -0.5597336341909703, -0.2995178688411242),
            ),
            (
                p(1.8736, 0.0175, -1.2139),
                p(1.460003264475843, 0.6960364056267249, -2.48611306651163),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[10]).unwrap(),
        "linear",
        vec![
            (
                p(1.4637, 0.7027, -2.4388),
                p(0.009262473932270865, 0.8013157577296064, -2.79219272077025),
            ),
            (
                p(1.4637, 0.7027, -2.4388),
                p(2.497628054457402, 1.3070193238693117, -3.3420114553104776),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[11]).unwrap(),
        "linear",
        vec![
            (
                p(0.7974, -0.5188, -0.4109),
                p(
                    -0.6246631366126646,
                    -0.3620925360609781,
                    -0.8616540418378413,
                ),
            ),
            (
                p(0.7974, -0.5188, -0.4109),
                p(1.1029975719280323, -1.2318729000478945, 0.8728979448686529),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[12]).unwrap(),
        "linear",
        vec![
            (
                p(0.0694, 0.6694, 0.1818),
                p(-1.1145299256047378, 0.45324866268629377, 1.0771147662331144),
            ),
            (
                p(0.0694, 0.6694, 0.1818),
                p(0.5322034015857405, 2.0625328425875695, -0.12641079537314223),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[13]).unwrap(),
        "linear",
        vec![
            (
                p(-1.0099, 0.3082, 0.9837),
                p(-1.838553491978018, 1.3676245485836003, 1.6477429324155721),
            ),
            (
                p(-1.0099, 0.3082, 0.9837),
                p(-1.3440652182512718, -1.1403750332818485, 1.1836099293771345),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[14]).unwrap(),
        "linear",
        vec![
            (
                p(-1.7495, 1.2879, 1.6234),
                p(-2.9132728826167718, 0.904147690922289, 2.4884819862668364),
            ),
            (
                p(-1.7495, 1.2879, 1.6234),
                p(-1.3871953085541886, 2.7333050766415155, 1.4515294819485003),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[15]).unwrap(),
        "cone",
        vec![(
            p(-1.4674, 2.5152, 1.5142),
            p(-1.1326372184619915, 3.971613902097121, 1.3846144071465774),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[16]).unwrap(),
        "linear",
        vec![
            (
                p(-2.9088, 0.8352, 2.4685),
                p(-3.751427279287107, 1.8456264145722994, 3.1889288507103086),
            ),
            (
                p(-2.9088, 0.8352, 2.4685),
                p(-3.222110545841037, -0.6266135045710152, 2.590720209943358),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[17]).unwrap(),
        "linear",
        vec![
            (
                p(-3.4718, -0.7782, 2.1463),
                p(-4.623072124663766, -1.3405665671099871, 2.926246369418504),
            ),
            (
                p(-3.4718, -0.7782, 2.1463),
                p(-2.8061742482599046, -1.60755964616719, 1.0884224381196725),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[18]).unwrap(),
        "cone",
        vec![(
            p(-4.5918, -1.0393, 3.1699),
            p(-5.683001261938577, -1.2936862941894303, 4.167180010464578),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[19]).unwrap(),
        "cone",
        vec![(
            p(-2.3879, -1.7833, 2.5431),
            p(-1.3256194189891346, -2.768352322145974, 2.9319854456546834),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[20]).unwrap(),
        "cone",
        vec![(
            p(-4.1461, -1.138, 0.6412),
            p(-4.745287773224306, -1.4577208376184274, -0.6962425589202206),
        )],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[21]).unwrap(),
        "linear",
        vec![
            (
                p(-5.5617, -1.5567, 0.7694),
                p(-6.349414355089712, -1.922220066539281, -0.4536703887114263),
            ),
            (
                p(-5.5617, -1.5567, 0.7694),
                p(-6.206993094507241, -1.6150286330554091, 2.122246847483911),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[22]).unwrap(),
        "linear",
        vec![
            (
                p(-6.4687, -0.7246, -0.0853),
                p(
                    -7.944523788859672,
                    -0.9926458206561319,
                    -0.09507661826312214,
                ),
            ),
            (
                p(-6.4687, -0.7246, -0.0853),
                p(-5.900774162761408, 0.37636961949174463, -0.9310695551112295),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[23]).unwrap(),
        "linear",
        vec![
            (
                p(-6.2039, 0.6293, 0.2187),
                p(-6.961916614834916, 1.717480085876177, -0.4822100600909719),
            ),
            (
                p(-6.2039, 0.6293, 0.2187),
                p(-5.164709720521175, 0.9787370763647079, 1.242407620709355),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[24]).unwrap(),
        "linear",
        vec![
            (
                p(-5.0043, 0.9715, -0.3642),
                p(-4.464274389378069, 2.3621805769055237, -0.20805944508190383),
            ),
            (
                p(-5.0043, 0.9715, -0.3642),
                p(-4.237497502650192, -0.04639249875619944, -1.155343976172453),
            ),
        ],
    );
    assert_vects(
        get_acceptor1_feat_vects(&conf, &[25]).unwrap(),
        "linear",
        vec![
            (
                p(-4.0185, -0.1979, -0.4235),
                p(
                    -2.6595441816425547,
                    -0.00332638191512169,
                    -1.027966865012306,
                ),
            ),
            (
                p(-4.0185, -0.1979, -0.4235),
                p(-4.411380198747151, -1.538475676892222, 0.12285336912785638),
            ),
        ],
    );
    assert!(get_acceptor2_feat_vects(&conf, &[0]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[1]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[2]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[3]).is_err());
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[4]).unwrap(),
        "linear",
        vec![(
            p(5.5768, 0.3927, -1.3977),
            p(6.598857830700517, 0.975871074713504, -2.3279200214579623),
        )],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[5]).unwrap(),
        "linear",
        vec![(
            p(4.236, 0.4482, -1.7047),
            p(3.8193878865977307, 1.1040104338374421, -2.987800550166058),
        )],
    );
    assert!(get_acceptor2_feat_vects(&conf, &[6]).is_err());
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[7]).unwrap(),
        "linear",
        vec![(
            p(3.6783, -0.6832, 0.3126),
            p(2.660411298147535, -1.269582117362117, 1.24536931932768),
        )],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[8]).unwrap(),
        "linear",
        vec![(
            p(5.0263, -0.7522, 0.6464),
            p(5.456783977285989, -1.406068927388016, 1.9259073157651274),
        )],
    );
    assert!(get_acceptor2_feat_vects(&conf, &[9]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[10]).is_err());
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[11]).unwrap(),
        "linear",
        vec![(
            p(0.7974, -0.5188, -0.4109),
            p(
                0.4669944529321134,
                -1.9674488444494114,
                -0.20535511453175148,
            ),
        )],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[12]).unwrap(),
        "linear",
        vec![(
            p(0.0694, 0.6694, 0.1818),
            p(0.47079044950862514, 2.092631931808385, -0.06978810648406478),
        )],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[13]).unwrap(),
        "linear",
        vec![
            (
                p(-1.0099, 0.3082, 0.9837),
                p(-1.962156849756822, -0.4602296081454991, 0.11610399412136063),
            ),
            (
                p(-1.0099, 0.3082, 0.9837),
                p(-0.4775689227526906, -0.6023727072008647, 2.0502275285231075),
            ),
        ],
    );
    assert!(get_acceptor2_feat_vects(&conf, &[14]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[15]).is_err());
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[16]).unwrap(),
        "linear",
        vec![(
            p(-2.9088, 0.8352, 2.4685),
            p(-3.533641345807814, 1.7117780338875894, 3.5130976465007224),
        )],
    );
    assert!(get_acceptor2_feat_vects(&conf, &[17]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[18]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[19]).is_err());
    assert!(get_acceptor2_feat_vects(&conf, &[20]).is_err());
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[21]).unwrap(),
        "linear",
        vec![(
            p(-5.5617, -1.5567, 0.7694),
            p(-6.030242593341148, -2.678087916090265, 1.6485797198934264),
        )],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[22]).unwrap(),
        "linear",
        vec![(
            p(-6.4687, -0.7246, -0.0853),
            p(-7.470783680228728, -1.234202064940725, -1.0783428154049535),
        )],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[23]).unwrap(),
        "linear",
        vec![
            (
                p(-6.2039, 0.6293, 0.2187),
                p(-7.297155188187929, 1.5029877183411575, -0.32117300758603995),
            ),
            (
                p(-6.2039, 0.6293, 0.2187),
                p(-6.117375800146415, 0.817441879094312, 1.7043366299235343),
            ),
        ],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[24]).unwrap(),
        "linear",
        vec![(
            p(-5.0043, 0.9715, -0.3642),
            p(-4.69639323260433, 2.3438913321117454, -0.8854825089493192),
        )],
    );
    assert_vects(
        get_acceptor2_feat_vects(&conf, &[25]).unwrap(),
        "linear",
        vec![(
            p(-4.0185, -0.1979, -0.4235),
            p(-2.999451482071919, -0.3434564347811764, -1.5145331994959519),
        )],
    );
    assert!(get_acceptor3_feat_vects(&conf, &[0]).is_err());
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[1]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[2]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[3]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[4]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[5]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[6]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[7]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[8]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[9]).unwrap(),
        "linear",
        vec![],
    );
    assert!(get_acceptor3_feat_vects(&conf, &[10]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[11]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[12]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[13]).is_err());
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[14]).unwrap(),
        "linear",
        vec![],
    );
    assert!(get_acceptor3_feat_vects(&conf, &[15]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[16]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[17]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[18]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[19]).is_err());
    assert_vects(
        get_acceptor3_feat_vects(&conf, &[20]).unwrap(),
        "linear",
        vec![(
            p(-4.1461, -1.138, 0.6412),
            p(-3.256898379409615, -2.2624769541402454, 0.19975684675497646),
        )],
    );
    assert!(get_acceptor3_feat_vects(&conf, &[21]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[22]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[23]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[24]).is_err());
    assert!(get_acceptor3_feat_vects(&conf, &[25]).is_err());
    assert_vects(
        get_donor1_feat_vects(&conf, &[0]).unwrap(),
        "cone",
        vec![(
            p(7.6856, -0.885, 1.3033),
            p(8.065798854319553, -1.556883945145039, 2.5893873980528843),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[1]).unwrap(),
        "cone",
        vec![(
            p(7.3613, -0.3119, 0.2063),
            p(6.981101145680446, 0.3599839451450389, -1.0797873980528843),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[2]).unwrap(),
        "cone",
        vec![(
            p(8.4168, 0.2289, -0.5967),
            p(9.522224795886176, 0.7952796585648919, -1.4376816305983882),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[3]).unwrap(),
        "cone",
        vec![(
            p(5.9647, -0.2134, -0.2106),
            p(
                4.5306438767079005,
                -0.11225827857348435,
                -0.6386810524133439,
            ),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[4]).unwrap(),
        "cone",
        vec![(
            p(5.5768, 0.3927, -1.3977),
            p(5.157652579897908, 1.0476245973804559, -2.6804272554864528),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[5]).unwrap(),
        "cone",
        vec![(
            p(4.236, 0.4482, -1.7047),
            p(2.775027011951434, 0.5086743442994446, -2.039215742341072),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[6]).unwrap(),
        "cone",
        vec![(
            p(3.2507, -0.0785, -0.8748),
            p(
                2.1887621985131407,
                -0.6461673500894439,
                0.019650605352628747,
            ),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[7]).unwrap(),
        "cone",
        vec![(
            p(3.6783, -0.6832, 0.3126),
            p(4.13662732307469, -1.3313537237213877, 1.5853265281077815),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[8]).unwrap(),
        "cone",
        vec![(
            p(5.0263, -0.7522, 0.6464),
            p(6.480529453706842, -0.8266375610576944, 1.0065051866820058),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[9]).unwrap(),
        "cone",
        vec![(
            p(1.8736, 0.0175, -1.2139),
            p(
                0.42043284849814744,
                0.11880277143575473,
                -1.5717309353527544,
            ),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[10]).unwrap(),
        "cone",
        vec![(
            p(1.4637, 0.7027, -2.4388),
            p(1.0431905283896725, 1.405635081598918, -3.695404176080728),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[11]).unwrap(),
        "cone",
        vec![(
            p(0.7974, -0.5188, -0.4109),
            p(-0.31906556468463265, -1.0751654361088723, 0.422143903030812),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[12]).unwrap(),
        "cone",
        vec![(
            p(0.0694, 0.6694, 0.1818),
            p(-0.6517265240189976, 1.8463815052738641, 0.7689039708599723),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[13]).unwrap(),
        "cone",
        vec![(
            p(-1.0099, 0.3082, 0.9837),
            p(-2.17271871022929, -0.08095048469824867, 1.8476528617927066),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[14]).unwrap(),
        "cone",
        vec![(
            p(-1.7495, 1.2879, 1.6234),
            p(-2.5509681911709605, 2.3495527675638046, 2.316611468215337),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[15]).unwrap(),
        "cone",
        vec![(
            p(-1.4674, 2.5152, 1.5142),
            p(-1.1326372184619915, 3.971613902097121, 1.3846144071465774),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[16]).unwrap(),
        "cone",
        vec![(
            p(-2.9088, 0.8352, 2.4685),
            p(-4.0647378251281445, 0.3838129100012843, 3.3111490606536664),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[17]).unwrap(),
        "cone",
        vec![(
            p(-3.4718, -0.7782, 2.1463),
            p(-3.957446372923671, -2.1699262132771775, 1.8683688075381761),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[18]).unwrap(),
        "cone",
        vec![(
            p(-4.5918, -1.0393, 3.1699),
            p(-5.683001261938577, -1.2936862941894303, 4.167180010464578),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[19]).unwrap(),
        "cone",
        vec![(
            p(-2.3879, -1.7833, 2.5431),
            p(-1.3256194189891346, -2.768352322145974, 2.9319854456546834),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[20]).unwrap(),
        "cone",
        vec![(
            p(-4.1461, -1.138, 0.6412),
            p(-4.745287773224306, -1.4577208376184274, -0.6962425589202206),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[21]).unwrap(),
        "cone",
        vec![(
            p(-5.5617, -1.5567, 0.7694),
            p(-6.9947074495969535, -1.9805486995946906, 0.8991764587724846),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[22]).unwrap(),
        "cone",
        vec![(
            p(-6.4687, -0.7246, -0.0853),
            p(-7.376597951621081, 0.10832379883561305, -0.9408461733743524),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[23]).unwrap(),
        "cone",
        vec![(
            p(-6.2039, 0.6293, 0.2187),
            p(-5.922726335356092, 2.0669171622408857, 0.5414975606183834),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[24]).unwrap(),
        "cone",
        vec![(
            p(-5.0043, 0.9715, -0.3642),
            p(-3.6974718920282608, 1.344288078149324, -0.999203421254357),
        )],
    );
    assert_vects(
        get_donor1_feat_vects(&conf, &[25]).unwrap(),
        "cone",
        vec![(
            p(-4.0185, -0.1979, -0.4235),
            p(-3.052424380389705, -1.3439020588073445, -0.4816134958844497),
        )],
    );
    assert!(get_donor2_feat_vects(&conf, &[0]).is_err());
    assert!(get_donor2_feat_vects(&conf, &[1]).is_err());
    assert!(get_donor2_feat_vects(&conf, &[2]).is_err());
    assert!(get_donor2_feat_vects(&conf, &[3]).is_err());
    assert_vects(
        get_donor2_feat_vects(&conf, &[4]).unwrap(),
        "linear",
        vec![(
            p(5.5768, 0.3927, -1.3977),
            p(6.658080850824108, 0.94946345358574, -2.275678461234883),
        )],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[5]).unwrap(),
        "linear",
        vec![(
            p(4.236, 0.4482, -1.7047),
            p(3.835816836967586, 1.106130561293414, -2.9919376674653),
        )],
    );
    assert!(get_donor2_feat_vects(&conf, &[6]).is_err());
    assert_vects(
        get_donor2_feat_vects(&conf, &[7]).unwrap(),
        "linear",
        vec![(
            p(3.6783, -0.6832, 0.3126),
            p(2.7177723820224564, -1.2846778853808958, 1.2952551014980345),
        )],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[8]).unwrap(),
        "linear",
        vec![(
            p(5.0263, -0.7522, 0.6464),
            p(5.411067390248432, -1.4153740804781898, 1.935645591182424),
        )],
    );
    assert!(get_donor2_feat_vects(&conf, &[9]).is_err());
    assert_vects(
        get_donor2_feat_vects(&conf, &[10]).unwrap(),
        "linear",
        vec![
            (
                p(1.4637, 0.7027, -2.4388),
                p(1.9723881663883418, 2.113689558247439, -2.457369218970211),
            ),
            (
                p(1.4637, 0.7027, -2.4388),
                p(2.1100302256486554, 0.05144814061003988, -3.6254458001670415),
            ),
            (
                p(1.4637, 0.7027, -2.4388),
                p(
                    -0.013856344626983264,
                    0.5967848803577922,
                    -2.6746161060736178,
                ),
            ),
        ],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[11]).unwrap(),
        "linear",
        vec![
            (
                p(0.7974, -0.5188, -0.4109),
                p(1.2334484342188072, -1.5602562000539866, 0.5766377189679579),
            ),
            (
                p(0.7974, -0.5188, -0.4109),
                p(
                    -0.08072882357092304,
                    -1.2399968203261649,
                    -1.3900654178765588,
                ),
            ),
        ],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[12]).unwrap(),
        "linear",
        vec![
            (
                p(0.0694, 0.6694, 0.1818),
                p(1.0685339736434496, 1.3318620529191263, 1.0833959910922561),
            ),
            (
                p(0.0694, 0.6694, 0.1818),
                p(
                    -0.38420083633518465,
                    1.5616508135350398,
                    -0.9353995197913756,
                ),
            ),
        ],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[13]).unwrap(),
        "linear",
        vec![(
            p(-1.0099, 0.3082, 0.9837),
            p(-1.37008724002419, -1.1379327483690003, 1.153891733696083),
        )],
    );
    assert!(get_donor2_feat_vects(&conf, &[14]).is_err());
    assert!(get_donor2_feat_vects(&conf, &[15]).is_err());
    assert_vects(
        get_donor2_feat_vects(&conf, &[16]).unwrap(),
        "linear",
        vec![
            (
                p(-2.9088, 0.8352, 2.4685),
                p(-2.511264409949158, 1.018319409152398, 3.9032239234900237),
            ),
            (
                p(-2.9088, 0.8352, 2.4685),
                p(-4.026972910484066, 1.8129575001493397, 2.259503317822251),
            ),
        ],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[17]).unwrap(),
        "linear",
        vec![],
    );
    assert!(get_donor2_feat_vects(&conf, &[18]).is_err());
    assert!(get_donor2_feat_vects(&conf, &[19]).is_err());
    assert!(get_donor2_feat_vects(&conf, &[20]).is_err());
    assert_vects(
        get_donor2_feat_vects(&conf, &[21]).unwrap(),
        "linear",
        vec![
            (
                p(-5.5617, -1.5567, 0.7694),
                p(-5.572748952617161, -2.965037228714254, 0.2531983234592129),
            ),
            (
                p(-5.5617, -1.5567, 0.7694),
                p(-5.983520676256519, -1.5882990948955482, 2.2085208476998295),
            ),
        ],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[22]).unwrap(),
        "linear",
        vec![
            (
                p(-6.4687, -0.7246, -0.0853),
                p(-6.127396889494448, -1.0072818305178866, -1.5183398352642858),
            ),
            (
                p(-6.4687, -0.7246, -0.0853),
                p(-7.923166806687549, -0.9297837102291364, 0.21871637011965306),
            ),
        ],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[23]).unwrap(),
        "linear",
        vec![(
            p(-6.2039, 0.6293, 0.2187),
            p(-7.077246713512679, 1.5360039947945077, 1.0342693617469232),
        )],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[24]).unwrap(),
        "linear",
        vec![
            (
                p(-5.0043, 0.9715, -0.3642),
                p(-4.407826694935622, 2.14700746370385, 0.35162246340983155),
            ),
            (
                p(-5.0043, 0.9715, -0.3642),
                p(-5.350217734363152, 1.3823267829164219, -1.7647578443932996),
            ),
        ],
    );
    assert_vects(
        get_donor2_feat_vects(&conf, &[25]).unwrap(),
        "linear",
        vec![
            (
                p(-4.0185, -0.1979, -0.4235),
                p(-4.269005857772371, -0.9003547680124745, -1.72496229838527),
            ),
            (
                p(-4.0185, -0.1979, -0.4235),
                p(-2.664339115319325, 0.41730302287095467, -0.6178541588217014),
            ),
        ],
    );
    assert!(get_donor3_feat_vects(&conf, &[0]).is_err());
    assert_vects(
        get_donor3_feat_vects(&conf, &[1]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[2]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[3]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[4]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[5]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[6]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[7]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[8]).unwrap(),
        "linear",
        vec![],
    );
    assert_vects(
        get_donor3_feat_vects(&conf, &[9]).unwrap(),
        "linear",
        vec![],
    );
    assert!(get_donor3_feat_vects(&conf, &[10]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[11]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[12]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[13]).is_err());
    assert_vects(
        get_donor3_feat_vects(&conf, &[14]).unwrap(),
        "linear",
        vec![],
    );
    assert!(get_donor3_feat_vects(&conf, &[15]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[16]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[17]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[18]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[19]).is_err());
    assert_vects(
        get_donor3_feat_vects(&conf, &[20]).unwrap(),
        "linear",
        vec![(
            p(-4.1461, -1.138, 0.6412),
            p(-3.256898379409615, -2.2624769541402454, 0.19975684675497646),
        )],
    );
    assert!(get_donor3_feat_vects(&conf, &[21]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[22]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[23]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[24]).is_err());
    assert!(get_donor3_feat_vects(&conf, &[25]).is_err());
    assert_vects(
        get_aromatic_feat_vects(
            &conf,
            &[3, 4, 5, 6, 7, 8],
            p(4.622133333333333, -0.14773333333333336, -0.5381333333333332),
        )
        .unwrap(),
        "linear",
        vec![
            (
                p(4.622133333333333, -0.14773333333333336, -0.5381333333333332),
                p(4.513806915106462, 1.170090641330508, 0.170107698318556),
            ),
            (
                p(4.622133333333333, -0.14773333333333336, -0.5381333333333332),
                p(4.730459751560204, -1.4655573079971747, -1.2463743649852224),
            ),
        ],
    );
    assert_close(
        arb_axis_rotation(30.0, p(0.0, 0.0, 1.0), p(1.0, 2.0, 3.0)),
        p(1.8660254037844386, 1.2320508075688774, 3.0),
    );
    let rotated = arb_axis_rotation_list(
        -75.5,
        p(0.6, 0.0, 0.8),
        &[p(1.0, 0.0, 0.0), p(0.3, -2.0, 5.0)],
    );
    assert_close(
        rotated[0],
        p(0.5202432025948425, 0.7745181123024862, 0.3598175980538681),
    );
    assert_close(
        rotated[1],
        p(3.5041971756527657, -3.1728474955524604, 2.5968521182604265),
    );
    assert_eq!(
        cross([1.5, 2.0, 3.0], [-4.0, 5.0, 6.25]),
        [-2.5, -21.375, 15.5]
    );
    let adj = rdkit_cxx::chem::get_adjacency_matrix(&mol).unwrap();
    assert_eq!(find_neighbors(1, &adj).unwrap(), [0, 2, 3]);
}

/// UnitTestFeatDirUtilsRD.py's molecule.
fn donor_molecule() -> Mol {
    let mut emol = RWMol::from(add_hs(&mol_from_smiles("C=CCONC").unwrap()).unwrap());
    let coords = [
        [-2.8272, -0.2716, 0.4130],
        [-1.7908, -0.1146, -0.4177],
        [-0.5452, 0.6287, -0.0653],
        [0.5603, -0.2584, -0.1671],
        [1.7601, 0.4902, 0.1811],
        [2.8427, -0.4743, 0.0559],
        [-3.7065, -0.8216, 0.0959],
        [-2.8190, 0.1408, 1.4159],
        [-1.8508, -0.5462, -1.4133],
        [-0.6006, 1.0355, 0.9521],
        [-0.4226, 1.4609, -0.7693],
        [1.8283, 1.1371, -0.6054],
        [2.8648, -0.9271, -0.9408],
        [2.7437, -1.2668, 0.8044],
        [3.8013, 0.0259, 0.2227],
    ];
    let mut conf = Conformer::new(15);
    for (i, c) in coords.iter().enumerate() {
        conf.set_atom_position(i as u32, p(c[0], c[1], c[2]))
            .unwrap();
    }
    emol.add_conformer(&conf).unwrap();
    emol.into()
}

fn assert_parallel(v: Point3D, reference: Point3D, same_way: bool) {
    let c = v.cross_product(&reference);
    assert!(c.x.abs() < 1e-7 && c.y.abs() < 1e-7 && c.z.abs() < 1e-7);
    assert_eq!(v.dot_product(&reference) > 0.0, same_way);
    assert!((v.length() - 1.5).abs() < 1e-7);
}

#[test]
fn donor2_matches_python_tests() {
    let mol = donor_molecule();
    let conf = mol.get_conformer().unwrap();
    let pos = |i: u32| conf.get_atom_position(i).unwrap();
    let unit = |i: u32, from: u32| {
        let mut v = pos(i) - pos(from);
        v.normalize().unwrap();
        v
    };
    // Case 1: two hydrogens.
    let (vects, _) = get_donor2_feat_vects_with_scale(&conf, &[2], 1.5).unwrap();
    assert_eq!(vects.len(), 2);
    assert_close(vects[0].0, pos(2));
    assert_close(vects[1].0, pos(2));
    assert_parallel(vects[0].1 - vects[0].0, pos(9) - pos(2), true);
    assert_parallel(vects[1].1 - vects[1].0, pos(10) - pos(2), true);
    // Case 2.1: one hydrogen, sp2.
    let (vects, _) = get_donor2_feat_vects(&conf, &[1]).unwrap();
    assert_eq!(vects.len(), 1);
    assert_parallel(vects[0].1 - vects[0].0, pos(8) - pos(1), true);
    // Case 2.2: one hydrogen, sp3.
    let (vects, _) = get_donor2_feat_vects(&conf, &[4]).unwrap();
    assert_eq!(vects.len(), 2);
    assert_parallel(vects[0].1 - vects[0].0, pos(11) - pos(4), true);
    assert_parallel(vects[1].1 - vects[1].0, unit(3, 4) + unit(5, 4), false);
    // Case 3: no hydrogens.
    let (vects, _) = get_donor2_feat_vects(&conf, &[3]).unwrap();
    assert_eq!(vects.len(), 1);
    assert_parallel(vects[0].1 - vects[0].0, unit(2, 3) + unit(4, 3), false);
}
