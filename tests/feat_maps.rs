// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.FeatMaps in RDKit 2025.09.6, and the tests
// of UnitTestFeatMap.py, UnitTestFeatMapParser.py, UnitTestFeatMapPoint.py
// (its doctests), and UnitTestFeatMapUtils.py. Python sets `featDirs` on
// free features as attributes; here those features are feature-map points.

use std::collections::BTreeMap;

use rdkit_cxx::Error;
use rdkit_cxx::chem::feat_maps::feat_map_parser::FeatMapParser;
use rdkit_cxx::chem::feat_maps::feat_map_point::{FeatMapPoint, FeatureLike};
use rdkit_cxx::chem::feat_maps::feat_map_utils::{
    CombineFeatMapsParams, DirMergeMode, MergeFeatPointsParams, MergeMethod, MergeMetric,
    combine_feat_maps, combine_feat_maps_with_params, families_match, feq, feq_with_tol,
    get_feat_feat_dist_matrix, merge_feat_points, merge_feat_points_with_params,
};
use rdkit_cxx::chem::feat_maps::feat_maps::{
    FeatDirScoreMode, FeatMap, FeatMapInitParams, FeatMapParams, FeatMapScoreMode, FeatProfile,
};
use rdkit_cxx::chem::rd_chemical_features::FreeChemicalFeature;
use rdkit_cxx::geometry::rd_geometry::Point3D;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-5
}

fn free(family: &str, x: f64, y: f64, z: f64) -> FreeChemicalFeature {
    FreeChemicalFeature::new(family, "Foo", Point3D::new(x, y, z))
}

fn with_dirs(family: &str, x: f64, y: f64, z: f64, dirs: &[Point3D]) -> FeatMapPoint {
    let mut p = FeatMapPoint::from(free(family, x, y, z));
    p.feat_dirs = dirs.to_vec();
    p
}

fn params(entries: &[(&str, f64, f64, FeatProfile)]) -> BTreeMap<String, FeatMapParams> {
    entries
        .iter()
        .map(|&(family, radius, width, feat_profile)| {
            (
                family.to_owned(),
                FeatMapParams {
                    radius,
                    width,
                    feat_profile,
                },
            )
        })
        .collect()
}

fn map_with(entries: &[(&str, f64, f64, FeatProfile)]) -> FeatMap {
    FeatMap::new_with_params(&FeatMapInitParams {
        params: params(entries),
        ..Default::default()
    })
    .unwrap()
}

#[test]
fn feat_map_point_doctests() {
    let mut sfeat = FeatMapPoint::from(FreeChemicalFeature::new(
        "Aromatic",
        "Foo",
        Point3D::new(0.0, 0.0, 0.0),
    ));
    let mut fmp = FeatMapPoint::default();
    fmp.init_from_feat(&*sfeat).unwrap();
    assert_eq!(fmp.get_family(), "Aromatic");
    assert_eq!(fmp.get_type(), "Foo");
    assert_eq!(fmp.get_pos(), Point3D::new(0.0, 0.0, 0.0));
    assert!(fmp.feat_dirs.is_empty());
    sfeat.feat_dirs = vec![Point3D::new(1.0, 0.0, 0.0)];
    fmp.init_from_feat(&sfeat).unwrap();
    assert_eq!(fmp.feat_dirs.len(), 1);

    let mut fmp = FeatMapPoint::default();
    fmp.init_from_feat(&*sfeat).unwrap();
    assert_eq!(fmp.get_dist2(&*sfeat).unwrap(), 0.0);
    sfeat.set_pos(Point3D::new(2.0, 0.0, 0.0));
    assert_eq!(fmp.get_dist2(&sfeat).unwrap(), 4.0);

    let mut sfeat = FeatMapPoint::from(free("Aromatic", 0.0, 0.0, 0.0));
    let mut fmp = FeatMapPoint::default();
    fmp.init_from_feat(&sfeat).unwrap();
    assert_eq!(fmp.get_dir_match(&sfeat), 1.0);
    sfeat.feat_dirs = vec![Point3D::new(0.0, 0.0, 1.0), Point3D::new(0.0, 0.0, -1.0)];
    fmp.feat_dirs = vec![Point3D::new(0.0, 0.0, 1.0), Point3D::new(1.0, 0.0, 0.0)];
    assert_eq!(fmp.get_dir_match(&sfeat), 1.0);
    assert_eq!(fmp.get_dir_match_with_use_best(&sfeat, true), 1.0);
    assert_eq!(fmp.get_dir_match_with_use_best(&sfeat, false), 0.0);
    sfeat.feat_dirs = vec![Point3D::new(0.0, 0.0, 1.0)];
    assert_eq!(fmp.get_dir_match_with_use_best(&sfeat, false), 0.5);
    fmp.feat_dirs = vec![Point3D::new(0.0, 0.0, -1.0)];
    assert_eq!(fmp.get_dir_match(&sfeat), -1.0);
    assert_eq!(fmp.get_dir_match_with_use_best(&sfeat, false), -1.0);
}

#[test]
fn feat_map_basics() {
    // UnitTestFeatMap.test1Basics.
    let fs = [
        free("Aromatic", 0.0, 0.0, 0.0),
        free("Acceptor", 2.0, 0.0, 0.0),
        free("Acceptor", 2.1, 0.0, 0.0),
    ];
    let mut fmap = map_with(&[
        ("Aromatic", 2.0, 1.0, FeatProfile::Gaussian),
        ("Acceptor", 1.0, 1.0, FeatProfile::Gaussian),
    ]);
    fmap.add_feature_with_weight(&fs[0], 1.0).unwrap();
    fmap.add_feature_with_weight(&fs[1], 1.0).unwrap();
    fmap.add_feature_with_weight(&fs[2], 2.0).unwrap();
    fmap.add_feature_with_weight(&fs[2], 2.0).unwrap();
    assert_eq!(fmap.get_num_features(), 4);
    assert_eq!(fmap.get_features().len(), 4);
    fmap.drop_feature(3).unwrap();
    assert_eq!(fmap.get_num_features(), 3);
    assert!(matches!(fmap.drop_feature(3), Err(Error::Index { .. })));
    for (i, family, weight) in [
        (0, "Aromatic", 1.0),
        (1, "Acceptor", 1.0),
        (2, "Acceptor", 2.0),
    ] {
        let f = fmap.get_feature(i).unwrap();
        assert_eq!((f.get_family().as_str(), f.weight), (family, weight));
    }
    assert!(matches!(
        fmap.add_feature(&free("Donor", 0.0, 0.0, 0.0)),
        Err(Error::Value { .. })
    ));

    // Python's feats and weights arguments.
    let points: Vec<FeatMapPoint> = fs.iter().cloned().map(FeatMapPoint::from).collect();
    let init = FeatMapInitParams {
        params: BTreeMap::new(),
        feats: points.clone(),
        weights: vec![1.0, 2.0, 3.0],
    };
    let fmap = FeatMap::new_with_params(&init).unwrap();
    assert_eq!(
        fmap.get_features()
            .iter()
            .map(|f| f.weight)
            .collect::<Vec<_>>(),
        [1.0, 2.0, 3.0]
    );
    let init = FeatMapInitParams {
        weights: vec![1.0],
        ..init
    };
    assert!(matches!(
        FeatMap::new_with_params(&init),
        Err(Error::Value { .. })
    ));
}

#[test]
fn feat_feat_scores() {
    // UnitTestFeatMap.test2FeatFeatScoreGauss to test4FeatFeatScoreBox.
    let score = |profile: FeatProfile,
                 radius: f64,
                 width: f64,
                 weight: f64,
                 other: FreeChemicalFeature,
                 type_match: bool| {
        let mut fmap = map_with(&[("Aromatic", radius, width, profile)]);
        fmap.add_feature_with_weight(&free("Aromatic", 0.0, 0.0, 0.0), weight)
            .unwrap();
        fmap.get_feat_feat_score_with_type_match(fmap.get_feature(0).unwrap(), &other, type_match)
            .unwrap()
    };
    let aro = |x: f64| FreeChemicalFeature::new("Aromatic", "", Point3D::new(x, 0.0, 0.0));
    let acc = |x: f64| FreeChemicalFeature::new("Acceptor", "", Point3D::new(x, 0.0, 0.0));
    use FeatProfile::*;
    assert!(close(
        score(Gaussian, 2.0, 1.0, 1.0, aro(1.0), true),
        (-1.0f64).exp()
    ));
    assert!(close(
        score(Gaussian, 2.0, 1.0, 1.0, aro(1.5), true),
        (-2.25f64).exp()
    ));
    assert!(close(score(Gaussian, 2.0, 1.0, 1.0, aro(0.0), true), 1.0));
    assert!(close(score(Gaussian, 2.0, 1.0, 1.0, aro(2.1), true), 0.0));
    assert!(close(score(Gaussian, 2.0, 1.0, 1.0, acc(1.0), true), 0.0));
    assert!(close(
        score(Gaussian, 2.0, 1.0, 1.0, acc(1.0), false),
        (-1.0f64).exp()
    ));
    let mut fmap = map_with(&[("Aromatic", 2.0, 1.0, Gaussian)]);
    fmap.add_feature_with_weight(&free("Aromatic", 0.0, 0.0, 0.0), 1.0)
        .unwrap();
    assert!(matches!(fmap.get_feature(1), Err(Error::Index { .. })));
    for (x, expected) in [(1.0, 0.5), (1.5, 0.25), (0.0, 1.0), (2.1, 0.0)] {
        assert!(close(
            score(Triangle, 3.0, 2.0, 1.0, aro(x), true),
            expected
        ));
    }
    for (x, expected) in [(1.0, 1.1), (1.5, 1.1), (0.0, 1.1), (2.1, 0.0)] {
        assert!(close(score(Box, 2.0, 1.0, 1.1, aro(x), true), expected));
    }
    // A family without parameters is Python's KeyError.
    let mut fmap = FeatMap::default();
    fmap.add_feature_with_weight(&free("Aromatic", 0.0, 0.0, 0.0), 1.0)
        .unwrap();
    assert!(matches!(
        fmap.get_feat_feat_score(fmap.get_feature(0).unwrap(), &aro(0.0)),
        Err(Error::Key { .. })
    ));
}

#[test]
fn score_feats() {
    // UnitTestFeatMap.test6ScoreFeats.
    let mut fmap = map_with(&[
        ("Aromatic", 2.0, 1.0, FeatProfile::Box),
        ("Acceptor", 1.0, 1.0, FeatProfile::Box),
    ]);
    fmap.add_feature_with_weight(&free("Aromatic", 0.0, 0.0, 0.0), 1.1)
        .unwrap();
    fmap.add_feature_with_weight(&free("Acceptor", 2.0, 0.0, 0.0), 1.1)
        .unwrap();
    fmap.add_feature_with_weight(&free("Acceptor", 2.1, 0.0, 0.0), 2.1)
        .unwrap();
    let fs = [
        FreeChemicalFeature::new("Aromatic", "", Point3D::new(0.0, 1.0, 0.0)),
        FreeChemicalFeature::new("Acceptor", "", Point3D::new(1.5, 0.0, 0.0)),
    ];
    let fs: Vec<&dyn FeatureLike> = fs.iter().map(|f| f as &dyn FeatureLike).collect();
    assert!(close(fmap.score_feats(&fs).unwrap(), 4.3));
    for (mode, total, idx) in [
        (FeatMapScoreMode::All, 4.3, vec![vec![0], vec![1, 2]]),
        (FeatMapScoreMode::Closest, 2.1, vec![vec![0], vec![1]]),
        (FeatMapScoreMode::Best, 3.2, vec![vec![0], vec![2]]),
    ] {
        fmap.score_mode = mode;
        let (sc, out) = fmap.score_feats_with_outputs(&fs).unwrap();
        assert!(close(sc, total), "{mode:?}");
        assert!(close(out.map_score_vect.iter().sum(), sc));
        assert!(close(out.feats_score_vect.iter().sum(), sc));
        assert_eq!(out.feats_to_feat_map_idx, idx);
    }

    // test7ScoreFeats.
    let mut fmap = map_with(&[
        ("Aromatic", 2.0, 1.0, FeatProfile::Box),
        ("Acceptor", 1.0, 1.0, FeatProfile::Box),
    ]);
    fmap.add_feature_with_weight(&free("Aromatic", 0.0, 0.0, 0.0), 1.1)
        .unwrap();
    for mode in [
        FeatMapScoreMode::All,
        FeatMapScoreMode::Closest,
        FeatMapScoreMode::Best,
    ] {
        fmap.score_mode = mode;
        let (sc, out) = fmap.score_feats_with_outputs(&fs).unwrap();
        assert!(close(sc, 1.1));
        assert!(close(out.map_score_vect.iter().sum(), sc));
        assert!(close(out.feats_score_vect.iter().sum(), sc));
        assert_eq!(out.feats_to_feat_map_idx, [vec![0], vec![]]);
    }

    // test8ScoreFeatDirs.
    let mut fmap = map_with(&[("Acceptor", 2.0, 1.0, FeatProfile::Box)]);
    fmap.dir_score_mode = FeatDirScoreMode::DotFullRange;
    fmap.add_feature_with_weight(
        &with_dirs("Acceptor", 0.0, 0.0, 0.0, &[Point3D::new(0.0, 1.0, 0.0)]),
        1.1,
    )
    .unwrap();
    let score = |fmap: &FeatMap, dir: Point3D| {
        let f = with_dirs("Acceptor", 1.0, 0.0, 0.0, &[dir]);
        fmap.score_feats(&[&f]).unwrap()
    };
    assert!(close(score(&fmap, Point3D::new(0.0, 1.0, 0.0)), 1.1));
    assert!(close(score(&fmap, Point3D::new(0.0, -1.0, 0.0)), -1.1));
    assert!(close(score(&fmap, Point3D::new(0.0, 0.0, 1.0)), 0.0));
    fmap.dir_score_mode = FeatDirScoreMode::DotPosRange;
    assert!(close(score(&fmap, Point3D::new(0.0, -1.0, 0.0)), 0.0));
    assert!(close(score(&fmap, Point3D::new(0.0, 1.0, 0.0)), 1.1));
}

#[test]
fn score_feats_match_python() {
    let mut p = FeatMapParser::from_data(
        "ScoreMode=All
BeginParams
  family=Aromatic radius=2.5 width=1.0 profile=Gaussian
  family=Acceptor radius=1.5 width=0.7 profile=Triangle
EndParams
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 5.0) weight=1.25 dir=(1, 1, 0)
  family=Aromatic pos=(0.0,1.0,0.0) weight=2.0 dir=(0,0,1) dir=(0,0,-1)
  family=Acceptor pos=(1.0,1.0,2.0) weight=1.25
  family=Aromatic pos=(-0.3,1.2,0.1) weight=0.5
EndPoints",
    );
    let mut fm = p.parse().unwrap();
    assert_eq!(
        fm.to_string(),
        "  1     Acceptor  1.0000  0.0000  5.0000  1.2500\n  2     Aromatic  0.0000  1.0000  0.0000  2.0000\n  3     Acceptor  1.0000  1.0000  2.0000  1.2500\n  4     Aromatic -0.3000  1.2000  0.1000  0.5000\n"
    );
    let fs = [
        FreeChemicalFeature::new("Aromatic", "", Point3D::new(0.1, 1.1, 0.2)),
        FreeChemicalFeature::new("Acceptor", "", Point3D::new(1.2, 0.9, 2.1)),
        FreeChemicalFeature::new("Acceptor", "", Point3D::new(1.0, 0.3, 4.6)),
    ];
    let fs: Vec<&dyn FeatureLike> = fs.iter().map(|f| f as &dyn FeatureLike).collect();
    for (mode, total, msv, fsv, idx) in [
        (
            FeatMapScoreMode::All,
            3.470898147377137,
            [
                0.3571428571428567,
                1.8835290671684974,
                0.8125911173601468,
                0.417635105705636,
            ],
            [2.3011641728741337, 0.8125911173601468, 0.3571428571428567],
            vec![vec![1, 3], vec![0, 2], vec![0, 2]],
        ),
        // Python's Closest pass scores each feature against the last point
        // it matched for any feature, not the closest one.
        (
            FeatMapScoreMode::Closest,
            0.8125911173601468,
            [0.0, 0.0, 0.8125911173601468, 0.0],
            [0.0, 0.8125911173601468, 0.0],
            vec![vec![1], vec![2], vec![0]],
        ),
        (
            FeatMapScoreMode::Best,
            3.0532630416715008,
            [
                0.3571428571428567,
                1.8835290671684974,
                0.8125911173601468,
                0.0,
            ],
            [1.8835290671684974, 0.8125911173601468, 0.3571428571428567],
            vec![vec![1], vec![2], vec![0]],
        ),
    ] {
        fm.score_mode = mode;
        let (sc, out) = fm.score_feats_with_outputs(&fs).unwrap();
        assert_eq!(sc, total, "{mode:?}");
        assert_eq!(out.map_score_vect, msv);
        assert_eq!(out.feats_score_vect, fsv);
        assert_eq!(out.feats_to_feat_map_idx, idx);
    }
}

const PARSER_DATA: &str = "
ScoreMode=Best
DirScoreMode=DotFullRange

BeginParams
  family=Aromatic radius=2.5 width=1.0 profile=Triangle
  family=Acceptor radius=1.5
EndParams

# optional
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 5.0) weight=1.25 dir=(1, 1, 0)
  family=Aromatic pos=(0.0,1.0,0.0) weight=2.0 dir=(0,0,1) dir=(0,0,-1)
  family=Acceptor pos=(1.0,1.0,2.0) weight=1.25
EndPoints
";

#[test]
fn parser_matches_python() {
    let mut p = FeatMapParser::default();
    p.set_data(PARSER_DATA);
    let fm = p.parse().unwrap();
    assert_eq!(fm.score_mode, FeatMapScoreMode::Best);
    assert_eq!(fm.dir_score_mode, FeatDirScoreMode::DotFullRange);
    assert_eq!(fm.get_num_features(), 3);
    let feats = fm.get_features();
    assert_eq!(
        feats.iter().map(|f| f.weight).collect::<Vec<_>>(),
        [1.25, 2.0, 1.25]
    );
    assert_eq!(
        feats.iter().map(|f| f.feat_dirs.len()).collect::<Vec<_>>(),
        [1, 2, 0]
    );
    assert_eq!(
        feats.iter().map(|f| f.get_family()).collect::<Vec<_>>(),
        ["Acceptor", "Aromatic", "Acceptor"]
    );
    assert_eq!(feats[0].get_pos(), Point3D::new(1.0, 0.0, 5.0));
    assert_eq!(
        feats[1].feat_dirs,
        [Point3D::new(0.0, 0.0, 1.0), Point3D::new(0.0, 0.0, -1.0)]
    );
    assert_eq!(
        fm.params["Aromatic"],
        FeatMapParams {
            radius: 2.5,
            width: 1.0,
            feat_profile: FeatProfile::Triangle
        }
    );
    assert_eq!(
        fm.params["Acceptor"],
        FeatMapParams {
            radius: 1.5,
            ..Default::default()
        }
    );

    let families = |mut p: FeatMapParser| -> Vec<String> {
        p.parse()
            .unwrap()
            .get_features()
            .iter()
            .map(|f| f.get_family())
            .collect()
    };
    assert_eq!(
        families(FeatMapParser::from_data(PARSER_DATA)),
        ["Acceptor", "Aromatic", "Acceptor"]
    );
    let mut lines = FeatMapParser::default();
    lines.set_data_lines(PARSER_DATA.split('\n').map(str::to_owned).collect());
    assert_eq!(families(lines), ["Acceptor", "Aromatic", "Acceptor"]);
    let file = FeatMapParser::from_file(&mut PARSER_DATA.as_bytes()).unwrap();
    assert_eq!(families(file), ["Acceptor", "Aromatic", "Acceptor"]);

    // Continued lines join with their backslash, which Python keeps in the
    // value before it.
    let mut p =
        FeatMapParser::from_data("BeginPoints\n  family=Acceptor \\\n pos=(1, 2, 3)\nEndPoints");
    let fm = p.parse().unwrap();
    assert_eq!(fm.get_feature(0).unwrap().get_family(), "Acceptor \\");
    assert_eq!(
        fm.get_feature(0).unwrap().get_pos(),
        Point3D::new(1.0, 2.0, 3.0)
    );

    for data in [
        "scoreMode = typo\nbeginParams\nfamily=Acceptor radius=1.5\nEndParams",
        "dirscoremode = typo\nbeginParams\nfamily=Acceptor radius=1.5\nEndParams",
        "typo = All\nbeginParams\nfamily=Acceptor radius=1.5\nEndParams",
        "beginTypo\nfamily=Acceptor radius=1.5\nEndParams",
        "beginParams\nfamily=Acceptor radius=1.5\nEndTypo",
        "beginParams\ntypo=Acceptor radius=1.5\nEndParams",
        "beginParams\nprofile=Typo\nEndParams",
        "BeginPoints\npos=(1.0, 0.0, 5.0) typo=Acceptor\nEndPoints",
    ] {
        let result = FeatMapParser::from_data(data).parse();
        assert!(
            matches!(result, Err(Error::Parse { .. })),
            "{data}: {result:?}"
        );
    }
    let result =
        FeatMapParser::from_data("BeginPoints\npos=(1.0, 0.0, 5.0, 4.0)\nEndPoints").parse();
    assert!(matches!(result, Err(Error::Value { .. })));
}

const PARAM_TXT: &str = "
BeginParams
  family=Acceptor radius=0.5 profile=Box
EndParams
";

fn parse_points(points: &str) -> FeatMap {
    FeatMapParser::from_data(&format!("{PARAM_TXT}\nBeginPoints\n{points}\nEndPoints"))
        .parse()
        .unwrap()
}

fn xs(fm: &FeatMap) -> Vec<f64> {
    fm.get_features().iter().map(|f| f.get_pos().x).collect()
}

fn assert_xs(fm: &FeatMap, expected: &[f64]) {
    let actual = xs(fm);
    assert_eq!(actual.len(), expected.len(), "{actual:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() < 1e-4, "{actual:?} != {expected:?}");
    }
    for f in fm.get_features() {
        assert_eq!((f.get_pos().y, f.get_pos().z), (0.0, 0.0));
    }
}

fn points(xs_weights: &[(f64, f64)]) -> String {
    xs_weights
        .iter()
        .map(|(x, w)| format!("  family=Acceptor pos=({x:?}, 0.0, 0.0) weight={w:?}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn distance(method: MergeMethod) -> MergeFeatPointsParams<'static> {
    MergeFeatPointsParams {
        merge_metric: MergeMetric::Distance,
        merge_method: method,
        ..Default::default()
    }
}

#[test]
fn merging_matches_python() {
    // UnitTestFeatMapUtils.test1Basics.
    let mut fm = parse_points(&points(&[(1.0, 1.0), (1.1, 1.0), (3.0, 1.0)]));
    assert!(!merge_feat_points(&mut fm).unwrap());
    assert!(
        merge_feat_points_with_params(&mut fm, &distance(MergeMethod::WeightedAverage)).unwrap()
    );
    assert_xs(&fm, &[1.05, 3.0]);
    let mut fm = parse_points(&points(&[(1.0, 1.0), (1.1, 1.0), (3.0, 1.0), (4.0, 1.0)]));
    assert!(
        merge_feat_points_with_params(&mut fm, &distance(MergeMethod::WeightedAverage)).unwrap()
    );
    assert_xs(&fm, &[1.05, 3.5]);
    let mut fm = parse_points(&points(&[(1.0, 1.0), (1.2, 1.0), (1.3, 1.0), (4.0, 1.0)]));
    assert!(
        merge_feat_points_with_params(&mut fm, &distance(MergeMethod::WeightedAverage)).unwrap()
    );
    assert_xs(&fm, &[1.0, 1.25, 4.0]);
    let heavy = points(&[(1.0, 1.0), (1.2, 3.0), (1.3, 1.0), (4.0, 1.0)]);
    for (method, x) in [
        (MergeMethod::Average, 1.25),
        (MergeMethod::WeightedAverage, 1.225),
        (MergeMethod::UseLarger, 1.2),
    ] {
        let mut fm = parse_points(&heavy);
        assert!(merge_feat_points_with_params(&mut fm, &distance(method)).unwrap());
        assert_xs(&fm, &[1.0, x, 4.0]);
    }
    let mut fm = parse_points(&points(&[(1.0, 1.0), (1.3, 1.0), (1.2, 3.0), (4.0, 1.0)]));
    assert!(merge_feat_points_with_params(&mut fm, &distance(MergeMethod::UseLarger)).unwrap());
    assert_xs(&fm, &[1.0, 1.2, 4.0]);

    // test1BasicsRepeated.
    let mut fm = parse_points(&points(&[
        (0.7, 1.0),
        (1.0, 1.0),
        (1.2, 1.0),
        (1.3, 1.0),
        (4.0, 1.0),
    ]));
    let params = distance(MergeMethod::WeightedAverage);
    assert!(merge_feat_points_with_params(&mut fm, &params).unwrap());
    assert_xs(&fm, &[0.7, 1.0, 1.25, 4.0]);
    assert!(merge_feat_points_with_params(&mut fm, &params).unwrap());
    assert_xs(&fm, &[0.7, 1.125, 4.0]);
    assert!(merge_feat_points_with_params(&mut fm, &params).unwrap());
    assert_xs(&fm, &[0.9125, 4.0]);

    // test2ScoreBasics.
    let overlap = MergeFeatPointsParams {
        merge_metric: MergeMetric::Overlap,
        merge_method: MergeMethod::Average,
        ..Default::default()
    };
    for (pts, expected) in [
        (vec![(1.0, 1.0), (1.2, 3.0), (4.0, 1.0)], vec![1.1, 4.0]),
        (
            vec![(1.0, 1.0), (1.1, 1.0), (1.3, 3.0), (4.0, 1.0)],
            vec![1.15, 1.1, 4.0],
        ),
        (
            vec![(1.0, 1.0), (1.2, 1.0), (1.6, 3.0), (4.0, 1.0)],
            vec![1.0, 1.4, 4.0],
        ),
    ] {
        let mut fm = parse_points(&points(&pts));
        assert!(merge_feat_points_with_params(&mut fm, &overlap).unwrap());
        assert_xs(&fm, &expected);
    }
}

#[test]
fn utils_match_python() {
    // test_CombineFeatMaps.
    let fm1 = parse_points(&points(&[(1.0, 1.0)]));
    let fm2 = parse_points(&points(&[(0.0, 1.0)]));
    assert_eq!(combine_feat_maps(&fm1, &fm2).unwrap().get_num_features(), 2);
    for (tol, n) in [(1.0, 2), (1.1, 1)] {
        let params = CombineFeatMapsParams {
            merge_metric: MergeMetric::Distance,
            merge_tol: tol,
            ..Default::default()
        };
        assert_eq!(
            combine_feat_maps_with_params(&fm1, &fm2, &params)
                .unwrap()
                .get_num_features(),
            n
        );
    }

    let fm = parse_points(&points(&[(1.0, 1.0), (1.1, 2.0), (3.0, 1.0)]));
    let matrix = get_feat_feat_dist_matrix(
        &fm,
        MergeMetric::Distance,
        1.5,
        DirMergeMode::NoMerge,
        |a, b| families_match(a, b),
    )
    .unwrap();
    assert_eq!(matrix[0][0], 1e8);
    assert_eq!(
        (matrix[0][1], matrix[1][0]),
        (0.010000000000000018, 0.010000000000000018)
    );
    assert_eq!(matrix[0][2], 1e8);
    let matrix = get_feat_feat_dist_matrix(
        &fm,
        MergeMetric::Overlap,
        1.5,
        DirMergeMode::NoMerge,
        |_, _| true,
    )
    .unwrap();
    // Python's -0.0 for points beyond the radius, below the tolerance.
    assert_eq!(matrix[0][1], -2.0);
    assert!(matrix[0][2] == 0.0 && matrix[0][2].is_sign_negative());
    assert_eq!(matrix[2][2], 1e8);
    let none = get_feat_feat_dist_matrix(
        &fm,
        MergeMetric::NoMerge,
        1.5,
        DirMergeMode::NoMerge,
        |_, _| true,
    )
    .unwrap();
    assert!(none.iter().flatten().all(|&d| d == 1e8));
    assert!(feq(1.0, 1.00005) && !feq(1.0, 1.0002) && feq_with_tol(1.0, 1.1, 0.2));
}
