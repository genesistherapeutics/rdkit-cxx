// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.Subshape in RDKit 2024.03.6 (built with
// Eigen3, as conda-forge's is), and the tests of UnitTestSubshape.py. The
// alignments' counts are UnitTestSubshape.py's; the shapes' points, moments,
// and distances come from embedding and eigenvectors, so they are compared
// with tolerances.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_transforms::canonicalize_conformer;
use rdkit_cxx::chem::subshape::builder_utils::*;
use rdkit_cxx::chem::subshape::subshape_aligner::*;
use rdkit_cxx::chem::subshape::subshape_builder::*;
use rdkit_cxx::chem::subshape::subshape_objects::*;
use rdkit_cxx::chem::{AddHsParams, Mol, SDMolSupplier, add_hs_with_params};
use rdkit_cxx::geometry::rd_geometry::Point3D;

fn molecules() -> Vec<Mol> {
    let path = format!(
        "{}/tests/data/Subshape_5ht3ligs.sdf",
        env!("CARGO_MANIFEST_DIR")
    );
    let mut suppl = SDMolSupplier::new(path).unwrap();
    let n = suppl.len().unwrap();
    let params = AddHsParams {
        add_coords: true,
        ..Default::default()
    };
    (0..n)
        .map(|i| {
            let mut m = add_hs_with_params(&suppl.get(i).unwrap(), &params).unwrap();
            canonicalize_conformer(&mut m.get_conformer_mut().unwrap()).unwrap();
            m
        })
        .collect()
}

fn builder() -> SubshapeBuilder {
    SubshapeBuilder {
        grid_dims: (20.0, 20.0, 10.0),
        grid_spacing: 0.5,
        win_rad: 4.0,
        ..Default::default()
    }
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[test]
fn subshapes_match_python() {
    let builder = builder();
    let mut ms = molecules();
    assert_eq!(ms.len(), 4);
    let terminal = GenerateSubshapeShapeParams {
        terminal_pts_only: true,
        ..Default::default()
    };
    let shapes: Vec<ShapeWithSkeleton> = ms
        .iter()
        .map(|m| builder.call_with_params(m, &terminal).unwrap())
        .collect();
    assert_eq!(
        shapes.iter().map(|s| s.skel_pts.len()).collect::<Vec<_>>(),
        [5, 5, 5, 5]
    );
    let expected_locations = [
        (3.161468, 1.112563, -0.144838),
        (-3.827586, 0.590582, 0.16356),
        (-0.708434, 0.126246, -0.127189),
        (2.171233, -1.690847, 0.3599),
        (-3.621544, -0.601959, -0.836694),
    ];
    let expected_moments = [
        (4.577346, 3.195432, 0.66169),
        (3.403519, 2.45853, 1.412388),
        (5.22222, 3.139929, 1.281316),
        (4.631228, 3.601128, 0.749204),
        (3.597608, 3.185625, 1.424715),
    ];
    for ((pt, loc), mom) in shapes[0]
        .skel_pts
        .iter()
        .zip(expected_locations)
        .zip(expected_moments)
    {
        let l = pt.location;
        assert!(
            close(l.x, loc.0, 1e-4) && close(l.y, loc.1, 1e-4) && close(l.z, loc.2, 1e-4),
            "{l:?}"
        );
        let m = pt.shape_moments;
        assert!(
            close(m[0], mom.0, 1e-4) && close(m[1], mom.1, 1e-4) && close(m[2], mom.2, 1e-4),
            "{m:?}"
        );
        // The directions are unit vectors, orthogonal to each other.
        let d = pt.shape_dirs.unwrap();
        for i in 0..3 {
            assert!(close(d[i].length(), 1.0, 1e-9));
            assert!(close(d[i].dot_product(&d[(i + 1) % 3]), 0.0, 1e-9));
        }
    }

    let mut ref_shape = builder.generate_subshape_shape(&ms[0]).unwrap();
    assert_eq!(ref_shape.skel_pts.len(), 15);
    let frac_vols = [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.379236, 0.483526, 0.220859, 0.342556, 0.47762, 0.321729,
        0.562948, 0.59108, 0.337116, 0.444515,
    ];
    for (pt, f) in ref_shape.skel_pts.iter().zip(frac_vols) {
        assert!(close(pt.frac_vol, f, 1e-4), "{}", pt.frac_vol);
    }

    let aligner = SubshapeAligner {
        shape_dist_tol: 0.30,
        ..Default::default()
    };
    let target_mol = ms[0].clone();
    let mut alg_store = vec![Vec::new()];
    let expected_stats = [
        (1113, vec![("coarseGrid", 127.0), ("direction", 981.0)]),
        (
            1713,
            vec![
                ("coarseGrid", 210.0),
                ("direction", 1480.0),
                ("fineGrid", 1.0),
                ("medGrid", 1.0),
            ],
        ),
        (661, vec![("coarseGrid", 103.0), ("direction", 558.0)]),
    ];
    for i in 1..4 {
        let n_matches = aligner
            .get_triangle_matches(&ref_shape, &shapes[i])
            .unwrap()
            .len();
        let (algs, stats) = aligner
            .get_subshape_alignments_with_outputs(
                &target_mol,
                &mut ref_shape,
                &mut ms[i],
                &shapes[i],
                &builder,
                &Default::default(),
            )
            .unwrap();
        let counts: Vec<(&str, f64)> = stats
            .iter()
            .filter(|(k, _)| !k.ends_with("_time"))
            .map(|(k, v)| (k.as_str(), *v))
            .collect();
        assert_eq!(
            (n_matches, counts),
            (expected_stats[i - 1].0, expected_stats[i - 1].1.clone())
        );
        alg_store.push(algs);
    }
    assert_eq!(
        alg_store.iter().map(Vec::len).collect::<Vec<_>>(),
        [0, 5, 21, 0]
    );
    let expected = [
        ((1, 2, 3), (1, 2, 4), 0.274702, 2.827916),
        ((1, 2, 6), (1, 2, 3), 0.285975, 2.84445),
        ((1, 2, 7), (1, 2, 3), 0.21817, 2.847784),
        ((1, 2, 8), (1, 2, 4), 0.257129, 2.931189),
        ((1, 2, 14), (1, 2, 4), 0.260279, 2.858114),
        ((1, 6, 7), (1, 2, 3), 0.248674, 2.773738),
        ((1, 6, 8), (1, 2, 4), 0.258621, 2.916515),
        ((1, 6, 8), (1, 3, 4), 0.268402, 2.806569),
        ((1, 6, 11), (1, 2, 4), 0.290782, 2.899487),
        ((1, 6, 14), (1, 2, 4), 0.258455, 2.92637),
        ((1, 6, 14), (1, 3, 4), 0.24685, 2.729016),
        ((1, 7, 8), (1, 3, 4), 0.217672, 2.901536),
        ((1, 7, 11), (1, 3, 4), 0.234416, 2.791708),
        ((1, 7, 14), (1, 3, 4), 0.214523, 2.862917),
        ((1, 11, 14), (1, 2, 4), 0.285809, 2.955373),
        ((2, 6, 8), (2, 3, 4), 0.288793, 2.85599),
        ((2, 6, 14), (2, 3, 4), 0.291943, 2.815237),
        ((2, 7, 8), (2, 3, 4), 0.242706, 2.960987),
        ((2, 7, 11), (2, 3, 4), 0.242208, 2.840687),
        ((2, 7, 14), (2, 3, 4), 0.231598, 2.932908),
        ((9, 11, 14), (1, 2, 4), 0.299901, 2.788602),
    ];
    for (a, (t, q, shape_dist, dir_match)) in alg_store[2].iter().zip(expected) {
        assert_eq!((a.target_tri, a.query_tri), (t, q));
        assert!(
            close(a.shape_dist, shape_dist, 1e-4) && close(a.dir_match, dir_match, 1e-4),
            "{a:?}"
        );
    }

    // The aligner's call finds the same alignments.
    for i in 1..4 {
        let algs = aligner
            .call(
                &target_mol,
                &mut ref_shape,
                &mut ms[i],
                &shapes[i],
                &builder,
            )
            .unwrap();
        assert_eq!(algs, alg_store[i]);
    }

    let pruned: Vec<usize> = (0..4)
        .map(|i| {
            let params = ClusterAlignmentsParams {
                neighbor_tol: 0.15,
                ..Default::default()
            };
            cluster_alignments_with_params(&mut ms[i], &alg_store[i], &builder, &params)
                .unwrap()
                .len()
        })
        .collect();
    assert_eq!(pruned, [0, 4, 11, 0]);
}

#[test]
fn subshape_objects_and_helpers() {
    // UnitTestSubshape.TestSubshapeObjects.
    assert!(SubshapeShape::default().shapes.is_empty());
    let pt = SkeletonPoint::from_location(Point3D::new(1.0, 2.0, 3.0));
    assert_eq!(
        (pt.shape_moments, pt.shape_dirs, pt.frac_vol),
        ([0.0; 3], None, 0.0)
    );

    let builder = builder();
    let ms = molecules();
    let no_skeleton = GenerateSubshapeShapeParams {
        add_skeleton: false,
        ..Default::default()
    };
    let s0 = builder
        .generate_subshape_shape_with_params(&ms[0], &no_skeleton)
        .unwrap();
    let s1 = builder
        .generate_subshape_shape_with_params(&ms[1], &no_skeleton)
        .unwrap();
    assert!(s0.skel_pts.is_empty());
    let total = |s: &ShapeWithSkeleton| {
        s.grid
            .as_ref()
            .unwrap()
            .get_occupancy_vect()
            .get_total_val()
    };
    assert_eq!((total(&s0), total(&s1)), (6032, 7427));
    let union = builder.combine_subshapes(&s0, &s1).unwrap();
    let sum = builder
        .combine_subshapes_with_operation(&s0, &s1, SubshapeCombineOperations::Sum)
        .unwrap();
    let inter = builder
        .combine_subshapes_with_operation(&s0, &s1, SubshapeCombineOperations::Intersect)
        .unwrap();
    // The grids hold 2-bit values, so sums saturate at 3.
    assert_eq!(
        (total(&union), total(&sum), total(&inter)),
        (9265, 9618, 4194)
    );
    let tanimoto = get_shape_shape_distance(&s0, &s1, SubshapeDistanceMetric::Tanimoto).unwrap();
    let protrude = get_shape_shape_distance(&s0, &s1, SubshapeDistanceMetric::Protrude).unwrap();
    assert!(
        close(tanimoto, 0.5473286562331354, 1e-12) && close(protrude, 0.3047082228116711, 1e-12)
    );
    let coarse = builder.sample_subshape(&s0, 1.0).unwrap();
    assert_eq!(
        (coarse.grid.as_ref().unwrap().get_spacing(), total(&coarse)),
        (1.0, 747)
    );
    // A finer grid of the same size reaches beyond the shape's, which
    // Python cannot store (OverflowError).
    assert!(matches!(
        builder.sample_subshape(&s0, 0.25),
        Err(Error::Value { .. })
    ));

    let grid = s0.grid.as_ref().unwrap();
    assert_eq!(compute_grid_indices(grid, 1.0).len(), 93);
    let (count, centroid) =
        compute_shape_grid_centroid(&Point3D::new(0.0, 0.0, 0.0), grid, 2.0).unwrap();
    // Python caches the offsets on the grid whatever the radius, so its
    // value here is for a grid not yet used with another radius.
    assert_eq!(count, 999);
    let expected = (
        0.01001001001001001,
        -0.12312312312312312,
        -0.2122122122122122,
    );
    assert!(
        close(centroid.x, expected.0, 1e-12)
            && close(centroid.y, expected.1, 1e-12)
            && close(centroid.z, expected.2, 1e-12)
    );
    let far = Point3D::new(100.0, 0.0, 0.0);
    assert!(matches!(
        compute_shape_grid_centroid(&far, grid, 1.0),
        Err(Error::Value { .. })
    ));
    let pts = [
        SkeletonPoint::from_location(Point3D::new(0.0, 0.0, 0.0)),
        SkeletonPoint::from_location(Point3D::new(1.0, 0.0, 0.0)),
        SkeletonPoint::from_location(Point3D::new(5.0, 0.0, 0.0)),
    ];
    let clustered = cluster_terminal_pts(&pts, 2.0, 0.75);
    assert_eq!(
        clustered.iter().map(|p| p.location.x).collect::<Vec<_>>(),
        [0.5, 5.0]
    );

    let mut m = ms[0].clone();
    let shift = vec![
        vec![1.0, 0.0, 0.0, 2.0],
        vec![0.0, 1.0, 0.0, 0.0],
        vec![0.0, 0.0, 1.0, 0.0],
        vec![0.0, 0.0, 0.0, 1.0],
    ];
    transform_mol(&mut m, &shift).unwrap();
    let a = m
        .get_conformer_with_id(-1)
        .unwrap()
        .get_atom_position(0)
        .unwrap();
    let b = m
        .get_conformer_with_id(100)
        .unwrap()
        .get_atom_position(0)
        .unwrap();
    assert!(close(b.x - a.x, 2.0, 1e-12) && b.y == a.y && b.z == a.z);
    // Only a 4x4 transform will do.
    for bad in [
        shift[..3].to_vec(),
        [shift.clone(), shift[..1].to_vec()].concat(),
    ] {
        assert!(matches!(
            transform_mol(&mut m, &bad),
            Err(rdkit_cxx::Error::Value { .. })
        ));
    }
    let mut ragged = shift.clone();
    ragged[3].push(0.0);
    assert!(matches!(
        transform_mol(&mut m, &ragged),
        Err(rdkit_cxx::Error::Value { .. })
    ));
}
