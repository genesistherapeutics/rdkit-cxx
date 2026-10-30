// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::f64::consts::{FRAC_PI_2, PI};

use rdkit_cxx::Error;
use rdkit_cxx::data_structs::DiscreteValueType;
use rdkit_cxx::geometry::*;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn point3d() {
    let p = Point3D::new(1.0, 2.0, 3.0);
    let q = Point3D::new(-1.0, 0.0, 2.0);
    assert_eq!(p.dot_product(&q), 5.0);
    assert_eq!(p.cross_product(&q), Point3D::new(4.0, -5.0, 2.0));
    assert!(close(p.angle_to(&q), 0.9302740141154721));
    assert!(close(p.signed_angle_to(&q), 0.9302740141154721));
    assert_eq!(p.distance(&q), 3.0);
    let dir = p.direction_vector(&q).unwrap();
    assert!(close(dir.x, -2.0 / 3.0) && close(dir.y, -2.0 / 3.0) && close(dir.z, -1.0 / 3.0));
    assert_eq!(p.length_sq(), 14.0);
    assert!(close(p.length(), 14f64.sqrt()));

    assert_eq!(p + q, Point3D::new(0.0, 2.0, 5.0));
    assert_eq!(p - q, Point3D::new(2.0, 2.0, 1.0));
    assert_eq!(p * 2.0, Point3D::new(2.0, 4.0, 6.0));
    assert_eq!(p / 2.0, Point3D::new(0.5, 1.0, 1.5));
    let mut r = p;
    r += q;
    r -= q;
    r *= 3.0;
    r /= 3.0;
    assert_eq!(r, p);
    assert_eq!(Point3D::default(), Point3D::new(0.0, 0.0, 0.0));

    let mut n = Point3D::new(3.0, 0.0, 4.0);
    n.normalize().unwrap();
    assert_eq!(n, Point3D::new(0.6, 0.0, 0.8));
    let mut zero = Point3D::default();
    assert!(matches!(zero.normalize(), Err(Error::Runtime { .. })));
    // Python's AngleTo is NaN for a zero vector, too.
    assert!(zero.angle_to(&p).is_nan());
    assert!(matches!(p.direction_vector(&p), Err(Error::Runtime { .. })));
}

#[test]
fn point2d() {
    let a = Point2D::new(1.0, 0.0);
    let b = Point2D::new(0.0, 1.0);
    assert!(close(a.angle_to(&b).unwrap(), FRAC_PI_2));
    assert!(close(a.signed_angle_to(&b).unwrap(), FRAC_PI_2));
    assert!(close(b.signed_angle_to(&a).unwrap(), 3.0 * FRAC_PI_2));
    assert_eq!(a.dot_product(&b), 0.0);
    assert_eq!((a + b) * 2.0, Point2D::new(2.0, 2.0));
    assert_eq!((a - b) / 2.0, Point2D::new(0.5, -0.5));
    assert!(close((a + b).length(), 2f64.sqrt()));
    assert_eq!((a + b).length_sq(), 2.0);
    let dir = a.direction_vector(&b).unwrap();
    assert!(close(dir.x, -dir.y) && close(dir.length(), 1.0));
    assert_eq!(
        Point2D::from(Point3D::new(1.0, 2.0, 3.0)),
        Point2D::new(1.0, 2.0)
    );
    // Point2D's AngleTo normalizes, so a zero vector raises.
    assert!(matches!(
        Point2D::default().angle_to(&a),
        Err(Error::Runtime { .. })
    ));
    let mut zero = Point2D::default();
    assert!(zero.normalize().is_err());
}

#[test]
fn dihedral_angles() {
    let (p1, p2, p3) = (
        Point3D::new(1.0, 0.0, 0.0),
        Point3D::default(),
        Point3D::new(0.0, 1.0, 0.0),
    );
    let up = Point3D::new(0.0, 1.0, 1.0);
    let down = Point3D::new(0.0, 1.0, -1.0);
    assert!(close(compute_dihedral_angle(&p1, &p2, &p3, &up), FRAC_PI_2));
    assert!(close(
        compute_signed_dihedral_angle(&p1, &p2, &p3, &up),
        -FRAC_PI_2
    ));
    assert!(close(
        compute_dihedral_angle(&p1, &p2, &p3, &down),
        FRAC_PI_2
    ));
    assert!(close(
        compute_signed_dihedral_angle(&p1, &p2, &p3, &down),
        FRAC_PI_2
    ));
}

#[test]
fn point_nd() {
    let mut n = PointND::new(4);
    let mut m = PointND::new(4);
    for i in 0..4 {
        n.set(i, f64::from(i) + 1.0).unwrap();
        m.set(i, 4.0 - f64::from(i)).unwrap();
    }
    assert_eq!(n.len(), 4);
    assert_eq!(n.coords(), [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(n.get(2).unwrap(), 3.0);
    assert_eq!(n.get(4).unwrap_err(), Error::Index { index: 4, len: 4 });
    assert_eq!(n.dot_product(&m).unwrap(), 20.0);
    assert!(close(n.length(), 30f64.sqrt()));
    assert_eq!(n.length_sq(), 30.0);
    assert!(close(n.angle_to(&m).unwrap(), 0.8410686705679303));
    let dir = n.direction_vector(&m).unwrap();
    assert!(close(dir.get(0).unwrap(), 0.6708203932499369));
    assert!(close(dir.get(3).unwrap(), -0.6708203932499369));
    assert!(matches!(
        n.dot_product(&PointND::new(3)),
        Err(Error::Value { .. })
    ));

    assert_eq!((&n + &m).coords(), [5.0; 4]);
    assert_eq!((&n - &m).coords(), [-3.0, -1.0, 1.0, 3.0]);
    assert_eq!((&n * 2.0).coords(), [2.0, 4.0, 6.0, 8.0]);
    assert_eq!((&n / 2.0).coords(), [0.5, 1.0, 1.5, 2.0]);
    let mut k = n.clone();
    k += &m;
    k -= &m;
    k *= 2.0;
    k /= 2.0;
    assert_eq!(k, n);
    k.normalize().unwrap();
    assert!(close(k.length(), 1.0));
    assert!(PointND::new(0).is_empty());
    assert_eq!(PointND::from_coords(vec![1.0, 2.0]).coords(), [1.0, 2.0]);
}

#[test]
#[should_panic(expected = "dimension mismatch")]
fn point_nd_operator_dimensions() {
    let _ = &PointND::new(2) + &PointND::new(3);
}

fn sphere_grid(center: Point3D) -> UniformGrid3D {
    let mut grid = UniformGrid3D::new(10.0, 10.0, 10.0).unwrap();
    grid.set_sphere_occupancy(&center, 1.5, 0.25).unwrap();
    grid
}

#[test]
fn uniform_grid3d() {
    let grid = UniformGrid3D::new(10.0, 10.0, 10.0).unwrap();
    assert_eq!(
        (grid.get_num_x(), grid.get_num_y(), grid.get_num_z()),
        (20, 20, 20)
    );
    assert_eq!(grid.get_size(), 8000);
    assert_eq!(grid.get_spacing(), 0.5);
    assert_eq!(grid.get_offset(), Point3D::new(-5.0, -5.0, -5.0));
    assert_eq!(grid.get_grid_point_index(&Point3D::default()), 4210);
    assert_eq!(grid.get_grid_index(1, 2, 3), 1241);
    assert_eq!(grid.get_grid_indices(321).unwrap(), (1, 16, 0));
    assert_eq!(
        grid.get_grid_point_loc(321).unwrap(),
        Point3D::new(-4.5, 3.0, -5.0)
    );
    assert_eq!(
        grid.get_grid_indices(100_000).unwrap_err(),
        Error::Index {
            index: 100_000,
            len: 8000
        }
    );
    assert!(grid.get_grid_point_loc(8000).is_err());
    assert!(grid.get_val(8000).is_err());

    // Points outside the grid, however far, have no index.
    let far = Point3D::new(100.0, 0.0, 0.0);
    assert_eq!(grid.get_grid_point_index(&far), -1);
    assert_eq!(grid.get_val_point(&far), -1);
    assert_eq!(grid.get_grid_index(100, 0, 0), -1);
    for far in [
        Point3D::new(1e300, 0.0, 0.0),
        Point3D::new(f64::NAN, 0.0, 0.0),
        Point3D::new(0.0, f64::NEG_INFINITY, 0.0),
    ] {
        assert_eq!(grid.get_grid_point_index(&far), -1);
        assert_eq!(grid.get_val_point(&far), -1);
    }

    let mut grid = sphere_grid(Point3D::default());
    assert_eq!(grid.get_occupancy_vect().get_total_val(), 523);
    assert_eq!(grid.get_val_point(&Point3D::default()), 3);
    assert_eq!(grid.get_val(4210).unwrap(), 3);
    grid.set_val(0, 2).unwrap();
    assert_eq!(grid.get_val(0).unwrap(), 2);
    grid.set_val_point(&Point3D::new(-5.0, -5.0, -5.0), 1)
        .unwrap();
    assert_eq!(grid.get_val(0).unwrap(), 1);
    // A point outside the grid sets nothing.
    grid.set_val_point(&Point3D::new(1e300, 0.0, 0.0), 1)
        .unwrap();
    assert!(matches!(grid.set_val(0, 4), Err(Error::Value { .. })));

    let copy = UniformGrid3D::from_pkl(&grid.to_binary()).unwrap();
    assert!(copy.compare_params(&grid));
    assert_eq!(copy.get_occupancy_vect(), grid.get_occupancy_vect());
    assert_eq!(grid.clone().to_binary(), grid.to_binary());
}

#[test]
fn uniform_grid3d_params() {
    let grid = UniformGrid3D::new_with_params(
        6.0,
        4.0,
        2.0,
        &UniformGrid3DParams {
            spacing: 1.0,
            val_type: DiscreteValueType::Fourbitvalue,
            off_set: Some(Point3D::new(1.0, 1.0, 1.0)),
        },
    )
    .unwrap();
    assert_eq!(
        (grid.get_num_x(), grid.get_num_y(), grid.get_num_z()),
        (6, 4, 2)
    );
    assert_eq!(grid.get_offset(), Point3D::new(1.0, 1.0, 1.0));
    assert_eq!(
        grid.get_occupancy_vect().get_value_type(),
        DiscreteValueType::Fourbitvalue
    );
    assert!(!grid.compare_params(&UniformGrid3D::new(6.0, 4.0, 2.0).unwrap()));

    assert!(matches!(
        UniformGrid3D::new(-1.0, 1.0, 1.0),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        UniformGrid3D::new(f64::NAN, 1.0, 1.0),
        Err(Error::Runtime { .. })
    ));
    // RDKit casts these sizes to unsigned int unchecked.
    for (dim, spacing) in [(f64::INFINITY, 0.5), (1e12, 0.5), (1e3, 1e-6)] {
        let params = UniformGrid3DParams {
            spacing,
            ..Default::default()
        };
        assert!(matches!(
            UniformGrid3D::new_with_params(dim, 1.0, 1.0, &params),
            Err(Error::Value { .. })
        ));
    }
    assert!(matches!(
        UniformRealValueGrid3D::new(f64::INFINITY, 1.0, 1.0),
        Err(Error::Value { .. })
    ));
}

#[test]
fn sphere_occupancy() {
    let mut grid = UniformGrid3D::new(10.0, 10.0, 10.0).unwrap();
    let outside = Point3D::new(100.0, 0.0, 0.0);
    grid.set_sphere_occupancy(&outside, 1.5, 0.25).unwrap();
    assert_eq!(grid.get_occupancy_vect().get_total_val(), 0);
    let strict = SetSphereOccupancyParams {
        ignore_out_of_bound: false,
        ..Default::default()
    };
    assert!(matches!(
        grid.set_sphere_occupancy_with_params(&outside, 1.5, 0.25, &strict),
        Err(Error::Runtime { .. })
    ));
    let center = Point3D::default();
    for (radius, step) in [(1.5, -0.25), (f64::NAN, 0.25), (1e300, 0.25), (1.5, 1e300)] {
        assert!(matches!(
            grid.set_sphere_occupancy(&center, radius, step),
            Err(Error::Value { .. })
        ));
    }
    let one_layer = SetSphereOccupancyParams {
        max_layers: 1,
        ..Default::default()
    };
    grid.set_sphere_occupancy_with_params(&center, 1.0, 0.5, &one_layer)
        .unwrap();
    assert_eq!(grid.get_val_point(&center), 3);
}

#[test]
fn grid_similarity() {
    let a = sphere_grid(Point3D::default());
    let b = sphere_grid(Point3D::new(1.0, 0.0, 0.0));
    assert!(close(tanimoto_distance(&a, &b).unwrap(), 0.572987721691678));
    assert!(close(
        protrude_distance(&a, &b).unwrap(),
        0.40152963671128106
    ));
    assert!(close(
        tversky_index(&a, &b, 0.5, 0.5).unwrap(),
        0.5984703632887189
    ));
    let small = UniformGrid3D::new(5.0, 5.0, 5.0).unwrap();
    assert!(matches!(
        tanimoto_distance(&a, &small),
        Err(Error::Value { .. })
    ));

    let (weight_sum, centroid) = compute_grid_centroid(&a, &Point3D::default(), 1.0).unwrap();
    assert_eq!(weight_sum, 279.0);
    assert!(centroid.length() < 1e-9);
    // RDKit uses only the point's grid index, -1 for any point outside.
    assert_eq!(
        compute_grid_centroid(&a, &Point3D::new(100.0, 0.0, 0.0), 1.0)
            .unwrap()
            .0,
        compute_grid_centroid(&a, &Point3D::new(f64::NAN, 0.0, 0.0), 1.0)
            .unwrap()
            .0
    );
    assert!(matches!(
        compute_grid_centroid(&a, &Point3D::default(), 1e9),
        Err(Error::Value { .. })
    ));

    let terminal = find_grid_terminal_points(&a, 1.0, 0.5).unwrap();
    assert_eq!(terminal.len(), 12);
    assert!(close(terminal[0].y, -0.656934306569343));
    assert!(matches!(
        find_grid_terminal_points(&a, f64::INFINITY, 0.5),
        Err(Error::Value { .. })
    ));

    let mut sum = UniformGrid3D::new(10.0, 10.0, 10.0).unwrap();
    sum += &a;
    sum += &b;
    assert_eq!(sum.get_occupancy_vect().get_total_val(), 781);
    sum -= &b;
    assert_eq!(sum.get_occupancy_vect().get_total_val(), 258);
    let mut both = a.clone();
    both &= &b;
    let mut either = a.clone();
    either |= &b;
    assert!(
        both.get_occupancy_vect().get_total_val() < either.get_occupancy_vect().get_total_val()
    );
}

#[test]
#[should_panic(expected = "parameters differ")]
fn grid_operator_params() {
    let mut a = UniformGrid3D::new(10.0, 10.0, 10.0).unwrap();
    a += &UniformGrid3D::new(5.0, 5.0, 5.0).unwrap();
}

#[test]
fn write_grid() {
    let grid = sphere_grid(Point3D::default());
    let dir = std::env::temp_dir().join(format!("rdkit-cxx-grid-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("grid.grd");
    write_grid_to_file(&grid, &path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.starts_with("Grid file representing a Shape"));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn uniform_real_value_grid3d() {
    let default = UniformRealValueGrid3D::default();
    assert_eq!(
        (
            default.get_num_x(),
            default.get_size(),
            default.get_spacing()
        ),
        (1, 1, 1.0)
    );
    assert_eq!(default.get_offset(), Point3D::default());

    let mut grid = UniformRealValueGrid3D::new(4.0, 4.0, 4.0).unwrap();
    assert_eq!(grid.get_num_x(), 8);
    assert_eq!(grid.get_size(), 512);
    assert_eq!(grid.get_offset(), Point3D::new(-2.0, -2.0, -2.0));
    grid.set_val(3, 2.5).unwrap();
    grid.set_val_point(&Point3D::default(), 1.25).unwrap();
    assert_eq!(grid.get_val(3).unwrap(), 2.5);
    assert_eq!(grid.get_val_point(&Point3D::default()), 1.25);
    assert_eq!(grid.get_val_point(&Point3D::new(10.0, 0.0, 0.0)), -1.0);
    assert_eq!(grid.get_val_point(&Point3D::new(1e300, 0.0, 0.0)), -1.0);
    assert_eq!(grid.get_occupancy_vect().get_total_val(), 3.75);
    assert_eq!(
        grid.get_val(1000).unwrap_err(),
        Error::Index {
            index: 1000,
            len: 512
        }
    );

    let mut copy = grid.clone();
    assert!(copy.compare_grids(&grid) && copy.compare_params(&grid) && copy.compare_vectors(&grid));
    copy.set_val(0, 7.0).unwrap();
    assert!(!copy.compare_grids(&grid));
    assert!(copy.compare_params(&grid));
    assert!(!copy.compare_vectors(&grid));
    copy += &grid;
    assert_eq!(copy.get_occupancy_vect().get_total_val(), 14.5);
    copy &= &grid;
    assert_eq!(copy.get_occupancy_vect().get_total_val(), 3.75);
    copy |= &grid;
    assert_eq!(copy.get_occupancy_vect().get_total_val(), 3.75);
    copy -= &grid;
    assert_eq!(copy.get_occupancy_vect().get_total_val(), 0.0);

    let pickled = UniformRealValueGrid3D::from_pkl(&grid.to_binary()).unwrap();
    assert!(pickled.compare_grids(&grid));
    let params = UniformRealValueGrid3DParams {
        spacing: 1.0,
        off_set: Some(Point3D::default()),
    };
    let other = UniformRealValueGrid3D::new_with_params(2.0, 2.0, 2.0, &params).unwrap();
    assert_eq!(other.get_size(), 8);
    assert!(!other.compare_params(&grid));
}

#[test]
fn point3d_signed_angle() {
    // Point3D's SignedAngleTo uses only the cross product's z.
    let p = Point3D::new(0.0, 1.0, 0.0);
    let q = Point3D::new(1.0, 0.0, 0.0);
    assert!(close(p.signed_angle_to(&q), 2.0 * PI - FRAC_PI_2));
}
