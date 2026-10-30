// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Expected values are from Python's rdDepictor. Coordinates depend on the
// platform's floating point, so the tests compare invariants (bond lengths,
// matched atoms at the reference's positions) with tolerances.

mod common;

use std::collections::BTreeMap;
use std::sync::{Mutex, MutexGuard};

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_depictor::*;
use rdkit_cxx::chem::{
    Atom, ChiralType, Conformer, GetPropNamesParams, Mol, RWMol, SmilesParserParams,
    mol_from_smarts, mol_from_smiles, mol_from_smiles_with_params,
};
use rdkit_cxx::geometry::{Point2D, Point3D};

// The CoordGen preference and the ring system templates are process-wide;
// tests that change them, or whose results depend on them, take this.
static SETTINGS: Mutex<()> = Mutex::new(());

fn settings() -> MutexGuard<'static, ()> {
    SETTINGS.lock().unwrap_or_else(|e| e.into_inner())
}

fn smiles(text: &str) -> Mol {
    mol_from_smiles(text).unwrap()
}

fn unsanitized(text: &str) -> Mol {
    let params = SmilesParserParams {
        sanitize: false,
        remove_hs: false,
        ..SmilesParserParams::default()
    };
    mol_from_smiles_with_params(text, &params).unwrap()
}

fn position(mol: &Mol, conf_id: i32, idx: u32) -> Point3D {
    mol.get_conformer_with_id(conf_id)
        .unwrap()
        .get_atom_position(idx)
        .unwrap()
}

fn distance(a: Point3D, b: Point3D) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

fn bond_lengths(mol: &Mol, conf_id: i32) -> Vec<f64> {
    mol.get_bonds()
        .map(|b| {
            distance(
                position(mol, conf_id, b.get_begin_atom_idx()),
                position(mol, conf_id, b.get_end_atom_idx()),
            )
        })
        .collect()
}

fn assert_bond_lengths(mol: &Mol, conf_id: i32, expected: f64, tol: f64) {
    for length in bond_lengths(mol, conf_id) {
        assert!(
            (length - expected).abs() < tol,
            "bond length {length}, expected {expected}"
        );
    }
}

fn set_position(mol: &mut Mol, idx: u32, p: Point3D) {
    mol.get_conformer_mut()
        .unwrap()
        .set_atom_position(idx, p)
        .unwrap();
}

fn is_value_error(result: Result<impl std::fmt::Debug, Error>, needle: &str) -> bool {
    matches!(result, Err(Error::Value { message }) if message.contains(needle))
}

#[test]
fn defaults_match_python() {
    let p = ConstrainedDepictionParams::default();
    assert!(!p.accept_failure && !p.force_rdkit && !p.allow_r_groups && !p.align_only);
    assert!(p.adjust_mol_block_wedging);
    assert_eq!((p.existing_conf_id, p.conf_id), (-1, -1));
    assert!(p.ref_patt.is_none());

    let c = Compute2DCoordsParams::default();
    assert!(c.canon_orient && c.clear_confs && c.coord_map.is_empty());
    assert_eq!((c.n_flips_per_sample, c.n_sample, c.sample_seed), (0, 0, 0));
    assert_eq!((c.permute_deg4_nodes, c.bond_length), (0, -1.0));
    assert!(!c.force_rdkit && !c.use_ring_templates);

    let m = Compute2DCoordsMimicDistmatParams::default();
    assert!(!m.canon_orient && m.clear_confs && m.permute_deg4_nodes && !m.force_rdkit);
    assert_eq!((m.weight_dist_mat, m.bond_length), (0.5, -1.0));
    assert_eq!(
        (m.n_flips_per_sample, m.n_sample, m.sample_seed),
        (3, 100, 100)
    );

    let g = GenerateDepictionMatching3DStructureParams::default();
    assert_eq!(g.conf_id, -1);
    assert!(g.ref_patt.is_none() && !g.accept_failure && !g.force_rdkit);

    assert_eq!(
        StraightenDepictionParams::default(),
        StraightenDepictionParams {
            conf_id: -1,
            minimize_rotation: false
        }
    );
    assert_eq!(
        NormalizeDepictionParams::default(),
        NormalizeDepictionParams {
            conf_id: -1,
            canonicalize: 1,
            scale_factor: -1.0
        }
    );
}

#[test]
fn compute_2d_coords_with_rdkit() {
    let _settings = settings();
    let mut mol = smiles("c1ccccc1CCO");
    assert_eq!(compute_2d_coords(&mut mol).unwrap(), 0);
    let conf = mol.get_conformer().unwrap();
    assert!(!conf.is_3d());
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);

    // bondLength holds for one call.
    let longer = Compute2DCoordsParams {
        bond_length: 2.0,
        ..Default::default()
    };
    assert_eq!(compute_2d_coords_with_params(&mut mol, &longer).unwrap(), 0);
    assert_bond_lengths(&mol, -1, 2.0, 1e-3);
    compute_2d_coords(&mut mol).unwrap();
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);

    let keep = Compute2DCoordsParams {
        clear_confs: false,
        ..Default::default()
    };
    assert_eq!(compute_2d_coords_with_params(&mut mol, &keep).unwrap(), 1);
    assert_eq!(mol.get_num_conformers(), 2);

    // Random sampling.
    let mut branched = smiles("CC(C)CC(C)CC(C)CO");
    let sampled = Compute2DCoordsParams {
        n_flips_per_sample: 2,
        n_sample: 20,
        sample_seed: 42,
        ..Default::default()
    };
    assert_eq!(
        compute_2d_coords_with_params(&mut branched, &sampled).unwrap(),
        0
    );
    assert_bond_lengths(&branched, -1, 1.5, 1e-3);

    let mut empty = Mol::default();
    assert_eq!(compute_2d_coords(&mut empty).unwrap(), 0);
    assert_eq!(empty.get_num_conformers(), 1);
}

#[test]
fn coordinate_maps_fix_atoms() {
    let _settings = settings();
    let mut mol = smiles("c1ccccc1CCO");
    let params = Compute2DCoordsParams {
        coord_map: BTreeMap::from([(0, Point2D::new(10.0, 10.0)), (1, Point2D::new(11.5, 10.0))]),
        ..Default::default()
    };
    compute_2d_coords_with_params(&mut mol, &params).unwrap();
    assert!(distance(position(&mol, -1, 0), Point3D::new(10.0, 10.0, 0.0)) < 1e-6);
    assert!(distance(position(&mol, -1, 1), Point3D::new(11.5, 10.0, 0.0)) < 1e-6);
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);

    // A single atom translates the depiction.
    let single = Compute2DCoordsParams {
        coord_map: BTreeMap::from([(7, Point2D::new(5.0, 5.0))]),
        ..Default::default()
    };
    compute_2d_coords_with_params(&mut mol, &single).unwrap();
    assert!(distance(position(&mol, -1, 7), Point3D::new(5.0, 5.0, 0.0)) < 1e-6);
}

#[test]
fn compute_2d_coords_refuses_bad_input() {
    let mut mol = smiles("CCO");
    // Python checks the index (as a ValueError); CoordGen would ignore it.
    let past_end = Compute2DCoordsParams {
        coord_map: BTreeMap::from([(3, Point2D::new(0.0, 0.0))]),
        ..Default::default()
    };
    assert_eq!(
        compute_2d_coords_with_params(&mut mol, &past_end).unwrap_err(),
        Error::Index { index: 3, len: 3 }
    );
    for p in [
        Point2D::new(f64::NAN, 0.0),
        Point2D::new(0.0, f64::INFINITY),
        Point2D::new(2e5, 0.0),
    ] {
        let bad = Compute2DCoordsParams {
            coord_map: BTreeMap::from([(0, p), (1, Point2D::new(1.5, 0.0))]),
            ..Default::default()
        };
        assert!(is_value_error(
            compute_2d_coords_with_params(&mut mol, &bad),
            "coordMap coordinates"
        ));
    }
    let infinite = Compute2DCoordsParams {
        bond_length: f64::INFINITY,
        ..Default::default()
    };
    assert!(is_value_error(
        compute_2d_coords_with_params(&mut mol, &infinite),
        "bondLength"
    ));
    assert_eq!(mol.get_num_conformers(), 0);
    // A bond length that is not positive, NaN included, keeps RDKit's.
    let _settings = settings();
    let nan = Compute2DCoordsParams {
        bond_length: f64::NAN,
        ..Default::default()
    };
    compute_2d_coords_with_params(&mut mol, &nan).unwrap();
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);
}

#[test]
fn rdkit_depictor_refuses_what_it_would_mishandle() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    // embedSquarePlanar reads the first neighbor of each square-planar atom
    // unchecked; Python segfaults here.
    let mut lone = smiles("[Pt]");
    lone.get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiSquareplanar);
    assert!(is_value_error(
        compute_2d_coords(&mut lone),
        "square-planar"
    ));
    let mut planar = unsanitized("[Pt](F)(F)(F)F");
    planar
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_chiral_tag(ChiralType::ChiSquareplanar);
    compute_2d_coords(&mut planar).unwrap();

    // embedTBP indexes a table of five positions by neighbor; Python reads
    // past it for six neighbors, or for four without a permutation.
    for text in ["[P](F)(F)(F)(F)(F)F", "[P](F)(F)(F)F"] {
        let mut tbp = unsanitized(text);
        tbp.get_atom_with_idx_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralType::ChiTrigonalbipyramidal);
        assert!(is_value_error(
            compute_2d_coords(&mut tbp),
            "trigonal-bipyramidal"
        ));
    }
    let mut tbp = smiles("[P@TB1](F)(F)(F)(F)F");
    compute_2d_coords(&mut tbp).unwrap();
    assert_bond_lengths(&tbp, -1, 1.5, 1e-3);

    // getDistanceMat computes n * n in int.
    let mut huge = unsanitized(&vec!["C"; 46341].join("."));
    assert!(is_value_error(compute_2d_coords(&mut huge), "46340"));
    // computeInitialCoords ranks atoms as rank * n + index in int, and a
    // hydrogen's rank is 100001.
    let mut hydrogens = unsanitized(&vec!["[H]"; 21475].join("."));
    assert!(is_value_error(
        compute_2d_coords(&mut hydrogens),
        "ranks atoms"
    ));
}

#[test]
fn mimic_distance_matrix() {
    let mut mol = smiles("CCCCC");
    // The topological distances times 1.5, lower triangle by rows.
    let n = 5;
    let mut dist = Vec::new();
    for i in 0..n {
        for j in 0..i {
            dist.push(1.5 * (i - j) as f64);
        }
    }
    assert_eq!(compute_2d_coords_mimic_distmat(&mut mol, &dist).unwrap(), 0);
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);
    assert!(is_value_error(
        compute_2d_coords_mimic_distmat(&mut mol, &dist[1..]),
        "array size"
    ));
    let mut empty = Mol::default();
    assert_eq!(compute_2d_coords_mimic_distmat(&mut empty, &[]).unwrap(), 0);

    // RDKit caches the distance matrix it uses on the molecule (U6); none is
    // left behind to go stale when the molecule grows.
    let names = mol
        .get_prop_names_with_params(&GetPropNamesParams {
            include_private: true,
            include_computed: true,
        })
        .unwrap();
    assert!(!names.iter().any(|name| name.contains("DistanceMatrix")));
    let mut grown = RWMol::from(mol);
    let idx = grown.add_atom(&Atom::new("O").unwrap());
    grown.add_bond(4, idx).unwrap();
    dist.extend([7.5, 6.0, 4.5, 3.0, 1.5]);
    let mut grown = Mol::from(grown);
    compute_2d_coords_mimic_distmat(&mut grown, &dist).unwrap();
    assert_bond_lengths(&grown, -1, 1.5, 1e-3);
}

fn benzene_reference() -> Mol {
    let mut reference = smiles("c1ccccc1");
    compute_2d_coords(&mut reference).unwrap();
    reference
}

#[test]
fn matching_2d_structure() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    let reference = benzene_reference();
    let mut mol = smiles("Oc1ccccc1CC");
    let matched = generate_depiction_matching_2d_structure(&mut mol, &reference).unwrap();
    assert_eq!(
        matched,
        vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]
    );
    assert_eq!(mol.get_num_conformers(), 1);
    for &(r, m) in &matched {
        assert!(distance(position(&mol, -1, m), position(&reference, -1, r)) < 1e-4);
    }

    let patterned = ConstrainedDepictionParams {
        ref_patt: Some(mol_from_smarts("c1ccccc1").unwrap()),
        ..Default::default()
    };
    let mut mol = smiles("Oc1ccccc1CC");
    assert_eq!(
        generate_depiction_matching_2d_structure_with_params(&mut mol, &reference, &patterned)
            .unwrap(),
        matched
    );

    // alignOnly keeps the molecule's own depiction and aligns it.
    let mut mol = smiles("Oc1ccccc1CC");
    compute_2d_coords(&mut mol).unwrap();
    let align = ConstrainedDepictionParams {
        align_only: true,
        ..Default::default()
    };
    let aligned =
        generate_depiction_matching_2d_structure_with_params(&mut mol, &reference, &align).unwrap();
    assert_eq!(aligned.len(), 6);
    for &(r, m) in &aligned {
        assert!(distance(position(&mol, -1, m), position(&reference, -1, r)) < 1e-3);
    }

    // Python's DepictException is a ValueError.
    let mut other = smiles("CCO");
    assert!(is_value_error(
        generate_depiction_matching_2d_structure(&mut other, &reference),
        "Substructure match with reference not found"
    ));
    let accept = ConstrainedDepictionParams {
        accept_failure: true,
        ..Default::default()
    };
    assert_eq!(
        generate_depiction_matching_2d_structure_with_params(&mut other, &reference, &accept)
            .unwrap(),
        Vec::<(u32, u32)>::new()
    );
    assert_eq!(other.get_num_conformers(), 1);
    let missing = ConstrainedDepictionParams {
        conf_id: 5,
        ..Default::default()
    };
    assert!(is_value_error(
        generate_depiction_matching_2d_structure_with_params(&mut mol, &reference, &missing),
        "conformation"
    ));
}

#[test]
fn matching_2d_structure_atom_map() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    let reference = benzene_reference();
    let mut mol = smiles("Oc1ccccc1CC");
    generate_depiction_matching_2d_structure_atom_map(&mut mol, &reference, &[(0, 1), (1, 2)])
        .unwrap();
    assert!(distance(position(&mol, -1, 1), position(&reference, -1, 0)) < 1e-4);
    assert!(distance(position(&mol, -1, 2), position(&reference, -1, 1)) < 1e-4);

    let forced = ConstrainedDepictionParams {
        force_rdkit: true,
        ..Default::default()
    };
    generate_depiction_matching_2d_structure_atom_map_with_params(
        &mut mol,
        &reference,
        &[(0, 1), (1, 2)],
        &forced,
    )
    .unwrap();

    // RDKit's own check lets an index equal to the atom count through, and
    // then sets a bit past the end of a dynamic_bitset.
    assert_eq!(
        generate_depiction_matching_2d_structure_atom_map(&mut mol, &reference, &[(6, 0)])
            .unwrap_err(),
        Error::Index { index: 6, len: 6 }
    );
    assert_eq!(
        generate_depiction_matching_2d_structure_atom_map(&mut mol, &reference, &[(0, 9)])
            .unwrap_err(),
        Error::Index { index: 9, len: 9 }
    );
    let too_long: Vec<(u32, u32)> = (0..7).map(|i| (i % 6, i)).collect();
    assert!(is_value_error(
        generate_depiction_matching_2d_structure_atom_map(&mut mol, &reference, &too_long),
        "refMatchVect"
    ));
}

#[test]
fn matching_refuses_bad_coordinates() {
    let _settings = settings();
    let mut reference = benzene_reference();
    set_position(&mut reference, 0, Point3D::new(f64::NAN, 0.0, 0.0));
    let mut mol = smiles("Oc1ccccc1CC");
    assert!(is_value_error(
        generate_depiction_matching_2d_structure(&mut mol, &reference),
        "the reference's coordinates"
    ));
    assert!(is_value_error(
        generate_depiction_matching_2d_structure_atom_map(&mut mol, &reference, &[(1, 1)]),
        "the reference's coordinates"
    ));
    let reference = benzene_reference();
    compute_2d_coords(&mut mol).unwrap();
    set_position(&mut mol, 3, Point3D::new(1e6, 0.0, 0.0));
    let align = ConstrainedDepictionParams {
        align_only: true,
        ..Default::default()
    };
    assert!(is_value_error(
        generate_depiction_matching_2d_structure_with_params(&mut mol, &reference, &align),
        "the molecule's coordinates"
    ));
}

fn propanol_3d() -> Mol {
    // Python's EmbedMolecule(AddHs(CCCO), randomSeed=42), without the Hs.
    let mut mol = smiles("CCCO");
    let mut conf = Conformer::new(4);
    for (i, (x, y, z)) in [
        (-1.2784, -0.0795, -0.4157),
        (-0.1858, 0.7166, 0.2689),
        (0.9073, -0.314, 0.5667),
        (1.2802, -0.8242, -0.6921),
    ]
    .into_iter()
    .enumerate()
    {
        conf.set_atom_position(i as u32, Point3D::new(x, y, z))
            .unwrap();
    }
    mol.add_conformer(&conf).unwrap();
    mol
}

#[test]
fn matching_3d_structure() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    let reference = propanol_3d();
    let mut mol = smiles("CCCO");
    generate_depiction_matching_3d_structure(&mut mol, &reference).unwrap();
    assert_eq!(mol.get_num_conformers(), 1);
    assert!(!mol.get_conformer().unwrap().is_3d());
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);

    let mut larger = smiles("CCCOC");
    assert!(is_value_error(
        generate_depiction_matching_3d_structure(&mut larger, &reference),
        "not compatible"
    ));
    let accept = GenerateDepictionMatching3DStructureParams {
        accept_failure: true,
        ..Default::default()
    };
    generate_depiction_matching_3d_structure_with_params(&mut larger, &reference, &accept).unwrap();
    assert_eq!(larger.get_num_conformers(), 1);

    let patterned = GenerateDepictionMatching3DStructureParams {
        ref_patt: Some(mol_from_smarts("CCO").unwrap()),
        ..Default::default()
    };
    let mut mol = smiles("CCCO");
    generate_depiction_matching_3d_structure_with_params(&mut mol, &reference, &patterned).unwrap();
    let missing = GenerateDepictionMatching3DStructureParams {
        conf_id: 3,
        ..Default::default()
    };
    assert!(is_value_error(
        generate_depiction_matching_3d_structure_with_params(&mut mol, &reference, &missing),
        "conformation"
    ));
    let mut infinite = propanol_3d();
    set_position(&mut infinite, 2, Point3D::new(f64::INFINITY, 0.0, 0.0));
    assert!(is_value_error(
        generate_depiction_matching_3d_structure(&mut mol, &infinite),
        "the reference's coordinates"
    ));
}

// The smaller of each bond's angle with the x axis, modulo 30 degrees, and
// 30 minus it.
fn off_30_degrees(mol: &Mol) -> Vec<f64> {
    mol.get_bonds()
        .map(|b| {
            let p = position(mol, -1, b.get_begin_atom_idx());
            let q = position(mol, -1, b.get_end_atom_idx());
            let angle = (p.y - q.y).atan2(p.x - q.x).to_degrees().rem_euclid(30.0);
            angle.min(30.0 - angle)
        })
        .collect()
}

#[test]
fn straighten_depiction_rotates_to_30_degrees() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    let mut mol = smiles("c1ccccc1CCO");
    compute_2d_coords(&mut mol).unwrap();
    // Rotate by 10 degrees first.
    let (sin, cos) = 10f64.to_radians().sin_cos();
    for i in 0..mol.get_num_atoms() {
        let p = position(&mol, -1, i);
        set_position(
            &mut mol,
            i,
            Point3D::new(p.x * cos - p.y * sin, p.x * sin + p.y * cos, 0.0),
        );
    }
    straighten_depiction(&mut mol).unwrap();
    assert!(off_30_degrees(&mol).iter().all(|&d| d < 1e-3));
    straighten_depiction_with_params(
        &mut mol,
        &StraightenDepictionParams {
            minimize_rotation: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(off_30_degrees(&mol).iter().all(|&d| d < 1e-3));
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);
}

#[test]
fn normalize_depiction_scales() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    let mut mol = smiles("c1ccccc1CCO");
    compute_2d_coords(&mut mol).unwrap();
    assert!((normalize_depiction(&mut mol).unwrap() - 1.0).abs() < 1e-9);
    let doubled = NormalizeDepictionParams {
        scale_factor: 2.0,
        ..Default::default()
    };
    assert!((normalize_depiction_with_params(&mut mol, &doubled).unwrap() - 2.0).abs() < 1e-9);
    assert_bond_lengths(&mol, -1, 3.0, 1e-3);
    // Back to RDKit's bond length, without reorienting.
    let unoriented = NormalizeDepictionParams {
        canonicalize: 0,
        ..Default::default()
    };
    let scale = normalize_depiction_with_params(&mut mol, &unoriented).unwrap();
    assert!((scale - 0.5).abs() < 1e-9);
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);
    let to_y = NormalizeDepictionParams {
        canonicalize: -1,
        ..Default::default()
    };
    normalize_depiction_with_params(&mut mol, &to_y).unwrap();
    assert_bond_lengths(&mol, -1, 1.5, 1e-3);
}

#[test]
fn straighten_and_normalize_refuse_bad_coordinates() {
    // Without bonds, RDKit returns before reading a conformer.
    let mut atom = smiles("C");
    assert_eq!(normalize_depiction(&mut atom).unwrap(), -1.0);
    straighten_depiction(&mut atom).unwrap();
    let mut mol = smiles("CCC");
    assert!(is_value_error(
        normalize_depiction(&mut mol),
        "No conformations"
    ));
    assert!(is_value_error(
        straighten_depiction(&mut mol),
        "No conformations"
    ));

    // NormalizeDepiction casts bond lengths times 10 to int, and
    // StraightenDepiction each bond's angle, which NaN makes undefined.
    for bad in [f64::NAN, f64::INFINITY, 1e12] {
        let mut mol = smiles("CCC");
        compute_2d_coords(&mut mol).unwrap();
        set_position(&mut mol, 0, Point3D::new(bad, 0.0, 0.0));
        assert!(is_value_error(
            normalize_depiction(&mut mol),
            "the coordinates"
        ));
        assert!(is_value_error(
            straighten_depiction(&mut mol),
            "the coordinates"
        ));
    }
}

#[test]
fn coord_gen_preference() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    assert!(!get_prefer_coord_gen());
    if !is_coord_gen_support_available() {
        set_prefer_coord_gen(true);
        assert!(!get_prefer_coord_gen());
        return;
    }
    {
        let _guard = UsingCoordGen::new(true);
        assert!(get_prefer_coord_gen());
        // CoordGen's bonds are 1.0 long; forceRDKit uses RDKit's 1.5.
        let mut mol = smiles("c1ccccc1CCO");
        assert_eq!(compute_2d_coords(&mut mol).unwrap(), 0);
        assert_bond_lengths(&mol, -1, 1.0, 1e-2);
        let forced = Compute2DCoordsParams {
            force_rdkit: true,
            ..Default::default()
        };
        compute_2d_coords_with_params(&mut mol, &forced).unwrap();
        assert_bond_lengths(&mol, -1, 1.5, 1e-3);
        {
            let _inner = UsingCoordGen::new(false);
            assert!(!get_prefer_coord_gen());
        }
        assert!(get_prefer_coord_gen());

        // CoordGen's canonical ordering sums Morgan scores in int, which a
        // 45-carbon chain overflows; RDKit's depictor is fine with it.
        let mut chain = smiles(&"C".repeat(45));
        assert!(is_value_error(compute_2d_coords(&mut chain), "Morgan"));
        compute_2d_coords_with_params(&mut chain, &forced).unwrap();
        let mut short = smiles(&"C".repeat(30));
        compute_2d_coords(&mut short).unwrap();
        assert_bond_lengths(&short, -1, 1.0, 1e-2);
    }
    assert!(!get_prefer_coord_gen());
}

#[test]
fn ring_system_templates() {
    let _settings = settings();
    set_prefer_coord_gen(false);
    // A bicyclo[2.2.1]heptane depiction from Python's MolToCXSmiles.
    let good = common::temp_path("templates.cxsmi");
    std::fs::write(
        &good,
        "C1CC2CCC1C2 |(1.55254,0.75,;1.55254,-0.75,;0.125954,-1.21353,;-1.30063,-0.75,;\
         -1.30063,0.75,;0.125954,1.21353,;-0.755724,0,)|\n",
    )
    .unwrap();
    set_ring_system_templates(&good).unwrap();
    add_ring_system_templates(&good).unwrap();
    let mut mol = smiles("C1CC2CCC1C2");
    let templated = Compute2DCoordsParams {
        use_ring_templates: true,
        ..Default::default()
    };
    compute_2d_coords_with_params(&mut mol, &templated).unwrap();

    // DepictException is a ValueError; RDKit's message names the file the
    // caller gave, not the checked copy RDKit reads.
    let chain = common::temp_path("chain.cxsmi");
    std::fs::write(&chain, "CCO |(0,0,;1,0,;2,0,)|\n").unwrap();
    let err = set_ring_system_templates(&chain).unwrap_err();
    assert!(
        matches!(&err, Error::Value { message } if message.contains("not a ring system")),
        "{err:?}"
    );
    let invalid = common::temp_path("invalid.cxsmi");
    std::fs::write(&invalid, "C1CC1(\n").unwrap();
    let message = format!("Could not load templates from {}", invalid.display());
    assert!(is_value_error(
        set_ring_system_templates(&invalid),
        &message
    ));
    let missing = common::temp_path("missing.cxsmi");
    let message = format!("Could not open file {}", missing.display());
    assert!(is_value_error(
        set_ring_system_templates(&missing),
        &message
    ));

    // RDKit copies template coordinates into depictions unchecked, and
    // CXSMILES coordinates are parsed with atof.
    let nan = common::temp_path("nan.cxsmi");
    std::fs::write(&nan, "C1CC1 |(nan,0,;1,0,;0,1,)|\n").unwrap();
    assert!(is_value_error(
        add_ring_system_templates(&nan),
        "ring system template coordinates"
    ));
    load_default_ring_system_templates();
    compute_2d_coords_with_params(&mut mol, &templated).unwrap();
    assert_bond_lengths(&mol, -1, 1.5, 0.1);
    for path in [good, chain, invalid, nan] {
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn depiction_settings_hold_for_a_whole_call_across_threads() {
    // One thread toggles the CoordGen preference while others depict with
    // RDKit's depictor at two bond lengths: each call sees its own bond
    // length throughout, and no call sees another's.
    let _settings = settings();
    let done = std::sync::atomic::AtomicBool::new(false);
    std::thread::scope(|scope| {
        scope.spawn(|| {
            let mut prefer = false;
            while !done.load(std::sync::atomic::Ordering::Relaxed) {
                prefer = !prefer;
                set_prefer_coord_gen(prefer);
                let _guard = UsingCoordGen::new(!prefer);
            }
            set_prefer_coord_gen(false);
        });
        let workers: Vec<_> = [1.5, 2.0, 1.5, 2.0]
            .into_iter()
            .map(|length| {
                scope.spawn(move || {
                    let params = Compute2DCoordsParams {
                        bond_length: length,
                        force_rdkit: true,
                        ..Default::default()
                    };
                    for text in ["c1ccccc1CCO", "CC(C)(C)C(=O)N", "C1CCC2CCCCC2C1"] {
                        for _ in 0..10 {
                            let mut mol = smiles(text);
                            compute_2d_coords_with_params(&mut mol, &params).unwrap();
                            assert_bond_lengths(&mol, -1, length, 1e-3);
                            // CoordGen or RDKit's depictor, whichever is
                            // preferred at the moment.
                            compute_2d_coords(&mut mol).unwrap();
                        }
                    }
                })
            })
            .collect();
        for worker in workers {
            worker.join().unwrap();
        }
        done.store(true, std::sync::atomic::Ordering::Relaxed);
    });
    assert!(!get_prefer_coord_gen());
}

#[test]
fn coord_gen_refuses_double_bond_stereo_without_stereo_atoms() {
    let _settings = settings();
    // RDKit 2024.09's CoordGen reads a Z double bond's two stereo atoms
    // from an empty list (RDKit 2025.03 skips such a bond).
    let mut mol = smiles("C=C1C=CC(=O)CC1");
    mol.get_bond_with_idx_mut(0)
        .unwrap()
        .set_stereo(rdkit_cxx::chem::BondStereo::Stereoz);
    set_prefer_coord_gen(true);
    let refused = compute_2d_coords(&mut mol);
    let mut plain = smiles("C=C1C=CC(=O)CC1");
    let fine = compute_2d_coords(&mut plain);
    set_prefer_coord_gen(false);
    assert!(is_value_error(refused, "fewer than two stereo atoms"));
    fine.unwrap();
}
