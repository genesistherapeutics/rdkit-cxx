// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;
use std::path::PathBuf;

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;
use rdkit_cxx::data_structs::*;

fn bits(size: u32, on: &[u32]) -> ExplicitBitVect {
    let mut bv = ExplicitBitVect::new(size);
    bv.set_bits_from_list(on).unwrap();
    bv
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-6
}

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

#[test]
fn explicit_bit_vects() {
    let mut a = bits(32, &[1, 3, 5, 7, 9]);
    assert_eq!(a.get_num_bits(), 32);
    assert_eq!(a.get_num_on_bits(), 5);
    assert_eq!(a.get_num_off_bits(), 27);
    assert_eq!(a.get_on_bits(), [1, 3, 5, 7, 9]);
    assert!(a.get_bit(3).unwrap());
    assert!(a.set_bit(3).unwrap());
    assert!(!a.un_set_bit(4).unwrap());
    assert_eq!(
        a.get_bit(32).unwrap_err(),
        Error::Index { index: 32, len: 32 }
    );
    // An out-of-range bit leaves the vector as it was.
    assert!(a.set_bits_from_list(&[0, 40]).is_err());
    assert!(!a.get_bit(0).unwrap());
    assert_eq!(a.to_list().iter().filter(|&&on| on).count(), 5);

    assert_eq!(
        ExplicitBitVect::from_pkl(&a.to_binary().unwrap()).unwrap(),
        a
    );
    assert_eq!(a.to_base64().unwrap(), "4P///yAAAAAFAAAAAgICAgIs");
    let mut copy = ExplicitBitVect::default();
    copy.from_base64(&a.to_base64().unwrap()).unwrap();
    assert_eq!(copy, a);
    assert_eq!(
        ExplicitBitVect::new_with_bits_set(4, true).get_num_on_bits(),
        4
    );

    let b = bits(32, &[1, 2, 3, 4]);
    assert_eq!((&a & &b).get_on_bits(), [1, 3]);
    assert_eq!((&a | &b).get_on_bits(), [1, 2, 3, 4, 5, 7, 9]);
    assert_eq!((&a ^ &b).get_on_bits(), [2, 4, 5, 7, 9]);
    assert_eq!((!&a).get_num_on_bits(), 27);
    assert_eq!((&a + &b).get_num_bits(), 64);
    let mut grown = a.clone();
    grown += &b;
    assert_eq!(grown, &a + &b);
}

#[test]
#[should_panic(expected = "BitVects must be same length")]
fn operators_on_vectors_of_different_sizes_panic() {
    let _ = &bits(32, &[1]) & &bits(64, &[1]);
}

#[test]
fn sparse_bit_vects() {
    let mut s = SparseBitVect::new(1000);
    s.set_bits_from_list(&[10, 500]).unwrap();
    assert_eq!(s.get_num_on_bits(), 2);
    assert_eq!(convert_to_explicit(&s).unwrap().get_on_bits(), [10, 500]);
    assert_eq!(SparseBitVect::from_pkl(&s.to_binary().unwrap()).unwrap(), s);
    let t = SparseBitVect::new(1000);
    assert_eq!((&s | &t), s);
    assert!(close(tanimoto_similarity(&s, &s).unwrap(), 1.0));
}

#[test]
fn similarity_metrics() {
    let a = bits(32, &[1, 3, 5, 7, 9]);
    let b = bits(32, &[1, 2, 3, 4]);
    let expected = [
        (tanimoto_similarity(&a, &b), 0.285714),
        (cosine_similarity(&a, &b), 0.447214),
        (kulczynski_similarity(&a, &b), 0.45),
        (dice_similarity(&a, &b), 0.444444),
        (sokal_similarity(&a, &b), 0.166667),
        (mc_connaughey_similarity(&a, &b), -0.1),
        (asymmetric_similarity(&a, &b), 0.5),
        (braun_blanquet_similarity(&a, &b), 0.4),
        (russel_similarity(&a, &b), 0.0625),
        (rogot_goldberg_similarity(&a, &b), 0.676768),
        (on_bit_similarity(&a, &b), 0.285714),
        (all_bit_similarity(&a, &b), 0.84375),
        (tversky_similarity(&a, &b, 0.5, 0.5), 0.444444),
    ];
    for (value, expected) in expected {
        assert!(close(value.unwrap(), expected), "{expected}");
    }
    assert!(close(
        tanimoto_similarity_with_return_distance(&a, &b, true).unwrap(),
        0.714286
    ));
    assert!(close(
        tanimoto_similarity_bytes(&a, &b.to_binary().unwrap()).unwrap(),
        0.285714
    ));
    let bulk = bulk_tanimoto_similarity(&a, [&a, &b]).unwrap();
    assert!(close(bulk[0], 1.0) && close(bulk[1], 0.285714));
    let distances =
        bulk_dice_similarity_with_return_distance(&a, &[a.clone(), b.clone()], true).unwrap();
    assert!(close(distances[0], 0.0) && close(distances[1], 0.555556));
    assert_eq!(
        tanimoto_similarity_neighbors([&a, &b], [&b, &a]).unwrap(),
        [(1, 1.0), (0, 1.0)]
    );
    assert!(matches!(
        tanimoto_similarity(&a, &bits(64, &[])),
        Err(Error::Value { .. })
    ));

    let on = on_bit_proj_similarity(&a, &b).unwrap();
    assert!(close(on[0], 0.4) && close(on[1], 0.5));
    assert_eq!(num_bits_in_common(&a, &b).unwrap(), 27);
    assert_eq!(on_bits_in_common(&a, &b).unwrap(), [1, 3]);
    assert_eq!(off_bits_in_common(&a, &b).unwrap().len(), 25);
    let folded = fold_fingerprint(&a).unwrap();
    assert_eq!(folded.get_num_bits(), 16);
    assert_eq!(folded.get_on_bits(), [1, 3, 5, 7, 9]);
    assert!(matches!(
        fold_fingerprint_with_fold_factor(&a, 0),
        Err(Error::Value { .. })
    ));
    assert!(!all_probe_bits_match(&b, &a).unwrap());
    assert!(all_probe_bits_match(&ExplicitBitVect::new(32), &a).unwrap());
    assert!(all_probe_bits_match_bytes(&bits(32, &[1]), &a.to_binary().unwrap()).unwrap());
}

#[test]
fn text_forms() {
    let b = bits(32, &[1, 2, 3, 4]);
    assert!(bit_vect_to_text(&b).starts_with("01111000"));
    assert_eq!(bit_vect_to_fps_text(&b), "1e000000");
    assert_eq!(bit_vect_to_binary_text(&b), b"\x1e\x00\x00\x00");
    assert_eq!(create_from_fps_text("1e000000").unwrap(), b);
    assert!(matches!(
        create_from_fps_text("1e0"),
        Err(Error::Value { .. })
    ));
    assert_eq!(create_from_binary_text(b"\x1e\x00\x00\x00"), b);
    assert_eq!(create_from_bit_string("0101").get_on_bits(), [1, 3]);
}

#[test]
fn sparse_int_vects() {
    let mut v = UIntSparseIntVect::new(10).unwrap();
    v.set_val(2, 3).unwrap();
    v.set_val(5, 1).unwrap();
    let mut w = UIntSparseIntVect::new(10).unwrap();
    w.set_val(2, 1).unwrap();
    w.set_val(7, 4).unwrap();
    assert_eq!(v.get_val(2).unwrap(), 3);
    assert_eq!(
        v.get_val(10).unwrap_err(),
        Error::Index { index: 10, len: 10 }
    );
    assert_eq!(v.get_total_val().unwrap(), 4);
    assert_eq!(v.get_length(), 10);
    assert_eq!(v.get_nonzero_elements(), BTreeMap::from([(2, 3), (5, 1)]));
    assert_eq!(UIntSparseIntVect::from_pkl(&v.to_binary()).unwrap(), v);

    assert!(close(
        dice_similarity_sparse_int_vect(&v, &w).unwrap(),
        0.222222
    ));
    assert!(close(
        tanimoto_similarity_sparse_int_vect(&v, &w).unwrap(),
        0.125
    ));
    assert!(close(
        tversky_similarity_sparse_int_vect(&v, &w, 0.5, 0.5).unwrap(),
        0.222222
    ));
    let params = SparseIntVectSimilarityParams {
        return_distance: true,
        ..Default::default()
    };
    assert!(close(
        tanimoto_similarity_sparse_int_vect_with_params(&v, &w, &params).unwrap(),
        0.875
    ));
    let bulk = bulk_dice_similarity_sparse_int_vect(&v, [&v, &w]).unwrap();
    assert!(close(bulk[0], 1.0) && close(bulk[1], 0.222222));

    assert_eq!((&v & &w).get_nonzero_elements(), BTreeMap::from([(2, 1)]));
    assert_eq!(
        (&v | &w).get_nonzero_elements(),
        BTreeMap::from([(2, 3), (5, 1), (7, 4)])
    );
    assert_eq!(
        (&v + &w).get_nonzero_elements(),
        BTreeMap::from([(2, 4), (5, 1), (7, 4)])
    );
    assert_eq!(
        (&v - &w).get_nonzero_elements(),
        BTreeMap::from([(2, 2), (5, 1), (7, -4)])
    );
    let mut sum = v.clone();
    sum += &w;
    assert_eq!(sum, &v + &w);

    let mut counts = IntSparseIntVect::new(5).unwrap();
    counts.update_from_sequence(&[1, 1, 3]).unwrap();
    counts *= 3;
    assert_eq!(counts.to_list(), [0, 6, 0, 3, 0]);
    counts /= 3;
    counts += 1;
    counts -= 2;
    assert_eq!(counts.to_list(), [0, 1, 0, 0, 0]);
    assert!(counts.update_from_sequence(&[-1]).is_err());
    assert!(matches!(
        IntSparseIntVect::new(-1),
        Err(Error::Value { .. })
    ));
    assert_eq!(LongSparseIntVect::default().get_length(), 0);
    let mut big = ULongSparseIntVect::new(1 << 40).unwrap();
    big.set_val(1 << 39, 7).unwrap();
    assert_eq!(big.get_nonzero_elements(), BTreeMap::from([(1 << 39, 7)]));
}

#[test]
#[should_panic(expected = "division by zero")]
fn dividing_by_zero_panics() {
    let mut v = IntSparseIntVect::new(2).unwrap();
    v.set_val(0, 1).unwrap();
    v /= 0;
}

#[test]
#[should_panic(expected = "overflow")]
fn overflow_panics() {
    let mut v = IntSparseIntVect::new(2).unwrap();
    v.set_val(0, i32::MAX).unwrap();
    v += 1;
}

#[test]
fn value_vects() {
    let mut d = DiscreteValueVect::new(DiscreteValueType::Twobitvalue, 5);
    d.set_val(0, 3).unwrap();
    d.set_val(1, 2).unwrap();
    let mut e = DiscreteValueVect::new(DiscreteValueType::Twobitvalue, 5);
    e.set_val(0, 1).unwrap();
    e.set_val(1, 3).unwrap();
    assert_eq!(d.get_value_type(), DiscreteValueType::Twobitvalue);
    assert_eq!(d.get_total_val(), 5);
    assert_eq!(d.len(), 5);
    assert_eq!(compute_l1_norm(&d, &e).unwrap(), 3);
    assert_eq!(convert_to_numpy_array(&(&d & &e)), [1, 2, 0, 0, 0]);
    assert_eq!(convert_to_numpy_array(&(&d + &e)), [3, 3, 0, 0, 0]);
    assert_eq!(convert_to_numpy_array(&(&d - &e)), [2, 0, 0, 0, 0]);
    assert!(matches!(d.set_val(0, 5), Err(Error::Value { .. })));
    assert_eq!(d.get_val(5).unwrap_err(), Error::Index { index: 5, len: 5 });
    assert_eq!(DiscreteValueVect::from_pkl(&d.to_binary()).unwrap(), d);
    assert_eq!(DiscreteValueVect::default().len(), 0);

    assert_eq!(convert_to_numpy_array(&bits(4, &[1])), [0, 1, 0, 0]);
}

#[test]
fn atom_bit_vect_props() {
    let mut atom = Atom::from_num(6);
    let bv = bits(8, &[2]);
    atom.set_explicit_bit_vect_prop("fp", &bv);
    assert_eq!(atom.get_explicit_bit_vect_prop("fp").unwrap(), bv);
    assert!(matches!(
        atom.get_explicit_bit_vect_prop("missing"),
        Err(Error::Key { .. })
    ));
    atom.set_int_prop("n", 1);
    assert!(matches!(
        atom.get_explicit_bit_vect_prop("n"),
        Err(Error::Value { .. })
    ));
}

#[test]
fn fpb_reader() {
    let mut fpbr = FPBReader::new(data("zim.head100.fpb")).unwrap();
    assert!(fpbr.len().is_err());
    fpbr.init().unwrap();
    assert_eq!(fpbr.len().unwrap(), 100);
    assert_eq!(fpbr.get_num_bits().unwrap(), 2048);
    assert_eq!(fpbr.get_id(0).unwrap(), "ZINC00902219");
    assert_eq!(fpbr.get_id(3).unwrap(), "ZINC04803506");
    let fp = fpbr.get_fp(0).unwrap();
    assert_eq!(fp.get_num_bits(), 2048);
    assert_eq!(fp.get_num_on_bits(), 17);
    let (fp0, name) = fpbr.get(0).unwrap();
    assert_eq!((fp0, name.as_str()), (fp, "ZINC00902219"));
    assert_eq!(
        fpbr.get_id(100).unwrap_err(),
        Error::Index {
            index: 100,
            len: 100
        }
    );

    let bv = fpbr.get_bytes(0).unwrap();
    assert_eq!(bv, bit_vect_to_binary_text(&fpbr.get_fp(0).unwrap()));
    assert!(close(fpbr.get_tanimoto(0, &bv).unwrap(), 1.0));
    assert!((fpbr.get_tanimoto(1, &bv).unwrap() - 0.3704).abs() < 1e-4);
    let nbrs = fpbr.get_tanimoto_neighbors(&bv).unwrap();
    assert_eq!(nbrs.len(), 1);
    assert_eq!(nbrs[0].1, 0);
    let nbrs = fpbr
        .get_tanimoto_neighbors_with_threshold(&bv, 0.3)
        .unwrap();
    assert_eq!(nbrs.len(), 5);
    assert_eq!(nbrs[1].1, 1);
    assert!((fpbr.get_tversky(1, &bv, 1.0, 1.0).unwrap() - 0.3704).abs() < 1e-4);
    assert_eq!(fpbr.get_tversky_neighbors(&bv, 1.0, 1.0).unwrap().len(), 1);
    assert_eq!(fpbr.get_containing_neighbors(&bv).unwrap(), [0]);
    let bv1 = fpbr.get_bytes(1).unwrap();
    assert_eq!(fpbr.get_containing_neighbors(&bv1).unwrap(), [1, 2, 3, 4]);

    // Fingerprints of the wrong length are refused, not read past their end.
    assert!(matches!(
        fpbr.get_tanimoto(0, &bv[..10]),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        FPBReader::new(data("missing.fpb")),
        Err(Error::Os { .. })
    ));
}

#[test]
fn multi_fpb_reader() {
    let mut mfpbr = MultiFPBReader::new();
    for (i, name) in ["zinc_random200.1.patt.fpb", "zinc_random200.2.patt.fpb"]
        .iter()
        .enumerate()
    {
        let reader = FPBReader::new(data(name)).unwrap();
        assert_eq!(mfpbr.add_reader(reader), i as u32 + 1);
    }
    mfpbr.init().unwrap();
    assert_eq!(mfpbr.get_num_bits().unwrap(), 1024);
    assert_eq!(mfpbr.len(), 2);
    assert_eq!(mfpbr.get_reader(0).unwrap().len().unwrap(), 200);
    assert!(mfpbr.get_reader(2).is_err());

    let query = mfpbr.get_reader(1).unwrap().get_bytes(5).unwrap();
    let params = MultiFPBReaderNeighborsParams {
        threshold: 0.99,
        num_threads: 2,
    };
    let nbrs = mfpbr
        .get_tanimoto_neighbors_with_params(&query, &params)
        .unwrap();
    assert!(nbrs.contains(&(1.0, 5, 1)), "{nbrs:?}");
    assert!(
        mfpbr
            .get_containing_neighbors(&query)
            .unwrap()
            .contains(&(5, 1))
    );
}
