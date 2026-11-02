// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use rdkit_cxx::Error;
use rdkit_cxx::data_structs::{ExplicitBitVect, SparseBitVect};
use rdkit_cxx::ml::info_theory::*;

// Expected values come from Python's rdInfoTheory.

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs().max(1.0)
}

fn assert_rows(actual: &[Vec<f64>], expected: &[&[f64]]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert_eq!(a.len(), e.len(), "{actual:?}");
        assert!(a.iter().zip(*e).all(|(x, y)| close(*x, *y)), "{actual:?}");
    }
}

fn ebv(size: u32, bits: &[u32]) -> ExplicitBitVect {
    let mut bv = ExplicitBitVect::new(size);
    bv.set_bits_from_list(bits).unwrap();
    bv
}

fn sbv(size: u32, bits: &[u32]) -> SparseBitVect {
    let mut bv = SparseBitVect::new(size);
    bv.set_bits_from_list(bits).unwrap();
    bv
}

// RDKit's docstring example: four fingerprints, two classes.
fn examples() -> Vec<(ExplicitBitVect, u32)> {
    vec![
        (ebv(4, &[3]), 0),
        (ebv(4, &[1, 3]), 0),
        (ebv(4, &[2]), 1),
        (ebv(4, &[0, 1, 2]), 1),
    ]
}

fn ranker(info_type: InfoType) -> InfoBitRanker {
    let mut ranker = InfoBitRanker::new_with_info_type(4, 2, info_type);
    if matches!(info_type, InfoType::Biasentropy | InfoType::Biaschisquare) {
        ranker.set_bias_list(&[1]).unwrap();
    }
    for (bv, label) in examples() {
        ranker.accumulate_votes(&bv, label).unwrap();
    }
    ranker
}

#[test]
fn statistics() {
    assert!(close(
        info_entropy(&[1.0, 2.0, 3.0]).unwrap(),
        1.459147917027245
    ));
    assert_eq!(info_entropy(&[0.0, 0.0]).unwrap(), 0.0);
    assert_eq!(info_entropy(&[5.0]).unwrap(), 0.0);
    let m = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![0.0, 5.0]];
    assert!(close(info_gain(&m).unwrap(), 0.19320844498095213));
    assert!(close(chi_square(&m).unwrap(), 2.824675324675326));
    let m = vec![vec![10.0, 0.0, 3.0], vec![2.0, 8.0, 1.0]];
    assert!(close(info_gain(&m).unwrap(), 0.5347605966186039));
    assert!(close(chi_square(&m).unwrap(), 14.265734265734267));
    let m = vec![vec![1.5, 2.25], vec![0.5, 3.75]];
    assert!(close(info_gain(&m).unwrap(), 0.07853536558575702));
    assert!(close(chi_square(&m).unwrap(), 0.7411764705882362));
    assert_eq!(info_gain(&[]).unwrap(), 0.0);
    assert_eq!(chi_square(&[]).unwrap(), 0.0);
    assert_eq!(info_gain(&[vec![], vec![]]).unwrap(), 0.0);
    assert!(chi_square(&[vec![], vec![]]).unwrap().is_nan());
}

#[test]
fn statistics_errors() {
    assert!(matches!(info_entropy(&[]), Err(Error::Runtime { .. })));
    let ragged = vec![vec![1.0], vec![1.0, 2.0]];
    assert!(matches!(info_gain(&ragged), Err(Error::Value { .. })));
    assert!(matches!(chi_square(&ragged), Err(Error::Value { .. })));
}

#[test]
fn statistics_guards() {
    // RDKit truncates the sums to int and adds those in int: garbage in
    // Python for these.
    for bad in [1e12, f64::NAN, f64::INFINITY, -3e9] {
        let m = vec![vec![bad, 2.0], vec![3.0, 4.0]];
        assert!(matches!(chi_square(&m), Err(Error::Value { .. })), "{bad}");
        assert!(matches!(info_gain(&m), Err(Error::Value { .. })), "{bad}");
    }
    // Each sum fits, but their sum does not.
    let m = vec![vec![2e9], vec![2e9]];
    assert!(matches!(chi_square(&m), Err(Error::Value { .. })));
    let m = vec![vec![2e9, 2e9]];
    assert!(matches!(info_gain(&m), Err(Error::Value { .. })));
    // Sums just inside int's range are taken.
    assert!(chi_square(&[vec![2147483647.5], vec![0.0]]).is_ok());
    assert!(info_gain(&[vec![-2147483648.5, 0.0]]).is_ok());
}

#[test]
fn ranking() {
    assert_rows(
        &ranker(InfoType::Entropy).get_top_n(3).unwrap(),
        &[
            &[3.0, 1.0, 2.0, 0.0],
            &[2.0, 1.0, 0.0, 2.0],
            &[0.0, 0.31127812445913294, 0.0, 1.0],
        ],
    );
    assert_rows(
        &ranker(InfoType::Biasentropy).get_top_n(3).unwrap(),
        &[
            &[2.0, 1.0, 0.0, 2.0],
            &[0.0, 0.31127812445913294, 0.0, 1.0],
            &[1.0, 0.0, 1.0, 1.0],
        ],
    );
    assert_rows(
        &ranker(InfoType::Chisquare).get_top_n(3).unwrap(),
        &[
            &[3.0, 4.0, 2.0, 0.0],
            &[2.0, 4.0, 0.0, 2.0],
            &[0.0, 1.333333333333333, 0.0, 1.0],
        ],
    );
    assert_rows(
        &ranker(InfoType::Biaschisquare).get_top_n(3).unwrap(),
        &[
            &[2.0, 4.0, 0.0, 2.0],
            &[0.0, 1.333333333333333, 0.0, 1.0],
            &[1.0, 0.0, 1.0, 1.0],
        ],
    );
    let mut sparse = InfoBitRanker::new(4, 2);
    for (bv, label) in examples() {
        sparse
            .accumulate_votes(&sbv(4, &bv.get_on_bits()), label)
            .unwrap();
    }
    assert_rows(
        &sparse.get_top_n(4).unwrap(),
        &[
            &[2.0, 1.0, 0.0, 2.0],
            &[3.0, 1.0, 2.0, 0.0],
            &[0.0, 0.31127812445913294, 0.0, 1.0],
            &[1.0, 0.0, 1.0, 1.0],
        ],
    );
    assert_rows(
        &InfoBitRanker::new(4, 2).get_top_n(2).unwrap(),
        &[&[1.0, 0.0, 0.0, 0.0], &[0.0, 0.0, 0.0, 0.0]],
    );
}

#[test]
fn masks() {
    let mut ranker = InfoBitRanker::new(4, 2);
    ranker.set_mask_bits(&[1, 2, 3]).unwrap();
    for (bv, label) in examples() {
        ranker.accumulate_votes(&bv, label).unwrap();
    }
    assert_rows(
        &ranker.get_top_n(3).unwrap(),
        &[
            &[3.0, 1.0, 2.0, 0.0],
            &[2.0, 1.0, 0.0, 2.0],
            &[1.0, 0.0, 1.0, 1.0],
        ],
    );
    let mut ranker = InfoBitRanker::new(6, 2);
    ranker.set_mask_bits(&[1, 4]).unwrap();
    for (bv, label) in examples() {
        ranker
            .accumulate_votes(&ebv(6, &bv.get_on_bits()), label)
            .unwrap();
    }
    assert_rows(
        &ranker.get_top_n(2).unwrap(),
        &[&[4.0, 0.0, 0.0, 0.0], &[1.0, 0.0, 1.0, 1.0]],
    );
    // As Python raises.
    assert!(matches!(ranker.get_top_n(7), Err(Error::Value { .. })));
    assert!(matches!(ranker.get_top_n(3), Err(Error::Runtime { .. })));
    assert!(matches!(
        ranker.accumulate_votes(&ebv(6, &[1]), 2),
        Err(Error::Runtime { .. })
    ));
    assert!(matches!(
        ranker.accumulate_votes(&ebv(5, &[1]), 0),
        Err(Error::Runtime { .. })
    ));
    // A bit past the ranker's, which leaves the mask as it was.
    assert!(matches!(
        ranker.set_mask_bits(&[0, 6]),
        Err(Error::Index { index: 6, len: 6 })
    ));
    assert!(ranker.get_top_n(2).is_ok());
}

#[test]
fn bias_lists() {
    let mut ranker = InfoBitRanker::new_with_info_type(8, 3, InfoType::Biasentropy);
    // As Python raises; RDKit's own check keeps the refused list, which a
    // biased ranking then reads past the classes with (a segfault in
    // Python).
    for list in [&[0, 2_000_000_000][..], &[0, 0], &[0, 1, 2], &[3]] {
        assert!(
            matches!(ranker.set_bias_list(list), Err(Error::Runtime { .. })),
            "{list:?}"
        );
    }
    ranker.accumulate_votes(&ebv(8, &[]), 0).unwrap();
    ranker.accumulate_votes(&ebv(8, &[1, 2]), 2).unwrap();
    let top = ranker.get_top_n(2).unwrap();
    assert_eq!(top.len(), 2);
    ranker.set_bias_list(&[2, 0]).unwrap();
    assert!(ranker.get_top_n(2).is_ok());
}

#[test]
fn top_n_guards() {
    // No bits to rank: RDKit reads its empty queue (a segfault in Python).
    assert!(matches!(
        InfoBitRanker::new(4, 2).get_top_n(0),
        Err(Error::Value { .. })
    ));
    // (The num * (2 + n_classes) bound needs a ranker of 4 GiB of counts.)
}

#[test]
#[should_panic(expected = "at most 32768 classes")]
fn too_many_classes() {
    InfoBitRanker::new(4, 32769);
}

#[test]
#[should_panic(expected = "at most i32::MAX bits")]
fn too_many_bits() {
    InfoBitRanker::new(1 << 31, 2);
}

#[test]
fn write_top_bits() {
    let mut ranker = InfoBitRanker::new(4, 2);
    let path = common::temp_path("top_bits.txt");
    ranker.write_top_bits_to_file(&path).unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "         Bit InfoContent     class0     class1\n"
    );
    for (bv, label) in examples() {
        ranker
            .accumulate_votes(&sbv(4, &bv.get_on_bits()), label)
            .unwrap();
    }
    ranker.get_top_n(4).unwrap();
    ranker.write_top_bits_to_file(&path).unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "         Bit InfoContent     class0     class1\n           2           1         0         2\n           3           1         2         0\n           0     0.31128         0         1\n           1           0         1         1\n"
    );
    std::fs::remove_file(&path).unwrap();
    assert!(matches!(
        ranker.write_top_bits_to_file("/nonexistent/dir/top.txt"),
        Err(Error::Os { .. })
    ));
}

#[cfg(target_os = "linux")]
#[test]
fn write_top_bits_failure() {
    // RDKit ignores a failed write.
    let ranker = InfoBitRanker::new(4, 2);
    assert!(matches!(
        ranker.write_top_bits_to_file("/dev/full"),
        Err(Error::Os { .. })
    ));
}

#[test]
fn correlations() {
    let mut generator = BitCorrMatGenerator::default();
    assert!(generator.get_corr_matrix().is_empty());
    generator.set_bit_list(&[0, 3, 5, 6]).unwrap();
    generator.collect_votes(&ebv(8, &[0, 3, 5])).unwrap();
    generator.collect_votes(&ebv(8, &[3, 5, 6])).unwrap();
    generator.collect_votes(&ebv(8, &[0, 6])).unwrap();
    generator.collect_votes(&sbv(8, &[0, 3, 5, 6])).unwrap();
    assert_eq!(generator.get_corr_matrix(), [2.0, 2.0, 3.0, 2.0, 2.0, 2.0]);
    // A listed bit past the fingerprint, as Python raises; nothing counts.
    assert!(matches!(
        generator.collect_votes(&ebv(4, &[0])),
        Err(Error::Index { index: 5, len: 4 })
    ));
    assert_eq!(generator.get_corr_matrix(), [2.0, 2.0, 3.0, 2.0, 2.0, 2.0]);
    generator.set_bit_list(&[1]).unwrap();
    assert!(generator.get_corr_matrix().is_empty());
}

#[test]
fn correlation_guards() {
    // RDKit sizes the matrix in int (bad_array_new_length in Python).
    let mut generator = BitCorrMatGenerator::default();
    let bits: Vec<u32> = (0..46342).collect();
    assert!(matches!(
        generator.set_bit_list(&bits),
        Err(Error::Value { .. })
    ));
    // The refused list leaves the generator as it was.
    assert!(generator.get_corr_matrix().is_empty());
    generator.set_bit_list(&bits[..1000]).unwrap();
    assert_eq!(generator.get_corr_matrix().len(), 1000 * 999 / 2);
}
