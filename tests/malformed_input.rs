// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Inputs RDKit misreads without checking: malformed pickles and text, and
//! sizes past what its `int` arithmetic holds. Each must fail cleanly.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_coord_gen::add_coords;
use rdkit_cxx::chem::rd_mol_descriptors::calc_mol_formula;
use rdkit_cxx::chem::rd_mol_hash::{HashFunction, mol_hash};
use rdkit_cxx::chem::{Atom, RWMol, ResonanceMolSupplier, get_formal_charge, mol_to_cml_block};
use rdkit_cxx::data_structs::*;
use rdkit_cxx::geometry::*;

fn is_value(result: Result<impl std::fmt::Debug, Error>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

fn with_i32(mut pkl: Vec<u8>, offset: usize, value: i32) -> Vec<u8> {
    pkl[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    pkl
}

fn with_u32(pkl: Vec<u8>, offset: usize, value: u32) -> Vec<u8> {
    with_i32(pkl, offset, value as i32)
}

fn with_f64(mut pkl: Vec<u8>, offset: usize, value: f64) -> Vec<u8> {
    pkl[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
    pkl
}

fn bits(size: u32, on: &[u32]) -> ExplicitBitVect {
    let mut bv = ExplicitBitVect::new(size);
    bv.set_bits_from_list(on).unwrap();
    bv
}

#[test]
fn bit_vect_pickles() {
    let pkl = bits(32, &[1, 5]).to_binary().unwrap();
    assert_eq!(ExplicitBitVect::from_pkl(&pkl).unwrap(), bits(32, &[1, 5]));
    // RDKit negates the version, which overflows for i32::MIN.
    assert!(is_value(ExplicitBitVect::from_pkl(&with_i32(
        pkl.clone(),
        0,
        i32::MIN
    ))));
    assert!(is_value(SparseBitVect::from_pkl(&with_i32(
        pkl.clone(),
        0,
        i32::MIN
    ))));
    // A negative size would be a huge unsigned one.
    assert!(is_value(ExplicitBitVect::from_pkl(&with_i32(
        pkl.clone(),
        4,
        -5
    ))));
    assert!(is_value(SparseBitVect::from_pkl(&with_i32(
        pkl.clone(),
        4,
        -5
    ))));
    assert!(ExplicitBitVect::from_pkl(&[]).is_err());
    assert!(ExplicitBitVect::from_pkl(&pkl[..6]).is_err());
    let bad = with_i32(pkl.clone(), 0, i32::MIN);
    assert!(is_value(tanimoto_similarity_bytes(&bits(32, &[1]), &bad)));
}

#[test]
fn base64() {
    let a = bits(32, &[1, 5]);
    let mut b = ExplicitBitVect::new(8);
    b.from_base64(&a.to_base64().unwrap()).unwrap();
    assert_eq!(b, a);
    // RDKit's decoder writes past its buffer for these lengths, and its
    // in-place parse leaves the vector broken when it fails.
    let before = bits(16, &[3]);
    for text in ["", "A", "AAA", "AAAAAAA", "\0AAA"] {
        let mut c = before.clone();
        assert!(c.from_base64(text).is_err(), "{text:?}");
        assert_eq!(c, before);
        assert_eq!(c.get_on_bits(), [3]);
    }
}

#[test]
fn probe_bits() {
    let probe = bits(16, &[1, 3]);
    let reference = bits(16, &[1, 2, 3]);
    assert!(all_probe_bits_match(&probe, &reference).unwrap());
    assert!(all_probe_bits_match_bytes(&probe, &reference.to_binary().unwrap()).unwrap());
    assert!(!all_probe_bits_match_bytes(&reference, &probe.to_binary().unwrap()).unwrap());
    // The pickle may be another size, as in RDKit.
    assert!(
        all_probe_bits_match_bytes(&probe, &bits(64, &[1, 3, 40]).to_binary().unwrap()).unwrap()
    );
    assert!(!all_probe_bits_match_bytes(&bits(64, &[40]), &probe.to_binary().unwrap()).unwrap());
    // RDKit reads the pickle through raw pointers, past short ones.
    assert!(all_probe_bits_match_bytes(&probe, &[1, 2]).is_err());
    assert!(all_probe_bits_match_bytes(&probe, &[]).is_err());
    // RDKit's ExplicitBitVect overloads of these have no size check.
    let other = bits(8, &[1]);
    assert!(is_value(all_probe_bits_match(&probe, &other)));
    assert!(is_value(num_bits_in_common(&probe, &other)));
}

#[test]
fn daylight_strings() {
    let mut bv = ExplicitBitVect::new(64);
    for text in ["", "\n", "....", "...2", "....1x", "....12"] {
        assert!(
            is_value(init_from_daylight_string(&mut bv, text)),
            "{text:?}"
        );
    }
    init_from_daylight_string(&mut bv, "....3\n").unwrap();
    assert_eq!(bv.get_num_on_bits(), 0);
    // A lone pad character is an empty fingerprint.
    init_from_daylight_string(&mut bv, "1").unwrap();
}

#[test]
#[should_panic(expected = "i32::MAX bits")]
fn explicit_bit_vect_size() {
    ExplicitBitVect::new(1 << 31);
}

#[test]
#[should_panic(expected = "i32::MAX bits")]
fn sparse_bit_vect_size() {
    SparseBitVect::new(u32::MAX);
}

#[test]
#[should_panic(expected = "i32::MAX values")]
fn discrete_value_vect_size() {
    DiscreteValueVect::new(DiscreteValueType::Onebitvalue, u32::MAX);
}

#[test]
fn discrete_value_vect_pickles() {
    let mut dvv = DiscreteValueVect::new(DiscreteValueType::Twobitvalue, 20);
    dvv.set_val(3, 2).unwrap();
    let pkl = dvv.to_binary();
    assert_eq!(DiscreteValueVect::from_pkl(&pkl).unwrap(), dvv);
    // Header: version, type, bits, mask, length, number of ints.
    for bad in [
        with_i32(pkl.clone(), 0, i32::MIN),
        with_u32(pkl.clone(), 4, 7),
        with_u32(pkl.clone(), 8, 0),
        with_u32(pkl.clone(), 8, 4),
        with_u32(pkl.clone(), 12, 0xff),
        with_u32(pkl.clone(), 16, 1 << 31),
        with_u32(pkl.clone(), 16, 100),
        with_u32(pkl.clone(), 20, 0),
        pkl[..pkl.len() - 1].to_vec(),
    ] {
        assert!(is_value(DiscreteValueVect::from_pkl(&bad)), "{bad:?}");
    }
    // RDKit reads a version it may not have, and divides by the bits per
    // value before reading the rest, so every truncated header is refused.
    for len in 0..24 {
        assert!(
            is_value(DiscreteValueVect::from_pkl(&pkl[..len])),
            "{len} bytes"
        );
    }
    assert!(is_value(DiscreteValueVect::from_pkl(&ZERO_BITS_PER_VALUE)));
}

/// A DiscreteValueVect pickle's version, a value type of 0, and 0 bits per
/// value, which RDKit divides by.
const ZERO_BITS_PER_VALUE: [u8; 12] = [0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0];

#[test]
fn real_value_vect_pickles() {
    let mut rvv = RealValueVect::new(4);
    rvv.set_val(1, 2.5).unwrap();
    let pkl = rvv.to_binary();
    assert_eq!(
        RealValueVect::from_pkl(&pkl).unwrap().get_val(1).unwrap(),
        2.5
    );
    assert!(is_value(RealValueVect::from_pkl(&with_i32(
        pkl,
        0,
        i32::MIN
    ))));
}

#[test]
fn sparse_int_vect_pickles() {
    let mut siv = IntSparseIntVect::new(10).unwrap();
    siv.set_val(3, 5).unwrap();
    // Version, index size, length, number of entries, then index and value.
    let pkl = siv.to_binary();
    assert_eq!(IntSparseIntVect::from_pkl(&pkl).unwrap(), siv);
    assert!(is_value(IntSparseIntVect::from_pkl(&with_i32(
        pkl.clone(),
        8,
        -1
    ))));
    assert!(is_value(IntSparseIntVect::from_pkl(&with_i32(
        pkl.clone(),
        16,
        50
    ))));
    assert!(is_value(IntSparseIntVect::from_pkl(&with_i32(
        pkl.clone(),
        16,
        -1
    ))));
    let mut long = ULongSparseIntVect::new(10).unwrap();
    long.set_val(9, 1).unwrap();
    let mut pkl = long.to_binary();
    let entry = pkl.len() - 12;
    pkl[entry..entry + 8].copy_from_slice(&10u64.to_le_bytes());
    assert!(is_value(ULongSparseIntVect::from_pkl(&pkl)));
}

#[test]
fn sparse_int_vect_sums() {
    // RDKit sums values, and their absolute values, in int.
    let mut siv = IntSparseIntVect::new(10).unwrap();
    siv.set_val(0, i32::MAX).unwrap();
    assert_eq!(siv.get_total_val().unwrap(), i32::MAX);
    siv.set_val(1, 1).unwrap();
    assert!(is_value(siv.get_total_val()));
    assert!(is_value(dice_similarity_sparse_int_vect(&siv, &siv)));
    assert!(is_value(tanimoto_similarity_sparse_int_vect(&siv, &siv)));
    assert!(is_value(tversky_similarity_sparse_int_vect(
        &siv, &siv, 0.5, 0.5
    )));
    let mut min = IntSparseIntVect::new(10).unwrap();
    min.set_val(0, i32::MIN).unwrap();
    assert!(is_value(min.get_total_val_with_use_abs(true)));
    assert!(is_value(dice_similarity_sparse_int_vect(&min, &min)));
}

#[test]
fn formal_charges_wrap_to_8_bits() {
    // RDKit stores formal charges in 8 bits, so no two atoms overflow a sum
    // in int (U9); formal_charge_sums rests on that.
    let mut atom = Atom::from_num(6);
    atom.set_formal_charge(i32::MAX);
    assert_eq!(atom.get_formal_charge(), -1);
    let mut mol = RWMol::default();
    mol.add_atom(&atom);
    mol.add_atom(&atom);
    assert_eq!(get_formal_charge(&mol).unwrap(), -2);
    atom.set_formal_charge(128);
    assert_eq!(atom.get_formal_charge(), -128);
}

#[test]
fn formal_charge_sums() {
    // RDKit sums formal charges in int (U9). 2^24 atoms of charge -128 sum to
    // i32::MIN, which fits but whose negation does not.
    let mut atom = Atom::from_num(6);
    atom.set_formal_charge(-128);
    let mut mol = RWMol::default();
    for _ in 0..1 << 24 {
        mol.add_atom(&atom);
    }
    assert_eq!(get_formal_charge(&mol).unwrap(), i32::MIN);
    assert!(is_value(calc_mol_formula(&mol)));
    for func in [
        HashFunction::MolFormula,
        HashFunction::Mesomer,
        HashFunction::RedoxPair,
        HashFunction::NetCharge,
        HashFunction::HetAtomTautomer,
        HashFunction::HetAtomTautomerv2,
        HashFunction::HetAtomProtomer,
        HashFunction::HetAtomProtomerv2,
    ] {
        assert!(is_value(mol_hash(&mol, func)), "{func:?}");
    }
    assert!(is_value(ResonanceMolSupplier::new(&mol).map(|_| ())));
    assert!(is_value(mol_to_cml_block(&mol)));
    // CoordGen sums each fragment's charges in int.
    assert!(is_value(add_coords(&mut mol)));
    mol.add_atom(&atom);
    assert!(is_value(get_formal_charge(&mol)));
}

#[test]
fn hydrogen_count_sums() {
    // The tautomer and protomer hashes also sum hydrogen counts in int, at
    // most 255 per atom, and subtract the charge from them (U9). 255 times
    // this many atoms passes i32::MAX, with no charges.
    let atom = Atom::from_num(6);
    let mut mol = RWMol::default();
    for _ in 0..i32::MAX as u32 / 255 + 1 {
        mol.add_atom(&atom);
    }
    for func in [
        HashFunction::HetAtomTautomer,
        HashFunction::HetAtomTautomerv2,
        HashFunction::HetAtomProtomer,
        HashFunction::HetAtomProtomerv2,
    ] {
        assert!(is_value(mol_hash(&mol, func)), "{func:?}");
    }
    // The other hashes reach RDKit, which needs the property cache.
    assert!(matches!(
        mol_hash(&mol, HashFunction::NetCharge),
        Err(Error::Runtime { .. })
    ));
}

fn grid_pickle() -> (UniformGrid3D, Vec<u8>) {
    let mut grid = UniformGrid3D::new(2.0, 2.0, 2.0).unwrap();
    grid.set_val(3, 2).unwrap();
    let pkl = grid.to_binary();
    (grid, pkl)
}

#[test]
fn uniform_grid3d_pickles() {
    let (grid, pkl) = grid_pickle();
    let copy = UniformGrid3D::from_pkl(&pkl).unwrap();
    assert_eq!(copy.get_val(3).unwrap(), 2);
    assert!(copy.compare_params(&grid));
    // Header: version, nx, ny, nz, spacing, offset, storage size; then the
    // storage, a DiscreteValueVect pickle.
    for bad in [
        with_i32(pkl.clone(), 0, i32::MIN),
        with_u32(pkl.clone(), 4, 5),
        with_u32(pkl.clone(), 4, 0),
        with_u32(pkl.clone(), 4, u32::MAX),
        with_f64(pkl.clone(), 16, f64::NAN),
        with_f64(pkl.clone(), 16, 0.0),
        with_f64(pkl.clone(), 16, -0.5),
        with_f64(pkl.clone(), 16, f64::INFINITY),
        with_f64(pkl.clone(), 24, f64::NAN),
        with_f64(pkl.clone(), 40, f64::INFINITY),
        with_u32(pkl.clone(), 48, u32::MAX),
        with_u32(pkl.clone(), 52 + 8, 0),
        pkl[..pkl.len() - 1].to_vec(),
    ] {
        assert!(is_value(UniformGrid3D::from_pkl(&bad)), "{bad:?}");
    }
    assert!(UniformGrid3D::from_pkl(&pkl[..30]).is_err());
    // The storage's header is checked as a DiscreteValueVect's is.
    let mut zero_bits = pkl[..48].to_vec();
    zero_bits.extend((ZERO_BITS_PER_VALUE.len() as u32).to_le_bytes());
    zero_bits.extend(ZERO_BITS_PER_VALUE);
    assert!(is_value(UniformGrid3D::from_pkl(&zero_bits)));
}

#[test]
fn uniform_real_value_grid3d_pickles() {
    let mut grid = UniformRealValueGrid3D::new(2.0, 2.0, 2.0).unwrap();
    grid.set_val(3, 2.5).unwrap();
    let pkl = grid.to_binary();
    let copy = UniformRealValueGrid3D::from_pkl(&pkl).unwrap();
    assert!(copy.compare_grids(&grid));
    assert!(copy.clone().compare_grids(&grid));
    for bad in [
        with_i32(pkl.clone(), 0, i32::MIN),
        with_u32(pkl.clone(), 4, 5),
        with_u32(pkl.clone(), 4, 0),
        with_f64(pkl.clone(), 16, f64::NAN),
        with_f64(pkl.clone(), 16, f64::INFINITY),
        with_f64(pkl.clone(), 32, f64::NAN),
        with_i32(pkl.clone(), 52, i32::MIN),
        with_u32(pkl.clone(), 56, 7),
    ] {
        assert!(is_value(UniformRealValueGrid3D::from_pkl(&bad)), "{bad:?}");
    }
}

#[test]
fn grid_constructors() {
    // A grid with no points along an axis can't be copied.
    assert!(is_value(UniformGrid3D::new(0.1, 1.0, 1.0)));
    assert!(is_value(UniformRealValueGrid3D::new(1.0, 0.1, 1.0)));
    let bad_offset = UniformGrid3DParams {
        off_set: Some(Point3D::new(f64::NAN, 0.0, 0.0)),
        ..Default::default()
    };
    assert!(is_value(UniformGrid3D::new_with_params(
        1.0,
        1.0,
        1.0,
        &bad_offset
    )));
    let bad_spacing = UniformRealValueGrid3DParams {
        spacing: f64::INFINITY,
        off_set: None,
    };
    assert!(is_value(UniformRealValueGrid3D::new_with_params(
        1.0,
        1.0,
        1.0,
        &bad_spacing
    )));
}

fn fpb_data(name: &str) -> Vec<u8> {
    std::fs::read(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name),
    )
    .unwrap()
}

/// The offset of a chunk's payload, and its size.
fn fpb_chunk(data: &[u8], tag: &[u8; 4]) -> (usize, usize) {
    let mut pos = 8;
    loop {
        let size = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap()) as usize;
        if &data[pos + 8..pos + 12] == tag {
            return (pos + 12, size);
        }
        pos += 12 + size;
    }
}

fn open_fpb(data: &[u8], name: &str) -> rdkit_cxx::Result<FPBReader> {
    let dir = std::env::temp_dir().join(format!("rdkit-cxx-fpb-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name);
    std::fs::write(&path, data).unwrap();
    let mut reader = FPBReader::new_with_lazy(&path, true)?;
    reader.init()?;
    std::fs::remove_file(&path).unwrap();
    Ok(reader)
}

#[test]
fn fpb_files() {
    let good = fpb_data("zim.head100.fpb");
    assert_eq!(open_fpb(&good, "good.fpb").unwrap().len().unwrap(), 100);
    let (arena, arena_size) = fpb_chunk(&good, b"AREN");
    let (ids, ids_size) = fpb_chunk(&good, b"FPID");
    let fend = fpb_chunk(&good, b"FEND").0 - 12;
    let set = |offset: usize, bytes: &[u8]| {
        let mut data = good.clone();
        data[offset..offset + bytes.len()].copy_from_slice(bytes);
        data
    };
    let table = ids + ids_size - 101 * 4;
    let cases = [
        ("truncated", good[..arena + 100].to_vec()),
        ("no FEND", good[..fend].to_vec()),
        (
            "chunk past the end",
            set(arena - 12, &u64::MAX.to_le_bytes()),
        ),
        ("zero storage", set(arena + 4, &0u32.to_le_bytes())),
        (
            "fingerprints past storage",
            set(arena, &257u32.to_le_bytes()),
        ),
        ("ragged arena", set(arena + 8, &[2])),
        ("huge fingerprints", set(arena, &u32::MAX.to_le_bytes())),
        ("8-byte ids", set(ids + 4, &1u32.to_le_bytes())),
        ("id table past its chunk", set(ids, &u32::MAX.to_le_bytes())),
        ("ids out of order", set(table + 8, &0u32.to_le_bytes())),
        (
            "id past its chunk",
            set(table + 400, &(ids_size as u32 + 1).to_le_bytes()),
        ),
    ];
    assert_eq!(arena_size, 25610);
    for (what, data) in cases {
        assert!(is_value(open_fpb(&data, "bad.fpb")), "{what}");
    }
    assert!(matches!(
        open_fpb(b"FPB2", "magic.fpb"),
        Err(Error::Os { .. })
    ));
}

#[test]
fn fpb_screens() {
    // RDKit's popcount screen indexes past its table when the threshold is
    // small enough, for any file with a popcount table; the wrapper searches
    // without the screen then.
    let reader = open_fpb(&fpb_data("zim.head100.fpb"), "screen.fpb").unwrap();
    let query = reader.get_bytes(0).unwrap();
    for threshold in [0.7, 0.01, 0.005, 1e-7, 0.0] {
        let expected: Vec<u32> = (0..100)
            .filter(|&i| reader.get_tanimoto(i, &query).unwrap() >= threshold)
            .collect();
        let mut found: Vec<u32> = reader
            .get_tanimoto_neighbors_with_threshold(&query, threshold)
            .unwrap()
            .into_iter()
            .map(|(_, i)| i)
            .collect();
        found.sort();
        assert_eq!(found, expected, "{threshold}");
        let tversky = reader
            .get_tversky_neighbors_with_threshold(&query, 1.0, 1.0, threshold)
            .unwrap();
        assert_eq!(tversky.len(), expected.len(), "{threshold}");
    }
    // Weights that make RDKit's bounds infinite or NaN.
    for (ca, cb) in [(0.0, 0.0), (1.0, 0.0), (-1.0, 2.0)] {
        reader
            .get_tversky_neighbors_with_threshold(&query, ca, cb, 1.0)
            .unwrap();
    }

    let mut multi = MultiFPBReader::new();
    multi.add_reader(open_fpb(&fpb_data("zinc_random200.1.patt.fpb"), "1.fpb").unwrap());
    multi.add_reader(open_fpb(&fpb_data("zinc_random200.2.patt.fpb"), "2.fpb").unwrap());
    multi.init().unwrap();
    let query = multi.get_reader(0).unwrap().get_bytes(3).unwrap();
    let params = MultiFPBReaderNeighborsParams {
        threshold: 0.001,
        num_threads: 2,
    };
    let nbrs = multi
        .get_tanimoto_neighbors_with_params(&query, &params)
        .unwrap();
    assert_eq!(nbrs[0], (1.0, 3, 0));
    let everything = (0..2)
        .map(|r| {
            let reader = multi.get_reader(r).unwrap();
            (0..200)
                .filter(|&i| reader.get_tanimoto(i, &query).unwrap() >= 0.001)
                .count()
        })
        .sum::<usize>();
    assert_eq!(nbrs.len(), everything);
    assert!(nbrs.windows(2).all(|w| w[0].0 >= w[1].0));
}
