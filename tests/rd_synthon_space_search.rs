// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use std::path::PathBuf;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_fingerprint_generator::{
    GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::chem::rd_generalized_substruct::{
    CreateExtendedQueryMolParams, create_extended_query_mol, create_extended_query_mol_with_params,
};
use rdkit_cxx::chem::rd_rascal_mces::RascalOptions;
use rdkit_cxx::chem::rd_synthon_space_search::*;
use rdkit_cxx::chem::*;

use common::{assert_close_to, temp_path};

// Expected values come from RDKit 2025.03.6's Python rdSynthonSpaceSearch,
// with the files in tests/data.

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn amide_space() -> SynthonSpace {
    let mut space = SynthonSpace::default();
    space.read_text_file(data("amide_space.txt")).unwrap();
    space
}

fn names(mols: &[Mol]) -> Vec<String> {
    let mut names: Vec<String> = mols.iter().map(|m| m.get_prop("_Name").unwrap()).collect();
    names.sort();
    names
}

fn result_names(result: &SubstructureResult) -> Vec<String> {
    names(&result.get_hit_molecules())
}

fn amides(ids: &[&str]) -> Vec<String> {
    ids.iter().map(|id| format!("{id};amide-1")).collect()
}

const BENZAMIDES: [&str; 6] = [
    "1-1;2-1", "1-1;2-2", "1-1;2-3", "1-3;2-1", "1-3;2-2", "1-3;2-3",
];

fn benzamide_query() -> Mol {
    mol_from_smarts("c1ccccc1C(=O)N").unwrap()
}

fn morgan() -> rdkit_cxx::chem::rd_fingerprint_generator::FingerprintGenerator64 {
    get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 2,
        fp_size: 2048,
        ..Default::default()
    })
    .unwrap()
}

fn write_text(name: &str, lines: &[&str]) -> PathBuf {
    let path = temp_path(name);
    std::fs::write(&path, lines.join("\n") + "\n").unwrap();
    path
}

const HEADER: &str = "SMILES\tsynton_id\tsynton#\treaction_id";

#[test]
fn search_params_defaults() {
    let p = SynthonSpaceSearchParams::default();
    assert_eq!(p.max_hits, 1000);
    assert_eq!(p.max_num_frags, 100000);
    assert_eq!(p.hit_start, 0);
    assert!(!p.random_sample);
    assert_eq!(p.random_seed, -1);
    assert!(p.build_hits);
    assert_eq!(p.num_random_sweeps, 10);
    assert_eq!(p.similarity_cutoff, 0.5);
    assert_eq!(p.frag_similarity_adjuster, 0.1);
    assert_eq!(p.approx_similarity_adjuster, 0.1);
    assert_eq!(p.min_hit_heavy_atoms, 0);
    assert_eq!(p.max_hit_heavy_atoms, -1);
    assert_eq!(p.min_hit_mol_wt, 0.0);
    assert_eq!(p.max_hit_mol_wt, 0.0);
    assert_eq!(p.min_hit_chiral_atoms, 0);
    assert_eq!(p.max_hit_chiral_atoms, -1);
    assert_eq!(p.time_out, 600);
    assert_eq!(p.num_threads, 1);
    let s = SubstructureSearchParams::default();
    assert_eq!(
        s.substruct_match_params,
        SubstructMatchParameters::default()
    );
    assert_eq!(s.params, p);
}

#[test]
fn read_text_files() {
    for (file, reactions, products) in [
        ("amide_space.txt", 1, 12),
        ("triazole_space.txt", 1, 8),
        ("extended_query.csv", 7, 52),
        ("doebner_miller_space.txt", 1, 12),
    ] {
        let mut space = SynthonSpace::default();
        space.read_text_file(data(file)).unwrap();
        assert_eq!(space.get_num_reactions(), reactions, "{file}");
        assert_eq!(space.get_num_products(), products, "{file}");
        assert_eq!(space.get_synthon_fingerprint_type(), "");
    }
    let empty = SynthonSpace::default();
    assert_eq!(empty.get_num_reactions(), 0);
    amide_space().summarise();
    assert!(matches!(
        SynthonSpace::default().read_text_file(temp_path("missing.txt")),
        Err(Error::Os { .. })
    ));
}

#[test]
fn substructure_search() {
    let mut space = amide_space();
    let result = space.substructure_search(&benzamide_query()).unwrap();
    assert_eq!(result_names(&result), amides(&BENZAMIDES));
    let mut smiles: Vec<String> = result
        .get_hit_molecules()
        .iter()
        .map(|m| mol_to_smiles(m).unwrap())
        .collect();
    smiles.sort();
    let mut expected = vec![
        "O=C(c1ccccc1)N1CCCC1",
        "O=C(NCC1CC1)c1ccccc1",
        "CCCCCNC(=O)c1ccccc1",
        "O=C(c1ccccc1Cl)N1CCCC1",
        "O=C(NCC1CC1)c1ccccc1Cl",
        "CCCCCNC(=O)c1ccccc1Cl",
    ];
    expected.sort();
    assert_eq!(smiles, expected);
    assert_eq!(result.get_max_num_results(), 6);
    assert!(!result.get_timed_out());
    assert!(!result.get_cancelled());

    let mut params = SubstructureSearchParams::default();
    params.params.max_hits = 2;
    let page = space
        .substructure_search_with_params(&benzamide_query(), &params)
        .unwrap();
    assert_eq!(result_names(&page), amides(&["1-1;2-1", "1-1;2-2"]));
    params.params.hit_start = 2;
    let page = space
        .substructure_search_with_params(&benzamide_query(), &params)
        .unwrap();
    assert_eq!(result_names(&page), amides(&["1-1;2-3", "1-3;2-1"]));

    let mut params = SubstructureSearchParams::default();
    params.params.build_hits = false;
    let counted = space
        .substructure_search_with_params(&benzamide_query(), &params)
        .unwrap();
    assert!(counted.get_hit_molecules().is_empty());
    assert_eq!(counted.get_max_num_results(), 6);

    let mut params = SubstructureSearchParams::default();
    params.params.max_hit_heavy_atoms = 13;
    let small = space
        .substructure_search_with_params(&benzamide_query(), &params)
        .unwrap();
    assert_eq!(result_names(&small), amides(&["1-1;2-1", "1-1;2-3"]));

    // RDKit's precondition.
    assert!(matches!(
        space.substructure_search(&Mol::default()),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn other_spaces() {
    let mut space = SynthonSpace::default();
    space
        .read_text_file(data("doebner_miller_space.txt"))
        .unwrap();
    let result = space
        .substructure_search(&mol_from_smarts("c1ccc2ncccc2c1").unwrap())
        .unwrap();
    let expected: Vec<String> = ["1-1", "1-2", "1-3", "1-4"]
        .iter()
        .flat_map(|a| ["2-1", "2-2", "2-3"].map(|b| format!("{a};{b};doebner-miller-quinoline")))
        .collect();
    assert_eq!(result_names(&result), expected);

    space.read_text_file(data("triazole_space.txt")).unwrap();
    assert_eq!(space.get_num_products(), 8);
    let result = space
        .substructure_search(&mol_from_smarts("OCc1ncnn1").unwrap())
        .unwrap();
    assert_eq!(result_names(&result).len(), 8);
    assert!(
        result_names(&result)
            .iter()
            .all(|n| n.ends_with(";triazole-1"))
    );
}

#[test]
fn extended_query_search() {
    let mut space = SynthonSpace::default();
    space.read_text_file(data("extended_query.csv")).unwrap();
    let query =
        create_extended_query_mol(&mol_from_smarts("[#6]-*.c1nc2cccnc2n1 |m:1:3.10|").unwrap())
            .unwrap();
    let result = space
        .substructure_search_extended_query_mol(&query)
        .unwrap();
    let expected = [
        "1-1;2-1;react_1",
        "1-1;2-2;react_1",
        "1-2;2-1;react_1",
        "1-2;2-2;react_1",
        "1-3;2-1;react_1",
        "1-3;2-2;react_1",
        "3-1;2-1;react_2",
        "3-1;2-2;react_2",
        "3-2;2-1;react_2",
        "3-2;2-2;react_2",
        "4-1;5-1;react_3",
        "4-2;5-1;react_3",
    ];
    assert_eq!(result_names(&result), expected);
}

#[test]
fn search_params_checks() {
    let mut space = amide_space();
    let search = |space: &mut SynthonSpace, f: &dyn Fn(&mut SynthonSpaceSearchParams)| {
        let mut params = SubstructureSearchParams::default();
        f(&mut params.params);
        space.substructure_search_with_params(&benzamide_query(), &params)
    };
    // RDKit erases from before the results' beginning for these (segfaults
    // in Python).
    assert!(matches!(
        search(&mut space, &|p| p.max_hits = -5),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        search(&mut space, &|p| p.hit_start = -5),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        search(&mut space, &|p| {
            p.max_hits = i64::MAX;
            p.hit_start = 1;
        }),
        Err(Error::Value { .. })
    ));
    assert!(matches!(
        search(&mut space, &|p| p.time_out = u64::MAX),
        Err(Error::Value { .. })
    ));
    // -1 is no limit, and maxHits plus hitStart at i64::MAX is fine.
    let all = search(&mut space, &|p| {
        p.max_hits = -1;
        p.hit_start = i64::MAX;
    })
    .unwrap();
    assert!(all.get_hit_molecules().is_empty());
    // RDKit's own check.
    assert!(matches!(
        search(&mut space, &|p| {
            p.random_sample = true;
            p.max_hits = -1;
        }),
        Err(Error::Runtime { .. })
    ));
    let sample = search(&mut space, &|p| {
        p.random_sample = true;
        p.random_seed = 42;
        p.max_hits = 3;
    })
    .unwrap();
    assert_eq!(sample.get_hit_molecules().len(), 3);
    // num_threads is not read.
    let threads = search(&mut space, &|p| p.num_threads = 8).unwrap();
    assert_eq!(result_names(&threads), amides(&BENZAMIDES));
}

#[test]
fn fingerprint_search() {
    let mut space = amide_space();
    let generator = morgan();
    let query = mol_from_smiles("c1ccccc1C(=O)N1CCCC1").unwrap();
    let params = SynthonSpaceSearchParams::default();
    let result = space
        .fingerprint_search_with_params(&query, &generator, &params)
        .unwrap();
    let hits = result.get_hit_molecules();
    assert_eq!(names(&hits), amides(&["1-1;2-1"]));
    assert_close_to(hits[0].get_double_prop("Similarity").unwrap(), 1.0, 1e-6);
    assert_eq!(result.get_max_num_results(), 5);
    assert_eq!(
        space.get_synthon_fingerprint_type(),
        generator.get_info_string().unwrap()
    );

    let low = SynthonSpaceSearchParams {
        similarity_cutoff: 0.3,
        ..Default::default()
    };
    let result = space
        .fingerprint_search_with_params(&query, &generator, &low)
        .unwrap();
    assert_eq!(
        result_names(&result),
        amides(&["1-1;2-1", "1-1;2-3", "1-3;2-1"])
    );

    // RDKit compares the fragment threshold as a double: below 0 it passes
    // every fragment.
    let negative = SynthonSpaceSearchParams {
        similarity_cutoff: 0.05,
        ..Default::default()
    };
    let result = space
        .fingerprint_search_with_params(&query, &generator, &negative)
        .unwrap();
    assert_eq!(result.get_hit_molecules().len(), 12);
    let zero = SynthonSpaceSearchParams {
        similarity_cutoff: 0.1,
        ..Default::default()
    };
    let result = space
        .fingerprint_search_with_params(&query, &generator, &zero)
        .unwrap();
    assert_eq!(result.get_hit_molecules().len(), 10);

    // Another generator's fingerprints are RDKit's error; building them
    // replaces the space's.
    let other = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius: 1,
        fp_size: 1024,
        ..Default::default()
    })
    .unwrap();
    assert!(matches!(
        space.fingerprint_search(&query, &other),
        Err(Error::Runtime { .. })
    ));
    space.build_synthon_fingerprints(&other).unwrap();
    assert_eq!(
        space.get_synthon_fingerprint_type(),
        other.get_info_string().unwrap()
    );
    let result = space.fingerprint_search(&query, &other).unwrap();
    assert_eq!(result_names(&result), amides(&["1-1;2-1", "1-3;2-1"]));
}

#[test]
fn rascal_search() {
    let mut space = amide_space();
    let query = mol_from_smiles("c1ccccc1C(=O)N1CCCC1").unwrap();
    let options = RascalOptions::default();
    let params = SynthonSpaceSearchParams {
        max_hits: 4,
        ..Default::default()
    };
    let result = space
        .rascal_search_with_params(&query, &options, &params)
        .unwrap();
    let hits = result.get_hit_molecules();
    assert_eq!(
        names(&hits),
        amides(&["1-1;2-1", "1-1;2-2", "1-1;2-3", "1-3;2-1"])
    );
    let mut sims: Vec<f64> = hits
        .iter()
        .map(|m| m.get_double_prop("Similarity").unwrap())
        .collect();
    sims.sort_by(|a, b| b.partial_cmp(a).unwrap());
    for (sim, expected) in sims.iter().zip([1.0, 0.931, 0.9273, 0.8942]) {
        assert_close_to(*sim, expected, 1e-4);
    }
    let all = space.rascal_search(&query, &options).unwrap();
    assert_eq!(result_names(&all), amides(&BENZAMIDES));

    // A fragment threshold of 0 is fine in RDKit 2025.03; RASCAL casts its
    // product with a NaN one to unsigned int.
    let zero = SynthonSpaceSearchParams {
        frag_similarity_adjuster: options.similarity_threshold,
        ..Default::default()
    };
    let result = space
        .rascal_search_with_params(&query, &options, &zero)
        .unwrap();
    assert_eq!(result_names(&result), amides(&BENZAMIDES));
    let nan = SynthonSpaceSearchParams {
        frag_similarity_adjuster: f64::NAN,
        ..Default::default()
    };
    assert!(matches!(
        space.rascal_search_with_params(&query, &options, &nan),
        Err(Error::Value { .. })
    ));
}

#[test]
fn database_files() {
    let space = amide_space();
    let file = temp_path("synthon-amide.spc");
    space.write_db_file(&file).unwrap();
    let mut read = SynthonSpace::default();
    // SAFETY: write_db_file wrote it, and nothing changes it.
    unsafe { read.read_db_file(&file) }.unwrap();
    assert_eq!(read.get_num_reactions(), 1);
    assert_eq!(read.get_num_products(), 12);
    let result = read.substructure_search(&benzamide_query()).unwrap();
    assert_eq!(result_names(&result), amides(&BENZAMIDES));

    // With fingerprints.
    let generator = morgan();
    let converted = temp_path("synthon-amide-fps.spc");
    convert_text_to_db_file_with_params(
        data("amide_space.txt"),
        &converted,
        &ConvertTextToDBFileParams {
            fp_gen: Some(&generator),
        },
    )
    .unwrap();
    let mut with_fps = SynthonSpace::default();
    // SAFETY: as above.
    unsafe { with_fps.read_db_file(&converted) }.unwrap();
    assert_eq!(
        with_fps.get_synthon_fingerprint_type(),
        generator.get_info_string().unwrap()
    );
    let query = mol_from_smiles("c1ccccc1C(=O)N1CCCC1").unwrap();
    let result = with_fps.fingerprint_search(&query, &generator).unwrap();
    assert_eq!(result_names(&result), amides(&["1-1;2-1"]));
    convert_text_to_db_file(data("amide_space.txt"), &converted).unwrap();
    // SAFETY: as above.
    unsafe { with_fps.read_db_file(&converted) }.unwrap();
    assert_eq!(with_fps.get_synthon_fingerprint_type(), "");
    std::fs::remove_file(&converted).unwrap();

    // A file cut short, its offsets' reads past the mapped file in RDKit.
    let bytes = std::fs::read(&file).unwrap();
    let short = temp_path("synthon-short.spc");
    for len in [10, 60, bytes.len() / 2] {
        std::fs::write(&short, &bytes[..len]).unwrap();
        // SAFETY: the wrapper refuses the file before RDKit reads it.
        let err = unsafe { read.read_db_file(&short) }.unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{len}: {err:?}");
    }
    // Cut inside the last record, which RDKit reads through a bounded
    // stream that throws.
    std::fs::write(&short, &bytes[..bytes.len() - 1]).unwrap();
    // SAFETY: the offsets lie within the file, and every record but the
    // last is RDKit's.
    assert!(unsafe { read.read_db_file(&short) }.is_err());
    // A later format, which RDKit would read as its own.
    let mut later = bytes.clone();
    later[8] = 3;
    std::fs::write(&short, &later).unwrap();
    // SAFETY: as above.
    assert!(matches!(
        unsafe { read.read_db_file(&short) },
        Err(Error::Value { .. })
    ));
    std::fs::remove_file(&short).unwrap();
    // The failed reads left the space as it was.
    assert_eq!(read.get_num_products(), 12);

    // An empty space writes a file RDKit cannot read: it reads the first
    // reaction's offset regardless (a segfault in Python).
    let empty = temp_path("synthon-empty.spc");
    SynthonSpace::default().write_db_file(&empty).unwrap();
    // SAFETY: as above.
    assert!(matches!(
        unsafe { read.read_db_file(&empty) },
        Err(Error::Value { .. })
    ));
    std::fs::remove_file(&empty).unwrap();

    // RDKit ignores a stream it cannot write.
    assert!(matches!(
        space.write_db_file(std::env::temp_dir()),
        Err(Error::Os { .. })
    ));
    // SAFETY: there is no file.
    assert!(unsafe { read.read_db_file(temp_path("synthon-missing.spc")) }.is_err());
    std::fs::remove_file(&file).unwrap();
}

#[test]
fn enumerated_file() {
    let space = amide_space();
    let file = temp_path("synthon-enum.smi");
    space.write_enumerated_file(&file).unwrap();
    let mut lines: Vec<String> = std::fs::read_to_string(&file)
        .unwrap()
        .lines()
        .map(String::from)
        .collect();
    lines.sort();
    let mut expected = vec![
        "O=C(c1ccccc1)N1CCCC1 1-1;2-1;amide-1",
        "CCCCCNC(=O)c1ccccc1 1-1;2-2;amide-1",
        "O=C(NCC1CC1)c1ccccc1 1-1;2-3;amide-1",
        "O=C(C1CCCCC1)N1CCCC1 1-2;2-1;amide-1",
        "CCCCCNC(=O)C1CCCCC1 1-2;2-2;amide-1",
        "O=C(NCC1CC1)C1CCCCC1 1-2;2-3;amide-1",
        "O=C(c1ccccc1Cl)N1CCCC1 1-3;2-1;amide-1",
        "CCCCCNC(=O)c1ccccc1Cl 1-3;2-2;amide-1",
        "O=C(NCC1CC1)c1ccccc1Cl 1-3;2-3;amide-1",
        "CCCCC(=O)N1CCCC1 1-4;2-1;amide-1",
        "CCCCCNC(=O)CCCC 1-4;2-2;amide-1",
        "CCCCC(=O)NCC1CC1 1-4;2-3;amide-1",
    ];
    expected.sort();
    assert_eq!(lines, expected);
    std::fs::remove_file(&file).unwrap();
    assert!(matches!(
        space.write_enumerated_file(std::env::temp_dir()),
        Err(Error::Os { .. })
    ));
}

#[test]
fn text_file_checks() {
    // A read RDKit refuses part way leaves its synthons without molecules,
    // which a search then dereferences (a segfault in Python); here the
    // space keeps what it held.
    let mut space = amide_space();
    let err = space
        .read_text_file(data("amide_space_error.txt"))
        .unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
    assert_eq!(space.get_num_products(), 12);
    let result = space.substructure_search(&benzamide_query()).unwrap();
    assert_eq!(result_names(&result).len(), 6);

    let refused = |name: &str, lines: &[&str]| {
        let path = write_text(name, &[&[HEADER], lines].concat());
        let mut space = SynthonSpace::default();
        let err = space.read_text_file(&path).unwrap_err();
        std::fs::remove_file(path).unwrap();
        err
    };
    // RDKit indexes the reaction's synthon sets by the number (a segfault in
    // Python for a negative one).
    for num in ["-3", "16", "2147483647"] {
        let line = format!("c1ccccc1C(=O)[1*]\t1-1\t{num}\tamide-1");
        let err = refused("synthon-num.txt", &[&line, "[1*]N1CCCC1\t2-1\t1\tamide-1"]);
        assert!(matches!(err, Error::Value { .. }), "{num}: {err:?}");
    }
    // Five synthon sets overrun RDKit's count of reactions by their sets.
    let err = refused(
        "synthon-five.txt",
        &[
            "[1*]C([2*])([3*])[4*]\ta\t0\tr",
            "[1*]N\tb\t1\tr",
            "[2*]O\tc\t2\tr",
            "[3*]S\td\t3\tr",
            "[4*]F\te\t4\tr",
        ],
    );
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    // A synthon without atoms (a segfault in Python), or with only dummies.
    for smiles in ["", "[1*][2*]", "[U][Np]"] {
        let line = format!("{smiles}\tb\t0\tr");
        let err = refused("synthon-dummies.txt", &[&line, "[1*]N\tc\t1\tr"]);
        assert!(matches!(err, Error::Value { .. }), "{smiles:?}: {err:?}");
    }
    // A run of separators RDKit's std::regex matches recursively.
    let spaced = format!("[1*]N{}c\t1\tr", " ".repeat(5000));
    let err = refused("synthon-spaces.txt", &["C[1*]\tb\t0\tr", &spaced]);
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    // A synthon SMILES that does not parse is RDKit's own error.
    let err = refused("synthon-smiles.txt", &["C1CC(\tb\t0\tr"]);
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    // And a line with too few fields.
    let err = refused("synthon-fields.txt", &["C[1*]\tb\t0"]);
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
    // A hydrogen synthon is fine.
    let path = write_text(
        "synthon-h.txt",
        &[HEADER, "[1*][H]\ta\t0\tr", "[1*]N\tc\t1\tr"],
    );
    let mut space = SynthonSpace::default();
    space.read_text_file(&path).unwrap();
    assert_eq!(space.get_num_products(), 1);
    std::fs::remove_file(path).unwrap();
}

/// Groups on atoms 1 and 2 whose PARENTs are `parents`, and one on atom 0,
/// which removing atom 0 drops: RDKit then follows the others' chains.
fn add_parents(mol: &mut Mol, parents: [u32; 2]) {
    for (i, (atom, parent)) in [(1, Some(parents[0])), (2, Some(parents[1])), (0, None)]
        .into_iter()
        .enumerate()
    {
        let mut group = create_mol_substance_group(mol, "DAT");
        group.add_atom_with_idx(atom).unwrap();
        group.set_unsigned_prop("index", i as u32 + 1);
        if let Some(parent) = parent {
            group.set_unsigned_prop("PARENT", parent);
        }
    }
}

// The reads and searches remove atoms and bonds from the query's fragments,
// the synthons, and their products: RDKit follows a cycle of substance
// groups' PARENTs forever (Python's searches hang) and misreads a bond's
// _MolFileBondEndPts.
#[test]
fn substance_group_parents() {
    let mut space = amide_space();
    let plain = result_names(&space.substructure_search(&benzamide_query()).unwrap());
    let mut chain = benzamide_query();
    add_parents(&mut chain, [2, 3]);
    let found = space.substructure_search(&chain).unwrap();
    assert_eq!(result_names(&found), plain);

    let mut cycle = benzamide_query();
    add_parents(&mut cycle, [2, 1]);
    let mut end_points = benzamide_query();
    end_points
        .get_bond_with_idx_mut(0)
        .unwrap()
        .set_prop("_MolFileBondEndPts", "()");
    for query in [&cycle, &end_points] {
        let err = space.substructure_search(query).unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{err:?}");
        // Enumerating or making tautomers removes atoms too.
        let err = create_extended_query_mol(query).unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{err:?}");
        let as_is = CreateExtendedQueryMolParams {
            do_enumeration: false,
            do_tautomers: false,
            ..Default::default()
        };
        let extended = create_extended_query_mol_with_params(query, &as_is).unwrap();
        let err = space
            .substructure_search_extended_query_mol(&extended)
            .unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{err:?}");
    }
    let mut cycle = mol_from_smiles("c1ccccc1C(=O)NC").unwrap();
    add_parents(&mut cycle, [2, 1]);
    let fp_gen = morgan();
    space.build_synthon_fingerprints(&fp_gen).unwrap();
    let err = space.fingerprint_search(&cycle, &fp_gen).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    let err = space
        .rascal_search(&cycle, &RascalOptions::default())
        .unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");

    // A synthon's own groups, as the reader removes atoms from it.
    let refused = |lines: &[&str]| {
        let path = write_text("synthon-parents.txt", &[&[HEADER], lines].concat());
        let mut space = SynthonSpace::default();
        let res = space.read_text_file(&path);
        std::fs::remove_file(path).unwrap();
        res
    };
    const CYCLE: &str = "|SgD:0:a:b::::,SgD:1:c:d::::,SgD:2:e:f::::,SgH:1:2,2:1|";
    let line = format!("c1ccccc1C(=O)[1*] {CYCLE}\t1-1\t1\tamide-1");
    let err = refused(&[&line, "[1*]N1CCCC1\t2-1\t2\tamide-1"]).unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    // A chain is fine.
    let line = "c1ccccc1C(=O)[1*] |SgD:0:a:b::::,SgD:1:c:d::::,SgH:1:0|\t1-1\t1\tamide-1";
    refused(&[line, "[1*]N1CCCC1\t2-1\t2\tamide-1"]).unwrap();
    // A reaction's synthons' chains are checked together, as a product
    // holds them all: the first synthon's group 2 names 3 as its PARENT, and
    // the second's group 3 names 2. In a product the second synthon's group
    // 2, which has no PARENT, hides the first's, so no product has a cycle
    // (Python reads and searches the space): the check refuses more than it
    // must. In two reactions they are fine.
    let first = "c1ccccc1C(=O)[1*] |SgD:0:a:b::::,SgD:1:c:d::::,SgD:2:e:f::::,SgH:2:1|";
    let second = "[1*]N1CCCC1 |SgD:1:a:b::::,SgD:2:c:d::::,SgD:3:e:f::::,SgH:1:2|";
    let err = refused(&[
        &format!("{first}\t1-1\t1\tamide-1"),
        &format!("{second}\t2-1\t2\tamide-1"),
    ])
    .unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
    refused(&[
        &format!("{first}\t1-1\t1\tamide-1"),
        &format!("{second}\t2-1\t2\tamide-2"),
    ])
    .unwrap();
}

#[test]
fn formatted_integers() {
    assert_eq!(formatted_integer_string(1234567), "1 234 567");
    assert_eq!(formatted_integer_string(-1234567), "-1 234 567");
    assert_eq!(formatted_integer_string(12), "12");
}

#[test]
fn space_is_send() {
    fn send<T: Send>() {}
    send::<SynthonSpace>();
    send::<SubstructureResult>();
}
