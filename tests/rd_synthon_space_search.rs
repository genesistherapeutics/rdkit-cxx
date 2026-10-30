// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use std::path::PathBuf;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_dist_geom::{embed_multiple_confs_with_params, etkdg_v3};
use rdkit_cxx::chem::rd_fingerprint_generator::{
    GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use rdkit_cxx::chem::rd_generalized_substruct::{
    CreateExtendedQueryMolParams, create_extended_query_mol, create_extended_query_mol_with_params,
};
use rdkit_cxx::chem::rd_rascal_mces::RascalOptions;
use rdkit_cxx::chem::rd_synthon_space_search::*;
use rdkit_cxx::chem::*;
use rdkit_cxx::geometry::rd_geometry::Point3D;

use common::{assert_close_to, temp_path};

// Expected values come from Python's rdSynthonSpaceSearch, with the files in
// tests/data.

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
    assert_eq!(p.similarity_cutoff, 0.5);
    assert_eq!(p.frag_similarity_adjuster, 0.1);
    assert_eq!(p.approx_similarity_adjuster, 0.1);
    assert_eq!(p.min_hit_heavy_atoms, 0);
    assert_eq!(p.max_hit_heavy_atoms, -1);
    assert_eq!(p.min_hit_mol_wt, 0.0);
    assert_eq!(p.max_hit_mol_wt, -1.0);
    assert_eq!(p.min_hit_chiral_atoms, 0);
    assert_eq!(p.max_hit_chiral_atoms, -1);
    assert_eq!(p.time_out, 600);
    assert_eq!(p.to_try_chunk_size, 2500000);
    assert_eq!(p.num_threads, 1);
    assert_eq!(p.possible_hits_file, None);
    assert_eq!(p.max_possible_hits_to_write, 10000000);
    assert!(!p.write_possible_hits_and_stop);
    let s = SubstructureSearchParams::default();
    assert_eq!(
        s.substruct_match_params,
        SubstructMatchParameters::default()
    );
    assert_eq!(s.params, p);
}

#[test]
fn read_text_files() {
    for (file, reactions, synthons, products) in [
        ("amide_space.txt", 1, 7, 12),
        ("triazole_space.txt", 1, 6, 8),
        ("extended_query.csv", 7, 28, 52),
        ("doebner_miller_space.txt", 1, 7, 12),
    ] {
        let mut space = SynthonSpace::default();
        space.read_text_file(data(file)).unwrap();
        assert_eq!(space.get_num_reactions(), reactions, "{file}");
        assert_eq!(space.get_num_synthons(), synthons, "{file}");
        assert_eq!(space.get_num_products(), products, "{file}");
        assert_eq!(space.get_num_synthons_with_shapes(), 0);
        assert_eq!(space.get_synthon_fingerprint_type(), "");
    }
    let empty = SynthonSpace::default();
    assert_eq!(empty.get_num_reactions(), 0);
    amide_space().summarise();
    amide_space().report_synthon_usage();
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
    assert!(result.get_best_hit().is_none());

    let mut params = SubstructureSearchParams::default();
    params.params.max_hits = 2;
    let page = space
        .substructure_search_with_params(&benzamide_query(), &params)
        .unwrap();
    assert_eq!(result_names(&page), amides(&["1-1;2-1", "1-1;2-3"]));
    params.params.hit_start = 2;
    let page = space
        .substructure_search_with_params(&benzamide_query(), &params)
        .unwrap();
    assert_eq!(result_names(&page), amides(&["1-1;2-2", "1-3;2-1"]));

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
    assert_eq!(space.get_num_synthons(), 6);
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

    let file = temp_path("synthon-ext-poss-hits.txt");
    let mut params = SubstructureSearchParams::default();
    params.params.max_hits = -1;
    params.params.possible_hits_file = Some(file.clone());
    params.params.write_possible_hits_and_stop = true;
    let written = space
        .substructure_search_extended_query_mol_with_params(&query, &params)
        .unwrap();
    assert!(written.get_hit_molecules().is_empty());
    let checked = space
        .substructure_search_extended_query_mol_lines_with_params(&query, 0, 10, &params)
        .unwrap();
    assert_eq!(checked.get_hit_molecules().len(), 10);
    std::fs::remove_file(file).unwrap();
}

#[test]
fn incremental_substructure_search() {
    let mut space = amide_space();
    let mut params = SubstructureSearchIncrementalParams::default();
    params.params.to_try_chunk_size = 2;
    let mut chunks = Vec::new();
    let mut hits = Vec::new();
    space
        .substructure_search_incremental_with_params(
            &benzamide_query(),
            |chunk| {
                chunks.push(chunk.len());
                hits.extend(chunk);
                Ok(false)
            },
            &params,
        )
        .unwrap();
    assert_eq!(chunks, [2, 2, 2]);
    assert_eq!(names(&hits), amides(&BENZAMIDES));

    // True stops the search after the first batch.
    let mut hits = Vec::new();
    space
        .substructure_search_incremental_with_params(
            &benzamide_query(),
            |chunk| {
                hits.extend(chunk);
                Ok(true)
            },
            &params,
        )
        .unwrap();
    assert_eq!(hits.len(), 2);

    // The closure's error ends the search and comes back as it was.
    let mut calls = 0;
    let err = space
        .substructure_search_incremental_with_params(
            &benzamide_query(),
            |_| {
                calls += 1;
                Err(Error::Value {
                    message: "stop here".into(),
                })
            },
            &params,
        )
        .unwrap_err();
    assert_eq!(calls, 1);
    assert!(matches!(err, Error::Value { ref message } if message == "stop here"));

    // A panic is caught at the boundary.
    let err = space
        .substructure_search_incremental(&benzamide_query(), |_| panic!("in the callback"))
        .unwrap_err();
    assert!(matches!(err, Error::Runtime { ref message } if message.contains("in the callback")));

    // RDKit's Ctrl-C handler is installed while the closure runs, so it
    // cannot match substructures.
    let pattern = mol_from_smarts("C=O").unwrap();
    let err = space
        .substructure_search_incremental(&benzamide_query(), |chunk| {
            chunk[0].has_substruct_match(&pattern)?;
            Ok(false)
        })
        .unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");
    // The lock is free again.
    assert!(benzamide_query().has_substruct_match(&pattern).unwrap());
}

#[test]
fn possible_hits_files() {
    let mut space = amide_space();
    let file = temp_path("synthon-poss-hits.txt");
    let mut params = SubstructureSearchParams::default();
    params.params.possible_hits_file = Some(file.clone());
    params.params.write_possible_hits_and_stop = true;
    let written = space
        .substructure_search_with_params(&benzamide_query(), &params)
        .unwrap();
    assert!(written.get_hit_molecules().is_empty());
    assert_eq!(written.get_max_num_results(), 6);
    let text = std::fs::read_to_string(&file).unwrap();
    assert_eq!(text.lines().count(), 6);
    assert!(text.starts_with("c1ccccc1C(=O)[1*].[1*]N1CCCC1 1-1;2-1;amide-1\n"));

    let all = space
        .substructure_search_lines_with_params(&benzamide_query(), 0, 1000, &params)
        .unwrap();
    assert_eq!(result_names(&all), amides(&BENZAMIDES));
    assert_eq!(all.get_max_num_results(), 6);
    let some = space
        .substructure_search_lines_with_params(&benzamide_query(), 1, 3, &params)
        .unwrap();
    assert_eq!(result_names(&some), amides(&["1-1;2-2", "1-1;2-3"]));
    assert_eq!(some.get_max_num_results(), 2);
    let none = space
        .substructure_search_lines_with_params(&benzamide_query(), 100, 200, &params)
        .unwrap();
    assert!(none.get_hit_molecules().is_empty());
    assert_eq!(none.get_max_num_results(), 0);
    // RDKit's precondition.
    assert!(matches!(
        space.substructure_search_lines_with_params(&benzamide_query(), 3, 3, &params),
        Err(Error::Runtime { .. })
    ));
    // With no file there is nothing to read.
    let empty = space
        .substructure_search_lines(&benzamide_query(), 0, 10)
        .unwrap();
    assert!(empty.get_hit_molecules().is_empty());

    // RDKit dereferences the null molecule of a line that does not parse (a
    // segfault in Python).
    std::fs::write(&file, "C1CC( 1-1;2-1;amide-1\n").unwrap();
    assert!(matches!(
        space.substructure_search_lines_with_params(&benzamide_query(), 0, 10, &params),
        Err(Error::Value { .. })
    ));
    // The bad line is outside these lines.
    std::fs::write(
        &file,
        "c1ccccc1C(=O)[1*].[1*]N1CCCC1 1-1;2-1;amide-1\nC1CC( x;amide-1\n",
    )
    .unwrap();
    let first = space
        .substructure_search_lines_with_params(&benzamide_query(), 0, 1, &params)
        .unwrap();
    assert_eq!(result_names(&first), amides(&["1-1;2-1"]));
    std::fs::remove_file(&file).unwrap();
    // A missing file gives no results, as RDKit's does.
    let missing = space
        .substructure_search_lines_with_params(&benzamide_query(), 0, 10, &params)
        .unwrap();
    assert!(missing.get_hit_molecules().is_empty());

    // RDKit ignores a file it cannot open for writing.
    params.params.possible_hits_file = Some(std::env::temp_dir());
    assert!(matches!(
        space.substructure_search_with_params(&benzamide_query(), &params),
        Err(Error::Os { .. })
    ));
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
    let best = result.get_best_hit().unwrap();
    assert_eq!(best.get_prop("_Name").unwrap(), "1-1;2-1;amide-1");
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

    let mut hits = Vec::new();
    space
        .fingerprint_search_incremental(&query, &generator, |chunk| {
            hits.extend(chunk);
            Ok(false)
        })
        .unwrap();
    assert_eq!(names(&hits), amides(&["1-1;2-1"]));

    // RDKit converts a negative fragment threshold times the bits to
    // unsigned int.
    let negative = SynthonSpaceSearchParams {
        similarity_cutoff: 0.05,
        ..Default::default()
    };
    assert!(matches!(
        space.fingerprint_search_with_params(&query, &generator, &negative),
        Err(Error::Value { .. })
    ));
    let nan = SynthonSpaceSearchParams {
        frag_similarity_adjuster: f64::NAN,
        ..Default::default()
    };
    assert!(matches!(
        space.fingerprint_search_with_params(&query, &generator, &nan),
        Err(Error::Value { .. })
    ));
    // A threshold of 0 is RDKit's own case.
    let zero = SynthonSpaceSearchParams {
        similarity_cutoff: 0.1,
        ..Default::default()
    };
    space
        .fingerprint_search_with_params(&query, &generator, &zero)
        .unwrap();
}

#[test]
fn fingerprint_possible_hits() {
    let mut space = amide_space();
    let generator = morgan();
    let query = mol_from_smiles("c1ccccc1C(=O)N1CCCC1").unwrap();
    let file = temp_path("synthon-fp-poss-hits.txt");
    let params = SynthonSpaceSearchParams {
        similarity_cutoff: 0.3,
        possible_hits_file: Some(file.clone()),
        write_possible_hits_and_stop: true,
        ..Default::default()
    };
    let written = space
        .fingerprint_search_with_params(&query, &generator, &params)
        .unwrap();
    assert!(written.get_hit_molecules().is_empty());
    let checked = space
        .fingerprint_search_lines(&query, &generator, &params, 0, 100)
        .unwrap();
    assert_eq!(
        result_names(&checked),
        amides(&["1-1;2-1", "1-1;2-3", "1-3;2-1"])
    );
    std::fs::remove_file(file).unwrap();
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

    let mut hits = Vec::new();
    space
        .rascal_search_incremental_with_params(
            &query,
            None,
            |chunk| {
                // The RASCAL search releases the lock for its closure.
                assert!(chunk[0].has_substruct_match(&benzamide_query())?);
                hits.extend(chunk);
                Ok(false)
            },
            &params,
        )
        .unwrap();
    assert_eq!(
        names(&hits),
        amides(&["1-1;2-1", "1-1;2-2", "1-1;2-3", "1-3;2-1"])
    );

    let file = temp_path("synthon-rascal-poss-hits.txt");
    let write = SynthonSpaceSearchParams {
        possible_hits_file: Some(file.clone()),
        write_possible_hits_and_stop: true,
        ..Default::default()
    };
    space
        .rascal_search_with_params(&query, &options, &write)
        .unwrap();
    let checked = space
        .rascal_search_lines(&query, &options, &write, 0, 3)
        .unwrap();
    assert_eq!(
        result_names(&checked),
        amides(&["1-1;2-1", "1-1;2-2", "1-1;2-3"])
    );
    // A reversed range is RDKit's length_error.
    assert!(matches!(
        space.rascal_search_lines(&query, &options, &write, 3, 1),
        Err(Error::Runtime { .. })
    ));
    std::fs::remove_file(file).unwrap();

    // RDKit divides by the fragment threshold, 0 here, and converts the
    // quotient to unsigned int.
    let zero = SynthonSpaceSearchParams {
        frag_similarity_adjuster: options.similarity_threshold,
        ..Default::default()
    };
    assert!(matches!(
        space.rascal_search_with_params(&query, &options, &zero),
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
    assert_eq!(read.get_num_synthons(), 7);
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
            ..Default::default()
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
    assert_eq!(read.get_num_synthons(), 7);

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
    assert_eq!(space.get_num_synthons(), 7);
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
        let err = space.substructure_search_lines(query, 0, 10).unwrap_err();
        assert!(matches!(err, Error::Value { .. }), "{err:?}");
        let err = space
            .substructure_search_incremental(query, |_| Ok(false))
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

fn quick_shapes() -> ShapeBuildParams {
    ShapeBuildParams {
        num_confs: 5,
        ..Default::default()
    }
}

/// Conformers from the crate's embedder, as Python's test makes them.
fn embedded(smiles: &str, num_confs: u32) -> Result<Option<Mol>, Error> {
    let mut mol = add_hs(&mol_from_smiles(smiles)?)?;
    let mut params = etkdg_v3();
    params.random_seed = 0xdac;
    embed_multiple_confs_with_params(&mut mol, num_confs, &mut params)?;
    Ok(Some(remove_hs(&mol)?))
}

#[test]
fn shape_build_params_defaults() {
    let p = ShapeBuildParams::default();
    assert_eq!(p.num_confs, 100);
    assert_eq!(p.rms_threshold, 0.5);
    assert_eq!(p.shape_sim_threshold, 0.95);
    assert_eq!(p.num_threads, 1);
    assert_eq!(p.random_seed, 0xdac);
    assert!(p.stereo_enum_opts.try_embedding);
    assert!(p.stereo_enum_opts.only_unassigned);
    assert!(!p.stereo_enum_opts.only_stereo_groups);
    assert!(p.stereo_enum_opts.unique);
    assert_eq!(p.stereo_enum_opts.max_isomers, 0);
    assert_eq!(p.stereo_enum_opts.random_seed, 0xdac);
    assert_eq!(p.max_synthon_atoms, 0);
    assert_eq!(p.max_embed_attempts, 10);
    assert_eq!(p.time_out, 60);
    assert_eq!(p.interim_file, None);
    assert_eq!(p.interim_writes, 1000);
}

#[test]
fn build_synthon_shapes() {
    let mut space = amide_space();
    space
        .build_synthon_shapes_with_py_params(&quick_shapes())
        .unwrap();
    assert_eq!(space.get_num_synthons_with_shapes(), 7);
    // Built already with this many conformers.
    space
        .build_synthon_shapes_with_py_params(&quick_shapes())
        .unwrap();
    assert_eq!(space.get_num_synthons_with_shapes(), 7);

    let file = temp_path("synthon-shapes.spc");
    space.write_db_file(&file).unwrap();
    let mut read = SynthonSpace::default();
    // SAFETY: write_db_file wrote it, and nothing changes it.
    unsafe { read.read_db_file(&file) }.unwrap();
    assert_eq!(read.get_num_synthons_with_shapes(), 7);
    let converted = temp_path("synthon-converted-shapes.spc");
    convert_text_to_db_file_with_params(
        data("amide_space.txt"),
        &converted,
        &ConvertTextToDBFileParams {
            py_shape_params: Some(quick_shapes()),
            ..Default::default()
        },
    )
    .unwrap();
    // SAFETY: as above.
    unsafe { read.read_db_file(&converted) }.unwrap();
    assert_eq!(read.get_num_synthons_with_shapes(), 7);
    std::fs::remove_file(&file).unwrap();
    std::fs::remove_file(&converted).unwrap();

    // The interim file holds the space as the build goes.
    let interim = temp_path("synthon-interim.spc");
    let params = ShapeBuildParams {
        interim_file: Some(interim.clone()),
        interim_writes: 3,
        ..quick_shapes()
    };
    let mut space = amide_space();
    space.build_synthon_shapes_with_py_params(&params).unwrap();
    // SAFETY: as above.
    unsafe { read.read_db_file(&interim) }.unwrap();
    assert_eq!(read.get_num_synthons(), 7);
    assert!(read.get_num_synthons_with_shapes() > 0);
    std::fs::remove_file(&interim).unwrap();

    // pruneShapes hands a NaN to the leader picker.
    let nan = ShapeBuildParams {
        shape_sim_threshold: f64::NAN,
        ..quick_shapes()
    };
    assert!(matches!(
        amide_space().build_synthon_shapes_with_py_params(&nan),
        Err(Error::Value { .. })
    ));

    // The embedder casts a sum over the atoms of a product with its
    // hydrogens to int.
    let iron = format!("[1*]{}", "[FeH90]".repeat(48));
    let path = write_text(
        "synthon-iron.txt",
        &[HEADER, &format!("{iron}\ta\t0\tr"), "[1*]C\tb\t1\tr"],
    );
    let mut heavy = SynthonSpace::default();
    heavy.read_text_file(&path).unwrap();
    std::fs::remove_file(path).unwrap();
    let err = heavy.build_synthon_shapes().unwrap_err();
    assert!(
        matches!(err, Error::Value { ref message } if message.contains("4295")),
        "{err:?}"
    );

    // A synthon without a connector would give a shape without volume.
    let path = write_text("synthon-loose.txt", &[HEADER, "CCO\ta\t0\tr"]);
    let mut loose = SynthonSpace::default();
    loose.read_text_file(&path).unwrap();
    std::fs::remove_file(path).unwrap();
    assert!(matches!(
        loose.build_synthon_shapes(),
        Err(Error::Value { .. })
    ));
}

#[test]
fn user_conformer_generator() {
    let mut params = quick_shapes();
    let calls = std::sync::Arc::new(std::sync::atomic::AtomicU32::new(0));
    let counted = calls.clone();
    params.set_user_conformer_generator(move |smiles, num_confs| {
        counted.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        // The build's locks are released: the closure embeds and matches.
        let mol = embedded(smiles, num_confs)?;
        let pattern = mol_from_smarts("[#6]")?;
        assert!(mol.as_ref().unwrap().has_substruct_match(&pattern)?);
        Ok(mol)
    });
    let mut space = amide_space();
    space.build_synthon_shapes_with_py_params(&params).unwrap();
    assert_eq!(space.get_num_synthons_with_shapes(), 7);
    assert!(calls.load(std::sync::atomic::Ordering::Relaxed) >= 7);
    assert_eq!(space.get_num_reactions(), 1);
    assert_eq!(space.get_num_products(), 12);

    // None drops the isomer (Python's null is dereferenced).
    params.set_user_conformer_generator(|_, _| Ok(None));
    let mut space = amide_space();
    space.build_synthon_shapes_with_py_params(&params).unwrap();
    assert_eq!(space.get_num_synthons_with_shapes(), 0);

    // An error, or a panic, ends the build.
    params.set_user_conformer_generator(|_, _| {
        Err(Error::Key {
            key: "no conformers".into(),
        })
    });
    let err = amide_space()
        .build_synthon_shapes_with_py_params(&params)
        .unwrap_err();
    assert!(matches!(err, Error::Key { ref key } if key == "no conformers"));
    params.set_user_conformer_generator(|_, _| panic!("no conformers"));
    let err = amide_space()
        .build_synthon_shapes_with_py_params(&params)
        .unwrap_err();
    assert!(matches!(err, Error::Runtime { .. }), "{err:?}");

    // A molecule RDKit removes bonds from.
    params.set_user_conformer_generator(|smiles, num_confs| {
        let mut mol = embedded(smiles, num_confs)?.unwrap();
        add_parents(&mut mol, [2, 1]);
        Ok(Some(mol))
    });
    let err = amide_space()
        .build_synthon_shapes_with_py_params(&params)
        .unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");

    // Coordinates the shapes' Gaussians cannot take.
    params.set_user_conformer_generator(|smiles, num_confs| {
        let mut mol = embedded(smiles, num_confs)?.unwrap();
        mol.get_conformer_mut()?
            .set_atom_position(0, Point3D::new(f64::NAN, 0.0, 0.0))?;
        Ok(Some(mol))
    });
    let err = amide_space()
        .build_synthon_shapes_with_py_params(&params)
        .unwrap_err();
    assert!(matches!(err, Error::Value { .. }), "{err:?}");
}
