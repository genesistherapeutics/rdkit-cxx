// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

mod common;

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use rdkit_cxx::Error;
use rdkit_cxx::chem::rdfragcatalog::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdfragcatalog in RDKit 2024.03.6.

const GROUPS: &str = "// functional groups\n-C(=O)O\t*-C(=O)[O;H1,-]\tcarboxylic acids\n-O\t*-[O;H1,-]\talcohols\n\n=O\t*=O\taldehydes\n";

/// Python's pickle of `catalog()`.
const CATALOG_PICKLE: &str = "efbeadde01000000000000000000000005000000050000003120332031652d30380a330a2d43283d4f294f092a2d43283d4f295b4f3b48312c2d5d0a2d4f092a2d5b4f3b48312c2d5d0a3d4f092a3d4f0aefbeadde000000001000000001000000000000000200000001000000800106006000000002020600680000000303010b000100170400000000000000160000000006000000433c2d4f3e430100000001000000000000000100000001000000efbeadde00000000100000000100000000000000020000000100000080010600680000000303010600680000000303010b0001280217040000000000000016010000000c000000433d433c2d43283d4f294f3e0100000001000000010000000100000000000000efbeadde000000001000000001000000000000000300000002000000800106006000000002020600680000000303010600680000000303010b00010001022802170400000000000000160200000011000000433c2d43283d4f294f3e3d43433c2d4f3e0200000002000000000000000100000001000000020000000100000000000000efbeadde0000000010000000010000000000000002000000010000008001060060000000020206006000000001030b00010017040000000000000016030000000200000043430100000000000000efbeadde00000000100000000100000000000000030000000200000080010600600000000202060060000000020206006000000001030b000100010200170400000000000000160400000007000000433c2d4f3e43430200000001000000000000000100000001000000020000000200000004000000010000000200000000000000010000000400000000000000";

fn hex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn groups_file(name: &str, text: &str) -> PathBuf {
    let path = common::temp_path(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// `GROUPS` in a file, written once: tests run on several threads.
fn groups() -> &'static Path {
    static PATH: OnceLock<PathBuf> = OnceLock::new();
    PATH.get_or_init(|| groups_file("groups.txt", GROUPS))
}

fn params(l_len: u32, u_len: u32) -> FragCatParams {
    FragCatParams::new(l_len, u_len, groups()).unwrap()
}

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn catalog() -> FragCatalog {
    let mut catalog = FragCatalog::new(&params(1, 3));
    let generator = FragCatGenerator;
    assert_eq!(
        generator
            .add_frags_from_mol(&mol("OCC=CC(=O)O"), &mut catalog)
            .unwrap(),
        3
    );
    assert_eq!(
        generator
            .add_frags_from_mol(&mol("OCC=CC(=O)O"), &mut catalog)
            .unwrap(),
        0
    );
    assert_eq!(
        generator
            .add_frags_from_mol(&mol("OCCC"), &mut catalog)
            .unwrap(),
        2
    );
    catalog
}

const DESCRIPTIONS: [&str; 5] = [
    "C<-O>C",
    "C=C<-C(=O)O>",
    "C<-C(=O)O>=CC<-O>",
    "CC",
    "C<-O>CC",
];

fn check_catalog(catalog: &FragCatalog) {
    assert_eq!(catalog.get_num_entries(), 5);
    assert_eq!(catalog.get_fp_length(), 5);
    let orders = [1, 1, 2, 1, 2];
    let groups: [&[i32]; 5] = [&[1], &[0], &[1, 0], &[], &[1]];
    let down: [&[u32]; 5] = [&[2, 4], &[2], &[], &[4], &[]];
    let discrims = [
        [1427608728.0, 1.0, 2.0],
        [1427608601.0, 1.0, 2.0],
        [4216696482.0, 2.0, 3.0],
        [2749589652.0, 1.0, 2.0],
        [1635518915.0, 2.0, 3.0],
    ];
    for i in 0..5 {
        let at = i as usize;
        assert_eq!(catalog.get_entry_description(i).unwrap(), DESCRIPTIONS[at]);
        assert_eq!(catalog.get_entry_order(i).unwrap(), orders[at]);
        assert_eq!(catalog.get_entry_func_group_ids(i).unwrap(), groups[at]);
        assert_eq!(catalog.get_entry_bit_id(i).unwrap(), Some(i));
        assert_eq!(catalog.get_entry_down_ids(i).unwrap(), down[at]);
        assert_eq!(catalog.get_bit_description(i).unwrap(), DESCRIPTIONS[at]);
        assert_eq!(catalog.get_bit_order(i).unwrap(), orders[at]);
        assert_eq!(catalog.get_bit_func_group_ids(i).unwrap(), groups[at]);
        assert_eq!(catalog.get_bit_entry_id(i).unwrap(), Some(i));
        assert_eq!(catalog.get_bit_discrims(i).unwrap(), discrims[at]);
    }
}

fn on_bits(catalog: &FragCatalog, smiles: &str) -> Vec<u32> {
    let fp = FragFPGenerator
        .get_fp_for_mol(&mol(smiles), catalog)
        .unwrap();
    assert_eq!(fp.get_num_bits(), catalog.get_fp_length());
    fp.get_on_bits()
}

fn is_index<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Index { .. }))
}

fn is_value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

#[test]
fn parameters() {
    let params = params(1, 3);
    assert_eq!(params.get_type_string(), "Fragment Catalog Parameters");
    assert_eq!(params.get_lower_frag_length(), 1);
    assert_eq!(params.get_upper_frag_length(), 3);
    assert!((params.get_tolerance() - 1e-8).abs() < 1e-20);
    assert_eq!(params.get_num_func_groups(), 3);
    let group = params.get_func_group(0).unwrap();
    assert_eq!(mol_to_smarts(&group).unwrap(), "*-C(=O)[O;H1,-]");
    assert_eq!(group.get_prop("_Name").unwrap(), "-C(=O)O");
    assert_eq!(group.get_prop("_fragSMARTS").unwrap(), "*-C(=O)[O;H1,-]");
    assert_eq!(params.get_func_group(2).unwrap().get_num_atoms(), 2);
    assert!(is_index(params.get_func_group(3)));
    assert_eq!(
        params.serialize(),
        b"1 3 1e-08\n3\n-C(=O)O\t*-C(=O)[O;H1,-]\n-O\t*-[O;H1,-]\n=O\t*=O\n"
    );

    let params = FragCatParams::new_with_tol(2, 4, groups(), 0.001).unwrap();
    assert!((params.get_tolerance() - 0.001).abs() < 1e-15);
    assert_eq!(&params.serialize()[..10], b"2 4 0.001\n");
}

#[test]
fn bad_parameters() {
    let path = groups();
    let runtime = |result: rdkit_cxx::Result<FragCatParams>| {
        assert!(matches!(result, Err(Error::Runtime { .. })), "{result:?}");
    };
    runtime(FragCatParams::new(3, 1, path));
    runtime(FragCatParams::new(
        1,
        3,
        groups_file("bad-smarts.txt", "name\tC(\n"),
    ));
    runtime(FragCatParams::new(
        1,
        3,
        groups_file("no-smarts.txt", "nameonly\n"),
    ));
    assert!(matches!(
        FragCatParams::new(1, 3, common::temp_path("missing.txt")),
        Err(Error::Os { .. })
    ));
    assert!(is_value(FragCatParams::new(1, 3, "groups\0.txt")));
}

#[test]
fn catalogs() {
    let catalog = catalog();
    check_catalog(&catalog);
    assert_eq!(catalog.serialize().unwrap(), hex(CATALOG_PICKLE));
    for pickle in [hex(CATALOG_PICKLE), catalog.serialize().unwrap()] {
        // SAFETY: a pickle RDKit wrote.
        let back = unsafe { FragCatalog::from_pickle(&pickle) }.unwrap();
        check_catalog(&back);
        assert_eq!(back.get_catalog_params().unwrap().get_num_func_groups(), 3);
    }
    let copy = catalog.get_catalog_params().unwrap();
    assert_eq!(copy.get_upper_frag_length(), 3);

    assert_eq!(on_bits(&catalog, "OCC=CC(=O)O"), [0, 1, 2]);
    assert_eq!(on_bits(&catalog, "OCCC"), [0, 3, 4]);
    assert_eq!(on_bits(&catalog, "CCCC"), [3]);
    assert!(on_bits(&catalog, "c1ccccc1").is_empty());

    // Fragments outside the lengths get no bits.
    let mut catalog = FragCatalog::new(&params(2, 4));
    assert_eq!(
        FragCatGenerator
            .add_frags_from_mol(&mol("OCC=CC(=O)O"), &mut catalog)
            .unwrap(),
        3
    );
    assert_eq!(catalog.get_fp_length(), 1);
    let bits: Vec<_> = (0..3)
        .map(|i| catalog.get_entry_bit_id(i).unwrap())
        .collect();
    assert_eq!(bits, [None, None, Some(0)]);
}

#[test]
fn catalogs_without_parameters() {
    let mut catalog = FragCatalog::default();
    assert!(catalog.get_catalog_params().is_none());
    assert!(matches!(catalog.serialize(), Err(Error::Runtime { .. })));
    // Python's AddFragsFromMol dereferences the missing parameters.
    assert!(is_value(
        FragCatGenerator.add_frags_from_mol(&mol("CCO"), &mut catalog)
    ));
    assert!(is_value(
        FragFPGenerator.get_fp_for_mol(&mol("CCO"), &catalog)
    ));
    // A catalog without fragments of one bond fingerprints nothing, as in
    // Python.
    let empty = FragCatalog::new(&params(1, 3));
    assert!(matches!(
        FragFPGenerator.get_fp_for_mol(&mol("CCO"), &empty),
        Err(Error::Runtime { .. })
    ));
}

// Python checks each index one past the end, where RDKit's own check
// raises; GetEntryDownIds has none, so it reads past the catalog's graph.
#[test]
fn indices_are_checked() {
    let catalog = catalog();
    assert!(is_index(catalog.get_bit_description(5)));
    assert!(is_index(catalog.get_bit_order(5)));
    assert!(is_index(catalog.get_bit_func_group_ids(5)));
    assert!(is_index(catalog.get_bit_entry_id(5)));
    assert!(is_index(catalog.get_bit_discrims(5)));
    assert!(is_index(catalog.get_entry_bit_id(5)));
    assert!(is_index(catalog.get_entry_description(5)));
    assert!(is_index(catalog.get_entry_order(5)));
    assert!(is_index(catalog.get_entry_func_group_ids(5)));
    assert!(is_index(catalog.get_entry_down_ids(5)));
    assert!(is_index(
        FragCatalog::new(&params(1, 3)).get_entry_down_ids(0)
    ));
}

// RDKit divides each path's product of primes by the prime of each of its
// bonds' indices, with no modulo, reading past its table of 1000 primes.
#[test]
fn large_molecules_are_checked() {
    let big = mol(&"C".repeat(1002));
    let mut catalog = FragCatalog::new(&params(1, 2));
    assert!(is_value(
        FragCatGenerator.add_frags_from_mol(&big, &mut catalog)
    ));
    assert_eq!(
        FragCatGenerator
            .add_frags_from_mol(&mol(&"C".repeat(1001)), &mut catalog)
            .unwrap(),
        2
    );
    assert!(is_value(FragFPGenerator.get_fp_for_mol(&big, &catalog)));
    let mut catalog = FragCatalog::new(&params(1, 1));
    assert_eq!(
        FragCatGenerator
            .add_frags_from_mol(&big, &mut catalog)
            .unwrap(),
        1
    );
    assert_eq!(on_bits(&catalog, &"C".repeat(1002)), [0]);
}

/// Where a description is in the pickle (its length, then its bytes).
fn description_at(pickle: &[u8], text: &str) -> usize {
    let mut pattern = (text.len() as i32).to_le_bytes().to_vec();
    pattern.extend(text.as_bytes());
    pickle
        .windows(pattern.len())
        .position(|w| w == pattern)
        .unwrap()
}

#[test]
fn description_lengths_are_checked() {
    let good = hex(CATALOG_PICKLE);
    let at = description_at(&good, "C<-O>CC");
    for length in [-1i32, 1 << 30, i32::MAX] {
        let mut bad = good.clone();
        bad[at..at + 4].copy_from_slice(&length.to_le_bytes());
        // SAFETY: a pickle RDKit wrote but for a description length, which the wrapper checks.
        assert!(
            is_value(unsafe { FragCatalog::from_pickle(&bad) }),
            "{length}"
        );
    }
}

#[test]
fn edges_are_checked() {
    let good = hex(CATALOG_PICKLE);
    // The fourth entry's edge, to entry 4.
    let at = good.len() - 8;
    assert_eq!(&good[at..at + 4], 4i32.to_le_bytes());
    for id in [5i32, -1] {
        let mut bad = good.clone();
        bad[at..at + 4].copy_from_slice(&id.to_le_bytes());
        // SAFETY: a pickle RDKit wrote but for an edge, which the wrapper checks.
        assert!(is_value(unsafe { FragCatalog::from_pickle(&bad) }), "{id}");
    }
}

#[test]
fn truncated_pickles() {
    let pickle = hex(CATALOG_PICKLE);
    for n in [0, 20, 40, 120, 300, pickle.len() - 1] {
        // SAFETY: a prefix of a pickle RDKit wrote, which the wrapper reads to the end first.
        assert!(
            unsafe { FragCatalog::from_pickle(&pickle[..n]) }.is_err(),
            "{n}"
        );
    }
    // RDKit cannot read back a catalog with no functional groups.
    let mut catalog =
        FragCatalog::new(&FragCatParams::new(1, 3, groups_file("none.txt", "// none\n")).unwrap());
    FragCatGenerator
        .add_frags_from_mol(&mol("CCO"), &mut catalog)
        .unwrap();
    assert!(matches!(
        // SAFETY: a pickle RDKit wrote.
        unsafe { FragCatalog::from_pickle(&catalog.serialize().unwrap()) },
        Err(Error::Runtime { .. })
    ));
}

// A pickle can give a catalog bits no entry has, which RDKit's lookups
// return as null entries and Python's dereference, and bits past the
// fingerprint, which RDKit's fingerprint refuses only after allocating it.
#[test]
fn bits_are_checked() {
    let good = hex(CATALOG_PICKLE);
    let mut pickle = good.clone();
    pickle[16..20].copy_from_slice(&7u32.to_le_bytes());
    // SAFETY: a pickle RDKit wrote.
    let catalog = unsafe { FragCatalog::from_pickle(&pickle) }.unwrap();
    assert_eq!(catalog.get_bit_description(4).unwrap(), "C<-O>CC");
    assert!(is_value(catalog.get_bit_description(5)));
    assert!(is_value(catalog.get_bit_order(6)));
    assert!(is_value(catalog.get_bit_func_group_ids(5)));
    assert!(is_value(catalog.get_bit_discrims(5)));
    assert_eq!(catalog.get_bit_entry_id(6).unwrap(), None);
    assert_eq!(on_bits(&catalog, "CCCC"), [3]);

    // The last entry's bit, before its description.
    let at = description_at(&good, "C<-O>CC") - 4;
    assert_eq!(&good[at..at + 4], 4i32.to_le_bytes());
    let mut pickle = good.clone();
    pickle[at..at + 4].copy_from_slice(&5i32.to_le_bytes());
    // SAFETY: a pickle RDKit wrote.
    let catalog = unsafe { FragCatalog::from_pickle(&pickle) }.unwrap();
    assert!(is_value(
        FragFPGenerator.get_fp_for_mol(&mol("CCCC"), &catalog)
    ));

    // More bits than a fingerprint can index (U3).
    let mut pickle = good.clone();
    pickle[16..20].copy_from_slice(&(1u32 << 31).to_le_bytes());
    // SAFETY: a pickle RDKit wrote.
    let catalog = unsafe { FragCatalog::from_pickle(&pickle) }.unwrap();
    assert!(is_value(
        FragFPGenerator.get_fp_for_mol(&mol("CCCC"), &catalog)
    ));
}

#[test]
fn molecules_are_checked() {
    let good = hex(CATALOG_PICKLE);
    // Replace the first fragment's molecule pickle with one whose atom has
    // plain monomer info, whose pickled type RDKit 2024.03 reads into info of
    // that type, as Python does: PDBRESIDUE gives PDB residue info, and an
    // unknown type is RDKit's error.
    let mut mol = mol("CC");
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_monomer_info(&AtomMonomerInfo::new_with_name(AtomMonomerType::Other, "MARK").unwrap());
    let mut marked = mol.to_binary();
    let at = marked
        .windows(12)
        .position(|w| w == b"\x04\0\0\0MARK\x02\0\0\0")
        .unwrap()
        + 8;
    // The first entry's molecule starts after the parameters' text and ends
    // before its bit id, which precedes its description.
    let start = good
        .windows(4)
        .skip(24)
        .position(|w| w == b"\xef\xbe\xad\xde")
        .unwrap()
        + 24;
    let end = description_at(&good, "C<-O>C") - 4;
    let splice = |molecule: &[u8]| [&good[..start], molecule, &good[end..]].concat();
    // SAFETY: a pickle RDKit wrote with another molecule pickle RDKit wrote spliced in.
    assert!(unsafe { FragCatalog::from_pickle(&splice(&marked)) }.is_ok());
    marked[at] = 1;
    // SAFETY: as above, but for a monomer type, which RDKit 2024.03 reads
    // into info of that type.
    assert!(unsafe { FragCatalog::from_pickle(&splice(&marked)) }.is_ok());
    marked[at] = 7;
    // SAFETY: as above; RDKit refuses the type.
    assert!(unsafe { FragCatalog::from_pickle(&splice(&marked)) }.is_err());
}

/// `smiles` with a substance group of `atoms` and their bonds, with `n`
/// brackets.
fn with_brackets(smiles: &str, atoms: u32, n: usize) -> Mol {
    let mut mol = mol(smiles);
    let mut sgroup = create_mol_substance_group(&mut mol, "GEN");
    for idx in 0..atoms {
        sgroup.add_atom_with_idx(idx).unwrap();
    }
    for idx in 0..atoms - 1 {
        sgroup.add_bond_with_idx(idx).unwrap();
    }
    let bracket = [
        rdkit_cxx::geometry::Point3D::default(),
        rdkit_cxx::geometry::Point3D::new(1.0, 0.0, 0.0),
    ];
    for _ in 0..n {
        sgroup.add_bracket(&bracket).unwrap();
    }
    mol
}

// RDKit 2024.03 makes the entries' fragments with Subgraphs::pathToSubmol,
// which keeps no substance or stereo groups, so a group a pickle could not
// hold does not reach them; the atoms keep their properties, which the
// pickle counts in two bytes.
#[test]
fn fragments_are_picklable() {
    let generator = FragCatGenerator;
    let mut catalog = FragCatalog::new(&params(1, 3));
    let small = with_brackets("OCC=CC(=O)O", 2, 256);
    assert!(generator.add_frags_from_mol(&small, &mut catalog).unwrap() > 0);
    // SAFETY: a pickle RDKit wrote.
    let read = unsafe { FragCatalog::from_pickle(&catalog.serialize().unwrap()) }.unwrap();
    assert_eq!(read.get_num_entries(), catalog.get_num_entries());

    let mut many = mol("OCC=CC(=O)O");
    {
        let mut atom = many.get_atom_with_idx_mut(0).unwrap();
        for i in 0..65536 {
            atom.set_prop(&format!("p{i}"), "x");
        }
    }
    let mut catalog = FragCatalog::new(&params(1, 3));
    let Err(Error::Value { message }) = generator.add_frags_from_mol(&many, &mut catalog) else {
        panic!("a fragment's pickle would wrap its atom's property count");
    };
    assert!(message.contains("properties"), "{message}");
    assert_eq!(catalog.get_num_entries(), 0);
}
