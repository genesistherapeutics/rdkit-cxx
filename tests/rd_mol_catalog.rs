// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::rd_mol_catalog::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdMolCatalog.

/// Python's pickle of `catalog()`.
const CATALOG_PICKLE: &str = "efbeadde0100000000000000000000000300000003000000efbeadde00000000100000000200000000000000030000000200000080010600600000000103060060000000020208006000000001010b000100010200420000000017040000000000000016000000000100000007000000657468616e6f6cefbeadde00000000100000000200000000000000040000000300000080010600600000000103060028000000030408002800000003020800680000000301010b0001000102280201032042000000001704000000000000001601000000020000000b0000006163657469632061636964efbeadde00000000100000000200000000000000060000000600000080010640680000000303010640680000000303010640680000000303010640680000000303010640680000000303010640680000000303010b0001680c0102680c0203680c0304680c0405680c0500680c4201000000060005040302011704000000000000001602000000020000000700000062656e7a656e650200000001000000020000000000000000000000";

/// Python's pickle of `ethanol()`.
const ENTRY_PICKLE: &str = "efbeadde00000000100000000200000000000000030000000200000080010600600000000103060060000000020208006000000001010b000100010200420000000017040000000000000016ffffffff0300000007000000657468616e6f6c";

fn hex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn entry(smiles: &str, description: &str, order: u32) -> MolCatalogEntry {
    let mut entry = MolCatalogEntry::default();
    entry.set_mol(&mol_from_smiles(smiles).unwrap());
    entry.set_description(description);
    entry.set_order(order);
    entry
}

fn ethanol() -> MolCatalogEntry {
    entry("CCO", "ethanol", 3)
}

fn catalog() -> MolCatalog {
    let mut catalog = create_mol_catalog();
    assert_eq!(catalog.add_entry(&entry("CCO", "ethanol", 1)), 0);
    assert_eq!(catalog.add_entry(&entry("CC(=O)O", "acetic acid", 2)), 1);
    assert_eq!(catalog.add_entry(&entry("c1ccccc1", "benzene", 2)), 2);
    catalog.add_edge(0, 1).unwrap();
    catalog.add_edge(0, 2).unwrap();
    catalog.add_edge(0, 1).unwrap();
    catalog
}

fn check_catalog(catalog: &MolCatalog) {
    assert_eq!(catalog.get_num_entries(), 3);
    assert_eq!(catalog.get_fp_length(), 3);
    let down: Vec<_> = (0..3)
        .map(|i| catalog.get_entry_down_ids(i).unwrap())
        .collect();
    assert_eq!(down, vec![vec![1, 2], vec![], vec![]]);
    let names = ["ethanol", "acetic acid", "benzene"];
    for (i, name) in (0..3).zip(names) {
        assert_eq!(catalog.get_bit_description(i).unwrap(), name);
        assert_eq!(catalog.get_entry_description(i).unwrap(), name);
        assert_eq!(catalog.get_bit_entry_id(i).unwrap(), Some(i));
        assert_eq!(catalog.get_entry_bit_id(i).unwrap(), Some(i));
    }
}

fn is_index<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Index { .. }))
}

fn is_value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

#[test]
fn catalogs() {
    let catalog = catalog();
    check_catalog(&catalog);
    assert_eq!(catalog.serialize().unwrap(), hex(CATALOG_PICKLE));
    // SAFETY: a pickle RDKit wrote.
    check_catalog(&unsafe { MolCatalog::from_pickle(&hex(CATALOG_PICKLE)) }.unwrap());
    // SAFETY: as above.
    check_catalog(&unsafe { MolCatalog::from_pickle(&catalog.serialize().unwrap()) }.unwrap());

    let empty = create_mol_catalog();
    assert_eq!((empty.get_num_entries(), empty.get_fp_length()), (0, 0));
    // SAFETY: a pickle RDKit wrote.
    let empty = unsafe { MolCatalog::from_pickle(&empty.serialize().unwrap()) }.unwrap();
    assert_eq!(empty.get_num_entries(), 0);
}

// Python checks each index one past the end, where RDKit's own check
// raises; GetEntryDownIds has none, so it reads past the catalog's graph.
#[test]
fn indices_are_checked() {
    let mut catalog = catalog();
    assert!(is_index(catalog.get_bit_description(3)));
    assert!(is_index(catalog.get_bit_entry_id(3)));
    assert!(is_index(catalog.get_entry_bit_id(3)));
    assert!(is_index(catalog.get_entry_description(3)));
    assert!(is_index(catalog.get_entry_down_ids(3)));
    assert!(is_index(catalog.add_edge(0, 3)));
    assert!(is_index(catalog.add_edge(3, 0)));
    assert!(is_index(create_mol_catalog().get_entry_down_ids(0)));
}

#[test]
fn entries() {
    let entry = MolCatalogEntry::default();
    assert_eq!(entry.get_description().unwrap(), "");
    assert_eq!(entry.get_order(), 0);
    // Python's GetMol dereferences null here.
    assert!(entry.get_mol().is_none());
    assert!(matches!(entry.serialize(), Err(Error::Runtime { .. })));

    let mut entry = ethanol();
    assert_eq!(entry.get_description().unwrap(), "ethanol");
    assert_eq!(entry.get_order(), 3);
    assert_eq!(mol_to_smiles(&entry.get_mol().unwrap()).unwrap(), "CCO");
    let pickle = entry.serialize().unwrap();
    assert_eq!(pickle, hex(ENTRY_PICKLE));
    for pickle in [pickle, hex(ENTRY_PICKLE)] {
        // SAFETY: a pickle RDKit wrote.
        let back = unsafe { MolCatalogEntry::from_pickle(&pickle) }.unwrap();
        assert_eq!(back.get_description().unwrap(), "ethanol");
        assert_eq!(back.get_order(), 3);
        assert_eq!(mol_to_smiles(&back.get_mol().unwrap()).unwrap(), "CCO");
    }

    // The entry holds a copy: changing the molecule changes neither.
    let mut mol = entry.get_mol().unwrap();
    mol.get_atom_with_idx_mut(0).unwrap().set_atomic_num(7);
    assert_eq!(mol_to_smiles(&entry.get_mol().unwrap()).unwrap(), "CCO");
    entry.set_mol(&mol);
    assert_eq!(mol_to_smiles(&entry.get_mol().unwrap()).unwrap(), "NCO");

    // The catalog holds a copy of the entry.
    let mut catalog = create_mol_catalog();
    catalog.add_entry(&entry);
    entry.set_description("changed");
    assert_eq!(catalog.get_entry_description(0).unwrap(), "ethanol");
    // An entry without a molecule can be added but not pickled.
    catalog.add_entry(&MolCatalogEntry::default());
    assert!(matches!(catalog.serialize(), Err(Error::Runtime { .. })));
}

/// Where the entry pickle's description length is.
fn description_length_at(pickle: &[u8]) -> usize {
    pickle.len() - "ethanol".len() - 4
}

// RDKit allocates a description's length plus one bytes and writes a NUL at
// its length: a negative length writes before the buffer, and one past the
// end leaves bytes unread.
#[test]
fn description_lengths_are_checked() {
    let good = hex(ENTRY_PICKLE);
    let at = description_length_at(&good);
    for length in [-1i32, -100_000, 8, 1 << 30, i32::MAX] {
        let mut bad = good.clone();
        bad[at..at + 4].copy_from_slice(&length.to_le_bytes());
        // SAFETY: a pickle RDKit wrote but for a description length, which the wrapper checks.
        assert!(
            is_value(unsafe { MolCatalogEntry::from_pickle(&bad) }),
            "{length}"
        );
    }
    // Shorter lengths stop early, as in Python, which ignores the rest.
    let mut short = good.clone();
    short[at..at + 4].copy_from_slice(&3i32.to_le_bytes());
    // SAFETY: as above, with a length the wrapper accepts.
    let entry = unsafe { MolCatalogEntry::from_pickle(&short) }.unwrap();
    assert_eq!(entry.get_description().unwrap(), "eth");

    // The catalog's entries are checked the same way: the last entry's
    // description ends 20 bytes before the edges' end.
    let catalog = hex(CATALOG_PICKLE);
    let at = catalog.len() - 20 - "benzene".len() - 4;
    assert_eq!(&catalog[at..at + 4], 7i32.to_le_bytes());
    let mut bad = catalog.clone();
    bad[at..at + 4].copy_from_slice(&(-1i32).to_le_bytes());
    // SAFETY: a pickle RDKit wrote but for a description length, which the wrapper checks.
    assert!(is_value(unsafe { MolCatalog::from_pickle(&bad) }));
}

#[test]
fn edges_are_checked() {
    let good = hex(CATALOG_PICKLE);
    // The first entry's second edge, to entry 2.
    let at = good.len() - 12;
    assert_eq!(&good[at..at + 4], 2i32.to_le_bytes());
    for id in [3i32, -1, i32::MAX] {
        let mut bad = good.clone();
        bad[at..at + 4].copy_from_slice(&id.to_le_bytes());
        // SAFETY: a pickle RDKit wrote but for an edge, which the wrapper checks.
        assert!(is_value(unsafe { MolCatalog::from_pickle(&bad) }), "{id}");
    }
}

#[test]
fn truncated_pickles() {
    let catalog = hex(CATALOG_PICKLE);
    let entry = hex(ENTRY_PICKLE);
    for n in [0, 10, 30, 100, catalog.len() - 1] {
        // SAFETY: a prefix of a pickle RDKit wrote, which the wrapper reads to the end first.
        assert!(
            unsafe { MolCatalog::from_pickle(&catalog[..n]) }.is_err(),
            "{n}"
        );
    }
    for n in [0, 10, 50, entry.len() - 8] {
        // SAFETY: as above.
        assert!(
            unsafe { MolCatalogEntry::from_pickle(&entry[..n]) }.is_err(),
            "{n}"
        );
    }
}

// A pickle can give the catalog bits no entry has; RDKit returns a null
// entry, which Python's GetBitDescription dereferences.
#[test]
fn bits_without_entries() {
    let mut pickle = hex(CATALOG_PICKLE);
    pickle[16..20].copy_from_slice(&5u32.to_le_bytes());
    // SAFETY: a pickle RDKit wrote.
    let catalog = unsafe { MolCatalog::from_pickle(&pickle) }.unwrap();
    assert_eq!(catalog.get_fp_length(), 5);
    assert_eq!(catalog.get_bit_description(2).unwrap(), "benzene");
    assert!(is_value(catalog.get_bit_description(3)));
    assert_eq!(catalog.get_bit_entry_id(4).unwrap(), None);
    assert!(is_index(catalog.get_bit_description(5)));
}

/// An entry pickle whose molecule's atom has monomer info of type
/// PDBRESIDUE but no residue fields, which RDKit 2024.09's pickler reads as
/// PDB residue info, reading the fields from the bytes after it (U10).
fn pdb_typed_entry_pickle() -> Vec<u8> {
    let mut mol = mol_from_smiles("C").unwrap();
    mol.get_atom_with_idx_mut(0)
        .unwrap()
        .set_monomer_info(&AtomMonomerInfo::new_with_name(AtomMonomerType::Other, "MARK").unwrap());
    let mut entry = MolCatalogEntry::default();
    entry.set_mol(&mol);
    let mut pickle = entry.serialize().unwrap();
    // The info's name (a u32 length, then its bytes), then its type.
    let at = pickle
        .windows(9)
        .position(|w| w == b"\x04\0\0\0MARK\x02")
        .unwrap();
    // SAFETY: a pickle RDKit wrote.
    assert!(unsafe { MolCatalogEntry::from_pickle(&pickle) }.is_ok());
    pickle[at + 8] = 1;
    pickle
}

// As Python's pickle.loads, the info comes back as PDB residue info.
#[test]
fn pdb_typed_molecules_have_pdb_residue_info() {
    let pickle = pdb_typed_entry_pickle();
    // SAFETY: a pickle RDKit wrote but for a monomer type byte, which RDKit
    // reads with its own checks.
    let entry = unsafe { MolCatalogEntry::from_pickle(&pickle) }.unwrap();
    let mol = entry.get_mol().unwrap();
    let atom = mol.get_atom_with_idx(0).unwrap();
    let info = atom.get_pdb_residue_info().unwrap().unwrap();
    assert_eq!(info.get_name().unwrap(), "MARK");
    // A catalog of one such entry, with no edges.
    let mut catalog = create_mol_catalog().serialize().unwrap();
    catalog[20..24].copy_from_slice(&1u32.to_le_bytes());
    catalog.extend(&pickle);
    catalog.extend(0u32.to_le_bytes());
    // SAFETY: a pickle RDKit wrote.
    let catalog = unsafe { MolCatalog::from_pickle(&catalog) }.unwrap();
    assert_eq!(catalog.get_num_entries(), 1);
}
