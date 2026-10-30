// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.MACCSkeys in RDKit 2025.09.6, and its
// doctests.

use rdkit_cxx::chem::maccs_keys::{maccs_keys, smarts_patts};
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_mol_descriptors::get_maccs_keys_fingerprint;

#[test]
fn smarts_patts_match_python() {
    let patts = smarts_patts();
    assert_eq!(patts.len(), 166);
    // The SHA-1 of Python's "key smarts count" lines, one per key in order.
    let text: Vec<String> = patts
        .iter()
        .map(|(key, (smarts, count))| format!("{key} {smarts} {count}"))
        .collect();
    let digest = sha1_smol::Sha1::from(text.join("\n")).digest().to_string();
    assert_eq!(digest, "61cc4a0974d8e7a42cfddca881f13688360bbea8");
    assert_eq!(patts[&118], ("[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]", 1));
    let unset: Vec<u32> = patts
        .iter()
        .filter(|(_, (s, _))| *s == "?")
        .map(|(k, _)| *k)
        .collect();
    assert_eq!(unset, [1, 125, 166]);
    assert_eq!(patts.values().filter(|(_, count)| *count > 0).count(), 14);
}

/// Python's `_pyGenMACCSKeys`: the keys from the patterns, with keys 125
/// (more than one aromatic ring) and 166 (more than one fragment) computed.
fn py_gen_maccs_keys(smiles: &str) -> Vec<u32> {
    let mol = mol_from_smiles(smiles).unwrap();
    let mut on = Vec::new();
    for (i, (patt, count)) in maccs_keys().iter().enumerate() {
        let key = i as u32 + 1;
        let set = match patt {
            Some(patt) if *count == 0 => mol.has_substruct_match(patt).unwrap(),
            Some(patt) => mol.get_substruct_matches(patt).unwrap().len() > *count as usize,
            None if key == 125 => {
                let rings = mol.get_ring_info().bond_rings();
                let aromatic = rings.iter().filter(|ring| {
                    ring.iter()
                        .all(|&b| mol.get_bond_with_idx(b).unwrap().get_is_aromatic())
                });
                aromatic.count() > 1
            }
            None if key == 166 => rdkit_cxx::chem::get_mol_frags(&mol).len() > 1,
            None => false,
        };
        if set {
            on.push(key);
        }
    }
    on
}

#[test]
fn maccs_keys_fingerprint_molecules() {
    let keys = maccs_keys();
    assert_eq!(keys.len(), 166);
    assert!(keys[0].0.is_none() && keys[124].0.is_none() && keys[165].0.is_none());
    assert_eq!(keys[117].1, 1);
    for (smiles, expected) in [
        (
            "CNO",
            vec![
                24, 68, 69, 71, 93, 94, 102, 124, 131, 139, 151, 158, 160, 161, 164,
            ],
        ),
        ("CCC", vec![74, 114, 149, 155, 160]),
    ] {
        assert_eq!(py_gen_maccs_keys(smiles), expected);
        let fp = get_maccs_keys_fingerprint(&mol_from_smiles(smiles).unwrap()).unwrap();
        assert_eq!(fp.get_on_bits(), expected);
    }
    // The patterns give rdMolDescriptors' keys, as Python's do.
    for smiles in [
        "c1ccccc1-c1ccncc1.[Na+]",
        "OC(=O)C1CCCCC1N",
        "FC(F)(F)c1ccc2[nH]ccc2c1",
    ] {
        let fp = get_maccs_keys_fingerprint(&mol_from_smiles(smiles).unwrap()).unwrap();
        assert_eq!(py_gen_maccs_keys(smiles), fp.get_on_bits(), "{smiles}");
    }
}
