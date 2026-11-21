// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The Cookbook's "Descriptor Calculations" chapter.

use std::collections::{BTreeMap, BTreeSet};

use rdkit_cxx::Result;
use rdkit_cxx::chem::rd_mol_hash::{HashFunction, mol_hash};
use rdkit_cxx::chem::{Mol, mol_from_smarts, mol_from_smiles};

use crate::Recipe;

pub const RECIPES: &[Recipe] = &[
    ("RDKitCB_21", "Molecule Hash Strings", molecule_hash_strings),
    ("RDKitCB_22", "Contiguous Rotatable Bonds", contiguous_rotatable_bonds),
];

/// Every hash function, by its name in Python's `HashFunction.names`.
const HASH_FUNCTIONS: [(&str, HashFunction); 19] = [
    ("AnonymousGraph", HashFunction::AnonymousGraph),
    ("ElementGraph", HashFunction::ElementGraph),
    ("CanonicalSmiles", HashFunction::CanonicalSmiles),
    ("MurckoScaffold", HashFunction::MurckoScaffold),
    ("ExtendedMurcko", HashFunction::ExtendedMurcko),
    ("MolFormula", HashFunction::MolFormula),
    ("AtomBondCounts", HashFunction::AtomBondCounts),
    ("DegreeVector", HashFunction::DegreeVector),
    ("Mesomer", HashFunction::Mesomer),
    ("HetAtomTautomer", HashFunction::HetAtomTautomer),
    ("HetAtomProtomer", HashFunction::HetAtomProtomer),
    ("RedoxPair", HashFunction::RedoxPair),
    ("Regioisomer", HashFunction::Regioisomer),
    ("NetCharge", HashFunction::NetCharge),
    ("SmallWorldIndexBR", HashFunction::SmallWorldIndexBr),
    ("SmallWorldIndexBRL", HashFunction::SmallWorldIndexBrl),
    ("ArthorSubstructureOrder", HashFunction::ArthorSubstructureOrder),
    ("HetAtomTautomerv2", HashFunction::HetAtomTautomerv2),
    ("HetAtomProtomerv2", HashFunction::HetAtomProtomerv2),
];

fn hash_all(smiles: &[&str], func: HashFunction) -> Result<Vec<String>> {
    smiles
        .iter()
        .map(|smi| mol_hash(&mol_from_smiles(smi)?, func))
        .collect()
}

/// Hashes molecules with NextMove's MolHash, and groups molecules by their
/// Murcko scaffold and regioisomer hashes. By Vincent Scalfani and Takayuki
/// Serizawa, with examples from O'Boyle and Sayle.
fn molecule_hash_strings() -> Result<()> {
    let s = mol_from_smiles("CC(C(C1=CC(=C(C=C1)O)O)O)N(C)C(=O)OCC2=CC=CC=C2")?;
    let expected = [
        "**1***(*(*)*(*)*(*)*(*)***2*****2)**1*",
        "CC(C(O)C1CCC(O)C(O)C1)N(C)C(O)OCC1CCCCC1",
        "CC(C(O)c1ccc(O)c(O)c1)N(C)C(=O)OCc1ccccc1",
        "c1ccc(CCNCOCc2ccccc2)cc1",
        "*c1ccc(C(*)C(*)N(*)C(=*)OCc2ccccc2)cc1*",
        "C18H21NO5",
        "24,25",
        "0,8,10,6",
        "CC(C(O)[C]1[CH][CH][C](O)[C](O)[CH]1)N(C)[C]([O])OC[C]1[CH][CH][CH][CH][CH]1_0",
        "CC(C([O])[C]1[CH][CH][C]([O])[C]([O])[CH]1)N(C)[C]([O])OC[C]1[CH][CH][CH][CH][CH]1_3_0",
        "CC(C([O])[C]1[CH][CH][C]([O])[C]([O])[CH]1)N(C)[C]([O])OC[C]1[CH][CH][CH][CH][CH]1_3",
        "CC(C(O)[C]1[CH][CH][C](O)[C](O)[CH]1)N(C)[C]([O])OC[C]1[CH][CH][CH][CH][CH]1",
        "*C.*CCC.*O.*O.*O.*OC(=O)N(*)*.C.c1ccccc1.c1ccccc1",
        "0",
        "B25R2",
        "B25R2L10",
        "00180019010012000600009b000000",
        "[CH3]-[CH](-[CH](-[OH])-[C]1:[C]:[C]:[C](:[O]):[C](:[O]):[C]:1)-[N](-[CH3])-[C](=[O])\
         -[O]-[CH2]-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1_5_0",
        "[CH3]-[CH](-[CH](-[OH])-[C]1:[C]:[C]:[C](:[O]):[C](:[O]):[C]:1)-[N](-[CH3])-[C](=[O])\
         -[O]-[CH2]-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1_5",
    ];
    for ((name, func), expected) in HASH_FUNCTIONS.iter().zip(expected) {
        let hash = mol_hash(&s, *func)?;
        println!("{name} {hash}");
        assert_eq!(hash, expected, "{name}");
    }

    // Murcko scaffold hashes (slide 16 of O'Boyle and Sayle).
    let m_list = [
        "CCC1CC(CCC1=O)C(=O)C1=CC=CC(C)=C1",
        "CCC1CC(CCC1=O)C(=O)C1=CC=CC=C1",
        "CC(=C)C(C1=CC=CC=C1)S(=O)CC(N)=O",
        "CC1=CC(=CC=C1)C(C1CCC(N)CC1)C(F)(F)F",
        "CNC1CCC(C2=CC(Cl)=C(Cl)C=C2)C2=CC=CC=C12",
        "CCCOC(C1CCCCC1)C1=CC=C(Cl)C=C1",
    ];
    let murcko_hashes = hash_all(&m_list, HashFunction::MurckoScaffold)?;
    println!("{murcko_hashes:?}");
    assert_eq!(
        murcko_hashes,
        [
            "c1ccc(CC2CCCCC2)cc1",
            "c1ccc(CC2CCCCC2)cc1",
            "c1ccccc1",
            "c1ccc(CC2CCCCC2)cc1",
            "c1ccc(C2CCCc3ccccc32)cc1",
            "c1ccc(CC2CCCCC2)cc1",
        ]
    );

    // The most frequent scaffold, and the atoms each molecule matches it on.
    let mut counts = BTreeMap::new();
    for hash in &murcko_hashes {
        *counts.entry(hash.as_str()).or_insert(0) += 1;
    }
    let most_frequent = counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(hash, _)| hash)
        .unwrap();
    println!("{most_frequent}");
    assert_eq!(most_frequent, "c1ccc(CC2CCCCC2)cc1");
    let scaffold = mol_from_smiles(most_frequent)?;
    for smi in m_list {
        let highlight = mol_from_smiles(smi)?.get_substruct_match(&scaffold)?;
        println!("{highlight:?}");
    }

    // Regioisomer hashes (slide 17): which molecules are regioisomers of r0?
    let r_list = [
        "CC1=CC2=C(C=C1)C(=CN2CCN1CCOCC1)C(=O)C1=CC=CC2=C1C=CC=C2",
        "CCCCCN1C=C(C2=CC=CC=C21)C(=O)C3=CC=CC4=CC=CC=C43",
        "CC1COCCN1CCN1C=C(C(=O)C2=CC=CC3=C2C=CC=C3)C2=C1C=CC=C2",
        "CC1=CC=C(C(=O)C2=CN(CCN3CCOCC3)C3=C2C=CC=C3)C2=C1C=CC=C2",
        "CC1=C(CCN2CCOCC2)C2=C(C=CC=C2)N1C(=O)C1=CC=CC2=CC=CC=C12",
        "CN1CCN(C(C1)CN2C=C(C3=CC=CC=C32)C(=O)C4=CC=CC5=CC=CC=C54)C",
    ];
    let regio_hashes = hash_all(&r_list, HashFunction::Regioisomer)?;
    let r0_hash = &regio_hashes[0];
    println!("{r0_hash}");
    assert_eq!(
        r0_hash,
        "*C.*C(*)=O.*CC*.C1COCCN1.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1"
    );
    println!("{regio_hashes:?}");
    let is_regioisomer: Vec<bool> = regio_hashes.iter().map(|hash| hash == r0_hash).collect();
    for (idx, matched) in is_regioisomer.iter().enumerate() {
        println!("r{idx}: Regioisomer: {}", if *matched { "True" } else { "False" });
    }
    assert_eq!(is_regioisomer, [true, false, true, true, true, false]);
    Ok(())
}

/// SMARTS for rotatable bonds: single, acyclic, between atoms that are
/// neither terminal nor triple-bonded. Python's `Lipinski.RotatableBondSmarts`.
const ROTATABLE_BOND_SMARTS: &str = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]";

/// Groups of contiguous rotatable bonds, largest first.
fn find_bond_groups(mol: &Mol) -> Result<Vec<BTreeSet<u32>>> {
    let rotatable = mol_from_smarts(ROTATABLE_BOND_SMARTS)?;
    let mut rot_bond_set = BTreeSet::new();
    for pair in mol.get_substruct_matches(&rotatable)? {
        let bond = mol.get_bond_between_atoms(pair[0], pair[1])?.unwrap();
        rot_bond_set.insert(bond.get_idx());
    }
    let mut rot_bond_groups = Vec::new();
    while let Some(first) = rot_bond_set.pop_first() {
        let mut connected = BTreeSet::from([first]);
        let mut stack = vec![first];
        while let Some(idx) = stack.pop() {
            let bond = mol.get_bond_with_idx(idx)?;
            for atom in [bond.get_begin_atom(), bond.get_end_atom()] {
                for neighbor in atom.get_bonds() {
                    let idx = neighbor.get_idx();
                    if rot_bond_set.remove(&idx) {
                        connected.insert(idx);
                        stack.push(idx);
                    }
                }
            }
        }
        rot_bond_groups.push(connected);
    }
    rot_bond_groups.sort_by_key(|group| std::cmp::Reverse(group.len()));
    Ok(rot_bond_groups)
}

/// Finds the longest run of contiguous rotatable bonds. By Paolo Tosco.
fn contiguous_rotatable_bonds() -> Result<()> {
    let mol = mol_from_smiles("CCC(CC(C)CC1CCC1)C(CC(=O)O)N")?;
    let bond_groups = find_bond_groups(&mol)?;
    // Groups come largest first, so the first is the longest run.
    let largest = bond_groups.first().map_or(0, BTreeSet::len);
    println!("{largest}");
    assert_eq!(largest, 8);
    println!("{bond_groups:?}");
    assert_eq!(bond_groups, [BTreeSet::from([1, 2, 3, 5, 6, 10, 11, 12])]);
    Ok(())
}
