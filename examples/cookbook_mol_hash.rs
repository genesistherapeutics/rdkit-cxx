// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Molecule Hash Strings (RDKitCB_21): hash molecules with NextMove's MolHash
//! functions, and find shared Murcko scaffolds and regioisomers.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#molecule-hash-strings),
//! by Vincent Scalfani and Takayuki Serizawa, with examples from O'Boyle and
//! Sayle, "Making a hash of it" (ACS, 2019). The recipe draws the molecules;
//! this prints their hashes.

use rdkit_cxx::Result;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_mol_hash::{HashFunction, mol_hash};

/// Python's `rdMolHash.HashFunction.names`, in its order.
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
    (
        "ArthorSubstructureOrder",
        HashFunction::ArthorSubstructureOrder,
    ),
    ("HetAtomTautomerv2", HashFunction::HetAtomTautomerv2),
    ("HetAtomProtomerv2", HashFunction::HetAtomProtomerv2),
];

/// The most frequent item. Python's `max(set(list), key=list.count)` breaks
/// ties in set order; this takes the first.
fn most_freq(list: &[String]) -> &str {
    list.iter()
        .rev()
        .max_by_key(|item| list.iter().filter(|x| x == item).count())
        .expect("the list is not empty")
}

fn main() -> Result<()> {
    let s = mol_from_smiles("CC(C(C1=CC(=C(C=C1)O)O)O)N(C)C(=O)OCC2=CC=CC=C2")?;

    // Hash 's' with every hash function.
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
        "[CH3]-[CH](-[CH](-[OH])-[C]1:[C]:[C]:[C](:[O]):[C](:[O]):[C]:1)-[N](-[CH3])-[C](=[O])-[O]-[CH2]-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1_5_0",
        "[CH3]-[CH](-[CH](-[OH])-[C]1:[C]:[C]:[C](:[O]):[C](:[O]):[C]:1)-[N](-[CH3])-[C](=[O])-[O]-[CH2]-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1_5",
    ];
    for ((name, func), expected) in HASH_FUNCTIONS.iter().zip(expected) {
        let hash = mol_hash(&s, *func)?;
        println!("{name} {hash}");
        assert_eq!(hash, expected, "{name}");
    }

    // Murcko scaffold hashes (from slide 16 of O'Boyle and Sayle)
    let m_list = [
        "CCC1CC(CCC1=O)C(=O)C1=CC=CC(C)=C1",
        "CCC1CC(CCC1=O)C(=O)C1=CC=CC=C1",
        "CC(=C)C(C1=CC=CC=C1)S(=O)CC(N)=O",
        "CC1=CC(=CC=C1)C(C1CCC(N)CC1)C(F)(F)F",
        "CNC1CCC(C2=CC(Cl)=C(Cl)C=C2)C2=CC=CC=C12",
        "CCCOC(C1CCCCC1)C1=CC=C(Cl)C=C1",
    ];
    let murcko_hash_list = m_list
        .iter()
        .map(|m| mol_hash(&mol_from_smiles(m)?, HashFunction::MurckoScaffold))
        .collect::<Result<Vec<_>>>()?;
    println!("{murcko_hash_list:?}");
    assert_eq!(
        murcko_hash_list,
        [
            "c1ccc(CC2CCCCC2)cc1",
            "c1ccc(CC2CCCCC2)cc1",
            "c1ccccc1",
            "c1ccc(CC2CCCCC2)cc1",
            "c1ccc(C2CCCc3ccccc32)cc1",
            "c1ccc(CC2CCCCC2)cc1",
        ]
    );

    // The most frequent Murcko scaffold hash
    let most_freq_murcko_hash = most_freq(&murcko_hash_list);
    println!("{most_freq_murcko_hash}");
    assert_eq!(most_freq_murcko_hash, "c1ccc(CC2CCCCC2)cc1");

    // Regioisomer hashes (from slide 17 of O'Boyle and Sayle): find the
    // regioisomers of r0.
    let r_list = [
        "CC1=CC2=C(C=C1)C(=CN2CCN1CCOCC1)C(=O)C1=CC=CC2=C1C=CC=C2",
        "CCCCCN1C=C(C2=CC=CC=C21)C(=O)C3=CC=CC4=CC=CC=C43",
        "CC1COCCN1CCN1C=C(C(=O)C2=CC=CC3=C2C=CC=C3)C2=C1C=CC=C2",
        "CC1=CC=C(C(=O)C2=CN(CCN3CCOCC3)C3=C2C=CC=C3)C2=C1C=CC=C2",
        "CC1=C(CCN2CCOCC2)C2=C(C=CC=C2)N1C(=O)C1=CC=CC2=CC=CC=C12",
        "CN1CCN(C(C1)CN2C=C(C3=CC=CC=C32)C(=O)C4=CC=CC5=CC=CC=C54)C",
    ];
    let regio_hash_list = r_list
        .iter()
        .map(|r| mol_hash(&mol_from_smiles(r)?, HashFunction::Regioisomer))
        .collect::<Result<Vec<_>>>()?;
    let r0_regio_hash = &regio_hash_list[0];
    println!("{r0_regio_hash}");
    assert_eq!(
        r0_regio_hash,
        "*C.*C(*)=O.*CC*.C1COCCN1.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1"
    );
    println!("{regio_hash_list:?}");
    assert_eq!(
        regio_hash_list,
        [
            "*C.*C(*)=O.*CC*.C1COCCN1.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1",
            "*C(*)=O.*CCCCC.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1",
            "*C.*C(*)=O.*CC*.C1COCCN1.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1",
            "*C.*C(*)=O.*CC*.C1COCCN1.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1",
            "*C.*C(*)=O.*CC*.C1COCCN1.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1",
            "*C.*C.*C(*)=O.*C*.C1CNCCN1.c1ccc2[nH]ccc2c1.c1ccc2ccccc2c1",
        ]
    );

    let rmatches: Vec<bool> = regio_hash_list
        .iter()
        .map(|regio_hash| regio_hash == r0_regio_hash)
        .collect();
    for (i, rmatch) in rmatches.iter().enumerate() {
        println!(
            "r{i}: Regioisomer: {}",
            if *rmatch { "True" } else { "False" }
        );
    }
    assert_eq!(rmatches, [true, false, true, true, true, false]);
    Ok(())
}
