// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::chem::rd_mol_hash::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdMolHash.

#[test]
fn hashes() {
    let mol = mol_from_smiles("C1CCCC(O)C1c1ccnc(OC)c1").unwrap();
    let expected = [
        (HashFunction::AnonymousGraph, "***1****(*2*****2*)*1"),
        (HashFunction::ElementGraph, "COC1CC(C2CCCCC2O)CCN1"),
        (HashFunction::CanonicalSmiles, "COc1cc(C2CCCCC2O)ccn1"),
        (HashFunction::MurckoScaffold, "c1cc(C2CCCCC2)ccn1"),
        (HashFunction::ExtendedMurcko, "*c1cc(C2CCCCC2*)ccn1"),
        (HashFunction::MolFormula, "C12H17NO2"),
        (HashFunction::AtomBondCounts, "15,16"),
        (HashFunction::DegreeVector, "0,4,9,2"),
        (
            HashFunction::Mesomer,
            "CO[C]1[CH][C](C2CCCCC2O)[CH][CH][N]1_0",
        ),
        (
            HashFunction::HetAtomTautomer,
            "CO[C]1[CH][C](C2CCCCC2[O])[CH][CH][N]1_1_0",
        ),
        (
            HashFunction::HetAtomProtomer,
            "CO[C]1[CH][C](C2CCCCC2[O])[CH][CH][N]1_1",
        ),
        (
            HashFunction::RedoxPair,
            "CO[C]1[CH][C](C2CCCCC2O)[CH][CH][N]1",
        ),
        (HashFunction::Regioisomer, "*O.*O*.C.C1CCCCC1.c1ccncc1"),
        (HashFunction::NetCharge, "0"),
        (HashFunction::SmallWorldIndexBr, "B16R2"),
        (HashFunction::SmallWorldIndexBrl, "B16R2L9"),
        (
            HashFunction::ArthorSubstructureOrder,
            "000f001001000c000300005f000000",
        ),
        (
            HashFunction::HetAtomTautomerv2,
            "[CH3]-[O]:[C]1:[C]:[C](-[CH]2-[CH2]-[CH2]-[CH2]-[CH2]-[CH]-2-[OH]):[C]:[C]:[N]:1_3_0",
        ),
        (
            HashFunction::HetAtomProtomerv2,
            "[CH3]-[O]:[C]1:[C]:[C](-[CH]2-[CH2]-[CH2]-[CH2]-[CH2]-[CH]-2-[OH]):[C]:[C]:[N]:1_3",
        ),
    ];
    for (func, hash) in expected {
        assert_eq!(mol_hash(&mol, func).unwrap(), hash, "{func:?}");
    }
    // The molecule is unchanged: MolHash edits a copy.
    assert_eq!(mol_to_smiles(&mol).unwrap(), "COc1cc(C2CCCCC2O)ccn1");
}

#[test]
fn cx_smiles() {
    let mol = mol_from_smiles("C[C@H](F)Cl |o1:1|").unwrap();
    let cx = MolHashParams {
        use_cx_smiles: true,
        ..Default::default()
    };
    assert_eq!(
        mol_hash_with_params(&mol, HashFunction::CanonicalSmiles, &cx).unwrap(),
        "C[C@H](F)Cl |o1:1|"
    );
    let skip_all = MolHashParams {
        use_cx_smiles: true,
        cx_flags_to_skip: 0x7fff_ffff,
    };
    assert_eq!(
        mol_hash_with_params(&mol, HashFunction::CanonicalSmiles, &skip_all).unwrap(),
        "C[C@H](F)Cl"
    );
}

#[test]
fn exclude_from_tautomerism() {
    assert_eq!(
        exclude_from_tautomerism_prop(),
        "_MolHashExcludeFromTautomerism"
    );
    let mut mol = mol_from_smiles("CC(=O)C=CN").unwrap();
    assert_eq!(
        mol_hash(&mol, HashFunction::HetAtomTautomerv2).unwrap(),
        "[C]:[C](:[O]):[C]:[C]:[N]_7_0"
    );
    for idx in 0..mol.get_num_atoms() {
        mol.get_atom_with_idx_mut(idx)
            .unwrap()
            .set_bool_prop(exclude_from_tautomerism_prop(), true);
    }
    // Python's values.
    assert_eq!(
        mol_hash(&mol, HashFunction::HetAtomTautomerv2).unwrap(),
        "[CH3]-[C](=[O])-[CH]=[CH]-[NH2]_0_0"
    );
    assert_eq!(
        mol_hash(&mol, HashFunction::HetAtomTautomer).unwrap(),
        "C[C]([O])[CH][CH][N]_2_0"
    );
}
