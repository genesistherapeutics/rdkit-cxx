// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

fn smiles(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn unsanitized(smiles: &str) -> Mol {
    let params = SmilesParserParams {
        sanitize: false,
        ..Default::default()
    };
    mol_from_smiles_with_params(smiles, &params).unwrap()
}

#[test]
fn rings() {
    let mut cubane = smiles("C12C3C4C1C5C2C3C45");
    assert_eq!(get_sssr(&mut cubane).unwrap().len(), 5);
    let symm = get_symm_sssr(&mut cubane).unwrap();
    assert_eq!(symm.len(), 6);
    assert_eq!(symm[0], [0, 1, 2, 3]);
    let rdl = GetSymmSSSRParams {
        algorithm: SymmetrizeSSSRAlgorithm::Rdl,
        ..Default::default()
    };
    assert_eq!(
        get_symm_sssr_with_params(&mut cubane, &rdl).unwrap().len(),
        6
    );
    find_ring_families(&mut cubane).unwrap();
    assert_eq!(cubane.get_ring_info().num_ring_families().unwrap(), 6);

    let mut benzoic = smiles("c1ccccc1C(=O)O");
    assert_eq!(get_sssr(&mut benzoic).unwrap(), [[0, 1, 2, 3, 4, 5]]);
    let mut fast = unsanitized("C1CC1C");
    fast_find_rings(&mut fast).unwrap();
    assert_eq!(fast.get_ring_info().num_rings().unwrap(), 1);
    assert!(!get_use_legacy_ring_finding());
}

#[test]
fn param_defaults() {
    assert_eq!(
        GetSSSRParams::default(),
        GetSSSRParams {
            include_dative_bonds: false,
            include_hydrogen_bonds: false
        }
    );
    let symm = GetSymmSSSRParams::default();
    assert_eq!(symm.algorithm, SymmetrizeSSSRAlgorithm::Default);
    assert!(symm.recalc_sssr && !symm.include_dative_bonds && !symm.include_hydrogen_bonds);
    assert_eq!(
        FindRingFamiliesParams::default(),
        FindRingFamiliesParams {
            include_dative_bonds: false,
            include_hydrogen_bonds: false
        }
    );
    let kekulize = KekulizeIfPossibleParams::default();
    assert!(!kekulize.clear_aromatic_flags && kekulize.canonical);
}

#[test]
fn sanitization_steps() {
    let mut nitro = unsanitized("CN(=O)=O");
    cleanup(&mut nitro).unwrap();
    assert_eq!(mol_to_smiles(&nitro).unwrap(), "C[N+](=O)[O-]");

    let mut benzene = unsanitized("C1=CC=CC=C1");
    benzene.update_property_cache().unwrap();
    fast_find_rings(&mut benzene).unwrap();
    set_aromaticity(&mut benzene).unwrap();
    assert!(benzene.get_atom_with_idx(0).unwrap().get_is_aromatic());
    // A custom model needs a C++ function, as in Python.
    assert!(matches!(
        set_aromaticity_with_model(&mut benzene, AromaticityModel::AromaticityCustom),
        Err(Error::Runtime { .. })
    ));
    set_aromaticity_with_model(&mut benzene, AromaticityModel::AromaticityMdl).unwrap();

    let mut kekule = smiles("c1ccccc1");
    kekulize_if_possible(&mut kekule).unwrap();
    assert_eq!(
        kekule.get_bond_with_idx(0).unwrap().get_bond_type(),
        BondType::Single
    );

    let mut methyl = unsanitized("[CH3]");
    methyl.update_property_cache_with_strict(false).unwrap();
    assert_eq!(
        methyl
            .get_atom_with_idx(0)
            .unwrap()
            .get_num_radical_electrons(),
        0
    );
    assign_radicals(&mut methyl).unwrap();
    assert_eq!(
        methyl
            .get_atom_with_idx(0)
            .unwrap()
            .get_num_radical_electrons(),
        1
    );

    let mut butadiene = unsanitized("C=CC=C");
    butadiene.update_property_cache().unwrap();
    set_conjugation(&mut butadiene).unwrap();
    set_hybridization(&mut butadiene).unwrap();
    assert!(butadiene.get_bond_with_idx(1).unwrap().get_is_conjugated());
    assert_eq!(
        butadiene.get_atom_with_idx(1).unwrap().get_hybridization(),
        HybridizationType::Sp2
    );
    cleanup_chirality(&mut butadiene).unwrap();
    cleanup_atropisomers(&mut butadiene).unwrap();
    cleanup_organometallics(&mut butadiene).unwrap();
}

#[test]
fn queries() {
    let benzoic = smiles("c1ccccc1C(=O)O");
    assert_eq!(
        get_formal_charge(&smiles("[NH4+].[O-]C(=O)C.[O-]C")).unwrap(),
        -1
    );
    assert_eq!(get_shortest_path(&benzoic, 0, 8).unwrap(), [0, 5, 6, 8]);
    assert!(get_shortest_path(&smiles("C.C"), 0, 1).unwrap().is_empty());
    assert_eq!(
        get_shortest_path(&benzoic, 0, 20).unwrap_err(),
        Error::Index { index: 20, len: 9 }
    );
    assert!(needs_hs(&benzoic).unwrap());
    assert!(!needs_hs(&add_hs(&benzoic).unwrap()).unwrap());

    let atom = benzoic.get_atom_with_idx(0).unwrap();
    assert_eq!(count_atom_elec(&atom).unwrap(), 1);
    assert!(atom_has_conjugated_bond(&atom).unwrap());
    let ethane = smiles("CC");
    assert!(!atom_has_conjugated_bond(&ethane.get_atom_with_idx(0).unwrap()).unwrap());
    // These need the atom's molecule.
    assert!(matches!(
        count_atom_elec(&Atom::new("C").unwrap()),
        Err(Error::Runtime { .. })
    ));
}

#[test]
fn haptic_bonds() {
    // Ferrocene with each ring drawn as a dummy atom and a dative bond.
    let haptic = mol_from_mol_block(FERROCENE_HAPTIC).unwrap();
    assert_eq!((haptic.get_num_atoms(), haptic.get_num_bonds()), (13, 12));
    let dative = haptic_bonds_to_dative(&haptic).unwrap();
    assert_eq!((dative.get_num_atoms(), dative.get_num_bonds()), (11, 20));
    let back = dative_bonds_to_haptic(&dative).unwrap();
    assert_eq!((back.get_num_atoms(), back.get_num_bonds()), (13, 12));
}

#[test]
fn terminal_atom_coords() {
    let mut mol = mol_from_mol_block(CHLOROETHANE).unwrap();
    set_terminal_atom_coords(&mut mol, 2, 1).unwrap();
    let conf = mol.get_conformer().unwrap();
    let h = conf.get_atom_position(2).unwrap();
    let c = conf.get_atom_position(1).unwrap();
    assert!((h.distance(&c) - 1.1).abs() < 1e-6);
    // The atom must be terminal.
    assert!(matches!(
        set_terminal_atom_coords(&mut mol, 1, 0),
        Err(Error::Runtime { .. })
    ));
    assert!(set_terminal_atom_coords(&mut mol, 9, 0).is_err());
}

const CHLOROETHANE: &str = "
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    5.0000    5.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
";

const FERROCENE_HAPTIC: &str = "
  Mrv2014

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 12 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Fe 0 0 0 0
M  V30 2 * 0 1.5 0 0
M  V30 3 C -1 2 0 0
M  V30 4 C -0.5 3 0 0
M  V30 5 C 0.5 3 0 0
M  V30 6 C 1 2 0 0
M  V30 7 C 0 1.2 0 0 CHG=-1
M  V30 8 * 0 -1.5 0 0
M  V30 9 C -1 -2 0 0
M  V30 10 C -0.5 -3 0 0
M  V30 11 C 0.5 -3 0 0
M  V30 12 C 1 -2 0 0
M  V30 13 C 0 -1.2 0 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 3 4
M  V30 2 1 4 5
M  V30 3 2 5 6
M  V30 4 1 6 7
M  V30 5 1 7 3
M  V30 6 2 9 10
M  V30 7 1 10 11
M  V30 8 2 11 12
M  V30 9 1 12 13
M  V30 10 1 13 9
M  V30 11 9 2 1 ENDPTS=(5 3 4 5 6 7) ATTACH=ALL
M  V30 12 9 8 1 ENDPTS=(5 9 10 11 12 13) ATTACH=ALL
M  V30 END BOND
M  V30 END CTAB
M  END
";

#[test]
fn added_rings_are_cycles() {
    // RDKit's ring algorithms and pickle reader take a ring to be a cycle.
    let mut m = smiles("CCCC");
    let value = |r: rdkit_cxx::Result<()>| matches!(r, Err(Error::Value { .. }));
    let mut info = m.get_ring_info_mut();
    assert!(value(info.add_ring(&[], &[])));
    assert!(value(info.add_ring(&[0, 3], &[0, 1])));
    assert!(value(info.add_ring(&[0, 1, 1], &[0, 1, 2])));
    // Bond 1 joins atoms 1 and 2, not 1 and 3.
    assert!(value(info.add_ring(&[0, 1, 3], &[0, 1, 2])));
    assert!(matches!(
        info.add_ring(&[0, 1, 9], &[0, 1, 2]),
        Err(Error::Index { .. })
    ));
    // A real cycle, in either direction around the ring.
    let mut ring = smiles("C1CCC1");
    let mut info = ring.get_ring_info_mut();
    info.add_ring(&[0, 1, 2, 3], &[0, 1, 2, 3]).unwrap();
    info.add_ring(&[3, 2, 1, 0], &[2, 1, 0, 3]).unwrap();
    assert_eq!(ring.get_ring_info().num_rings().unwrap(), 3);
    // A molecule with an added ring pickles and reads back.
    let pkl = ring.to_binary();
    // SAFETY: a pickle RDKit wrote.
    assert_eq!(
        unsafe { Mol::from_pkl_string(&pkl) }
            .unwrap()
            .get_num_atoms(),
        4
    );
}
