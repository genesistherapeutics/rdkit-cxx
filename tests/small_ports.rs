// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The small Python modules: MolStandardize's ReorderTautomers,
// ChemicalFeatures, Graphs, and FilterCatalog. Values from Python's modules
// in RDKit 2025.09.6.

use rdkit_cxx::chem::chemical_features::{
    GetMolFeatureParams, build_feature_factory_from_string, mcff_get_features_for_mol,
    mcff_get_features_for_mol_with_params,
};
use rdkit_cxx::chem::filter_catalog::{CustomFilterMatcher, FilterMatch, FilterMatcher};
use rdkit_cxx::chem::graphs::characteristic_polynomial;
use rdkit_cxx::chem::mol_standardize::reorder_tautomers;
use rdkit_cxx::chem::rd_depictor::compute_2d_coords;
use rdkit_cxx::chem::{
    Conformer, Mol, get_adjacency_matrix, get_adjacency_matrix_use_bo, mol_from_smiles,
    mol_to_smiles,
};
use rdkit_cxx::geometry::rd_geometry::Point3D;
use rdkit_cxx::{Error, Result};

#[test]
fn reordered_tautomers() -> Result<()> {
    for (smiles, expected) in [
        (
            "Oc1c(cccc3)c3nc2ccncc12",
            &[
                "O=c1c2ccccc2[nH]c2ccncc12",
                "O=c1c2c[nH]ccc-2nc2ccccc12",
                "Oc1c2ccccc2nc2ccncc12",
            ][..],
        ),
        (
            "OC(C)=CC(=O)C",
            &[
                "CC(=O)CC(C)=O",
                "C=C(O)C=C(C)O",
                "C=C(O)CC(=C)O",
                "C=C(O)CC(C)=O",
                "CC(=O)C=C(C)O",
            ],
        ),
        ("c1ccccc1", &["c1ccccc1"]),
    ] {
        let tauts = reorder_tautomers(&mol_from_smiles(smiles)?)?;
        let smis: Vec<String> = tauts.iter().map(mol_to_smiles).collect::<Result<_>>()?;
        assert_eq!(smis, expected, "{smiles}");
    }
    Ok(())
}

#[test]
fn characteristic_polynomials() -> Result<()> {
    let m = mol_from_smiles("C1CCC1C")?;
    let adj: Vec<Vec<f64>> = get_adjacency_matrix(&m)?
        .into_iter()
        .map(|row| row.into_iter().map(f64::from).collect())
        .collect();
    assert_eq!(
        characteristic_polynomial(&m, &adj)?,
        [1.0, -0.0, -5.0, -0.0, 2.0, -0.0]
    );
    let m = mol_from_smiles("c1ccccc1O")?;
    assert_eq!(
        characteristic_polynomial(&m, &get_adjacency_matrix_use_bo(&m)?)?,
        [1.0, -0.0, -14.5, -0.0, 54.5625, -0.0, -60.75, -0.0]
    );
    assert!(matches!(
        characteristic_polynomial(&m, &adj),
        Err(Error::Value { .. })
    ));
    Ok(())
}

const FDEF: &str = "DefineFeature HDonor1 [N,O;!H0]
  Family HBondDonor
  Weights 1.0
EndFeature
DefineFeature HAcceptor1 [N,O]
  Family HBondAcceptor
  Weights 1.0
EndFeature
DefineFeature Arom6 a1aaaaa1
  Family Aromatic
  Weights 1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
";

#[test]
fn features_for_mol() -> Result<()> {
    let factory = build_feature_factory_from_string(FDEF)?;
    let mut m = mol_from_smiles("OCCc1ccccc1N")?;
    compute_2d_coords(&mut m)?;
    let mut conf = Conformer::new(m.get_num_atoms());
    for i in 0..m.get_num_atoms() {
        conf.set_atom_position(i, Point3D::new(f64::from(i), 2.0 * f64::from(i), 0.5))?;
    }
    conf.set_id(7)?;
    m.add_conformer(&conf)?;
    let feats: Vec<(String, String, Vec<u32>)> = mcff_get_features_for_mol(&factory, &m)?
        .iter()
        .map(|f| (f.get_family(), f.get_type(), f.get_atom_ids()))
        .collect();
    let expected = [
        ("HBondDonor", "HDonor1", vec![0]),
        ("HBondDonor", "HDonor1", vec![9]),
        ("HBondAcceptor", "HAcceptor1", vec![0]),
        ("HBondAcceptor", "HAcceptor1", vec![9]),
        ("Aromatic", "Arom6", vec![3, 4, 5, 6, 7, 8]),
    ]
    .map(|(f, t, a)| (f.to_owned(), t.to_owned(), a));
    assert_eq!(feats, expected);
    let params = GetMolFeatureParams {
        include_only: "HBondAcceptor".into(),
        conf_id: 7,
    };
    let feats = mcff_get_features_for_mol_with_params(&factory, &m, &params)?;
    assert_eq!(feats.len(), 2);
    assert_eq!(feats[1].get_atom_ids(), [9]);
    let p = feats[1].get_pos()?;
    assert_eq!((p.x, p.y, p.z), (9.0, 18.0, 0.5));
    let params = GetMolFeatureParams {
        include_only: "Aromatic".into(),
        conf_id: 7,
    };
    let feats = mcff_get_features_for_mol_with_params(&factory, &m, &params)?;
    let p = feats[0].get_pos()?;
    assert!((p.x - 5.5).abs() < 1e-12 && (p.y - 11.0).abs() < 1e-12 && (p.z - 0.5).abs() < 1e-12);
    Ok(())
}

/// The module docstring's example: a filter for molecular weights outside
/// a range, here counting heavy atoms instead.
struct HeavyAtomFilter {
    max: u32,
}

impl CustomFilterMatcher for HeavyAtomFilter {
    fn has_match(&self, mol: &Mol) -> bool {
        mol.get_num_heavy_atoms() > self.max
    }

    fn get_matches(&self, mol: &Mol, _matches: &mut Vec<FilterMatch>) -> bool {
        self.has_match(mol)
    }
}

#[test]
fn custom_filter_matchers() -> Result<()> {
    let matcher = FilterMatcher::new("too big", HeavyAtomFilter { max: 3 });
    assert_eq!(matcher.get_name()?, "too big");
    assert!(matcher.is_valid());
    assert!(matcher.has_match(&mol_from_smiles("CCCCC")?)?);
    assert!(!matcher.has_match(&mol_from_smiles("CC")?)?);
    assert!(matcher.get_matches(&mol_from_smiles("CC")?)?.is_empty());
    Ok(())
}

// TemplateAlign's doctests.
#[test]
fn templates_align_in_2d() -> Result<()> {
    use rdkit_cxx::chem::template_align::*;
    let mut patt = mol_from_smiles("C1CC1")?;
    compute_2d_coords(&mut patt)?;
    let mut mol = mol_from_smiles("OC1CC1CC1CCC1")?;
    compute_2d_coords(&mut mol)?;
    // Whether each template atom's x is where the molecule's atom `offset`
    // further on is.
    let aligned = |mol: &Mol, offset: u32| -> Result<Vec<bool>> {
        let (pc, mc) = (patt.get_conformer()?, mol.get_conformer()?);
        (0..3)
            .map(|i| {
                Ok((pc.get_atom_position(i)?.x - mc.get_atom_position(i + offset)?.x).abs() < 1e-4)
            })
            .collect()
    };
    assert_eq!(aligned(&mol, 1)?, [false, false, false]);
    let params = AlignMolToTemplate2DParams {
        clear_confs: true,
        ..Default::default()
    };
    assert_eq!(
        align_mol_to_template_2d_with_params(&mut mol, &patt, &params)?,
        0
    );
    assert_eq!(mol.get_num_conformers(), 1);
    assert_eq!(aligned(&mol, 1)?, [true, true, true]);
    let params = AlignMolToTemplate2DParams {
        match_: Some(vec![5, 6, 7]),
        clear_confs: true,
        ..Default::default()
    };
    assert_eq!(
        align_mol_to_template_2d_with_params(&mut mol, &patt, &params)?,
        0
    );
    assert_eq!(aligned(&mol, 1)?, [false, false, false]);
    assert_eq!(aligned(&mol, 5)?, [true, true, true]);
    assert_eq!(align_mol_to_template_2d(&mut mol, &patt)?, 1);
    assert_eq!(mol.get_num_conformers(), 2);
    match align_mol_to_template_2d(&mut mol_from_smiles("CCC")?, &patt) {
        Err(Error::Value { message }) => assert_eq!(message, "no match between mol and template"),
        other => panic!("{other:?}"),
    }
    Ok(())
}

// rdkit/Chem/UnitTestSATIS.py, and other neighbor counts.
#[test]
fn satis_codes() -> Result<()> {
    use rdkit_cxx::chem::satis::*;
    for (smiles, expected) in [
        (
            "CC(=O)NC",
            &[
                "0601010106",
                "0606070895",
                "0806999995",
                "0701060699",
                "0601010107",
            ][..],
        ),
        (
            "O=CC(=O)C",
            &[
                "0806999993",
                "0601060893",
                "0606060894",
                "0806999994",
                "0601010106",
            ],
        ),
        (
            "C(=O)OCC(=O)O",
            &[
                "0601080896",
                "0806999996",
                "0806069996",
                "0601010608",
                "0606080898",
                "0806999998",
                "0801069998",
            ],
        ),
        ("C(=O)[O-]", &["0601080897", "0806999997", "0806999997"]),
        (
            "CP(F)(Cl)(Br)(O)",
            &[
                "0601010115",
                "1508091735",
                "0915999999",
                "1715999999",
                "3515999999",
                "0801159999",
            ],
        ),
        ("C=C", &["0601010699", "0601010699"]),
        (
            "c1ccccc1[N+](=O)[O-]",
            &[
                "0601060699",
                "0601060699",
                "0601060699",
                "0601060699",
                "0601060699",
                "0606060799",
                "0706080899",
                "0807999999",
                "0807999999",
            ],
        ),
    ] {
        assert_eq!(
            satis_types(&mol_from_smiles(smiles)?)?,
            expected,
            "{smiles}"
        );
    }
    assert_eq!(
        satis_types_with_neighbors_to_include(&mol_from_smiles("CC(=O)NC")?, 2)?,
        ["060106", "060708", "080695", "070606", "060107"]
    );
    assert!(matches!(
        satis_types_with_neighbors_to_include(&mol_from_smiles("CC")?, 0),
        Err(Error::Index { .. })
    ));
    assert_eq!(
        satis_types_with_neighbors_to_include(&mol_from_smiles("[Na+]")?, 0)?,
        ["11"]
    );
    assert_eq!(special_cases().len(), 6);
    Ok(())
}

// FragmentMatcher's doctests.
#[test]
fn fragment_matchers() -> Result<()> {
    use rdkit_cxx::chem::BondType;
    use rdkit_cxx::chem::fragment_matcher::FragmentMatcher;
    let empty = FragmentMatcher::default();
    assert!(!empty.has_match(&mol_from_smiles("CC")?)?);
    assert_eq!(empty.get_match(&mol_from_smiles("CC")?)?, None);
    assert!(empty.get_bond(0)?.is_none());
    let mut p = FragmentMatcher::default();
    p.init("C=O")?;
    assert!(p.has_match(&mol_from_smiles("CC(=O)O")?)?);
    assert!(!p.has_match(&mol_from_smiles("CC(C)C")?)?);
    assert_eq!(p.get_matches(&mol_from_smiles("CC=O")?)?.unwrap().len(), 1);
    assert_eq!(
        p.get_matches(&mol_from_smiles("O=CC=O")?)?.unwrap().len(),
        2
    );
    p.add_exclusion("c1ccccc1")?;
    assert!(p.has_match(&mol_from_smiles("CC=O")?)?);
    assert!(!p.has_match(&mol_from_smiles("c1ccccc1CC=O")?)?);
    assert_eq!(
        p.get_match(&mol_from_smiles("c1ccccc1CC=O")?)?,
        Some(vec![7, 8])
    );
    assert_eq!(
        p.get_matches_with_uniquify(&mol_from_smiles("OC(=O)C=O")?, false)?,
        Some(vec![vec![1, 2], vec![3, 4]])
    );
    assert_eq!(p.get_match(&mol_from_smiles("CC")?)?, Some(vec![]));
    assert_eq!(p.get_bond(0)?.unwrap().get_bond_type(), BondType::Double);
    assert!(matches!(p.add_exclusion("C1"), Err(Error::Parse { .. })));
    Ok(())
}

// Randomize, with Python's generator seeded as `random.seed(42)` and
// `random.seed(5)` seed it.
#[test]
fn randomized_atom_orders() -> Result<()> {
    use rdkit_cxx::chem::mol_to_mol_block;
    use rdkit_cxx::chem::randomize::{randomize_mol, randomize_mol_block};
    use rdkit_cxx::rd_random::Random;
    let bonds = |m: &Mol| -> Vec<(u32, u32)> {
        m.get_bonds()
            .map(|b| (b.get_begin_atom_idx(), b.get_end_atom_idx()))
            .collect()
    };
    let mut random = Random::new(42);
    for (smiles, expected) in [
        ("CC(=O)[O-]", &[(0, 2), (2, 3), (2, 1)][..]),
        ("[NH3+]CC(=O)[O-]", &[(3, 0), (0, 2), (2, 4), (2, 1)]),
        ("[O-][O-]", &[(1, 0)]),
        (
            "c1ccccc1O",
            &[(0, 2), (2, 4), (4, 6), (6, 3), (3, 5), (5, 1), (5, 0)],
        ),
    ] {
        let m = randomize_mol(&mol_from_smiles(smiles)?, &mut random)?;
        assert_eq!(bonds(&m), expected, "{smiles}");
    }
    assert_eq!(
        randomize_mol(&mol_from_smiles("C")?, &mut random)?.get_num_atoms(),
        1
    );
    let block = mol_to_mol_block(&mol_from_smiles("OCC=O")?)?;
    assert_eq!(
        randomize_mol_block(&block, &mut Random::new(5))?,
        "
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  4  3  1  0
  3  1  1  0
  1  2  2  0
M  END
"
    );
    Ok(())
}
