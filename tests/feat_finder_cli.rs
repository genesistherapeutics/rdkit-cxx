// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Values from Python's rdkit.Chem.FeatFinderCLI in RDKit 2026.03.6, with
// RDKit's BaseFeatures.fdef.

use rdkit_cxx::chem::feat_finder_cli::get_atom_feat_info;
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_mol_chemical_features::build_feature_factory_from_string;

const BASE_FEATURES: &str = include_str!("../wrapper/third_party/rdkit/Data/BaseFeatures.fdef");

fn names(names: &[&str]) -> Option<Vec<String>> {
    Some(names.iter().map(|s| s.to_string()).collect())
}

#[test]
fn atom_feat_info_matches_python() {
    let factory = build_feature_factory_from_string(BASE_FEATURES).unwrap();
    let info =
        |smiles: &str| get_atom_feat_info(&factory, &mol_from_smiles(smiles).unwrap()).unwrap();
    let donor_acceptor = ["Donor-SingleAtomDonor", "Acceptor-SingleAtomAcceptor"];
    assert_eq!(
        info("OCC(=O)O"),
        [
            names(&donor_acceptor),
            names(&["ZnBinder-ZnBinder5"]),
            names(&["NegIonizable-AcidicGroup", "ZnBinder-ZnBinder5"]),
            names(&[
                "Acceptor-SingleAtomAcceptor",
                "NegIonizable-AcidicGroup",
                "ZnBinder-ZnBinder5"
            ]),
            names(&[
                "Donor-SingleAtomDonor",
                "Acceptor-SingleAtomAcceptor",
                "NegIonizable-AcidicGroup",
                "ZnBinder-ZnBinder5"
            ]),
        ]
    );
    let ring = names(&["Aromatic-Arom6", "LumpedHydrophobe-RH6_6"]);
    let mut expected = vec![ring; 6];
    expected.push(names(&["Donor-SingleAtomDonor"]));
    assert_eq!(info("c1ccccc1N"), expected);
    assert_eq!(info("CC"), [None, None]);
}
