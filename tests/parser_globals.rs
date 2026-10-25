// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The rdmolfiles functions that write RDKit's SMILES and SMARTS parser debug
//! levels (CONVENTIONS.md U8). This binary has one test, so no other thread
//! parses while they run. Expected values come from RDKit's Python wrappers.

use rdkit_cxx::Error;
use rdkit_cxx::chem::*;

const SVG: &str = "<svg><metadata><rdkit:mol xmlns:rdkit=\"http://www.rdkit.org/xml\" \
                   version=\"0.9\"><rdkit:atom idx=\"1\" atom-smiles=\"[CH3]\" x=\"0\" y=\"0\"/>\
                   <rdkit:atom idx=\"2\" atom-smiles=\"[OH]\" x=\"1\" y=\"0\"/>\
                   <rdkit:bond idx=\"1\" begin-atom-idx=\"1\" end-atom-idx=\"2\" \
                   bond-smiles=\"-\"/></rdkit:mol></metadata></svg>";

#[test]
fn atoms_bonds_and_svg() {
    // SAFETY: this is the only test in this binary, and nothing else in it
    // parses SMILES or SMARTS.
    unsafe {
        let atom = atom_from_smiles("[13CH3-:4]").unwrap();
        assert_eq!(
            (
                atom.get_atomic_num(),
                atom.get_isotope(),
                atom.get_formal_charge(),
                atom.get_atom_map_num()
            ),
            (6, 13, -1, 4)
        );
        assert!(matches!(atom_from_smiles("["), Err(Error::Parse { .. })));
        assert_eq!(
            atom_from_smarts("[C,N]").unwrap().get_smarts().unwrap(),
            "[C,N]"
        );
        assert_eq!(
            bond_from_smiles("=").unwrap().get_bond_type(),
            BondType::Double
        );
        assert!(matches!(bond_from_smiles("x"), Err(Error::Parse { .. })));
        assert_eq!(bond_from_smarts("~").unwrap().get_smarts().unwrap(), "~");

        let mol = mol_from_rdkit_svg(SVG).unwrap();
        assert_eq!(mol_to_smiles(&mol).unwrap(), "CO");
        let unsanitized = MolFromRDKitSVGParams {
            sanitize: false,
            remove_hs: false,
        };
        assert_eq!(
            mol_from_rdkit_svg_with_params(SVG, &unsanitized)
                .unwrap()
                .get_num_atoms(),
            2
        );
        // RDKit subtracts 1 from each bond's atom indices in int.
        let bad = SVG.replace("begin-atom-idx=\"1\"", "begin-atom-idx=\"-2147483648\"");
        assert!(matches!(mol_from_rdkit_svg(&bad), Err(Error::Value { .. })));
        assert!(matches!(
            mol_from_rdkit_svg("<svg/>"),
            Err(Error::Parse { .. })
        ));
    }
}
