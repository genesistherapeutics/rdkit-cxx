// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ports testVectorDescriptors and testVectorDescriptorsInDescList of RDKit's
// rdkit/Chem/UnitTestDescriptors.py. setup_autocorr_descriptors changes the
// descriptor list for the whole process, as Python's does, so these run in a
// test binary of their own.

use std::collections::BTreeSet;

use rdkit_cxx::chem::descriptors::{calc_mol_descriptors, desc_list, setup_autocorr_descriptors};
use rdkit_cxx::chem::mol_from_smiles;
use rdkit_cxx::chem::rd_mol_descriptors::calc_autocorr_2d;

#[test]
fn autocorr_descriptors() {
    let before: Vec<&str> = desc_list().iter().map(|(name, _)| *name).collect();
    assert_eq!(before.len(), 210);
    setup_autocorr_descriptors();
    let after: Vec<&str> = desc_list().iter().map(|(name, _)| *name).collect();
    // As in Python: the AUTOCORR2D descriptors, descList, and then the
    // descriptors from GraphDescriptors on a second time.
    assert_eq!(after.len(), 586);
    assert_eq!(after[0], "AUTOCORR2D_1");
    assert_eq!(after[191], "AUTOCORR2D_192");
    assert_eq!(after[192..402], before[..]);
    assert_eq!(after[402], "AvgIpc");
    assert_eq!(after[402..], before[26..]);
    let names: BTreeSet<&str> = after.iter().copied().collect();
    assert!(before.iter().all(|name| names.contains(name)));
    assert_eq!(names.len(), 402);

    let m = mol_from_smiles("CCCc1ccccc1").unwrap();
    let results = calc_autocorr_2d(&m).unwrap();
    let list = desc_list();
    for (i, expected) in results.iter().enumerate() {
        let name = format!("AUTOCORR2D_{}", i + 1);
        let (_, f) = list.iter().find(|(n, _)| *n == name).unwrap();
        assert_eq!(f(&m).unwrap(), *expected);
    }
    let descs = calc_mol_descriptors(&m);
    // Python's dictionary keeps each name where it first appears.
    assert_eq!(descs.len(), 402);
    let names: Vec<&str> = descs.iter().map(|(name, _)| *name).collect();
    assert_eq!(names, after[..402]);
    assert_eq!(descs[0], ("AUTOCORR2D_1", Some(results[0])));
}
