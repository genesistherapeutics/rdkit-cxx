// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolHash::HashFunction`; cxx checks the values.
    #[repr(i32)]
    enum HashFunction {
        AnonymousGraph = 1,
        ElementGraph = 2,
        CanonicalSmiles = 3,
        MurckoScaffold = 4,
        ExtendedMurcko = 5,
        MolFormula = 6,
        AtomBondCounts = 7,
        DegreeVector = 8,
        Mesomer = 9,
        HetAtomTautomer = 10,
        HetAtomProtomer = 11,
        RedoxPair = 12,
        Regioisomer = 13,
        NetCharge = 14,
        SmallWorldIndexBR = 15,
        SmallWorldIndexBRL = 16,
        ArthorSubstructureOrder = 17,
        HetAtomTautomerv2 = 18,
        HetAtomProtomerv2 = 19,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol_hash.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type HashFunction;

        fn mol_hash(
            mol: &RWMol,
            func: HashFunction,
            use_cx_smiles: bool,
            cx_flags_to_skip: u32,
        ) -> Result<UniquePtr<CxxString>>;
    }
}
