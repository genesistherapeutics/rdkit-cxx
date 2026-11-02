// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKFingerprint`'s and `UnfoldedRDKFingerprintCountBased`'s
    /// arguments; the optional lists are `None` when their `has_` is false.
    struct RDKFingerprintData {
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        n_bits_per_hash: u32,
        use_hs: bool,
        tgt_density: f64,
        min_size: u32,
        branched_paths: bool,
        use_bond_order: bool,
        has_atom_invariants: bool,
        atom_invariants: Vec<u32>,
        has_from_atoms: bool,
        from_atoms: Vec<u32>,
    }

    /// `LayeredFingerprint`'s arguments.
    struct LayeredFingerprintData {
        layer_flags: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        has_atom_counts: bool,
        atom_counts: Vec<u32>,
        has_set_only_bits: bool,
        branched_paths: bool,
        has_from_atoms: bool,
        from_atoms: Vec<u32>,
    }

    /// `PatternFingerprint`'s arguments.
    struct PatternFingerprintData {
        fp_size: u32,
        has_atom_counts: bool,
        atom_counts: Vec<u32>,
        has_set_only_bits: bool,
        tautomer_fingerprints: bool,
    }

    /// A list of bit ids.
    struct U64List {
        values: Vec<u64>,
    }

    /// What Python's output arguments receive: each atom's bits, each bit's
    /// paths (flattened, with a count per bit), and the updated atom counts.
    struct FingerprintOutputs {
        atom_bits: Vec<U64List>,
        bit_info_bits: Vec<u64>,
        bit_info_num_paths: Vec<u64>,
        bit_info_paths: Vec<IndexList>,
        atom_counts: Vec<u32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/fingerprints.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolBundle = crate::ffi::mol_bundle::ffi::MolBundle;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type SparseIntVect = crate::ffi::data_structs::ffi::SparseIntVect;
        type IndexList = crate::ffi::geometry::ffi::IndexList;

        fn rdk_fingerprint(
            mol: &RWMol,
            params: &RDKFingerprintData,
            with_outputs: bool,
            outputs: &mut FingerprintOutputs,
        ) -> Result<UniquePtr<BitVect>>;
        fn unfolded_rdk_fingerprint_count_based(
            mol: &RWMol,
            params: &RDKFingerprintData,
            with_outputs: bool,
            outputs: &mut FingerprintOutputs,
        ) -> Result<UniquePtr<SparseIntVect>>;
        fn layered_fingerprint(
            mol: &RWMol,
            params: &LayeredFingerprintData,
            set_only_bits: &BitVect,
            outputs: &mut FingerprintOutputs,
        ) -> Result<UniquePtr<BitVect>>;
        fn pattern_fingerprint(
            mol: &RWMol,
            params: &PatternFingerprintData,
            set_only_bits: &BitVect,
            outputs: &mut FingerprintOutputs,
        ) -> Result<UniquePtr<BitVect>>;
        fn pattern_fingerprint_mol_bundle(
            bundle: &MolBundle,
            params: &PatternFingerprintData,
            set_only_bits: &BitVect,
        ) -> Result<UniquePtr<BitVect>>;
        fn layered_fingerprint_substruct_layers() -> Result<u32>;
    }
}
