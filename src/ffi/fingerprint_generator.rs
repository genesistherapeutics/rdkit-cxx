// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::FPType`; cxx checks the values.
    #[repr(i32)]
    enum FPType {
        AtomPairFP,
        MorganFP,
        RDKitFP,
        TopologicalTorsionFP,
    }

    /// Which `FingerprintGenerator` method a call stands for.
    #[repr(u8)]
    enum FingerprintOutput {
        Bits,
        SparseBits,
        Counts,
        SparseCounts,
    }

    /// `FingerprintArguments` and the fields of its subclass for `kind`.
    struct FingerprintOptionsData {
        kind: FPType,
        count_simulation: bool,
        include_chirality: bool,
        fp_size: u32,
        num_bits_per_feature: u32,
        count_bounds: Vec<u32>,
        radius: u32,
        only_nonzero_invariants: bool,
        include_redundant_environments: bool,
        min_distance: u32,
        max_distance: u32,
        use_2d: bool,
        min_path: u32,
        max_path: u32,
        use_hs: bool,
        branched_paths: bool,
        use_bond_order: bool,
        torsion_atom_count: u32,
        only_shortest_paths: bool,
    }

    /// `GetMorganGenerator`'s arguments but the invariant generators.
    struct MorganGeneratorData {
        radius: u32,
        count_simulation: bool,
        include_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
        include_ring_membership: bool,
        has_count_bounds: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
        include_redundant_environments: bool,
    }

    /// `GetAtomPairGenerator`'s arguments but the invariant generator.
    struct AtomPairGeneratorData {
        min_distance: u32,
        max_distance: u32,
        include_chirality: bool,
        use_2d: bool,
        count_simulation: bool,
        has_count_bounds: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
    }

    /// `GetRDKitFPGenerator`'s arguments but the invariant generator.
    struct RDKitFPGeneratorData {
        min_path: u32,
        max_path: u32,
        use_hs: bool,
        branched_paths: bool,
        use_bond_order: bool,
        count_simulation: bool,
        has_count_bounds: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
        num_bits_per_feature: u32,
    }

    /// `GetTopologicalTorsionGenerator`'s arguments but the invariant
    /// generator.
    struct TopologicalTorsionGeneratorData {
        include_chirality: bool,
        torsion_atom_count: u32,
        count_simulation: bool,
        has_count_bounds: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
    }

    /// The per-call arguments of `FingerprintGenerator.GetFingerprint` and
    /// its siblings.
    struct FingerprintArgsData {
        has_from_atoms: bool,
        from_atoms: Vec<u32>,
        has_ignore_atoms: bool,
        ignore_atoms: Vec<u32>,
        conf_id: i32,
        has_custom_atom_invariants: bool,
        custom_atom_invariants: Vec<u32>,
        has_custom_bond_invariants: bool,
        custom_bond_invariants: Vec<u32>,
    }

    /// `AdditionalOutput`: what to collect, and what was. Maps cross as keys,
    /// the number of values for each key, and the values in order.
    struct AdditionalOutputData {
        collect_atom_to_bits: bool,
        collect_bit_info_map: bool,
        collect_bit_paths: bool,
        collect_atom_counts: bool,
        collect_atoms_per_bit: bool,
        atom_to_bits: Vec<U64List>,
        bit_info_map_keys: Vec<u64>,
        bit_info_map_counts: Vec<u32>,
        bit_info_map_values: Vec<IndexPair>,
        bit_paths_keys: Vec<u64>,
        bit_paths_counts: Vec<u32>,
        bit_paths_values: Vec<IndexList>,
        atom_counts: Vec<u32>,
        atoms_per_bit_keys: Vec<u64>,
        atoms_per_bit_counts: Vec<u32>,
        atoms_per_bit_values: Vec<IndexList>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/fingerprint_generator.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");
        include!("rdkit-cxx/ffi/fingerprints.rs.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type SparseIntVect = crate::ffi::data_structs::ffi::SparseIntVect;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type U64List = crate::ffi::fingerprints::ffi::U64List;
        type FPType;
        type FpGenerator;
        type AtomInvGen;
        type BondInvGen;
        type FingerprintList;

        fn get_morgan_generator(
            data: &MorganGeneratorData,
            atom_inv_gen: UniquePtr<AtomInvGen>,
            bond_inv_gen: UniquePtr<BondInvGen>,
        ) -> Result<UniquePtr<FpGenerator>>;
        fn get_atom_pair_generator(
            data: &AtomPairGeneratorData,
            atom_inv_gen: UniquePtr<AtomInvGen>,
        ) -> Result<UniquePtr<FpGenerator>>;
        fn get_rdkit_fp_generator(
            data: &RDKitFPGeneratorData,
            atom_inv_gen: UniquePtr<AtomInvGen>,
        ) -> Result<UniquePtr<FpGenerator>>;
        fn get_topological_torsion_generator(
            data: &TopologicalTorsionGeneratorData,
            atom_inv_gen: UniquePtr<AtomInvGen>,
        ) -> Result<UniquePtr<FpGenerator>>;
        fn fingerprint_generator_from_json(json: &str) -> Result<UniquePtr<FpGenerator>>;

        fn get_morgan_atom_inv_gen(include_ring_membership: bool) -> Result<UniquePtr<AtomInvGen>>;
        fn get_morgan_feature_atom_inv_gen(
            has_patterns: bool,
            patterns: Pin<&mut BorrowedMols>,
        ) -> Result<UniquePtr<AtomInvGen>>;
        fn get_morgan_bond_inv_gen(
            use_bond_types: bool,
            use_chirality: bool,
        ) -> Result<UniquePtr<BondInvGen>>;
        fn get_atom_pair_atom_inv_gen(include_chirality: bool) -> Result<UniquePtr<AtomInvGen>>;
        fn get_rdkit_atom_inv_gen() -> Result<UniquePtr<AtomInvGen>>;
        fn atom_inv_gen_clone(generator: &AtomInvGen) -> Result<UniquePtr<AtomInvGen>>;
        fn bond_inv_gen_clone(generator: &BondInvGen) -> Result<UniquePtr<BondInvGen>>;
        fn atom_inv_gen_info_string(generator: &AtomInvGen) -> Result<UniquePtr<CxxString>>;
        fn bond_inv_gen_info_string(generator: &BondInvGen) -> Result<UniquePtr<CxxString>>;

        fn fp_generator_get_options(generator: &FpGenerator) -> Result<FingerprintOptionsData>;
        fn fp_generator_set_options(
            generator: Pin<&mut FpGenerator>,
            data: &FingerprintOptionsData,
        ) -> Result<()>;
        fn fp_generator_info_string(generator: &FpGenerator) -> Result<UniquePtr<CxxString>>;
        fn fp_generator_to_json(generator: &FpGenerator) -> Result<UniquePtr<CxxString>>;

        fn fp_generator_bits(
            generator: &FpGenerator,
            mol: &RWMol,
            output: FingerprintOutput,
            args: &FingerprintArgsData,
            outputs: &mut AdditionalOutputData,
        ) -> Result<UniquePtr<BitVect>>;
        fn fp_generator_counts(
            generator: &FpGenerator,
            mol: &RWMol,
            output: FingerprintOutput,
            args: &FingerprintArgsData,
            outputs: &mut AdditionalOutputData,
        ) -> Result<UniquePtr<SparseIntVect>>;
        fn fp_generator_values(
            generator: &FpGenerator,
            mol: &RWMol,
            output: FingerprintOutput,
            args: &FingerprintArgsData,
            outputs: &mut AdditionalOutputData,
        ) -> Result<Vec<u32>>;

        fn fp_generator_bulk(
            generator: &FpGenerator,
            mols: Pin<&mut BorrowedMols>,
            output: FingerprintOutput,
            num_threads: i32,
        ) -> Result<UniquePtr<FingerprintList>>;
        fn fp_type_bulk(
            mols: Pin<&mut BorrowedMols>,
            fp_type: FPType,
            output: FingerprintOutput,
        ) -> Result<UniquePtr<FingerprintList>>;
        fn fingerprint_list_size(list: &FingerprintList) -> Result<usize>;
        fn fingerprint_list_take_bits(
            list: Pin<&mut FingerprintList>,
            i: usize,
        ) -> Result<UniquePtr<BitVect>>;
        fn fingerprint_list_take_counts(
            list: Pin<&mut FingerprintList>,
            i: usize,
        ) -> Result<UniquePtr<SparseIntVect>>;
    }
}

// SAFETY: a generator owns its options and invariant generators and has no
// thread affinity. It is not Sync; only the crate's bulk calls share one
// across threads, as RDKit's own multithreaded calls do.
unsafe impl Send for ffi::FpGenerator {}
// SAFETY: an invariant generator owns its data (the feature generator copies
// its patterns).
unsafe impl Send for ffi::AtomInvGen {}
// SAFETY: as for AtomInvGen.
unsafe impl Send for ffi::BondInvGen {}
