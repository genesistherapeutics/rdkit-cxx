// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::FingerprintType`; cxx checks the values.
    #[repr(i32)]
    #[allow(clippy::upper_case_acronyms)]
    enum FingerprintType {
        AtomPairFP = 1,
        TopologicalTorsion,
        MorganFP,
        RDKitFP,
        PatternFP,
    }

    /// `RDKit::RxnOps::SanitizeRxnFlags`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum SanitizeRxnFlags {
        SANITIZE_NONE = 0x0,
        SANITIZE_RGROUP_NAMES = 0x1,
        SANITIZE_ATOM_MAPS = 0x2,
        SANITIZE_ADJUST_REACTANTS = 0x4,
        SANITIZE_MERGEHS = 0x8,
        SANITIZE_ADJUST_PRODUCTS = 0x10,
        SANITIZE_ALL = 0xFFFFFFFF,
    }

    /// `ReactionFingerprintParams`.
    struct ReactionFingerprintParamsData {
        include_agents: bool,
        bit_ratio_agents: f64,
        non_agent_weight: u32,
        agent_weight: i32,
        fp_size: u32,
        fp_type: FingerprintType,
    }

    /// One label `addRecursiveQueriesToReaction` reports: the reactant
    /// template, the atom, and the query's name.
    struct ReactantLabel {
        reactant: u32,
        atom_idx: u32,
        label: UniquePtr<CxxString>,
    }

    /// `PreprocessReaction`'s results.
    struct PreprocessReactionData {
        num_warnings: u32,
        num_errors: u32,
        num_reactants: u32,
        num_products: u32,
        /// The number of reactants labeled: zero when validation failed.
        num_label_lists: u32,
        labels: Vec<ReactantLabel>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/chem_reactions.h");
        include!("rdkit-cxx/ffi/geometry.rs.h");
        include!("rdkit-cxx/ffi/molfiles.rs.h");
        include!("rdkit-cxx/ffi/substruct.rs.h");
        include!("rdkit-cxx/ffi/adjust_query.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type RDProps = crate::ffi::props::ffi::RDProps;
        type MolList = crate::ffi::molops::ffi::MolList;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type SparseIntVect = crate::ffi::data_structs::ffi::SparseIntVect;
        type IndexList = crate::ffi::geometry::ffi::IndexList;
        type SmilesWriteParamsData = crate::ffi::molfiles::ffi::SmilesWriteParamsData;
        type SubstructMatchParametersData =
            crate::ffi::substruct::ffi::SubstructMatchParametersData;
        type AdjustQueryParametersData = crate::ffi::adjust_query::ffi::AdjustQueryParametersData;
        type ChemicalReaction;
        type ProductSets;
        type FingerprintType;
        type SanitizeRxnFlags;

        fn product_sets_size(sets: &ProductSets) -> Result<usize>;
        fn product_sets_set_size(sets: &ProductSets, i: usize) -> Result<usize>;
        fn product_sets_take(
            sets: Pin<&mut ProductSets>,
            i: usize,
            j: usize,
        ) -> Result<UniquePtr<RWMol>>;

        fn chemical_reaction_new() -> Result<UniquePtr<ChemicalReaction>>;
        fn chemical_reaction_clone(rxn: &ChemicalReaction) -> Result<UniquePtr<ChemicalReaction>>;
        fn chemical_reaction_from_bin_str(bin_str: &[u8]) -> Result<UniquePtr<ChemicalReaction>>;
        fn chemical_reaction_to_binary(
            rxn: &ChemicalReaction,
            has_flags: bool,
            flags: u32,
        ) -> Result<Vec<u8>>;
        fn chemical_reaction_as_props(rxn: &ChemicalReaction) -> Result<&RDProps>;
        fn chemical_reaction_get_num_templates(rxn: &ChemicalReaction, kind: u8) -> Result<u32>;
        fn chemical_reaction_get_template(
            rxn: &ChemicalReaction,
            kind: u8,
            which: u32,
        ) -> Result<&RWMol>;
        fn chemical_reaction_add_template(
            rxn: Pin<&mut ChemicalReaction>,
            kind: u8,
            mol: &RWMol,
        ) -> Result<u32>;
        fn chemical_reaction_remove_unmapped_templates(
            rxn: Pin<&mut ChemicalReaction>,
            kind: u8,
            threshold_unmapped_atoms: f64,
            move_to_agent_templates: bool,
        ) -> Result<UniquePtr<MolList>>;
        fn chemical_reaction_remove_agent_templates(
            rxn: Pin<&mut ChemicalReaction>,
        ) -> Result<UniquePtr<MolList>>;
        fn chemical_reaction_initialize(
            rxn: Pin<&mut ChemicalReaction>,
            silent: bool,
        ) -> Result<()>;
        fn chemical_reaction_is_initialized(rxn: &ChemicalReaction) -> Result<bool>;
        fn chemical_reaction_validate(
            rxn: Pin<&mut ChemicalReaction>,
            silent: bool,
        ) -> Result<Vec<u32>>;
        fn chemical_reaction_run_reactants(
            rxn: Pin<&mut ChemicalReaction>,
            reactants: Pin<&mut BorrowedMols>,
            max_products: u32,
        ) -> Result<UniquePtr<ProductSets>>;
        fn chemical_reaction_run_reactant(
            rxn: Pin<&mut ChemicalReaction>,
            reactant: &RWMol,
            reaction_idx: u32,
        ) -> Result<UniquePtr<ProductSets>>;
        fn chemical_reaction_run_reactant_in_place(
            rxn: Pin<&mut ChemicalReaction>,
            reactant: Pin<&mut RWMol>,
            remove_unmatched_atoms: bool,
        ) -> Result<bool>;
        fn chemical_reaction_is_molecule_template(
            rxn: &ChemicalReaction,
            kind: u8,
            mol: &RWMol,
        ) -> Result<bool>;
        fn chemical_reaction_get_reacting_atoms(
            rxn: &ChemicalReaction,
            mapped_atoms_only: bool,
        ) -> Result<Vec<IndexList>>;
        fn chemical_reaction_add_recursive_queries(
            rxn: Pin<&mut ChemicalReaction>,
            queries: Pin<&mut BorrowedMols>,
            prop_name: &str,
            get_labels: bool,
        ) -> Result<Vec<ReactantLabel>>;
        fn chemical_reaction_get_substruct_params(
            rxn: &ChemicalReaction,
        ) -> Result<SubstructMatchParametersData>;
        fn chemical_reaction_set_substruct_params(
            rxn: Pin<&mut ChemicalReaction>,
            params: &SubstructMatchParametersData,
        ) -> Result<()>;

        fn reaction_from_smarts(
            smarts: &str,
            replacement_keys: &Vec<String>,
            replacement_values: &Vec<String>,
            use_smiles: bool,
        ) -> Result<UniquePtr<ChemicalReaction>>;
        fn reaction_to_smarts(
            rxn: &ChemicalReaction,
            params: &SmilesWriteParamsData,
        ) -> Result<UniquePtr<CxxString>>;
        fn reaction_to_smiles(
            rxn: &ChemicalReaction,
            params: &SmilesWriteParamsData,
        ) -> Result<UniquePtr<CxxString>>;
        fn reaction_from_rxn_block(
            rxn_block: &str,
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<UniquePtr<ChemicalReaction>>;
        fn reaction_from_rxn_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
            strict_parsing: bool,
        ) -> Result<UniquePtr<ChemicalReaction>>;
        fn reaction_from_mrv_block(
            mrv: &str,
            sanitize: bool,
            remove_hs: bool,
        ) -> Result<UniquePtr<ChemicalReaction>>;
        fn reaction_from_mrv_file(
            path: &[u8],
            sanitize: bool,
            remove_hs: bool,
        ) -> Result<UniquePtr<ChemicalReaction>>;
        fn mrv_block_is_reaction(mrv: &str) -> Result<bool>;
        fn mrv_file_is_reaction(path: &[u8]) -> Result<bool>;
        fn reaction_to_rxn_block(
            rxn: &ChemicalReaction,
            separate_agents: bool,
            force_v3000: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn reaction_to_v3k_rxn_block(
            rxn: &ChemicalReaction,
            separate_agents: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn reaction_to_mrv_block(
            rxn: &ChemicalReaction,
            pretty_print: bool,
        ) -> Result<UniquePtr<CxxString>>;
        fn reaction_to_mrv_file(
            rxn: &ChemicalReaction,
            path: &[u8],
            pretty_print: bool,
        ) -> Result<()>;
        fn reaction_from_molecule(mol: &RWMol) -> Result<UniquePtr<ChemicalReaction>>;
        fn reaction_to_molecule(rxn: &ChemicalReaction) -> Result<UniquePtr<RWMol>>;

        fn reaction_from_png_string(data: &[u8]) -> Result<UniquePtr<ChemicalReaction>>;
        fn reaction_from_png_file(path: &[u8]) -> Result<UniquePtr<ChemicalReaction>>;
        fn reaction_metadata_to_png_string(
            rxn: &ChemicalReaction,
            png: &[u8],
            include_pkl: bool,
            include_smiles: bool,
            include_smarts: bool,
            include_rxn: bool,
        ) -> Result<Vec<u8>>;
        fn reaction_metadata_to_png_file(
            rxn: &ChemicalReaction,
            path: &[u8],
            include_pkl: bool,
            include_smiles: bool,
            include_smarts: bool,
            include_rxn: bool,
        ) -> Result<Vec<u8>>;

        #[allow(clippy::too_many_arguments)]
        fn compute_2d_coords_for_reaction(
            rxn: Pin<&mut ChemicalReaction>,
            spacing: f64,
            update_props: bool,
            canon_orient: bool,
            n_flips_per_sample: u32,
            n_sample: u32,
            sample_seed: i32,
            permute_deg4_nodes: bool,
            bond_length: f64,
        ) -> Result<()>;

        fn reaction_fingerprint_params_default() -> Result<ReactionFingerprintParamsData>;
        fn default_difference_fp_params() -> Result<ReactionFingerprintParamsData>;
        fn default_structural_fp_params() -> Result<ReactionFingerprintParamsData>;
        fn create_difference_fingerprint_for_reaction(
            rxn: &ChemicalReaction,
            params: &ReactionFingerprintParamsData,
        ) -> Result<UniquePtr<SparseIntVect>>;
        fn create_structural_fingerprint_for_reaction(
            rxn: &ChemicalReaction,
            params: &ReactionFingerprintParamsData,
        ) -> Result<UniquePtr<BitVect>>;

        fn is_reaction_template_molecule_agent(mol: &RWMol, agent_threshold: f64) -> Result<bool>;
        fn has_reaction_atom_mapping(rxn: &ChemicalReaction) -> Result<bool>;
        fn has_reaction_substruct_match(
            rxn: &ChemicalReaction,
            query_rxn: &ChemicalReaction,
            include_agents: bool,
        ) -> Result<bool>;
        fn has_template_substruct_match(
            rxn: &ChemicalReaction,
            query_rxn: &ChemicalReaction,
            kind: u8,
        ) -> Result<bool>;
        fn update_products_stereochemistry(rxn: Pin<&mut ChemicalReaction>) -> Result<()>;
        fn reduce_product_to_side_chains(
            product: &RWMol,
            add_dummy_atoms: bool,
        ) -> Result<UniquePtr<RWMol>>;
        fn remove_mapping_numbers_from_reactions(rxn: Pin<&mut ChemicalReaction>) -> Result<()>;
        fn preprocess_reaction(
            rxn: Pin<&mut ChemicalReaction>,
            queries: Pin<&mut BorrowedMols>,
            prop_name: &str,
        ) -> Result<PreprocessReactionData>;

        fn sanitize_rxn(
            rxn: Pin<&mut ChemicalReaction>,
            sanitize_ops: u32,
            params: &AdjustQueryParametersData,
            catch_errors: bool,
        ) -> Result<u32>;
        fn sanitize_rxn_as_mols(rxn: Pin<&mut ChemicalReaction>, sanitize_ops: u32) -> Result<()>;
        fn get_default_adjust_params() -> Result<AdjustQueryParametersData>;
        fn get_chem_draw_rxn_adjust_params() -> Result<AdjustQueryParametersData>;
        fn match_only_at_rgroups_adjust_params() -> Result<AdjustQueryParametersData>;
    }
}

// SAFETY: a ChemicalReaction owns its templates: the wrappers copy every
// molecule they add, normalize_templates gives every reaction RDKit makes
// templates no other reaction shares, and the template views borrow the
// reaction, so moving one to another thread moves all it owns. It is not
// Sync: RDKit writes through a const reaction (validation sets product
// template properties, the SMILES writers set computed properties on the
// templates, ReactionToMolecule sets atom roles), as it does through a const
// molecule (CONVENTIONS.md S1).
unsafe impl Send for ffi::ChemicalReaction {}
