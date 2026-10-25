// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::Descriptors::NumRotatableBondsOptions`; cxx checks the values.
    #[repr(i32)]
    enum NumRotatableBondsOptions {
        Default = -1,
        NonStrict = 0,
        Strict = 1,
        StrictLinkages = 2,
    }

    /// `AtomPairs.h`'s constants: Python's `AtomPairsParameters`.
    struct AtomPairsConstants {
        num_type_bits: u32,
        num_pi_bits: u32,
        num_branch_bits: u32,
        num_chiral_bits: u32,
        code_size: u32,
        num_path_bits: u32,
        num_atom_pair_fingerprint_bits: u32,
        atom_types: Vec<u32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/descriptors.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Atom = crate::ffi::atom::ffi::Atom;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type NumRotatableBondsOptions;

        fn calc_crippen_descriptors(mol: &RWMol, include_hs: bool) -> Result<Vec<f64>>;
        fn calc_labute_asa(mol: &RWMol, include_hs: bool) -> Result<f64>;
        fn calc_tpsa(mol: &RWMol, include_s_and_p: bool) -> Result<f64>;
        fn calc_exact_mol_wt(mol: &RWMol, only_heavy: bool) -> Result<f64>;
        fn calc_mol_formula(
            mol: &RWMol,
            separate_isotopes: bool,
            abbreviate_h_isotopes: bool,
        ) -> Result<UniquePtr<CxxString>>;

        fn calc_num_lipinski_hbd(mol: &RWMol) -> Result<u32>;
        fn calc_num_lipinski_hba(mol: &RWMol) -> Result<u32>;
        fn calc_num_hbd(mol: &RWMol) -> Result<u32>;
        fn calc_num_hba(mol: &RWMol) -> Result<u32>;
        fn calc_num_rotatable_bonds(mol: &RWMol, strict: NumRotatableBondsOptions) -> Result<u32>;
        fn calc_num_rings(mol: &RWMol) -> Result<u32>;
        fn calc_num_aromatic_rings(mol: &RWMol) -> Result<u32>;
        fn calc_num_saturated_rings(mol: &RWMol) -> Result<u32>;
        fn calc_num_heterocycles(mol: &RWMol) -> Result<u32>;
        fn calc_num_aromatic_heterocycles(mol: &RWMol) -> Result<u32>;
        fn calc_num_aromatic_carbocycles(mol: &RWMol) -> Result<u32>;
        fn calc_num_saturated_heterocycles(mol: &RWMol) -> Result<u32>;
        fn calc_num_saturated_carbocycles(mol: &RWMol) -> Result<u32>;
        fn calc_num_aliphatic_rings(mol: &RWMol) -> Result<u32>;
        fn calc_num_aliphatic_heterocycles(mol: &RWMol) -> Result<u32>;
        fn calc_num_aliphatic_carbocycles(mol: &RWMol) -> Result<u32>;
        fn calc_num_heavy_atoms(mol: &RWMol) -> Result<u32>;
        fn calc_num_atoms(mol: &RWMol) -> Result<u32>;
        fn calc_num_heteroatoms(mol: &RWMol) -> Result<u32>;
        fn calc_num_amide_bonds(mol: &RWMol) -> Result<u32>;
        fn calc_fraction_csp3(mol: &RWMol) -> Result<f64>;
        fn calc_num_atom_stereo_centers(mol: &RWMol) -> Result<u32>;
        fn calc_num_unspecified_atom_stereo_centers(mol: &RWMol) -> Result<u32>;
        fn calc_num_spiro_atoms(mol: &RWMol, atoms: &mut Vec<u32>) -> Result<u32>;
        fn calc_num_bridgehead_atoms(mol: &RWMol, atoms: &mut Vec<u32>) -> Result<u32>;

        fn slogp_vsa(mol: &RWMol, bins: &[f64]) -> Result<Vec<f64>>;
        fn smr_vsa(mol: &RWMol, bins: &[f64]) -> Result<Vec<f64>>;
        fn peoe_vsa(mol: &RWMol, bins: &[f64]) -> Result<Vec<f64>>;
        fn custom_prop_vsa(mol: &RWMol, custom_prop_name: &str, bins: &[f64]) -> Result<Vec<f64>>;
        fn mqns(mol: &RWMol) -> Result<Vec<u32>>;

        fn calc_chi0v(mol: &RWMol) -> Result<f64>;
        fn calc_chi1v(mol: &RWMol) -> Result<f64>;
        fn calc_chi2v(mol: &RWMol) -> Result<f64>;
        fn calc_chi3v(mol: &RWMol) -> Result<f64>;
        fn calc_chi4v(mol: &RWMol) -> Result<f64>;
        fn calc_chi_nv(mol: &RWMol, n: u32) -> Result<f64>;
        fn calc_chi0n(mol: &RWMol) -> Result<f64>;
        fn calc_chi1n(mol: &RWMol) -> Result<f64>;
        fn calc_chi2n(mol: &RWMol) -> Result<f64>;
        fn calc_chi3n(mol: &RWMol) -> Result<f64>;
        fn calc_chi4n(mol: &RWMol) -> Result<f64>;
        fn calc_chi_nn(mol: &RWMol, n: u32) -> Result<f64>;
        fn calc_hall_kier_alpha(mol: &RWMol, atom_contribs: &mut Vec<f64>) -> Result<f64>;
        fn calc_kappa1(mol: &RWMol) -> Result<f64>;
        fn calc_kappa2(mol: &RWMol) -> Result<f64>;
        fn calc_kappa3(mol: &RWMol) -> Result<f64>;
        fn calc_phi(mol: &RWMol) -> Result<f64>;

        fn get_atom_features(mol: &RWMol, atomid: u32, addchiral: bool) -> Result<Vec<f64>>;
        fn calc_oxidation_numbers(mol: Pin<&mut RWMol>) -> Result<()>;
        fn calc_autocorr_2d(mol: &RWMol, custom_atom_property: &str) -> Result<Vec<f64>>;
        fn bcut_2d(mol: &RWMol) -> Result<Vec<f64>>;
        fn bcut_2d_atom_props(mol: &RWMol, atom_props: &[f64]) -> Result<Vec<f64>>;
        fn bcut_2d_atom_propname(mol: &RWMol, atom_propname: &str) -> Result<Vec<f64>>;

        fn get_maccs_keys_fingerprint(mol: &RWMol) -> Result<UniquePtr<BitVect>>;
        fn get_atom_pair_atom_code(
            atom: &Atom,
            branch_subtract: u32,
            include_chirality: bool,
        ) -> Result<u32>;
        fn get_atom_pair_code(
            atom1_code: u32,
            atom2_code: u32,
            distance: u32,
            include_chirality: bool,
        ) -> Result<u32>;
        fn get_connectivity_invariants(
            mol: &RWMol,
            include_ring_membership: bool,
        ) -> Result<Vec<u32>>;
        fn get_feature_invariants(mol: &RWMol) -> Result<Vec<u32>>;
        fn atom_pairs_constants() -> Result<AtomPairsConstants>;
        fn atom_pairs_version() -> Result<UniquePtr<CxxString>>;
    }
}
