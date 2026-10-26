// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's 2D descriptors: counts, surface areas, partition
//! coefficients, connectivity indices, and BCUT and autocorrelation vectors.

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::error::{expect_ok, to_string};
use crate::ffi::descriptors::ffi;

rdkit_enum! {
    /// The rotatable-bond definitions. Mirrors
    /// `rdMolDescriptors.NumRotatableBondsOptions`.
    NumRotatableBondsOptions: ffi::NumRotatableBondsOptions {
        NonStrict = NonStrict,
        Strict = Strict,
        StrictLinkages = StrictLinkages,
        Default = Default,
    }
}

/// Python's `bool` overload: `true` is `Strict` and `false` `NonStrict`, as
/// in RDKit.
impl From<bool> for NumRotatableBondsOptions {
    fn from(strict: bool) -> NumRotatableBondsOptions {
        if strict {
            NumRotatableBondsOptions::Strict
        } else {
            NumRotatableBondsOptions::NonStrict
        }
    }
}

fn pair(values: Vec<f64>) -> (f64, f64) {
    (values[0], values[1])
}

/// The Wildman-Crippen logP and MR, with hydrogens. Mirrors
/// `rdMolDescriptors.CalcCrippenDescriptors`. Python's `force`, which only
/// steers RDKit's descriptor caches, is not bound: the crate never keeps
/// them (CONVENTIONS.md U8), so every result is computed afresh.
pub fn calc_crippen_descriptors(mol: &Mol) -> Result<(f64, f64)> {
    calc_crippen_descriptors_with_include_hs(mol, true)
}

/// Mirrors `rdMolDescriptors.CalcCrippenDescriptors(mol, includeHs)`.
pub fn calc_crippen_descriptors_with_include_hs(mol: &Mol, include_hs: bool) -> Result<(f64, f64)> {
    Ok(pair(ffi::calc_crippen_descriptors(
        mol.as_ffi(),
        include_hs,
    )?))
}

/// Labute's approximate surface area, with hydrogens. Mirrors
/// `rdMolDescriptors.CalcLabuteASA`; `force` is not bound (U8).
pub fn calc_labute_asa(mol: &Mol) -> Result<f64> {
    calc_labute_asa_with_include_hs(mol, true)
}

/// Mirrors `rdMolDescriptors.CalcLabuteASA(mol, includeHs)`.
pub fn calc_labute_asa_with_include_hs(mol: &Mol, include_hs: bool) -> Result<f64> {
    Ok(ffi::calc_labute_asa(mol.as_ffi(), include_hs)?)
}

/// The topological polar surface area from nitrogens and oxygens. Mirrors
/// `rdMolDescriptors.CalcTPSA`; `force` is not bound (U8).
pub fn calc_tpsa(mol: &Mol) -> Result<f64> {
    calc_tpsa_with_include_s_and_p(mol, false)
}

/// As `calc_tpsa`, with sulfurs and phosphorus when `include_s_and_p`.
/// Mirrors `rdMolDescriptors.CalcTPSA(mol, includeSandP=...)`.
pub fn calc_tpsa_with_include_s_and_p(mol: &Mol, include_s_and_p: bool) -> Result<f64> {
    Ok(ffi::calc_tpsa(mol.as_ffi(), include_s_and_p)?)
}

/// The monoisotopic molecular weight. Mirrors
/// `rdMolDescriptors.CalcExactMolWt`.
pub fn calc_exact_mol_wt(mol: &Mol) -> Result<f64> {
    calc_exact_mol_wt_with_only_heavy(mol, false)
}

/// Without hydrogens when `only_heavy`. Mirrors
/// `rdMolDescriptors.CalcExactMolWt(mol, onlyHeavy)`.
pub fn calc_exact_mol_wt_with_only_heavy(mol: &Mol, only_heavy: bool) -> Result<f64> {
    Ok(ffi::calc_exact_mol_wt(mol.as_ffi(), only_heavy)?)
}

/// Parameters for `calc_mol_formula_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct CalcMolFormulaParams {
    /// Writes each isotope separately (`C[13C]H6`).
    pub separate_isotopes: bool,
    /// With `separate_isotopes`, writes deuterium and tritium as D and T.
    pub abbreviate_h_isotopes: bool,
}

impl Default for CalcMolFormulaParams {
    fn default() -> CalcMolFormulaParams {
        CalcMolFormulaParams {
            separate_isotopes: false,
            abbreviate_h_isotopes: true,
        }
    }
}

/// The molecular formula, in Hill order. RDKit sums the formal charges in
/// `int` and negates it; a molecule whose charges could overflow it (only
/// possible with 2^24 or more atoms) is `Error::Value`. Mirrors
/// `rdMolDescriptors.CalcMolFormula`.
pub fn calc_mol_formula(mol: &Mol) -> Result<String> {
    calc_mol_formula_with_params(mol, &CalcMolFormulaParams::default())
}

/// Mirrors `rdMolDescriptors.CalcMolFormula(mol, separateIsotopes,
/// abbreviateHIsotopes)`.
pub fn calc_mol_formula_with_params(mol: &Mol, params: &CalcMolFormulaParams) -> Result<String> {
    to_string(ffi::calc_mol_formula(
        mol.as_ffi(),
        params.separate_isotopes,
        params.abbreviate_h_isotopes,
    )?)
}

/// The average molecular weight, of the heavy atoms alone when
/// `only_heavy`. Python's private `rdMolDescriptors._CalcMolWt`, behind
/// `Descriptors.MolWt`.
pub(crate) fn calc_mol_wt(mol: &Mol, only_heavy: bool) -> Result<f64> {
    Ok(ffi::calc_amw(mol.as_ffi(), only_heavy)?)
}

/// Each atom's contribution to Labute's surface area, after the implicit
/// hydrogens' (which is 0 unless `include_hs`). Python's private
/// `rdMolDescriptors._CalcLabuteASAContribs`, as `MolSurf._LabuteHelper`
/// arranges it.
pub(crate) fn calc_labute_asa_contribs(mol: &Mol, include_hs: bool) -> Result<Vec<f64>> {
    Ok(ffi::calc_labute_asa_contribs(mol.as_ffi(), include_hs)?)
}

/// Each atom's Wildman-Crippen (logP, MR) contributions. Python's private
/// `rdMolDescriptors._CalcCrippenContribs`.
pub(crate) fn calc_crippen_contribs(mol: &Mol) -> Result<Vec<(f64, f64)>> {
    Ok(ffi::calc_crippen_contribs(mol.as_ffi())?
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect())
}

/// Defines functions of a molecule that return a count and can fail.
macro_rules! fallible_counts {
    ($($(#[$meta:meta])* $name:ident;)*) => {
        $(
            $(#[$meta])*
            pub fn $name(mol: &Mol) -> Result<u32> {
                Ok(ffi::$name(mol.as_ffi())?)
            }
        )*
    };
}

/// Defines functions of a molecule that return a count and cannot fail:
/// they read only atomic numbers and the ring info's lists, which are empty
/// before rings are found.
macro_rules! counts {
    ($($(#[$meta:meta])* $name:ident;)*) => {
        $(
            $(#[$meta])*
            pub fn $name(mol: &Mol) -> u32 {
                expect_ok(ffi::$name(mol.as_ffi()))
            }
        )*
    };
}

fallible_counts! {
    /// The hydrogens on nitrogens and oxygens. Mirrors
    /// `rdMolDescriptors.CalcNumLipinskiHBD`.
    calc_num_lipinski_hbd;
    /// The hydrogen-bond donors. Mirrors `rdMolDescriptors.CalcNumHBD`.
    calc_num_hbd;
    /// The hydrogen-bond acceptors. Mirrors `rdMolDescriptors.CalcNumHBA`.
    calc_num_hba;
    /// The rings in the ring info; `Error::Runtime` before rings are found.
    /// Mirrors `rdMolDescriptors.CalcNumRings`.
    calc_num_rings;
    /// The atoms, with hydrogens. Mirrors `rdMolDescriptors.CalcNumAtoms`.
    calc_num_atoms;
    /// The atoms other than carbon and hydrogen. Mirrors
    /// `rdMolDescriptors.CalcNumHeteroatoms`.
    calc_num_heteroatoms;
    /// The amide bonds. Mirrors `rdMolDescriptors.CalcNumAmideBonds`.
    calc_num_amide_bonds;
    /// The atomic stereocenters, specified or not. Mirrors
    /// `rdMolDescriptors.CalcNumAtomStereoCenters`.
    calc_num_atom_stereo_centers;
    /// The unspecified atomic stereocenters. Mirrors
    /// `rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters`.
    calc_num_unspecified_atom_stereo_centers;
}

counts! {
    /// The nitrogens and oxygens. Mirrors
    /// `rdMolDescriptors.CalcNumLipinskiHBA`.
    calc_num_lipinski_hba;
    /// The aromatic rings. Mirrors `rdMolDescriptors.CalcNumAromaticRings`.
    calc_num_aromatic_rings;
    /// The rings of single, non-aromatic bonds. Mirrors
    /// `rdMolDescriptors.CalcNumSaturatedRings`.
    calc_num_saturated_rings;
    /// The rings with an atom other than carbon. Mirrors
    /// `rdMolDescriptors.CalcNumHeterocycles`.
    calc_num_heterocycles;
    /// Mirrors `rdMolDescriptors.CalcNumAromaticHeterocycles`.
    calc_num_aromatic_heterocycles;
    /// Mirrors `rdMolDescriptors.CalcNumAromaticCarbocycles`.
    calc_num_aromatic_carbocycles;
    /// Mirrors `rdMolDescriptors.CalcNumSaturatedHeterocycles`.
    calc_num_saturated_heterocycles;
    /// Mirrors `rdMolDescriptors.CalcNumSaturatedCarbocycles`.
    calc_num_saturated_carbocycles;
    /// The rings with a non-aromatic bond. Mirrors
    /// `rdMolDescriptors.CalcNumAliphaticRings`.
    calc_num_aliphatic_rings;
    /// Mirrors `rdMolDescriptors.CalcNumAliphaticHeterocycles`.
    calc_num_aliphatic_heterocycles;
    /// Mirrors `rdMolDescriptors.CalcNumAliphaticCarbocycles`.
    calc_num_aliphatic_carbocycles;
    /// The atoms other than hydrogen. Mirrors
    /// `rdMolDescriptors.CalcNumHeavyAtoms`.
    calc_num_heavy_atoms;
}

/// The rotatable bonds, by RDKit's default definition: `Strict`, unless
/// RDKit was built without `RDK_USE_STRICT_ROTOR_DEFINITION`. Mirrors
/// `rdMolDescriptors.CalcNumRotatableBonds`.
pub fn calc_num_rotatable_bonds(mol: &Mol) -> Result<u32> {
    calc_num_rotatable_bonds_with_strict(mol, NumRotatableBondsOptions::Default)
}

/// By the given definition, or `bool` (`Strict` or `NonStrict`). Mirrors
/// `rdMolDescriptors.CalcNumRotatableBonds(mol, strict)`, both overloads.
pub fn calc_num_rotatable_bonds_with_strict(
    mol: &Mol,
    strict: impl Into<NumRotatableBondsOptions>,
) -> Result<u32> {
    Ok(ffi::calc_num_rotatable_bonds(
        mol.as_ffi(),
        strict.into().to_ffi(),
    )?)
}

/// The fraction of carbons that are sp3 (four neighbors, hydrogens
/// included); 0 without carbons. Mirrors `rdMolDescriptors.CalcFractionCSP3`.
pub fn calc_fraction_csp3(mol: &Mol) -> Result<f64> {
    Ok(ffi::calc_fraction_csp3(mol.as_ffi())?)
}

/// What `calc_num_spiro_atoms_with_outputs` fills in (A9).
/// A Rust addition, for Python's output argument (A9).
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct CalcNumSpiroAtomsOutputs {
    /// The spiro atoms. Python's `atoms`.
    pub atoms: Vec<u32>,
}

/// The atoms shared by two rings that share only them. Mirrors
/// `rdMolDescriptors.CalcNumSpiroAtoms`.
pub fn calc_num_spiro_atoms(mol: &Mol) -> Result<u32> {
    Ok(calc_num_spiro_atoms_with_outputs(mol)?.0)
}

/// With the atoms. Mirrors `rdMolDescriptors.CalcNumSpiroAtoms(mol, atoms)`.
pub fn calc_num_spiro_atoms_with_outputs(mol: &Mol) -> Result<(u32, CalcNumSpiroAtomsOutputs)> {
    let mut atoms = Vec::new();
    let count = ffi::calc_num_spiro_atoms(mol.as_ffi(), &mut atoms)?;
    Ok((count, CalcNumSpiroAtomsOutputs { atoms }))
}

/// What `calc_num_bridgehead_atoms_with_outputs` fills in (A9).
/// A Rust addition, for Python's output argument (A9).
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct CalcNumBridgeheadAtomsOutputs {
    /// The bridgehead atoms. Python's `atoms`.
    pub atoms: Vec<u32>,
}

/// The atoms shared by rings that share at least two bonds. Mirrors
/// `rdMolDescriptors.CalcNumBridgeheadAtoms`.
pub fn calc_num_bridgehead_atoms(mol: &Mol) -> Result<u32> {
    Ok(calc_num_bridgehead_atoms_with_outputs(mol)?.0)
}

/// With the atoms. Mirrors `rdMolDescriptors.CalcNumBridgeheadAtoms(mol,
/// atoms)`.
pub fn calc_num_bridgehead_atoms_with_outputs(
    mol: &Mol,
) -> Result<(u32, CalcNumBridgeheadAtomsOutputs)> {
    let mut atoms = Vec::new();
    let count = ffi::calc_num_bridgehead_atoms(mol.as_ffi(), &mut atoms)?;
    Ok((count, CalcNumBridgeheadAtomsOutputs { atoms }))
}

/// The Labute surface area binned by Crippen logP contribution, with
/// RDKit's 11 default bins. Mirrors `rdMolDescriptors.SlogP_VSA_`; `force`
/// is not bound (U8).
pub fn slogp_vsa(mol: &Mol) -> Result<Vec<f64>> {
    slogp_vsa_with_bins(mol, &[])
}

/// With the given bin boundaries, which must be sorted, without NaN; RDKit
/// bins with `std::upper_bound`, which is undefined otherwise, so others are
/// `Error::Value`. Empty means the defaults, as Python's empty list does.
/// Mirrors `rdMolDescriptors.SlogP_VSA_(mol, bins)`.
pub fn slogp_vsa_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    Ok(ffi::slogp_vsa(mol.as_ffi(), bins)?)
}

/// The Labute surface area binned by Crippen MR contribution, with RDKit's
/// 9 default bins. Mirrors `rdMolDescriptors.SMR_VSA_`; `force` is not
/// bound (U8).
pub fn smr_vsa(mol: &Mol) -> Result<Vec<f64>> {
    smr_vsa_with_bins(mol, &[])
}

/// With the given bins, as `slogp_vsa_with_bins`. Mirrors
/// `rdMolDescriptors.SMR_VSA_(mol, bins)`.
pub fn smr_vsa_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    Ok(ffi::smr_vsa(mol.as_ffi(), bins)?)
}

/// The Labute surface area binned by Gasteiger charge, with RDKit's 13
/// default bins. Mirrors `rdMolDescriptors.PEOE_VSA_`; `force` is not bound
/// (U8).
pub fn peoe_vsa(mol: &Mol) -> Result<Vec<f64>> {
    peoe_vsa_with_bins(mol, &[])
}

/// With the given bins, as `slogp_vsa_with_bins`. Mirrors
/// `rdMolDescriptors.PEOE_VSA_(mol, bins)`.
pub fn peoe_vsa_with_bins(mol: &Mol, bins: &[f64]) -> Result<Vec<f64>> {
    Ok(ffi::peoe_vsa(mol.as_ffi(), bins)?)
}

/// The Labute surface area binned by each atom's `f64` property
/// `custom_prop_name` (1 where an atom lacks it). The bins must be sorted,
/// without NaN, as in `slogp_vsa_with_bins`. Mirrors
/// `rdMolDescriptors.CustomProp_VSA_`; `force` is not bound (U8).
pub fn custom_prop_vsa(mol: &Mol, custom_prop_name: &str, bins: &[f64]) -> Result<Vec<f64>> {
    Ok(ffi::custom_prop_vsa(mol.as_ffi(), custom_prop_name, bins)?)
}

/// The 42 molecular quantum numbers. Mirrors `rdMolDescriptors.MQNs_`,
/// whose `force` RDKit ignores.
pub fn mqns(mol: &Mol) -> Result<Vec<u32>> {
    Ok(ffi::mqns(mol.as_ffi())?)
}

/// Defines functions of a molecule that return an `f64`.
macro_rules! real_descriptors {
    ($($(#[$meta:meta])* $name:ident;)*) => {
        $(
            $(#[$meta])*
            pub fn $name(mol: &Mol) -> Result<f64> {
                Ok(ffi::$name(mol.as_ffi())?)
            }
        )*
    };
}

real_descriptors! {
    /// Kier and Hall's valence connectivity index of order 0. Mirrors
    /// `rdMolDescriptors.CalcChi0v`; `force` is not bound (U8).
    calc_chi0v;
    /// Of order 1. Mirrors `rdMolDescriptors.CalcChi1v`.
    calc_chi1v;
    /// Of order 2. Mirrors `rdMolDescriptors.CalcChi2v`.
    calc_chi2v;
    /// Of order 3. Mirrors `rdMolDescriptors.CalcChi3v`.
    calc_chi3v;
    /// Of order 4. Mirrors `rdMolDescriptors.CalcChi4v`.
    calc_chi4v;
    /// As `calc_chi0v`, from valence electrons rather than valence. Mirrors
    /// `rdMolDescriptors.CalcChi0n`; `force` is not bound (U8).
    calc_chi0n;
    /// Of order 1. Mirrors `rdMolDescriptors.CalcChi1n`.
    calc_chi1n;
    /// Of order 2. Mirrors `rdMolDescriptors.CalcChi2n`.
    calc_chi2n;
    /// Of order 3. Mirrors `rdMolDescriptors.CalcChi3n`.
    calc_chi3n;
    /// Of order 4. Mirrors `rdMolDescriptors.CalcChi4n`.
    calc_chi4n;
    /// Kier's first shape index. Mirrors `rdMolDescriptors.CalcKappa1`.
    calc_kappa1;
    /// The second. Mirrors `rdMolDescriptors.CalcKappa2`.
    calc_kappa2;
    /// The third. Mirrors `rdMolDescriptors.CalcKappa3`.
    calc_kappa3;
    /// Kier's flexibility index. Mirrors `rdMolDescriptors.CalcPhi`.
    calc_phi;
}

/// The valence connectivity index of order `n`; 0 when `n` exceeds the
/// number of atoms, which no path is long enough for. Mirrors
/// `rdMolDescriptors.CalcChiNv`; `force` is not bound (U8).
pub fn calc_chi_nv(mol: &Mol, n: u32) -> Result<f64> {
    Ok(ffi::calc_chi_nv(mol.as_ffi(), n)?)
}

/// As `calc_chi_nv`, from valence electrons. Mirrors
/// `rdMolDescriptors.CalcChiNn`; `force` is not bound (U8).
pub fn calc_chi_nn(mol: &Mol, n: u32) -> Result<f64> {
    Ok(ffi::calc_chi_nn(mol.as_ffi(), n)?)
}

/// What `calc_hall_kier_alpha_with_outputs` fills in (A9).
/// A Rust addition, for Python's output argument (A9).
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct CalcHallKierAlphaOutputs {
    /// Each atom's alpha; 0 for dummy atoms, which RDKit skips. Python's
    /// `atomContribs`.
    pub atom_contribs: Vec<f64>,
}

/// Hall and Kier's alpha, from atom radii and hybridizations. Mirrors
/// `rdMolDescriptors.CalcHallKierAlpha`.
pub fn calc_hall_kier_alpha(mol: &Mol) -> Result<f64> {
    Ok(calc_hall_kier_alpha_with_outputs(mol)?.0)
}

/// With each atom's contribution. Mirrors
/// `rdMolDescriptors.CalcHallKierAlpha(mol, atomContribs)`.
pub fn calc_hall_kier_alpha_with_outputs(mol: &Mol) -> Result<(f64, CalcHallKierAlphaOutputs)> {
    let mut atom_contribs = Vec::new();
    let alpha = ffi::calc_hall_kier_alpha(mol.as_ffi(), &mut atom_contribs)?;
    Ok((alpha, CalcHallKierAlphaOutputs { atom_contribs }))
}

/// An atom's 49 one-hot features: element, degree, hybridization, implicit
/// valence, charge, ring sizes, aromaticity, hydrogens, and 1 / the number
/// of atoms. `Error::Index` for an atom out of range. Mirrors
/// `rdMolDescriptors.GetAtomFeatures`.
pub fn get_atom_features(mol: &Mol, atomid: u32) -> Result<Vec<f64>> {
    get_atom_features_with_addchiral(mol, atomid, false)
}

/// With three more for the chiral tag when `addchiral`. Mirrors
/// `rdMolDescriptors.GetAtomFeatures(mol, atomid, addchiral)`.
pub fn get_atom_features_with_addchiral(
    mol: &Mol,
    atomid: u32,
    addchiral: bool,
) -> Result<Vec<f64>> {
    Ok(ffi::get_atom_features(mol.as_ffi(), atomid, addchiral)?)
}

/// Sets each atom's `i32` property `OxidationNumber` from Pauling
/// electronegativities, on a kekulized copy. Mirrors
/// `rdMolDescriptors.CalcOxidationNumbers`.
pub fn calc_oxidation_numbers(mol: &mut Mol) -> Result<()> {
    Ok(ffi::calc_oxidation_numbers(mol.pin_mut())?)
}

/// The 192 Moreau-Broto autocorrelations of six atom properties over
/// topological distances 1 to 8. RDKit reads the properties from tables
/// indexed by atomic number that cover 1 to 110, and past their ends for
/// other atoms, so dummy atoms and heavier ones are `Error::Value`; so are
/// molecules above 46340 atoms (U3). Mirrors `rdMolDescriptors.CalcAUTOCORR2D`.
pub fn calc_autocorr_2d(mol: &Mol) -> Result<Vec<f64>> {
    calc_autocorr_2d_with_custom_atom_property(mol, "")
}

/// The 32 autocorrelations of each atom's `f64` property
/// `custom_atom_property` (1 where an atom lacks it), which takes any atoms;
/// empty means the six table properties. Mirrors
/// `rdMolDescriptors.CalcAUTOCORR2D(mol, CustomAtomProperty)`.
pub fn calc_autocorr_2d_with_custom_atom_property(
    mol: &Mol,
    custom_atom_property: &str,
) -> Result<Vec<f64>> {
    Ok(ffi::calc_autocorr_2d(mol.as_ffi(), custom_atom_property)?)
}

/// The eight BCUT descriptors: the highest and lowest eigenvalues of the
/// Burden matrix with mass, Gasteiger charge, Crippen logP, and Crippen MR on
/// the diagonal, for the molecule without hydrogens. Mirrors
/// `rdMolDescriptors.BCUT2D(mol)`.
pub fn bcut_2d(mol: &Mol) -> Result<Vec<f64>> {
    Ok(ffi::bcut_2d(mol.as_ffi())?)
}

/// The (highest, lowest) eigenvalues with `atom_props`, one per atom, on the
/// diagonal; `Error::Runtime` for another number of values. Mirrors
/// `rdMolDescriptors.BCUT2D(mol, atom_props)`, for a list or a tuple.
pub fn bcut_2d_atom_props(mol: &Mol, atom_props: &[f64]) -> Result<(f64, f64)> {
    Ok(pair(ffi::bcut_2d_atom_props(mol.as_ffi(), atom_props)?))
}

/// With each atom's `f64` property `atom_propname` on the diagonal;
/// `Error::Key` when an atom lacks it. Mirrors
/// `rdMolDescriptors.BCUT2D(mol, atom_propname)`.
pub fn bcut_2d_atom_propname(mol: &Mol, atom_propname: &str) -> Result<(f64, f64)> {
    Ok(pair(ffi::bcut_2d_atom_propname(
        mol.as_ffi(),
        atom_propname,
    )?))
}
