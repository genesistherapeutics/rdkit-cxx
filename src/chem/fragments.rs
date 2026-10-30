// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Fragments`, a Python module: counts of functional
//! groups and ring systems, each the number of matches of a SMARTS pattern
//! from RDKit's `Data/FragmentDescriptors.csv`.
//!
//! Python reads the file from RDKit's data directory on import and makes a
//! function per pattern; the crate embeds the file (vendored in
//! `wrapper/third_party/rdkit`) and defines the same functions. Each
//! parses its pattern anew, since a `Mol` cannot be shared between threads
//! (CONVENTIONS.md S1). Python's functions take the pattern as a defaulted
//! argument; here each function is its pattern.

use crate::Result;
use crate::chem::{Mol, SubstructMatchParameters, mol_from_smarts};

/// `Data/FragmentDescriptors.csv` from RDKit 2024.03.6: a code, a
/// description, a SMARTS pattern, and notes per tab-separated line.
const FRAGMENT_DESCRIPTORS: &str =
    include_str!("../../wrapper/third_party/rdkit/Data/FragmentDescriptors.csv");

/// The SMARTS of the pattern whose code is `code`, as Python's
/// `_LoadPatterns` reads it.
fn smarts(code: &str) -> &'static str {
    FRAGMENT_DESCRIPTORS
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|fields| fields.len() >= 3 && fields[0] == code)
        .map(|fields| fields[2])
        .expect("the fragment code is in FragmentDescriptors.csv")
}

/// The number of matches of the pattern whose code is `code`, unique ones
/// only when `count_unique`. Python's `_CountMatches`.
fn count_matches(mol: &Mol, code: &str, count_unique: bool) -> Result<u32> {
    let pattern = mol_from_smarts(smarts(code)).expect("the fragment SMARTS parse");
    let params = SubstructMatchParameters {
        uniquify: count_unique,
        ..SubstructMatchParameters::default()
    };
    Ok(mol
        .get_substruct_matches_with_params(&pattern, &params)?
        .len() as u32)
}

/// A fragment count: a function of a molecule.
/// A Rust addition, naming the functions of `fns`.
pub type FragmentFn = fn(&Mol) -> Result<u32>;

/// Defines a function per pattern, and its form taking `count_unique`.
macro_rules! fragments {
    ($($name:ident, $with_count_unique:ident, $code:literal, $python:literal, $doc:literal;)*) => {
        $(
            #[doc = concat!(
                $doc, ": the unique matches of the pattern `", $code,
                "` in `FragmentDescriptors.csv`. Mirrors `Fragments.", $python, "`."
            )]
            pub fn $name(mol: &Mol) -> Result<u32> {
                count_matches(mol, $code, true)
            }

            #[doc = concat!(
                "Every match when `count_unique` is false. Mirrors `Fragments.", $python,
                "(mol, countUnique)`."
            )]
            pub fn $with_count_unique(mol: &Mol, count_unique: bool) -> Result<u32> {
                count_matches(mol, $code, count_unique)
            }
        )*
    };
}

fragments! {
    fr_c_o, fr_c_o_with_count_unique, "fr_C=O", "fr_C_O", "Number of carbonyl O";
    fr_c_o_no_coo, fr_c_o_no_coo_with_count_unique, "fr_C=O_noCOO", "fr_C_O_noCOO", "Number of carbonyl O, excluding COOH";
    fr_al_oh, fr_al_oh_with_count_unique, "fr_Al_OH", "fr_Al_OH", "Number of aliphatic hydroxyl groups";
    fr_ar_oh, fr_ar_oh_with_count_unique, "fr_Ar_OH", "fr_Ar_OH", "Number of aromatic hydroxyl groups";
    fr_methoxy, fr_methoxy_with_count_unique, "fr_methoxy", "fr_methoxy", "Number of methoxy groups -OCH3";
    fr_oxime, fr_oxime_with_count_unique, "fr_oxime", "fr_oxime", "Number of oxime groups";
    fr_ester, fr_ester_with_count_unique, "fr_ester", "fr_ester", "Number of esters";
    fr_al_coo, fr_al_coo_with_count_unique, "fr_Al_COO", "fr_Al_COO", "Number of aliphatic carboxylic acids";
    fr_ar_coo, fr_ar_coo_with_count_unique, "fr_Ar_COO", "fr_Ar_COO", "Number of Aromatic carboxylic acide";
    fr_coo, fr_coo_with_count_unique, "fr_COO", "fr_COO", "Number of carboxylic acids";
    fr_coo2, fr_coo2_with_count_unique, "fr_COO2", "fr_COO2", "Number of carboxylic acids";
    fr_ketone, fr_ketone_with_count_unique, "fr_ketone", "fr_ketone", "Number of ketones";
    fr_ether, fr_ether_with_count_unique, "fr_ether", "fr_ether", "Number of ether oxygens (including phenoxy)";
    fr_phenol, fr_phenol_with_count_unique, "fr_phenol", "fr_phenol", "Number of phenols";
    fr_aldehyde, fr_aldehyde_with_count_unique, "fr_aldehyde", "fr_aldehyde", "Number of aldehydes";
    fr_quat_n, fr_quat_n_with_count_unique, "fr_quatN", "fr_quatN", "Number of quaternary nitrogens";
    fr_nh2, fr_nh2_with_count_unique, "fr_NH2", "fr_NH2", "Number of Primary amines";
    fr_nh1, fr_nh1_with_count_unique, "fr_NH1", "fr_NH1", "Number of Secondary amines";
    fr_nh0, fr_nh0_with_count_unique, "fr_NH0", "fr_NH0", "Number of Tertiary amines";
    fr_ar_n, fr_ar_n_with_count_unique, "fr_Ar_N", "fr_Ar_N", "Number of aromatic nitrogens";
    fr_ar_nh, fr_ar_nh_with_count_unique, "fr_Ar_NH", "fr_Ar_NH", "Number of aromatic amines";
    fr_aniline, fr_aniline_with_count_unique, "fr_aniline", "fr_aniline", "Number of anilines";
    fr_imine, fr_imine_with_count_unique, "fr_Imine", "fr_Imine", "Number of Imines";
    fr_nitrile, fr_nitrile_with_count_unique, "fr_nitrile", "fr_nitrile", "Number of nitriles";
    fr_hdrzine, fr_hdrzine_with_count_unique, "fr_hdrzine", "fr_hdrzine", "Number of hydrazine groups";
    fr_hdrzone, fr_hdrzone_with_count_unique, "fr_hdrzone", "fr_hdrzone", "Number of hydrazone groups";
    fr_nitroso, fr_nitroso_with_count_unique, "fr_nitroso", "fr_nitroso", "Number of nitroso groups, excluding NO2";
    fr_n_o, fr_n_o_with_count_unique, "fr_N-O", "fr_N_O", "Number of hydroxylamine groups";
    fr_nitro, fr_nitro_with_count_unique, "fr_nitro", "fr_nitro", "Number of nitro groups";
    fr_azo, fr_azo_with_count_unique, "fr_azo", "fr_azo", "Number of azo groups";
    fr_diazo, fr_diazo_with_count_unique, "fr_diazo", "fr_diazo", "Number of diazo groups";
    fr_azide, fr_azide_with_count_unique, "fr_azide", "fr_azide", "Number of azide groups";
    fr_amide, fr_amide_with_count_unique, "fr_amide", "fr_amide", "Number of amides";
    fr_priamide, fr_priamide_with_count_unique, "fr_priamide", "fr_priamide", "Number of primary amides";
    fr_amidine, fr_amidine_with_count_unique, "fr_amidine", "fr_amidine", "Number of amidine groups";
    fr_guanido, fr_guanido_with_count_unique, "fr_guanido", "fr_guanido", "Number of guanidine groups";
    fr_nhpyrrole, fr_nhpyrrole_with_count_unique, "fr_Nhpyrrole", "fr_Nhpyrrole", "Number of H-pyrrole nitrogens";
    fr_imide, fr_imide_with_count_unique, "fr_imide", "fr_imide", "Number of imide groups";
    fr_isocyan, fr_isocyan_with_count_unique, "fr_isocyan", "fr_isocyan", "Number of isocyanates";
    fr_isothiocyan, fr_isothiocyan_with_count_unique, "fr_isothiocyan", "fr_isothiocyan", "Number of isothiocyanates";
    fr_thiocyan, fr_thiocyan_with_count_unique, "fr_thiocyan", "fr_thiocyan", "Number of thiocyanates";
    fr_halogen, fr_halogen_with_count_unique, "fr_halogen", "fr_halogen", "Number of halogens";
    fr_alkyl_halide, fr_alkyl_halide_with_count_unique, "fr_alkyl_halide", "fr_alkyl_halide", "Number of alkyl halides";
    fr_sulfide, fr_sulfide_with_count_unique, "fr_sulfide", "fr_sulfide", "Number of thioether";
    fr_sh, fr_sh_with_count_unique, "fr_SH", "fr_SH", "Number of thiol groups";
    fr_c_s, fr_c_s_with_count_unique, "fr_C=S", "fr_C_S", "Number of thiocarbonyl";
    fr_sulfone, fr_sulfone_with_count_unique, "fr_sulfone", "fr_sulfone", "Number of sulfone groups";
    fr_sulfonamd, fr_sulfonamd_with_count_unique, "fr_sulfonamd", "fr_sulfonamd", "Number of sulfonamides";
    fr_prisulfonamd, fr_prisulfonamd_with_count_unique, "fr_prisulfonamd", "fr_prisulfonamd", "Number of primary sulfonamides";
    fr_barbitur, fr_barbitur_with_count_unique, "fr_barbitur", "fr_barbitur", "Number of barbiturate groups";
    fr_urea, fr_urea_with_count_unique, "fr_urea", "fr_urea", "Number of urea groups";
    fr_term_acetylene, fr_term_acetylene_with_count_unique, "fr_term_acetylene", "fr_term_acetylene", "Number of terminal acetylenes";
    fr_imidazole, fr_imidazole_with_count_unique, "fr_imidazole", "fr_imidazole", "Number of imidazole rings";
    fr_furan, fr_furan_with_count_unique, "fr_furan", "fr_furan", "Number of furan rings";
    fr_thiophene, fr_thiophene_with_count_unique, "fr_thiophene", "fr_thiophene", "Number of thiophene rings";
    fr_thiazole, fr_thiazole_with_count_unique, "fr_thiazole", "fr_thiazole", "Number of thiazole rings";
    fr_oxazole, fr_oxazole_with_count_unique, "fr_oxazole", "fr_oxazole", "Number of oxazole rings";
    fr_pyridine, fr_pyridine_with_count_unique, "fr_pyridine", "fr_pyridine", "Number of pyridine rings";
    fr_piperdine, fr_piperdine_with_count_unique, "fr_piperdine", "fr_piperdine", "Number of piperdine rings";
    fr_piperzine, fr_piperzine_with_count_unique, "fr_piperzine", "fr_piperzine", "Number of piperzine rings";
    fr_morpholine, fr_morpholine_with_count_unique, "fr_morpholine", "fr_morpholine", "Number of morpholine rings";
    fr_lactam, fr_lactam_with_count_unique, "fr_lactam", "fr_lactam", "Number of beta lactams";
    fr_lactone, fr_lactone_with_count_unique, "fr_lactone", "fr_lactone", "Number of cyclic esters (lactones)";
    fr_tetrazole, fr_tetrazole_with_count_unique, "fr_tetrazole", "fr_tetrazole", "Number of tetrazole rings";
    fr_epoxide, fr_epoxide_with_count_unique, "fr_epoxide", "fr_epoxide", "Number of epoxide rings";
    fr_unbrch_alkane, fr_unbrch_alkane_with_count_unique, "fr_unbrch_alkane", "fr_unbrch_alkane", "Number of unbranched alkanes  of at least 4 members (excludes halogenated alkanes)";
    fr_bicyclic, fr_bicyclic_with_count_unique, "fr_bicyclic", "fr_bicyclic", "Bicyclic";
    fr_benzene, fr_benzene_with_count_unique, "fr_benzene", "fr_benzene", "Number of benzene rings";
    fr_phos_acid, fr_phos_acid_with_count_unique, "fr_phos_acid", "fr_phos_acid", "Number of phosphoric acid groups";
    fr_phos_ester, fr_phos_ester_with_count_unique, "fr_phos_ester", "fr_phos_ester", "Number of phosphoric ester groups";
    fr_nitro_arom, fr_nitro_arom_with_count_unique, "fr_nitro_arom", "fr_nitro_arom", "Number of nitro benzene ring substituents";
    fr_nitro_arom_nonortho, fr_nitro_arom_nonortho_with_count_unique, "fr_nitro_arom_nonortho", "fr_nitro_arom_nonortho", "Number of non-ortho nitro benzene ring substituents";
    fr_dihydropyridine, fr_dihydropyridine_with_count_unique, "fr_dihydropyridine", "fr_dihydropyridine", "Number of dihydropyridines";
    fr_phenol_no_ortho_hbond, fr_phenol_no_ortho_hbond_with_count_unique, "fr_phenol_noOrthoHbond", "fr_phenol_noOrthoHbond", "Number of phenolic OH excluding ortho intramolecular Hbond substituents";
    fr_al_oh_no_tert, fr_al_oh_no_tert_with_count_unique, "fr_Al_OH_noTert", "fr_Al_OH_noTert", "Number of aliphatic hydroxyl groups excluding tert-OH";
    fr_benzodiazepine, fr_benzodiazepine_with_count_unique, "fr_benzodiazepine", "fr_benzodiazepine", "Number of benzodiazepines with no additional fused rings";
    fr_para_hydroxylation, fr_para_hydroxylation_with_count_unique, "fr_para_hydroxylation", "fr_para_hydroxylation", "Number of para-hydroxylation sites";
    fr_allylic_oxid, fr_allylic_oxid_with_count_unique, "fr_allylic_oxid", "fr_allylic_oxid", "Number of allylic oxidation sites excluding steroid dienone";
    fr_aryl_methyl, fr_aryl_methyl_with_count_unique, "fr_aryl_methyl", "fr_aryl_methyl", "Number of aryl methyl sites for hydroxylation";
    fr_ndealkylation1, fr_ndealkylation1_with_count_unique, "fr_Ndealkylation1", "fr_Ndealkylation1", "Number of XCCNR groups";
    fr_ndealkylation2, fr_ndealkylation2_with_count_unique, "fr_Ndealkylation2", "fr_Ndealkylation2", "Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)";
    fr_alkyl_carbamate, fr_alkyl_carbamate_with_count_unique, "fr_alkyl_carbamate", "fr_alkyl_carbamate", "Number of alkyl carbamates (subject to hydrolysis)";
    fr_ketone_topliss, fr_ketone_topliss_with_count_unique, "fr_ketone_Topliss", "fr_ketone_Topliss", "Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha";
    fr_arn, fr_arn_with_count_unique, "fr_ArN", "fr_ArN", "Number of N functional groups attached to aromatics";
    fr_hoccn, fr_hoccn_with_count_unique, "fr_HOCCN", "fr_HOCCN", "Number of C(OH)CCN-Ctert-alkyl or  C(OH)CCNcyclic";
}

/// The fragment counts' names and functions, in the file's order. Mirrors
/// `Fragments.fns`.
pub fn fns() -> Vec<(&'static str, FragmentFn)> {
    vec![
        ("fr_C_O", fr_c_o),
        ("fr_C_O_noCOO", fr_c_o_no_coo),
        ("fr_Al_OH", fr_al_oh),
        ("fr_Ar_OH", fr_ar_oh),
        ("fr_methoxy", fr_methoxy),
        ("fr_oxime", fr_oxime),
        ("fr_ester", fr_ester),
        ("fr_Al_COO", fr_al_coo),
        ("fr_Ar_COO", fr_ar_coo),
        ("fr_COO", fr_coo),
        ("fr_COO2", fr_coo2),
        ("fr_ketone", fr_ketone),
        ("fr_ether", fr_ether),
        ("fr_phenol", fr_phenol),
        ("fr_aldehyde", fr_aldehyde),
        ("fr_quatN", fr_quat_n),
        ("fr_NH2", fr_nh2),
        ("fr_NH1", fr_nh1),
        ("fr_NH0", fr_nh0),
        ("fr_Ar_N", fr_ar_n),
        ("fr_Ar_NH", fr_ar_nh),
        ("fr_aniline", fr_aniline),
        ("fr_Imine", fr_imine),
        ("fr_nitrile", fr_nitrile),
        ("fr_hdrzine", fr_hdrzine),
        ("fr_hdrzone", fr_hdrzone),
        ("fr_nitroso", fr_nitroso),
        ("fr_N_O", fr_n_o),
        ("fr_nitro", fr_nitro),
        ("fr_azo", fr_azo),
        ("fr_diazo", fr_diazo),
        ("fr_azide", fr_azide),
        ("fr_amide", fr_amide),
        ("fr_priamide", fr_priamide),
        ("fr_amidine", fr_amidine),
        ("fr_guanido", fr_guanido),
        ("fr_Nhpyrrole", fr_nhpyrrole),
        ("fr_imide", fr_imide),
        ("fr_isocyan", fr_isocyan),
        ("fr_isothiocyan", fr_isothiocyan),
        ("fr_thiocyan", fr_thiocyan),
        ("fr_halogen", fr_halogen),
        ("fr_alkyl_halide", fr_alkyl_halide),
        ("fr_sulfide", fr_sulfide),
        ("fr_SH", fr_sh),
        ("fr_C_S", fr_c_s),
        ("fr_sulfone", fr_sulfone),
        ("fr_sulfonamd", fr_sulfonamd),
        ("fr_prisulfonamd", fr_prisulfonamd),
        ("fr_barbitur", fr_barbitur),
        ("fr_urea", fr_urea),
        ("fr_term_acetylene", fr_term_acetylene),
        ("fr_imidazole", fr_imidazole),
        ("fr_furan", fr_furan),
        ("fr_thiophene", fr_thiophene),
        ("fr_thiazole", fr_thiazole),
        ("fr_oxazole", fr_oxazole),
        ("fr_pyridine", fr_pyridine),
        ("fr_piperdine", fr_piperdine),
        ("fr_piperzine", fr_piperzine),
        ("fr_morpholine", fr_morpholine),
        ("fr_lactam", fr_lactam),
        ("fr_lactone", fr_lactone),
        ("fr_tetrazole", fr_tetrazole),
        ("fr_epoxide", fr_epoxide),
        ("fr_unbrch_alkane", fr_unbrch_alkane),
        ("fr_bicyclic", fr_bicyclic),
        ("fr_benzene", fr_benzene),
        ("fr_phos_acid", fr_phos_acid),
        ("fr_phos_ester", fr_phos_ester),
        ("fr_nitro_arom", fr_nitro_arom),
        ("fr_nitro_arom_nonortho", fr_nitro_arom_nonortho),
        ("fr_dihydropyridine", fr_dihydropyridine),
        ("fr_phenol_noOrthoHbond", fr_phenol_no_ortho_hbond),
        ("fr_Al_OH_noTert", fr_al_oh_no_tert),
        ("fr_benzodiazepine", fr_benzodiazepine),
        ("fr_para_hydroxylation", fr_para_hydroxylation),
        ("fr_allylic_oxid", fr_allylic_oxid),
        ("fr_aryl_methyl", fr_aryl_methyl),
        ("fr_Ndealkylation1", fr_ndealkylation1),
        ("fr_Ndealkylation2", fr_ndealkylation2),
        ("fr_alkyl_carbamate", fr_alkyl_carbamate),
        ("fr_ketone_Topliss", fr_ketone_topliss),
        ("fr_ArN", fr_arn),
        ("fr_HOCCN", fr_hoccn),
    ]
}
