// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Descriptors`, a Python module: molecular weights,
//! electron counts, partial-charge extremes, fingerprint densities, and the
//! list of every 2D descriptor RDKit's Python offers (`descList`), from
//! `Lipinski`, `Crippen`, `MolSurf`, `GraphDescriptors`, `EState`,
//! `Fragments`, `QED`, `SpacialScore`, and `rdMolDescriptors`.
//!
//! Python gathers `descList` from those modules' namespaces on import; here
//! it is a table of the ported functions, in Python's order, under Python's
//! names. Python's `force` arguments and the results it caches on the
//! molecule are not ported: every call computes its result afresh.

use std::collections::HashSet;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::Result;
use crate::chem::e_state::e_state_vsa;
use crate::chem::rd_fingerprint_generator::{
    GetMorganGeneratorParams, get_morgan_generator_with_params,
};
use crate::chem::rd_mol_descriptors as rd_descriptors;
use crate::chem::{
    Mol, crippen, e_state, fragments, get_periodic_table, graph_descriptors, lipinski, mol_surf,
    qed, spacial_score,
};
use crate::py::{py_max, py_min};

/// The average molecular weight. Mirrors `Descriptors.MolWt`, which calls
/// `rdMolDescriptors._CalcMolWt`.
pub fn mol_wt(mol: &Mol) -> Result<f64> {
    mol_wt_with_only_heavy(mol, false)
}

/// Of the heavy atoms alone when `only_heavy`. Mirrors
/// `Descriptors.MolWt(mol, onlyHeavy)`.
pub fn mol_wt_with_only_heavy(mol: &Mol, only_heavy: bool) -> Result<f64> {
    rd_descriptors::calc_mol_wt(mol, only_heavy)
}

/// The average molecular weight ignoring hydrogens. Mirrors
/// `Descriptors.HeavyAtomMolWt`.
pub fn heavy_atom_mol_wt(mol: &Mol) -> Result<f64> {
    mol_wt_with_only_heavy(mol, true)
}

/// The monoisotopic molecular weight. Mirrors `Descriptors.ExactMolWt`,
/// which calls `rdMolDescriptors.CalcExactMolWt`.
pub fn exact_mol_wt(mol: &Mol) -> Result<f64> {
    rd_descriptors::calc_exact_mol_wt(mol)
}

/// Of the heavy atoms alone when `only_heavy`. Mirrors
/// `Descriptors.ExactMolWt(mol, onlyHeavy)`.
pub fn exact_mol_wt_with_only_heavy(mol: &Mol, only_heavy: bool) -> Result<f64> {
    rd_descriptors::calc_exact_mol_wt_with_only_heavy(mol, only_heavy)
}

/// An atomic number as the periodic table takes it; a negative one, which
/// it lacks, is refused there.
pub(crate) fn atomic_number(num: i32) -> u32 {
    u32::try_from(num).unwrap_or(u32::MAX)
}

/// The number of valence electrons: each atom's outer-shell electrons, less
/// its charge, plus its hydrogens. Mirrors `Descriptors.NumValenceElectrons`.
pub fn num_valence_electrons(mol: &Mol) -> Result<i64> {
    let tbl = get_periodic_table();
    let mut total = 0;
    for atom in mol.get_atoms() {
        total += i64::from(tbl.get_n_outer_elecs(atomic_number(atom.get_atomic_num()))?)
            - i64::from(atom.get_formal_charge())
            + i64::from(atom.get_total_num_hs()?);
    }
    Ok(total)
}

/// The number of radical electrons, which says nothing of the spin state.
/// Mirrors `Descriptors.NumRadicalElectrons`.
pub fn num_radical_electrons(mol: &Mol) -> u64 {
    mol.get_atoms()
        .map(|atom| u64::from(atom.get_num_radical_electrons()))
        .sum()
}

/// The smallest and largest Gasteiger charges, from 500 and −500, as
/// Python's `_ChargeDescriptors` compares them (so NaN charges give
/// Python's answers).
fn charge_descriptors(mol: &Mol) -> Result<(f64, f64)> {
    let (mut min_chg, mut max_chg) = (500.0, -500.0);
    for chg in mol_surf::gasteiger_charges(mol)? {
        // Python's min(chg, minChg) and max(chg, maxChg).
        min_chg = if min_chg < chg { min_chg } else { chg };
        max_chg = if max_chg > chg { max_chg } else { chg };
    }
    Ok((min_chg, max_chg))
}

/// The largest Gasteiger charge, computed on a copy of the molecule.
/// Mirrors `Descriptors.MaxPartialCharge`.
pub fn max_partial_charge(mol: &Mol) -> Result<f64> {
    Ok(charge_descriptors(mol)?.1)
}

/// The smallest Gasteiger charge. Mirrors `Descriptors.MinPartialCharge`.
pub fn min_partial_charge(mol: &Mol) -> Result<f64> {
    Ok(charge_descriptors(mol)?.0)
}

/// The larger magnitude of the smallest and largest Gasteiger charges.
/// Mirrors `Descriptors.MaxAbsPartialCharge`.
pub fn max_abs_partial_charge(mol: &Mol) -> Result<f64> {
    let (v1, v2) = charge_descriptors(mol)?;
    py_max([v1.abs(), v2.abs()])
}

/// The smaller magnitude of the smallest and largest Gasteiger charges.
/// Mirrors `Descriptors.MinAbsPartialCharge`.
pub fn min_abs_partial_charge(mol: &Mol) -> Result<f64> {
    let (v1, v2) = charge_descriptors(mol)?;
    py_min([v1.abs(), v2.abs()])
}

/// The number of distinct Morgan environments of the given radius per heavy
/// atom; 0 without heavy atoms. Python's `_FingerprintDensity` of
/// `_getMorganCountFingerprint`.
fn fp_density_morgan(mol: &Mol, radius: u32) -> Result<f64> {
    let fpg = get_morgan_generator_with_params(&GetMorganGeneratorParams {
        radius,
        ..GetMorganGeneratorParams::default()
    })?;
    let val = fpg
        .get_sparse_count_fingerprint(mol)?
        .get_nonzero_elements()
        .len();
    match mol.get_num_heavy_atoms() {
        0 => Ok(0.0),
        num_heavy_atoms => Ok(val as f64 / f64::from(num_heavy_atoms)),
    }
}

/// The number of distinct Morgan environments of radius 1 per heavy atom.
/// Mirrors `Descriptors.FpDensityMorgan1`.
pub fn fp_density_morgan1(mol: &Mol) -> Result<f64> {
    fp_density_morgan(mol, 1)
}

/// Of radius 2. Mirrors `Descriptors.FpDensityMorgan2`.
pub fn fp_density_morgan2(mol: &Mol) -> Result<f64> {
    fp_density_morgan(mol, 2)
}

/// Of radius 3. Mirrors `Descriptors.FpDensityMorgan3`.
pub fn fp_density_morgan3(mol: &Mol) -> Result<f64> {
    fp_density_morgan(mol, 3)
}

/// The descriptors Python makes from `rdMolDescriptors` functions that
/// return vectors, with `ChemUtils.DescriptorUtilities.VectorDescriptorWrapper`:
/// each is one entry, or NaN when the function fails.
mod vector {
    use super::*;

    macro_rules! entries {
        ($function:path: $($name:ident = $index:literal,)*) => {
            $(
                pub(super) fn $name(mol: &Mol) -> Result<f64> {
                    Ok($function(mol).map_or(f64::NAN, |values| values[$index]))
                }
            )*
        };
    }

    entries! {
        rd_descriptors::bcut_2d:
        bcut2d_0 = 0, bcut2d_1 = 1, bcut2d_2 = 2, bcut2d_3 = 3, bcut2d_4 = 4, bcut2d_5 = 5,
        bcut2d_6 = 6, bcut2d_7 = 7,
    }

    entries! {
        rd_descriptors::calc_autocorr_2d:
        autocorr2d_0 = 0,
        autocorr2d_1 = 1,
        autocorr2d_2 = 2,
        autocorr2d_3 = 3,
        autocorr2d_4 = 4,
        autocorr2d_5 = 5,
        autocorr2d_6 = 6,
        autocorr2d_7 = 7,
        autocorr2d_8 = 8,
        autocorr2d_9 = 9,
        autocorr2d_10 = 10,
        autocorr2d_11 = 11,
        autocorr2d_12 = 12,
        autocorr2d_13 = 13,
        autocorr2d_14 = 14,
        autocorr2d_15 = 15,
        autocorr2d_16 = 16,
        autocorr2d_17 = 17,
        autocorr2d_18 = 18,
        autocorr2d_19 = 19,
        autocorr2d_20 = 20,
        autocorr2d_21 = 21,
        autocorr2d_22 = 22,
        autocorr2d_23 = 23,
        autocorr2d_24 = 24,
        autocorr2d_25 = 25,
        autocorr2d_26 = 26,
        autocorr2d_27 = 27,
        autocorr2d_28 = 28,
        autocorr2d_29 = 29,
        autocorr2d_30 = 30,
        autocorr2d_31 = 31,
        autocorr2d_32 = 32,
        autocorr2d_33 = 33,
        autocorr2d_34 = 34,
        autocorr2d_35 = 35,
        autocorr2d_36 = 36,
        autocorr2d_37 = 37,
        autocorr2d_38 = 38,
        autocorr2d_39 = 39,
        autocorr2d_40 = 40,
        autocorr2d_41 = 41,
        autocorr2d_42 = 42,
        autocorr2d_43 = 43,
        autocorr2d_44 = 44,
        autocorr2d_45 = 45,
        autocorr2d_46 = 46,
        autocorr2d_47 = 47,
        autocorr2d_48 = 48,
        autocorr2d_49 = 49,
        autocorr2d_50 = 50,
        autocorr2d_51 = 51,
        autocorr2d_52 = 52,
        autocorr2d_53 = 53,
        autocorr2d_54 = 54,
        autocorr2d_55 = 55,
        autocorr2d_56 = 56,
        autocorr2d_57 = 57,
        autocorr2d_58 = 58,
        autocorr2d_59 = 59,
        autocorr2d_60 = 60,
        autocorr2d_61 = 61,
        autocorr2d_62 = 62,
        autocorr2d_63 = 63,
        autocorr2d_64 = 64,
        autocorr2d_65 = 65,
        autocorr2d_66 = 66,
        autocorr2d_67 = 67,
        autocorr2d_68 = 68,
        autocorr2d_69 = 69,
        autocorr2d_70 = 70,
        autocorr2d_71 = 71,
        autocorr2d_72 = 72,
        autocorr2d_73 = 73,
        autocorr2d_74 = 74,
        autocorr2d_75 = 75,
        autocorr2d_76 = 76,
        autocorr2d_77 = 77,
        autocorr2d_78 = 78,
        autocorr2d_79 = 79,
        autocorr2d_80 = 80,
        autocorr2d_81 = 81,
        autocorr2d_82 = 82,
        autocorr2d_83 = 83,
        autocorr2d_84 = 84,
        autocorr2d_85 = 85,
        autocorr2d_86 = 86,
        autocorr2d_87 = 87,
        autocorr2d_88 = 88,
        autocorr2d_89 = 89,
        autocorr2d_90 = 90,
        autocorr2d_91 = 91,
        autocorr2d_92 = 92,
        autocorr2d_93 = 93,
        autocorr2d_94 = 94,
        autocorr2d_95 = 95,
        autocorr2d_96 = 96,
        autocorr2d_97 = 97,
        autocorr2d_98 = 98,
        autocorr2d_99 = 99,
        autocorr2d_100 = 100,
        autocorr2d_101 = 101,
        autocorr2d_102 = 102,
        autocorr2d_103 = 103,
        autocorr2d_104 = 104,
        autocorr2d_105 = 105,
        autocorr2d_106 = 106,
        autocorr2d_107 = 107,
        autocorr2d_108 = 108,
        autocorr2d_109 = 109,
        autocorr2d_110 = 110,
        autocorr2d_111 = 111,
        autocorr2d_112 = 112,
        autocorr2d_113 = 113,
        autocorr2d_114 = 114,
        autocorr2d_115 = 115,
        autocorr2d_116 = 116,
        autocorr2d_117 = 117,
        autocorr2d_118 = 118,
        autocorr2d_119 = 119,
        autocorr2d_120 = 120,
        autocorr2d_121 = 121,
        autocorr2d_122 = 122,
        autocorr2d_123 = 123,
        autocorr2d_124 = 124,
        autocorr2d_125 = 125,
        autocorr2d_126 = 126,
        autocorr2d_127 = 127,
        autocorr2d_128 = 128,
        autocorr2d_129 = 129,
        autocorr2d_130 = 130,
        autocorr2d_131 = 131,
        autocorr2d_132 = 132,
        autocorr2d_133 = 133,
        autocorr2d_134 = 134,
        autocorr2d_135 = 135,
        autocorr2d_136 = 136,
        autocorr2d_137 = 137,
        autocorr2d_138 = 138,
        autocorr2d_139 = 139,
        autocorr2d_140 = 140,
        autocorr2d_141 = 141,
        autocorr2d_142 = 142,
        autocorr2d_143 = 143,
        autocorr2d_144 = 144,
        autocorr2d_145 = 145,
        autocorr2d_146 = 146,
        autocorr2d_147 = 147,
        autocorr2d_148 = 148,
        autocorr2d_149 = 149,
        autocorr2d_150 = 150,
        autocorr2d_151 = 151,
        autocorr2d_152 = 152,
        autocorr2d_153 = 153,
        autocorr2d_154 = 154,
        autocorr2d_155 = 155,
        autocorr2d_156 = 156,
        autocorr2d_157 = 157,
        autocorr2d_158 = 158,
        autocorr2d_159 = 159,
        autocorr2d_160 = 160,
        autocorr2d_161 = 161,
        autocorr2d_162 = 162,
        autocorr2d_163 = 163,
        autocorr2d_164 = 164,
        autocorr2d_165 = 165,
        autocorr2d_166 = 166,
        autocorr2d_167 = 167,
        autocorr2d_168 = 168,
        autocorr2d_169 = 169,
        autocorr2d_170 = 170,
        autocorr2d_171 = 171,
        autocorr2d_172 = 172,
        autocorr2d_173 = 173,
        autocorr2d_174 = 174,
        autocorr2d_175 = 175,
        autocorr2d_176 = 176,
        autocorr2d_177 = 177,
        autocorr2d_178 = 178,
        autocorr2d_179 = 179,
        autocorr2d_180 = 180,
        autocorr2d_181 = 181,
        autocorr2d_182 = 182,
        autocorr2d_183 = 183,
        autocorr2d_184 = 184,
        autocorr2d_185 = 185,
        autocorr2d_186 = 186,
        autocorr2d_187 = 187,
        autocorr2d_188 = 188,
        autocorr2d_189 = 189,
        autocorr2d_190 = 190,
        autocorr2d_191 = 191,
    }
}

/// A descriptor: a function of a molecule. Python's descriptors return
/// `int`s or `float`s; here every value is an `f64`.
pub type DescriptorFn = fn(&Mol) -> Result<f64>;

/// A descriptor's Python name, its version (Python's `version` attribute),
/// the summary `MolecularDescriptorCalculator.GetDescriptorSummaries`
/// gives it, and its function.
pub(crate) type DescriptorEntry = (&'static str, &'static str, &'static str, DescriptorFn);

/// Python's `descList`, with the 192 `AUTOCORR2D_*` descriptors after it,
/// which Python defines in the module but lists only after
/// `setupAUTOCorrDescriptors()`. The summaries are those of conda-forge's
/// RDKit 2024.03.6, whose Boost.Python wrappers give the functions Lipinski
/// takes from `rdMolDescriptors` their signatures as docstrings.
pub(crate) static DESCRIPTORS: [DescriptorEntry; 402] = [
    (
        "MaxAbsEStateIndex",
        "1.0.0",
        "N/A",
        e_state::max_abs_e_state_index,
    ),
    ("MaxEStateIndex", "1.0.0", "N/A", e_state::max_e_state_index),
    (
        "MinAbsEStateIndex",
        "1.0.0",
        "N/A",
        e_state::min_abs_e_state_index,
    ),
    ("MinEStateIndex", "1.0.0", "N/A", e_state::min_e_state_index),
    (
        "qed",
        "1.1.0",
        "Calculate the weighted sum of ADS mapped properties",
        qed::qed,
    ),
    (
        "SPS",
        "1.0.0",
        "Calculates the SpacialScore descriptor. By default, the score is normalized by the number of heavy atoms (nSPS) resulting in a float value, otherwise (normalize=False) the absolute score is returned as an integer.",
        spacial_score::sps,
    ),
    (
        "MolWt",
        "1.0.0",
        "The average molecular weight of the molecule",
        mol_wt,
    ),
    (
        "HeavyAtomMolWt",
        "1.0.0",
        "The average molecular weight of the molecule ignoring hydrogens",
        heavy_atom_mol_wt,
    ),
    (
        "ExactMolWt",
        "1.0.0",
        "The exact molecular weight of the molecule",
        exact_mol_wt,
    ),
    (
        "NumValenceElectrons",
        "1.1.0",
        "The number of valence electrons the molecule has",
        |m| Ok(num_valence_electrons(m)? as f64),
    ),
    (
        "NumRadicalElectrons",
        "1.1.0",
        "The number of radical electrons the molecule has (says nothing about spin state)",
        |m| Ok(num_radical_electrons(m) as f64),
    ),
    ("MaxPartialCharge", "1.0.0", "N/A", max_partial_charge),
    ("MinPartialCharge", "1.0.0", "N/A", min_partial_charge),
    (
        "MaxAbsPartialCharge",
        "1.0.0",
        "N/A",
        max_abs_partial_charge,
    ),
    (
        "MinAbsPartialCharge",
        "1.0.0",
        "N/A",
        min_abs_partial_charge,
    ),
    ("FpDensityMorgan1", "1.0.0", "N/A", fp_density_morgan1),
    ("FpDensityMorgan2", "1.0.0", "N/A", fp_density_morgan2),
    ("FpDensityMorgan3", "1.0.0", "N/A", fp_density_morgan3),
    ("BCUT2D_MWHI", "1.0.0", "N/A", vector::bcut2d_0),
    ("BCUT2D_MWLOW", "1.0.0", "N/A", vector::bcut2d_1),
    ("BCUT2D_CHGHI", "1.0.0", "N/A", vector::bcut2d_2),
    ("BCUT2D_CHGLO", "1.0.0", "N/A", vector::bcut2d_3),
    ("BCUT2D_LOGPHI", "1.0.0", "N/A", vector::bcut2d_4),
    ("BCUT2D_LOGPLOW", "1.0.0", "N/A", vector::bcut2d_5),
    ("BCUT2D_MRHI", "1.0.0", "N/A", vector::bcut2d_6),
    ("BCUT2D_MRLOW", "1.0.0", "N/A", vector::bcut2d_7),
    (
        "AvgIpc",
        "1.0.0",
        "This returns the average information content of the coefficients of the characteristic polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.",
        graph_descriptors::avg_ipc,
    ),
    (
        "BalabanJ",
        "1.0.0",
        "Calculate Balaban's J value for a molecule",
        graph_descriptors::balaban_j,
    ),
    (
        "BertzCT",
        "2.0.0",
        "A topological index meant to quantify \"complexity\" of molecules.",
        graph_descriptors::bertz_ct,
    ),
    (
        "Chi0",
        "1.0.0",
        "From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)",
        |m| Ok(graph_descriptors::chi0(m)),
    ),
    ("Chi0n", "1.2.0", "N/A", graph_descriptors::chi0n),
    ("Chi0v", "1.2.0", "N/A", graph_descriptors::chi0v),
    (
        "Chi1",
        "1.0.0",
        "From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)",
        |m| Ok(graph_descriptors::chi1(m)),
    ),
    ("Chi1n", "1.2.0", "N/A", graph_descriptors::chi1n),
    ("Chi1v", "1.2.0", "N/A", graph_descriptors::chi1v),
    ("Chi2n", "1.2.0", "N/A", graph_descriptors::chi2n),
    ("Chi2v", "1.2.0", "N/A", graph_descriptors::chi2v),
    ("Chi3n", "1.2.0", "N/A", graph_descriptors::chi3n),
    ("Chi3v", "1.2.0", "N/A", graph_descriptors::chi3v),
    ("Chi4n", "1.2.0", "N/A", graph_descriptors::chi4n),
    ("Chi4v", "1.2.0", "N/A", graph_descriptors::chi4v),
    (
        "HallKierAlpha",
        "1.2.0",
        "N/A",
        graph_descriptors::hall_kier_alpha,
    ),
    (
        "Ipc",
        "1.0.0",
        "This returns the information content of the coefficients of the characteristic polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.",
        graph_descriptors::ipc,
    ),
    ("Kappa1", "1.1.0", "N/A", graph_descriptors::kappa1),
    ("Kappa2", "1.1.0", "N/A", graph_descriptors::kappa2),
    ("Kappa3", "1.1.0", "N/A", graph_descriptors::kappa3),
    ("LabuteASA", "1.0.2", "N/A", mol_surf::labute_asa),
    (
        "PEOE_VSA1",
        "1.0.1",
        "MOE Charge VSA Descriptor 1 (-inf < x < -0.30)",
        mol_surf::peoe_vsa1,
    ),
    (
        "PEOE_VSA10",
        "1.0.1",
        "MOE Charge VSA Descriptor 10 ( 0.10 <= x <  0.15)",
        mol_surf::peoe_vsa10,
    ),
    (
        "PEOE_VSA11",
        "1.0.1",
        "MOE Charge VSA Descriptor 11 ( 0.15 <= x <  0.20)",
        mol_surf::peoe_vsa11,
    ),
    (
        "PEOE_VSA12",
        "1.0.1",
        "MOE Charge VSA Descriptor 12 ( 0.20 <= x <  0.25)",
        mol_surf::peoe_vsa12,
    ),
    (
        "PEOE_VSA13",
        "1.0.1",
        "MOE Charge VSA Descriptor 13 ( 0.25 <= x <  0.30)",
        mol_surf::peoe_vsa13,
    ),
    (
        "PEOE_VSA14",
        "1.0.1",
        "MOE Charge VSA Descriptor 14 ( 0.30 <= x < inf)",
        mol_surf::peoe_vsa14,
    ),
    (
        "PEOE_VSA2",
        "1.0.1",
        "MOE Charge VSA Descriptor 2 (-0.30 <= x < -0.25)",
        mol_surf::peoe_vsa2,
    ),
    (
        "PEOE_VSA3",
        "1.0.1",
        "MOE Charge VSA Descriptor 3 (-0.25 <= x < -0.20)",
        mol_surf::peoe_vsa3,
    ),
    (
        "PEOE_VSA4",
        "1.0.1",
        "MOE Charge VSA Descriptor 4 (-0.20 <= x < -0.15)",
        mol_surf::peoe_vsa4,
    ),
    (
        "PEOE_VSA5",
        "1.0.1",
        "MOE Charge VSA Descriptor 5 (-0.15 <= x < -0.10)",
        mol_surf::peoe_vsa5,
    ),
    (
        "PEOE_VSA6",
        "1.0.1",
        "MOE Charge VSA Descriptor 6 (-0.10 <= x < -0.05)",
        mol_surf::peoe_vsa6,
    ),
    (
        "PEOE_VSA7",
        "1.0.1",
        "MOE Charge VSA Descriptor 7 (-0.05 <= x <  0.00)",
        mol_surf::peoe_vsa7,
    ),
    (
        "PEOE_VSA8",
        "1.0.1",
        "MOE Charge VSA Descriptor 8 ( 0.00 <= x <  0.05)",
        mol_surf::peoe_vsa8,
    ),
    (
        "PEOE_VSA9",
        "1.0.1",
        "MOE Charge VSA Descriptor 9 ( 0.05 <= x <  0.10)",
        mol_surf::peoe_vsa9,
    ),
    (
        "SMR_VSA1",
        "1.0.1",
        "MOE MR VSA Descriptor 1 (-inf < x <  1.29)",
        mol_surf::smr_vsa1,
    ),
    (
        "SMR_VSA10",
        "1.0.1",
        "MOE MR VSA Descriptor 10 ( 4.00 <= x < inf)",
        mol_surf::smr_vsa10,
    ),
    (
        "SMR_VSA2",
        "1.0.1",
        "MOE MR VSA Descriptor 2 ( 1.29 <= x <  1.82)",
        mol_surf::smr_vsa2,
    ),
    (
        "SMR_VSA3",
        "1.0.1",
        "MOE MR VSA Descriptor 3 ( 1.82 <= x <  2.24)",
        mol_surf::smr_vsa3,
    ),
    (
        "SMR_VSA4",
        "1.0.1",
        "MOE MR VSA Descriptor 4 ( 2.24 <= x <  2.45)",
        mol_surf::smr_vsa4,
    ),
    (
        "SMR_VSA5",
        "1.0.1",
        "MOE MR VSA Descriptor 5 ( 2.45 <= x <  2.75)",
        mol_surf::smr_vsa5,
    ),
    (
        "SMR_VSA6",
        "1.0.1",
        "MOE MR VSA Descriptor 6 ( 2.75 <= x <  3.05)",
        mol_surf::smr_vsa6,
    ),
    (
        "SMR_VSA7",
        "1.0.1",
        "MOE MR VSA Descriptor 7 ( 3.05 <= x <  3.63)",
        mol_surf::smr_vsa7,
    ),
    (
        "SMR_VSA8",
        "1.0.1",
        "MOE MR VSA Descriptor 8 ( 3.63 <= x <  3.80)",
        mol_surf::smr_vsa8,
    ),
    (
        "SMR_VSA9",
        "1.0.1",
        "MOE MR VSA Descriptor 9 ( 3.80 <= x <  4.00)",
        mol_surf::smr_vsa9,
    ),
    (
        "SlogP_VSA1",
        "1.0.1",
        "MOE logP VSA Descriptor 1 (-inf < x < -0.40)",
        mol_surf::slogp_vsa1,
    ),
    (
        "SlogP_VSA10",
        "1.0.1",
        "MOE logP VSA Descriptor 10 ( 0.40 <= x <  0.50)",
        mol_surf::slogp_vsa10,
    ),
    (
        "SlogP_VSA11",
        "1.0.1",
        "MOE logP VSA Descriptor 11 ( 0.50 <= x <  0.60)",
        mol_surf::slogp_vsa11,
    ),
    (
        "SlogP_VSA12",
        "1.0.1",
        "MOE logP VSA Descriptor 12 ( 0.60 <= x < inf)",
        mol_surf::slogp_vsa12,
    ),
    (
        "SlogP_VSA2",
        "1.0.1",
        "MOE logP VSA Descriptor 2 (-0.40 <= x < -0.20)",
        mol_surf::slogp_vsa2,
    ),
    (
        "SlogP_VSA3",
        "1.0.1",
        "MOE logP VSA Descriptor 3 (-0.20 <= x <  0.00)",
        mol_surf::slogp_vsa3,
    ),
    (
        "SlogP_VSA4",
        "1.0.1",
        "MOE logP VSA Descriptor 4 ( 0.00 <= x <  0.10)",
        mol_surf::slogp_vsa4,
    ),
    (
        "SlogP_VSA5",
        "1.0.1",
        "MOE logP VSA Descriptor 5 ( 0.10 <= x <  0.15)",
        mol_surf::slogp_vsa5,
    ),
    (
        "SlogP_VSA6",
        "1.0.1",
        "MOE logP VSA Descriptor 6 ( 0.15 <= x <  0.20)",
        mol_surf::slogp_vsa6,
    ),
    (
        "SlogP_VSA7",
        "1.0.1",
        "MOE logP VSA Descriptor 7 ( 0.20 <= x <  0.25)",
        mol_surf::slogp_vsa7,
    ),
    (
        "SlogP_VSA8",
        "1.0.1",
        "MOE logP VSA Descriptor 8 ( 0.25 <= x <  0.30)",
        mol_surf::slogp_vsa8,
    ),
    (
        "SlogP_VSA9",
        "1.0.1",
        "MOE logP VSA Descriptor 9 ( 0.30 <= x <  0.40)",
        mol_surf::slogp_vsa9,
    ),
    ("TPSA", "2.0.0", "N/A", mol_surf::tpsa),
    (
        "EState_VSA1",
        "1.0.1",
        "EState VSA Descriptor 1 (-inf < x <  -0.39)",
        e_state_vsa::e_state_vsa1,
    ),
    (
        "EState_VSA10",
        "1.0.1",
        "EState VSA Descriptor 10 ( 9.17 <= x <  15.00)",
        e_state_vsa::e_state_vsa10,
    ),
    (
        "EState_VSA11",
        "1.0.1",
        "EState VSA Descriptor 11 ( 15.00 <= x < inf)",
        e_state_vsa::e_state_vsa11,
    ),
    (
        "EState_VSA2",
        "1.0.1",
        "EState VSA Descriptor 2 ( -0.39 <= x <  0.29)",
        e_state_vsa::e_state_vsa2,
    ),
    (
        "EState_VSA3",
        "1.0.1",
        "EState VSA Descriptor 3 ( 0.29 <= x <  0.72)",
        e_state_vsa::e_state_vsa3,
    ),
    (
        "EState_VSA4",
        "1.0.1",
        "EState VSA Descriptor 4 ( 0.72 <= x <  1.17)",
        e_state_vsa::e_state_vsa4,
    ),
    (
        "EState_VSA5",
        "1.0.1",
        "EState VSA Descriptor 5 ( 1.17 <= x <  1.54)",
        e_state_vsa::e_state_vsa5,
    ),
    (
        "EState_VSA6",
        "1.0.1",
        "EState VSA Descriptor 6 ( 1.54 <= x <  1.81)",
        e_state_vsa::e_state_vsa6,
    ),
    (
        "EState_VSA7",
        "1.0.1",
        "EState VSA Descriptor 7 ( 1.81 <= x <  2.05)",
        e_state_vsa::e_state_vsa7,
    ),
    (
        "EState_VSA8",
        "1.0.1",
        "EState VSA Descriptor 8 ( 2.05 <= x <  4.69)",
        e_state_vsa::e_state_vsa8,
    ),
    (
        "EState_VSA9",
        "1.0.1",
        "EState VSA Descriptor 9 ( 4.69 <= x <  9.17)",
        e_state_vsa::e_state_vsa9,
    ),
    (
        "VSA_EState1",
        "1.0.0",
        "VSA EState Descriptor 1 (-inf < x <  4.78)",
        e_state_vsa::vsa_e_state1,
    ),
    (
        "VSA_EState10",
        "1.0.0",
        "VSA EState Descriptor 10 ( 11.00 <= x < inf)",
        e_state_vsa::vsa_e_state10,
    ),
    (
        "VSA_EState2",
        "1.0.0",
        "VSA EState Descriptor 2 ( 4.78 <= x <  5.00)",
        e_state_vsa::vsa_e_state2,
    ),
    (
        "VSA_EState3",
        "1.0.0",
        "VSA EState Descriptor 3 ( 5.00 <= x <  5.41)",
        e_state_vsa::vsa_e_state3,
    ),
    (
        "VSA_EState4",
        "1.0.0",
        "VSA EState Descriptor 4 ( 5.41 <= x <  5.74)",
        e_state_vsa::vsa_e_state4,
    ),
    (
        "VSA_EState5",
        "1.0.0",
        "VSA EState Descriptor 5 ( 5.74 <= x <  6.00)",
        e_state_vsa::vsa_e_state5,
    ),
    (
        "VSA_EState6",
        "1.0.0",
        "VSA EState Descriptor 6 ( 6.00 <= x <  6.07)",
        e_state_vsa::vsa_e_state6,
    ),
    (
        "VSA_EState7",
        "1.0.0",
        "VSA EState Descriptor 7 ( 6.07 <= x <  6.45)",
        e_state_vsa::vsa_e_state7,
    ),
    (
        "VSA_EState8",
        "1.0.0",
        "VSA EState Descriptor 8 ( 6.45 <= x <  7.00)",
        e_state_vsa::vsa_e_state8,
    ),
    (
        "VSA_EState9",
        "1.0.0",
        "VSA EState Descriptor 9 ( 7.00 <= x <  11.00)",
        e_state_vsa::vsa_e_state9,
    ),
    (
        "FractionCSP3",
        "1.0.0",
        "CalcFractionCSP3( (Mol)mol) -> float : returns the fraction of C atoms that are SP3 hybridized",
        lipinski::fraction_csp3,
    ),
    (
        "HeavyAtomCount",
        "1.0.1",
        "Number of heavy atoms a molecule.",
        |m| Ok(f64::from(lipinski::heavy_atom_count(m))),
    ),
    ("NHOHCount", "2.0.0", "Number of NHs or OHs", |m| {
        Ok(f64::from(lipinski::nhoh_count(m)?))
    }),
    ("NOCount", "1.0.0", "Number of Nitrogens and Oxygens", |m| {
        Ok(f64::from(lipinski::no_count(m)))
    }),
    (
        "NumAliphaticCarbocycles",
        "1.0.0",
        "CalcNumAliphaticCarbocycles( (Mol)mol) -> int : returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule",
        |m| Ok(f64::from(lipinski::num_aliphatic_carbocycles(m))),
    ),
    (
        "NumAliphaticHeterocycles",
        "1.0.0",
        "CalcNumAliphaticHeterocycles( (Mol)mol) -> int : returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule",
        |m| Ok(f64::from(lipinski::num_aliphatic_heterocycles(m))),
    ),
    (
        "NumAliphaticRings",
        "1.0.0",
        "CalcNumAliphaticRings( (Mol)mol) -> int : returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule",
        |m| Ok(f64::from(lipinski::num_aliphatic_rings(m))),
    ),
    (
        "NumAromaticCarbocycles",
        "1.0.0",
        "CalcNumAromaticCarbocycles( (Mol)mol) -> int : returns the number of aromatic carbocycles for a molecule",
        |m| Ok(f64::from(lipinski::num_aromatic_carbocycles(m))),
    ),
    (
        "NumAromaticHeterocycles",
        "1.0.0",
        "CalcNumAromaticHeterocycles( (Mol)mol) -> int : returns the number of aromatic heterocycles for a molecule",
        |m| Ok(f64::from(lipinski::num_aromatic_heterocycles(m))),
    ),
    (
        "NumAromaticRings",
        "1.0.0",
        "CalcNumAromaticRings( (Mol)mol) -> int : returns the number of aromatic rings for a molecule",
        |m| Ok(f64::from(lipinski::num_aromatic_rings(m))),
    ),
    (
        "NumHAcceptors",
        "2.0.0",
        "Number of Hydrogen Bond Acceptors",
        |m| Ok(f64::from(lipinski::num_h_acceptors(m)?)),
    ),
    (
        "NumHDonors",
        "1.0.0",
        "Number of Hydrogen Bond Donors",
        |m| Ok(f64::from(lipinski::num_h_donors(m)?)),
    ),
    ("NumHeteroatoms", "1.0.0", "Number of Heteroatoms", |m| {
        Ok(f64::from(lipinski::num_heteroatoms(m)?))
    }),
    (
        "NumRotatableBonds",
        "1.0.0",
        "Number of Rotatable Bonds",
        |m| Ok(f64::from(lipinski::num_rotatable_bonds(m)?)),
    ),
    (
        "NumSaturatedCarbocycles",
        "1.0.0",
        "CalcNumSaturatedCarbocycles( (Mol)mol) -> int : returns the number of saturated carbocycles for a molecule",
        |m| Ok(f64::from(lipinski::num_saturated_carbocycles(m))),
    ),
    (
        "NumSaturatedHeterocycles",
        "1.0.0",
        "CalcNumSaturatedHeterocycles( (Mol)mol) -> int : returns the number of saturated heterocycles for a molecule",
        |m| Ok(f64::from(lipinski::num_saturated_heterocycles(m))),
    ),
    (
        "NumSaturatedRings",
        "1.0.0",
        "CalcNumSaturatedRings( (Mol)mol) -> int : returns the number of saturated rings for a molecule",
        |m| Ok(f64::from(lipinski::num_saturated_rings(m))),
    ),
    ("RingCount", "1.0.0", "N/A", |m| {
        Ok(f64::from(lipinski::ring_count(m)?))
    }),
    (
        "MolLogP",
        "1.2.0",
        "Wildman-Crippen LogP value",
        crippen::mol_log_p,
    ),
    (
        "MolMR",
        "1.2.0",
        "Wildman-Crippen MR value",
        crippen::mol_mr,
    ),
    (
        "fr_Al_COO",
        "N/A",
        "Number of aliphatic carboxylic acids",
        |m| Ok(f64::from(fragments::fr_al_coo(m)?)),
    ),
    (
        "fr_Al_OH",
        "N/A",
        "Number of aliphatic hydroxyl groups",
        |m| Ok(f64::from(fragments::fr_al_oh(m)?)),
    ),
    (
        "fr_Al_OH_noTert",
        "N/A",
        "Number of aliphatic hydroxyl groups excluding tert-OH",
        |m| Ok(f64::from(fragments::fr_al_oh_no_tert(m)?)),
    ),
    (
        "fr_ArN",
        "N/A",
        "Number of N functional groups attached to aromatics",
        |m| Ok(f64::from(fragments::fr_arn(m)?)),
    ),
    (
        "fr_Ar_COO",
        "N/A",
        "Number of Aromatic carboxylic acide",
        |m| Ok(f64::from(fragments::fr_ar_coo(m)?)),
    ),
    ("fr_Ar_N", "N/A", "Number of aromatic nitrogens", |m| {
        Ok(f64::from(fragments::fr_ar_n(m)?))
    }),
    ("fr_Ar_NH", "N/A", "Number of aromatic amines", |m| {
        Ok(f64::from(fragments::fr_ar_nh(m)?))
    }),
    (
        "fr_Ar_OH",
        "N/A",
        "Number of aromatic hydroxyl groups",
        |m| Ok(f64::from(fragments::fr_ar_oh(m)?)),
    ),
    ("fr_COO", "N/A", "Number of carboxylic acids", |m| {
        Ok(f64::from(fragments::fr_coo(m)?))
    }),
    ("fr_COO2", "N/A", "Number of carboxylic acids", |m| {
        Ok(f64::from(fragments::fr_coo2(m)?))
    }),
    ("fr_C_O", "N/A", "Number of carbonyl O", |m| {
        Ok(f64::from(fragments::fr_c_o(m)?))
    }),
    (
        "fr_C_O_noCOO",
        "N/A",
        "Number of carbonyl O, excluding COOH",
        |m| Ok(f64::from(fragments::fr_c_o_no_coo(m)?)),
    ),
    ("fr_C_S", "N/A", "Number of thiocarbonyl", |m| {
        Ok(f64::from(fragments::fr_c_s(m)?))
    }),
    (
        "fr_HOCCN",
        "N/A",
        "Number of C(OH)CCN-Ctert-alkyl or  C(OH)CCNcyclic",
        |m| Ok(f64::from(fragments::fr_hoccn(m)?)),
    ),
    ("fr_Imine", "N/A", "Number of Imines", |m| {
        Ok(f64::from(fragments::fr_imine(m)?))
    }),
    ("fr_NH0", "N/A", "Number of Tertiary amines", |m| {
        Ok(f64::from(fragments::fr_nh0(m)?))
    }),
    ("fr_NH1", "N/A", "Number of Secondary amines", |m| {
        Ok(f64::from(fragments::fr_nh1(m)?))
    }),
    ("fr_NH2", "N/A", "Number of Primary amines", |m| {
        Ok(f64::from(fragments::fr_nh2(m)?))
    }),
    ("fr_N_O", "N/A", "Number of hydroxylamine groups", |m| {
        Ok(f64::from(fragments::fr_n_o(m)?))
    }),
    ("fr_Ndealkylation1", "N/A", "Number of XCCNR groups", |m| {
        Ok(f64::from(fragments::fr_ndealkylation1(m)?))
    }),
    (
        "fr_Ndealkylation2",
        "N/A",
        "Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)",
        |m| Ok(f64::from(fragments::fr_ndealkylation2(m)?)),
    ),
    (
        "fr_Nhpyrrole",
        "N/A",
        "Number of H-pyrrole nitrogens",
        |m| Ok(f64::from(fragments::fr_nhpyrrole(m)?)),
    ),
    ("fr_SH", "N/A", "Number of thiol groups", |m| {
        Ok(f64::from(fragments::fr_sh(m)?))
    }),
    ("fr_aldehyde", "N/A", "Number of aldehydes", |m| {
        Ok(f64::from(fragments::fr_aldehyde(m)?))
    }),
    (
        "fr_alkyl_carbamate",
        "N/A",
        "Number of alkyl carbamates (subject to hydrolysis)",
        |m| Ok(f64::from(fragments::fr_alkyl_carbamate(m)?)),
    ),
    ("fr_alkyl_halide", "N/A", "Number of alkyl halides", |m| {
        Ok(f64::from(fragments::fr_alkyl_halide(m)?))
    }),
    (
        "fr_allylic_oxid",
        "N/A",
        "Number of allylic oxidation sites excluding steroid dienone",
        |m| Ok(f64::from(fragments::fr_allylic_oxid(m)?)),
    ),
    ("fr_amide", "N/A", "Number of amides", |m| {
        Ok(f64::from(fragments::fr_amide(m)?))
    }),
    ("fr_amidine", "N/A", "Number of amidine groups", |m| {
        Ok(f64::from(fragments::fr_amidine(m)?))
    }),
    ("fr_aniline", "N/A", "Number of anilines", |m| {
        Ok(f64::from(fragments::fr_aniline(m)?))
    }),
    (
        "fr_aryl_methyl",
        "N/A",
        "Number of aryl methyl sites for hydroxylation",
        |m| Ok(f64::from(fragments::fr_aryl_methyl(m)?)),
    ),
    ("fr_azide", "N/A", "Number of azide groups", |m| {
        Ok(f64::from(fragments::fr_azide(m)?))
    }),
    ("fr_azo", "N/A", "Number of azo groups", |m| {
        Ok(f64::from(fragments::fr_azo(m)?))
    }),
    ("fr_barbitur", "N/A", "Number of barbiturate groups", |m| {
        Ok(f64::from(fragments::fr_barbitur(m)?))
    }),
    ("fr_benzene", "N/A", "Number of benzene rings", |m| {
        Ok(f64::from(fragments::fr_benzene(m)?))
    }),
    (
        "fr_benzodiazepine",
        "N/A",
        "Number of benzodiazepines with no additional fused rings",
        |m| Ok(f64::from(fragments::fr_benzodiazepine(m)?)),
    ),
    ("fr_bicyclic", "N/A", "Bicyclic", |m| {
        Ok(f64::from(fragments::fr_bicyclic(m)?))
    }),
    ("fr_diazo", "N/A", "Number of diazo groups", |m| {
        Ok(f64::from(fragments::fr_diazo(m)?))
    }),
    (
        "fr_dihydropyridine",
        "N/A",
        "Number of dihydropyridines",
        |m| Ok(f64::from(fragments::fr_dihydropyridine(m)?)),
    ),
    ("fr_epoxide", "N/A", "Number of epoxide rings", |m| {
        Ok(f64::from(fragments::fr_epoxide(m)?))
    }),
    ("fr_ester", "N/A", "Number of esters", |m| {
        Ok(f64::from(fragments::fr_ester(m)?))
    }),
    (
        "fr_ether",
        "N/A",
        "Number of ether oxygens (including phenoxy)",
        |m| Ok(f64::from(fragments::fr_ether(m)?)),
    ),
    ("fr_furan", "N/A", "Number of furan rings", |m| {
        Ok(f64::from(fragments::fr_furan(m)?))
    }),
    ("fr_guanido", "N/A", "Number of guanidine groups", |m| {
        Ok(f64::from(fragments::fr_guanido(m)?))
    }),
    ("fr_halogen", "N/A", "Number of halogens", |m| {
        Ok(f64::from(fragments::fr_halogen(m)?))
    }),
    ("fr_hdrzine", "N/A", "Number of hydrazine groups", |m| {
        Ok(f64::from(fragments::fr_hdrzine(m)?))
    }),
    ("fr_hdrzone", "N/A", "Number of hydrazone groups", |m| {
        Ok(f64::from(fragments::fr_hdrzone(m)?))
    }),
    ("fr_imidazole", "N/A", "Number of imidazole rings", |m| {
        Ok(f64::from(fragments::fr_imidazole(m)?))
    }),
    ("fr_imide", "N/A", "Number of imide groups", |m| {
        Ok(f64::from(fragments::fr_imide(m)?))
    }),
    ("fr_isocyan", "N/A", "Number of isocyanates", |m| {
        Ok(f64::from(fragments::fr_isocyan(m)?))
    }),
    ("fr_isothiocyan", "N/A", "Number of isothiocyanates", |m| {
        Ok(f64::from(fragments::fr_isothiocyan(m)?))
    }),
    ("fr_ketone", "N/A", "Number of ketones", |m| {
        Ok(f64::from(fragments::fr_ketone(m)?))
    }),
    (
        "fr_ketone_Topliss",
        "N/A",
        "Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha",
        |m| Ok(f64::from(fragments::fr_ketone_topliss(m)?)),
    ),
    ("fr_lactam", "N/A", "Number of beta lactams", |m| {
        Ok(f64::from(fragments::fr_lactam(m)?))
    }),
    (
        "fr_lactone",
        "N/A",
        "Number of cyclic esters (lactones)",
        |m| Ok(f64::from(fragments::fr_lactone(m)?)),
    ),
    ("fr_methoxy", "N/A", "Number of methoxy groups -OCH3", |m| {
        Ok(f64::from(fragments::fr_methoxy(m)?))
    }),
    ("fr_morpholine", "N/A", "Number of morpholine rings", |m| {
        Ok(f64::from(fragments::fr_morpholine(m)?))
    }),
    ("fr_nitrile", "N/A", "Number of nitriles", |m| {
        Ok(f64::from(fragments::fr_nitrile(m)?))
    }),
    ("fr_nitro", "N/A", "Number of nitro groups", |m| {
        Ok(f64::from(fragments::fr_nitro(m)?))
    }),
    (
        "fr_nitro_arom",
        "N/A",
        "Number of nitro benzene ring substituents",
        |m| Ok(f64::from(fragments::fr_nitro_arom(m)?)),
    ),
    (
        "fr_nitro_arom_nonortho",
        "N/A",
        "Number of non-ortho nitro benzene ring substituents",
        |m| Ok(f64::from(fragments::fr_nitro_arom_nonortho(m)?)),
    ),
    (
        "fr_nitroso",
        "N/A",
        "Number of nitroso groups, excluding NO2",
        |m| Ok(f64::from(fragments::fr_nitroso(m)?)),
    ),
    ("fr_oxazole", "N/A", "Number of oxazole rings", |m| {
        Ok(f64::from(fragments::fr_oxazole(m)?))
    }),
    ("fr_oxime", "N/A", "Number of oxime groups", |m| {
        Ok(f64::from(fragments::fr_oxime(m)?))
    }),
    (
        "fr_para_hydroxylation",
        "N/A",
        "Number of para-hydroxylation sites",
        |m| Ok(f64::from(fragments::fr_para_hydroxylation(m)?)),
    ),
    ("fr_phenol", "N/A", "Number of phenols", |m| {
        Ok(f64::from(fragments::fr_phenol(m)?))
    }),
    (
        "fr_phenol_noOrthoHbond",
        "N/A",
        "Number of phenolic OH excluding ortho intramolecular Hbond substituents",
        |m| Ok(f64::from(fragments::fr_phenol_no_ortho_hbond(m)?)),
    ),
    (
        "fr_phos_acid",
        "N/A",
        "Number of phosphoric acid groups",
        |m| Ok(f64::from(fragments::fr_phos_acid(m)?)),
    ),
    (
        "fr_phos_ester",
        "N/A",
        "Number of phosphoric ester groups",
        |m| Ok(f64::from(fragments::fr_phos_ester(m)?)),
    ),
    ("fr_piperdine", "N/A", "Number of piperdine rings", |m| {
        Ok(f64::from(fragments::fr_piperdine(m)?))
    }),
    ("fr_piperzine", "N/A", "Number of piperzine rings", |m| {
        Ok(f64::from(fragments::fr_piperzine(m)?))
    }),
    ("fr_priamide", "N/A", "Number of primary amides", |m| {
        Ok(f64::from(fragments::fr_priamide(m)?))
    }),
    (
        "fr_prisulfonamd",
        "N/A",
        "Number of primary sulfonamides",
        |m| Ok(f64::from(fragments::fr_prisulfonamd(m)?)),
    ),
    ("fr_pyridine", "N/A", "Number of pyridine rings", |m| {
        Ok(f64::from(fragments::fr_pyridine(m)?))
    }),
    ("fr_quatN", "N/A", "Number of quaternary nitrogens", |m| {
        Ok(f64::from(fragments::fr_quat_n(m)?))
    }),
    ("fr_sulfide", "N/A", "Number of thioether", |m| {
        Ok(f64::from(fragments::fr_sulfide(m)?))
    }),
    ("fr_sulfonamd", "N/A", "Number of sulfonamides", |m| {
        Ok(f64::from(fragments::fr_sulfonamd(m)?))
    }),
    ("fr_sulfone", "N/A", "Number of sulfone groups", |m| {
        Ok(f64::from(fragments::fr_sulfone(m)?))
    }),
    (
        "fr_term_acetylene",
        "N/A",
        "Number of terminal acetylenes",
        |m| Ok(f64::from(fragments::fr_term_acetylene(m)?)),
    ),
    ("fr_tetrazole", "N/A", "Number of tetrazole rings", |m| {
        Ok(f64::from(fragments::fr_tetrazole(m)?))
    }),
    ("fr_thiazole", "N/A", "Number of thiazole rings", |m| {
        Ok(f64::from(fragments::fr_thiazole(m)?))
    }),
    ("fr_thiocyan", "N/A", "Number of thiocyanates", |m| {
        Ok(f64::from(fragments::fr_thiocyan(m)?))
    }),
    ("fr_thiophene", "N/A", "Number of thiophene rings", |m| {
        Ok(f64::from(fragments::fr_thiophene(m)?))
    }),
    (
        "fr_unbrch_alkane",
        "N/A",
        "Number of unbranched alkanes  of at least 4 members (excludes halogenated alkanes)",
        |m| Ok(f64::from(fragments::fr_unbrch_alkane(m)?)),
    ),
    ("fr_urea", "N/A", "Number of urea groups", |m| {
        Ok(f64::from(fragments::fr_urea(m)?))
    }),
    ("AUTOCORR2D_1", "1.0.0", "N/A", vector::autocorr2d_0),
    ("AUTOCORR2D_2", "1.0.0", "N/A", vector::autocorr2d_1),
    ("AUTOCORR2D_3", "1.0.0", "N/A", vector::autocorr2d_2),
    ("AUTOCORR2D_4", "1.0.0", "N/A", vector::autocorr2d_3),
    ("AUTOCORR2D_5", "1.0.0", "N/A", vector::autocorr2d_4),
    ("AUTOCORR2D_6", "1.0.0", "N/A", vector::autocorr2d_5),
    ("AUTOCORR2D_7", "1.0.0", "N/A", vector::autocorr2d_6),
    ("AUTOCORR2D_8", "1.0.0", "N/A", vector::autocorr2d_7),
    ("AUTOCORR2D_9", "1.0.0", "N/A", vector::autocorr2d_8),
    ("AUTOCORR2D_10", "1.0.0", "N/A", vector::autocorr2d_9),
    ("AUTOCORR2D_11", "1.0.0", "N/A", vector::autocorr2d_10),
    ("AUTOCORR2D_12", "1.0.0", "N/A", vector::autocorr2d_11),
    ("AUTOCORR2D_13", "1.0.0", "N/A", vector::autocorr2d_12),
    ("AUTOCORR2D_14", "1.0.0", "N/A", vector::autocorr2d_13),
    ("AUTOCORR2D_15", "1.0.0", "N/A", vector::autocorr2d_14),
    ("AUTOCORR2D_16", "1.0.0", "N/A", vector::autocorr2d_15),
    ("AUTOCORR2D_17", "1.0.0", "N/A", vector::autocorr2d_16),
    ("AUTOCORR2D_18", "1.0.0", "N/A", vector::autocorr2d_17),
    ("AUTOCORR2D_19", "1.0.0", "N/A", vector::autocorr2d_18),
    ("AUTOCORR2D_20", "1.0.0", "N/A", vector::autocorr2d_19),
    ("AUTOCORR2D_21", "1.0.0", "N/A", vector::autocorr2d_20),
    ("AUTOCORR2D_22", "1.0.0", "N/A", vector::autocorr2d_21),
    ("AUTOCORR2D_23", "1.0.0", "N/A", vector::autocorr2d_22),
    ("AUTOCORR2D_24", "1.0.0", "N/A", vector::autocorr2d_23),
    ("AUTOCORR2D_25", "1.0.0", "N/A", vector::autocorr2d_24),
    ("AUTOCORR2D_26", "1.0.0", "N/A", vector::autocorr2d_25),
    ("AUTOCORR2D_27", "1.0.0", "N/A", vector::autocorr2d_26),
    ("AUTOCORR2D_28", "1.0.0", "N/A", vector::autocorr2d_27),
    ("AUTOCORR2D_29", "1.0.0", "N/A", vector::autocorr2d_28),
    ("AUTOCORR2D_30", "1.0.0", "N/A", vector::autocorr2d_29),
    ("AUTOCORR2D_31", "1.0.0", "N/A", vector::autocorr2d_30),
    ("AUTOCORR2D_32", "1.0.0", "N/A", vector::autocorr2d_31),
    ("AUTOCORR2D_33", "1.0.0", "N/A", vector::autocorr2d_32),
    ("AUTOCORR2D_34", "1.0.0", "N/A", vector::autocorr2d_33),
    ("AUTOCORR2D_35", "1.0.0", "N/A", vector::autocorr2d_34),
    ("AUTOCORR2D_36", "1.0.0", "N/A", vector::autocorr2d_35),
    ("AUTOCORR2D_37", "1.0.0", "N/A", vector::autocorr2d_36),
    ("AUTOCORR2D_38", "1.0.0", "N/A", vector::autocorr2d_37),
    ("AUTOCORR2D_39", "1.0.0", "N/A", vector::autocorr2d_38),
    ("AUTOCORR2D_40", "1.0.0", "N/A", vector::autocorr2d_39),
    ("AUTOCORR2D_41", "1.0.0", "N/A", vector::autocorr2d_40),
    ("AUTOCORR2D_42", "1.0.0", "N/A", vector::autocorr2d_41),
    ("AUTOCORR2D_43", "1.0.0", "N/A", vector::autocorr2d_42),
    ("AUTOCORR2D_44", "1.0.0", "N/A", vector::autocorr2d_43),
    ("AUTOCORR2D_45", "1.0.0", "N/A", vector::autocorr2d_44),
    ("AUTOCORR2D_46", "1.0.0", "N/A", vector::autocorr2d_45),
    ("AUTOCORR2D_47", "1.0.0", "N/A", vector::autocorr2d_46),
    ("AUTOCORR2D_48", "1.0.0", "N/A", vector::autocorr2d_47),
    ("AUTOCORR2D_49", "1.0.0", "N/A", vector::autocorr2d_48),
    ("AUTOCORR2D_50", "1.0.0", "N/A", vector::autocorr2d_49),
    ("AUTOCORR2D_51", "1.0.0", "N/A", vector::autocorr2d_50),
    ("AUTOCORR2D_52", "1.0.0", "N/A", vector::autocorr2d_51),
    ("AUTOCORR2D_53", "1.0.0", "N/A", vector::autocorr2d_52),
    ("AUTOCORR2D_54", "1.0.0", "N/A", vector::autocorr2d_53),
    ("AUTOCORR2D_55", "1.0.0", "N/A", vector::autocorr2d_54),
    ("AUTOCORR2D_56", "1.0.0", "N/A", vector::autocorr2d_55),
    ("AUTOCORR2D_57", "1.0.0", "N/A", vector::autocorr2d_56),
    ("AUTOCORR2D_58", "1.0.0", "N/A", vector::autocorr2d_57),
    ("AUTOCORR2D_59", "1.0.0", "N/A", vector::autocorr2d_58),
    ("AUTOCORR2D_60", "1.0.0", "N/A", vector::autocorr2d_59),
    ("AUTOCORR2D_61", "1.0.0", "N/A", vector::autocorr2d_60),
    ("AUTOCORR2D_62", "1.0.0", "N/A", vector::autocorr2d_61),
    ("AUTOCORR2D_63", "1.0.0", "N/A", vector::autocorr2d_62),
    ("AUTOCORR2D_64", "1.0.0", "N/A", vector::autocorr2d_63),
    ("AUTOCORR2D_65", "1.0.0", "N/A", vector::autocorr2d_64),
    ("AUTOCORR2D_66", "1.0.0", "N/A", vector::autocorr2d_65),
    ("AUTOCORR2D_67", "1.0.0", "N/A", vector::autocorr2d_66),
    ("AUTOCORR2D_68", "1.0.0", "N/A", vector::autocorr2d_67),
    ("AUTOCORR2D_69", "1.0.0", "N/A", vector::autocorr2d_68),
    ("AUTOCORR2D_70", "1.0.0", "N/A", vector::autocorr2d_69),
    ("AUTOCORR2D_71", "1.0.0", "N/A", vector::autocorr2d_70),
    ("AUTOCORR2D_72", "1.0.0", "N/A", vector::autocorr2d_71),
    ("AUTOCORR2D_73", "1.0.0", "N/A", vector::autocorr2d_72),
    ("AUTOCORR2D_74", "1.0.0", "N/A", vector::autocorr2d_73),
    ("AUTOCORR2D_75", "1.0.0", "N/A", vector::autocorr2d_74),
    ("AUTOCORR2D_76", "1.0.0", "N/A", vector::autocorr2d_75),
    ("AUTOCORR2D_77", "1.0.0", "N/A", vector::autocorr2d_76),
    ("AUTOCORR2D_78", "1.0.0", "N/A", vector::autocorr2d_77),
    ("AUTOCORR2D_79", "1.0.0", "N/A", vector::autocorr2d_78),
    ("AUTOCORR2D_80", "1.0.0", "N/A", vector::autocorr2d_79),
    ("AUTOCORR2D_81", "1.0.0", "N/A", vector::autocorr2d_80),
    ("AUTOCORR2D_82", "1.0.0", "N/A", vector::autocorr2d_81),
    ("AUTOCORR2D_83", "1.0.0", "N/A", vector::autocorr2d_82),
    ("AUTOCORR2D_84", "1.0.0", "N/A", vector::autocorr2d_83),
    ("AUTOCORR2D_85", "1.0.0", "N/A", vector::autocorr2d_84),
    ("AUTOCORR2D_86", "1.0.0", "N/A", vector::autocorr2d_85),
    ("AUTOCORR2D_87", "1.0.0", "N/A", vector::autocorr2d_86),
    ("AUTOCORR2D_88", "1.0.0", "N/A", vector::autocorr2d_87),
    ("AUTOCORR2D_89", "1.0.0", "N/A", vector::autocorr2d_88),
    ("AUTOCORR2D_90", "1.0.0", "N/A", vector::autocorr2d_89),
    ("AUTOCORR2D_91", "1.0.0", "N/A", vector::autocorr2d_90),
    ("AUTOCORR2D_92", "1.0.0", "N/A", vector::autocorr2d_91),
    ("AUTOCORR2D_93", "1.0.0", "N/A", vector::autocorr2d_92),
    ("AUTOCORR2D_94", "1.0.0", "N/A", vector::autocorr2d_93),
    ("AUTOCORR2D_95", "1.0.0", "N/A", vector::autocorr2d_94),
    ("AUTOCORR2D_96", "1.0.0", "N/A", vector::autocorr2d_95),
    ("AUTOCORR2D_97", "1.0.0", "N/A", vector::autocorr2d_96),
    ("AUTOCORR2D_98", "1.0.0", "N/A", vector::autocorr2d_97),
    ("AUTOCORR2D_99", "1.0.0", "N/A", vector::autocorr2d_98),
    ("AUTOCORR2D_100", "1.0.0", "N/A", vector::autocorr2d_99),
    ("AUTOCORR2D_101", "1.0.0", "N/A", vector::autocorr2d_100),
    ("AUTOCORR2D_102", "1.0.0", "N/A", vector::autocorr2d_101),
    ("AUTOCORR2D_103", "1.0.0", "N/A", vector::autocorr2d_102),
    ("AUTOCORR2D_104", "1.0.0", "N/A", vector::autocorr2d_103),
    ("AUTOCORR2D_105", "1.0.0", "N/A", vector::autocorr2d_104),
    ("AUTOCORR2D_106", "1.0.0", "N/A", vector::autocorr2d_105),
    ("AUTOCORR2D_107", "1.0.0", "N/A", vector::autocorr2d_106),
    ("AUTOCORR2D_108", "1.0.0", "N/A", vector::autocorr2d_107),
    ("AUTOCORR2D_109", "1.0.0", "N/A", vector::autocorr2d_108),
    ("AUTOCORR2D_110", "1.0.0", "N/A", vector::autocorr2d_109),
    ("AUTOCORR2D_111", "1.0.0", "N/A", vector::autocorr2d_110),
    ("AUTOCORR2D_112", "1.0.0", "N/A", vector::autocorr2d_111),
    ("AUTOCORR2D_113", "1.0.0", "N/A", vector::autocorr2d_112),
    ("AUTOCORR2D_114", "1.0.0", "N/A", vector::autocorr2d_113),
    ("AUTOCORR2D_115", "1.0.0", "N/A", vector::autocorr2d_114),
    ("AUTOCORR2D_116", "1.0.0", "N/A", vector::autocorr2d_115),
    ("AUTOCORR2D_117", "1.0.0", "N/A", vector::autocorr2d_116),
    ("AUTOCORR2D_118", "1.0.0", "N/A", vector::autocorr2d_117),
    ("AUTOCORR2D_119", "1.0.0", "N/A", vector::autocorr2d_118),
    ("AUTOCORR2D_120", "1.0.0", "N/A", vector::autocorr2d_119),
    ("AUTOCORR2D_121", "1.0.0", "N/A", vector::autocorr2d_120),
    ("AUTOCORR2D_122", "1.0.0", "N/A", vector::autocorr2d_121),
    ("AUTOCORR2D_123", "1.0.0", "N/A", vector::autocorr2d_122),
    ("AUTOCORR2D_124", "1.0.0", "N/A", vector::autocorr2d_123),
    ("AUTOCORR2D_125", "1.0.0", "N/A", vector::autocorr2d_124),
    ("AUTOCORR2D_126", "1.0.0", "N/A", vector::autocorr2d_125),
    ("AUTOCORR2D_127", "1.0.0", "N/A", vector::autocorr2d_126),
    ("AUTOCORR2D_128", "1.0.0", "N/A", vector::autocorr2d_127),
    ("AUTOCORR2D_129", "1.0.0", "N/A", vector::autocorr2d_128),
    ("AUTOCORR2D_130", "1.0.0", "N/A", vector::autocorr2d_129),
    ("AUTOCORR2D_131", "1.0.0", "N/A", vector::autocorr2d_130),
    ("AUTOCORR2D_132", "1.0.0", "N/A", vector::autocorr2d_131),
    ("AUTOCORR2D_133", "1.0.0", "N/A", vector::autocorr2d_132),
    ("AUTOCORR2D_134", "1.0.0", "N/A", vector::autocorr2d_133),
    ("AUTOCORR2D_135", "1.0.0", "N/A", vector::autocorr2d_134),
    ("AUTOCORR2D_136", "1.0.0", "N/A", vector::autocorr2d_135),
    ("AUTOCORR2D_137", "1.0.0", "N/A", vector::autocorr2d_136),
    ("AUTOCORR2D_138", "1.0.0", "N/A", vector::autocorr2d_137),
    ("AUTOCORR2D_139", "1.0.0", "N/A", vector::autocorr2d_138),
    ("AUTOCORR2D_140", "1.0.0", "N/A", vector::autocorr2d_139),
    ("AUTOCORR2D_141", "1.0.0", "N/A", vector::autocorr2d_140),
    ("AUTOCORR2D_142", "1.0.0", "N/A", vector::autocorr2d_141),
    ("AUTOCORR2D_143", "1.0.0", "N/A", vector::autocorr2d_142),
    ("AUTOCORR2D_144", "1.0.0", "N/A", vector::autocorr2d_143),
    ("AUTOCORR2D_145", "1.0.0", "N/A", vector::autocorr2d_144),
    ("AUTOCORR2D_146", "1.0.0", "N/A", vector::autocorr2d_145),
    ("AUTOCORR2D_147", "1.0.0", "N/A", vector::autocorr2d_146),
    ("AUTOCORR2D_148", "1.0.0", "N/A", vector::autocorr2d_147),
    ("AUTOCORR2D_149", "1.0.0", "N/A", vector::autocorr2d_148),
    ("AUTOCORR2D_150", "1.0.0", "N/A", vector::autocorr2d_149),
    ("AUTOCORR2D_151", "1.0.0", "N/A", vector::autocorr2d_150),
    ("AUTOCORR2D_152", "1.0.0", "N/A", vector::autocorr2d_151),
    ("AUTOCORR2D_153", "1.0.0", "N/A", vector::autocorr2d_152),
    ("AUTOCORR2D_154", "1.0.0", "N/A", vector::autocorr2d_153),
    ("AUTOCORR2D_155", "1.0.0", "N/A", vector::autocorr2d_154),
    ("AUTOCORR2D_156", "1.0.0", "N/A", vector::autocorr2d_155),
    ("AUTOCORR2D_157", "1.0.0", "N/A", vector::autocorr2d_156),
    ("AUTOCORR2D_158", "1.0.0", "N/A", vector::autocorr2d_157),
    ("AUTOCORR2D_159", "1.0.0", "N/A", vector::autocorr2d_158),
    ("AUTOCORR2D_160", "1.0.0", "N/A", vector::autocorr2d_159),
    ("AUTOCORR2D_161", "1.0.0", "N/A", vector::autocorr2d_160),
    ("AUTOCORR2D_162", "1.0.0", "N/A", vector::autocorr2d_161),
    ("AUTOCORR2D_163", "1.0.0", "N/A", vector::autocorr2d_162),
    ("AUTOCORR2D_164", "1.0.0", "N/A", vector::autocorr2d_163),
    ("AUTOCORR2D_165", "1.0.0", "N/A", vector::autocorr2d_164),
    ("AUTOCORR2D_166", "1.0.0", "N/A", vector::autocorr2d_165),
    ("AUTOCORR2D_167", "1.0.0", "N/A", vector::autocorr2d_166),
    ("AUTOCORR2D_168", "1.0.0", "N/A", vector::autocorr2d_167),
    ("AUTOCORR2D_169", "1.0.0", "N/A", vector::autocorr2d_168),
    ("AUTOCORR2D_170", "1.0.0", "N/A", vector::autocorr2d_169),
    ("AUTOCORR2D_171", "1.0.0", "N/A", vector::autocorr2d_170),
    ("AUTOCORR2D_172", "1.0.0", "N/A", vector::autocorr2d_171),
    ("AUTOCORR2D_173", "1.0.0", "N/A", vector::autocorr2d_172),
    ("AUTOCORR2D_174", "1.0.0", "N/A", vector::autocorr2d_173),
    ("AUTOCORR2D_175", "1.0.0", "N/A", vector::autocorr2d_174),
    ("AUTOCORR2D_176", "1.0.0", "N/A", vector::autocorr2d_175),
    ("AUTOCORR2D_177", "1.0.0", "N/A", vector::autocorr2d_176),
    ("AUTOCORR2D_178", "1.0.0", "N/A", vector::autocorr2d_177),
    ("AUTOCORR2D_179", "1.0.0", "N/A", vector::autocorr2d_178),
    ("AUTOCORR2D_180", "1.0.0", "N/A", vector::autocorr2d_179),
    ("AUTOCORR2D_181", "1.0.0", "N/A", vector::autocorr2d_180),
    ("AUTOCORR2D_182", "1.0.0", "N/A", vector::autocorr2d_181),
    ("AUTOCORR2D_183", "1.0.0", "N/A", vector::autocorr2d_182),
    ("AUTOCORR2D_184", "1.0.0", "N/A", vector::autocorr2d_183),
    ("AUTOCORR2D_185", "1.0.0", "N/A", vector::autocorr2d_184),
    ("AUTOCORR2D_186", "1.0.0", "N/A", vector::autocorr2d_185),
    ("AUTOCORR2D_187", "1.0.0", "N/A", vector::autocorr2d_186),
    ("AUTOCORR2D_188", "1.0.0", "N/A", vector::autocorr2d_187),
    ("AUTOCORR2D_189", "1.0.0", "N/A", vector::autocorr2d_188),
    ("AUTOCORR2D_190", "1.0.0", "N/A", vector::autocorr2d_189),
    ("AUTOCORR2D_191", "1.0.0", "N/A", vector::autocorr2d_190),
    ("AUTOCORR2D_192", "1.0.0", "N/A", vector::autocorr2d_191),
];

/// The number of descriptors in `descList` before
/// `setup_autocorr_descriptors`.
const NUM_DESC_LIST: usize = 217;

/// The index of `AvgIpc`, from which Python's `setupAUTOCorrDescriptors`
/// lists the descriptors of the modules it gathers a second time.
const FIRST_MODULE_DESCRIPTOR: usize = 26;

static AUTOCORR: AtomicBool = AtomicBool::new(false);

/// The descriptor named `name` in Python's `Descriptors` module: one of
/// `descList`, or an `AUTOCORR2D_*` descriptor.
pub(crate) fn lookup(name: &str) -> Option<&'static DescriptorEntry> {
    DESCRIPTORS.iter().find(|entry| entry.0 == name)
}

/// Every descriptor's name and function, in Python's order. After
/// [`setup_autocorr_descriptors`], as in Python, the 192 `AUTOCORR2D_*`
/// descriptors come first, and those from `GraphDescriptors` on appear
/// twice. Mirrors `Descriptors.descList`.
pub fn desc_list() -> Vec<(&'static str, DescriptorFn)> {
    let entries: Vec<&DescriptorEntry> = if AUTOCORR.load(Ordering::Relaxed) {
        DESCRIPTORS[NUM_DESC_LIST..]
            .iter()
            .chain(&DESCRIPTORS[..NUM_DESC_LIST])
            .chain(&DESCRIPTORS[FIRST_MODULE_DESCRIPTOR..NUM_DESC_LIST])
            .collect()
    } else {
        DESCRIPTORS[..NUM_DESC_LIST].iter().collect()
    };
    entries
        .into_iter()
        .map(|&(name, _, _, function)| (name, function))
        .collect()
}

/// Adds the 192 `AUTOCORR2D_*` descriptors to [`desc_list`], and so to
/// [`calc_mol_descriptors`], as Python's does. Mirrors
/// `Descriptors.setupAUTOCorrDescriptors`.
///
/// **This changes process-wide state.** Every later call of `desc_list` or
/// `calc_mol_descriptors`, on any thread and from any crate in the process,
/// sees the longer list, and nothing undoes it. A library should not call
/// it on its callers' behalf; to compute the `AUTOCORR2D_*` descriptors
/// alone, call [`calc_autocorr_2d`](crate::chem::rd_mol_descriptors::calc_autocorr_2d)
/// instead.
pub fn setup_autocorr_descriptors() {
    AUTOCORR.store(true, Ordering::Relaxed);
}

/// Every descriptor of [`desc_list`], by name, with `None` for one that
/// fails, in the order of Python's dictionary: `desc_list`'s, each name
/// once, where it first appears. Mirrors `Descriptors.CalcMolDescriptors`.
pub fn calc_mol_descriptors(mol: &Mol) -> Vec<(&'static str, Option<f64>)> {
    let mut seen = HashSet::new();
    desc_list()
        .into_iter()
        .filter(|(name, _)| seen.insert(*name))
        .map(|(name, function)| (name, function(mol).ok()))
        .collect()
}

/// With `missing_val` for a descriptor that fails. Python's `silent`, which
/// prints the failures' tracebacks when false, is not bound. Mirrors
/// `Descriptors.CalcMolDescriptors(mol, missingVal)`.
pub fn calc_mol_descriptors_with_missing_val(
    mol: &Mol,
    missing_val: f64,
) -> Vec<(&'static str, f64)> {
    calc_mol_descriptors(mol)
        .into_iter()
        .map(|(name, value)| (name, value.unwrap_or(missing_val)))
        .collect()
}
