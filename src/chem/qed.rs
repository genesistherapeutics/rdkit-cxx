// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.QED`, a Python module: the quantitative estimate
//! of drug-likeness (Bickerton, Paolini, Besnard, Muresan, and Hopkins,
//! Nature Chemistry 4, 90-98 (2012)), as RDKit's Biscu-it implementation
//! computes it from eight properties.
//!
//! Python builds its patterns once, on import; each pattern function here
//! parses new molecules, since a `Mol` cannot be shared between threads
//! (CONVENTIONS.md S1).

use std::collections::BTreeMap;

use crate::Result;
use crate::chem::py_math::py_sum;
use crate::chem::rd_mol_descriptors::{
    self as descriptors, NumRotatableBondsOptions, calc_num_hbd,
    calc_num_rotatable_bonds_with_strict,
};
use crate::chem::{
    Mol, crippen, delete_substructs, get_sssr, mol_from_smarts, mol_surf, remove_hs,
};

/// The eight properties QED weighs, or their weights. Python's namedtuple
/// `QEDproperties`; the derived `Debug` and `Clone` stand in for its
/// `__repr__` and `__getnewargs__`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QEDproperties {
    /// The average molecular weight. Python's `MW`.
    pub mw: f64,
    /// Wildman and Crippen's logP. Python's `ALOGP`.
    pub alogp: f64,
    /// The hydrogen-bond acceptors, by [`acceptors`]. Python's `HBA`.
    pub hba: f64,
    /// The hydrogen-bond donors. Python's `HBD`.
    pub hbd: f64,
    /// The topological polar surface area. Python's `PSA`.
    pub psa: f64,
    /// The rotatable bonds, by the strict definition. Python's `ROTB`.
    pub rotb: f64,
    /// The aromatic rings: the SSSR once aliphatic ring atoms are deleted.
    /// Python's `AROM`.
    pub arom: f64,
    /// The structural alerts that match. Python's `ALERTS`.
    pub alerts: f64,
}

impl QEDproperties {
    /// The values in Python's field order.
    fn values(&self) -> [f64; 8] {
        [
            self.mw,
            self.alogp,
            self.hba,
            self.hbd,
            self.psa,
            self.rotb,
            self.arom,
            self.alerts,
        ]
    }
}

/// The parameters of an asymmetric double sigmoidal (ADS) function.
/// Python's namedtuple `ADSparameter`; the derived `Debug` and `Clone`
/// stand in for its `__repr__` and `__getnewargs__`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ADSparameter {
    /// Python's `A`.
    pub a: f64,
    /// Python's `B`.
    pub b: f64,
    /// Python's `C`.
    pub c: f64,
    /// Python's `D`.
    pub d: f64,
    /// Python's `E`.
    pub e: f64,
    /// Python's `F`.
    pub f: f64,
    /// Python's `DMAX`.
    pub dmax: f64,
}

const WEIGHT_MAX: QEDproperties = QEDproperties {
    mw: 0.50,
    alogp: 0.25,
    hba: 0.00,
    hbd: 0.50,
    psa: 0.00,
    rotb: 0.50,
    arom: 0.25,
    alerts: 1.00,
};

const WEIGHT_MEAN: QEDproperties = QEDproperties {
    mw: 0.66,
    alogp: 0.46,
    hba: 0.05,
    hbd: 0.61,
    psa: 0.06,
    rotb: 0.65,
    arom: 0.48,
    alerts: 0.95,
};

const WEIGHT_NONE: QEDproperties = QEDproperties {
    mw: 1.00,
    alogp: 1.00,
    hba: 1.00,
    hbd: 1.00,
    psa: 1.00,
    rotb: 1.00,
    arom: 1.00,
    alerts: 1.00,
};

/// The weights that maximize QED's information content. Mirrors
/// `QED.WEIGHT_MAX`.
pub fn weight_max() -> QEDproperties {
    WEIGHT_MAX
}

/// The mean of the most informative weightings in Bickerton et al., QED's
/// default. Mirrors `QED.WEIGHT_MEAN`.
pub fn weight_mean() -> QEDproperties {
    WEIGHT_MEAN
}

/// Equal weights. Mirrors `QED.WEIGHT_NONE`.
pub fn weight_none() -> QEDproperties {
    WEIGHT_NONE
}

/// Aliphatic ring atoms, which `arom` deletes before counting rings.
/// Mirrors `QED.AliphaticRings`.
pub fn aliphatic_rings() -> Mol {
    mol_from_smarts("[$([A;R][!a])]").expect("the QED SMARTS parse")
}

const ACCEPTOR_SMARTS: [&str; 11] = [
    "[oH0;X2]",
    "[OH1;X2;v2]",
    "[OH0;X2;v2]",
    "[OH0;X1;v2]",
    "[O-;X1]",
    "[SH0;X2;v2]",
    "[SH0;X1;v2]",
    "[S-;X1]",
    "[nH0;X2]",
    "[NH0;X1;v3]",
    "[$([N;+0;X3;v3]);!$(N[C,S]=O)]",
];

/// The hydrogen-bond acceptor SMARTS. Mirrors `QED.AcceptorSmarts`.
pub fn acceptor_smarts() -> &'static [&'static str] {
    &ACCEPTOR_SMARTS
}

/// The hydrogen-bond acceptor patterns. Mirrors `QED.Acceptors`.
pub fn acceptors() -> Vec<Mol> {
    patterns(&ACCEPTOR_SMARTS)
}

const STRUCTURAL_ALERT_SMARTS: [&str; 116] = [
    "*1[O,S,N]*1",
    "[S,C](=[O,S])[F,Br,Cl,I]",
    "[CX4][Cl,Br,I]",
    "[#6]S(=O)(=O)O[#6]",
    "[$([CH]),$(CC)]#CC(=O)[#6]",
    "[$([CH]),$(CC)]#CC(=O)O[#6]",
    "n[OH]",
    "[$([CH]),$(CC)]#CS(=O)(=O)[#6]",
    "C=C(C=O)C=O",
    "n1c([F,Cl,Br,I])cccc1",
    "[CH1](=O)",
    "[#8][#8]",
    "[C;!R]=[N;!R]",
    "[N!R]=[N!R]",
    "[#6](=O)[#6](=O)",
    "[#16][#16]",
    "[#7][NH2]",
    "C(=O)N[NH2]",
    "[#6]=S",
    "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
    "C1(=[O,N])C=CC(=[O,N])C=C1",
    "C1(=[O,N])C(=[O,N])C=CC=C1",
    "a21aa3a(aa1aaaa2)aaaa3",
    "a31a(a2a(aa1)aaaa2)aaaa3",
    "a1aa2a3a(a1)A=AA=A3=AA=A2",
    "c1cc([NH2])ccc1",
    "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
    "I",
    "OS(=O)(=O)[O-]",
    "[N+](=O)[O-]",
    "C(=O)N[OH]",
    "C1NC(=O)NC(=O)1",
    "[SH]",
    "[S-]",
    "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
    "[CR1]1[CR1][CR1]cc[CR1][CR1]1",
    "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
    "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "C#C",
    "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]",
    "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
    "[#6]=N[OH]",
    "[#6]=NOC=O",
    "[#6](=O)[CX4,CR0X3,O][#6](=O)",
    "c1ccc2c(c1)ccc(=O)o2",
    "[O+,o+,S+,s+]",
    "N=C=O",
    "[NX3,NX4][F,Cl,Br,I]",
    "c1ccccc1OC(=O)[#6]",
    "[CR0]=[CR0][CR0]=[CR0]",
    "[C+,c+,C-,c-]",
    "N=[N+]=[N-]",
    "C12C(NC(N1)=O)CSC2",
    "c1c([OH])c([OH,NH2,NH])ccc1",
    "P",
    "[N,O,S]C#N",
    "C=C=O",
    "[Si][F,Cl,Br,I]",
    "[SX2]O",
    "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
    "O1CCCCC1OC2CCC3CCCCC3C2",
    "N=[CR0][N,n,O,S]",
    "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
    "C=[C!r]C#N",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])",
    "[OH]c1ccc([OH,NH2,NH])cc1",
    "c1ccccc1OC(=O)O",
    "[SX2H0][N]",
    "c12ccccc1(SC(S)=N2)",
    "c12ccccc1(SC(=S)N2)",
    "c1nnnn1C=O",
    "s1c(S)nnc1NC=O",
    "S1C=CSC1=S",
    "C(=O)Onnn",
    "OS(=O)(=O)C(F)(F)F",
    "N#CC[OH]",
    "N#CC(=O)",
    "S(=O)(=O)C#N",
    "N[CH2]C#N",
    "C1(=O)NCC1",
    "S(=O)(=O)[O-,OH]",
    "NC[F,Cl,Br,I]",
    "C=[C!r]O",
    "[NX2+0]=[O+0]",
    "[OR0,NR0][OR0,NR0]",
    "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]",
    "[CX2R0][NX3R0]",
    "c1ccccc1[C;!R]=[C;!R]c2ccccc2",
    "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
    "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
    "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]",
    "[*]=[N+]=[*]",
    "[SX3](=O)[O-,OH]",
    "N#N",
    "F.F.F.F",
    "[R0;D2][R0;D2][R0;D2][R0;D2]",
    "[cR,CR]~C(=O)NC(=O)~[cR,CR]",
    "C=!@CC=[O,S]",
    "[#6,#8,#16][#6](=O)O[#6]",
    "c[C;R0](=[O,S])[#6]",
    "c[SX2][C;!R]",
    "C=C=C",
    "c1nc([F,Cl,Br,I,S])ncc1",
    "c1ncnc([F,Cl,Br,I,S])c1",
    "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
    "[#6]S(=O)(=O)c1ccc(cc1)F",
    "[15N]",
    "[13C]",
    "[18O]",
    "[34S]",
];

/// The SMARTS of the unwanted functionalities. Mirrors
/// `QED.StructuralAlertSmarts`.
pub fn structural_alert_smarts() -> &'static [&'static str] {
    &STRUCTURAL_ALERT_SMARTS
}

/// The unwanted functionalities' patterns. Mirrors `QED.StructuralAlerts`.
pub fn structural_alerts() -> Vec<Mol> {
    patterns(&STRUCTURAL_ALERT_SMARTS)
}

fn patterns(smarts: &[&str]) -> Vec<Mol> {
    smarts
        .iter()
        .map(|sma| mol_from_smarts(sma).expect("the QED SMARTS parse"))
        .collect()
}

/// The ADS function's parameters for each property, by Python's name for
/// it. Mirrors `QED.adsParameters`.
pub fn ads_parameters() -> BTreeMap<String, ADSparameter> {
    [
        (
            "MW",
            ADSparameter {
                a: 2.817065973,
                b: 392.5754953,
                c: 290.7489764,
                d: 2.419764353,
                e: 49.22325677,
                f: 65.37051707,
                dmax: 104.9805561,
            },
        ),
        (
            "ALOGP",
            ADSparameter {
                a: 3.172690585,
                b: 137.8624751,
                c: 2.534937431,
                d: 4.581497897,
                e: 0.822739154,
                f: 0.576295591,
                dmax: 131.3186604,
            },
        ),
        (
            "HBA",
            ADSparameter {
                a: 2.948620388,
                b: 160.4605972,
                c: 3.615294657,
                d: 4.435986202,
                e: 0.290141953,
                f: 1.300669958,
                dmax: 148.7763046,
            },
        ),
        (
            "HBD",
            ADSparameter {
                a: 1.618662227,
                b: 1010.051101,
                c: 0.985094388,
                d: 0.000000001,
                e: 0.713820843,
                f: 0.920922555,
                dmax: 258.1632616,
            },
        ),
        (
            "PSA",
            ADSparameter {
                a: 1.876861559,
                b: 125.2232657,
                c: 62.90773554,
                d: 87.83366614,
                e: 12.01999824,
                f: 28.51324732,
                dmax: 104.5686167,
            },
        ),
        (
            "ROTB",
            ADSparameter {
                a: 0.010000000,
                b: 272.4121427,
                c: 2.558379970,
                d: 1.565547684,
                e: 1.271567166,
                f: 2.758063707,
                dmax: 105.4420403,
            },
        ),
        (
            "AROM",
            ADSparameter {
                a: 3.217788970,
                b: 957.7374108,
                c: 2.274627939,
                d: 0.000000001,
                e: 1.317690384,
                f: 0.375760881,
                dmax: 312.3372610,
            },
        ),
        (
            "ALERTS",
            ADSparameter {
                a: 0.010000000,
                b: 1199.094025,
                c: -0.09002883,
                d: 0.000000001,
                e: 0.185904477,
                f: 0.875193782,
                dmax: 417.7253140,
            },
        ),
    ]
    .into_iter()
    .map(|(name, p)| (name.to_owned(), p))
    .collect()
}

/// The parameters in [`QEDproperties`]' field order.
fn ads_parameters_in_order() -> [ADSparameter; 8] {
    let params = ads_parameters();
    ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"].map(|name| params[name])
}

/// The asymmetric double sigmoidal function, scaled by its maximum. Mirrors
/// `QED.ads`.
pub fn ads(x: f64, ads_parameter: &ADSparameter) -> f64 {
    let p = ads_parameter;
    let exp1 = 1.0 + (-(x - p.c + p.d / 2.0) / p.e).exp();
    let exp2 = 1.0 + (-(x - p.c - p.d / 2.0) / p.f).exp();
    let dx = p.a + p.b / exp1 * (1.0 - 1.0 / exp2);
    dx / p.dmax
}

/// The properties QED weighs, of the molecule without its explicit
/// hydrogens. Mirrors `QED.properties`.
pub fn properties(mol: &Mol) -> Result<QEDproperties> {
    let mol = remove_hs(mol)?;
    let mut hba = 0;
    for pattern in acceptors() {
        if mol.has_substruct_match(&pattern)? {
            hba += mol.get_substruct_matches(&pattern)?.len();
        }
    }
    let mut alerts = 0;
    for alert in structural_alerts() {
        if mol.has_substruct_match(&alert)? {
            alerts += 1;
        }
    }
    let mut aromatic = delete_substructs(&mol, &aliphatic_rings())?;
    Ok(QEDproperties {
        mw: descriptors::calc_mol_wt(&mol, false)?,
        alogp: crippen::mol_log_p(&mol)?,
        hba: hba as f64,
        hbd: f64::from(calc_num_hbd(&mol)?),
        psa: mol_surf::tpsa(&mol)?,
        rotb: f64::from(calc_num_rotatable_bonds_with_strict(
            &mol,
            NumRotatableBondsOptions::Strict,
        )?),
        arom: get_sssr(&mut aromatic)?.len() as f64,
        alerts: f64::from(alerts),
    })
}

/// Parameters for `qed_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct QedParams {
    /// The properties' weights.
    pub w: QEDproperties,
    /// The properties to weigh; the molecule's when `None`.
    pub qed_properties: Option<QEDproperties>,
}

impl Default for QedParams {
    fn default() -> QedParams {
        QedParams {
            w: WEIGHT_MEAN,
            qed_properties: None,
        }
    }
}

/// QED with the mean weights: the weighted geometric mean of the
/// properties' desirabilities. Mirrors `QED.qed`.
pub fn qed(mol: &Mol) -> Result<f64> {
    qed_with_params(mol, &QedParams::default())
}

/// Mirrors `QED.qed(mol, w, qedProperties)`.
pub fn qed_with_params(mol: &Mol, params: &QedParams) -> Result<f64> {
    let qed_properties = match params.qed_properties {
        Some(qed_properties) => qed_properties,
        None => properties(mol)?,
    };
    let w = params.w.values();
    let d: Vec<f64> = qed_properties
        .values()
        .iter()
        .zip(ads_parameters_in_order())
        .map(|(&pi, p)| ads(pi, &p))
        .collect();
    let t = py_sum(w.iter().zip(&d).map(|(wi, di)| wi * di.ln()));
    Ok((t / py_sum(w)).exp())
}

/// QED with [`weight_max`]. Mirrors `QED.weights_max`.
pub fn weights_max(mol: &Mol) -> Result<f64> {
    qed_with_params(
        mol,
        &QedParams {
            w: WEIGHT_MAX,
            ..QedParams::default()
        },
    )
}

/// QED with [`weight_mean`]. Mirrors `QED.weights_mean`.
pub fn weights_mean(mol: &Mol) -> Result<f64> {
    qed(mol)
}

/// QED with [`weight_none`]. Mirrors `QED.weights_none`.
pub fn weights_none(mol: &Mol) -> Result<f64> {
    qed_with_params(
        mol,
        &QedParams {
            w: WEIGHT_NONE,
            ..QedParams::default()
        },
    )
}

/// QED with the mean weights. Mirrors `QED.default`.
pub fn default(mol: &Mol) -> Result<f64> {
    weights_mean(mol)
}
