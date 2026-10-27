// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's 3D descriptors, of a conformer: principal moments and
//! the shape indices from them, the plane of best fit, autocorrelations,
//! RDF, MoRSE, WHIM, Coulomb matrices, and EEM charges.
//!
//! GETAWAY (`rdMolDescriptors.CalcGETAWAY`) is not bound: RDKit 2024.09's
//! `GetGeodesicMatrix` (GETAWAY.cpp) reserves a vector and writes past its
//! size on every call (fixed in RDKit 2025.03.4, rdkit/rdkit#8601), so no
//! input check could make it sound (U1).
//!
//! Each takes a conformer id, −1 for the default conformer. A molecule
//! without conformers is `Error::Runtime`, as in Python, and a missing id
//! `Error::Value`.
//!
//! The moments and the indices from them, the plane of best fit, WHIM, and
//! EEM solve matrices of the coordinates, and Eigen leaves its
//! output uninitialized when they are not finite, so for them coordinates
//! beyond 1e9 or not finite are `Error::Value`, as are WHIM's custom weights.
//! EEM divides by the atoms' distances, so two atoms within 1e-6 of each
//! other are `Error::Value` there too, as embedding `[Na+].[Cl-]`
//! leaves them.

use crate::Result;
use crate::chem::Mol;
use crate::ffi::descriptors_3d::ffi;

/// Defines a moment descriptor: a parameters struct with Python's defaults
/// (A2) and the two functions.
macro_rules! moment_descriptor {
    ($params:ident, $name:ident, $with_params:ident, $python:literal, $what:literal) => {
        #[doc = concat!("Parameters for `", stringify!($with_params), "`, with Python's defaults (A2).")]
        #[derive(Clone, Debug, PartialEq)]
        pub struct $params {
            /// The conformer; −1 is the default conformer.
            pub conf_id: i32,
            /// Weights the atoms by mass; otherwise equally.
            pub use_atomic_masses: bool,
        }

        impl Default for $params {
            fn default() -> $params {
                $params {
                    conf_id: -1,
                    use_atomic_masses: true,
                }
            }
        }

        #[doc = concat!(
                    $what, ", mass-weighted, for the default conformer. Mirrors `rdMolDescriptors.",
                    $python, "`; `force` is not bound (CONVENTIONS.md U8)."
                )]
        pub fn $name(mol: &Mol) -> Result<f64> {
            $with_params(mol, &$params::default())
        }

        #[doc = concat!(
                    "Mirrors `rdMolDescriptors.", $python, "(mol, confId, useAtomicMasses)`."
                )]
        pub fn $with_params(mol: &Mol, params: &$params) -> Result<f64> {
            Ok(ffi::$name(mol.as_ffi(), params.conf_id, params.use_atomic_masses)?)
        }
    };
}

moment_descriptor!(
    CalcPMI1Params,
    calc_pmi1,
    calc_pmi1_with_params,
    "CalcPMI1",
    "The smallest principal moment of inertia"
);
moment_descriptor!(
    CalcPMI2Params,
    calc_pmi2,
    calc_pmi2_with_params,
    "CalcPMI2",
    "The middle principal moment of inertia"
);
moment_descriptor!(
    CalcPMI3Params,
    calc_pmi3,
    calc_pmi3_with_params,
    "CalcPMI3",
    "The largest principal moment of inertia"
);
moment_descriptor!(
    CalcNPR1Params,
    calc_npr1,
    calc_npr1_with_params,
    "CalcNPR1",
    "The normalized principal moments ratio PMI1 / PMI3"
);
moment_descriptor!(
    CalcNPR2Params,
    calc_npr2,
    calc_npr2_with_params,
    "CalcNPR2",
    "The normalized principal moments ratio PMI2 / PMI3"
);
moment_descriptor!(
    CalcRadiusOfGyrationParams,
    calc_radius_of_gyration,
    calc_radius_of_gyration_with_params,
    "CalcRadiusOfGyration",
    "The radius of gyration"
);
moment_descriptor!(
    CalcInertialShapeFactorParams,
    calc_inertial_shape_factor,
    calc_inertial_shape_factor_with_params,
    "CalcInertialShapeFactor",
    "The inertial shape factor, PMI2 / (PMI1 PMI3)"
);
moment_descriptor!(
    CalcEccentricityParams,
    calc_eccentricity,
    calc_eccentricity_with_params,
    "CalcEccentricity",
    "The eccentricity"
);
moment_descriptor!(
    CalcAsphericityParams,
    calc_asphericity,
    calc_asphericity_with_params,
    "CalcAsphericity",
    "The asphericity"
);

/// The spherocity index, from the unweighted gyration moments, for the
/// default conformer. Mirrors `rdMolDescriptors.CalcSpherocityIndex`;
/// `force` is not bound (U8).
pub fn calc_spherocity_index(mol: &Mol) -> Result<f64> {
    calc_spherocity_index_with_conf_id(mol, -1)
}

/// Mirrors `rdMolDescriptors.CalcSpherocityIndex(mol, confId)`.
pub fn calc_spherocity_index_with_conf_id(mol: &Mol, conf_id: i32) -> Result<f64> {
    Ok(ffi::calc_spherocity_index(mol.as_ffi(), conf_id)?)
}

/// The plane of best fit: the atoms' mean distance from it, 0 below four
/// atoms or for a 2D conformer. Mirrors `rdMolDescriptors.CalcPBF`, which
/// in Python reads the plane cached for whichever conformer came first; here
/// it is each conformer's own (U8).
pub fn calc_pbf(mol: &Mol) -> Result<f64> {
    calc_pbf_with_conf_id(mol, -1)
}

/// Mirrors `rdMolDescriptors.CalcPBF(mol, confId)`.
pub fn calc_pbf_with_conf_id(mol: &Mol, conf_id: i32) -> Result<f64> {
    Ok(ffi::calc_pbf(mol.as_ffi(), conf_id)?)
}

/// Defines the parameters of a descriptor that takes a conformer and a
/// custom atom property.
macro_rules! custom_property_params {
    ($(#[$meta:meta])* $params:ident { $($(#[$field_meta:meta])* $field:ident: $ty:ty = $default:expr,)* }) => {
        $(#[$meta])*
        #[derive(Clone, Debug, PartialEq)]
        pub struct $params {
            /// The conformer; −1 is the default conformer.
            pub conf_id: i32,
            $($(#[$field_meta])* pub $field: $ty,)*
            /// Each atom's `f64` property to use (1 where an atom lacks it),
            /// which takes any atoms; empty means RDKit's property tables,
            /// which take atomic numbers 1 to 110 (`Error::Value` for
            /// others, which RDKit would read past the tables' ends).
            pub custom_atom_property: String,
        }

        impl Default for $params {
            fn default() -> $params {
                $params {
                    conf_id: -1,
                    $($field: $default,)*
                    custom_atom_property: String::new(),
                }
            }
        }
    };
}

custom_property_params! {
    /// Parameters for `calc_autocorr_3d_with_params`, with Python's defaults
    /// (A2).
    CalcAUTOCORR3DParams {}
}

/// The 80 3D autocorrelations of eight atom properties over topological
/// distances 1 to 10. `Error::Value` for a molecule without atoms, where
/// RDKit takes the address of a first element that does not exist, and
/// above 46340 atoms (U3). Mirrors `rdMolDescriptors.CalcAUTOCORR3D`.
pub fn calc_autocorr_3d(mol: &Mol) -> Result<Vec<f64>> {
    calc_autocorr_3d_with_params(mol, &CalcAUTOCORR3DParams::default())
}

/// With a custom property, 10 values. Mirrors
/// `rdMolDescriptors.CalcAUTOCORR3D(mol, confId, CustomAtomProperty)`.
pub fn calc_autocorr_3d_with_params(mol: &Mol, params: &CalcAUTOCORR3DParams) -> Result<Vec<f64>> {
    Ok(ffi::calc_autocorr_3d(
        mol.as_ffi(),
        params.conf_id,
        &params.custom_atom_property,
    )?)
}

custom_property_params! {
    /// Parameters for `calc_rdf_with_params`, with Python's defaults (A2).
    CalcRDFParams {}
}

/// The 210 radial distribution function descriptors: seven atom properties
/// at 30 radii. `Error::Value` above 46340 atoms (U3). Mirrors
/// `rdMolDescriptors.CalcRDF`.
pub fn calc_rdf(mol: &Mol) -> Result<Vec<f64>> {
    calc_rdf_with_params(mol, &CalcRDFParams::default())
}

/// With a custom property, 30 values. Mirrors
/// `rdMolDescriptors.CalcRDF(mol, confId, CustomAtomProperty)`.
pub fn calc_rdf_with_params(mol: &Mol, params: &CalcRDFParams) -> Result<Vec<f64>> {
    Ok(ffi::calc_rdf(
        mol.as_ffi(),
        params.conf_id,
        &params.custom_atom_property,
    )?)
}

custom_property_params! {
    /// Parameters for `calc_morse_with_params`, with Python's defaults (A2).
    CalcMORSEParams {}
}

/// The 224 MoRSE descriptors: seven atom properties at 32 scattering
/// parameters. `Error::Value` above 46340 atoms (U3). Mirrors
/// `rdMolDescriptors.CalcMORSE`.
pub fn calc_morse(mol: &Mol) -> Result<Vec<f64>> {
    calc_morse_with_params(mol, &CalcMORSEParams::default())
}

/// With a custom property, 32 values. Mirrors
/// `rdMolDescriptors.CalcMORSE(mol, confId, CustomAtomProperty)`.
pub fn calc_morse_with_params(mol: &Mol, params: &CalcMORSEParams) -> Result<Vec<f64>> {
    Ok(ffi::calc_morse(
        mol.as_ffi(),
        params.conf_id,
        &params.custom_atom_property,
    )?)
}

custom_property_params! {
    /// Parameters for `calc_whim_with_params`, with Python's defaults (A2).
    CalcWHIMParams {
        /// The tolerance for symmetric atoms.
        thresh: f64 = 0.001,
    }
}

/// The 114 WHIM descriptors. `Error::Value` for a molecule without atoms,
/// where RDKit takes the address of a first element that does not exist.
/// Mirrors `rdMolDescriptors.CalcWHIM`.
pub fn calc_whim(mol: &Mol) -> Result<Vec<f64>> {
    calc_whim_with_params(mol, &CalcWHIMParams::default())
}

/// With a custom property, 17 values. Mirrors
/// `rdMolDescriptors.CalcWHIM(mol, confId, thresh, CustomAtomProperty)`.
pub fn calc_whim_with_params(mol: &Mol, params: &CalcWHIMParams) -> Result<Vec<f64>> {
    Ok(ffi::calc_whim(
        mol.as_ffi(),
        params.conf_id,
        params.thresh,
        &params.custom_atom_property,
    )?)
}

/// The Coulomb matrix, as rows: Z_i Z_j / r_ij off the diagonal and
/// Z_i^2.4 / 2 on it. Mirrors `rdMolDescriptors.CalcCoulombMat`.
pub fn calc_coulomb_mat(mol: &Mol) -> Result<Vec<Vec<f64>>> {
    calc_coulomb_mat_with_conf_id(mol, -1)
}

/// Mirrors `rdMolDescriptors.CalcCoulombMat(mol, confId)`.
pub fn calc_coulomb_mat_with_conf_id(mol: &Mol, conf_id: i32) -> Result<Vec<Vec<f64>>> {
    let values = ffi::calc_coulomb_mat(mol.as_ffi(), conf_id)?;
    let n = mol.get_num_atoms() as usize;
    Ok(if n == 0 {
        Vec::new()
    } else {
        values.chunks(n).map(<[f64]>::to_vec).collect()
    })
}

/// Each atom's partial charge by the electronegativity equalization method.
/// RDKit's parameter tables cover atomic numbers below 36 and bond orders up
/// to 3 (aromatic bonds are kekulized first); it reads past the tables, or
/// uninitialized memory, for others, which are `Error::Value`, as are
/// molecules above 46339 atoms, whose matrix size overflows `int`. The
/// molecule is unchanged, as in Python. Mirrors
/// `rdMolDescriptors.CalcEEMcharges`.
pub fn calc_eem_charges(mol: &Mol) -> Result<Vec<f64>> {
    calc_eem_charges_with_conf_id(mol, -1)
}

/// Mirrors `rdMolDescriptors.CalcEEMcharges(mol, confId)`.
pub fn calc_eem_charges_with_conf_id(mol: &Mol, conf_id: i32) -> Result<Vec<f64>> {
    Ok(ffi::calc_eem_charges(mol.as_ffi(), conf_id)?)
}
