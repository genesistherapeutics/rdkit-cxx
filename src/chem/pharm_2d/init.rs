// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports the Python code in `rdkit/Chem/Pharm2D/__init__.py`, re-exported
//! from `rdkit_cxx::chem::pharm_2d` as Python's is from
//! `rdkit.Chem.Pharm2D`.

use std::path::PathBuf;

use crate::Result;
use crate::chem::pharm_2d::sig_factory::{SigFactory, SigFactoryParams};
use crate::chem::rd_mol_chemical_features::{
    build_feature_factory, build_feature_factory_from_string,
};

/// RDKit's `Data/BaseFeatures.fdef`, which Python reads from
/// `RDConfig.RDDataDir`.
const BASE_FEATURES: &str =
    include_str!("../../../wrapper/third_party/rdkit/Data/BaseFeatures.fdef");

/// Parameters for [`default_sig_factory_with_params`], with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct DefaultSigFactoryParams {
    /// A feature definition file; `None` for RDKit's `BaseFeatures.fdef`,
    /// which the crate embeds.
    pub fdef_file: Option<PathBuf>,
    /// The fewest points in a pharmacophore.
    pub min_point_count: usize,
    /// The most points in a pharmacophore.
    pub max_point_count: usize,
    /// The distance bins.
    pub bins: Vec<(u32, u32)>,
}

impl Default for DefaultSigFactoryParams {
    fn default() -> DefaultSigFactoryParams {
        DefaultSigFactoryParams {
            fdef_file: None,
            min_point_count: 2,
            max_point_count: 3,
            bins: vec![(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 100)],
        }
    }
}

/// A signature factory over RDKit's base features, less `ZnBinder` and
/// `LumpedHydrophobe`, without triangle pruning, initialized. Mirrors
/// `Pharm2D.DefaultSigFactory`.
pub fn default_sig_factory() -> Result<SigFactory> {
    default_sig_factory_with_params(&DefaultSigFactoryParams::default())
}

/// Mirrors `Pharm2D.DefaultSigFactory(fdefFile, minPointCount,
/// maxPointCount, bins)`.
pub fn default_sig_factory_with_params(params: &DefaultSigFactoryParams) -> Result<SigFactory> {
    let feat_factory = match &params.fdef_file {
        Some(path) => build_feature_factory(path)?,
        None => build_feature_factory_from_string(BASE_FEATURES)?,
    };
    let factory_params = SigFactoryParams {
        skip_feats: vec!["ZnBinder".to_owned(), "LumpedHydrophobe".to_owned()],
        min_point_count: params.min_point_count,
        max_point_count: params.max_point_count,
        triangle_prune_bins: false,
        ..Default::default()
    };
    let mut factory = SigFactory::new_with_params(feat_factory, &factory_params);
    factory.set_bins(&params.bins)?;
    factory.init()?;
    Ok(factory)
}
