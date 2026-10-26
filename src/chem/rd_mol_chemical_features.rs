// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdMolChemicalFeatures`: feature factories built from
//! feature definitions, and the features they find in molecules.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::path::Path;

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::mol_chemical_features::ffi;
use crate::geometry::rd_geometry::Point3D;
use crate::{Error, Result};

// Every string is UTF-8: the definitions came from a &str.
fn text(value: std::result::Result<UniquePtr<cxx::CxxString>, cxx::Exception>) -> String {
    to_string(expect_ok(value)).expect("a feature definition string is not UTF-8")
}

/// Finds chemical features (donors, acceptors, rings, and so on) in
/// molecules, by matching the SMARTS patterns of its feature definitions.
/// Mirrors `rdMolChemicalFeatures.MolChemicalFeatureFactory`.
///
/// It is `Send` and `Sync`: its matches hold the crate's Ctrl-C lock, so
/// threads that share one match one at a time.
pub struct MolChemicalFeatureFactory {
    /// Never null.
    ptr: UniquePtr<ffi::MolChemicalFeatureFactory>,
}

impl fmt::Debug for MolChemicalFeatureFactory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MolChemicalFeatureFactory")
            .field("feature_defs", &self.get_feature_defs())
            .finish()
    }
}

/// A factory from the feature definitions in a file (an `.fdef` file).
/// `Error::Os` if the file cannot be read, `Error::Value` if it is not
/// UTF-8, and otherwise as `build_feature_factory_from_string`. Mirrors
/// `rdMolChemicalFeatures.BuildFeatureFactory`.
pub fn build_feature_factory(file_name: impl AsRef<Path>) -> Result<MolChemicalFeatureFactory> {
    let path = file_name.as_ref();
    let bytes = std::fs::read(path).map_err(|err| Error::Os {
        message: format!("File: {} could not be opened: {err}", path.display()),
    })?;
    let text = String::from_utf8(bytes).map_err(|err| Error::Value {
        message: format!("the feature definitions are not UTF-8: {err}"),
    })?;
    build_feature_factory_from_string(&text)
}

/// A factory from feature definitions: `AtomType` lines, and
/// `DefineFeature` blocks with a SMARTS pattern, a `Family`, optional
/// `Weights`, and `EndFeature`. `Error::Parse` for definitions RDKit cannot
/// parse (Python's `FeatureFileParseException`), and `Error::Value` for
/// weights that do not match a pattern's atoms. Mirrors
/// `rdMolChemicalFeatures.BuildFeatureFactoryFromString`.
pub fn build_feature_factory_from_string(fdef_string: &str) -> Result<MolChemicalFeatureFactory> {
    let ptr = ffi::build_feature_factory_from_string(fdef_string)?;
    assert!(!ptr.is_null(), "RDKit returned a null feature factory");
    Ok(MolChemicalFeatureFactory { ptr })
}

/// Parameters for `MolChemicalFeatureFactory::get_mol_feature_with_params`,
/// with Python's defaults (A2). Python's `recompute`, which only chooses
/// between computing the features and reading a list that Python's binding
/// keeps from the last call (for any molecule, from any thread), is not
/// bound: the features are computed every time.
#[derive(Clone, Debug, PartialEq)]
pub struct GetMolFeatureParams {
    /// Only features of this family; empty for all.
    pub include_only: String,
    /// The conformer the feature's position comes from; -1 for the first.
    pub conf_id: i32,
}

impl Default for GetMolFeatureParams {
    fn default() -> GetMolFeatureParams {
        GetMolFeatureParams {
            include_only: String::new(),
            conf_id: -1,
        }
    }
}

impl MolChemicalFeatureFactory {
    /// The number of feature definitions. Mirrors
    /// `MolChemicalFeatureFactory.GetNumFeatureDefs`.
    pub fn get_num_feature_defs(&self) -> i32 {
        expect_ok(ffi::factory_get_num_feature_defs(&self.ptr))
    }

    fn num_defs(&self) -> u32 {
        u32::try_from(self.get_num_feature_defs())
            .expect("a negative number of feature definitions")
    }

    /// The definitions' families, each once, in the order they first
    /// appear. Mirrors `MolChemicalFeatureFactory.GetFeatureFamilies`.
    pub fn get_feature_families(&self) -> Vec<String> {
        let mut families: Vec<String> = Vec::new();
        for i in 0..self.num_defs() {
            let family = text(ffi::factory_feature_def_family(&self.ptr, i));
            if !families.contains(&family) {
                families.push(family);
            }
        }
        families
    }

    /// Each definition's SMARTS, keyed by `"family.type"`; a later
    /// definition with the same key replaces an earlier one. Mirrors
    /// `MolChemicalFeatureFactory.GetFeatureDefs`.
    pub fn get_feature_defs(&self) -> BTreeMap<String, String> {
        (0..self.num_defs())
            .map(|i| {
                let family = text(ffi::factory_feature_def_family(&self.ptr, i));
                let type_ = text(ffi::factory_feature_def_type(&self.ptr, i));
                let smarts = text(ffi::factory_feature_def_smarts(&self.ptr, i));
                (format!("{family}.{type_}"), smarts)
            })
            .collect()
    }

    /// The number of features the molecule has. Mirrors
    /// `MolChemicalFeatureFactory.GetNumMolFeatures`.
    pub fn get_num_mol_features(&self, mol: &Mol) -> Result<i32> {
        self.get_num_mol_features_with_include_only(mol, "")
    }

    /// The number of features of family `include_only` (all for an empty
    /// one, and RDKit reads it only up to a NUL). Mirrors
    /// `MolChemicalFeatureFactory.GetNumMolFeatures(mol, includeOnly)`.
    pub fn get_num_mol_features_with_include_only(
        &self,
        mol: &Mol,
        include_only: &str,
    ) -> Result<i32> {
        Ok(ffi::factory_get_num_mol_features(
            &self.ptr,
            mol.as_ffi(),
            include_only,
        )?)
    }

    /// The molecule's feature `idx`, in the order of the definitions and
    /// then of their matches; `Error::Index` past the last. The feature
    /// borrows the molecule and the factory. Mirrors
    /// `MolChemicalFeatureFactory.GetMolFeature`.
    pub fn get_mol_feature<'a>(&'a self, mol: &'a Mol, idx: i32) -> Result<MolChemicalFeature<'a>> {
        self.get_mol_feature_with_params(mol, idx, &GetMolFeatureParams::default())
    }

    /// Mirrors `MolChemicalFeatureFactory.GetMolFeature(mol, idx,
    /// includeOnly, recompute, confId)`; see `GetMolFeatureParams`.
    pub fn get_mol_feature_with_params<'a>(
        &'a self,
        mol: &'a Mol,
        idx: i32,
        params: &GetMolFeatureParams,
    ) -> Result<MolChemicalFeature<'a>> {
        let ptr = ffi::factory_get_mol_feature(
            &self.ptr,
            mol.as_ffi(),
            idx,
            &params.include_only,
            params.conf_id,
        )?;
        assert!(!ptr.is_null(), "RDKit returned a null feature");
        Ok(MolChemicalFeature {
            ptr,
            mol,
            factory: self,
        })
    }
}

/// A chemical feature found in a molecule: a family and type from the
/// factory's definition, the matched atoms, and a position from the
/// molecule's coordinates. It borrows the molecule and the factory, as
/// Python's keeps them alive. Mirrors
/// `rdMolChemicalFeatures.MolChemicalFeature`.
pub struct MolChemicalFeature<'a> {
    // Declared first, so it is dropped before what it points into.
    ptr: UniquePtr<ffi::MolChemicalFeature>,
    mol: &'a Mol,
    factory: &'a MolChemicalFeatureFactory,
}

impl fmt::Debug for MolChemicalFeature<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MolChemicalFeature")
            .field("id", &self.get_id())
            .field("family", &self.get_family())
            .field("type", &self.get_type())
            .field("atom_ids", &self.get_atom_ids())
            .finish()
    }
}

impl<'a> MolChemicalFeature<'a> {
    /// The feature's id: its position among the molecule's features,
    /// starting from 1. Mirrors `MolChemicalFeature.GetId`.
    pub fn get_id(&self) -> i32 {
        expect_ok(ffi::mol_chemical_feature_get_id(&self.ptr))
    }

    /// The feature's family. Mirrors `MolChemicalFeature.GetFamily`.
    pub fn get_family(&self) -> String {
        text(ffi::mol_chemical_feature_get_family(&self.ptr))
    }

    /// The feature's specific type. Mirrors `MolChemicalFeature.GetType`.
    pub fn get_type(&self) -> String {
        text(ffi::mol_chemical_feature_get_type(&self.ptr))
    }

    /// The feature's position in the active conformer: its atom's position,
    /// or its atoms' positions averaged with the definition's weights.
    /// `Error::Runtime` for a molecule without conformers or a feature of
    /// several atoms whose definition has no weights, and `Error::Value` for
    /// an active conformer the molecule lacks. Mirrors
    /// `MolChemicalFeature.GetPos()`.
    pub fn get_pos(&self) -> Result<Point3D> {
        Ok(ffi::mol_chemical_feature_get_pos(&self.ptr)?)
    }

    /// The feature's position in conformer `conf_id` (-1 for the first), as
    /// `get_pos`. Mirrors `MolChemicalFeature.GetPos(confId)`.
    pub fn get_pos_with_conf_id(&self, conf_id: i32) -> Result<Point3D> {
        Ok(ffi::mol_chemical_feature_get_pos_with_conf_id(
            &self.ptr, conf_id,
        )?)
    }

    /// The indices of the feature's atoms, in the order of the definition's
    /// pattern. Mirrors `MolChemicalFeature.GetAtomIds`.
    pub fn get_atom_ids(&self) -> Vec<u32> {
        expect_ok(ffi::mol_chemical_feature_get_atom_ids(&self.ptr))
    }

    /// The molecule the feature was found in. Mirrors
    /// `MolChemicalFeature.GetMol`.
    pub fn get_mol(&self) -> &'a Mol {
        self.mol
    }

    /// The factory that found the feature. Mirrors
    /// `MolChemicalFeature.GetFactory`.
    pub fn get_factory(&self) -> &'a MolChemicalFeatureFactory {
        self.factory
    }

    /// The conformer `get_pos` uses; -1 for the first. Mirrors
    /// `MolChemicalFeature.GetActiveConformer`.
    pub fn get_active_conformer(&self) -> i32 {
        expect_ok(ffi::mol_chemical_feature_get_active_conformer(&self.ptr))
    }

    /// Mirrors `MolChemicalFeature.SetActiveConformer`.
    pub fn set_active_conformer(&mut self, conf_id: i32) {
        expect_ok(ffi::mol_chemical_feature_set_active_conformer(
            self.ptr.pin_mut(),
            conf_id,
        ));
    }

    /// Forgets the positions computed so far. Mirrors
    /// `MolChemicalFeature.ClearCache`.
    pub fn clear_cache(&mut self) {
        expect_ok(ffi::mol_chemical_feature_clear_cache(self.ptr.pin_mut()));
    }
}

/// Each feature's atom indices, or an empty list if two of the features share
/// an atom. Mirrors `rdMolChemicalFeatures.GetAtomMatch`.
pub fn get_atom_match(feat_match: &[MolChemicalFeature<'_>]) -> Result<Vec<Vec<u32>>> {
    get_atom_match_with_max_ats(feat_match, 1024)
}

/// As `get_atom_match` for atom indices below `max_ats`: RDKit marks the
/// atoms it has seen in a vector of `maxAts` flags, which it indexes without
/// a check, so an atom at `max_ats` or beyond is `Error::Index`, and a
/// negative `max_ats` `Error::Value`. Mirrors
/// `rdMolChemicalFeatures.GetAtomMatch(featMatch, maxAts)`.
pub fn get_atom_match_with_max_ats(
    feat_match: &[MolChemicalFeature<'_>],
    max_ats: i32,
) -> Result<Vec<Vec<u32>>> {
    let len = u32::try_from(max_ats).map_err(|_| Error::Value {
        message: format!("maxAts must not be negative, not {max_ats}"),
    })?;
    let mut seen = BTreeSet::new();
    let mut result = Vec::with_capacity(feat_match.len());
    for feat in feat_match {
        let ids = feat.get_atom_ids();
        for &idx in &ids {
            if idx >= len {
                return Err(Error::index(idx, len as usize));
            }
            if !seen.insert(idx) {
                return Ok(Vec::new());
            }
        }
        result.push(ids);
    }
    Ok(result)
}
