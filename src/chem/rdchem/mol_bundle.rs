// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::ops::Deref;

use cxx::UniquePtr;

use crate::chem::rdchem::{Mol, SubstructMatchParameters};
use crate::chem::rdmolops::mol_list_into_vec;
use crate::error::expect_ok;
use crate::ffi::geometry::ffi::IndexList;
use crate::ffi::mol_bundle::ffi;
use crate::ffi::substruct::ffi as substruct;
use crate::{Error, Result};

/// A group of related molecules. Mirrors `rdchem.MolBundle`.
///
/// The bundle owns copies of its molecules, and `get_mol` borrows them, where
/// Python's `GetMol` returns a reference that can edit the bundle's molecule.
#[derive(Clone, Debug, Default)]
pub struct MolBundle {
    mols: Vec<Mol>,
}

impl MolBundle {
    /// The bundle a pickle holds, as `to_binary` writes it. Mirrors
    /// `MolBundle(pklString)`.
    ///
    /// # Safety
    ///
    /// `pkl_string` must be a bundle pickle RDKit wrote (this crate's
    /// `to_binary` or RDKit's), unmodified: it holds molecule pickles, which
    /// RDKit reads as `Mol::from_pkl_string` does.
    pub unsafe fn from_pkl_string(pkl_string: &[u8]) -> Result<MolBundle> {
        Ok(MolBundle {
            mols: mol_list_into_vec(ffi::mol_bundle_from_pkl_string(pkl_string)?),
        })
    }

    /// A bundle of molecules RDKit returned (W8).
    pub(crate) fn from_mols(mols: Vec<Mol>) -> MolBundle {
        MolBundle { mols }
    }

    /// Calls `f` with a C++ bundle that borrows the molecules.
    pub(crate) fn with_ffi<T>(&self, f: impl FnOnce(&ffi::MolBundle) -> T) -> T {
        let mut bundle: UniquePtr<ffi::MolBundle> = expect_ok(ffi::mol_bundle_new());
        for mol in &self.mols {
            expect_ok(ffi::mol_bundle_add_borrowed(bundle.pin_mut(), mol.as_ffi()));
        }
        // The bundle drops here, before the borrow of the molecules ends.
        f(&bundle)
    }

    /// A pickle of the molecules, without the bundle's properties. Mirrors
    /// `MolBundle.ToBinary`; `Error::Value` if RDKit was built without
    /// Boost.Serialization (`mol_bundle_can_serialize`), or for a molecule
    /// `Mol::to_binary` refuses.
    pub fn to_binary(&self) -> Result<Vec<u8>> {
        Ok(self.with_ffi(ffi::mol_bundle_to_binary)?)
    }

    /// Adds a copy of `nmol`, returning the new number of molecules. Mirrors
    /// `MolBundle.AddMol`.
    pub fn add_mol(&mut self, nmol: &Mol) -> usize {
        self.mols.push(nmol.clone());
        self.mols.len()
    }

    /// Mirrors `MolBundle.GetMol`.
    pub fn get_mol(&self, idx: usize) -> Result<&Mol> {
        self.mols
            .get(idx)
            .ok_or_else(|| Error::index(idx, self.mols.len()))
    }

    /// The number of molecules. Mirrors `MolBundle.Size`.
    pub fn size(&self) -> usize {
        self.mols.len()
    }
}

/// Defines a substructure query of a bundle by a molecule, and its
/// `_with_params` form (A2, A3).
macro_rules! bundle_matches {
    ($(
        $(#[$meta:meta])*
        $name:ident, $name_with_params:ident, $ffi:ident -> $out:ty $(, $convert:path)?;
    )*) => {
        impl MolBundle {
            $(
                $(#[$meta])*
                pub fn $name(&self, query: &Mol) -> Result<$out> {
                    self.$name_with_params(query, &SubstructMatchParameters::default())
                }

                #[doc = concat!(
                    "As [`", stringify!($name), "`](Self::", stringify!($name),
                    "), with parameters; also covers the keyword form."
                )]
                pub fn $name_with_params(
                    &self,
                    query: &Mol,
                    params: &SubstructMatchParameters,
                ) -> Result<$out> {
                    let params = params.to_ffi();
                    let result =
                        self.with_ffi(|bundle| substruct::$ffi(bundle, query.as_ffi(), &params))?;
                    Ok($($convert)?(result))
                }
            )*
        }
    };
}

bundle_matches! {
    /// Whether `query` matches part of any molecule. Mirrors
    /// `MolBundle.HasSubstructMatch`.
    has_substruct_match, has_substruct_match_with_params,
        mol_bundle_has_substruct_match -> bool;
    /// The first match in the first molecule that matches. Mirrors
    /// `MolBundle.GetSubstructMatch`.
    get_substruct_match, get_substruct_match_with_params,
        mol_bundle_get_substruct_match -> Vec<u32>;
    /// Every match in the first molecule that matches. Mirrors
    /// `MolBundle.GetSubstructMatches`.
    get_substruct_matches, get_substruct_matches_with_params,
        mol_bundle_get_substruct_matches -> Vec<Vec<u32>>, IndexList::into_nested;
}

impl MolBundle {
    /// Whether any molecule in `query` matches part of any molecule here.
    /// Mirrors `MolBundle.HasSubstructMatch(query: MolBundle)` (A5).
    pub fn has_substruct_match_mol_bundle(&self, query: &MolBundle) -> Result<bool> {
        self.has_substruct_match_mol_bundle_with_params(query, &SubstructMatchParameters::default())
    }

    /// Mirrors `MolBundle.HasSubstructMatch(query: MolBundle, params)` and its
    /// keyword form.
    pub fn has_substruct_match_mol_bundle_with_params(
        &self,
        query: &MolBundle,
        params: &SubstructMatchParameters,
    ) -> Result<bool> {
        let params = params.to_ffi();
        Ok(self.with_ffi(|target| {
            query.with_ffi(|query| {
                substruct::mol_bundle_has_substruct_match_mol_bundle(target, query, &params)
            })
        })?)
    }

    /// The first match of a molecule in `query` in a molecule here. Mirrors
    /// `MolBundle.GetSubstructMatch(query: MolBundle)` (A5).
    pub fn get_substruct_match_mol_bundle(&self, query: &MolBundle) -> Result<Vec<u32>> {
        self.get_substruct_match_mol_bundle_with_params(query, &SubstructMatchParameters::default())
    }

    /// Mirrors `MolBundle.GetSubstructMatch(query: MolBundle, params)` and its
    /// keyword form.
    pub fn get_substruct_match_mol_bundle_with_params(
        &self,
        query: &MolBundle,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<u32>> {
        let params = params.to_ffi();
        Ok(self.with_ffi(|target| {
            query.with_ffi(|query| {
                substruct::mol_bundle_get_substruct_match_mol_bundle(target, query, &params)
            })
        })?)
    }

    /// Every match of the first molecule in `query` that matches, in the
    /// first molecule here that it matches. Mirrors
    /// `MolBundle.GetSubstructMatches(query: MolBundle)` (A5).
    pub fn get_substruct_matches_mol_bundle(&self, query: &MolBundle) -> Result<Vec<Vec<u32>>> {
        self.get_substruct_matches_mol_bundle_with_params(
            query,
            &SubstructMatchParameters::default(),
        )
    }

    /// Mirrors `MolBundle.GetSubstructMatches(query: MolBundle, params)` and
    /// its keyword form.
    pub fn get_substruct_matches_mol_bundle_with_params(
        &self,
        query: &MolBundle,
        params: &SubstructMatchParameters,
    ) -> Result<Vec<Vec<u32>>> {
        let params = params.to_ffi();
        Ok(IndexList::into_nested(self.with_ffi(|target| {
            query.with_ffi(|query| {
                substruct::mol_bundle_get_substruct_matches_mol_bundle(target, query, &params)
            })
        })?))
    }
}

/// A bundle whose molecules all have the same numbers of atoms and bonds.
/// Mirrors `rdchem.FixedMolSizeMolBundle`.
///
/// It dereferences to `MolBundle` for everything but `add_mol`, which checks
/// the sizes.
#[derive(Clone, Debug, Default)]
pub struct FixedMolSizeMolBundle {
    bundle: MolBundle,
}

impl FixedMolSizeMolBundle {
    /// Adds a copy of `nmol`, returning the new number of molecules; a
    /// molecule of another size is `Error::Value`. Mirrors
    /// `FixedMolSizeMolBundle.AddMol`.
    pub fn add_mol(&mut self, nmol: &Mol) -> Result<usize> {
        if let Some(first) = self.bundle.mols.first() {
            ffi::fixed_mol_size_mol_bundle_check(first.as_ffi(), nmol.as_ffi())?;
        }
        Ok(self.bundle.add_mol(nmol))
    }
}

impl Deref for FixedMolSizeMolBundle {
    type Target = MolBundle;

    fn deref(&self) -> &MolBundle {
        &self.bundle
    }
}

/// Whether `MolBundle::to_binary` works: RDKit was built with
/// Boost.Serialization. Mirrors `rdchem.MolBundleCanSerialize`.
pub fn mol_bundle_can_serialize() -> bool {
    expect_ok(ffi::mol_bundle_can_serialize())
}
