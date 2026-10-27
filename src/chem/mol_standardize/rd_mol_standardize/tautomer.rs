// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Enumerating, scoring, and canonicalizing tautomers.

use std::fmt;

use cxx::{CxxString, CxxVector, UniquePtr};

use super::CleanupParameters;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::error::{expect_ok, to_string};
use crate::ffi::standardize_tautomer::{ScoreFunction, TautomerCallback, ffi};
use crate::{Error, Result};

rdkit_enum! {
    /// Whether an enumeration completed, and if not, why. Mirrors
    /// `rdMolStandardize.TautomerEnumeratorStatus`.
    TautomerEnumeratorStatus: ffi::TautomerEnumeratorStatus {
        Completed = Completed,
        MaxTautomersReached = MaxTautomersReached,
        MaxTransformsReached = MaxTransformsReached,
        Canceled = Canceled,
    }
}

fn strings(values: UniquePtr<CxxVector<CxxString>>) -> Result<Vec<String>> {
    values
        .iter()
        .map(|value| {
            value
                .to_str()
                .map(str::to_owned)
                .map_err(|err| Error::Value {
                    message: format!("invalid UTF-8 in a string from RDKit: {err}"),
                })
        })
        .collect()
}

/// The tautomers an enumeration found, by SMILES. Mirrors
/// `rdMolStandardize.TautomerEnumeratorResult`.
///
/// Accessors return copies. A result is `Send`, but not `Sync`, and has no
/// `Clone`: the one a callback is lent shares its tautomers with the
/// enumeration under way.
#[non_exhaustive]
pub struct TautomerEnumeratorResult {
    /// Never null.
    ptr: UniquePtr<ffi::TautomerEnumeratorResult>,
}

impl TautomerEnumeratorResult {
    pub(crate) fn from_ffi(
        ptr: UniquePtr<ffi::TautomerEnumeratorResult>,
    ) -> TautomerEnumeratorResult {
        assert!(
            !ptr.is_null(),
            "RDKit returned a null TautomerEnumeratorResult"
        );
        TautomerEnumeratorResult { ptr }
    }

    fn raw(&self) -> &ffi::TautomerEnumeratorResult {
        self.ptr.as_ref().unwrap()
    }

    /// The number of tautomers. Mirrors `TautomerEnumeratorResult.__len__`.
    pub fn len(&self) -> usize {
        expect_ok(ffi::tautomer_result_size(self.raw()))
    }

    /// Whether there are no tautomers.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The tautomer at a position, in SMILES order; a negative position
    /// counts from the end, as in Python. `Error::Index` past either end.
    /// Mirrors `TautomerEnumeratorResult.__getitem__`.
    pub fn get(&self, pos: i32) -> Result<Mol> {
        let len = self.len();
        let index = if pos < 0 {
            len as i64 + i64::from(pos)
        } else {
            i64::from(pos)
        };
        if index < 0 || index >= len as i64 {
            return Err(Error::index(pos, len));
        }
        Ok(Mol::from_ffi(ffi::tautomer_result_tautomer(
            self.raw(),
            index as usize,
        )?))
    }

    /// The tautomers, in SMILES order. Mirrors
    /// `TautomerEnumeratorResult.tautomers`, which `__call__` and
    /// `__iter__` also give.
    pub fn tautomers(&self) -> Vec<Mol> {
        (0..self.len())
            .map(|i| Mol::from_ffi(expect_ok(ffi::tautomer_result_tautomer(self.raw(), i))))
            .collect()
    }

    /// The tautomers' canonical SMILES, in order. Mirrors
    /// `TautomerEnumeratorResult.smiles`.
    pub fn smiles(&self) -> Result<Vec<String>> {
        strings(expect_ok(ffi::tautomer_result_smiles(self.raw())))
    }

    /// The tautomers by SMILES. Mirrors
    /// `TautomerEnumeratorResult.smilesTautomerMap`.
    pub fn smiles_tautomer_map(&self) -> Result<SmilesTautomerMap> {
        let smiles = self.smiles()?;
        let entries = smiles
            .into_iter()
            .enumerate()
            .map(|(i, smiles)| {
                let tautomer = ResultTautomer {
                    result: self,
                    index: i,
                };
                (smiles, tautomer.owned())
            })
            .collect();
        Ok(SmilesTautomerMap { entries })
    }

    /// Whether the enumeration completed. Mirrors
    /// `TautomerEnumeratorResult.status`.
    pub fn status(&self) -> TautomerEnumeratorStatus {
        TautomerEnumeratorStatus::from_ffi(expect_ok(ffi::tautomer_result_status(self.raw())))
    }

    /// The atoms the transforms changed: the ends of each transform's path.
    /// Mirrors `TautomerEnumeratorResult.modifiedAtoms`.
    pub fn modified_atoms(&self) -> Vec<u32> {
        expect_ok(ffi::tautomer_result_modified_atoms(self.raw()))
    }

    /// The bonds whose order the transforms changed. Mirrors
    /// `TautomerEnumeratorResult.modifiedBonds`.
    pub fn modified_bonds(&self) -> Vec<u32> {
        expect_ok(ffi::tautomer_result_modified_bonds(self.raw()))
    }
}

impl fmt::Debug for TautomerEnumeratorResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TautomerEnumeratorResult")
            .field("smiles", &self.smiles())
            .field("status", &self.status())
            .finish()
    }
}

/// A result's tautomer, before it is copied.
struct ResultTautomer<'a> {
    result: &'a TautomerEnumeratorResult,
    index: usize,
}

impl ResultTautomer<'_> {
    fn owned(&self) -> Tautomer {
        let tautomer = Mol::from_ffi(expect_ok(ffi::tautomer_result_tautomer(
            self.result.raw(),
            self.index,
        )));
        let kekulized = ffi::tautomer_result_kekulized(self.result.raw(), self.index)
            .map(Mol::from_ffi)
            .map_err(Error::from);
        Tautomer {
            tautomer,
            kekulized,
        }
    }
}

/// A tautomer, aromatic and kekulized. Mirrors `rdMolStandardize.Tautomer`.
#[derive(Clone, Debug)]
pub struct Tautomer {
    tautomer: Mol,
    kekulized: Result<Mol>,
}

impl Tautomer {
    /// The aromatic tautomer. Mirrors `Tautomer.tautomer`.
    pub fn tautomer(&self) -> Mol {
        self.tautomer.clone()
    }

    /// The kekulized tautomer: the form RDKit matched transforms against, or
    /// one made now when it had not been needed. Its kekulization's error,
    /// which Python raises here, for one that does not kekulize. Mirrors
    /// `Tautomer.kekulized`.
    pub fn kekulized(&self) -> Result<Mol> {
        self.kekulized.clone()
    }
}

/// A result's tautomers by SMILES, in SMILES order. Mirrors
/// `rdMolStandardize.SmilesTautomerMap`.
#[derive(Clone, Debug)]
pub struct SmilesTautomerMap {
    entries: Vec<(String, Tautomer)>,
}

impl SmilesTautomerMap {
    /// The SMILES. Mirrors `SmilesTautomerMap.keys`.
    pub fn keys(&self) -> Vec<String> {
        self.entries
            .iter()
            .map(|(smiles, _)| smiles.clone())
            .collect()
    }

    /// The tautomers. Mirrors `SmilesTautomerMap.values`.
    pub fn values(&self) -> Vec<&Tautomer> {
        self.entries.iter().map(|(_, tautomer)| tautomer).collect()
    }

    /// The SMILES and tautomers. Mirrors `SmilesTautomerMap.items`.
    pub fn items(&self) -> Vec<(&str, &Tautomer)> {
        self.entries
            .iter()
            .map(|(smiles, tautomer)| (smiles.as_str(), tautomer))
            .collect()
    }

    /// The tautomer with this SMILES. Mirrors `SmilesTautomerMap[smiles]`,
    /// which raises `KeyError` where this returns `None`.
    pub fn get(&self, smiles: &str) -> Option<&Tautomer> {
        self.entries
            .iter()
            .find(|(key, _)| key == smiles)
            .map(|(_, tautomer)| tautomer)
    }

    /// Whether there is a tautomer with this SMILES. Mirrors
    /// `smiles in SmilesTautomerMap`.
    pub fn contains_key(&self, smiles: &str) -> bool {
        self.get(smiles).is_some()
    }

    /// The number of tautomers. Mirrors `SmilesTautomerMap.__len__`.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether there are no tautomers.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Monitors, and can cancel, an enumeration: RDKit calls it before each
/// transform it applies, with the molecule and the result so far, and stops
/// when it returns false. Mirrors
/// `rdMolStandardize.TautomerEnumeratorCallback`, the class Python
/// subclasses to define one. Clones share it.
pub struct TautomerEnumeratorCallback {
    /// Never null.
    ptr: UniquePtr<ffi::TautomerCallbackHandle>,
}

impl fmt::Debug for TautomerEnumeratorCallback {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TautomerEnumeratorCallback")
            .finish_non_exhaustive()
    }
}

impl Clone for TautomerEnumeratorCallback {
    fn clone(&self) -> TautomerEnumeratorCallback {
        TautomerEnumeratorCallback {
            ptr: expect_ok(ffi::tautomer_callback_clone(self.raw())),
        }
    }
}

impl TautomerEnumeratorCallback {
    fn raw(&self) -> &ffi::TautomerCallbackHandle {
        self.ptr.as_ref().unwrap()
    }

    /// A callback that calls `call` with a copy of the molecule and the
    /// result so far. Python subclasses the class and overrides `__call__`;
    /// here the closure is the argument. Enumerators copied from one another
    /// share it, so it must be `Send + Sync`. A panic in it is caught, and
    /// the enumeration that made it returns `Error::Runtime`. Mirrors
    /// `TautomerEnumeratorCallback()`.
    pub fn new<F>(call: F) -> TautomerEnumeratorCallback
    where
        F: Fn(&Mol, &TautomerEnumeratorResult) -> bool + Send + Sync + 'static,
    {
        let callback = Box::new(TautomerCallback(Box::new(call)));
        TautomerEnumeratorCallback {
            ptr: expect_ok(ffi::tautomer_callback_new(callback)),
        }
    }

    /// Calls the closure, as an enumeration does: whether to go on. Mirrors
    /// `TautomerEnumeratorCallback.__call__`.
    pub fn call(&self, mol: &Mol, res: &TautomerEnumeratorResult) -> Result<bool> {
        Ok(ffi::tautomer_callback_call(
            self.raw(),
            mol.as_ffi(),
            res.raw(),
        )?)
    }
}

/// Enumerates a molecule's tautomers and picks a canonical one. Mirrors
/// `rdMolStandardize.TautomerEnumerator`.
///
/// `Send`, not `Sync`. Clones share the transforms and the callback.
pub struct TautomerEnumerator {
    /// Never null.
    ptr: UniquePtr<ffi::TautomerEnumerator>,
    /// The callback RDKit holds, which its getter returns only as a pointer.
    callback: Option<TautomerEnumeratorCallback>,
}

/// RDKit's default transforms and options. Mirrors `TautomerEnumerator()`.
impl Default for TautomerEnumerator {
    fn default() -> TautomerEnumerator {
        TautomerEnumerator::from_ffi(expect_ok(ffi::tautomer_enumerator_default()))
    }
}

impl fmt::Debug for TautomerEnumerator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TautomerEnumerator")
            .field("max_tautomers", &self.get_max_tautomers())
            .field("max_transforms", &self.get_max_transforms())
            .field("remove_sp3_stereo", &self.get_remove_sp3_stereo())
            .field("remove_bond_stereo", &self.get_remove_bond_stereo())
            .field("reassign_stereo", &self.get_reassign_stereo())
            .field("callback", &self.callback.is_some())
            .finish_non_exhaustive()
    }
}

/// Mirrors `TautomerEnumerator(other)`.
impl Clone for TautomerEnumerator {
    fn clone(&self) -> TautomerEnumerator {
        TautomerEnumerator {
            ptr: expect_ok(ffi::tautomer_enumerator_clone(self.raw())),
            callback: self.callback.clone(),
        }
    }
}

impl TautomerEnumerator {
    fn from_ffi(ptr: UniquePtr<ffi::TautomerEnumerator>) -> TautomerEnumerator {
        assert!(!ptr.is_null(), "RDKit returned a null TautomerEnumerator");
        TautomerEnumerator {
            ptr,
            callback: None,
        }
    }

    fn raw(&self) -> &ffi::TautomerEnumerator {
        self.ptr.as_ref().unwrap()
    }

    /// The transforms and options the parameters describe. `Error::Os` for
    /// a transforms file that cannot be read, `Error::Value` for SMARTS that
    /// do not parse, and `Error::Value` for a transform with fewer bonds or
    /// charges than its SMARTS has bonds or atoms, which RDKit would read
    /// past their ends. Mirrors `TautomerEnumerator(params)`.
    pub fn new(params: &CleanupParameters) -> Result<TautomerEnumerator> {
        Ok(TautomerEnumerator::from_ffi(ffi::tautomer_enumerator_new(
            &params.to_ffi(),
        )?))
    }

    /// The tautomers of a molecule. Atoms with the `_protected` property
    /// are not changed. A callback's panic is `Error::Runtime`. Mirrors
    /// `TautomerEnumerator.Enumerate`.
    pub fn enumerate(&self, mol: &Mol) -> Result<TautomerEnumeratorResult> {
        Ok(TautomerEnumeratorResult::from_ffi(
            ffi::tautomer_enumerator_enumerate(self.raw(), mol.as_ffi())?,
        ))
    }

    /// The canonical tautomer, by RDKit's default score: the highest
    /// scoring, then the one with the most stereo, then the first SMILES.
    /// `Error::Value` for a molecule whose scores could overflow RDKit's
    /// `int`s (millions of atoms, or of rings). Mirrors
    /// `TautomerEnumerator.Canonicalize`.
    pub fn canonicalize(&self, mol: &Mol) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::tautomer_enumerator_canonicalize(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// By `score_func`, which gets a copy of each tautomer. A score of
    /// `i32::MIN` counts as `i32::MIN + 1`: RDKit dereferences null when no
    /// tautomer scores above `i32::MIN`. A panic in it is `Error::Runtime`.
    /// Mirrors `TautomerEnumerator.Canonicalize(mol, scoreFunc)`.
    pub fn canonicalize_with_score_func<F>(&self, mol: &Mol, score_func: F) -> Result<Mol>
    where
        F: Fn(&Mol) -> i32,
    {
        let score = ScoreFunction(&score_func);
        Ok(Mol::from_ffi(ffi::tautomer_enumerator_canonicalize_with(
            self.raw(),
            mol.as_ffi(),
            &score,
        )?))
    }

    /// The canonical tautomer of these, by RDKit's default score.
    /// `Error::Value` for none, where RDKit dereferences null, and as
    /// `canonicalize`. RDKit gives each molecule the symmetrized SSSR when
    /// it lacks it. Mirrors `TautomerEnumerator.PickCanonical(iterable)` for
    /// molecules.
    pub fn pick_canonical(&self, mols: &[Mol]) -> Result<Mol> {
        Ok(
            with_borrowed_mols_shared(mols.iter().map(|mol| ("", mol)), |list| {
                ffi::tautomer_enumerator_pick_canonical(self.raw(), list)
            })
            .map(Mol::from_ffi)?,
        )
    }

    /// By `score_func`, as `canonicalize_with_score_func`. Mirrors
    /// `TautomerEnumerator.PickCanonical(iterable, scoreFunc)` for
    /// molecules.
    pub fn pick_canonical_with_score_func<F>(&self, mols: &[Mol], score_func: F) -> Result<Mol>
    where
        F: Fn(&Mol) -> i32,
    {
        let score = ScoreFunction(&score_func);
        Ok(
            with_borrowed_mols_shared(mols.iter().map(|mol| ("", mol)), |list| {
                ffi::tautomer_enumerator_pick_canonical_with(self.raw(), list, &score)
            })
            .map(Mol::from_ffi)?,
        )
    }

    /// The canonical tautomer of a result, by RDKit's default score, as
    /// `pick_canonical`. Mirrors `TautomerEnumerator.PickCanonical(iterable)`
    /// for a `TautomerEnumeratorResult` (A5).
    pub fn pick_canonical_tautomer_enumerator_result(
        &self,
        result: &TautomerEnumeratorResult,
    ) -> Result<Mol> {
        Ok(Mol::from_ffi(
            ffi::tautomer_enumerator_pick_canonical_result(self.raw(), result.raw())?,
        ))
    }

    /// By `score_func`, as `canonicalize_with_score_func`. Mirrors
    /// `TautomerEnumerator.PickCanonical(iterable, scoreFunc)` for a
    /// `TautomerEnumeratorResult`.
    pub fn pick_canonical_tautomer_enumerator_result_with_score_func<F>(
        &self,
        result: &TautomerEnumeratorResult,
        score_func: F,
    ) -> Result<Mol>
    where
        F: Fn(&Mol) -> i32,
    {
        let score = ScoreFunction(&score_func);
        Ok(Mol::from_ffi(
            ffi::tautomer_enumerator_pick_canonical_result_with(self.raw(), result.raw(), &score)?,
        ))
    }

    /// RDKit's default score of a tautomer: `score_rings`, plus
    /// `score_substructs`, plus `score_hetero_hs`. `Error::Value` as
    /// `canonicalize`. Mirrors `TautomerEnumerator.ScoreTautomer`.
    pub fn score_tautomer(mol: &Mol) -> Result<i32> {
        Ok(ffi::score_tautomer(mol.as_ffi())?)
    }

    /// The version of the default score. Mirrors
    /// `TautomerEnumerator.tautomerScoreVersion`.
    pub fn tautomer_score_version() -> String {
        to_string(expect_ok(ffi::tautomer_score_version())).expect("the version is UTF-8")
    }

    /// The most tautomers to enumerate. Mirrors
    /// `TautomerEnumerator.GetMaxTautomers`.
    pub fn get_max_tautomers(&self) -> u32 {
        expect_ok(ffi::tautomer_enumerator_get_max_tautomers(self.raw()))
    }

    /// Mirrors `TautomerEnumerator.SetMaxTautomers`.
    pub fn set_max_tautomers(&mut self, max_tautomers: u32) {
        expect_ok(ffi::tautomer_enumerator_set_max_tautomers(
            self.ptr.pin_mut(),
            max_tautomers,
        ))
    }

    /// The most transforms to apply. Mirrors
    /// `TautomerEnumerator.GetMaxTransforms`.
    pub fn get_max_transforms(&self) -> u32 {
        expect_ok(ffi::tautomer_enumerator_get_max_transforms(self.raw()))
    }

    /// Mirrors `TautomerEnumerator.SetMaxTransforms`.
    pub fn set_max_transforms(&mut self, max_transforms: u32) {
        expect_ok(ffi::tautomer_enumerator_set_max_transforms(
            self.ptr.pin_mut(),
            max_transforms,
        ))
    }

    /// Whether sp3 centers involved in tautomerism lose their stereo.
    /// Mirrors `TautomerEnumerator.GetRemoveSp3Stereo`.
    pub fn get_remove_sp3_stereo(&self) -> bool {
        expect_ok(ffi::tautomer_enumerator_get_remove_sp3_stereo(self.raw()))
    }

    /// Mirrors `TautomerEnumerator.SetRemoveSp3Stereo`.
    pub fn set_remove_sp3_stereo(&mut self, remove_sp3_stereo: bool) {
        expect_ok(ffi::tautomer_enumerator_set_remove_sp3_stereo(
            self.ptr.pin_mut(),
            remove_sp3_stereo,
        ))
    }

    /// Whether double bonds involved in tautomerism lose their stereo.
    /// Mirrors `TautomerEnumerator.GetRemoveBondStereo`.
    pub fn get_remove_bond_stereo(&self) -> bool {
        expect_ok(ffi::tautomer_enumerator_get_remove_bond_stereo(self.raw()))
    }

    /// Mirrors `TautomerEnumerator.SetRemoveBondStereo`.
    pub fn set_remove_bond_stereo(&mut self, remove_bond_stereo: bool) {
        expect_ok(ffi::tautomer_enumerator_set_remove_bond_stereo(
            self.ptr.pin_mut(),
            remove_bond_stereo,
        ))
    }

    /// Whether `enumerate` assigns stereochemistry to each tautomer.
    /// Mirrors `TautomerEnumerator.GetReassignStereo`.
    pub fn get_reassign_stereo(&self) -> bool {
        expect_ok(ffi::tautomer_enumerator_get_reassign_stereo(self.raw()))
    }

    /// Mirrors `TautomerEnumerator.SetReassignStereo`.
    pub fn set_reassign_stereo(&mut self, reassign_stereo: bool) {
        expect_ok(ffi::tautomer_enumerator_set_reassign_stereo(
            self.ptr.pin_mut(),
            reassign_stereo,
        ))
    }

    /// Sets the callback, shared, or with `None` clears it. Mirrors
    /// `TautomerEnumerator.SetCallback`.
    pub fn set_callback(&mut self, callback: Option<&TautomerEnumeratorCallback>) {
        match callback {
            Some(callback) => expect_ok(ffi::tautomer_enumerator_set_callback(
                self.ptr.pin_mut(),
                callback.raw(),
            )),
            None => expect_ok(ffi::tautomer_enumerator_clear_callback(self.ptr.pin_mut())),
        }
        // The old callback is dropped after RDKit lets it go.
        self.callback = callback.cloned();
    }

    /// The callback, shared, if one is set. Mirrors
    /// `TautomerEnumerator.GetCallback`.
    pub fn get_callback(&self) -> Option<TautomerEnumeratorCallback> {
        self.callback.clone()
    }
}

/// An enumerator with the first version of RDKit's transforms. Mirrors
/// `rdMolStandardize.GetV1TautomerEnumerator`.
pub fn get_v1_tautomer_enumerator() -> TautomerEnumerator {
    TautomerEnumerator::from_ffi(expect_ok(ffi::tautomer_enumerator_v1()))
}
