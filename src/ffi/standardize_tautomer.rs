// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::panic::{AssertUnwindSafe, catch_unwind};

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::chem::mol_standardize::rd_mol_standardize::TautomerEnumeratorResult;
use crate::ffi::mol::ffi::RWMol;
use crate::ffi::standardize_validate::panic_message;

/// A tautomer enumerator callback's closure: the molecule and the result so
/// far, to whether to go on.
pub(crate) type TautomerCallbackFn = dyn Fn(&Mol, &TautomerEnumeratorResult) -> bool + Send + Sync;

/// The closure behind a tautomer enumerator callback made in Rust.
/// Enumerators copied from one another share it, so threads holding them
/// may call it at once.
pub(crate) struct TautomerCallback(pub(crate) Box<TautomerCallbackFn>);

/// A scoring function, lent to C++ for one call.
pub(crate) struct ScoreFunction<'a>(pub(crate) &'a dyn Fn(&Mol) -> i32);

/// Called by the C++ callback with its own copy of the molecule and of the
/// result so far. A panic is caught here, at the boundary, and becomes the
/// exception C++ sees.
fn call_tautomer_callback(
    callback: &TautomerCallback,
    mol: UniquePtr<RWMol>,
    result: UniquePtr<ffi::TautomerEnumeratorResult>,
) -> Result<bool, String> {
    let mol = Mol::from_ffi(mol);
    let result = TautomerEnumeratorResult::from_ffi(result);
    catch_unwind(AssertUnwindSafe(|| (callback.0)(&mol, &result)))
        .map_err(|payload| panic_message(payload, "a tautomer enumerator callback"))
}

/// Called by the C++ scoring function with its own copy of the tautomer.
///
/// # Safety
///
/// cxx makes a function whose signature names a lifetime `unsafe` for C++:
/// the C++ caller must not keep `score` past the call that lent it, which
/// the wrapper's scorer does not.
unsafe fn call_score_function(
    score: &ScoreFunction<'_>,
    mol: UniquePtr<RWMol>,
) -> Result<i32, String> {
    let mol = Mol::from_ffi(mol);
    catch_unwind(AssertUnwindSafe(|| (score.0)(&mol)))
        .map_err(|payload| panic_message(payload, "a tautomer scoring function"))
}

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::MolStandardize::TautomerEnumeratorStatus`; cxx checks the
    /// values.
    #[repr(i32)]
    enum TautomerEnumeratorStatus {
        Completed = 0,
        MaxTautomersReached,
        MaxTransformsReached,
        Canceled,
    }

    extern "Rust" {
        type TautomerCallback;
        type ScoreFunction<'a>;

        fn call_tautomer_callback(
            callback: &TautomerCallback,
            mol: UniquePtr<RWMol>,
            result: UniquePtr<TautomerEnumeratorResult>,
        ) -> Result<bool>;
        unsafe fn call_score_function<'a>(
            score: &ScoreFunction<'a>,
            mol: UniquePtr<RWMol>,
        ) -> Result<i32>;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/standardize_tautomer.h");
        include!("rdkit-cxx/ffi/standardize.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type CleanupParametersData = crate::ffi::standardize::ffi::CleanupParametersData;
        type TautomerEnumerator;
        type TautomerEnumeratorResult;
        type TautomerEnumeratorStatus;
        type TautomerCallbackHandle;

        fn tautomer_enumerator_default() -> Result<UniquePtr<TautomerEnumerator>>;
        fn tautomer_enumerator_new(
            params: &CleanupParametersData,
        ) -> Result<UniquePtr<TautomerEnumerator>>;
        fn tautomer_enumerator_v1() -> Result<UniquePtr<TautomerEnumerator>>;
        fn tautomer_enumerator_clone(
            enumerator: &TautomerEnumerator,
        ) -> Result<UniquePtr<TautomerEnumerator>>;

        fn tautomer_enumerator_get_max_tautomers(enumerator: &TautomerEnumerator) -> Result<u32>;
        fn tautomer_enumerator_set_max_tautomers(
            enumerator: Pin<&mut TautomerEnumerator>,
            max_tautomers: u32,
        ) -> Result<()>;
        fn tautomer_enumerator_get_max_transforms(enumerator: &TautomerEnumerator) -> Result<u32>;
        fn tautomer_enumerator_set_max_transforms(
            enumerator: Pin<&mut TautomerEnumerator>,
            max_transforms: u32,
        ) -> Result<()>;
        fn tautomer_enumerator_get_remove_sp3_stereo(
            enumerator: &TautomerEnumerator,
        ) -> Result<bool>;
        fn tautomer_enumerator_set_remove_sp3_stereo(
            enumerator: Pin<&mut TautomerEnumerator>,
            value: bool,
        ) -> Result<()>;
        fn tautomer_enumerator_get_remove_bond_stereo(
            enumerator: &TautomerEnumerator,
        ) -> Result<bool>;
        fn tautomer_enumerator_set_remove_bond_stereo(
            enumerator: Pin<&mut TautomerEnumerator>,
            value: bool,
        ) -> Result<()>;
        fn tautomer_enumerator_get_reassign_stereo(enumerator: &TautomerEnumerator)
        -> Result<bool>;
        fn tautomer_enumerator_set_reassign_stereo(
            enumerator: Pin<&mut TautomerEnumerator>,
            value: bool,
        ) -> Result<()>;
        fn tautomer_enumerator_set_callback(
            enumerator: Pin<&mut TautomerEnumerator>,
            callback: &TautomerCallbackHandle,
        ) -> Result<()>;
        fn tautomer_enumerator_clear_callback(
            enumerator: Pin<&mut TautomerEnumerator>,
        ) -> Result<()>;

        fn tautomer_enumerator_enumerate(
            enumerator: &TautomerEnumerator,
            mol: &RWMol,
        ) -> Result<UniquePtr<TautomerEnumeratorResult>>;
        fn tautomer_enumerator_canonicalize(
            enumerator: &TautomerEnumerator,
            mol: &RWMol,
        ) -> Result<UniquePtr<RWMol>>;
        fn tautomer_enumerator_canonicalize_with<'a>(
            enumerator: &TautomerEnumerator,
            mol: &RWMol,
            score: &ScoreFunction<'a>,
        ) -> Result<UniquePtr<RWMol>>;
        fn tautomer_enumerator_pick_canonical(
            enumerator: &TautomerEnumerator,
            mols: Pin<&mut BorrowedMols>,
        ) -> Result<UniquePtr<RWMol>>;
        fn tautomer_enumerator_pick_canonical_with<'a>(
            enumerator: &TautomerEnumerator,
            mols: Pin<&mut BorrowedMols>,
            score: &ScoreFunction<'a>,
        ) -> Result<UniquePtr<RWMol>>;
        fn tautomer_enumerator_pick_canonical_result(
            enumerator: &TautomerEnumerator,
            result: &TautomerEnumeratorResult,
        ) -> Result<UniquePtr<RWMol>>;
        fn tautomer_enumerator_pick_canonical_result_with<'a>(
            enumerator: &TautomerEnumerator,
            result: &TautomerEnumeratorResult,
            score: &ScoreFunction<'a>,
        ) -> Result<UniquePtr<RWMol>>;

        fn tautomer_callback_new(
            callback: Box<TautomerCallback>,
        ) -> Result<UniquePtr<TautomerCallbackHandle>>;
        fn tautomer_callback_clone(
            handle: &TautomerCallbackHandle,
        ) -> Result<UniquePtr<TautomerCallbackHandle>>;
        fn tautomer_callback_call(
            handle: &TautomerCallbackHandle,
            mol: &RWMol,
            result: &TautomerEnumeratorResult,
        ) -> Result<bool>;

        fn tautomer_result_size(result: &TautomerEnumeratorResult) -> Result<usize>;
        fn tautomer_result_tautomer(
            result: &TautomerEnumeratorResult,
            i: usize,
        ) -> Result<UniquePtr<RWMol>>;
        fn tautomer_result_kekulized(
            result: &TautomerEnumeratorResult,
            i: usize,
        ) -> Result<UniquePtr<RWMol>>;
        fn tautomer_result_smiles(
            result: &TautomerEnumeratorResult,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn tautomer_result_status(
            result: &TautomerEnumeratorResult,
        ) -> Result<TautomerEnumeratorStatus>;
        fn tautomer_result_modified_atoms(result: &TautomerEnumeratorResult) -> Result<Vec<u32>>;
        fn tautomer_result_modified_bonds(result: &TautomerEnumeratorResult) -> Result<Vec<u32>>;

        fn score_tautomer(mol: &RWMol) -> Result<i32>;
        fn tautomer_score_version() -> Result<UniquePtr<CxxString>>;
    }
}

// SAFETY: an enumerator shares its catalog, which nothing changes after it
// is made, and its callback, a Send + Sync closure, with its copies, through
// shared_ptrs (atomic counts); the enumeration only reads them.
unsafe impl Send for ffi::TautomerEnumerator {}
// SAFETY: a result from enumerate is the only owner of its tautomers.
// Callbacks are lent a copy that shares them with RDKit's enumeration, but
// only by reference, for the call, on the enumerating thread (the Rust type
// is not Sync and has no Clone).
unsafe impl Send for ffi::TautomerEnumeratorResult {}
// SAFETY: a handle shares a callback whose closure is Send + Sync.
unsafe impl Send for ffi::TautomerCallbackHandle {}
unsafe impl Sync for ffi::TautomerCallbackHandle {}
