// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::any::Any;
use std::cell::{Cell, RefCell};
use std::panic::{AssertUnwindSafe, catch_unwind};

use cxx::CxxString;

use crate::chem::Mol;
use crate::chem::rd_fmcs::{
    AtomTyper, BondTyper, MCSAtomCompareParameters, MCSBondCompareParameters, MCSParameters,
    MCSProgressData,
};
use crate::error::{Error, Result, decode};
use crate::ffi::chirality::ffi::IndexPair;

/// The Rust side of one `FindMCS` call: the molecules RDKit compares, the
/// parameters whose closures the C++ trampolines call, and the first error
/// a closure returned (or a panic it raised). The trampolines check
/// `failed` before each call: once a closure has failed, the rest answer
/// false without running, the progress callback cancels the search, and
/// `find_mcs_with_params` returns the error.
pub(crate) struct McsCallbacks {
    /// The molecules of the call, which outlive this value: it is made and
    /// dropped inside `with_callbacks`, while the slice is borrowed.
    mols: *const [Mol],
    params: MCSParameters,
    failed: Cell<bool>,
    error: RefCell<Option<Error>>,
}

impl McsCallbacks {
    /// Runs `f` with the callbacks for a call on `mols`, then returns the
    /// first error a callback recorded, if any, in place of `f`'s value.
    pub(crate) fn with_callbacks<T>(
        mols: &[Mol],
        params: &MCSParameters,
        f: impl FnOnce(&McsCallbacks) -> T,
    ) -> std::result::Result<T, Error> {
        let callbacks = McsCallbacks {
            mols: mols as *const [Mol],
            params: params.clone(),
            failed: Cell::new(false),
            error: RefCell::new(None),
        };
        let value = f(&callbacks);
        if !callbacks.failed.get() {
            return Ok(value);
        }
        Err(callbacks.error.into_inner().unwrap_or(Error::Runtime {
            message: "a FindMCS callback failed".to_string(),
        }))
    }

    fn mols(&self) -> &[Mol] {
        // SAFETY: `with_callbacks` makes this value from a borrowed slice and
        // drops it before the borrow ends, and the trampolines only call in
        // while it lives (during `find_mcs`).
        unsafe { &*self.mols }
    }

    fn mol(&self, index: usize) -> Result<&Mol> {
        self.mols().get(index).ok_or_else(|| Error::Runtime {
            message: format!("FindMCS called back with molecule {index}, which it was not given"),
        })
    }

    fn fail(&self, error: Error) {
        if self.failed.replace(true) {
            return;
        }
        if let Ok(mut slot) = self.error.try_borrow_mut() {
            *slot = Some(error);
        }
    }

    /// Runs one callback, catching a panic at the boundary: this is called
    /// from C++, so it must not unwind.
    fn run(&self, what: &str, f: impl FnOnce() -> Result<bool>) -> bool {
        if self.failed.get() {
            return false;
        }
        match catch_unwind(AssertUnwindSafe(f)) {
            Ok(Ok(value)) => value,
            Ok(Err(error)) => {
                self.fail(error);
                false
            }
            Err(payload) => {
                let message = panic_message(payload.as_ref());
                self.fail(Error::Runtime {
                    message: format!("{what} panicked: {message}"),
                });
                false
            }
        }
    }
}

fn panic_message(payload: &(dyn Any + Send)) -> String {
    payload
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "a panic without a message".to_string())
}

fn pairs(pairs: &[IndexPair]) -> Vec<(u32, u32)> {
    pairs.iter().map(|pair| (pair.first, pair.second)).collect()
}

fn mcs_callbacks_failed(callbacks: &McsCallbacks) -> bool {
    callbacks.failed.get()
}

/// Records an exception a trampoline caught, encoded as `errors.h` encodes
/// any exception.
fn mcs_callbacks_fail(callbacks: &McsCallbacks, what: &CxxString) {
    let error = catch_unwind(AssertUnwindSafe(|| decode(&what.to_string_lossy()))).unwrap_or(
        Error::Runtime {
            message: "a FindMCS callback failed".to_string(),
        },
    );
    callbacks.fail(error);
}

fn mcs_call_atom_typer(
    callbacks: &McsCallbacks,
    params: &ffi::MCSAtomCompareParametersData,
    mol1: usize,
    atom1: u32,
    mol2: usize,
    atom2: u32,
) -> bool {
    callbacks.run("an MCSAtomCompare", || {
        let AtomTyper::Custom(typer) = &callbacks.params.atom_typer else {
            return Err(Error::Runtime {
                message: "FindMCS called an atom typer Rust did not set".to_string(),
            });
        };
        let params = MCSAtomCompareParameters::from_ffi(params);
        typer.call(
            &params,
            callbacks.mol(mol1)?,
            atom1,
            callbacks.mol(mol2)?,
            atom2,
        )
    })
}

fn mcs_call_bond_typer(
    callbacks: &McsCallbacks,
    params: &ffi::MCSBondCompareParametersData,
    mol1: usize,
    bond1: u32,
    mol2: usize,
    bond2: u32,
) -> bool {
    callbacks.run("an MCSBondCompare", || {
        let BondTyper::Custom(typer) = &callbacks.params.bond_typer else {
            return Err(Error::Runtime {
                message: "FindMCS called a bond typer Rust did not set".to_string(),
            });
        };
        let params = MCSBondCompareParameters::from_ffi(params);
        typer.call(
            &params,
            callbacks.mol(mol1)?,
            bond1,
            callbacks.mol(mol2)?,
            bond2,
        )
    })
}

fn mcs_call_progress(
    callbacks: &McsCallbacks,
    num_atoms: u32,
    num_bonds: u32,
    seed_processed: u32,
    params: &ffi::MCSParametersData,
) -> bool {
    callbacks.run("an MCSProgress", || {
        let Some(progress) = &callbacks.params.progress_callback else {
            return Err(Error::Runtime {
                message: "FindMCS called a progress callback Rust did not set".to_string(),
            });
        };
        let stat = MCSProgressData::new(num_atoms, num_bonds, seed_processed);
        progress.call(
            &stat,
            &MCSParameters::from_callback(params, &callbacks.params),
        )
    })
}

fn mcs_call_final_match_checker(
    callbacks: &McsCallbacks,
    mol1: usize,
    mol2: usize,
    atom_idx_match: &[IndexPair],
    bond_idx_match: &[IndexPair],
    params: &ffi::MCSParametersData,
) -> bool {
    callbacks.run("an MCSFinalMatchCheck", || {
        let Some(check) = &callbacks.params.final_match_checker else {
            return Err(Error::Runtime {
                message: "FindMCS called a final match checker Rust did not set".to_string(),
            });
        };
        check.call(
            callbacks.mol(mol1)?,
            callbacks.mol(mol2)?,
            &pairs(atom_idx_match),
            &pairs(bond_idx_match),
            &MCSParameters::from_callback(params, &callbacks.params),
        )
    })
}

fn mcs_call_should_accept_mcs(
    callbacks: &McsCallbacks,
    query: usize,
    target: usize,
    atom_idx_match: &[IndexPair],
    bond_idx_match: &[IndexPair],
    params: &ffi::MCSParametersData,
) -> bool {
    callbacks.run("an MCSAcceptance", || {
        let Some(accept) = &callbacks.params.should_accept_mcs else {
            return Err(Error::Runtime {
                message: "FindMCS called an acceptance callback Rust did not set".to_string(),
            });
        };
        accept.call(
            callbacks.mol(query)?,
            callbacks.mol(target)?,
            &pairs(atom_idx_match),
            &pairs(bond_idx_match),
            &MCSParameters::from_callback(params, &callbacks.params),
        )
    })
}

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::AtomComparator`; cxx checks the values.
    #[repr(u32)]
    enum AtomComparator {
        AtomCompareAny = 0,
        AtomCompareElements,
        AtomCompareIsotopes,
        AtomCompareAnyHeavyAtom,
    }

    /// `RDKit::BondComparator`; cxx checks the values.
    #[repr(u32)]
    enum BondComparator {
        BondCompareAny = 0,
        BondCompareOrder,
        BondCompareOrderExact,
    }

    /// `RDKit::MCSAtomCompareParameters`.
    struct MCSAtomCompareParametersData {
        match_valences: bool,
        match_chiral_tag: bool,
        match_formal_charge: bool,
        ring_matches_ring_only: bool,
        complete_rings_only: bool,
        match_isotope: bool,
        max_distance: f64,
    }

    /// `RDKit::MCSBondCompareParameters`.
    struct MCSBondCompareParametersData {
        ring_matches_ring_only: bool,
        complete_rings_only: bool,
        match_fused_rings: bool,
        match_fused_rings_strict: bool,
        match_stereo: bool,
    }

    /// `RDKit::MCSParameters`. Each callback is a flag: set, the C++ side
    /// installs a trampoline that calls the closure through `McsCallbacks`.
    /// From C++ (defaults, and the parameters a callback sees), a typer
    /// flag is set when the typer is not one of RDKit's.
    struct MCSParametersData {
        maximize_bonds: bool,
        threshold: f64,
        timeout: u32,
        verbose: bool,
        atom_compare_parameters: MCSAtomCompareParametersData,
        bond_compare_parameters: MCSBondCompareParametersData,
        atom_typer: AtomComparator,
        custom_atom_typer: bool,
        bond_typer: BondComparator,
        custom_bond_typer: bool,
        progress_callback: bool,
        final_match_checker: bool,
        should_accept_mcs: bool,
        /// Comes from Rust, which RDKit copies but never changes, so it is
        /// UTF-8 on the way back too.
        initial_seed: String,
        store_all: bool,
    }

    /// `RDKit::MCSResult`, with the molecules copied out.
    struct MCSResultData {
        num_atoms: u32,
        num_bonds: u32,
        smarts_string: UniquePtr<CxxString>,
        canceled: bool,
        /// Null when the MCS is empty or `StoreAll` is set.
        query_mol: UniquePtr<RWMol>,
        /// Each degenerate MCS's query molecule, named by its SMARTS.
        degenerate_smarts_query_mol_dict: UniquePtr<MolList>,
    }

    extern "Rust" {
        type McsCallbacks;

        fn mcs_callbacks_failed(callbacks: &McsCallbacks) -> bool;
        fn mcs_callbacks_fail(callbacks: &McsCallbacks, what: &CxxString);
        fn mcs_call_atom_typer(
            callbacks: &McsCallbacks,
            params: &MCSAtomCompareParametersData,
            mol1: usize,
            atom1: u32,
            mol2: usize,
            atom2: u32,
        ) -> bool;
        fn mcs_call_bond_typer(
            callbacks: &McsCallbacks,
            params: &MCSBondCompareParametersData,
            mol1: usize,
            bond1: u32,
            mol2: usize,
            bond2: u32,
        ) -> bool;
        fn mcs_call_progress(
            callbacks: &McsCallbacks,
            num_atoms: u32,
            num_bonds: u32,
            seed_processed: u32,
            params: &MCSParametersData,
        ) -> bool;
        fn mcs_call_final_match_checker(
            callbacks: &McsCallbacks,
            mol1: usize,
            mol2: usize,
            atom_idx_match: &[IndexPair],
            bond_idx_match: &[IndexPair],
            params: &MCSParametersData,
        ) -> bool;
        fn mcs_call_should_accept_mcs(
            callbacks: &McsCallbacks,
            query: usize,
            target: usize,
            atom_idx_match: &[IndexPair],
            bond_idx_match: &[IndexPair],
            params: &MCSParametersData,
        ) -> bool;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/fmcs.h");
        include!("rdkit-cxx/ffi/chirality.rs.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type MolList = crate::ffi::molops::ffi::MolList;
        type BorrowedMols = crate::ffi::molops::ffi::BorrowedMols;
        type IndexPair = crate::ffi::chirality::ffi::IndexPair;
        type AtomComparator;
        type BondComparator;

        fn mcs_atom_compare_parameters_default() -> Result<MCSAtomCompareParametersData>;
        fn mcs_bond_compare_parameters_default() -> Result<MCSBondCompareParametersData>;
        fn mcs_parameters_default() -> Result<MCSParametersData>;

        fn find_mcs(
            mols: Pin<&mut BorrowedMols>,
            params: &MCSParametersData,
            callbacks: &McsCallbacks,
        ) -> Result<MCSResultData>;

        fn mcs_atom_compare_check_atom_ring_match(
            params: &MCSAtomCompareParametersData,
            mol1: &RWMol,
            atom1: u32,
            mol2: &RWMol,
            atom2: u32,
        ) -> Result<bool>;
        fn mcs_atom_compare_check_atom_charge(
            params: &MCSAtomCompareParametersData,
            mol1: &RWMol,
            atom1: u32,
            mol2: &RWMol,
            atom2: u32,
        ) -> Result<bool>;
        fn mcs_atom_compare_check_atom_chirality(
            params: &MCSAtomCompareParametersData,
            mol1: &RWMol,
            atom1: u32,
            mol2: &RWMol,
            atom2: u32,
        ) -> Result<bool>;
        fn mcs_bond_compare_check_bond_stereo(
            params: &MCSBondCompareParametersData,
            mol1: &RWMol,
            bond1: u32,
            mol2: &RWMol,
            bond2: u32,
        ) -> Result<bool>;
        fn mcs_bond_compare_check_bond_ring_match(
            params: &MCSBondCompareParametersData,
            mol1: &RWMol,
            bond1: u32,
            mol2: &RWMol,
            bond2: u32,
        ) -> Result<bool>;
    }
}
