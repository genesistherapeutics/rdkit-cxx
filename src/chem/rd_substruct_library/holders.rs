// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The holders a library keeps its molecules, fingerprints, and keys in.
//!
//! Each holder class of Python's is a type that dereferences to its base
//! class's type, which a library holds. `downcast_ref` and `downcast_mut`
//! find the class again, as Python does for the holders it returns.

use std::fmt;
use std::ops::Deref;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::data_structs::c_data_structs::ExplicitBitVect;
use crate::error::{expect_ok, to_string};
use crate::ffi::substruct_library::ffi;

mod sealed {
    pub trait Sealed {}
}

/// A molecule holder's class: `MolHolder`, `CachedMolHolder`,
/// `CachedSmilesMolHolder`, or `CachedTrustedSmilesMolHolder`.
/// A Rust addition, for Python's duck typing; sealed.
pub trait MolHolderType: sealed::Sealed {
    #[doc(hidden)]
    const KIND: u8;
}

/// A fingerprint holder's class: `PatternHolder` or `TautomerPatternHolder`.
/// A Rust addition, for Python's duck typing; sealed.
pub trait FPHolderType: sealed::Sealed {
    #[doc(hidden)]
    const KIND: u8;
}

/// A key holder's class: `KeyFromPropHolder`.
/// A Rust addition, for Python's duck typing; sealed.
pub trait KeyHolderType: sealed::Sealed {}

/// Molecules for a library, of any of the holder classes. Mirrors
/// `rdSubstructLibrary.MolHolderBase`.
///
/// Holders are `Send` but not `Sync`: searching a library writes the ring
/// information of a `MolHolder`'s molecules (CONVENTIONS.md S1).
pub struct MolHolderBase {
    /// Never null.
    ptr: UniquePtr<ffi::MolHolderHandle>,
}

impl MolHolderBase {
    fn with_kind(kind: u8) -> MolHolderBase {
        MolHolderBase::from_ffi(expect_ok(ffi::mol_holder_new(kind)))
    }

    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::MolHolderHandle>) -> MolHolderBase {
        assert!(!ptr.is_null(), "RDKit returned a null holder");
        MolHolderBase { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::MolHolderHandle {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::MolHolderHandle> {
        self.ptr.pin_mut()
    }

    fn kind(&self) -> u8 {
        expect_ok(ffi::mol_holder_kind(self.raw()))
    }

    /// The number of molecules. Mirrors `MolHolderBase.__len__`.
    pub fn len(&self) -> u32 {
        expect_ok(ffi::mol_holder_len(self.raw()))
    }

    /// Whether the holder has no molecules.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Adds a molecule (a copy, a pickle, or SMILES, by the holder's class)
    /// and returns its index. A holder of `i32::MAX` molecules is full
    /// (`Error::Value`): RDKit counts a library's matches in `int`. A
    /// `CachedMolHolder` refuses a molecule `Mol::to_binary` refuses
    /// (`Error::Value`), whose pickle it could not read back. Mirrors
    /// `MolHolderBase.AddMol`.
    pub fn add_mol(&mut self, m: &Mol) -> Result<u32> {
        Ok(ffi::mol_holder_add_mol(self.raw_mut(), m.as_ffi())?)
    }

    /// A copy of the molecule at `idx`. A stored SMILES that does not parse
    /// is `Error::Parse` (Python returns `None`), and one that does not
    /// sanitize the sanitization variant; a stored pickle whose plain
    /// monomer info has the PDB residue type is `Error::Value`
    /// (CONVENTIONS.md U10). The trusted SMILES holder's molecules are not
    /// sanitized and have no ring information. Mirrors
    /// `MolHolderBase.GetMol`.
    pub fn get_mol(&self, idx: u32) -> Result<Mol> {
        Ok(Mol::from_ffi(ffi::mol_holder_get_mol(self.raw(), idx)?))
    }

    /// The holder as its class, if it is a `T`.
    pub fn downcast_ref<T: MolHolderType>(&self) -> Option<&T> {
        if self.kind() != T::KIND {
            return None;
        }
        // SAFETY: every MolHolderType (the trait is sealed) is a
        // #[repr(transparent)] wrapper of a MolHolderBase.
        Some(unsafe { &*(self as *const MolHolderBase).cast::<T>() })
    }

    /// The holder as its class, if it is a `T`.
    pub fn downcast_mut<T: MolHolderType>(&mut self) -> Option<&mut T> {
        if self.kind() != T::KIND {
            return None;
        }
        // SAFETY: as in downcast_ref.
        Some(unsafe { &mut *(self as *mut MolHolderBase).cast::<T>() })
    }
}

impl fmt::Debug for MolHolderBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let class = [
            "MolHolder",
            "CachedMolHolder",
            "CachedSmilesMolHolder",
            "CachedTrustedSmilesMolHolder",
        ][usize::from(self.kind())];
        f.debug_struct(class).field("len", &self.len()).finish()
    }
}

macro_rules! mol_holder {
    ($(#[$meta:meta])* $name:ident, $kind:expr) => {
        $(#[$meta])*
        #[repr(transparent)]
        pub struct $name {
            base: MolHolderBase,
        }

        #[doc = concat!("An empty holder. Mirrors `", stringify!($name), "()`.")]
        impl Default for $name {
            fn default() -> $name {
                $name {
                    base: MolHolderBase::with_kind($kind),
                }
            }
        }

        impl $name {
            /// Mirrors `MolHolderBase.AddMol`, as `MolHolderBase::add_mol`.
            pub fn add_mol(&mut self, m: &Mol) -> Result<u32> {
                self.base.add_mol(m)
            }
        }

        impl Deref for $name {
            type Target = MolHolderBase;

            fn deref(&self) -> &MolHolderBase {
                &self.base
            }
        }

        impl From<$name> for MolHolderBase {
            fn from(holder: $name) -> MolHolderBase {
                holder.base
            }
        }

        impl sealed::Sealed for $name {}

        impl MolHolderType for $name {
            const KIND: u8 = $kind;
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.base.fmt(f)
            }
        }
    };
}

mol_holder!(
    /// Molecules held as copies, the fastest to search and the largest.
    /// Mirrors `rdSubstructLibrary.MolHolder`.
    MolHolder,
    0
);
mol_holder!(
    /// Molecules held as pickles, parsed for each search. Mirrors
    /// `rdSubstructLibrary.CachedMolHolder`.
    CachedMolHolder,
    1
);
mol_holder!(
    /// Molecules held as SMILES, parsed and sanitized for each search.
    /// Mirrors `rdSubstructLibrary.CachedSmilesMolHolder`.
    CachedSmilesMolHolder,
    2
);
mol_holder!(
    /// Molecules held as SMILES RDKit wrote, parsed without sanitization
    /// for each search: only their property caches are updated, and their
    /// rings are found when a query needs them. Other SMILES give wrong
    /// matches, or errors, but nothing worse. Mirrors
    /// `rdSubstructLibrary.CachedTrustedSmilesMolHolder`.
    CachedTrustedSmilesMolHolder,
    3
);

impl CachedMolHolder {
    /// Adds a molecule pickle as it is, unchecked until a search or
    /// `get_mol` parses it, and returns its index. Mirrors
    /// `CachedMolHolder.AddBinary`.
    ///
    /// # Safety
    ///
    /// `pickle` must be a molecule pickle RDKit wrote, unmodified, as for
    /// `Mol::from_pkl_string`: the holder's `get_mol` and every search that
    /// reads it, safe as they are, read it that way.
    pub unsafe fn add_binary(&mut self, pickle: &[u8]) -> Result<u32> {
        Ok(ffi::mol_holder_add_binary(self.base.raw_mut(), pickle)?)
    }
}

impl CachedSmilesMolHolder {
    /// Adds SMILES as they are, unchecked until a search or `get_mol`
    /// parses them, and returns their index. Mirrors
    /// `CachedSmilesMolHolder.AddSmiles`.
    pub fn add_smiles(&mut self, smiles: &str) -> Result<u32> {
        Ok(ffi::mol_holder_add_smiles(self.base.raw_mut(), smiles)?)
    }
}

impl CachedTrustedSmilesMolHolder {
    /// Adds SMILES as they are, which RDKit should have written, and
    /// returns their index. Mirrors `CachedTrustedSmilesMolHolder.AddSmiles`.
    pub fn add_smiles(&mut self, smiles: &str) -> Result<u32> {
        Ok(ffi::mol_holder_add_smiles(self.base.raw_mut(), smiles)?)
    }
}

/// Pattern fingerprints that screen a library's molecules before they are
/// matched, of any of the holder classes. Every fingerprint has the
/// holder's number of bits. Mirrors `rdSubstructLibrary.FPHolderBase`.
pub struct FPHolderBase {
    /// Never null.
    ptr: UniquePtr<ffi::FpHolderHandle>,
}

impl FPHolderBase {
    fn with_kind(kind: u8, num_bits: u32) -> FPHolderBase {
        assert!(
            num_bits <= i32::MAX as u32,
            "a pattern holder has at most i32::MAX bits"
        );
        FPHolderBase::from_ffi(expect_ok(ffi::fp_holder_new(kind, num_bits)))
    }

    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::FpHolderHandle>) -> FPHolderBase {
        assert!(!ptr.is_null(), "RDKit returned a null holder");
        FPHolderBase { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::FpHolderHandle {
        self.ptr.as_ref().unwrap()
    }

    pub(crate) fn raw_mut(&mut self) -> Pin<&mut ffi::FpHolderHandle> {
        self.ptr.pin_mut()
    }

    fn kind(&self) -> u8 {
        expect_ok(ffi::fp_holder_kind(self.raw()))
    }

    /// The fingerprints' number of bits.
    pub fn get_num_bits(&self) -> u32 {
        expect_ok(ffi::fp_holder_num_bits(self.raw()))
    }

    /// The number of fingerprints. Mirrors `FPHolderBase.__len__`.
    pub fn len(&self) -> u32 {
        expect_ok(ffi::fp_holder_len(self.raw()))
    }

    /// Whether the holder has no fingerprints.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Adds a molecule's fingerprint and returns its index; a holder of no
    /// bits is `Error::Runtime`, and a full one `Error::Value`. Mirrors
    /// `FPHolderBase.AddMol`.
    pub fn add_mol(&mut self, m: &Mol) -> Result<u32> {
        Ok(ffi::fp_holder_add_mol(self.raw_mut(), m.as_ffi())?)
    }

    /// Adds a fingerprint (made by `make_fingerprint`, here or elsewhere)
    /// and returns its index. One of another size is `Error::Value`:
    /// RDKit's filter would read past the end of a shorter one. Mirrors
    /// `FPHolderBase.AddFingerprint`.
    pub fn add_fingerprint(&mut self, v: &ExplicitBitVect) -> Result<u32> {
        Ok(ffi::fp_holder_add_fingerprint(self.raw_mut(), v.raw())?)
    }

    /// A copy of the fingerprint at `idx`. Mirrors
    /// `FPHolderBase.GetFingerprint`, which Python returns as a view.
    pub fn get_fingerprint(&self, idx: u32) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(ffi::fp_holder_get_fingerprint(
            self.raw(),
            idx,
        )?))
    }

    /// Whether every bit of `query` is set in the fingerprint at `idx`,
    /// which a match needs. A query of another size is `Error::Value`, where
    /// RDKit reads past the end of the shorter one. Mirrors
    /// `FPHolderBase.PassesFilter`.
    pub fn passes_filter(&self, idx: u32, query: &ExplicitBitVect) -> Result<bool> {
        Ok(ffi::fp_holder_passes_filter(self.raw(), idx, query.raw())?)
    }

    /// The fingerprint the holder makes for a molecule or query. A holder
    /// of no bits is `Error::Runtime`. Mirrors `FPHolderBase.MakeFingerprint`.
    pub fn make_fingerprint(&self, mol: &Mol) -> Result<ExplicitBitVect> {
        Ok(ExplicitBitVect::from_ffi(ffi::fp_holder_make_fingerprint(
            self.raw(),
            mol.as_ffi(),
        )?))
    }

    /// The holder as its class, if it is a `T`.
    pub fn downcast_ref<T: FPHolderType>(&self) -> Option<&T> {
        if self.kind() != T::KIND {
            return None;
        }
        // SAFETY: every FPHolderType (the trait is sealed) is a
        // #[repr(transparent)] wrapper of an FPHolderBase.
        Some(unsafe { &*(self as *const FPHolderBase).cast::<T>() })
    }

    /// The holder as its class, if it is a `T`.
    pub fn downcast_mut<T: FPHolderType>(&mut self) -> Option<&mut T> {
        if self.kind() != T::KIND {
            return None;
        }
        // SAFETY: as in downcast_ref.
        Some(unsafe { &mut *(self as *mut FPHolderBase).cast::<T>() })
    }
}

impl fmt::Debug for FPHolderBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let class = ["PatternHolder", "TautomerPatternHolder"][usize::from(self.kind())];
        f.debug_struct(class)
            .field("num_bits", &self.get_num_bits())
            .field("len", &self.len())
            .finish()
    }
}

macro_rules! fp_holder {
    ($(#[$meta:meta])* $name:ident, $kind:expr) => {
        $(#[$meta])*
        #[repr(transparent)]
        pub struct $name {
            base: FPHolderBase,
        }

        #[doc = concat!("A holder of 2048-bit fingerprints. Mirrors `", stringify!($name), "()`.")]
        impl Default for $name {
            fn default() -> $name {
                $name::new(2048)
            }
        }

        impl $name {
            #[doc = concat!("A holder of fingerprints with `num_bits` bits. Mirrors `", stringify!($name), "(numBits)`.")]
            ///
            /// # Panics
            ///
            /// For more than `i32::MAX` bits (CONVENTIONS.md U3).
            pub fn new(num_bits: u32) -> $name {
                $name {
                    base: FPHolderBase::with_kind($kind, num_bits),
                }
            }

            /// Mirrors `FPHolderBase.AddMol`, as `FPHolderBase::add_mol`.
            pub fn add_mol(&mut self, m: &Mol) -> Result<u32> {
                self.base.add_mol(m)
            }

            /// Mirrors `FPHolderBase.AddFingerprint`, as
            /// `FPHolderBase::add_fingerprint`.
            pub fn add_fingerprint(&mut self, v: &ExplicitBitVect) -> Result<u32> {
                self.base.add_fingerprint(v)
            }
        }

        impl Deref for $name {
            type Target = FPHolderBase;

            fn deref(&self) -> &FPHolderBase {
                &self.base
            }
        }

        impl From<$name> for FPHolderBase {
            fn from(holder: $name) -> FPHolderBase {
                holder.base
            }
        }

        impl sealed::Sealed for $name {}

        impl FPHolderType for $name {
            const KIND: u8 = $kind;
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.base.fmt(f)
            }
        }
    };
}

fp_holder!(
    /// Pattern fingerprints (`pattern_fingerprint`). Mirrors
    /// `rdSubstructLibrary.PatternHolder`.
    PatternHolder,
    0
);
fp_holder!(
    /// Tautomeric pattern fingerprints, for tautomer queries. Mirrors
    /// `rdSubstructLibrary.TautomerPatternHolder`.
    TautomerPatternHolder,
    1
);

/// Keys for a library's molecules, which name them outside it. Mirrors
/// `rdSubstructLibrary.KeyHolderBase`.
pub struct KeyHolderBase {
    /// Never null.
    ptr: UniquePtr<ffi::KeyHolderHandle>,
}

// Keys come from Rust or from molecules' properties, which may hold invalid
// UTF-8 read from files (T2).
fn text(value: UniquePtr<cxx::CxxString>) -> Result<String> {
    to_string(value)
}

impl KeyHolderBase {
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::KeyHolderHandle>) -> KeyHolderBase {
        assert!(!ptr.is_null(), "RDKit returned a null holder");
        KeyHolderBase { ptr }
    }

    pub(crate) fn raw(&self) -> &ffi::KeyHolderHandle {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::KeyHolderHandle> {
        self.ptr.pin_mut()
    }

    /// The number of keys. Mirrors `KeyHolderBase.__len__`.
    pub fn len(&self) -> u32 {
        expect_ok(ffi::key_holder_len(self.raw()))
    }

    /// Whether the holder has no keys.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Adds a molecule's key, its property's value or `LIBIDX-<index>`
    /// without one, and returns its index. Mirrors `KeyHolderBase.AddMol`.
    pub fn add_mol(&mut self, m: &Mol) -> Result<u32> {
        Ok(ffi::key_holder_add_mol(self.raw_mut(), m.as_ffi())?)
    }

    /// Adds a key, which must come in the molecules' order, and returns its
    /// index. Mirrors `KeyHolderBase.AddKey`.
    pub fn add_key(&mut self, key: &str) -> Result<u32> {
        Ok(ffi::key_holder_add_key(self.raw_mut(), key)?)
    }

    /// The key at `idx`; invalid UTF-8 is `Error::Value`. Mirrors
    /// `KeyHolderBase.GetKey`.
    pub fn get_key(&self, idx: u32) -> Result<String> {
        text(ffi::key_holder_get_key(self.raw(), idx)?)
    }

    /// The keys at the indices, such as a search's matches. Mirrors
    /// `KeyHolderBase.GetKeys`.
    pub fn get_keys(&self, indices: &[u32]) -> Result<Vec<String>> {
        ffi::key_holder_get_keys(self.raw(), indices)?
            .iter()
            .map(|key| {
                key.to_str()
                    .map(str::to_owned)
                    .map_err(|err| crate::Error::Value {
                        message: format!("invalid UTF-8 in a string from RDKit: {err}"),
                    })
            })
            .collect()
    }

    /// The holder as its class, if it is a `T`.
    pub fn downcast_ref<T: KeyHolderType>(&self) -> Option<&T> {
        // SAFETY: KeyFromPropHolder, the only KeyHolderType (the trait is
        // sealed) and the only key holder, is a #[repr(transparent)] wrapper
        // of a KeyHolderBase.
        Some(unsafe { &*(self as *const KeyHolderBase).cast::<T>() })
    }

    /// The holder as its class, if it is a `T`.
    pub fn downcast_mut<T: KeyHolderType>(&mut self) -> Option<&mut T> {
        // SAFETY: as in downcast_ref.
        Some(unsafe { &mut *(self as *mut KeyHolderBase).cast::<T>() })
    }
}

impl fmt::Debug for KeyHolderBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("KeyFromPropHolder")
            .field(
                "prop_name",
                &text(expect_ok(ffi::key_holder_get_prop_name(self.raw()))),
            )
            .field("len", &self.len())
            .finish()
    }
}

/// Keys read from a property of each molecule, `_Name` by default. Mirrors
/// `rdSubstructLibrary.KeyFromPropHolder`.
#[repr(transparent)]
pub struct KeyFromPropHolder {
    base: KeyHolderBase,
}

/// Keys from `_Name`. Mirrors `KeyFromPropHolder()`.
impl Default for KeyFromPropHolder {
    fn default() -> KeyFromPropHolder {
        KeyFromPropHolder::new("_Name")
    }
}

impl KeyFromPropHolder {
    /// Keys from the property `propname`. Mirrors
    /// `KeyFromPropHolder(propname)`.
    pub fn new(propname: &str) -> KeyFromPropHolder {
        KeyFromPropHolder {
            base: KeyHolderBase::from_ffi(expect_ok(ffi::key_holder_new(propname))),
        }
    }

    /// The property the keys come from; invalid UTF-8 (in a library read
    /// from a pickle) is `Error::Value`. Mirrors
    /// `KeyFromPropHolder.GetPropName`.
    pub fn get_prop_name(&self) -> Result<String> {
        text(expect_ok(ffi::key_holder_get_prop_name(self.base.raw())))
    }

    /// Mirrors `KeyHolderBase.AddMol`, as `KeyHolderBase::add_mol`.
    pub fn add_mol(&mut self, m: &Mol) -> Result<u32> {
        self.base.add_mol(m)
    }

    /// Mirrors `KeyHolderBase.AddKey`, as `KeyHolderBase::add_key`.
    pub fn add_key(&mut self, key: &str) -> Result<u32> {
        self.base.add_key(key)
    }
}

impl Deref for KeyFromPropHolder {
    type Target = KeyHolderBase;

    fn deref(&self) -> &KeyHolderBase {
        &self.base
    }
}

impl From<KeyFromPropHolder> for KeyHolderBase {
    fn from(holder: KeyFromPropHolder) -> KeyHolderBase {
        holder.base
    }
}

impl sealed::Sealed for KeyFromPropHolder {}

impl KeyHolderType for KeyFromPropHolder {}

impl fmt::Debug for KeyFromPropHolder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.base.fmt(f)
    }
}
