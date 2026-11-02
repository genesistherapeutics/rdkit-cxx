// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::collections::BTreeMap;
use std::fmt;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::error::{expect_ok, op};
use crate::ffi::data_structs::ffi;
use crate::{Error, Result};

pub(crate) mod sealed {
    use crate::ffi::data_structs::ffi;

    pub trait SparseIntVectLike {
        fn sparse_int_vect(&self) -> &ffi::SparseIntVect;
    }
}

/// A sparse integer vector the module functions take: `IntSparseIntVect`,
/// `LongSparseIntVect`, `UIntSparseIntVect`, or `ULongSparseIntVect`.
/// Generic functions over it cover Python's overloads for all four
/// (CONVENTIONS.md A5).
pub trait SparseIntVectType: sealed::SparseIntVectLike {}

/// An index as the bridge carries it; negative indices are out of range.
fn index<T: TryInto<u64> + TryInto<i64> + Copy>(idx: T, len: u64) -> Result<u64> {
    idx.try_into()
        .map_err(|_| Error::index(idx, usize::try_from(len).unwrap_or(usize::MAX)))
}

/// Defines a sparse integer vector type with index type `$idx`, the
/// `SparseIntVect` variant `$kind` in `wrapper/include/rdkit-cxx/data_structs.h`.
macro_rules! sparse_int_vect_type {
    ($(#[$meta:meta])* $name:ident, $idx:ty, $kind:expr) => {
        $(#[$meta])*
        ///
        /// Values are `i32`s. Arithmetic that would overflow one, or divide
        /// by zero, panics in the operators.
        pub struct $name {
            /// Never null.
            ptr: UniquePtr<ffi::SparseIntVect>,
        }

        impl sealed::SparseIntVectLike for $name {
            fn sparse_int_vect(&self) -> &ffi::SparseIntVect {
                self.raw()
            }
        }

        impl SparseIntVectType for $name {}

        impl $name {
            pub(crate) fn raw(&self) -> &ffi::SparseIntVect {
                self.ptr.as_ref().unwrap()
            }

            fn raw_mut(&mut self) -> Pin<&mut ffi::SparseIntVect> {
                self.ptr.pin_mut()
            }

            pub(crate) fn from_ffi(ptr: UniquePtr<ffi::SparseIntVect>) -> $name {
                assert!(!ptr.is_null(), "RDKit returned a null vector");
                $name { ptr }
            }

            /// A vector of `length` zeros; a negative length is
            /// `Error::Value`. Mirrors the class's `size` constructor.
            pub fn new(length: $idx) -> Result<$name> {
                let length = u64::try_from(length).map_err(|_| Error::Value {
                    message: format!("negative SparseIntVect length {length}"),
                })?;
                Ok($name::from_ffi(ffi::sparse_int_vect_new($kind, length)?))
            }

            /// A vector from a pickle (`to_binary`). Mirrors the class's
            /// `pkl` constructor.
            pub fn from_pkl(pkl: &[u8]) -> Result<$name> {
                Ok($name::from_ffi(ffi::sparse_int_vect_from_pkl($kind, pkl)?))
            }

            /// Mirrors `vect[idx]`.
            pub fn get_val(&self, idx: $idx) -> Result<i32> {
                let idx = index(idx, self.length_u64())?;
                Ok(ffi::sparse_int_vect_get_val(self.raw(), idx)?)
            }

            /// Mirrors `vect[idx] = val`.
            pub fn set_val(&mut self, idx: $idx, val: i32) -> Result<()> {
                let idx = index(idx, self.length_u64())?;
                Ok(ffi::sparse_int_vect_set_val(self.raw_mut(), idx, val)?)
            }

            fn length_u64(&self) -> u64 {
                expect_ok(ffi::sparse_int_vect_get_length(self.raw()))
            }

            /// Mirrors `GetLength`.
            pub fn get_length(&self) -> $idx {
                <$idx>::try_from(self.length_u64()).expect("the length fits its type")
            }

            /// The sum of the values. RDKit sums in `i32`, so values whose
            /// absolute values sum past `i32::MAX` are `Error::Value`.
            /// Mirrors `GetTotalVal`.
            pub fn get_total_val(&self) -> Result<i32> {
                self.get_total_val_with_use_abs(false)
            }

            /// The sum of the absolute values when `use_abs`. Mirrors
            /// `GetTotalVal(useAbs)`.
            pub fn get_total_val_with_use_abs(&self, use_abs: bool) -> Result<i32> {
                Ok(ffi::sparse_int_vect_get_total_val(self.raw(), use_abs)?)
            }

            /// RDKit's binary form, which `from_pkl` reads. Mirrors
            /// `ToBinary`.
            pub fn to_binary(&self) -> Vec<u8> {
                expect_ok(ffi::sparse_int_vect_to_binary(self.raw()))
            }

            /// Adds one at each index in `seq`; nothing changes if one is
            /// out of range. Mirrors `UpdateFromSequence`.
            pub fn update_from_sequence(&mut self, seq: &[$idx]) -> Result<()> {
                let len = self.length_u64();
                let seq = seq
                    .iter()
                    .map(|&idx| index(idx, len))
                    .collect::<Result<Vec<u64>>>()?;
                Ok(ffi::sparse_int_vect_update_from_sequence(
                    self.raw_mut(),
                    &seq,
                )?)
            }

            /// Mirrors `GetNonzeroElements`.
            pub fn get_nonzero_elements(&self) -> BTreeMap<$idx, i32> {
                let elements = expect_ok(ffi::sparse_int_vect_get_nonzero_elements(self.raw()));
                elements
                    .indices
                    .into_iter()
                    .map(|idx| <$idx>::try_from(idx).expect("the index fits its type"))
                    .zip(elements.values)
                    .collect()
            }

            /// Every value, zeros included. Mirrors `ToList`.
            pub fn to_list(&self) -> Vec<i32> {
                expect_ok(ffi::sparse_int_vect_to_list(self.raw()))
            }
        }

        /// A vector of length 0, a Rust addition (Python's constructors take a
        /// length).
        impl Default for $name {
            fn default() -> $name {
                $name::new(0).expect("0 is a length")
            }
        }

        impl Clone for $name {
            fn clone(&self) -> $name {
                $name::from_ffi(expect_ok(ffi::sparse_int_vect_clone(self.raw())))
            }
        }

        impl PartialEq for $name {
            fn eq(&self, other: &$name) -> bool {
                expect_ok(ffi::sparse_int_vect_eq(self.raw(), other.raw()))
            }
        }

        impl Eq for $name {}

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("length", &self.get_length())
                    .field("nonzero_elements", &self.get_nonzero_elements())
                    .finish()
            }
        }

        /// The smaller value at each index. Mirrors `&`.
        ///
        /// # Panics
        ///
        /// If the vectors differ in length, where Python raises `ValueError`.
        impl std::ops::BitAnd for &$name {
            type Output = $name;

            fn bitand(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::sparse_int_vect_binary_op(self.raw(), rhs.raw(), b'&')))
            }
        }

        /// The larger value at each index. Mirrors `|`.
        ///
        /// # Panics
        ///
        /// If the vectors differ in length, where Python raises `ValueError`.
        impl std::ops::BitOr for &$name {
            type Output = $name;

            fn bitor(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::sparse_int_vect_binary_op(self.raw(), rhs.raw(), b'|')))
            }
        }

        /// Mirrors `+`.
        ///
        /// # Panics
        ///
        /// If the vectors differ in length, or a value would overflow an `i32`.
        impl std::ops::Add for &$name {
            type Output = $name;

            fn add(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::sparse_int_vect_binary_op(self.raw(), rhs.raw(), b'+')))
            }
        }

        /// Mirrors `-`.
        ///
        /// # Panics
        ///
        /// If the vectors differ in length, or a value would overflow an `i32`.
        impl std::ops::Sub for &$name {
            type Output = $name;

            fn sub(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::sparse_int_vect_binary_op(self.raw(), rhs.raw(), b'-')))
            }
        }

        /// Mirrors `+=`.
        ///
        /// # Panics
        ///
        /// If the vectors differ in length, or a value would overflow an `i32`.
        impl std::ops::AddAssign<&$name> for $name {
            fn add_assign(&mut self, rhs: &$name) {
                op(ffi::sparse_int_vect_assign_op(self.raw_mut(), rhs.raw(), b'+'))
            }
        }

        /// Mirrors `-=`.
        ///
        /// # Panics
        ///
        /// If the vectors differ in length, or a value would overflow an `i32`.
        impl std::ops::SubAssign<&$name> for $name {
            fn sub_assign(&mut self, rhs: &$name) {
                op(ffi::sparse_int_vect_assign_op(self.raw_mut(), rhs.raw(), b'-'))
            }
        }

        /// Adds to each nonzero value, as RDKit does. Mirrors `+= int`.
        ///
        /// # Panics
        ///
        /// If a value would overflow an `i32`.
        impl std::ops::AddAssign<i32> for $name {
            fn add_assign(&mut self, rhs: i32) {
                op(ffi::sparse_int_vect_scalar_op(self.raw_mut(), b'+', rhs))
            }
        }

        /// Subtracts from each nonzero value. Mirrors `-= int`.
        ///
        /// # Panics
        ///
        /// If a value would overflow an `i32`.
        impl std::ops::SubAssign<i32> for $name {
            fn sub_assign(&mut self, rhs: i32) {
                op(ffi::sparse_int_vect_scalar_op(self.raw_mut(), b'-', rhs))
            }
        }

        /// Mirrors `*= int`.
        ///
        /// # Panics
        ///
        /// If a value would overflow an `i32`.
        impl std::ops::MulAssign<i32> for $name {
            fn mul_assign(&mut self, rhs: i32) {
                op(ffi::sparse_int_vect_scalar_op(self.raw_mut(), b'*', rhs))
            }
        }

        /// Integer division. Mirrors `/= int`.
        ///
        /// # Panics
        ///
        /// If `rhs` is 0, where Python raises, or a value would overflow an `i32` (`i32::MIN / -1`).
        impl std::ops::DivAssign<i32> for $name {
            fn div_assign(&mut self, rhs: i32) {
                op(ffi::sparse_int_vect_scalar_op(self.raw_mut(), b'/', rhs))
            }
        }
    };
}

sparse_int_vect_type!(
    /// A sparse vector of integers with `i32` indices. Mirrors
    /// `cDataStructs.IntSparseIntVect`.
    IntSparseIntVect,
    i32,
    0
);

sparse_int_vect_type!(
    /// A sparse vector of integers with `i64` indices. Mirrors
    /// `cDataStructs.LongSparseIntVect`.
    LongSparseIntVect,
    i64,
    1
);

sparse_int_vect_type!(
    /// A sparse vector of integers with `u32` indices. Mirrors
    /// `cDataStructs.UIntSparseIntVect`.
    UIntSparseIntVect,
    u32,
    2
);

sparse_int_vect_type!(
    /// A sparse vector of integers with `u64` indices. Mirrors
    /// `cDataStructs.ULongSparseIntVect`.
    ULongSparseIntVect,
    u64,
    3
);
