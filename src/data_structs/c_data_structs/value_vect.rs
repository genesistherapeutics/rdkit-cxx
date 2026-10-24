// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::data_structs::c_data_structs::{
    ExplicitBitVect, IntSparseIntVect, LongSparseIntVect, UIntSparseIntVect, ULongSparseIntVect,
};
use crate::error::{expect_ok, op};
use crate::ffi::data_structs::ffi;

rdkit_enum! {
    /// The bits per value of a `DiscreteValueVect`. Mirrors
    /// `cDataStructs.DiscreteValueType`.
    DiscreteValueType: ffi::DiscreteValueType {
        Onebitvalue = ONEBITVALUE,
        Twobitvalue = TWOBITVALUE,
        Fourbitvalue = FOURBITVALUE,
        Eightbitvalue = EIGHTBITVALUE,
        Sixteenbitvalue = SIXTEENBITVALUE,
    }
}

/// Defines the methods a value vector shares, over the `$prefix` wrappers
/// (W3).
macro_rules! value_vect_type {
    (
        $(#[$meta:meta])* $name:ident, $ffi:ident, $value:ty,
        $clone:ident, $from_pkl:ident, $get_val:ident, $set_val:ident, $get_length:ident,
        $get_total_val:ident, $to_binary:ident, $binary_op:ident, $assign_op:ident
    ) => {
        $(#[$meta])*
        pub struct $name {
            /// Never null.
            ptr: UniquePtr<ffi::$ffi>,
        }

        impl $name {
            fn raw(&self) -> &ffi::$ffi {
                self.ptr.as_ref().unwrap()
            }

            fn raw_mut(&mut self) -> Pin<&mut ffi::$ffi> {
                self.ptr.pin_mut()
            }

            pub(crate) fn from_ffi(ptr: UniquePtr<ffi::$ffi>) -> $name {
                assert!(!ptr.is_null(), "RDKit returned a null vector");
                $name { ptr }
            }

            /// A vector from a pickle (`to_binary`). Mirrors the class's
            /// `pkl` constructor.
            pub fn from_pkl(pkl: &[u8]) -> Result<$name> {
                Ok($name::from_ffi(ffi::$from_pkl(pkl)?))
            }

            /// Mirrors `vect[i]`.
            pub fn get_val(&self, i: u32) -> Result<$value> {
                Ok(ffi::$get_val(self.raw(), i)?)
            }

            /// Mirrors `vect[i] = val`.
            pub fn set_val(&mut self, i: u32, val: $value) -> Result<()> {
                Ok(ffi::$set_val(self.raw_mut(), i, val)?)
            }

            /// The number of values. Mirrors `len(vect)`.
            pub fn len(&self) -> u32 {
                expect_ok(ffi::$get_length(self.raw()))
            }

            /// Whether the vector has no values.
            pub fn is_empty(&self) -> bool {
                self.len() == 0
            }

            /// The sum of the values. Mirrors `GetTotalVal`.
            pub fn get_total_val(&self) -> $value {
                expect_ok(ffi::$get_total_val(self.raw()))
            }

            /// RDKit's binary form, which `from_pkl` reads; Python pickles
            /// it.
            pub fn to_binary(&self) -> Vec<u8> {
                expect_ok(ffi::$to_binary(self.raw()))
            }
        }

        impl Clone for $name {
            fn clone(&self) -> $name {
                $name::from_ffi(expect_ok(ffi::$clone(self.raw())))
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("len", &self.len())
                    .finish()
            }
        }

        /// The smaller value at each index. Mirrors `&`.
        impl std::ops::BitAnd for &$name {
            type Output = $name;

            fn bitand(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::$binary_op(self.raw(), rhs.raw(), b'&')))
            }
        }

        /// The larger value at each index. Mirrors `|`.
        impl std::ops::BitOr for &$name {
            type Output = $name;

            fn bitor(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::$binary_op(self.raw(), rhs.raw(), b'|')))
            }
        }

        /// Mirrors `+`.
        impl std::ops::Add for &$name {
            type Output = $name;

            fn add(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::$binary_op(self.raw(), rhs.raw(), b'+')))
            }
        }

        /// Mirrors `-`.
        impl std::ops::Sub for &$name {
            type Output = $name;

            fn sub(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::$binary_op(self.raw(), rhs.raw(), b'-')))
            }
        }

        /// Mirrors `+=`.
        impl std::ops::AddAssign<&$name> for $name {
            fn add_assign(&mut self, rhs: &$name) {
                op(ffi::$assign_op(self.raw_mut(), rhs.raw(), b'+'))
            }
        }

        /// Mirrors `-=`.
        impl std::ops::SubAssign<&$name> for $name {
            fn sub_assign(&mut self, rhs: &$name) {
                op(ffi::$assign_op(self.raw_mut(), rhs.raw(), b'-'))
            }
        }
    };
}

value_vect_type!(
    /// A vector of small unsigned integers, packed by their bits per value.
    /// Mirrors `cDataStructs.DiscreteValueVect`.
    ///
    /// Sums saturate, and differences stop at zero.
    DiscreteValueVect, DiscreteValueVect, u32,
    discrete_value_vect_clone, discrete_value_vect_from_pkl, discrete_value_vect_get_val,
    discrete_value_vect_set_val, discrete_value_vect_get_length,
    discrete_value_vect_get_total_val, discrete_value_vect_to_binary,
    discrete_value_vect_binary_op, discrete_value_vect_assign_op
);

value_vect_type!(
    /// A vector of real numbers. Mirrors `cDataStructs.RealValueVect`.
    RealValueVect, RealValueVect, f64,
    real_value_vect_clone, real_value_vect_from_pkl, real_value_vect_get_val,
    real_value_vect_set_val, real_value_vect_get_length, real_value_vect_get_total_val,
    real_value_vect_to_binary, real_value_vect_binary_op, real_value_vect_assign_op
);

impl DiscreteValueVect {
    /// `length` zeros of `val_type`. Mirrors `DiscreteValueVect(valType,
    /// length)`.
    pub fn new(val_type: DiscreteValueType, length: u32) -> DiscreteValueVect {
        DiscreteValueVect::from_ffi(expect_ok(ffi::discrete_value_vect_new(
            val_type.to_ffi(),
            length,
        )))
    }

    /// Mirrors `GetValueType`.
    pub fn get_value_type(&self) -> DiscreteValueType {
        DiscreteValueType::from_ffi(expect_ok(ffi::discrete_value_vect_get_value_type(
            self.raw(),
        )))
    }
}

/// An empty vector of one-bit values. Mirrors `DiscreteValueVect()`.
impl Default for DiscreteValueVect {
    fn default() -> DiscreteValueVect {
        DiscreteValueVect::new(DiscreteValueType::Onebitvalue, 0)
    }
}

impl PartialEq for DiscreteValueVect {
    fn eq(&self, other: &DiscreteValueVect) -> bool {
        self.to_binary() == other.to_binary()
    }
}

impl RealValueVect {
    /// `length` zeros. Mirrors `RealValueVect(length)`.
    pub fn new(length: u32) -> RealValueVect {
        RealValueVect::from_ffi(expect_ok(ffi::real_value_vect_new(length)))
    }
}

/// An empty vector. Mirrors `RealValueVect()`.
impl Default for RealValueVect {
    fn default() -> RealValueVect {
        RealValueVect::new(0)
    }
}

/// The sum of the absolute differences, for vectors of one length and value
/// type. Mirrors `cDataStructs.ComputeL1Norm`.
pub fn compute_l1_norm(v1: &DiscreteValueVect, v2: &DiscreteValueVect) -> Result<u32> {
    Ok(ffi::discrete_value_vect_compute_l1_norm(
        v1.raw(),
        v2.raw(),
    )?)
}

/// As [`compute_l1_norm`], for real vectors. Mirrors
/// `cDataStructs.ComputeL1Norm(RealValueVect, RealValueVect)` (A5).
pub fn compute_l1_norm_real_value_vect(v1: &RealValueVect, v2: &RealValueVect) -> Result<f64> {
    Ok(ffi::real_value_vect_compute_l1_norm(v1.raw(), v2.raw())?)
}

/// A vector `convert_to_numpy_array` takes, with its values' type.
pub trait ConvertToNumpyArray {
    type Value;

    #[doc(hidden)]
    fn __values(&self) -> Vec<Self::Value>;
}

impl ConvertToNumpyArray for DiscreteValueVect {
    type Value = u32;

    fn __values(&self) -> Vec<u32> {
        (0..self.len())
            .map(|i| self.get_val(i).expect("the index is in range"))
            .collect()
    }
}

impl ConvertToNumpyArray for RealValueVect {
    type Value = f64;

    fn __values(&self) -> Vec<f64> {
        (0..self.len())
            .map(|i| self.get_val(i).expect("the index is in range"))
            .collect()
    }
}

impl ConvertToNumpyArray for ExplicitBitVect {
    type Value = i32;

    fn __values(&self) -> Vec<i32> {
        self.to_list().into_iter().map(i32::from).collect()
    }
}

macro_rules! sparse_int_vect_numpy {
    ($($name:ident),*) => {
        $(
            impl ConvertToNumpyArray for $name {
                type Value = i32;

                fn __values(&self) -> Vec<i32> {
                    self.to_list()
                }
            }
        )*
    };
}

sparse_int_vect_numpy!(
    IntSparseIntVect,
    LongSparseIntVect,
    UIntSparseIntVect,
    ULongSparseIntVect
);

/// Every value of the vector, as Python writes them into `destArray`.
/// Mirrors `cDataStructs.ConvertToNumpyArray`.
pub fn convert_to_numpy_array<V: ConvertToNumpyArray>(v: &V) -> Vec<V::Value> {
    v.__values()
}
