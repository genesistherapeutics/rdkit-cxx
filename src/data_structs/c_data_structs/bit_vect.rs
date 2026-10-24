// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::fmt;
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::error::{expect_ok, to_string};
use crate::ffi::data_structs::ffi;

pub(crate) mod sealed {
    use std::pin::Pin;

    use cxx::UniquePtr;

    use crate::ffi::data_structs::ffi;

    pub trait BitVectLike: Sized {
        fn bit_vect(&self) -> &ffi::BitVect;
        fn bit_vect_mut(&mut self) -> Pin<&mut ffi::BitVect>;
        fn from_ffi(ptr: UniquePtr<ffi::BitVect>) -> Self;
    }
}

/// A bit vector the module functions take: `ExplicitBitVect` or
/// `SparseBitVect`. Generic functions over it cover Python's overloads for
/// both (CONVENTIONS.md A5).
pub trait BitVect: sealed::BitVectLike {}

/// Unwraps an operator's result; operators panic where Python raises.
pub(crate) fn op<T>(result: std::result::Result<T, cxx::Exception>) -> T {
    result.unwrap_or_else(|err| panic!("{}", crate::Error::from(err)))
}

/// Defines a bit vector type over `ffi::BitVect` (W3).
macro_rules! bit_vect_type {
    ($(#[$meta:meta])* $name:ident, sparse = $sparse:expr) => {
        $(#[$meta])*
        pub struct $name {
            /// Never null.
            ptr: UniquePtr<ffi::BitVect>,
        }

        impl sealed::BitVectLike for $name {
            fn bit_vect(&self) -> &ffi::BitVect {
                self.ptr.as_ref().unwrap()
            }

            fn bit_vect_mut(&mut self) -> Pin<&mut ffi::BitVect> {
                self.ptr.pin_mut()
            }

            fn from_ffi(ptr: UniquePtr<ffi::BitVect>) -> $name {
                assert!(!ptr.is_null(), "RDKit returned a null bit vector");
                $name { ptr }
            }
        }

        impl BitVect for $name {}

        impl $name {
            pub(crate) fn raw(&self) -> &ffi::BitVect {
                sealed::BitVectLike::bit_vect(self)
            }

            fn raw_mut(&mut self) -> Pin<&mut ffi::BitVect> {
                sealed::BitVectLike::bit_vect_mut(self)
            }

            pub(crate) fn from_ffi(ptr: UniquePtr<ffi::BitVect>) -> $name {
                sealed::BitVectLike::from_ffi(ptr)
            }

            /// A vector from a pickle (`to_binary`). Mirrors the class's
            /// `pkl` constructor.
            pub fn from_pkl(pkl: &[u8]) -> Result<$name> {
                Ok($name::from_ffi(ffi::bit_vect_from_pkl($sparse, pkl)?))
            }

            /// Turns a bit on, returning its old value. Mirrors `SetBit`.
            pub fn set_bit(&mut self, which: u32) -> Result<bool> {
                Ok(ffi::bit_vect_set_bit(self.raw_mut(), which)?)
            }

            /// Turns bits on; none changes if any is out of range. Mirrors
            /// `SetBitsFromList`.
            pub fn set_bits_from_list(&mut self, on_bit_list: &[u32]) -> Result<()> {
                Ok(ffi::bit_vect_set_bits_from_list(self.raw_mut(), on_bit_list)?)
            }

            /// Turns a bit off, returning its old value. Mirrors `UnSetBit`.
            pub fn un_set_bit(&mut self, which: u32) -> Result<bool> {
                Ok(ffi::bit_vect_un_set_bit(self.raw_mut(), which)?)
            }

            /// Turns bits off; none changes if any is out of range. Mirrors
            /// `UnSetBitsFromList`.
            pub fn un_set_bits_from_list(&mut self, off_bit_list: &[u32]) -> Result<()> {
                Ok(ffi::bit_vect_un_set_bits_from_list(
                    self.raw_mut(),
                    off_bit_list,
                )?)
            }

            /// Mirrors `GetBit`.
            pub fn get_bit(&self, which: u32) -> Result<bool> {
                Ok(ffi::bit_vect_get_bit(self.raw(), which)?)
            }

            /// The vector's size. Mirrors `GetNumBits`.
            pub fn get_num_bits(&self) -> u32 {
                expect_ok(ffi::bit_vect_get_num_bits(self.raw()))
            }

            /// Mirrors `GetNumOnBits`.
            pub fn get_num_on_bits(&self) -> u32 {
                expect_ok(ffi::bit_vect_get_num_on_bits(self.raw()))
            }

            /// Mirrors `GetNumOffBits`.
            pub fn get_num_off_bits(&self) -> u32 {
                expect_ok(ffi::bit_vect_get_num_off_bits(self.raw()))
            }

            /// The on bits' indices, in order. Mirrors `GetOnBits`.
            pub fn get_on_bits(&self) -> Vec<u32> {
                expect_ok(ffi::bit_vect_get_on_bits(self.raw()))
            }

            /// RDKit's binary form, which `from_pkl` reads. Mirrors
            /// `ToBinary`.
            pub fn to_binary(&self) -> Vec<u8> {
                expect_ok(ffi::bit_vect_to_binary(self.raw()))
            }

            /// Replaces the vector with one from base64 text (`to_base64`).
            /// Mirrors `FromBase64`.
            pub fn from_base64(&mut self, in_d: &str) -> Result<()> {
                Ok(ffi::bit_vect_from_base64(self.raw_mut(), in_d)?)
            }

            /// The binary form in base64. Mirrors `ToBase64`.
            pub fn to_base64(&self) -> String {
                to_string(expect_ok(ffi::bit_vect_to_base64(self.raw())))
                    .expect("base64 is ASCII")
            }

            /// Every bit's value. Mirrors `ToList`.
            pub fn to_list(&self) -> Vec<bool> {
                expect_ok(ffi::bit_vect_to_list(self.raw()))
            }
        }

        impl Clone for $name {
            fn clone(&self) -> $name {
                $name::from_ffi(expect_ok(ffi::bit_vect_clone(self.raw())))
            }
        }

        impl PartialEq for $name {
            fn eq(&self, other: &$name) -> bool {
                expect_ok(ffi::bit_vect_eq(self.raw(), other.raw()))
            }
        }

        impl Eq for $name {}

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("num_bits", &self.get_num_bits())
                    .field("on_bits", &self.get_on_bits())
                    .finish()
            }
        }

        /// The bitwise AND of two vectors of one size. Mirrors `&`.
        impl std::ops::BitAnd for &$name {
            type Output = $name;

            fn bitand(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::bit_vect_binary_op(self.raw(), rhs.raw(), b'&')))
            }
        }

        /// The bitwise OR of two vectors of one size. Mirrors `|`.
        impl std::ops::BitOr for &$name {
            type Output = $name;

            fn bitor(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::bit_vect_binary_op(self.raw(), rhs.raw(), b'|')))
            }
        }

        /// The bitwise XOR of two vectors of one size. Mirrors `^`.
        impl std::ops::BitXor for &$name {
            type Output = $name;

            fn bitxor(self, rhs: &$name) -> $name {
                $name::from_ffi(op(ffi::bit_vect_binary_op(self.raw(), rhs.raw(), b'^')))
            }
        }

        /// The bitwise negation. Mirrors `~`.
        impl std::ops::Not for &$name {
            type Output = $name;

            fn not(self) -> $name {
                $name::from_ffi(expect_ok(ffi::bit_vect_invert(self.raw())))
            }
        }
    };
}

bit_vect_type!(
    /// A bit vector stored bit by bit, for sizes up to tens of thousands.
    /// Mirrors `cDataStructs.ExplicitBitVect`.
    ///
    /// `&a + &b` concatenates, as Python's `+` does.
    ExplicitBitVect,
    sparse = false
);

bit_vect_type!(
    /// A bit vector that stores its on bits, for large vectors with few of
    /// them. Mirrors `cDataStructs.SparseBitVect`.
    SparseBitVect,
    sparse = true
);

impl ExplicitBitVect {
    /// All bits off. Mirrors `ExplicitBitVect(size)`.
    pub fn new(size: u32) -> ExplicitBitVect {
        ExplicitBitVect::new_with_bits_set(size, false)
    }

    /// All bits on when `bits_set`. Mirrors `ExplicitBitVect(size, bitsSet)`.
    pub fn new_with_bits_set(size: u32, bits_set: bool) -> ExplicitBitVect {
        ExplicitBitVect::from_ffi(expect_ok(ffi::bit_vect_new_explicit(size, bits_set)))
    }
}

/// An empty vector. Mirrors `ExplicitBitVect()`.
impl Default for ExplicitBitVect {
    fn default() -> ExplicitBitVect {
        ExplicitBitVect::new(0)
    }
}

/// The concatenation. Mirrors `+`.
impl std::ops::Add for &ExplicitBitVect {
    type Output = ExplicitBitVect;

    fn add(self, rhs: &ExplicitBitVect) -> ExplicitBitVect {
        ExplicitBitVect::from_ffi(op(ffi::bit_vect_binary_op(self.raw(), rhs.raw(), b'+')))
    }
}

/// Appends `rhs`. Mirrors `+=`.
impl std::ops::AddAssign<&ExplicitBitVect> for ExplicitBitVect {
    fn add_assign(&mut self, rhs: &ExplicitBitVect) {
        op(ffi::bit_vect_append(self.raw_mut(), rhs.raw()))
    }
}

impl SparseBitVect {
    /// All bits off. Mirrors `SparseBitVect(size)`.
    pub fn new(size: u32) -> SparseBitVect {
        SparseBitVect::from_ffi(expect_ok(ffi::bit_vect_new_sparse(size)))
    }
}

/// An empty vector. Mirrors `SparseBitVect()`.
impl Default for SparseBitVect {
    fn default() -> SparseBitVect {
        SparseBitVect::new(0)
    }
}

/// The vector as 0s and 1s. Mirrors `cDataStructs.BitVectToText`.
pub fn bit_vect_to_text<B: BitVect>(bv1: &B) -> String {
    to_string(expect_ok(ffi::bit_vect_to_text(bv1.bit_vect()))).expect("the text is ASCII")
}

/// The vector in the FPS format's hex. Mirrors
/// `cDataStructs.BitVectToFPSText`.
pub fn bit_vect_to_fps_text<B: BitVect>(bv1: &B) -> String {
    to_string(expect_ok(ffi::bit_vect_to_fps_text(bv1.bit_vect()))).expect("the text is ASCII")
}

/// The vector's bits packed into bytes. Mirrors
/// `cDataStructs.BitVectToBinaryText`.
pub fn bit_vect_to_binary_text<B: BitVect>(bv: &B) -> Vec<u8> {
    expect_ok(ffi::bit_vect_to_binary_text(bv.bit_vect()))
}

/// Sets bits from Daylight's fingerprint encoding. Mirrors
/// `cDataStructs.InitFromDaylightString`.
pub fn init_from_daylight_string<B: BitVect>(sbv: &mut B, s: &str) -> Result<()> {
    Ok(ffi::bit_vect_init_from_daylight_string(
        sbv.bit_vect_mut(),
        s,
    )?)
}

/// Mirrors `cDataStructs.ConvertToExplicit`.
pub fn convert_to_explicit(sbv: &SparseBitVect) -> ExplicitBitVect {
    ExplicitBitVect::from_ffi(expect_ok(ffi::convert_to_explicit(sbv.raw())))
}

/// A vector with a bit for each character, on for `1`. Mirrors
/// `cDataStructs.CreateFromBitString`.
pub fn create_from_bit_string(bits: &str) -> ExplicitBitVect {
    ExplicitBitVect::from_ffi(expect_ok(ffi::create_from_bit_string(bits)))
}

/// A vector from FPS hex, which must have an even length. Mirrors
/// `cDataStructs.CreateFromFPSText`.
pub fn create_from_fps_text(fps: &str) -> Result<ExplicitBitVect> {
    Ok(ExplicitBitVect::from_ffi(ffi::create_from_fps_text(fps)?))
}

/// A vector from packed bytes (`bit_vect_to_binary_text`). Mirrors
/// `cDataStructs.CreateFromBinaryText`.
pub fn create_from_binary_text(fps: &[u8]) -> ExplicitBitVect {
    ExplicitBitVect::from_ffi(expect_ok(ffi::create_from_binary_text(fps)))
}
