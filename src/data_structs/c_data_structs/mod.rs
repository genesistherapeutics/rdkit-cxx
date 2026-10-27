// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.DataStructs.cDataStructs`: bit vectors, sparse integer
//! vectors, value vectors, their similarity metrics, and FPB files.
//!
//! The binary operators (`&a & &b` and the like) panic where Python raises:
//! on vectors of different sizes, and on integer overflow or division by
//! zero.

mod bit_vect;
mod fpb;
mod similarity;
mod sparse_int_vect;
mod value_vect;

pub use bit_vect::{
    BitVect, ExplicitBitVect, SparseBitVect, bit_vect_to_binary_text, bit_vect_to_fps_text,
    bit_vect_to_text, convert_to_explicit, create_from_binary_text, create_from_bit_string,
    create_from_fps_text, init_from_daylight_string,
};
pub use fpb::{FPBReader, MultiFPBReader, MultiFPBReaderNeighborsParams};
pub use similarity::*;
pub use sparse_int_vect::{
    IntSparseIntVect, LongSparseIntVect, SparseIntVectType, UIntSparseIntVect, ULongSparseIntVect,
};
pub use value_vect::{
    ConvertToNumpyArray, DiscreteValueType, DiscreteValueVect, compute_l1_norm, convert_to_numpy_array,
};
