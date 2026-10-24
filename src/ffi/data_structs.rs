// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::DiscreteValueVect::DiscreteValueType`; cxx checks the values.
    #[repr(u32)]
    #[allow(clippy::upper_case_acronyms)]
    enum DiscreteValueType {
        ONEBITVALUE = 0,
        TWOBITVALUE,
        FOURBITVALUE,
        EIGHTBITVALUE,
        SIXTEENBITVALUE,
    }

    /// A sparse vector's nonzero entries.
    struct SparseIntVectElements {
        indices: Vec<u64>,
        values: Vec<i32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/data_structs.h");

        type BitVect;
        type SparseIntVect;
        type DiscreteValueVect;
        type DiscreteValueType;
        type RealValueVect;

        fn bit_vect_new_explicit(size: u32, bits_set: bool) -> Result<UniquePtr<BitVect>>;
        fn bit_vect_new_sparse(size: u32) -> Result<UniquePtr<BitVect>>;
        fn bit_vect_from_pkl(sparse: bool, pkl: &[u8]) -> Result<UniquePtr<BitVect>>;
        fn bit_vect_clone(bv: &BitVect) -> Result<UniquePtr<BitVect>>;
        fn bit_vect_set_bit(bv: Pin<&mut BitVect>, which: u32) -> Result<bool>;
        fn bit_vect_un_set_bit(bv: Pin<&mut BitVect>, which: u32) -> Result<bool>;
        fn bit_vect_set_bits_from_list(bv: Pin<&mut BitVect>, on_bit_list: &[u32]) -> Result<()>;
        fn bit_vect_un_set_bits_from_list(
            bv: Pin<&mut BitVect>,
            off_bit_list: &[u32],
        ) -> Result<()>;
        fn bit_vect_get_bit(bv: &BitVect, which: u32) -> Result<bool>;
        fn bit_vect_get_num_bits(bv: &BitVect) -> Result<u32>;
        fn bit_vect_get_num_on_bits(bv: &BitVect) -> Result<u32>;
        fn bit_vect_get_num_off_bits(bv: &BitVect) -> Result<u32>;
        fn bit_vect_get_on_bits(bv: &BitVect) -> Result<Vec<u32>>;
        fn bit_vect_to_binary(bv: &BitVect) -> Result<Vec<u8>>;
        fn bit_vect_from_base64(bv: Pin<&mut BitVect>, in_d: &str) -> Result<()>;
        fn bit_vect_to_base64(bv: &BitVect) -> Result<UniquePtr<CxxString>>;
        fn bit_vect_to_list(bv: &BitVect) -> Result<Vec<bool>>;
        fn bit_vect_binary_op(a: &BitVect, b: &BitVect, op: u8) -> Result<UniquePtr<BitVect>>;
        fn bit_vect_append(bv: Pin<&mut BitVect>, other: &BitVect) -> Result<()>;
        fn bit_vect_invert(bv: &BitVect) -> Result<UniquePtr<BitVect>>;
        fn bit_vect_eq(a: &BitVect, b: &BitVect) -> Result<bool>;
        fn bit_vect_similarity(
            a: &BitVect,
            b: &BitVect,
            metric: u8,
            alpha: f64,
            beta: f64,
            return_distance: bool,
        ) -> Result<f64>;
        fn bit_vect_similarity_pkl(
            a: &BitVect,
            pkl: &[u8],
            metric: u8,
            alpha: f64,
            beta: f64,
            return_distance: bool,
        ) -> Result<f64>;
        fn bit_vect_proj_similarity(a: &BitVect, b: &BitVect, on_bits: bool) -> Result<Vec<f64>>;
        fn bit_vect_num_bits_in_common(a: &BitVect, b: &BitVect) -> Result<i32>;
        fn bit_vect_bits_in_common(a: &BitVect, b: &BitVect, on_bits: bool) -> Result<Vec<u32>>;
        fn bit_vect_fold(bv: &BitVect, fold_factor: u32) -> Result<UniquePtr<BitVect>>;
        fn bit_vect_all_probe_bits_match(probe: &BitVect, reference: &BitVect) -> Result<bool>;
        fn bit_vect_all_probe_bits_match_pkl(probe: &BitVect, reference: &[u8]) -> Result<bool>;
        fn bit_vect_to_text(bv: &BitVect) -> Result<UniquePtr<CxxString>>;
        fn bit_vect_to_fps_text(bv: &BitVect) -> Result<UniquePtr<CxxString>>;
        fn bit_vect_to_binary_text(bv: &BitVect) -> Result<Vec<u8>>;
        fn bit_vect_init_from_daylight_string(bv: Pin<&mut BitVect>, s: &str) -> Result<()>;
        fn convert_to_explicit(sbv: &BitVect) -> Result<UniquePtr<BitVect>>;
        fn create_from_bit_string(bits: &str) -> Result<UniquePtr<BitVect>>;
        fn create_from_fps_text(fps: &str) -> Result<UniquePtr<BitVect>>;
        fn create_from_binary_text(fps: &[u8]) -> Result<UniquePtr<BitVect>>;

        fn sparse_int_vect_new(kind: u8, length: u64) -> Result<UniquePtr<SparseIntVect>>;
        fn sparse_int_vect_from_pkl(kind: u8, pkl: &[u8]) -> Result<UniquePtr<SparseIntVect>>;
        fn sparse_int_vect_clone(siv: &SparseIntVect) -> Result<UniquePtr<SparseIntVect>>;
        fn sparse_int_vect_get_val(siv: &SparseIntVect, idx: u64) -> Result<i32>;
        fn sparse_int_vect_set_val(siv: Pin<&mut SparseIntVect>, idx: u64, val: i32) -> Result<()>;
        fn sparse_int_vect_get_length(siv: &SparseIntVect) -> Result<u64>;
        fn sparse_int_vect_get_total_val(siv: &SparseIntVect, use_abs: bool) -> Result<i32>;
        fn sparse_int_vect_to_binary(siv: &SparseIntVect) -> Result<Vec<u8>>;
        fn sparse_int_vect_update_from_sequence(
            siv: Pin<&mut SparseIntVect>,
            seq: &[u64],
        ) -> Result<()>;
        fn sparse_int_vect_get_nonzero_elements(
            siv: &SparseIntVect,
        ) -> Result<SparseIntVectElements>;
        fn sparse_int_vect_to_list(siv: &SparseIntVect) -> Result<Vec<i32>>;
        fn sparse_int_vect_binary_op(
            a: &SparseIntVect,
            b: &SparseIntVect,
            op: u8,
        ) -> Result<UniquePtr<SparseIntVect>>;
        fn sparse_int_vect_assign_op(
            a: Pin<&mut SparseIntVect>,
            b: &SparseIntVect,
            op: u8,
        ) -> Result<()>;
        fn sparse_int_vect_scalar_op(siv: Pin<&mut SparseIntVect>, op: u8, val: i32) -> Result<()>;
        fn sparse_int_vect_eq(a: &SparseIntVect, b: &SparseIntVect) -> Result<bool>;
        #[allow(clippy::too_many_arguments)]
        fn sparse_int_vect_similarity(
            a: &SparseIntVect,
            b: &SparseIntVect,
            metric: u8,
            alpha: f64,
            beta: f64,
            return_distance: bool,
            bounds: f64,
        ) -> Result<f64>;

        fn discrete_value_vect_new(
            val_type: DiscreteValueType,
            length: u32,
        ) -> Result<UniquePtr<DiscreteValueVect>>;
        fn discrete_value_vect_from_pkl(pkl: &[u8]) -> Result<UniquePtr<DiscreteValueVect>>;
        fn discrete_value_vect_clone(
            dvv: &DiscreteValueVect,
        ) -> Result<UniquePtr<DiscreteValueVect>>;
        fn discrete_value_vect_get_val(dvv: &DiscreteValueVect, i: u32) -> Result<u32>;
        fn discrete_value_vect_set_val(
            dvv: Pin<&mut DiscreteValueVect>,
            i: u32,
            val: u32,
        ) -> Result<()>;
        fn discrete_value_vect_get_length(dvv: &DiscreteValueVect) -> Result<u32>;
        fn discrete_value_vect_get_value_type(dvv: &DiscreteValueVect)
        -> Result<DiscreteValueType>;
        fn discrete_value_vect_get_total_val(dvv: &DiscreteValueVect) -> Result<u32>;
        fn discrete_value_vect_to_binary(dvv: &DiscreteValueVect) -> Result<Vec<u8>>;
        fn discrete_value_vect_binary_op(
            a: &DiscreteValueVect,
            b: &DiscreteValueVect,
            op: u8,
        ) -> Result<UniquePtr<DiscreteValueVect>>;
        fn discrete_value_vect_assign_op(
            a: Pin<&mut DiscreteValueVect>,
            b: &DiscreteValueVect,
            op: u8,
        ) -> Result<()>;
        fn discrete_value_vect_compute_l1_norm(
            a: &DiscreteValueVect,
            b: &DiscreteValueVect,
        ) -> Result<u32>;

        fn real_value_vect_new(length: u32) -> Result<UniquePtr<RealValueVect>>;
        fn real_value_vect_from_pkl(pkl: &[u8]) -> Result<UniquePtr<RealValueVect>>;
        fn real_value_vect_clone(rvv: &RealValueVect) -> Result<UniquePtr<RealValueVect>>;
        fn real_value_vect_get_val(rvv: &RealValueVect, i: u32) -> Result<f64>;
        fn real_value_vect_set_val(rvv: Pin<&mut RealValueVect>, i: u32, val: f64) -> Result<()>;
        fn real_value_vect_get_length(rvv: &RealValueVect) -> Result<u32>;
        fn real_value_vect_get_total_val(rvv: &RealValueVect) -> Result<f64>;
        fn real_value_vect_to_binary(rvv: &RealValueVect) -> Result<Vec<u8>>;
        fn real_value_vect_binary_op(
            a: &RealValueVect,
            b: &RealValueVect,
            op: u8,
        ) -> Result<UniquePtr<RealValueVect>>;
        fn real_value_vect_assign_op(
            a: Pin<&mut RealValueVect>,
            b: &RealValueVect,
            op: u8,
        ) -> Result<()>;
        fn real_value_vect_compute_l1_norm(a: &RealValueVect, b: &RealValueVect) -> Result<f64>;
    }
}

// SAFETY: the vectors are plain data with no thread affinity or caches their
// const methods write.
unsafe impl Send for ffi::BitVect {}
unsafe impl Sync for ffi::BitVect {}
unsafe impl Send for ffi::SparseIntVect {}
unsafe impl Sync for ffi::SparseIntVect {}
unsafe impl Send for ffi::DiscreteValueVect {}
unsafe impl Sync for ffi::DiscreteValueVect {}
unsafe impl Send for ffi::RealValueVect {}
unsafe impl Sync for ffi::RealValueVect {}
