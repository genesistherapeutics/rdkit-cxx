// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/metric_matrix_calc.h");

        type BitVect = crate::ffi::data_structs::ffi::BitVect;
        type BitVectList;

        fn bit_vect_list_new() -> Result<UniquePtr<BitVectList>>;
        fn bit_vect_list_push(list: Pin<&mut BitVectList>, bv: &BitVect) -> Result<()>;

        fn get_euclidean_dist_mat(values: &[f64], nrows: u32, ncols: u32) -> Result<Vec<f64>>;
        fn get_tanimoto_dist_mat(list: &BitVectList) -> Result<Vec<f64>>;
        fn get_tanimoto_sim_mat(list: &BitVectList) -> Result<Vec<f64>>;
    }
}
