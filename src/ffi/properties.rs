// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::panic::{AssertUnwindSafe, catch_unwind};

use cxx::UniquePtr;

use crate::chem::Mol;
use crate::ffi::mol::ffi::RWMol;

/// The closure behind a property functor made in Rust. RDKit's registry is
/// global, so any thread may call it.
pub(crate) struct PropertyCallback(pub(crate) Box<dyn Fn(&Mol) -> f64 + Send + Sync>);

/// Called by the C++ functor with its own copy of the molecule. A panic is
/// caught here, at the boundary, and becomes the exception C++ sees.
fn call_property_callback(
    callback: &PropertyCallback,
    mol: UniquePtr<RWMol>,
) -> Result<f64, String> {
    let mol = Mol::from_ffi(mol);
    catch_unwind(AssertUnwindSafe(|| (callback.0)(&mol))).map_err(|payload| {
        let message = payload
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| payload.downcast_ref::<String>().cloned())
            .unwrap_or_else(|| "a property functor panicked".to_string());
        format!("property functor panicked: {message}")
    })
}

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    extern "Rust" {
        type PropertyCallback;

        fn call_property_callback(
            callback: &PropertyCallback,
            mol: UniquePtr<RWMol>,
        ) -> Result<f64>;
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/properties.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;
        type PropertyFunctorHandle;
        type Properties;
        type PropertyRangeQuery;

        fn property_functor_new(
            name: &str,
            version: &str,
            callback: Box<PropertyCallback>,
        ) -> Result<UniquePtr<PropertyFunctorHandle>>;
        fn property_functor_get_name(
            functor: &PropertyFunctorHandle,
        ) -> Result<UniquePtr<CxxString>>;
        fn property_functor_get_version(
            functor: &PropertyFunctorHandle,
        ) -> Result<UniquePtr<CxxString>>;
        fn property_functor_call(functor: &PropertyFunctorHandle, mol: &RWMol) -> Result<f64>;

        fn properties_new() -> Result<UniquePtr<Properties>>;
        fn properties_from_prop_names(prop_names: &[&str]) -> Result<UniquePtr<Properties>>;
        fn properties_get_property_names(
            properties: &Properties,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn properties_compute_properties(properties: &Properties, mol: &RWMol) -> Result<Vec<f64>>;
        fn properties_compute_properties_annotate(
            properties: &Properties,
            mol: Pin<&mut RWMol>,
        ) -> Result<Vec<f64>>;
        fn properties_annotate_properties(
            properties: &Properties,
            mol: Pin<&mut RWMol>,
        ) -> Result<()>;
        fn properties_get_available_properties() -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn properties_get_property(name: &str) -> Result<UniquePtr<PropertyFunctorHandle>>;
        fn properties_register_property(functor: &PropertyFunctorHandle) -> Result<i32>;

        fn make_property_range_query(
            name: &str,
            min: f64,
            max: f64,
        ) -> Result<UniquePtr<PropertyRangeQuery>>;
        fn property_range_query_match(query: &PropertyRangeQuery, mol: &RWMol) -> Result<bool>;
    }
}
