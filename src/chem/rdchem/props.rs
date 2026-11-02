// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Typed properties, shared by every class that has them (CONVENTIONS.md
//! P1).

use std::collections::BTreeMap;

use crate::error::to_string;
use crate::ffi::props::ffi;
use crate::{Error, Result};

/// A property value, as `GetPropsAsDict` and `GetProp(autoConvert=True)`
/// convert it. A Rust addition, for the Python values they return.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum PropValue {
    /// An `int`.
    Int(i32),
    /// An `unsigned int`.
    Unsigned(u32),
    /// Also holds C++ `float` values.
    Double(f64),
    /// A `bool`.
    Bool(bool),
    /// A string.
    String(String),
    /// A vector of `int`s.
    IntVec(Vec<i32>),
    /// A vector of `unsigned int`s.
    UnsignedVec(Vec<u32>),
    /// Also holds C++ `float` vectors.
    DoubleVec(Vec<f64>),
    /// A vector of strings.
    StringVec(Vec<String>),
    /// An empty value. Python's `None`.
    None,
}

impl PropValue {
    fn from_ffi(data: ffi::PropValueData) -> Result<PropValue> {
        Ok(match data.tag {
            0 => PropValue::None,
            1 => PropValue::Int(data.int_value),
            2 => PropValue::Unsigned(data.unsigned_value),
            3 => PropValue::Double(data.double_value),
            4 => PropValue::Bool(data.bool_value),
            5 => PropValue::String(to_string(data.string_value)?),
            6 => PropValue::IntVec(data.int_vec),
            7 => PropValue::UnsignedVec(data.unsigned_vec),
            8 => PropValue::DoubleVec(data.double_vec),
            9 => PropValue::StringVec(strings(data.string_vec)?),
            tag => panic!("unknown property tag {tag}"),
        })
    }
}

/// Parameters for `get_prop_names_with_params`, with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GetPropNamesParams {
    /// Includes names starting with `_`.
    pub include_private: bool,
    /// Includes properties RDKit computed.
    pub include_computed: bool,
}

/// Parameters for `get_props_as_dict_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetPropsAsDictParams {
    /// Includes names starting with `_`.
    pub include_private: bool,
    /// Includes properties RDKit computed.
    pub include_computed: bool,
    /// Converts strings that parse as integers or doubles.
    pub auto_convert_strings: bool,
}

impl Default for GetPropsAsDictParams {
    fn default() -> GetPropsAsDictParams {
        GetPropsAsDictParams {
            include_private: false,
            include_computed: false,
            auto_convert_strings: true,
        }
    }
}

pub(crate) fn strings(list: cxx::UniquePtr<cxx::CxxVector<cxx::CxxString>>) -> Result<Vec<String>> {
    let list = list.as_ref().expect("RDKit returned a null list");
    list.iter()
        .map(|text| {
            text.to_str()
                .map(str::to_owned)
                .map_err(|err| Error::Value {
                    message: format!("invalid UTF-8 in a string from RDKit: {err}"),
                })
        })
        .collect()
}

pub(crate) fn get_auto_convert(obj: &ffi::RDProps, key: &str) -> Result<PropValue> {
    PropValue::from_ffi(ffi::props_get_prop_auto_convert(obj, key)?)
}

pub(crate) fn get_names(obj: &ffi::RDProps, params: &GetPropNamesParams) -> Result<Vec<String>> {
    strings(ffi::props_get_prop_names(
        obj,
        params.include_private,
        params.include_computed,
    )?)
}

pub(crate) fn get_dict(
    obj: &ffi::RDProps,
    params: &GetPropsAsDictParams,
) -> Result<BTreeMap<String, PropValue>> {
    ffi::props_get_props_as_dict(
        obj,
        params.include_private,
        params.include_computed,
        params.auto_convert_strings,
    )?
    .into_iter()
    .map(|mut data| {
        let key = to_string(std::mem::replace(&mut data.key, cxx::UniquePtr::null()))?;
        Ok((key, PropValue::from_ffi(data)?))
    })
    .collect()
}

/// Stamps the property getters onto a type with
/// `fn props(&self) -> &ffi::RDProps` (W3, P1). The `dict_defaults` form
/// sets the parameters `get_props_as_dict` uses, for classes whose Python
/// `GetPropsAsDict` has other defaults.
macro_rules! prop_getters {
    ($type:ident $(<$lifetime:lifetime>)?) => {
        $crate::chem::rdchem::props::prop_getters!(
            $type $(<$lifetime>)?,
            dict_defaults = $crate::chem::rdchem::GetPropsAsDictParams::default(),
            dict_doc = "Public, non-computed properties"
        );
    };
    (
        $type:ident $(<$lifetime:lifetime>)?,
        dict_defaults = $dict_defaults:expr,
        dict_doc = $dict_doc:literal
    ) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Mirrors `HasProp`.
            pub fn has_prop(&self, key: &str) -> bool {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_has_prop(self.props(), key))
            }

            /// The property as a string, converting other types. A missing key
            /// is `Error::Key`. Mirrors `GetProp`.
            pub fn get_prop(&self, key: &str) -> $crate::Result<String> {
                $crate::error::to_string($crate::ffi::props::ffi::props_get_prop(self.props(), key)?)
            }

            /// The property, with strings that parse as numbers converted.
            /// Mirrors `GetProp(key, autoConvert=True)`.
            pub fn get_prop_auto_convert(
                &self,
                key: &str,
            ) -> $crate::Result<$crate::chem::rdchem::PropValue> {
                $crate::chem::rdchem::props::get_auto_convert(self.props(), key)
            }

            /// Mirrors `GetIntProp`.
            pub fn get_int_prop(&self, key: &str) -> $crate::Result<i32> {
                Ok($crate::ffi::props::ffi::props_get_int_prop(self.props(), key)?)
            }

            /// Mirrors `GetUnsignedProp`.
            pub fn get_unsigned_prop(&self, key: &str) -> $crate::Result<u32> {
                Ok($crate::ffi::props::ffi::props_get_unsigned_prop(self.props(), key)?)
            }

            /// Mirrors `GetDoubleProp`.
            pub fn get_double_prop(&self, key: &str) -> $crate::Result<f64> {
                Ok($crate::ffi::props::ffi::props_get_double_prop(self.props(), key)?)
            }

            /// Mirrors `GetBoolProp`.
            pub fn get_bool_prop(&self, key: &str) -> $crate::Result<bool> {
                Ok($crate::ffi::props::ffi::props_get_bool_prop(self.props(), key)?)
            }

            /// Public, non-computed property names. Mirrors `GetPropNames`.
            pub fn get_prop_names(&self) -> $crate::Result<Vec<String>> {
                self.get_prop_names_with_params(&Default::default())
            }

            /// Mirrors `GetPropNames(includePrivate, includeComputed)`.
            pub fn get_prop_names_with_params(
                &self,
                params: &$crate::chem::rdchem::GetPropNamesParams,
            ) -> $crate::Result<Vec<String>> {
                $crate::chem::rdchem::props::get_names(self.props(), params)
            }

            #[doc = concat!(
                $dict_doc,
                ", with strings that parse as numbers converted. Mirrors `GetPropsAsDict`."
            )]
            pub fn get_props_as_dict(
                &self,
            ) -> $crate::Result<
                std::collections::BTreeMap<String, $crate::chem::rdchem::PropValue>,
            > {
                self.get_props_as_dict_with_params(&$dict_defaults)
            }

            /// Mirrors `GetPropsAsDict(includePrivate, includeComputed,
            /// autoConvertStrings)`.
            pub fn get_props_as_dict_with_params(
                &self,
                params: &$crate::chem::rdchem::GetPropsAsDictParams,
            ) -> $crate::Result<
                std::collections::BTreeMap<String, $crate::chem::rdchem::PropValue>,
            > {
                $crate::chem::rdchem::props::get_dict(self.props(), params)
            }
        }
    };
}

/// Stamps the property setters onto a type with
/// `fn props(&self) -> &ffi::RDProps`. The `clear_computed` form adds
/// `clear_computed_props`, for classes whose Python counterpart has it (`Mol`
/// defines its own, with `include_rings`) (W3, P1).
macro_rules! prop_setters {
    ($type:ident $(<$lifetime:lifetime>)?, clear_computed) => {
        $crate::chem::rdchem::props::prop_setters!($type $(<$lifetime>)?);

        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Removes the properties RDKit computed. Mirrors
            /// `ClearComputedProps`.
            pub fn clear_computed_props(&mut self) {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_clear_computed_props(
                    self.props(),
                ))
            }
        }
    };
    ($type:ident $(<$lifetime:lifetime>)?) => {
        impl$(<$lifetime>)? $type$(<$lifetime>)? {
            /// Mirrors `SetProp`.
            pub fn set_prop(&mut self, key: &str, val: &str) {
                self.set_prop_with_computed(key, val, false)
            }

            /// Marks the property as computed when `computed`. Mirrors
            /// `SetProp(key, val, computed)`.
            pub fn set_prop_with_computed(&mut self, key: &str, val: &str, computed: bool) {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_set_prop(
                    self.props(),
                    key,
                    val,
                    computed,
                ))
            }

            /// Mirrors `SetIntProp`.
            pub fn set_int_prop(&mut self, key: &str, val: i32) {
                self.set_int_prop_with_computed(key, val, false)
            }

            /// Mirrors `SetIntProp(key, val, computed)`.
            pub fn set_int_prop_with_computed(&mut self, key: &str, val: i32, computed: bool) {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_set_int_prop(
                    self.props(),
                    key,
                    val,
                    computed,
                ))
            }

            /// Mirrors `SetUnsignedProp`.
            pub fn set_unsigned_prop(&mut self, key: &str, val: u32) {
                self.set_unsigned_prop_with_computed(key, val, false)
            }

            /// Mirrors `SetUnsignedProp(key, val, computed)`.
            pub fn set_unsigned_prop_with_computed(&mut self, key: &str, val: u32, computed: bool) {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_set_unsigned_prop(
                    self.props(),
                    key,
                    val,
                    computed,
                ))
            }

            /// Mirrors `SetDoubleProp`.
            pub fn set_double_prop(&mut self, key: &str, val: f64) {
                self.set_double_prop_with_computed(key, val, false)
            }

            /// Mirrors `SetDoubleProp(key, val, computed)`.
            pub fn set_double_prop_with_computed(&mut self, key: &str, val: f64, computed: bool) {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_set_double_prop(
                    self.props(),
                    key,
                    val,
                    computed,
                ))
            }

            /// Mirrors `SetBoolProp`.
            pub fn set_bool_prop(&mut self, key: &str, val: bool) {
                self.set_bool_prop_with_computed(key, val, false)
            }

            /// Mirrors `SetBoolProp(key, val, computed)`.
            pub fn set_bool_prop_with_computed(&mut self, key: &str, val: bool, computed: bool) {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_set_bool_prop(
                    self.props(),
                    key,
                    val,
                    computed,
                ))
            }

            /// Removes a property; a missing key is not an error. Mirrors
            /// `ClearProp`.
            pub fn clear_prop(&mut self, key: &str) {
                $crate::error::expect_ok($crate::ffi::props::ffi::props_clear_prop(
                    self.props(),
                    key,
                ))
            }
        }
    };
}

pub(crate) use {prop_getters, prop_setters};
