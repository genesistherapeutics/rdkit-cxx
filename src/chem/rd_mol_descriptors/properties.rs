// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! rdMolDescriptors's property registry: named descriptors to compute
//! together, and range queries on them.

use cxx::{CxxString, CxxVector, UniquePtr};

use crate::Result;
use crate::chem::Mol;
use crate::error::{expect_ok, to_string};
use crate::ffi::properties::{PropertyCallback, ffi};

/// Names cross as C++ strings (T2); every one came from Rust or RDKit.
fn names(names: UniquePtr<CxxVector<CxxString>>) -> Vec<String> {
    names
        .iter()
        .map(|name| name.to_str().expect("property names are UTF-8").to_owned())
        .collect()
}

/// A named descriptor in RDKit's property registry: one of RDKit's, or one
/// made from a Rust closure. Mirrors `rdMolDescriptors.PythonPropertyFunctor`,
/// the class Python subclasses to define one.
pub struct PythonPropertyFunctor {
    /// Never null.
    ptr: UniquePtr<ffi::PropertyFunctorHandle>,
}

impl std::fmt::Debug for PythonPropertyFunctor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PythonPropertyFunctor")
            .finish_non_exhaustive()
    }
}

impl PythonPropertyFunctor {
    /// A property computed by `call`. Python subclasses the class and
    /// overrides `__call__`; here the closure is the third argument. RDKit's
    /// registry is global, so any thread may call it: it is `Send + Sync`.
    /// It gets a copy of the molecule. A panic in it is caught, and the call
    /// that made it returns `Error::Runtime`. Mirrors
    /// `PythonPropertyFunctor(name, version)`.
    pub fn new<F>(name: &str, version: &str, call: F) -> PythonPropertyFunctor
    where
        F: Fn(&Mol) -> f64 + Send + Sync + 'static,
    {
        let callback = Box::new(PropertyCallback(Box::new(call)));
        PythonPropertyFunctor {
            ptr: expect_ok(ffi::property_functor_new(name, version, callback)),
        }
    }

    fn raw(&self) -> &ffi::PropertyFunctorHandle {
        self.ptr.as_ref().unwrap()
    }

    /// The property's name. Mirrors `GetName`.
    pub fn get_name(&self) -> String {
        to_string(expect_ok(ffi::property_functor_get_name(self.raw()))).expect("names are UTF-8")
    }

    /// The property's version. Mirrors `GetVersion`.
    pub fn get_version(&self) -> String {
        to_string(expect_ok(ffi::property_functor_get_version(self.raw())))
            .expect("versions are UTF-8")
    }

    /// The property's value for a molecule. Mirrors `__call__`.
    pub fn call(&self, mol: &Mol) -> Result<f64> {
        Ok(ffi::property_functor_call(self.raw(), mol.as_ffi())?)
    }
}

/// A set of registered properties to compute together. Mirrors
/// `rdMolDescriptors.Properties`.
///
/// The registry is process-wide. RDKit writes it without a lock, so the
/// wrappers hold one around every use (CONVENTIONS.md U4).
pub struct Properties {
    /// Never null.
    ptr: UniquePtr<ffi::Properties>,
}

impl std::fmt::Debug for Properties {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Properties").finish_non_exhaustive()
    }
}

/// Every registered property. Mirrors `Properties()`.
impl Default for Properties {
    fn default() -> Properties {
        Properties {
            ptr: expect_ok(ffi::properties_new()),
        }
    }
}

impl Properties {
    fn raw(&self) -> &ffi::Properties {
        self.ptr.as_ref().unwrap()
    }

    /// The named properties; `Error::Key` for one not registered. Mirrors
    /// `Properties(propNames)`.
    pub fn new(prop_names: &[&str]) -> Result<Properties> {
        Ok(Properties {
            ptr: ffi::properties_from_prop_names(prop_names)?,
        })
    }

    /// The properties' names, in order. Mirrors `GetPropertyNames`.
    pub fn get_property_names(&self) -> Vec<String> {
        names(expect_ok(ffi::properties_get_property_names(self.raw())))
    }

    /// The properties' values for a molecule, in order. Mirrors
    /// `ComputeProperties`.
    pub fn compute_properties(&self, mol: &Mol) -> Result<Vec<f64>> {
        Ok(ffi::properties_compute_properties(
            self.raw(),
            mol.as_ffi(),
        )?)
    }

    /// With `annotate_mol`, also sets each as an `f64` property of the
    /// molecule. Mirrors `ComputeProperties(mol, annotateMol)`.
    pub fn compute_properties_with_annotate_mol(
        &self,
        mol: &mut Mol,
        annotate_mol: bool,
    ) -> Result<Vec<f64>> {
        if !annotate_mol {
            return self.compute_properties(mol);
        }
        Ok(ffi::properties_compute_properties_annotate(
            self.raw(),
            mol.pin_mut(),
        )?)
    }

    /// Sets each property as an `f64` property of the molecule. Mirrors
    /// `AnnotateProperties`.
    pub fn annotate_properties(&self, mol: &mut Mol) -> Result<()> {
        Ok(ffi::properties_annotate_properties(
            self.raw(),
            mol.pin_mut(),
        )?)
    }

    /// The names of every registered property, RDKit's 43 first. Mirrors
    /// `Properties.GetAvailableProperties`.
    pub fn get_available_properties() -> Vec<String> {
        names(expect_ok(ffi::properties_get_available_properties()))
    }

    /// The registered property named `prop_name`; `Error::Key` for none.
    /// Mirrors `Properties.GetProperty`.
    pub fn get_property(prop_name: &str) -> Result<PythonPropertyFunctor> {
        Ok(PythonPropertyFunctor {
            ptr: ffi::properties_get_property(prop_name)?,
        })
    }

    /// Registers a property, in place of any of the same name, and returns
    /// its index in the registry. `Properties` made before keep the ones
    /// they had. Mirrors `Properties.RegisterProperty`.
    pub fn register_property(property_functor: &PythonPropertyFunctor) -> i32 {
        expect_ok(ffi::properties_register_property(property_functor.raw()))
    }
}

/// Whether a registered property of a molecule lies in a range. Mirrors
/// `rdMolDescriptors.PropertyRangeQuery`.
pub struct PropertyRangeQuery {
    /// Never null.
    ptr: UniquePtr<ffi::PropertyRangeQuery>,
}

impl std::fmt::Debug for PropertyRangeQuery {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PropertyRangeQuery").finish_non_exhaustive()
    }
}

impl PropertyRangeQuery {
    /// Whether the molecule's property is in range. `Error::Runtime` for a
    /// property made in Rust: the query calls RDKit's function for the
    /// property, which only RDKit's own have. Mirrors `Match`.
    pub fn match_(&self, what: &Mol) -> Result<bool> {
        Ok(ffi::property_range_query_match(
            self.ptr.as_ref().unwrap(),
            what.as_ffi(),
        )?)
    }
}

/// A query for the registered property `name` between `min` and `max`;
/// `Error::Key` for a name not registered. Mirrors
/// `rdMolDescriptors.MakePropertyRangeQuery`.
pub fn make_property_range_query(name: &str, min: f64, max: f64) -> Result<PropertyRangeQuery> {
    Ok(PropertyRangeQuery {
        ptr: ffi::make_property_range_query(name, min, max)?,
    })
}
