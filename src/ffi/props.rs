// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A property value and its key. `tag` selects the field that holds the
    /// value; see `PropValue::from_ffi`.
    struct PropValueData {
        key: UniquePtr<CxxString>,
        tag: u8,
        int_value: i32,
        unsigned_value: u32,
        double_value: f64,
        bool_value: bool,
        string_value: UniquePtr<CxxString>,
        int_vec: Vec<i32>,
        unsigned_vec: Vec<u32>,
        double_vec: Vec<f64>,
        string_vec: UniquePtr<CxxVector<CxxString>>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/props.h");

        type RDProps;
        type RWMol = crate::ffi::mol::ffi::RWMol;
        type Atom = crate::ffi::atom::ffi::Atom;
        type Bond = crate::ffi::bond::ffi::Bond;
        type Conformer = crate::ffi::conformer::ffi::Conformer;
        type SubstanceGroup = crate::ffi::substance_group::ffi::SubstanceGroup;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;

        fn mol_as_props(mol: &RWMol) -> Result<&RDProps>;
        fn atom_as_props(atom: &Atom) -> Result<&RDProps>;
        fn bond_as_props(bond: &Bond) -> Result<&RDProps>;
        fn conformer_as_props(conf: &Conformer) -> Result<&RDProps>;
        fn substance_group_as_props(sgroup: &SubstanceGroup) -> Result<&RDProps>;

        fn props_has_prop(obj: &RDProps, key: &str) -> Result<bool>;
        fn props_get_prop(obj: &RDProps, key: &str) -> Result<UniquePtr<CxxString>>;
        fn props_get_prop_auto_convert(obj: &RDProps, key: &str) -> Result<PropValueData>;
        fn props_get_int_prop(obj: &RDProps, key: &str) -> Result<i32>;
        fn props_get_unsigned_prop(obj: &RDProps, key: &str) -> Result<u32>;
        fn props_get_double_prop(obj: &RDProps, key: &str) -> Result<f64>;
        fn props_get_bool_prop(obj: &RDProps, key: &str) -> Result<bool>;
        fn props_get_explicit_bit_vect_prop(obj: &RDProps, key: &str)
        -> Result<UniquePtr<BitVect>>;
        fn props_get_unsigned_vect_prop(obj: &RDProps, key: &str) -> Result<Vec<u32>>;
        fn props_get_string_vect_prop(
            obj: &RDProps,
            key: &str,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn props_get_prop_names(
            obj: &RDProps,
            include_private: bool,
            include_computed: bool,
        ) -> Result<UniquePtr<CxxVector<CxxString>>>;
        fn props_get_props_as_dict(
            obj: &RDProps,
            include_private: bool,
            include_computed: bool,
            auto_convert_strings: bool,
        ) -> Result<Vec<PropValueData>>;

        fn props_set_prop(obj: &RDProps, key: &str, val: &str, computed: bool) -> Result<()>;
        fn props_set_int_prop(obj: &RDProps, key: &str, val: i32, computed: bool) -> Result<()>;
        fn props_set_unsigned_prop(
            obj: &RDProps,
            key: &str,
            val: u32,
            computed: bool,
        ) -> Result<()>;
        fn props_set_double_prop(obj: &RDProps, key: &str, val: f64, computed: bool) -> Result<()>;
        fn props_set_bool_prop(obj: &RDProps, key: &str, val: bool, computed: bool) -> Result<()>;
        fn props_set_explicit_bit_vect_prop(obj: &RDProps, key: &str, val: &BitVect) -> Result<()>;
        fn props_clear_prop(obj: &RDProps, key: &str) -> Result<()>;
        fn props_clear_computed_props(obj: &RDProps) -> Result<()>;
        fn mol_clear_computed_props(mol: &RWMol, include_rings: bool) -> Result<()>;
    }
}
