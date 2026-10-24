// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// The atom value an `rdqueries` value query compares.
    #[repr(u8)]
    enum AtomValueQuery {
        AtomNum,
        ExplicitValence,
        TotalValence,
        ExplicitDegree,
        TotalDegree,
        NonHydrogenDegree,
        HCount,
        Mass,
        Isotope,
        FormalCharge,
        Hybridization,
        InNRings,
        MinRingSize,
        RingBondCount,
        NumRadicalElectrons,
        NumHeteroatomNeighbors,
        NumAliphaticHeteroatomNeighbors,
    }

    /// How a value query compares; named as Python names the functions.
    #[repr(u8)]
    enum QueryComparison {
        Equals,
        Less,
        Greater,
    }

    /// The `rdqueries` queries that take no value.
    #[repr(u8)]
    enum AtomFlagQuery {
        IsUnsaturated,
        IsAromatic,
        IsAliphatic,
        IsInRing,
        HasChiralTag,
        MissingChiralTag,
        IsBridgehead,
        AAtom,
        AHAtom,
        QAtom,
        QHAtom,
        XAtom,
        XHAtom,
        MAtom,
        MHAtom,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/rdqueries.h");

        type Atom = crate::ffi::atom::ffi::Atom;
        type Bond = crate::ffi::bond::ffi::Bond;
        type RWMol = crate::ffi::mol::ffi::RWMol;
        type BitVect = crate::ffi::data_structs::ffi::BitVect;

        fn value_query_atom(
            kind: AtomValueQuery,
            comparison: QueryComparison,
            val: i32,
            negate: bool,
        ) -> Result<UniquePtr<Atom>>;
        fn flag_query_atom(kind: AtomFlagQuery, negate: bool) -> Result<UniquePtr<Atom>>;

        fn has_prop_query_atom(propname: &str, negate: bool) -> Result<UniquePtr<Atom>>;
        fn has_int_prop_with_value_query_atom(
            propname: &str,
            val: i32,
            negate: bool,
            tolerance: i32,
        ) -> Result<UniquePtr<Atom>>;
        fn has_bool_prop_with_value_query_atom(
            propname: &str,
            val: bool,
            negate: bool,
        ) -> Result<UniquePtr<Atom>>;
        fn has_string_prop_with_value_query_atom(
            propname: &str,
            val: &str,
            negate: bool,
        ) -> Result<UniquePtr<Atom>>;
        fn has_double_prop_with_value_query_atom(
            propname: &str,
            val: f64,
            negate: bool,
            tolerance: f64,
        ) -> Result<UniquePtr<Atom>>;
        fn has_bit_vect_prop_with_value_query_atom(
            propname: &str,
            val: &BitVect,
            negate: bool,
            tolerance: f64,
        ) -> Result<UniquePtr<Atom>>;

        fn has_prop_query_bond(propname: &str, negate: bool) -> Result<UniquePtr<Bond>>;
        fn has_int_prop_with_value_query_bond(
            propname: &str,
            val: i32,
            negate: bool,
            tolerance: i32,
        ) -> Result<UniquePtr<Bond>>;
        fn has_bool_prop_with_value_query_bond(
            propname: &str,
            val: bool,
            negate: bool,
        ) -> Result<UniquePtr<Bond>>;
        fn has_string_prop_with_value_query_bond(
            propname: &str,
            val: &str,
            negate: bool,
        ) -> Result<UniquePtr<Bond>>;
        fn has_double_prop_with_value_query_bond(
            propname: &str,
            val: f64,
            negate: bool,
            tolerance: f64,
        ) -> Result<UniquePtr<Bond>>;

        fn replace_atom_with_query_atom(mol: Pin<&mut RWMol>, idx: u32) -> Result<()>;
    }
}
