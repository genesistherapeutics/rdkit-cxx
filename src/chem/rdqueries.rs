// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdqueries`: query atoms and bonds.
//!
//! Each function returns an `Atom` or `Bond` holding the query (Python's
//! `QueryAtom` and `QueryBond`), to match against with `match_`, add to a
//! molecule, or combine with `expand_query`. As in Python, `..._less_...`
//! matches atoms whose value is less than `val`, and `..._greater_...` those
//! whose value is greater.

use crate::Result;
use crate::chem::rdchem::Atom;
use crate::chem::rdchem::{AtomMut, Bond, Mol};
use crate::data_structs::ExplicitBitVect;
use crate::ffi::rdqueries::ffi;

/// Defines the equals, less, and greater query functions for an atom value,
/// each with a `_with_negate` form (A2).
macro_rules! value_queries {
    ($(
        $kind:ident, $what:literal:
        $equals:ident, $equals_negate:ident,
        $less:ident, $less_negate:ident,
        $greater:ident, $greater_negate:ident;
    )*) => {
        $(
            value_query!($kind, Equals, $what, "equals", $equals, $equals_negate);
            value_query!($kind, Less, $what, "is less than", $less, $less_negate);
            value_query!($kind, Greater, $what, "is greater than", $greater, $greater_negate);
        )*
    };
}

macro_rules! value_query {
    ($kind:ident, $comparison:ident, $what:literal, $relation:literal, $name:ident, $negate:ident) => {
        #[doc = concat!(
                    "A query atom that matches atoms whose ", $what, " ", $relation,
                    " `val`. Mirrors `rdqueries.", stringify!($kind), stringify!($comparison),
                    "QueryAtom`. `Error::Value` for `val` beyond 2^30 (2^30 / 1000 for masses), ",
                    "where RDKit's `int` comparison overflows (CONVENTIONS.md U7)."
                )]
        pub fn $name(val: i32) -> Result<Atom> {
            $negate(val, false)
        }

        #[doc = concat!(
                    "Mirrors `rdqueries.", stringify!($kind), stringify!($comparison),
                    "QueryAtom(val, negate)`; with `negate`, the query matches the other atoms."
                )]
        pub fn $negate(val: i32, negate: bool) -> Result<Atom> {
            Ok(Atom::from_ffi(ffi::value_query_atom(
                ffi::AtomValueQuery::$kind,
                ffi::QueryComparison::$comparison,
                val,
                negate,
            )?))
        }
    };
}

value_queries! {
    AtomNum, "atomic number":
        atom_num_equals_query_atom, atom_num_equals_query_atom_with_negate,
        atom_num_less_query_atom, atom_num_less_query_atom_with_negate,
        atom_num_greater_query_atom, atom_num_greater_query_atom_with_negate;
    ExplicitValence, "explicit valence":
        explicit_valence_equals_query_atom, explicit_valence_equals_query_atom_with_negate,
        explicit_valence_less_query_atom, explicit_valence_less_query_atom_with_negate,
        explicit_valence_greater_query_atom, explicit_valence_greater_query_atom_with_negate;
    TotalValence, "total valence":
        total_valence_equals_query_atom, total_valence_equals_query_atom_with_negate,
        total_valence_less_query_atom, total_valence_less_query_atom_with_negate,
        total_valence_greater_query_atom, total_valence_greater_query_atom_with_negate;
    ExplicitDegree, "explicit degree":
        explicit_degree_equals_query_atom, explicit_degree_equals_query_atom_with_negate,
        explicit_degree_less_query_atom, explicit_degree_less_query_atom_with_negate,
        explicit_degree_greater_query_atom, explicit_degree_greater_query_atom_with_negate;
    TotalDegree, "total degree":
        total_degree_equals_query_atom, total_degree_equals_query_atom_with_negate,
        total_degree_less_query_atom, total_degree_less_query_atom_with_negate,
        total_degree_greater_query_atom, total_degree_greater_query_atom_with_negate;
    NonHydrogenDegree, "number of heavy-atom neighbors":
        non_hydrogen_degree_equals_query_atom, non_hydrogen_degree_equals_query_atom_with_negate,
        non_hydrogen_degree_less_query_atom, non_hydrogen_degree_less_query_atom_with_negate,
        non_hydrogen_degree_greater_query_atom, non_hydrogen_degree_greater_query_atom_with_negate;
    HCount, "total hydrogen count":
        h_count_equals_query_atom, h_count_equals_query_atom_with_negate,
        h_count_less_query_atom, h_count_less_query_atom_with_negate,
        h_count_greater_query_atom, h_count_greater_query_atom_with_negate;
    Mass, "mass, rounded to an integer,":
        mass_equals_query_atom, mass_equals_query_atom_with_negate,
        mass_less_query_atom, mass_less_query_atom_with_negate,
        mass_greater_query_atom, mass_greater_query_atom_with_negate;
    Isotope, "isotope":
        isotope_equals_query_atom, isotope_equals_query_atom_with_negate,
        isotope_less_query_atom, isotope_less_query_atom_with_negate,
        isotope_greater_query_atom, isotope_greater_query_atom_with_negate;
    FormalCharge, "formal charge":
        formal_charge_equals_query_atom, formal_charge_equals_query_atom_with_negate,
        formal_charge_less_query_atom, formal_charge_less_query_atom_with_negate,
        formal_charge_greater_query_atom, formal_charge_greater_query_atom_with_negate;
    Hybridization, "hybridization, as its `HybridizationType` value,":
        hybridization_equals_query_atom, hybridization_equals_query_atom_with_negate,
        hybridization_less_query_atom, hybridization_less_query_atom_with_negate,
        hybridization_greater_query_atom, hybridization_greater_query_atom_with_negate;
    InNRings, "number of rings":
        in_n_rings_equals_query_atom, in_n_rings_equals_query_atom_with_negate,
        in_n_rings_less_query_atom, in_n_rings_less_query_atom_with_negate,
        in_n_rings_greater_query_atom, in_n_rings_greater_query_atom_with_negate;
    MinRingSize, "smallest ring size":
        min_ring_size_equals_query_atom, min_ring_size_equals_query_atom_with_negate,
        min_ring_size_less_query_atom, min_ring_size_less_query_atom_with_negate,
        min_ring_size_greater_query_atom, min_ring_size_greater_query_atom_with_negate;
    RingBondCount, "number of ring bonds":
        ring_bond_count_equals_query_atom, ring_bond_count_equals_query_atom_with_negate,
        ring_bond_count_less_query_atom, ring_bond_count_less_query_atom_with_negate,
        ring_bond_count_greater_query_atom, ring_bond_count_greater_query_atom_with_negate;
    NumRadicalElectrons, "number of radical electrons":
        num_radical_electrons_equals_query_atom, num_radical_electrons_equals_query_atom_with_negate,
        num_radical_electrons_less_query_atom, num_radical_electrons_less_query_atom_with_negate,
        num_radical_electrons_greater_query_atom, num_radical_electrons_greater_query_atom_with_negate;
    NumHeteroatomNeighbors, "number of heteroatom neighbors":
        num_heteroatom_neighbors_equals_query_atom,
        num_heteroatom_neighbors_equals_query_atom_with_negate,
        num_heteroatom_neighbors_less_query_atom,
        num_heteroatom_neighbors_less_query_atom_with_negate,
        num_heteroatom_neighbors_greater_query_atom,
        num_heteroatom_neighbors_greater_query_atom_with_negate;
    NumAliphaticHeteroatomNeighbors, "number of aliphatic heteroatom neighbors":
        num_aliphatic_heteroatom_neighbors_equals_query_atom,
        num_aliphatic_heteroatom_neighbors_equals_query_atom_with_negate,
        num_aliphatic_heteroatom_neighbors_less_query_atom,
        num_aliphatic_heteroatom_neighbors_less_query_atom_with_negate,
        num_aliphatic_heteroatom_neighbors_greater_query_atom,
        num_aliphatic_heteroatom_neighbors_greater_query_atom_with_negate;
}

/// Defines a query function that takes no value, with a `_with_negate`
/// form (A2).
macro_rules! flag_queries {
    ($($kind:ident, $python:literal, $what:literal: $name:ident, $negate:ident;)*) => {
        $(
            #[doc = concat!(
                "A query atom that matches ", $what, ". Mirrors `rdqueries.", $python, "`."
            )]
            pub fn $name() -> Atom {
                $negate(false)
            }

            #[doc = concat!(
                "Mirrors `rdqueries.", $python,
                "(negate)`; with `negate`, the query matches the other atoms."
            )]
            pub fn $negate(negate: bool) -> Atom {
                Atom::from_ffi(crate::error::expect_ok(ffi::flag_query_atom(
                    ffi::AtomFlagQuery::$kind,
                    negate,
                )))
            }
        )*
    };
}

flag_queries! {
    IsUnsaturated, "IsUnsaturatedQueryAtom", "atoms with a multiple bond":
        is_unsaturated_query_atom, is_unsaturated_query_atom_with_negate;
    IsAromatic, "IsAromaticQueryAtom", "aromatic atoms":
        is_aromatic_query_atom, is_aromatic_query_atom_with_negate;
    IsAliphatic, "IsAliphaticQueryAtom", "aliphatic atoms":
        is_aliphatic_query_atom, is_aliphatic_query_atom_with_negate;
    IsInRing, "IsInRingQueryAtom", "ring atoms":
        is_in_ring_query_atom, is_in_ring_query_atom_with_negate;
    HasChiralTag, "HasChiralTagQueryAtom", "atoms with a chiral tag":
        has_chiral_tag_query_atom, has_chiral_tag_query_atom_with_negate;
    MissingChiralTag, "MissingChiralTagQueryAtom",
        "potential stereocenters without a chiral tag":
        missing_chiral_tag_query_atom, missing_chiral_tag_query_atom_with_negate;
    IsBridgehead, "IsBridgeheadQueryAtom", "bridgehead atoms":
        is_bridgehead_query_atom, is_bridgehead_query_atom_with_negate;
    AAtom, "AAtomQueryAtom", "any heavy atom (the `A` query)":
        a_atom_query_atom, a_atom_query_atom_with_negate;
    AHAtom, "AHAtomQueryAtom", "any atom (the `AH` query)":
        ah_atom_query_atom, ah_atom_query_atom_with_negate;
    QAtom, "QAtomQueryAtom", "any heteroatom but hydrogen (the `Q` query)":
        q_atom_query_atom, q_atom_query_atom_with_negate;
    QHAtom, "QHAtomQueryAtom", "any heteroatom (the `QH` query)":
        qh_atom_query_atom, qh_atom_query_atom_with_negate;
    XAtom, "XAtomQueryAtom", "halogens (the `X` query)":
        x_atom_query_atom, x_atom_query_atom_with_negate;
    XHAtom, "XHAtomQueryAtom", "halogens and hydrogen (the `XH` query)":
        xh_atom_query_atom, xh_atom_query_atom_with_negate;
    MAtom, "MAtomQueryAtom", "metals (the `M` query)":
        m_atom_query_atom, m_atom_query_atom_with_negate;
    MHAtom, "MHAtomQueryAtom", "metals and hydrogen (the `MH` query)":
        mh_atom_query_atom, mh_atom_query_atom_with_negate;
}

/// A query atom that matches atoms with the property. Mirrors
/// `rdqueries.HasPropQueryAtom`.
pub fn has_prop_query_atom(propname: &str) -> Atom {
    has_prop_query_atom_with_negate(propname, false)
}

/// Mirrors `rdqueries.HasPropQueryAtom(propname, negate)`.
pub fn has_prop_query_atom_with_negate(propname: &str, negate: bool) -> Atom {
    Atom::from_ffi(crate::error::expect_ok(ffi::has_prop_query_atom(
        propname, negate,
    )))
}

/// A query bond that matches bonds with the property. Mirrors
/// `rdqueries.HasPropQueryBond`.
pub fn has_prop_query_bond(propname: &str) -> Bond {
    has_prop_query_bond_with_negate(propname, false)
}

/// Mirrors `rdqueries.HasPropQueryBond(propname, negate)`.
pub fn has_prop_query_bond_with_negate(propname: &str, negate: bool) -> Bond {
    Bond::from_ffi(crate::error::expect_ok(ffi::has_prop_query_bond(
        propname, negate,
    )))
}

/// Parameters for `has_int_prop_with_value_query_atom_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HasIntPropWithValueQueryAtomParams {
    /// Matches the other atoms.
    pub negate: bool,
    /// How far the property may be from `val`.
    pub tolerance: i32,
}

/// Parameters for `has_int_prop_with_value_query_bond_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HasIntPropWithValueQueryBondParams {
    /// Matches the other bonds.
    pub negate: bool,
    /// How far the property may be from `val`.
    pub tolerance: i32,
}

/// Parameters for `has_double_prop_with_value_query_atom_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HasDoublePropWithValueQueryAtomParams {
    /// Matches the other atoms.
    pub negate: bool,
    /// How far the property may be from `val`.
    pub tolerance: f64,
}

/// Parameters for `has_double_prop_with_value_query_bond_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HasDoublePropWithValueQueryBondParams {
    /// Matches the other bonds.
    pub negate: bool,
    /// How far the property may be from `val`.
    pub tolerance: f64,
}

/// Parameters for `has_bit_vect_prop_with_value_query_atom_with_params`,
/// with Python's defaults (A2).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HasBitVectPropWithValueQueryAtomParams {
    /// Matches the other atoms.
    pub negate: bool,
    /// The largest Tanimoto distance from `val` that matches.
    pub tolerance: f64,
}

/// A query atom that matches atoms whose int property is `val`. Mirrors
/// `rdqueries.HasIntPropWithValueQueryAtom`. The comparison is made in 64
/// bits, where RDKit's own would overflow `int` for distant values.
pub fn has_int_prop_with_value_query_atom(propname: &str, val: i32) -> Atom {
    has_int_prop_with_value_query_atom_with_params(
        propname,
        val,
        &HasIntPropWithValueQueryAtomParams::default(),
    )
}

/// Mirrors `rdqueries.HasIntPropWithValueQueryAtom(propname, val, negate,
/// tolerance)`.
pub fn has_int_prop_with_value_query_atom_with_params(
    propname: &str,
    val: i32,
    params: &HasIntPropWithValueQueryAtomParams,
) -> Atom {
    Atom::from_ffi(crate::error::expect_ok(
        ffi::has_int_prop_with_value_query_atom(propname, val, params.negate, params.tolerance),
    ))
}

/// A query bond that matches bonds whose int property is `val`. Mirrors
/// `rdqueries.HasIntPropWithValueQueryBond`, compared in 64 bits as
/// `has_int_prop_with_value_query_atom` is.
pub fn has_int_prop_with_value_query_bond(propname: &str, val: i32) -> Bond {
    has_int_prop_with_value_query_bond_with_params(
        propname,
        val,
        &HasIntPropWithValueQueryBondParams::default(),
    )
}

/// Mirrors `rdqueries.HasIntPropWithValueQueryBond(propname, val, negate,
/// tolerance)`.
pub fn has_int_prop_with_value_query_bond_with_params(
    propname: &str,
    val: i32,
    params: &HasIntPropWithValueQueryBondParams,
) -> Bond {
    Bond::from_ffi(crate::error::expect_ok(
        ffi::has_int_prop_with_value_query_bond(propname, val, params.negate, params.tolerance),
    ))
}

/// A query atom that matches atoms whose bool property is `val`. Mirrors
/// `rdqueries.HasBoolPropWithValueQueryAtom`.
pub fn has_bool_prop_with_value_query_atom(propname: &str, val: bool) -> Atom {
    has_bool_prop_with_value_query_atom_with_negate(propname, val, false)
}

/// Mirrors `rdqueries.HasBoolPropWithValueQueryAtom(propname, val, negate)`.
pub fn has_bool_prop_with_value_query_atom_with_negate(
    propname: &str,
    val: bool,
    negate: bool,
) -> Atom {
    Atom::from_ffi(crate::error::expect_ok(
        ffi::has_bool_prop_with_value_query_atom(propname, val, negate),
    ))
}

/// A query bond that matches bonds whose bool property is `val`. Mirrors
/// `rdqueries.HasBoolPropWithValueQueryBond`.
pub fn has_bool_prop_with_value_query_bond(propname: &str, val: bool) -> Bond {
    has_bool_prop_with_value_query_bond_with_negate(propname, val, false)
}

/// Mirrors `rdqueries.HasBoolPropWithValueQueryBond(propname, val, negate)`.
pub fn has_bool_prop_with_value_query_bond_with_negate(
    propname: &str,
    val: bool,
    negate: bool,
) -> Bond {
    Bond::from_ffi(crate::error::expect_ok(
        ffi::has_bool_prop_with_value_query_bond(propname, val, negate),
    ))
}

/// A query atom that matches atoms whose string property is `val`. Mirrors
/// `rdqueries.HasStringPropWithValueQueryAtom`.
pub fn has_string_prop_with_value_query_atom(propname: &str, val: &str) -> Atom {
    has_string_prop_with_value_query_atom_with_negate(propname, val, false)
}

/// Mirrors `rdqueries.HasStringPropWithValueQueryAtom(propname, val,
/// negate)`.
pub fn has_string_prop_with_value_query_atom_with_negate(
    propname: &str,
    val: &str,
    negate: bool,
) -> Atom {
    Atom::from_ffi(crate::error::expect_ok(
        ffi::has_string_prop_with_value_query_atom(propname, val, negate),
    ))
}

/// A query bond that matches bonds whose string property is `val`. Mirrors
/// `rdqueries.HasStringPropWithValueQueryBond`.
pub fn has_string_prop_with_value_query_bond(propname: &str, val: &str) -> Bond {
    has_string_prop_with_value_query_bond_with_negate(propname, val, false)
}

/// Mirrors `rdqueries.HasStringPropWithValueQueryBond(propname, val,
/// negate)`.
pub fn has_string_prop_with_value_query_bond_with_negate(
    propname: &str,
    val: &str,
    negate: bool,
) -> Bond {
    Bond::from_ffi(crate::error::expect_ok(
        ffi::has_string_prop_with_value_query_bond(propname, val, negate),
    ))
}

/// A query atom that matches atoms whose double property is `val`. Mirrors
/// `rdqueries.HasDoublePropWithValueQueryAtom`.
pub fn has_double_prop_with_value_query_atom(propname: &str, val: f64) -> Atom {
    has_double_prop_with_value_query_atom_with_params(
        propname,
        val,
        &HasDoublePropWithValueQueryAtomParams::default(),
    )
}

/// Mirrors `rdqueries.HasDoublePropWithValueQueryAtom(propname, val, negate,
/// tolerance)`.
pub fn has_double_prop_with_value_query_atom_with_params(
    propname: &str,
    val: f64,
    params: &HasDoublePropWithValueQueryAtomParams,
) -> Atom {
    Atom::from_ffi(crate::error::expect_ok(
        ffi::has_double_prop_with_value_query_atom(propname, val, params.negate, params.tolerance),
    ))
}

/// A query bond that matches bonds whose double property is `val`. Mirrors
/// `rdqueries.HasDoublePropWithValueQueryBond`.
pub fn has_double_prop_with_value_query_bond(propname: &str, val: f64) -> Bond {
    has_double_prop_with_value_query_bond_with_params(
        propname,
        val,
        &HasDoublePropWithValueQueryBondParams::default(),
    )
}

/// Mirrors `rdqueries.HasDoublePropWithValueQueryBond(propname, val, negate,
/// tolerance)`.
pub fn has_double_prop_with_value_query_bond_with_params(
    propname: &str,
    val: f64,
    params: &HasDoublePropWithValueQueryBondParams,
) -> Bond {
    Bond::from_ffi(crate::error::expect_ok(
        ffi::has_double_prop_with_value_query_bond(propname, val, params.negate, params.tolerance),
    ))
}

/// A query atom that matches atoms whose bit-vector property is `val`.
/// Mirrors `rdqueries.HasBitVectPropWithValueQueryAtom`.
pub fn has_bit_vect_prop_with_value_query_atom(propname: &str, val: &ExplicitBitVect) -> Atom {
    has_bit_vect_prop_with_value_query_atom_with_params(
        propname,
        val,
        &HasBitVectPropWithValueQueryAtomParams::default(),
    )
}

/// Mirrors `rdqueries.HasBitVectPropWithValueQueryAtom(propname, val, negate,
/// tolerance)`.
pub fn has_bit_vect_prop_with_value_query_atom_with_params(
    propname: &str,
    val: &ExplicitBitVect,
    params: &HasBitVectPropWithValueQueryAtomParams,
) -> Atom {
    Atom::from_ffi(crate::error::expect_ok(
        ffi::has_bit_vect_prop_with_value_query_atom(
            propname,
            val.raw(),
            params.negate,
            params.tolerance,
        ),
    ))
}

/// Makes the atom at `idx` a query atom that matches what it is, and
/// returns it to take more queries; an atom that has a query already is
/// left alone. Mirrors `rdqueries.ReplaceAtomWithQueryAtom`, with the atom
/// given by index. `Error::Index` for an index past the end.
pub fn replace_atom_with_query_atom(mol: &mut Mol, idx: u32) -> Result<AtomMut<'_>> {
    ffi::replace_atom_with_query_atom(mol.pin_mut(), idx)?;
    mol.get_atom_with_idx_mut(idx)
}
