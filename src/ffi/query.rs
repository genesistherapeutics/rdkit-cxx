// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/query.h");

        type Atom = crate::ffi::atom::ffi::Atom;
        type Bond = crate::ffi::bond::ffi::Bond;
        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn atom_describe_query(atom: &Atom) -> Result<UniquePtr<CxxString>>;
        fn atom_get_query_type(atom: &Atom) -> Result<UniquePtr<CxxString>>;
        fn atom_match(atom: &Atom, other: &Atom) -> Result<bool>;
        fn atom_expand_query(
            atom: Pin<&mut Atom>,
            other: &Atom,
            how: i32,
            maintain_order: bool,
        ) -> Result<()>;
        fn atom_set_query(atom: Pin<&mut Atom>, other: &Atom) -> Result<()>;

        fn bond_describe_query(bond: &Bond) -> Result<UniquePtr<CxxString>>;
        fn bond_match(bond: &Bond, other: &Bond) -> Result<bool>;
        fn bond_expand_query(
            bond: Pin<&mut Bond>,
            other: &Bond,
            how: i32,
            maintain_order: bool,
        ) -> Result<()>;
        fn bond_set_query(bond: Pin<&mut Bond>, other: &Bond) -> Result<()>;

        // Only the crate's unit tests call these.
        #[allow(dead_code)]
        fn query_atom_without_query_for_tests() -> Result<UniquePtr<Atom>>;
        #[allow(dead_code)]
        fn query_bond_without_query_for_tests() -> Result<UniquePtr<Bond>>;

        fn mol_get_atoms_matching_query(mol: &RWMol, query: &Atom) -> Result<Vec<u32>>;
    }
}
