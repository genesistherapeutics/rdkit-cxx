// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdchem`.

mod atom;
mod bond;
mod conformer;
mod enums;
mod mol;
mod monomer_info;
mod periodic_table;
pub(crate) mod props;
mod ring_info;
mod rwmol;
mod substruct;

pub use atom::{__noatom_for_tests, Atom, AtomMut, AtomRef, GetSmartsParams};
pub use bond::{BondMut, BondRef};
pub use conformer::{Conformer, ConformerMut, ConformerRef};
pub use enums::{BondDir, BondStereo, BondType, ChiralType, HybridizationType, ValenceType};
pub use mol::{Atoms, Bonds, Conformers, Mol};
pub use monomer_info::{
    __monomer_type_round_trip_for_tests, AtomMonomerInfo, AtomMonomerInfoParams, AtomMonomerType,
    AtomPDBResidueInfo, AtomPDBResidueInfoParams,
};
pub use periodic_table::{PeriodicTable, get_periodic_table};
pub use props::{GetPropNamesParams, GetPropsAsDictParams, PropValue};
pub use ring_info::{RingInfoMut, RingInfoRef};
pub use rwmol::{RWMol, ReplaceAtomParams, ReplaceBondParams};
pub use substruct::SubstructMatchParameters;
