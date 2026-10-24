// Copyright (c) 2026 Steven Kearnes All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdchem`.

mod atom;
mod atom_labels;
mod bond;
mod conformer;
mod editable_mol;
mod enums;
mod mol;
mod mol_bundle;
mod monomer_info;
mod periodic_table;
mod pickle;
pub(crate) mod props;
mod ring_info;
mod rwmol;
mod stereo;
mod substance_group;
mod substruct;

pub use atom::{__noatom_for_tests, Atom, AtomMut, AtomRef, GetSmartsParams};
pub use atom_labels::{
    AsAtom, AsAtomMut, get_atom_alias, get_atom_r_label, get_atom_value, get_num_pi_electrons,
    get_supplemental_smiles_label, set_atom_alias, set_atom_r_label, set_atom_value,
    set_supplemental_smiles_label,
};
pub use bond::{BondMut, BondRef};
pub use conformer::{Conformer, ConformerMut, ConformerRef};
pub use editable_mol::EditableMol;
pub use enums::{
    BondDir, BondStereo, BondType, ChiralType, HybridizationType, StereoDescriptor,
    StereoGroupType, StereoSpecified, StereoType, ValenceType,
};
pub use mol::{Atoms, Bonds, Conformers, Mol};
pub use mol_bundle::{FixedMolSizeMolBundle, MolBundle, mol_bundle_can_serialize};
pub use monomer_info::{
    __monomer_type_round_trip_for_tests, AtomMonomerInfo, AtomMonomerInfoParams, AtomMonomerType,
    AtomPDBResidueInfo, AtomPDBResidueInfoParams,
};
pub use periodic_table::{PeriodicTable, get_periodic_table};
pub use pickle::{
    __property_pickle_option_values_for_tests, PropertyPickleOptions,
    get_default_pickle_properties, set_default_pickle_properties,
};
pub use props::{GetPropNamesParams, GetPropsAsDictParams, PropValue};
pub use ring_info::{RingInfoMut, RingInfoRef};
pub use rwmol::{RWMol, ReplaceAtomParams, ReplaceBondParams};
pub use stereo::{
    CreateStereoGroupParams, StereoGroup, StereoInfo, create_stereo_group,
    create_stereo_group_with_params, forward_stereo_group_ids,
};
pub use substance_group::{
    SubstanceGroupAttach, SubstanceGroupCState, SubstanceGroupMut, SubstanceGroupRef,
    add_mol_substance_group, clear_mol_substance_groups, create_mol_data_substance_group,
    create_mol_substance_group, get_mol_substance_group_with_idx,
    get_mol_substance_group_with_idx_mut, get_mol_substance_groups,
};
pub use substruct::SubstructMatchParameters;
