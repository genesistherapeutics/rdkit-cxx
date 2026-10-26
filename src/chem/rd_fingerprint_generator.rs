// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Mirrors `rdkit.Chem.rdFingerprintGenerator`: Morgan, atom-pair, RDKit,
//! and topological-torsion fingerprint generators.
//!
//! The generators are 64-bit (`FingerprintGenerator64`), as every Python
//! factory returns. The wrapper checks what RDKit trusts (CONVENTIONS.md U1):
//! atom indices and custom invariants per call, and per generator a
//! fingerprint size within `i32::MAX`, count bounds when count simulation is
//! on (RDKit divides by their number), and a torsion length whose codes fit
//! in 64 bits.

use std::collections::BTreeMap;
use std::fmt;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::enums::rdkit_enum;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::data_structs::{ExplicitBitVect, SparseBitVect, UIntSparseIntVect, ULongSparseIntVect};
use crate::error::{expect_ok, to_string};
use crate::ffi::fingerprint_generator::ffi;
use crate::ffi::geometry::ffi::IndexList;

rdkit_enum! {
    /// A fingerprint type. Mirrors `rdFingerprintGenerator.FPType`.
    FPType: ffi::FPType {
        AtomPairFp = AtomPairFP,
        MorganFp = MorganFP,
        RdkitFp = RDKitFP,
        TopologicalTorsionFp = TopologicalTorsionFP,
    }
}

/// Computes atom invariants for a fingerprint generator. Mirrors
/// `rdFingerprintGenerator.AtomInvariantsGenerator`.
pub struct AtomInvariantsGenerator {
    /// Never null.
    ptr: UniquePtr<ffi::AtomInvGen>,
}

impl AtomInvariantsGenerator {
    fn from_ffi(ptr: UniquePtr<ffi::AtomInvGen>) -> AtomInvariantsGenerator {
        assert!(!ptr.is_null(), "RDKit returned a null invariants generator");
        AtomInvariantsGenerator { ptr }
    }

    fn clone_ffi(&self) -> UniquePtr<ffi::AtomInvGen> {
        expect_ok(ffi::atom_inv_gen_clone(&self.ptr))
    }
}

/// A copy; generators take copies, as Python's do.
impl Clone for AtomInvariantsGenerator {
    fn clone(&self) -> AtomInvariantsGenerator {
        AtomInvariantsGenerator::from_ffi(self.clone_ffi())
    }
}

impl fmt::Debug for AtomInvariantsGenerator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let info = to_string(expect_ok(ffi::atom_inv_gen_info_string(&self.ptr)));
        f.debug_tuple("AtomInvariantsGenerator")
            .field(&info)
            .finish()
    }
}

/// Computes bond invariants for a fingerprint generator. Mirrors
/// `rdFingerprintGenerator.BondInvariantsGenerator`.
pub struct BondInvariantsGenerator {
    /// Never null.
    ptr: UniquePtr<ffi::BondInvGen>,
}

impl BondInvariantsGenerator {
    fn from_ffi(ptr: UniquePtr<ffi::BondInvGen>) -> BondInvariantsGenerator {
        assert!(!ptr.is_null(), "RDKit returned a null invariants generator");
        BondInvariantsGenerator { ptr }
    }

    fn clone_ffi(&self) -> UniquePtr<ffi::BondInvGen> {
        expect_ok(ffi::bond_inv_gen_clone(&self.ptr))
    }
}

/// A copy; generators take copies, as Python's do.
impl Clone for BondInvariantsGenerator {
    fn clone(&self) -> BondInvariantsGenerator {
        BondInvariantsGenerator::from_ffi(self.clone_ffi())
    }
}

impl fmt::Debug for BondInvariantsGenerator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let info = to_string(expect_ok(ffi::bond_inv_gen_info_string(&self.ptr)));
        f.debug_tuple("BondInvariantsGenerator")
            .field(&info)
            .finish()
    }
}

fn atom_inv_gen_ffi(generator: &Option<AtomInvariantsGenerator>) -> UniquePtr<ffi::AtomInvGen> {
    generator
        .as_ref()
        .map_or_else(UniquePtr::null, AtomInvariantsGenerator::clone_ffi)
}

/// Morgan's atom invariants: the ECFP connectivity invariants. Mirrors
/// `rdFingerprintGenerator.GetMorganAtomInvGen`.
pub fn get_morgan_atom_inv_gen(include_ring_membership: bool) -> AtomInvariantsGenerator {
    AtomInvariantsGenerator::from_ffi(expect_ok(ffi::get_morgan_atom_inv_gen(
        include_ring_membership,
    )))
}

/// Morgan's feature invariants (FCFP), from RDKit's default feature
/// patterns. Mirrors `rdFingerprintGenerator.GetMorganFeatureAtomInvGen`.
pub fn get_morgan_feature_atom_inv_gen() -> AtomInvariantsGenerator {
    let mut none = expect_ok(crate::ffi::molops::ffi::borrowed_mols_new());
    AtomInvariantsGenerator::from_ffi(expect_ok(ffi::get_morgan_feature_atom_inv_gen(
        false,
        none.pin_mut(),
    )))
}

/// Morgan's feature invariants, one bit for each query in `patterns`, which
/// the generator copies. Mirrors
/// `rdFingerprintGenerator.GetMorganFeatureAtomInvGen(patterns)`.
pub fn get_morgan_feature_atom_inv_gen_with_patterns(patterns: &[Mol]) -> AtomInvariantsGenerator {
    AtomInvariantsGenerator::from_ffi(expect_ok(with_borrowed_mols_shared(
        patterns.iter().map(|m| ("", m)),
        |list| ffi::get_morgan_feature_atom_inv_gen(true, list),
    )))
}

/// Parameters for `get_morgan_bond_inv_gen_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct GetMorganBondInvGenParams {
    /// Includes bond types in the invariants.
    pub use_bond_types: bool,
    /// Includes bond stereo in the invariants.
    pub use_chirality: bool,
}

impl Default for GetMorganBondInvGenParams {
    fn default() -> GetMorganBondInvGenParams {
        GetMorganBondInvGenParams {
            use_bond_types: true,
            use_chirality: false,
        }
    }
}

/// Morgan's bond invariants. Mirrors
/// `rdFingerprintGenerator.GetMorganBondInvGen`.
pub fn get_morgan_bond_inv_gen() -> BondInvariantsGenerator {
    get_morgan_bond_inv_gen_with_params(&GetMorganBondInvGenParams::default())
}

/// Mirrors `rdFingerprintGenerator.GetMorganBondInvGen(useBondTypes,
/// useChirality)`.
pub fn get_morgan_bond_inv_gen_with_params(
    params: &GetMorganBondInvGenParams,
) -> BondInvariantsGenerator {
    BondInvariantsGenerator::from_ffi(expect_ok(ffi::get_morgan_bond_inv_gen(
        params.use_bond_types,
        params.use_chirality,
    )))
}

/// The atom-pair atom invariants. Mirrors
/// `rdFingerprintGenerator.GetAtomPairAtomInvGen`.
pub fn get_atom_pair_atom_inv_gen() -> AtomInvariantsGenerator {
    get_atom_pair_atom_inv_gen_with_include_chirality(false)
}

/// Mirrors `rdFingerprintGenerator.GetAtomPairAtomInvGen(includeChirality)`.
pub fn get_atom_pair_atom_inv_gen_with_include_chirality(
    include_chirality: bool,
) -> AtomInvariantsGenerator {
    AtomInvariantsGenerator::from_ffi(expect_ok(ffi::get_atom_pair_atom_inv_gen(
        include_chirality,
    )))
}

/// The RDKit fingerprint's atom invariants. Mirrors
/// `rdFingerprintGenerator.GetRDKitAtomInvGen`.
pub fn get_rdkit_atom_inv_gen() -> AtomInvariantsGenerator {
    AtomInvariantsGenerator::from_ffi(expect_ok(ffi::get_rdkit_atom_inv_gen()))
}

/// The Morgan-specific options. Mirrors
/// `rdFingerprintGenerator.MorganFingerprintOptions`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MorganFingerprintOptions {
    /// The radius of the environments.
    pub radius: u32,
    /// Leaves out atoms whose invariant is zero.
    pub only_nonzero_invariants: bool,
    /// Keeps environments that repeat a smaller one.
    pub include_redundant_environments: bool,
}

/// The atom-pair-specific options. Mirrors
/// `rdFingerprintGenerator.AtomPairFingerprintOptions`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AtomPairFingerprintOptions {
    /// The shortest distance between a pair's atoms.
    pub min_distance: u32,
    /// The longest distance between a pair's atoms.
    pub max_distance: u32,
    /// Uses topological distances; otherwise 3D ones.
    pub use_2d: bool,
}

/// The RDKit-fingerprint-specific options. Mirrors
/// `rdFingerprintGenerator.RDKitFingerprintOptions`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RDKitFingerprintOptions {
    /// The fewest bonds in a path.
    pub min_path: u32,
    /// The most bonds in a path.
    pub max_path: u32,
    /// Includes explicit hydrogens in the paths.
    pub use_hs: bool,
    /// Uses branched subgraphs, not only linear paths.
    pub branched_paths: bool,
    /// Includes bond orders in the path hashes.
    pub use_bond_order: bool,
}

/// The topological-torsion-specific options. Mirrors
/// `rdFingerprintGenerator.TopologicalTorsionFingerprintOptions`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TopologicalTorsionFingerprintOptions {
    /// The number of atoms in each torsion.
    pub torsion_atom_count: u32,
    /// Only torsions that are shortest paths between their ends.
    pub only_shortest_paths: bool,
}

/// A generator's type-specific options: Python's subclasses of
/// `FingerprintOptions`.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum FingerprintKindOptions {
    /// Mirrors `MorganFingerprintOptions`.
    Morgan(MorganFingerprintOptions),
    /// Mirrors `AtomPairFingerprintOptions`.
    AtomPair(AtomPairFingerprintOptions),
    /// Mirrors `RDKitFingerprintOptions`.
    Rdkit(RDKitFingerprintOptions),
    /// Mirrors `TopologicalTorsionFingerprintOptions`.
    TopologicalTorsion(TopologicalTorsionFingerprintOptions),
}

/// A generator's options. Mirrors `rdFingerprintGenerator.FingerprintOptions`
/// and, in `kind`, its subclass's fields.
#[derive(Clone, Debug, PartialEq)]
pub struct FingerprintOptions {
    /// Simulates counts in bit fingerprints, one bit per count bound.
    pub count_simulation: bool,
    /// Includes chirality in the atom invariants.
    pub include_chirality: bool,
    /// The number of bits in folded fingerprints.
    pub fp_size: u32,
    /// The bits each feature sets.
    pub num_bits_per_feature: u32,
    /// The count bounds for count simulation.
    pub count_bounds: Vec<u32>,
    /// The type-specific options.
    pub kind: FingerprintKindOptions,
}

impl FingerprintOptions {
    /// Sets the count bounds. Mirrors `FingerprintOptions.SetCountBounds`;
    /// pass the options to `FingerprintGenerator64::set_options`.
    pub fn set_count_bounds(&mut self, bounds: &[u32]) {
        self.count_bounds = bounds.to_vec();
    }

    fn from_ffi(data: ffi::FingerprintOptionsData) -> FingerprintOptions {
        let kind = match FPType::from_ffi(data.kind) {
            FPType::MorganFp => FingerprintKindOptions::Morgan(MorganFingerprintOptions {
                radius: data.radius,
                only_nonzero_invariants: data.only_nonzero_invariants,
                include_redundant_environments: data.include_redundant_environments,
            }),
            FPType::AtomPairFp => FingerprintKindOptions::AtomPair(AtomPairFingerprintOptions {
                min_distance: data.min_distance,
                max_distance: data.max_distance,
                use_2d: data.use_2d,
            }),
            FPType::RdkitFp => FingerprintKindOptions::Rdkit(RDKitFingerprintOptions {
                min_path: data.min_path,
                max_path: data.max_path,
                use_hs: data.use_hs,
                branched_paths: data.branched_paths,
                use_bond_order: data.use_bond_order,
            }),
            FPType::TopologicalTorsionFp => {
                FingerprintKindOptions::TopologicalTorsion(TopologicalTorsionFingerprintOptions {
                    torsion_atom_count: data.torsion_atom_count,
                    only_shortest_paths: data.only_shortest_paths,
                })
            }
        };
        FingerprintOptions {
            count_simulation: data.count_simulation,
            include_chirality: data.include_chirality,
            fp_size: data.fp_size,
            num_bits_per_feature: data.num_bits_per_feature,
            count_bounds: data.count_bounds,
            kind,
        }
    }

    fn to_ffi(&self) -> ffi::FingerprintOptionsData {
        let mut data = ffi::FingerprintOptionsData {
            kind: ffi::FPType::MorganFP,
            count_simulation: self.count_simulation,
            include_chirality: self.include_chirality,
            fp_size: self.fp_size,
            num_bits_per_feature: self.num_bits_per_feature,
            count_bounds: self.count_bounds.clone(),
            radius: 0,
            only_nonzero_invariants: false,
            include_redundant_environments: false,
            min_distance: 0,
            max_distance: 0,
            use_2d: false,
            min_path: 0,
            max_path: 0,
            use_hs: false,
            branched_paths: false,
            use_bond_order: false,
            torsion_atom_count: 0,
            only_shortest_paths: false,
        };
        match &self.kind {
            FingerprintKindOptions::Morgan(m) => {
                data.radius = m.radius;
                data.only_nonzero_invariants = m.only_nonzero_invariants;
                data.include_redundant_environments = m.include_redundant_environments;
            }
            FingerprintKindOptions::AtomPair(a) => {
                data.kind = ffi::FPType::AtomPairFP;
                data.min_distance = a.min_distance;
                data.max_distance = a.max_distance;
                data.use_2d = a.use_2d;
            }
            FingerprintKindOptions::Rdkit(r) => {
                data.kind = ffi::FPType::RDKitFP;
                data.min_path = r.min_path;
                data.max_path = r.max_path;
                data.use_hs = r.use_hs;
                data.branched_paths = r.branched_paths;
                data.use_bond_order = r.use_bond_order;
            }
            FingerprintKindOptions::TopologicalTorsion(t) => {
                data.kind = ffi::FPType::TopologicalTorsionFP;
                data.torsion_atom_count = t.torsion_atom_count;
                data.only_shortest_paths = t.only_shortest_paths;
            }
        }
        data
    }
}

/// Information about how each bit was set. Mirrors
/// `rdFingerprintGenerator.AdditionalOutput`: choose what to collect, pass
/// it in `GetFingerprintParams::additional_output`, and read it from the
/// `_with_outputs` result (A9).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AdditionalOutput {
    atom_to_bits: Option<Vec<Vec<u64>>>,
    bit_info_map: Option<BTreeMap<u64, Vec<(u32, u32)>>>,
    bit_paths: Option<BTreeMap<u64, Vec<Vec<u32>>>>,
    atom_counts: Option<Vec<u32>>,
    atoms_per_bit: Option<BTreeMap<u64, Vec<Vec<u32>>>>,
}

impl AdditionalOutput {
    /// Collects the bits each atom sets. Mirrors
    /// `AdditionalOutput.CollectAtomToBits`.
    pub fn collect_atom_to_bits(&mut self) {
        self.atom_to_bits = Some(Vec::new());
    }

    /// Collects each bit's (atom, radius) pairs for Morgan, or atom pairs for
    /// atom-pair fingerprints. Mirrors `AdditionalOutput.CollectBitInfoMap`.
    pub fn collect_bit_info_map(&mut self) {
        self.bit_info_map = Some(BTreeMap::new());
    }

    /// Collects each bit's bond paths (RDKit fingerprints) or atom paths
    /// (torsions). Mirrors `AdditionalOutput.CollectBitPaths`.
    pub fn collect_bit_paths(&mut self) {
        self.bit_paths = Some(BTreeMap::new());
    }

    /// Collects the number of bits each atom sets. Mirrors
    /// `AdditionalOutput.CollectAtomCounts`.
    pub fn collect_atom_counts(&mut self) {
        self.atom_counts = Some(Vec::new());
    }

    /// Collects the atoms involved in setting each bit. Mirrors
    /// `AdditionalOutput.CollectAtomsPerBit`.
    pub fn collect_atoms_per_bit(&mut self) {
        self.atoms_per_bit = Some(BTreeMap::new());
    }

    /// Mirrors `AdditionalOutput.AllocateAtomToBits`, a synonym of
    /// `collect_atom_to_bits`.
    pub fn allocate_atom_to_bits(&mut self) {
        self.collect_atom_to_bits();
    }

    /// Mirrors `AdditionalOutput.AllocateBitInfoMap`, a synonym of
    /// `collect_bit_info_map`.
    pub fn allocate_bit_info_map(&mut self) {
        self.collect_bit_info_map();
    }

    /// Mirrors `AdditionalOutput.AllocateBitPaths`, a synonym of
    /// `collect_bit_paths`.
    pub fn allocate_bit_paths(&mut self) {
        self.collect_bit_paths();
    }

    /// Mirrors `AdditionalOutput.AllocateAtomCounts`, a synonym of
    /// `collect_atom_counts`.
    pub fn allocate_atom_counts(&mut self) {
        self.collect_atom_counts();
    }

    /// Mirrors `AdditionalOutput.AllocateAtomsPerBit`, a synonym of
    /// `collect_atoms_per_bit`.
    pub fn allocate_atoms_per_bit(&mut self) {
        self.collect_atoms_per_bit();
    }

    /// The bits each atom set, if collected. Mirrors
    /// `AdditionalOutput.GetAtomToBits`.
    pub fn get_atom_to_bits(&self) -> Option<&[Vec<u64>]> {
        self.atom_to_bits.as_deref()
    }

    /// Each bit's (atom, radius) or atom pairs, if collected. Mirrors
    /// `AdditionalOutput.GetBitInfoMap`.
    pub fn get_bit_info_map(&self) -> Option<&BTreeMap<u64, Vec<(u32, u32)>>> {
        self.bit_info_map.as_ref()
    }

    /// Each bit's paths, if collected. Mirrors `AdditionalOutput.GetBitPaths`.
    pub fn get_bit_paths(&self) -> Option<&BTreeMap<u64, Vec<Vec<u32>>>> {
        self.bit_paths.as_ref()
    }

    /// The bits each atom set, counted, if collected. Mirrors
    /// `AdditionalOutput.GetAtomCounts`.
    pub fn get_atom_counts(&self) -> Option<&[u32]> {
        self.atom_counts.as_deref()
    }

    /// The atoms behind each bit, if collected. Mirrors
    /// `AdditionalOutput.GetAtomsPerBit`.
    pub fn get_atoms_per_bit(&self) -> Option<&BTreeMap<u64, Vec<Vec<u32>>>> {
        self.atoms_per_bit.as_ref()
    }

    fn ffi(&self) -> ffi::AdditionalOutputData {
        ffi::AdditionalOutputData {
            collect_atom_to_bits: self.atom_to_bits.is_some(),
            collect_bit_info_map: self.bit_info_map.is_some(),
            collect_bit_paths: self.bit_paths.is_some(),
            collect_atom_counts: self.atom_counts.is_some(),
            collect_atoms_per_bit: self.atoms_per_bit.is_some(),
            atom_to_bits: Vec::new(),
            bit_info_map_keys: Vec::new(),
            bit_info_map_counts: Vec::new(),
            bit_info_map_values: Vec::new(),
            bit_paths_keys: Vec::new(),
            bit_paths_counts: Vec::new(),
            bit_paths_values: Vec::new(),
            atom_counts: Vec::new(),
            atoms_per_bit_keys: Vec::new(),
            atoms_per_bit_counts: Vec::new(),
            atoms_per_bit_values: Vec::new(),
        }
    }

    fn from_ffi(data: ffi::AdditionalOutputData) -> AdditionalOutput {
        fn lists(
            keys: Vec<u64>,
            counts: Vec<u32>,
            values: Vec<IndexList>,
        ) -> BTreeMap<u64, Vec<Vec<u32>>> {
            let mut values = values.into_iter().map(|l| l.indices);
            keys.into_iter()
                .zip(counts)
                .map(|(k, n)| (k, values.by_ref().take(n as usize).collect()))
                .collect()
        }
        let mut pairs = data
            .bit_info_map_values
            .into_iter()
            .map(|p| (p.first, p.second));
        AdditionalOutput {
            atom_to_bits: data
                .collect_atom_to_bits
                .then(|| data.atom_to_bits.into_iter().map(|b| b.values).collect()),
            bit_info_map: data.collect_bit_info_map.then(|| {
                data.bit_info_map_keys
                    .into_iter()
                    .zip(data.bit_info_map_counts)
                    .map(|(k, n)| (k, pairs.by_ref().take(n as usize).collect()))
                    .collect()
            }),
            bit_paths: data.collect_bit_paths.then(|| {
                lists(
                    data.bit_paths_keys,
                    data.bit_paths_counts,
                    data.bit_paths_values,
                )
            }),
            atom_counts: data.collect_atom_counts.then_some(data.atom_counts),
            atoms_per_bit: data.collect_atoms_per_bit.then(|| {
                lists(
                    data.atoms_per_bit_keys,
                    data.atoms_per_bit_counts,
                    data.atoms_per_bit_values,
                )
            }),
        }
    }
}

/// Parameters for the `FingerprintGenerator64::get_*_fingerprint*_with_params`
/// methods, which share Python's parameter list, with Python's defaults (A2).
/// As in Python, an empty list is the same as `None`.
#[derive(Clone, Debug, PartialEq)]
pub struct GetFingerprintParams {
    /// Only environments that start at or are centered on these atoms.
    pub from_atoms: Option<Vec<u32>>,
    /// Leaves out environments that include these atoms.
    pub ignore_atoms: Option<Vec<u32>>,
    /// The conformer, for atom pairs from 3D distances; −1 is the default.
    pub conf_id: i32,
    /// Atom invariants to use instead of the generator's, one per atom.
    pub custom_atom_invariants: Option<Vec<u32>>,
    /// Bond invariants to use instead of the generator's, one per bond.
    pub custom_bond_invariants: Option<Vec<u32>>,
    /// What to collect about the bits (an in-out argument, A9).
    pub additional_output: Option<AdditionalOutput>,
}

impl Default for GetFingerprintParams {
    fn default() -> GetFingerprintParams {
        GetFingerprintParams {
            from_atoms: None,
            ignore_atoms: None,
            conf_id: -1,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            additional_output: None,
        }
    }
}

impl GetFingerprintParams {
    fn to_ffi(&self) -> ffi::FingerprintArgsData {
        fn list(values: &Option<Vec<u32>>) -> (bool, Vec<u32>) {
            (values.is_some(), values.clone().unwrap_or_default())
        }
        let (has_from_atoms, from_atoms) = list(&self.from_atoms);
        let (has_ignore_atoms, ignore_atoms) = list(&self.ignore_atoms);
        let (has_custom_atom_invariants, custom_atom_invariants) =
            list(&self.custom_atom_invariants);
        let (has_custom_bond_invariants, custom_bond_invariants) =
            list(&self.custom_bond_invariants);
        ffi::FingerprintArgsData {
            has_from_atoms,
            from_atoms,
            has_ignore_atoms,
            ignore_atoms,
            conf_id: self.conf_id,
            has_custom_atom_invariants,
            custom_atom_invariants,
            has_custom_bond_invariants,
            custom_bond_invariants,
        }
    }

    fn outputs(&self) -> ffi::AdditionalOutputData {
        self.additional_output.clone().unwrap_or_default().ffi()
    }
}

/// A fingerprint generator. Mirrors
/// `rdFingerprintGenerator.FingerprintGenerator64`.
///
/// The `get_*` calls return `Error::Index` for an atom index past the end,
/// `Error::Value` for custom invariants shorter than the atoms or bonds, and
/// for atom pairs from 3D distances, `Error::Value` for a missing conformer
/// or coordinates beyond 1e9 or not finite (RDKit casts each distance to
/// `unsigned`), and `Error::Value` above 46340 atoms (U3).
pub struct FingerprintGenerator64 {
    /// Never null.
    ptr: UniquePtr<ffi::FpGenerator>,
}

impl fmt::Debug for FingerprintGenerator64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("FingerprintGenerator64")
            .field(&self.get_info_string().unwrap_or_default())
            .finish()
    }
}

macro_rules! fingerprint_methods {
    (
        $kind:ident, $out:ty, $call:ident, $python:literal, $what:literal:
        $name:ident, $with_params:ident, $with_outputs:ident
    ) => {
        #[doc = concat!("The molecule's ", $what, ". Mirrors `FingerprintGenerator64.", $python, "`.")]
        pub fn $name(&self, mol: &Mol) -> Result<$out> {
            self.$with_params(mol, &GetFingerprintParams::default())
        }

        #[doc = concat!(
            "Mirrors `FingerprintGenerator64.", $python,
            "(mol, fromAtoms, ignoreAtoms, confId, customAtomInvariants, customBondInvariants)`."
        )]
        pub fn $with_params(&self, mol: &Mol, params: &GetFingerprintParams) -> Result<$out> {
            Ok(self.$with_outputs(mol, params)?.0)
        }

        #[doc = concat!(
            "The ", $what, " and the `AdditionalOutput` it filled (A9). Mirrors `FingerprintGenerator64.",
            $python, "(..., additionalOutput)`."
        )]
        pub fn $with_outputs(
            &self,
            mol: &Mol,
            params: &GetFingerprintParams,
        ) -> Result<($out, AdditionalOutput)> {
            let mut outputs = params.outputs();
            let fp = ffi::$call(
                &self.ptr,
                mol.as_ffi(),
                ffi::FingerprintOutput::$kind,
                &params.to_ffi(),
                &mut outputs,
            )?;
            Ok((<$out>::from_ffi(fp), AdditionalOutput::from_ffi(outputs)))
        }
    };
}

macro_rules! bulk_methods {
    (
        $kind:ident, $out:ty, $take:ident, $python:literal, $what:literal:
        $name:ident, $with_num_threads:ident
    ) => {
        #[doc = concat!("Each molecule's ", $what, ". Mirrors `FingerprintGenerator64.", $python, "`.")]
        pub fn $name(&self, mols: &[Mol]) -> Result<Vec<$out>> {
            self.$with_num_threads(mols, 1)
        }

        #[doc = concat!(
            "Mirrors `FingerprintGenerator64.", $python, "(mols, numThreads)`; 0 or less uses all ",
            "but that many of the cores. Each molecule is read by one thread."
        )]
        pub fn $with_num_threads(&self, mols: &[Mol], num_threads: i32) -> Result<Vec<$out>> {
            let list = with_borrowed_mols_shared(mols.iter().map(|m| ("", m)), |list| {
                ffi::fp_generator_bulk(&self.ptr, list, ffi::FingerprintOutput::$kind, num_threads)
            })?;
            Ok(take_all(list, |l, i| <$out>::from_ffi(expect_ok(ffi::$take(l, i)))))
        }
    };
}

fn take_all<T>(
    mut list: UniquePtr<ffi::FingerprintList>,
    mut take: impl FnMut(std::pin::Pin<&mut ffi::FingerprintList>, usize) -> T,
) -> Vec<T> {
    let size = expect_ok(ffi::fingerprint_list_size(&list));
    (0..size).map(|i| take(list.pin_mut(), i)).collect()
}

impl FingerprintGenerator64 {
    fn from_ffi(ptr: UniquePtr<ffi::FpGenerator>) -> FingerprintGenerator64 {
        assert!(!ptr.is_null(), "RDKit returned a null generator");
        FingerprintGenerator64 { ptr }
    }

    pub(crate) fn as_ffi(&self) -> &ffi::FpGenerator {
        self.ptr.as_ref().unwrap()
    }

    fingerprint_methods!(
        Bits, ExplicitBitVect, fp_generator_bits, "GetFingerprint", "fingerprint":
        get_fingerprint, get_fingerprint_with_params, get_fingerprint_with_outputs
    );
    fingerprint_methods!(
        SparseBits, SparseBitVect, fp_generator_bits, "GetSparseFingerprint", "unfolded fingerprint":
        get_sparse_fingerprint, get_sparse_fingerprint_with_params, get_sparse_fingerprint_with_outputs
    );
    fingerprint_methods!(
        Counts, UIntSparseIntVect, fp_generator_counts, "GetCountFingerprint", "count fingerprint":
        get_count_fingerprint, get_count_fingerprint_with_params, get_count_fingerprint_with_outputs
    );
    fingerprint_methods!(
        SparseCounts, ULongSparseIntVect, fp_generator_counts, "GetSparseCountFingerprint",
        "unfolded count fingerprint":
        get_sparse_count_fingerprint, get_sparse_count_fingerprint_with_params,
        get_sparse_count_fingerprint_with_outputs
    );

    /// The fingerprint as one 0 or 1 per bit. Mirrors
    /// `FingerprintGenerator64.GetFingerprintAsNumPy`.
    pub fn get_fingerprint_as_num_py(&self, mol: &Mol) -> Result<Vec<u8>> {
        self.get_fingerprint_as_num_py_with_params(mol, &GetFingerprintParams::default())
    }

    /// Mirrors `FingerprintGenerator64.GetFingerprintAsNumPy(mol, fromAtoms,
    /// ignoreAtoms, confId, customAtomInvariants, customBondInvariants,
    /// additionalOutput)`; the additional output is collected and dropped.
    pub fn get_fingerprint_as_num_py_with_params(
        &self,
        mol: &Mol,
        params: &GetFingerprintParams,
    ) -> Result<Vec<u8>> {
        let values = ffi::fp_generator_values(
            &self.ptr,
            mol.as_ffi(),
            ffi::FingerprintOutput::Bits,
            &params.to_ffi(),
            &mut params.outputs(),
        )?;
        Ok(values.into_iter().map(|v| v as u8).collect())
    }

    /// The count fingerprint as one count per bit. Mirrors
    /// `FingerprintGenerator64.GetCountFingerprintAsNumPy`.
    pub fn get_count_fingerprint_as_num_py(&self, mol: &Mol) -> Result<Vec<u32>> {
        self.get_count_fingerprint_as_num_py_with_params(mol, &GetFingerprintParams::default())
    }

    /// Mirrors `FingerprintGenerator64.GetCountFingerprintAsNumPy(mol, ...)`;
    /// the additional output is collected and dropped.
    pub fn get_count_fingerprint_as_num_py_with_params(
        &self,
        mol: &Mol,
        params: &GetFingerprintParams,
    ) -> Result<Vec<u32>> {
        Ok(ffi::fp_generator_values(
            &self.ptr,
            mol.as_ffi(),
            ffi::FingerprintOutput::Counts,
            &params.to_ffi(),
            &mut params.outputs(),
        )?)
    }

    bulk_methods!(
        Bits, ExplicitBitVect, fingerprint_list_take_bits, "GetFingerprints", "fingerprint":
        get_fingerprints, get_fingerprints_with_num_threads
    );
    bulk_methods!(
        SparseBits, SparseBitVect, fingerprint_list_take_bits, "GetSparseFingerprints",
        "unfolded fingerprint":
        get_sparse_fingerprints, get_sparse_fingerprints_with_num_threads
    );
    bulk_methods!(
        Counts, UIntSparseIntVect, fingerprint_list_take_counts, "GetCountFingerprints",
        "count fingerprint":
        get_count_fingerprints, get_count_fingerprints_with_num_threads
    );
    bulk_methods!(
        SparseCounts, ULongSparseIntVect, fingerprint_list_take_counts,
        "GetSparseCountFingerprints", "unfolded count fingerprint":
        get_sparse_count_fingerprints, get_sparse_count_fingerprints_with_num_threads
    );

    /// A description of the generator and its options. Mirrors
    /// `FingerprintGenerator64.GetInfoString`.
    pub fn get_info_string(&self) -> Result<String> {
        to_string(expect_ok(ffi::fp_generator_info_string(&self.ptr)))
    }

    /// A copy of the generator's options. Mirrors
    /// `FingerprintGenerator64.GetOptions`, which returns the options to
    /// change in place; here `set_options` changes them.
    pub fn get_options(&self) -> FingerprintOptions {
        FingerprintOptions::from_ffi(expect_ok(ffi::fp_generator_get_options(&self.ptr)))
    }

    /// Replaces the generator's options with `options`, which must be of the
    /// generator's kind. `Error::Value` for another kind, and for options
    /// the generator could not use safely (see the module documentation).
    pub fn set_options(&mut self, options: &FingerprintOptions) -> Result<()> {
        Ok(ffi::fp_generator_set_options(
            self.ptr.pin_mut(),
            &options.to_ffi(),
        )?)
    }

    /// The generator as JSON. Mirrors `FingerprintGenerator64.ToJSON`.
    pub fn to_json(&self) -> Result<String> {
        to_string(ffi::fp_generator_to_json(&self.ptr)?)
    }
}

/// A generator from `FingerprintGenerator64::to_json`'s JSON. Mirrors
/// `rdFingerprintGenerator.FingerprintGeneratorFromJSON`. `Error::Value` for
/// JSON that leaves out the fingerprint arguments, the environment
/// generator, or the atom-invariants generator, or whose arguments and
/// environment generator are of different kinds (RDKit would dereference
/// null), and for options the generator could not use safely;
/// `Error::Runtime` for text that is not JSON.
pub fn fingerprint_generator_from_json(json_string: &str) -> Result<FingerprintGenerator64> {
    Ok(FingerprintGenerator64::from_ffi(
        ffi::fingerprint_generator_from_json(json_string)?,
    ))
}

fn count_bounds_ffi(bounds: &Option<Vec<u32>>) -> (bool, Vec<u32>) {
    (bounds.is_some(), bounds.clone().unwrap_or_default())
}

/// Parameters for `get_morgan_generator_with_params`, with Python's defaults
/// (A2).
#[derive(Clone, Debug)]
pub struct GetMorganGeneratorParams {
    /// The radius of the environments.
    pub radius: u32,
    /// Simulates counts in bit fingerprints.
    pub count_simulation: bool,
    /// Includes chirality in the atom invariants.
    pub include_chirality: bool,
    /// Includes bond types in the bond invariants.
    pub use_bond_types: bool,
    /// Leaves out atoms whose invariant is zero.
    pub only_nonzero_invariants: bool,
    /// Python accepts and ignores this; pass
    /// `get_morgan_atom_inv_gen(false)` as the atom invariants generator
    /// instead.
    pub include_ring_membership: bool,
    /// The count bounds for count simulation; `None` for 1, 2, 4, 8.
    pub count_bounds: Option<Vec<u32>>,
    /// The number of bits in folded fingerprints.
    pub fp_size: u32,
    /// Atom invariants to use instead of Morgan's; the generator copies it.
    pub atom_invariants_generator: Option<AtomInvariantsGenerator>,
    /// Bond invariants to use instead of Morgan's; the generator copies it.
    pub bond_invariants_generator: Option<BondInvariantsGenerator>,
    /// Keeps environments that repeat a smaller one.
    pub include_redundant_environments: bool,
}

impl Default for GetMorganGeneratorParams {
    fn default() -> GetMorganGeneratorParams {
        GetMorganGeneratorParams {
            radius: 3,
            count_simulation: false,
            include_chirality: false,
            use_bond_types: true,
            only_nonzero_invariants: false,
            include_ring_membership: true,
            count_bounds: None,
            fp_size: 2048,
            atom_invariants_generator: None,
            bond_invariants_generator: None,
            include_redundant_environments: false,
        }
    }
}

/// A Morgan (circular) fingerprint generator. Mirrors
/// `rdFingerprintGenerator.GetMorganGenerator`.
pub fn get_morgan_generator() -> FingerprintGenerator64 {
    expect_ok_generator(get_morgan_generator_with_params(
        &GetMorganGeneratorParams::default(),
    ))
}

/// Mirrors `rdFingerprintGenerator.GetMorganGenerator(radius, ...)`.
/// `Error::Value` for an `fp_size` above `i32::MAX`, or count simulation
/// without count bounds.
pub fn get_morgan_generator_with_params(
    params: &GetMorganGeneratorParams,
) -> Result<FingerprintGenerator64> {
    let (has_count_bounds, count_bounds) = count_bounds_ffi(&params.count_bounds);
    let data = ffi::MorganGeneratorData {
        radius: params.radius,
        count_simulation: params.count_simulation,
        include_chirality: params.include_chirality,
        use_bond_types: params.use_bond_types,
        only_nonzero_invariants: params.only_nonzero_invariants,
        include_ring_membership: params.include_ring_membership,
        has_count_bounds,
        count_bounds,
        fp_size: params.fp_size,
        include_redundant_environments: params.include_redundant_environments,
    };
    Ok(FingerprintGenerator64::from_ffi(ffi::get_morgan_generator(
        &data,
        atom_inv_gen_ffi(&params.atom_invariants_generator),
        params
            .bond_invariants_generator
            .as_ref()
            .map_or_else(UniquePtr::null, BondInvariantsGenerator::clone_ffi),
    )?))
}

/// Parameters for `get_atom_pair_generator_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug)]
pub struct GetAtomPairGeneratorParams {
    /// The shortest distance between a pair's atoms.
    pub min_distance: u32,
    /// The longest distance between a pair's atoms.
    pub max_distance: u32,
    /// Includes chirality in the atom invariants.
    pub include_chirality: bool,
    /// Uses topological distances; otherwise 3D ones.
    pub use_2d: bool,
    /// Simulates counts in bit fingerprints.
    pub count_simulation: bool,
    /// The count bounds for count simulation; `None` for 1, 2, 4, 8.
    pub count_bounds: Option<Vec<u32>>,
    /// The number of bits in folded fingerprints.
    pub fp_size: u32,
    /// Atom invariants to use instead of the atom-pair ones; the generator
    /// copies it.
    pub atom_invariants_generator: Option<AtomInvariantsGenerator>,
}

impl Default for GetAtomPairGeneratorParams {
    fn default() -> GetAtomPairGeneratorParams {
        GetAtomPairGeneratorParams {
            min_distance: 1,
            max_distance: 30,
            include_chirality: false,
            use_2d: true,
            count_simulation: true,
            count_bounds: None,
            fp_size: 2048,
            atom_invariants_generator: None,
        }
    }
}

/// An atom-pair fingerprint generator. Mirrors
/// `rdFingerprintGenerator.GetAtomPairGenerator`.
pub fn get_atom_pair_generator() -> FingerprintGenerator64 {
    expect_ok_generator(get_atom_pair_generator_with_params(
        &GetAtomPairGeneratorParams::default(),
    ))
}

/// Mirrors `rdFingerprintGenerator.GetAtomPairGenerator(minDistance, ...)`.
/// `Error::Value` as `get_morgan_generator_with_params`.
pub fn get_atom_pair_generator_with_params(
    params: &GetAtomPairGeneratorParams,
) -> Result<FingerprintGenerator64> {
    let (has_count_bounds, count_bounds) = count_bounds_ffi(&params.count_bounds);
    let data = ffi::AtomPairGeneratorData {
        min_distance: params.min_distance,
        max_distance: params.max_distance,
        include_chirality: params.include_chirality,
        use_2d: params.use_2d,
        count_simulation: params.count_simulation,
        has_count_bounds,
        count_bounds,
        fp_size: params.fp_size,
    };
    Ok(FingerprintGenerator64::from_ffi(
        ffi::get_atom_pair_generator(&data, atom_inv_gen_ffi(&params.atom_invariants_generator))?,
    ))
}

/// Parameters for `get_rdkit_fp_generator_with_params`, with Python's
/// defaults (A2).
#[derive(Clone, Debug)]
pub struct GetRDKitFPGeneratorParams {
    /// The fewest bonds in a path.
    pub min_path: u32,
    /// The most bonds in a path.
    pub max_path: u32,
    /// Includes explicit hydrogens in the paths.
    pub use_hs: bool,
    /// Uses branched subgraphs, not only linear paths.
    pub branched_paths: bool,
    /// Includes bond orders in the path hashes.
    pub use_bond_order: bool,
    /// Simulates counts in bit fingerprints.
    pub count_simulation: bool,
    /// The count bounds for count simulation; `None` for 1, 2, 4, 8.
    pub count_bounds: Option<Vec<u32>>,
    /// The number of bits in folded fingerprints.
    pub fp_size: u32,
    /// The bits each path sets.
    pub num_bits_per_feature: u32,
    /// Atom invariants to use instead of the RDKit ones; the generator
    /// copies it.
    pub atom_invariants_generator: Option<AtomInvariantsGenerator>,
}

impl Default for GetRDKitFPGeneratorParams {
    fn default() -> GetRDKitFPGeneratorParams {
        GetRDKitFPGeneratorParams {
            min_path: 1,
            max_path: 7,
            use_hs: true,
            branched_paths: true,
            use_bond_order: true,
            count_simulation: false,
            count_bounds: None,
            fp_size: 2048,
            num_bits_per_feature: 2,
            atom_invariants_generator: None,
        }
    }
}

/// An RDKit (path-based) fingerprint generator. Mirrors
/// `rdFingerprintGenerator.GetRDKitFPGenerator`.
pub fn get_rdkit_fp_generator() -> FingerprintGenerator64 {
    expect_ok_generator(get_rdkit_fp_generator_with_params(
        &GetRDKitFPGeneratorParams::default(),
    ))
}

/// Mirrors `rdFingerprintGenerator.GetRDKitFPGenerator(minPath, ...)`.
/// `Error::Value` as `get_morgan_generator_with_params`; `Error::Runtime`
/// for `num_bits_per_feature` 0.
pub fn get_rdkit_fp_generator_with_params(
    params: &GetRDKitFPGeneratorParams,
) -> Result<FingerprintGenerator64> {
    let (has_count_bounds, count_bounds) = count_bounds_ffi(&params.count_bounds);
    let data = ffi::RDKitFPGeneratorData {
        min_path: params.min_path,
        max_path: params.max_path,
        use_hs: params.use_hs,
        branched_paths: params.branched_paths,
        use_bond_order: params.use_bond_order,
        count_simulation: params.count_simulation,
        has_count_bounds,
        count_bounds,
        fp_size: params.fp_size,
        num_bits_per_feature: params.num_bits_per_feature,
    };
    Ok(FingerprintGenerator64::from_ffi(
        ffi::get_rdkit_fp_generator(&data, atom_inv_gen_ffi(&params.atom_invariants_generator))?,
    ))
}

/// Parameters for `get_topological_torsion_generator_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug)]
pub struct GetTopologicalTorsionGeneratorParams {
    /// Includes chirality in the atom invariants.
    pub include_chirality: bool,
    /// The number of atoms in each torsion.
    pub torsion_atom_count: u32,
    /// Simulates counts in bit fingerprints.
    pub count_simulation: bool,
    /// The count bounds for count simulation; `None` for 1, 2, 4, 8.
    pub count_bounds: Option<Vec<u32>>,
    /// The number of bits in folded fingerprints.
    pub fp_size: u32,
    /// Atom invariants to use instead of the torsion ones; the generator
    /// copies it.
    pub atom_invariants_generator: Option<AtomInvariantsGenerator>,
}

impl Default for GetTopologicalTorsionGeneratorParams {
    fn default() -> GetTopologicalTorsionGeneratorParams {
        GetTopologicalTorsionGeneratorParams {
            include_chirality: false,
            torsion_atom_count: 4,
            count_simulation: true,
            count_bounds: None,
            fp_size: 2048,
            atom_invariants_generator: None,
        }
    }
}

/// A topological-torsion fingerprint generator. Mirrors
/// `rdFingerprintGenerator.GetTopologicalTorsionGenerator`.
pub fn get_topological_torsion_generator() -> FingerprintGenerator64 {
    expect_ok_generator(get_topological_torsion_generator_with_params(
        &GetTopologicalTorsionGeneratorParams::default(),
    ))
}

/// Mirrors `rdFingerprintGenerator.GetTopologicalTorsionGenerator(
/// includeChirality, torsionAtomCount, ...)`. `Error::Value` as
/// `get_morgan_generator_with_params`, and for a `torsion_atom_count` above
/// 7 (5 with chirality), whose codes RDKit would shift past 64 bits.
pub fn get_topological_torsion_generator_with_params(
    params: &GetTopologicalTorsionGeneratorParams,
) -> Result<FingerprintGenerator64> {
    let (has_count_bounds, count_bounds) = count_bounds_ffi(&params.count_bounds);
    let data = ffi::TopologicalTorsionGeneratorData {
        include_chirality: params.include_chirality,
        torsion_atom_count: params.torsion_atom_count,
        count_simulation: params.count_simulation,
        has_count_bounds,
        count_bounds,
        fp_size: params.fp_size,
    };
    Ok(FingerprintGenerator64::from_ffi(
        ffi::get_topological_torsion_generator(
            &data,
            atom_inv_gen_ffi(&params.atom_invariants_generator),
        )?,
    ))
}

fn expect_ok_generator(result: Result<FingerprintGenerator64>) -> FingerprintGenerator64 {
    result.unwrap_or_else(|err| panic!("the default generator failed: {err}"))
}

macro_rules! type_bulk {
    ($kind:ident, $out:ty, $take:ident, $python:literal, $what:literal: $name:ident, $with_fp_type:ident) => {
        #[doc = concat!(
            "Each molecule's Morgan ", $what, ", from a default radius-2 generator. Mirrors `rdFingerprintGenerator.",
            $python, "`; Python's `molecules=None` is an empty slice."
        )]
        pub fn $name(molecules: &[Mol]) -> Result<Vec<$out>> {
            $with_fp_type(molecules, FPType::MorganFp)
        }

        #[doc = concat!(
            "Each molecule's ", $what, " from a default generator of the type (Morgan's radius is 2). ",
            "Mirrors `rdFingerprintGenerator.", $python, "(molecules, fpType)`."
        )]
        pub fn $with_fp_type(molecules: &[Mol], fp_type: FPType) -> Result<Vec<$out>> {
            let list = with_borrowed_mols_shared(molecules.iter().map(|m| ("", m)), |list| {
                ffi::fp_type_bulk(list, fp_type.to_ffi(), ffi::FingerprintOutput::$kind)
            })?;
            Ok(take_all(list, |l, i| <$out>::from_ffi(expect_ok(ffi::$take(l, i)))))
        }
    };
}

type_bulk!(Bits, ExplicitBitVect, fingerprint_list_take_bits, "GetFPs", "fingerprint": get_fps, get_fps_with_fp_type);
type_bulk!(
    SparseBits, SparseBitVect, fingerprint_list_take_bits, "GetSparseFPs", "unfolded fingerprint":
    get_sparse_fps, get_sparse_fps_with_fp_type
);
type_bulk!(
    Counts, UIntSparseIntVect, fingerprint_list_take_counts, "GetCountFPs", "count fingerprint":
    get_count_fps, get_count_fps_with_fp_type
);
type_bulk!(
    SparseCounts, ULongSparseIntVect, fingerprint_list_take_counts, "GetSparseCountFPs",
    "unfolded count fingerprint":
    get_sparse_count_fps, get_sparse_count_fps_with_fp_type
);
