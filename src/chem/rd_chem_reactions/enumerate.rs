// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Enumerating a reaction's products over sets of building blocks.

use std::fmt;
use std::ops::{Deref, DerefMut};
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::rd_chem_reactions::ChemicalReaction;
use crate::chem::rd_chem_reactions::reaction::product_sets_into_vec;
use crate::chem::rdmolops::with_borrowed_mols_shared;
use crate::chem::{Mol, mol_to_smiles};
use crate::error::{expect_ok, to_string};
use crate::ffi::enumerate::ffi;

/// Calls `f` with the building blocks lent as one list and their groups'
/// sizes.
fn with_building_blocks<T>(
    groups: &[Vec<Mol>],
    f: impl FnOnce(Pin<&mut crate::ffi::molops::ffi::BorrowedMols>, &[u32]) -> T,
) -> T {
    let sizes: Vec<u32> = groups
        .iter()
        .map(|group| u32::try_from(group.len()).expect("a group has more than u32::MAX molecules"))
        .collect();
    with_borrowed_mols_shared(groups.iter().flatten().map(|mol| ("", mol)), |list| {
        f(list, &sizes)
    })
}

/// A way to walk the combinations of building blocks, one index per group.
/// Mirrors `rdChemReactions.EnumerationStrategyBase`, which holds any of the
/// strategies; each strategy type dereferences to it.
pub struct EnumerationStrategyBase {
    /// Never null.
    ptr: UniquePtr<ffi::EnumerationStrategyBase>,
}

impl EnumerationStrategyBase {
    fn from_ffi(ptr: UniquePtr<ffi::EnumerationStrategyBase>) -> EnumerationStrategyBase {
        assert!(!ptr.is_null(), "RDKit returned a null strategy");
        EnumerationStrategyBase { ptr }
    }

    fn new(kind: u8) -> EnumerationStrategyBase {
        EnumerationStrategyBase::from_ffi(expect_ok(ffi::enumeration_strategy_new(kind)))
    }

    fn raw(&self) -> &ffi::EnumerationStrategyBase {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::EnumerationStrategyBase> {
        self.ptr.pin_mut()
    }

    /// The number of combinations, `u64::MAX` when there are more than a
    /// `u64` counts. Mirrors `EnumerationStrategyBase.GetNumPermutations`.
    pub fn get_num_permutations(&self) -> u64 {
        expect_ok(ffi::enumeration_strategy_get_num_permutations(self.raw()))
    }

    /// The current indices into the groups. Mirrors
    /// `EnumerationStrategyBase.GetPosition`.
    pub fn get_position(&self) -> Vec<u64> {
        expect_ok(ffi::enumeration_strategy_get_position(self.raw()))
    }

    /// Starts over for groups of building blocks, of which only the sizes
    /// matter. The random strategies need between 1 and 2^31 building
    /// blocks in each group (and `RandomSampleAllBBsStrategy` a group), and
    /// `EvenSamplePairsStrategy` at most 2^63 combinations; others are
    /// `Error::Value`. Mirrors `EnumerationStrategyBase.Initialize`.
    pub fn initialize(&mut self, rxn: &ChemicalReaction, ob: &[Vec<Mol>]) -> Result<()> {
        Ok(with_building_blocks(ob, |list, sizes| {
            ffi::enumeration_strategy_initialize(self.raw_mut(), rxn.as_ffi(), list, sizes)
        })?)
    }

    /// Steps past `skip_count` combinations, one at a time. Mirrors
    /// `EnumerationStrategyBase.Skip`.
    pub fn skip(&mut self, skip_count: u64) -> Result<bool> {
        Ok(ffi::enumeration_strategy_skip(self.raw_mut(), skip_count)?)
    }

    /// The strategy's class name. Mirrors `EnumerationStrategyBase.Type`.
    pub fn type_(&self) -> String {
        to_string(expect_ok(ffi::enumeration_strategy_type(self.raw())))
            .expect("a strategy's name is ASCII")
    }

    /// Whether combinations remain; the random strategies never run out.
    /// Mirrors `EnumerationStrategyBase.__bool__`.
    pub fn has_next(&self) -> bool {
        expect_ok(ffi::enumeration_strategy_has_next(self.raw()))
    }

    /// Steps to the next combination and returns its indices.
    /// `EvenSamplePairsStrategy` past its last is `Error::Runtime`. Mirrors
    /// `EnumerationStrategyBase.next`.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Vec<u64>> {
        Ok(ffi::enumeration_strategy_next(self.raw_mut())?)
    }
}

/// A copy, state included. Mirrors `EnumerationStrategyBase.__copy__`.
impl Clone for EnumerationStrategyBase {
    fn clone(&self) -> EnumerationStrategyBase {
        EnumerationStrategyBase::from_ffi(expect_ok(ffi::enumeration_strategy_clone(self.raw())))
    }
}

impl fmt::Debug for EnumerationStrategyBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EnumerationStrategyBase")
            .field("type", &self.type_())
            .field("position", &self.get_position())
            .finish()
    }
}

macro_rules! strategy {
    ($(#[$meta:meta])* $name:ident, $kind:expr) => {
        $(#[$meta])*
        #[derive(Clone, Debug)]
        pub struct $name {
            base: EnumerationStrategyBase,
        }

        #[doc = concat!("A new strategy. Mirrors `", stringify!($name), "()`.")]
        impl Default for $name {
            fn default() -> $name {
                $name {
                    base: EnumerationStrategyBase::new($kind),
                }
            }
        }

        impl Deref for $name {
            type Target = EnumerationStrategyBase;

            fn deref(&self) -> &EnumerationStrategyBase {
                &self.base
            }
        }

        impl DerefMut for $name {
            fn deref_mut(&mut self) -> &mut EnumerationStrategyBase {
                &mut self.base
            }
        }

        impl From<$name> for EnumerationStrategyBase {
            fn from(strategy: $name) -> EnumerationStrategyBase {
                strategy.base
            }
        }
    };
}

strategy!(
    /// Every combination in turn, the first group's index fastest:
    /// (0,0,0), (1,0,0), (2,0,0), ... Mirrors
    /// `rdChemReactions.CartesianProductStrategy`; `clone` mirrors its
    /// `__copy__`.
    CartesianProductStrategy,
    0
);
strategy!(
    /// Combinations drawn at random, without end and with repeats. Mirrors
    /// `rdChemReactions.RandomSampleStrategy`; `clone` mirrors its `__copy__`.
    RandomSampleStrategy,
    1
);
strategy!(
    /// Combinations drawn at random so that every building block is used
    /// early, without end. Mirrors `rdChemReactions.RandomSampleAllBBsStrategy`;
    /// `clone` mirrors its `__copy__`.
    RandomSampleAllBBsStrategy,
    2
);
strategy!(
    /// Combinations drawn so that pairs of building blocks are used evenly.
    /// Mirrors `rdChemReactions.EvenSamplePairsStrategy`; `clone` mirrors its
    /// `__copy__`.
    EvenSamplePairsStrategy,
    3
);

impl EvenSamplePairsStrategy {
    /// How often each building block and pair has been used. Mirrors
    /// `EvenSamplePairsStrategy.Stats`.
    ///
    /// # Panics
    ///
    /// If the strategy behind it was swapped for another kind through
    /// `DerefMut`.
    pub fn stats(&self) -> String {
        to_string(expect_ok(ffi::even_sample_pairs_strategy_stats(
            self.base.raw(),
        )))
        .expect("the statistics are ASCII")
    }
}

/// How `EnumerateLibrary` filters building blocks. Mirrors
/// `rdChemReactions.EnumerationParams`.
#[derive(Clone, Debug, PartialEq)]
pub struct EnumerationParams {
    /// Drops building blocks a reactant template matches more often.
    pub reagent_max_match_count: i32,
    /// Drops building blocks whose products with the product templates fail
    /// sanitization.
    pub sane_partial_products: bool,
}

/// RDKit's defaults. Mirrors `EnumerationParams()`.
impl Default for EnumerationParams {
    fn default() -> EnumerationParams {
        let data = expect_ok(ffi::enumeration_params_default());
        EnumerationParams {
            reagent_max_match_count: data.reagent_max_match_count,
            sane_partial_products: data.sane_partial_products,
        }
    }
}

impl EnumerationParams {
    fn to_ffi(&self) -> ffi::EnumerationParamsData {
        ffi::EnumerationParamsData {
            reagent_max_match_count: self.reagent_max_match_count,
            sane_partial_products: self.sane_partial_products,
        }
    }
}

/// A reaction's products over groups of building blocks, one group per
/// reactant template, in the order a strategy picks them. Mirrors
/// `rdChemReactions.EnumerateLibrary` and its base class,
/// `EnumerateLibraryBase`.
///
/// The library owns copies of its building blocks; Python's shares Python's
/// molecules without keeping them alive. It is an iterator over product
/// sets, as in Python; `has_next` mirrors its `__bool__`.
///
/// `EnumerateLibrary` is `Send` but not `Sync` (CONVENTIONS.md S4).
///
/// ```no_run
/// use rdkit_cxx::chem::mol_from_smiles;
/// use rdkit_cxx::chem::rd_chem_reactions::{EnumerateLibrary, reaction_from_smarts};
///
/// let rxn = reaction_from_smarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]")?;
/// let acids = vec![mol_from_smiles("CC(=O)O")?, mol_from_smiles("CCC(=O)O")?];
/// let amines = vec![mol_from_smiles("CN")?];
/// for products in EnumerateLibrary::new(&rxn, &[acids, amines])? {
///     let products = products?;
/// }
/// # Ok::<(), rdkit_cxx::Error>(())
/// ```
pub struct EnumerateLibrary {
    /// Never null.
    ptr: UniquePtr<ffi::EnumerateLibrary>,
}

impl EnumerateLibrary {
    fn from_ffi(ptr: UniquePtr<ffi::EnumerateLibrary>) -> EnumerateLibrary {
        assert!(!ptr.is_null(), "RDKit returned a null library");
        EnumerateLibrary { ptr }
    }

    fn raw(&self) -> &ffi::EnumerateLibrary {
        self.ptr.as_ref().unwrap()
    }

    fn raw_mut(&mut self) -> Pin<&mut ffi::EnumerateLibrary> {
        self.ptr.pin_mut()
    }

    /// Every combination of the building blocks each reactant template
    /// matches. Mirrors `EnumerateLibrary(rxn, reagents)`.
    pub fn new(rxn: &ChemicalReaction, reagents: &[Vec<Mol>]) -> Result<EnumerateLibrary> {
        EnumerateLibrary::new_with_params(rxn, reagents, &EnumerationParams::default())
    }

    /// Mirrors `EnumerateLibrary(rxn, reagents, params)`.
    pub fn new_with_params(
        rxn: &ChemicalReaction,
        reagents: &[Vec<Mol>],
        params: &EnumerationParams,
    ) -> Result<EnumerateLibrary> {
        let ptr = with_building_blocks(reagents, |list, sizes| {
            ffi::enumerate_library_new_with_params(rxn.as_ffi(), list, sizes, &params.to_ffi())
        })?;
        Ok(EnumerateLibrary::from_ffi(ptr))
    }

    /// The combinations a copy of `enumerator` picks, which must suit the
    /// building blocks the templates match (see
    /// [`EnumerationStrategyBase::initialize`]). Mirrors
    /// `EnumerateLibrary(rxn, reagents, enumerator)`.
    pub fn new_with_enumerator(
        rxn: &ChemicalReaction,
        reagents: &[Vec<Mol>],
        enumerator: &EnumerationStrategyBase,
    ) -> Result<EnumerateLibrary> {
        EnumerateLibrary::new_with_enumerator_and_params(
            rxn,
            reagents,
            enumerator,
            &EnumerationParams::default(),
        )
    }

    /// Mirrors `EnumerateLibrary(rxn, reagents, enumerator, params)`.
    pub fn new_with_enumerator_and_params(
        rxn: &ChemicalReaction,
        reagents: &[Vec<Mol>],
        enumerator: &EnumerationStrategyBase,
        params: &EnumerationParams,
    ) -> Result<EnumerateLibrary> {
        let ptr = with_building_blocks(reagents, |list, sizes| {
            ffi::enumerate_library_new_with_enumerator(
                rxn.as_ffi(),
                list,
                sizes,
                enumerator.raw(),
                &params.to_ffi(),
            )
        })?;
        Ok(EnumerateLibrary::from_ffi(ptr))
    }

    /// Copies of the building blocks the templates match, by group. Mirrors
    /// `EnumerateLibrary.GetReagents`.
    pub fn get_reagents(&self) -> Vec<Vec<Mol>> {
        product_sets_into_vec(expect_ok(ffi::enumerate_library_get_reagents(self.raw())))
    }

    /// A copy of the library's reaction. Mirrors
    /// `EnumerateLibraryBase.GetReaction`, which Python returns as a view.
    pub fn get_reaction(&self) -> ChemicalReaction {
        ChemicalReaction::from_ffi(expect_ok(ffi::enumerate_library_get_reaction(self.raw())))
    }

    /// A copy of the strategy, at its position. A library made by
    /// `default` has none: `Error::Value`. Mirrors
    /// `EnumerateLibraryBase.GetEnumerator`, which Python returns as a view.
    pub fn get_enumerator(&self) -> Result<EnumerationStrategyBase> {
        Ok(EnumerationStrategyBase::from_ffi(
            ffi::enumerate_library_get_enumerator(self.raw())?,
        ))
    }

    /// The strategy's indices into `get_reagents`. A library without a
    /// strategy is `Error::Value` (Python segfaults). Mirrors
    /// `EnumerateLibraryBase.GetPosition`.
    pub fn get_position(&self) -> Result<Vec<u64>> {
        Ok(ffi::enumerate_library_get_position(self.raw())?)
    }

    /// Whether products remain. A library without a strategy is
    /// `Error::Value`. Mirrors `EnumerateLibraryBase.__bool__`.
    pub fn has_next(&self) -> Result<bool> {
        Ok(ffi::enumerate_library_has_next(self.raw())?)
    }

    fn next_products(&mut self) -> Result<Vec<Vec<Mol>>> {
        Ok(product_sets_into_vec(ffi::enumerate_library_next(
            self.raw_mut(),
        )?))
    }

    /// The next product sets as SMILES. A library with no products left is
    /// `Error::Runtime`. Mirrors `EnumerateLibraryBase.nextSmiles`.
    pub fn next_smiles(&mut self) -> Result<Vec<Vec<String>>> {
        self.next_products()?
            .iter()
            .map(|set| set.iter().map(mol_to_smiles).collect())
            .collect()
    }

    /// The strategy's state, for `set_state`. Mirrors
    /// `EnumerateLibraryBase.GetState`.
    pub fn get_state(&self) -> Result<String> {
        to_string(ffi::enumerate_library_get_state(self.raw())?)
    }

    /// Restores a state from `get_state`. A state that does not fit the
    /// library's building blocks (one from another library) is
    /// `Error::Value`; Python's crashes on one from a larger library.
    /// Mirrors `EnumerateLibraryBase.SetState`.
    ///
    /// # Safety
    ///
    /// `state` must be text `get_state` (this crate's or RDKit's) returned.
    /// RDKit reads it as a Boost.Serialization archive, which trusts its
    /// class and object ids; other text may make it read out of bounds
    /// (CONVENTIONS.md U2).
    pub unsafe fn set_state(&mut self, state: &str) -> Result<()> {
        Ok(ffi::enumerate_library_set_state(self.raw_mut(), state)?)
    }

    /// Returns the strategy to where it started. Mirrors
    /// `EnumerateLibraryBase.ResetState`.
    pub fn reset_state(&mut self) -> Result<()> {
        Ok(ffi::enumerate_library_reset_state(self.raw_mut())?)
    }

    /// The library as bytes, position included. Mirrors
    /// `EnumerateLibraryBase.Serialize`.
    pub fn serialize(&self) -> Result<Vec<u8>> {
        Ok(ffi::enumerate_library_serialize(self.raw())?)
    }

    /// Fills a library made by `default` from `serialize`'s bytes; a library
    /// with a reaction is `Error::Value`, as is a strategy state that does
    /// not fit the building blocks, or a molecule whose plain monomer info
    /// has the PDB residue type (CONVENTIONS.md U10). Mirrors
    /// `EnumerateLibraryBase.InitFromString`.
    ///
    /// # Safety
    ///
    /// `data` must be bytes `serialize` (this crate's or RDKit's) returned,
    /// for the reasons `set_state` gives.
    pub unsafe fn init_from_string(&mut self, data: &[u8]) -> Result<()> {
        if !expect_ok(ffi::enumerate_library_is_empty(self.raw())) {
            return Err(crate::Error::Value {
                message: "EnumerateBase already created from archive.".to_owned(),
            });
        }
        self.ptr = ffi::enumerate_library_from_string(data)?;
        Ok(())
    }
}

/// An empty library, for `init_from_string`. Mirrors `EnumerateLibrary()`.
impl Default for EnumerateLibrary {
    fn default() -> EnumerateLibrary {
        EnumerateLibrary::from_ffi(expect_ok(ffi::enumerate_library_new()))
    }
}

/// The product sets, one list per product template, for each combination
/// in turn, until the strategy runs out (never, for the random strategies).
/// Mirrors `EnumerateLibraryBase.next` and `__next__`.
impl Iterator for EnumerateLibrary {
    type Item = Result<Vec<Vec<Mol>>>;

    fn next(&mut self) -> Option<Result<Vec<Vec<Mol>>>> {
        match self.has_next() {
            Ok(true) => Some(self.next_products()),
            _ => None,
        }
    }
}

impl fmt::Debug for EnumerateLibrary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EnumerateLibrary")
            .field("position", &self.get_position().ok())
            .finish()
    }
}

/// Whether libraries can be serialized (RDKit built with Boost
/// serialization). Mirrors `rdChemReactions.EnumerateLibraryCanSerialize`.
pub fn enumerate_library_can_serialize() -> bool {
    expect_ok(ffi::enumerate_library_can_serialize())
}
