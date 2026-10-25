// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The `ChemicalReaction` class and its template views.

use std::collections::BTreeMap;
use std::fmt;
use std::marker::PhantomData;
use std::mem::ManuallyDrop;
use std::ops::{Deref, Range};
use std::pin::Pin;

use cxx::UniquePtr;

use crate::Result;
use crate::chem::Mol;
use crate::chem::rdchem::props::{prop_getters, prop_setters};
use crate::chem::rdchem::{PropertyPickleOptions, SubstructMatchParameters};
use crate::chem::rdmolops::{mol_list_into_vec, with_borrowed_mols_shared};
use crate::error::{expect_ok, to_string};
use crate::ffi::chem_reactions::ffi;
use crate::ffi::mol::ffi as mol_ffi;
use crate::ffi::props::ffi::RDProps;

pub(crate) const REACTANT: u8 = 0;
pub(crate) const PRODUCT: u8 = 1;
pub(crate) const AGENT: u8 = 2;

/// Each reactant template's labels: the atom index and the query's name.
/// Mirrors the tuples `AddRecursiveQueriesToReaction(getLabels=True)` and
/// `PreprocessReaction` return.
pub type ReactantLabels = Vec<Vec<(u32, String)>>;

/// A chemical reaction: reactant, product, and agent templates, with
/// properties. Mirrors `rdChemReactions.ChemicalReaction`.
///
/// The reaction owns its templates. The methods that add templates copy the
/// molecule, where Python's reaction keeps Python's molecule, and the
/// template getters lend views that borrow the reaction (CONVENTIONS.md W2).
///
/// Running a reaction takes `&mut self`: like Python's runs, these
/// initialize the reaction first when it is not, and initialization
/// validates the reaction, which writes properties on its product templates
/// (CONVENTIONS.md S4).
///
/// `ChemicalReaction` is `Send` but not `Sync`, like `Mol`: RDKit writes
/// through a read-only reaction too (the SMILES writers set computed
/// properties on its templates), so parallel code gives each thread its own
/// reaction (S1).
///
/// ```compile_fail
/// fn assert_sync<T: Sync>() {}
/// assert_sync::<rdkit_cxx::chem::rd_chem_reactions::ChemicalReaction>();
/// ```
pub struct ChemicalReaction {
    /// Never null.
    ptr: UniquePtr<ffi::ChemicalReaction>,
}

/// A template of a reaction, borrowed from it: one of the molecules
/// `ChemicalReaction.GetReactantTemplate`, `GetProductTemplate`, and
/// `GetAgentTemplate` return, which Python also lends as views into the
/// reaction. It dereferences to [`Mol`] for reading; it is not `Send`, since
/// it shares the reaction.
pub struct TemplateRef<'r> {
    /// A molecule over the reaction's template, never dropped, since the
    /// reaction owns the template.
    mol: ManuallyDrop<Mol>,
    reaction: PhantomData<&'r ChemicalReaction>,
}

impl<'r> TemplateRef<'r> {
    fn new(template: &'r mol_ffi::RWMol) -> TemplateRef<'r> {
        let raw = template as *const mol_ffi::RWMol as *mut mol_ffi::RWMol;
        // SAFETY: `raw` points to an RWMol the reaction owns for 'r (every
        // template is one, normalize_templates in chem_reactions.cpp). The
        // ManuallyDrop never frees it, and TemplateRef lends only `&Mol`,
        // through which nothing can move the pointer out or edit the template
        // (Mol is not Copy, and Clone copies the molecule).
        let ptr = unsafe { UniquePtr::from_raw(raw) };
        TemplateRef {
            mol: ManuallyDrop::new(Mol::from_ffi(ptr)),
            reaction: PhantomData,
        }
    }
}

impl Deref for TemplateRef<'_> {
    type Target = Mol;

    fn deref(&self) -> &Mol {
        &self.mol
    }
}

impl fmt::Debug for TemplateRef<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("TemplateRef").field(&*self.mol).finish()
    }
}

/// The templates of one kind, as views. Mirrors the lists
/// `ChemicalReaction.GetReactants`, `GetProducts`, and `GetAgents` return
/// (CONVENTIONS.md W5).
#[derive(Clone, Debug)]
pub struct Templates<'r> {
    rxn: &'r ChemicalReaction,
    kind: u8,
    range: Range<u32>,
}

impl<'r> Iterator for Templates<'r> {
    type Item = TemplateRef<'r>;

    fn next(&mut self) -> Option<TemplateRef<'r>> {
        let idx = self.range.next()?;
        Some(self.rxn.template(self.kind, idx))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.range.size_hint()
    }
}

impl DoubleEndedIterator for Templates<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let idx = self.range.next_back()?;
        Some(self.rxn.template(self.kind, idx))
    }
}

impl ExactSizeIterator for Templates<'_> {}

/// Parameters for `remove_unmapped_reactant_templates_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct RemoveUnmappedReactantTemplatesParams {
    /// Templates whose fraction of mapped heavy atoms is below this are
    /// removed.
    pub threshold_unmapped_atoms: f64,
    /// Moves the removed templates to the agents.
    pub move_to_agent_templates: bool,
}

impl Default for RemoveUnmappedReactantTemplatesParams {
    fn default() -> RemoveUnmappedReactantTemplatesParams {
        RemoveUnmappedReactantTemplatesParams {
            threshold_unmapped_atoms: 0.2,
            move_to_agent_templates: true,
        }
    }
}

/// Parameters for `remove_unmapped_product_templates_with_params`, with
/// Python's defaults (A2): those of the reactant form.
pub type RemoveUnmappedProductTemplatesParams = RemoveUnmappedReactantTemplatesParams;

/// Parameters for `add_recursive_queries_to_reaction_with_params`, with
/// Python's defaults (A2).
#[derive(Clone, Debug)]
pub struct AddRecursiveQueriesToReactionParams {
    /// The queries, by the names the atoms' labels use.
    pub queries: BTreeMap<String, Mol>,
    /// The atom property holding the labels.
    pub prop_name: String,
}

impl Default for AddRecursiveQueriesToReactionParams {
    fn default() -> AddRecursiveQueriesToReactionParams {
        AddRecursiveQueriesToReactionParams {
            queries: BTreeMap::new(),
            prop_name: "molFileValue".to_owned(),
        }
    }
}

/// Groups the labels by reactant template.
pub(crate) fn group_labels(
    labels: Vec<ffi::ReactantLabel>,
    num_lists: u32,
) -> Result<ReactantLabels> {
    let mut res = vec![Vec::new(); num_lists as usize];
    for label in labels {
        res[label.reactant as usize].push((label.atom_idx, to_string(label.label)?));
    }
    Ok(res)
}

/// Moves RDKit's product sets out, one list per set (W6).
pub(crate) fn product_sets_into_vec(mut sets: UniquePtr<ffi::ProductSets>) -> Vec<Vec<Mol>> {
    let num_sets = expect_ok(ffi::product_sets_size(
        sets.as_ref().expect("RDKit returned null products"),
    ));
    (0..num_sets)
        .map(|i| {
            let size = expect_ok(ffi::product_sets_set_size(&sets, i));
            (0..size)
                .map(|j| Mol::from_ffi(expect_ok(ffi::product_sets_take(sets.pin_mut(), i, j))))
                .collect()
        })
        .collect()
}

impl ChemicalReaction {
    /// Takes ownership of a reaction from the C++ wrappers.
    pub(crate) fn from_ffi(ptr: UniquePtr<ffi::ChemicalReaction>) -> ChemicalReaction {
        assert!(!ptr.is_null(), "RDKit returned a null reaction");
        ChemicalReaction { ptr }
    }

    pub(crate) fn as_ffi(&self) -> &ffi::ChemicalReaction {
        self.ptr.as_ref().unwrap()
    }

    pub(crate) fn pin_mut(&mut self) -> Pin<&mut ffi::ChemicalReaction> {
        self.ptr.pin_mut()
    }

    fn props(&self) -> &RDProps {
        expect_ok(ffi::chemical_reaction_as_props(self.as_ffi()))
    }

    fn num_templates(&self, kind: u8) -> u32 {
        expect_ok(ffi::chemical_reaction_get_num_templates(
            self.as_ffi(),
            kind,
        ))
    }

    fn template(&self, kind: u8, which: u32) -> TemplateRef<'_> {
        expect_ok_index(self.get_template(kind, which))
    }

    fn get_template(&self, kind: u8, which: u32) -> Result<TemplateRef<'_>> {
        let mol = ffi::chemical_reaction_get_template(self.as_ffi(), kind, which)?;
        Ok(TemplateRef::new(mol))
    }

    fn templates(&self, kind: u8) -> Templates<'_> {
        Templates {
            rxn: self,
            kind,
            range: 0..self.num_templates(kind),
        }
    }

    fn add_template(&mut self, kind: u8, mol: &Mol) -> u32 {
        expect_ok(ffi::chemical_reaction_add_template(
            self.pin_mut(),
            kind,
            mol.as_ffi(),
        ))
    }

    fn is_molecule_template(&self, kind: u8, mol: &Mol) -> Result<bool> {
        Ok(ffi::chemical_reaction_is_molecule_template(
            self.as_ffi(),
            kind,
            mol.as_ffi(),
        )?)
    }

    /// A reaction from a pickle (`to_binary`). A template whose plain
    /// monomer info has the PDB residue type is `Error::Value`, as for a
    /// molecule pickle (CONVENTIONS.md U10). Mirrors
    /// `ChemicalReaction(binStr)`.
    pub fn from_bin_str(bin_str: &[u8]) -> Result<ChemicalReaction> {
        Ok(ChemicalReaction::from_ffi(
            ffi::chemical_reaction_from_bin_str(bin_str)?,
        ))
    }

    /// The reaction as a pickle, with the default properties. Mirrors
    /// `ChemicalReaction.ToBinary`.
    pub fn to_binary(&self) -> Vec<u8> {
        expect_ok(ffi::chemical_reaction_to_binary(self.as_ffi(), false, 0))
    }

    /// Mirrors `ChemicalReaction.ToBinary(propertyFlags)`.
    pub fn to_binary_with_property_flags(&self, property_flags: PropertyPickleOptions) -> Vec<u8> {
        expect_ok(ffi::chemical_reaction_to_binary(
            self.as_ffi(),
            true,
            property_flags.bits(),
        ))
    }

    /// Mirrors `ChemicalReaction.GetNumReactantTemplates`.
    pub fn get_num_reactant_templates(&self) -> u32 {
        self.num_templates(REACTANT)
    }

    /// Mirrors `ChemicalReaction.GetNumProductTemplates`.
    pub fn get_num_product_templates(&self) -> u32 {
        self.num_templates(PRODUCT)
    }

    /// Mirrors `ChemicalReaction.GetNumAgentTemplates`.
    pub fn get_num_agent_templates(&self) -> u32 {
        self.num_templates(AGENT)
    }

    /// Adds a copy of `mol` as a reactant template, returning the number of
    /// reactant templates. Mirrors `ChemicalReaction.AddReactantTemplate`.
    pub fn add_reactant_template(&mut self, mol: &Mol) -> u32 {
        self.add_template(REACTANT, mol)
    }

    /// Adds a copy of `mol` as a product template. Mirrors
    /// `ChemicalReaction.AddProductTemplate`.
    pub fn add_product_template(&mut self, mol: &Mol) -> u32 {
        self.add_template(PRODUCT, mol)
    }

    /// Adds a copy of `mol` as an agent template. Mirrors
    /// `ChemicalReaction.AddAgentTemplate`.
    pub fn add_agent_template(&mut self, mol: &Mol) -> u32 {
        self.add_template(AGENT, mol)
    }

    /// Mirrors `ChemicalReaction.GetReactantTemplate`.
    pub fn get_reactant_template(&self, which: u32) -> Result<TemplateRef<'_>> {
        self.get_template(REACTANT, which)
    }

    /// Mirrors `ChemicalReaction.GetProductTemplate`.
    pub fn get_product_template(&self, which: u32) -> Result<TemplateRef<'_>> {
        self.get_template(PRODUCT, which)
    }

    /// Mirrors `ChemicalReaction.GetAgentTemplate`.
    pub fn get_agent_template(&self, which: u32) -> Result<TemplateRef<'_>> {
        self.get_template(AGENT, which)
    }

    /// Mirrors `ChemicalReaction.GetReactants`.
    pub fn get_reactants(&self) -> Templates<'_> {
        self.templates(REACTANT)
    }

    /// Mirrors `ChemicalReaction.GetProducts`.
    pub fn get_products(&self) -> Templates<'_> {
        self.templates(PRODUCT)
    }

    /// Mirrors `ChemicalReaction.GetAgents`.
    pub fn get_agents(&self) -> Templates<'_> {
        self.templates(AGENT)
    }

    /// Removes the reactant templates with less than a fifth of their heavy
    /// atoms mapped, moving them to the agents, and returns copies of them
    /// (Python's `targetList`). Mirrors
    /// `ChemicalReaction.RemoveUnmappedReactantTemplates`.
    pub fn remove_unmapped_reactant_templates(&mut self) -> Vec<Mol> {
        self.remove_unmapped_reactant_templates_with_params(&Default::default())
    }

    /// Mirrors `ChemicalReaction.RemoveUnmappedReactantTemplates(
    /// thresholdUnmappedAtoms, moveToAgentTemplates)`.
    pub fn remove_unmapped_reactant_templates_with_params(
        &mut self,
        params: &RemoveUnmappedReactantTemplatesParams,
    ) -> Vec<Mol> {
        mol_list_into_vec(expect_ok(ffi::chemical_reaction_remove_unmapped_templates(
            self.pin_mut(),
            REACTANT,
            params.threshold_unmapped_atoms,
            params.move_to_agent_templates,
        )))
    }

    /// Removes the product templates with less than a fifth of their heavy
    /// atoms mapped, moving them to the agents, and returns copies of them.
    /// Mirrors `ChemicalReaction.RemoveUnmappedProductTemplates`.
    pub fn remove_unmapped_product_templates(&mut self) -> Vec<Mol> {
        self.remove_unmapped_product_templates_with_params(&Default::default())
    }

    /// Mirrors `ChemicalReaction.RemoveUnmappedProductTemplates(
    /// thresholdUnmappedAtoms, moveToAgentTemplates)`.
    pub fn remove_unmapped_product_templates_with_params(
        &mut self,
        params: &RemoveUnmappedProductTemplatesParams,
    ) -> Vec<Mol> {
        mol_list_into_vec(expect_ok(ffi::chemical_reaction_remove_unmapped_templates(
            self.pin_mut(),
            PRODUCT,
            params.threshold_unmapped_atoms,
            params.move_to_agent_templates,
        )))
    }

    /// Removes the agent templates and returns them. Mirrors
    /// `ChemicalReaction.RemoveAgentTemplates`.
    pub fn remove_agent_templates(&mut self) -> Vec<Mol> {
        mol_list_into_vec(expect_ok(ffi::chemical_reaction_remove_agent_templates(
            self.pin_mut(),
        )))
    }

    /// Runs the reaction on one reactant per reactant template, returning a
    /// list of products (one per product template) for each way the
    /// templates match, at most 1000 lists. Initializes the reaction first
    /// if needed; a reaction that fails validation, or a wrong number of
    /// reactants, is `Error::Value`. The reactants are not changed. Mirrors
    /// `ChemicalReaction.RunReactants`.
    pub fn run_reactants(&mut self, reactants: &[Mol]) -> Result<Vec<Vec<Mol>>> {
        self.run_reactants_with_max_products(reactants, 1000)
    }

    /// Stops after `max_products` lists, or never when it is 0. Mirrors
    /// `ChemicalReaction.RunReactants(reactants, maxProducts)`.
    pub fn run_reactants_with_max_products(
        &mut self,
        reactants: &[Mol],
        max_products: u32,
    ) -> Result<Vec<Vec<Mol>>> {
        let products = with_borrowed_mols_shared(reactants.iter().map(|mol| ("", mol)), |list| {
            ffi::chemical_reaction_run_reactants(self.pin_mut(), list, max_products)
        })?;
        Ok(product_sets_into_vec(products))
    }

    /// Runs the reaction with one reactant, for the reactant template
    /// `reaction_idx` (`Error::Index` past the templates). Mirrors
    /// `ChemicalReaction.RunReactant`.
    pub fn run_reactant(&mut self, reactant: &Mol, reaction_idx: u32) -> Result<Vec<Vec<Mol>>> {
        Ok(product_sets_into_vec(ffi::chemical_reaction_run_reactant(
            self.pin_mut(),
            reactant.as_ffi(),
            reaction_idx,
        )?))
    }

    /// Runs a reaction with one reactant template and one product template on
    /// the reactant itself, removing reactant atoms the product template
    /// lacks, and returns whether the reactant changed. Other reactions are
    /// `Error::Value`. Mirrors `ChemicalReaction.RunReactantInPlace`.
    pub fn run_reactant_in_place(&mut self, reactant: &mut Mol) -> Result<bool> {
        self.run_reactant_in_place_with_remove_unmatched_atoms(reactant, true)
    }

    /// Keeps the atoms the template did not match unless
    /// `remove_unmatched_atoms`. Mirrors
    /// `ChemicalReaction.RunReactantInPlace(reactant, removeUnmatchedAtoms)`.
    pub fn run_reactant_in_place_with_remove_unmatched_atoms(
        &mut self,
        reactant: &mut Mol,
        remove_unmatched_atoms: bool,
    ) -> Result<bool> {
        Ok(ffi::chemical_reaction_run_reactant_in_place(
            self.pin_mut(),
            reactant.pin_mut(),
            remove_unmatched_atoms,
        )?)
    }

    /// Validates the reaction and readies it to run, logging its problems.
    /// A reaction that fails validation stays uninitialized. Mirrors
    /// `ChemicalReaction.Initialize`.
    pub fn initialize(&mut self) -> Result<()> {
        self.initialize_with_silent(false)
    }

    /// Logs nothing when `silent`. Mirrors
    /// `ChemicalReaction.Initialize(silent)`.
    pub fn initialize_with_silent(&mut self, silent: bool) -> Result<()> {
        Ok(ffi::chemical_reaction_initialize(self.pin_mut(), silent)?)
    }

    /// Mirrors `ChemicalReaction.IsInitialized`.
    pub fn is_initialized(&self) -> bool {
        expect_ok(ffi::chemical_reaction_is_initialized(self.as_ffi()))
    }

    /// Checks the reaction for problems, logging them, and returns the
    /// numbers of warnings and errors. It writes properties on the product
    /// templates. Mirrors `ChemicalReaction.Validate`.
    pub fn validate(&mut self) -> Result<(u32, u32)> {
        self.validate_with_silent(false)
    }

    /// Logs nothing when `silent`. Mirrors
    /// `ChemicalReaction.Validate(silent)`.
    pub fn validate_with_silent(&mut self, silent: bool) -> Result<(u32, u32)> {
        let counts = ffi::chemical_reaction_validate(self.pin_mut(), silent)?;
        Ok((counts[0], counts[1]))
    }

    /// Whether a reactant template matches part of `mol`. An uninitialized
    /// reaction is `Error::Value`. Mirrors
    /// `ChemicalReaction.IsMoleculeReactant`.
    pub fn is_molecule_reactant(&self, mol: &Mol) -> Result<bool> {
        self.is_molecule_template(REACTANT, mol)
    }

    /// Whether a product template matches part of `mol`. Mirrors
    /// `ChemicalReaction.IsMoleculeProduct`.
    pub fn is_molecule_product(&self, mol: &Mol) -> Result<bool> {
        self.is_molecule_template(PRODUCT, mol)
    }

    /// Whether an agent template matches `mol` whole. Mirrors
    /// `ChemicalReaction.IsMoleculeAgent`.
    pub fn is_molecule_agent(&self, mol: &Mol) -> Result<bool> {
        self.is_molecule_template(AGENT, mol)
    }

    /// The atoms of each reactant template that the reaction changes,
    /// unmapped atoms included. An uninitialized reaction is `Error::Value`.
    /// Mirrors `ChemicalReaction.GetReactingAtoms`.
    pub fn get_reacting_atoms(&self) -> Result<Vec<Vec<u32>>> {
        self.get_reacting_atoms_with_mapped_atoms_only(false)
    }

    /// Leaves out unmapped atoms when `mapped_atoms_only`. Mirrors
    /// `ChemicalReaction.GetReactingAtoms(mappedAtomsOnly)`.
    pub fn get_reacting_atoms_with_mapped_atoms_only(
        &self,
        mapped_atoms_only: bool,
    ) -> Result<Vec<Vec<u32>>> {
        Ok(
            ffi::chemical_reaction_get_reacting_atoms(self.as_ffi(), mapped_atoms_only)?
                .into_iter()
                .map(|list| list.indices)
                .collect(),
        )
    }

    fn add_recursive_queries(
        &mut self,
        params: &AddRecursiveQueriesToReactionParams,
        get_labels: bool,
    ) -> Result<Vec<ffi::ReactantLabel>> {
        Ok(with_borrowed_mols_shared(
            params
                .queries
                .iter()
                .map(|(name, mol)| (name.as_str(), mol)),
            |list| {
                ffi::chemical_reaction_add_recursive_queries(
                    self.pin_mut(),
                    list,
                    &params.prop_name,
                    get_labels,
                )
            },
        )?)
    }

    /// Replaces the reactant atoms labeled in their `molFileValue` property
    /// with the named recursive queries, of which there are none. A label
    /// with no query is `Error::Key`; an uninitialized reaction is
    /// `Error::Value`. Mirrors `ChemicalReaction.AddRecursiveQueriesToReaction`.
    pub fn add_recursive_queries_to_reaction(&mut self) -> Result<()> {
        self.add_recursive_queries_to_reaction_with_params(&Default::default())
    }

    /// Mirrors `ChemicalReaction.AddRecursiveQueriesToReaction(queries,
    /// propName)`.
    pub fn add_recursive_queries_to_reaction_with_params(
        &mut self,
        params: &AddRecursiveQueriesToReactionParams,
    ) -> Result<()> {
        self.add_recursive_queries(params, false).map(drop)
    }

    /// Also returns each reactant template's labels. Mirrors
    /// `ChemicalReaction.AddRecursiveQueriesToReaction(getLabels=True)`.
    pub fn add_recursive_queries_to_reaction_get_labels(&mut self) -> Result<ReactantLabels> {
        self.add_recursive_queries_to_reaction_get_labels_with_params(&Default::default())
    }

    /// Mirrors `ChemicalReaction.AddRecursiveQueriesToReaction(queries,
    /// propName, getLabels=True)`.
    pub fn add_recursive_queries_to_reaction_get_labels_with_params(
        &mut self,
        params: &AddRecursiveQueriesToReactionParams,
    ) -> Result<ReactantLabels> {
        let labels = self.add_recursive_queries(params, true)?;
        group_labels(labels, self.get_num_reactant_templates())
    }

    /// A copy of the parameters for matching the reactant templates. Mirrors
    /// `ChemicalReaction.GetSubstructParams`, which Python returns for
    /// changing in place; `set_substruct_params` writes a changed copy back.
    pub fn get_substruct_params(&self) -> SubstructMatchParameters {
        SubstructMatchParameters::from_ffi(expect_ok(ffi::chemical_reaction_get_substruct_params(
            self.as_ffi(),
        )))
    }

    /// Sets the parameters for matching the reactant templates. Python
    /// changes the object `GetSubstructParams` returns instead.
    pub fn set_substruct_params(&mut self, params: &SubstructMatchParameters) {
        expect_ok(ffi::chemical_reaction_set_substruct_params(
            self.pin_mut(),
            &params.to_ffi(),
        ))
    }
}

/// Unwraps an index the caller checked.
fn expect_ok_index<T>(result: Result<T>) -> T {
    result.unwrap_or_else(|err| panic!("a checked index was out of range: {err}"))
}

prop_getters!(ChemicalReaction);
prop_setters!(ChemicalReaction, clear_computed);

/// An empty reaction. Mirrors `ChemicalReaction()`.
impl Default for ChemicalReaction {
    fn default() -> ChemicalReaction {
        ChemicalReaction::from_ffi(expect_ok(ffi::chemical_reaction_new()))
    }
}

/// A deep copy, templates included. Mirrors `ChemicalReaction(other)`.
impl Clone for ChemicalReaction {
    fn clone(&self) -> ChemicalReaction {
        ChemicalReaction::from_ffi(expect_ok(ffi::chemical_reaction_clone(self.as_ffi())))
    }
}

impl fmt::Debug for ChemicalReaction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ChemicalReaction")
            .field("num_reactant_templates", &self.get_num_reactant_templates())
            .field("num_product_templates", &self.get_num_product_templates())
            .field("num_agent_templates", &self.get_num_agent_templates())
            .finish()
    }
}
