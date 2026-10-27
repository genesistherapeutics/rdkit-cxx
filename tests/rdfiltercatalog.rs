// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use rdkit_cxx::Error;
use rdkit_cxx::chem::rdfiltercatalog::filter_match_ops::{And, Not, Or};
use rdkit_cxx::chem::rdfiltercatalog::*;
use rdkit_cxx::chem::*;

// Expected values come from Python's rdfiltercatalog.

/// Python's archive of a catalog of one entry, "root", whose matcher is a
/// FilterHierarchyMatcher without a pattern.
const EMPTY_NODE_ARCHIVE: &str = "32322073657269616c697a6174696f6e3a3a617263686976652032302030203020312031203020312032203120310a30203020312038203120310a312030203020323520556e6e616d65642046696c7465724d617463686572426173652030203020302031202d312030203020322030203131206465736372697074696f6e203420726f6f74";

/// Python's archive of a catalog of one entry, "cycle", whose hierarchy's
/// child holds a copy of the hierarchy.
const CYCLE_ARCHIVE: &str = "32322073657269616c697a6174696f6e3a3a617263686976652032302030203020312031203020312032203120310a30203020312038203120310a312030203020323520556e6e616d65642046696c7465724d6174636865724261736520302030203120312030203120380a3220323520556e6e616d65642046696c7465724d617463686572426173652031203120380a3320323520556e6e616d65642046696c7465724d6174636865724261736520312031203820322036203120310a342031204320373920efbeadde000000001000000004000000000000000100000000000000800106100800000000190d00000041746f6d41746f6d69634e756d211a06000000000000001c002b0b170400000000000000162031203432393439363732393520360a352031204f20373920efbeadde000000001000000004000000000000000100000000000000800108100800000000190d00000041746f6d41746f6d69634e756d211a08000000000000001c002b0b1704000000000000001620312034323934393637323935203620342030203020322030203131206465736372697074696f6e2035206379636c65";

/// Python's archive of `small_catalog()`.
const SMALL_ARCHIVE: &str = "32322073657269616c697a6174696f6e3a3a617263686976652032302030203020322031203020312032203120310a30203020312036203120310a31203020302031204320373920efbeadde000000001000000004000000000000000100000000000000800106100800000000190d00000041746f6d41746f6d69634e756d211a06000000000000001c002b0b17040000000000000016203220332030203020342030203131206465736372697074696f6e203620636172626f6e20352053636f7065203520746573747320320a322033203120300a33203320416e6420360a342031204320373920efbeadde000000001000000004000000000000000100000000000000800106100800000000190d00000041746f6d41746f6d69634e756d211a06000000000000001c002b0b17040000000000000016203120343239343936373239352035203120300a352033204e6f7420360a362031204e20373920efbeadde000000001000000004000000000000000100000000000000800107100800000000190d00000041746f6d41746f6d69634e756d211a07000000000000001c002b0b170400000000000000162031203432393439363732393520322030203131206465736372697074696f6e203320616e64";

const HZONE: &str = "O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2";

fn hex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn mol(smiles: &str) -> Mol {
    mol_from_smiles(smiles).unwrap()
}

fn smarts(name: &str, smarts: &str) -> SmartsMatcher {
    SmartsMatcher::from_smarts(name, smarts)
}

fn names(matches: &[FilterMatch]) -> Vec<(String, Vec<(u32, u32)>)> {
    matches
        .iter()
        .map(|m| (m.filter_match.get_name().unwrap(), m.atom_pairs.clone()))
        .collect()
}

fn descriptions(entries: &[FilterCatalogEntry]) -> Vec<String> {
    entries
        .iter()
        .map(|e| e.get_description().unwrap())
        .collect()
}

fn is_value<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Value { .. }))
}

fn is_runtime<T>(result: rdkit_cxx::Result<T>) -> bool {
    matches!(result, Err(Error::Runtime { .. }))
}

/// Entries "carbon", which matches two or three carbons (with property
/// Scope = tests), and "and", carbon and not nitrogen.
fn small_catalog() -> FilterCatalog {
    let mut catalog = FilterCatalog::default();
    let counts = SmartsMatcherParams {
        min_count: 2,
        max_count: 3,
    };
    let mut carbon = FilterCatalogEntry::new(
        "carbon",
        &SmartsMatcher::from_smarts_with_params("C", "[#6]", &counts),
    )
    .unwrap();
    carbon.set_prop("Scope", "tests");
    catalog.add_entry(&carbon);
    let not_n = Not::new(&smarts("N", "[#7]")).unwrap();
    let and = And::new(&smarts("C", "[#6]"), &not_n).unwrap();
    catalog.add_entry(&FilterCatalogEntry::new("and", &and).unwrap());
    catalog
}

fn check_small_catalog(catalog: &FilterCatalog) {
    let ethanol = mol("CCO");
    assert_eq!(catalog.get_num_entries(), 2);
    assert_eq!(
        descriptions(&catalog.get_matches(&ethanol).unwrap()),
        ["carbon", "and"]
    );
    assert_eq!(
        catalog.get_entry(0).unwrap().get_prop("Scope").unwrap(),
        "tests"
    );
    // RDKit's And replaces the matches gathered so far with its own.
    assert_eq!(
        names(&catalog.get_filter_matches(&ethanol).unwrap()),
        [("C".to_string(), vec![(0, 0)])]
    );
    let carbon = catalog.get_entry(0).unwrap();
    assert_eq!(
        names(&carbon.get_filter_matches(&ethanol).unwrap()),
        [
            ("C".to_string(), vec![(0, 0)]),
            ("C".to_string(), vec![(0, 1)])
        ]
    );
    assert!(catalog.has_match(&mol("C")).unwrap());
    assert!(!catalog.has_match(&mol("N")).unwrap());
    assert!(catalog.get_first_match(&mol("N")).unwrap().is_none());
}

#[test]
fn built_in_catalogs() {
    let sizes = [
        (FilterCatalogs::PAINS_A, 16),
        (FilterCatalogs::PAINS_B, 55),
        (FilterCatalogs::PAINS_C, 409),
        (FilterCatalogs::PAINS, 480),
        (FilterCatalogs::BRENK, 105),
        (FilterCatalogs::NIH, 180),
        (FilterCatalogs::ZINC, 50),
        (FilterCatalogs::CHEMBL_Glaxo, 55),
        (FilterCatalogs::CHEMBL_Dundee, 105),
        (FilterCatalogs::CHEMBL_BMS, 180),
        (FilterCatalogs::CHEMBL_SureChEMBL, 166),
        (FilterCatalogs::CHEMBL_MLSMR, 116),
        (FilterCatalogs::CHEMBL_Inpharmatica, 91),
        (FilterCatalogs::CHEMBL_LINT, 57),
        (FilterCatalogs::CHEMBL, 770),
    ];
    for (catalogs, size) in sizes {
        assert_eq!(
            FilterCatalog::from_catalogs(catalogs).get_num_entries(),
            size,
            "{catalogs:?}"
        );
    }
    // RDKit casts the bits to its enum, where bits outside its range are
    // undefined; they are dropped first, as RDKit ignores them.
    let all = FilterCatalog::from_catalogs(FilterCatalogs::from_bits_retain(u32::MAX));
    assert_eq!(all.get_num_entries(), 1585);

    let mut params = FilterCatalogParams::default();
    assert!(!params.add_catalog(FilterCatalogs::from_bits_retain(1 | (1 << 20))));
    assert!(params.add_catalog(FilterCatalogs::PAINS_A));
    assert!(params.add_catalog(FilterCatalogs::BRENK));
    let catalog = FilterCatalog::new(&params);
    assert_eq!(catalog.get_num_entries(), 121);
    let empty = FilterCatalogParams::new(FilterCatalogs::from_bits_retain(1 << 20));
    assert_eq!(FilterCatalog::new(&empty).get_num_entries(), 0);

    let hzone = mol(HZONE);
    let first = catalog.get_first_match(&hzone).unwrap().unwrap();
    assert_eq!(first.get_description().unwrap(), "hzone_phenol_A(479)");
    assert_eq!(first.get_prop("Scope").unwrap(), "PAINS filters (family A)");
    let mut props = first.get_prop_list().unwrap();
    props.sort();
    assert_eq!(props, ["FilterSet", "Reference", "Scope", "description"]);
    assert_eq!(
        descriptions(&catalog.get_matches(&hzone).unwrap()),
        [
            "hzone_phenol_A(479)",
            "imine_1",
            "Oxygen-nitrogen_single_bond"
        ]
    );
    assert_eq!(
        names(&catalog.get_filter_matches(&hzone).unwrap()),
        [
            (
                "hzone_phenol_A(479)".to_string(),
                vec![
                    (0, 23),
                    (1, 22),
                    (2, 20),
                    (3, 19),
                    (4, 25),
                    (5, 24),
                    (6, 18),
                    (7, 17),
                    (8, 16),
                    (9, 21)
                ]
            ),
            ("imine_1".to_string(), vec![(0, 18), (1, 17)]),
            (
                "Oxygen-nitrogen_single_bond".to_string(),
                vec![(0, 16), (1, 17)]
            ),
        ]
    );
    assert!(catalog.has_match(&hzone).unwrap());
    assert!(!catalog.has_match(&mol("CC")).unwrap());
    assert!(catalog.get_first_match(&mol("CC")).unwrap().is_none());

    let entry = FilterCatalog::from_catalogs(FilterCatalogs::PAINS_A)
        .get_entry(0)
        .unwrap();
    assert_eq!(entry.get_description().unwrap(), "ene_six_het_A(483)");
    assert_eq!(entry.get_prop("FilterSet").unwrap(), "PAINS_A");
    assert!(matches!(entry.get_prop("missing"), Err(Error::Key { .. })));
}

#[test]
fn entries() {
    let entry = FilterCatalogEntry::default();
    assert!(!entry.is_valid());
    assert_eq!(entry.get_description().unwrap(), "");
    assert!(!entry.has_filter_match(&mol("C")).unwrap());
    assert!(entry.get_filter_matches(&mol("C")).unwrap().is_empty());

    let mut entry = FilterCatalogEntry::new("carbon", &smarts("C", "[#6]")).unwrap();
    assert!(entry.is_valid());
    assert_eq!(entry.get_description().unwrap(), "carbon");
    entry.set_prop("Scope", "tests");
    assert_eq!(entry.get_prop_list().unwrap(), ["description", "Scope"]);
    entry.set_description("C");
    assert_eq!(entry.get_prop("description").unwrap(), "C");
    entry.clear_prop("Scope");
    entry.clear_prop("missing");
    assert_eq!(entry.get_prop_list().unwrap(), ["description"]);
    assert!(entry.has_filter_match(&mol("C")).unwrap());
    assert!(!entry.has_filter_match(&mol("N")).unwrap());
    assert!(!entry.serialize().unwrap().is_empty());

    // An invalid matcher makes an invalid entry, which matches nothing.
    let invalid = FilterCatalogEntry::new("bad", &smarts("bad", "C(")).unwrap();
    assert!(!invalid.is_valid());
    assert!(!invalid.has_filter_match(&mol("C")).unwrap());
}

#[test]
fn catalogs() {
    let catalog = small_catalog();
    check_small_catalog(&catalog);
    assert!(filter_catalog_can_serialize());
    let archive = catalog.serialize().unwrap();
    // SAFETY: bytes serialize wrote, and bytes Python's Serialize wrote.
    let back = unsafe { FilterCatalog::from_pickle(&archive) }.unwrap();
    check_small_catalog(&back);
    let python = unsafe { FilterCatalog::from_pickle(&hex(SMALL_ARCHIVE)) }.unwrap();
    check_small_catalog(&python);

    let catalog = FilterCatalog::from_catalogs(FilterCatalogs::PAINS_A);
    let archive = catalog.serialize().unwrap();
    // SAFETY: bytes serialize wrote.
    let back = unsafe { FilterCatalog::from_pickle(&archive) }.unwrap();
    assert_eq!(back.get_num_entries(), 16);
    assert_eq!(
        descriptions(&back.get_matches(&mol(HZONE)).unwrap()),
        ["hzone_phenol_A(479)"]
    );

    assert!(matches!(catalog.get_entry(16), Err(Error::Index { .. })));
    assert!(matches!(
        catalog.get_entry_with_idx(16),
        Err(Error::Index { .. })
    ));
    assert_eq!(
        catalog
            .get_entry_with_idx(1)
            .unwrap()
            .get_description()
            .unwrap(),
        catalog.get_entry(1).unwrap().get_description().unwrap()
    );
}

// Python's entries are the catalog's own, so changing one changes the
// catalog, from any thread that holds it. Here a change to a shared entry
// changes a copy; RemoveEntry finds entries by identity, as in Python.
#[test]
fn entries_are_shared_until_changed() {
    let mut catalog = small_catalog();
    let mut entry = catalog.get_entry(0).unwrap();
    assert!(entry.is_same_entry(&catalog.get_entry(0).unwrap()));
    let copy = entry.clone();
    assert!(copy.is_same_entry(&entry));
    entry.set_description("changed");
    assert!(!entry.is_same_entry(&copy));
    assert_eq!(copy.get_description().unwrap(), "carbon");
    assert_eq!(
        catalog.get_entry(0).unwrap().get_description().unwrap(),
        "carbon"
    );

    // AddEntry adds a copy, which the original is not.
    let added = FilterCatalogEntry::new("oxygen", &smarts("O", "[#8]")).unwrap();
    catalog.add_entry(&added);
    assert!(!catalog.remove_entry_filter_catalog_entry(&added));
    assert!(!catalog.remove_entry_filter_catalog_entry(&entry));
    let third = catalog.get_entry(2).unwrap();
    assert!(catalog.remove_entry_filter_catalog_entry(&third));
    assert_eq!(catalog.get_num_entries(), 2);
    // The removed entry lives on in its handle.
    assert_eq!(third.get_description().unwrap(), "oxygen");
    assert!(catalog.remove_entry(1));
    assert!(!catalog.remove_entry(1));
    assert_eq!(catalog.get_num_entries(), 1);
}

#[test]
fn smarts_matchers() {
    let mut matcher = SmartsMatcher::from_smarts("C", "[#6]");
    assert_eq!(matcher.get_name().unwrap(), "C");
    assert_eq!(matcher.to_string(), "C");
    assert_eq!(
        (matcher.get_min_count(), matcher.get_max_count()),
        (1, u32::MAX)
    );
    assert!(matcher.is_valid());
    assert_eq!(
        mol_to_smarts(&matcher.get_pattern().unwrap()).unwrap(),
        "[#6]"
    );
    assert_eq!(
        names(&matcher.get_matches(&mol("CCO")).unwrap()),
        [("C".to_string(), vec![(0, 0)])]
    );
    matcher.set_min_count(2);
    assert_eq!(
        names(&matcher.get_matches(&mol("CCO")).unwrap()),
        [
            ("C".to_string(), vec![(0, 0)]),
            ("C".to_string(), vec![(0, 1)])
        ]
    );
    assert!(!matcher.has_match(&mol("C")).unwrap());
    matcher.set_max_count(3);
    assert!(matcher.has_match(&mol("CCC")).unwrap());
    assert!(!matcher.has_match(&mol("CCCC")).unwrap());

    let counts = SmartsMatcherParams {
        min_count: 1,
        max_count: 1,
    };
    let mut once = SmartsMatcher::from_smarts_with_params("O", "[#8]", &counts);
    assert!(once.get_matches(&mol("OCO")).unwrap().is_empty());
    assert!(!once.has_match(&mol("OCO")).unwrap());
    assert!(once.has_match(&mol("CO")).unwrap());
    // SetPattern keeps the counts.
    once.set_pattern_str("[#6]");
    assert_eq!((once.get_min_count(), once.get_max_count()), (1, 1));
    assert!(once.has_match(&mol("CO")).unwrap());
    once.set_pattern_str("C(");
    assert!(!once.is_valid());
    assert!(once.get_pattern().is_none());
    // RDKit's precondition, as in Python.
    assert!(is_runtime(once.has_match(&mol("CO"))));
    once.set_pattern(&mol_from_smarts("[#7]").unwrap());
    assert!(once.has_match(&mol("CN")).unwrap());

    let invalid = SmartsMatcher::from_smarts("bad", "C(");
    assert!(!invalid.is_valid());
    assert!(invalid.get_pattern().is_none());

    let rhs = SmartsMatcher::from_rhs(&mol_from_smarts("[#7]").unwrap());
    assert_eq!(rhs.get_name().unwrap(), "Unnamed SmartsMatcher");
    assert_eq!((rhs.get_min_count(), rhs.get_max_count()), (1, u32::MAX));
    let named = SmartsMatcher::new("n4");
    assert_eq!(named.get_name().unwrap(), "n4");
    assert!(!named.is_valid());
    assert_eq!(
        (named.get_min_count(), named.get_max_count()),
        (0, u32::MAX)
    );
    let from_mol = SmartsMatcher::from_mol("n5", &mol_from_smarts("[#7]").unwrap());
    assert_eq!(
        (from_mol.get_min_count(), from_mol.get_max_count()),
        (1, u32::MAX)
    );
    assert_eq!(
        SmartsMatcherParams::default(),
        SmartsMatcherParams {
            min_count: 1,
            max_count: u32::MAX
        }
    );

    // Clones are copies.
    let mut clone = matcher.clone();
    clone.set_min_count(5);
    assert_eq!(matcher.get_min_count(), 2);
}

#[test]
fn composite_matchers() {
    let c = smarts("C", "[#6]");
    let o = smarts("O", "[#8]");
    let and = And::new(&c, &o).unwrap();
    assert_eq!(and.get_name().unwrap(), "(C And O)");
    assert!(and.has_match(&mol("CO")).unwrap());
    assert!(!and.has_match(&mol("CC")).unwrap());
    assert_eq!(
        names(&and.get_matches(&mol("CO")).unwrap()),
        [
            ("C".to_string(), vec![(0, 0)]),
            ("O".to_string(), vec![(0, 1)])
        ]
    );
    let or = Or::new(&c, &o).unwrap();
    assert_eq!(or.get_name().unwrap(), "(C Or O)");
    assert_eq!(names(&or.get_matches(&mol("CO")).unwrap()).len(), 2);
    let not = Not::new(&c).unwrap();
    assert_eq!(not.get_name().unwrap(), "(Not C)");
    assert!(not.has_match(&mol("O")).unwrap());
    assert!(not.get_matches(&mol("O")).unwrap().is_empty());
    let nested = Not::new(&And::new(&c, &ExclusionList::default()).unwrap()).unwrap();
    assert_eq!(nested.get_name().unwrap(), "(Not (C And (Not any of)))");
    // An invalid argument makes an invalid matcher, which RDKit refuses.
    let invalid = And::new(&smarts("bad", "C("), &c).unwrap();
    assert!(!invalid.is_valid());
    assert!(is_runtime(invalid.has_match(&mol("C"))));

    let mut list = ExclusionList::default();
    assert_eq!(list.get_name().unwrap(), "(Not any of)");
    assert!(list.is_valid());
    assert!(list.has_match(&mol("C")).unwrap());
    list.add_pattern(&o).unwrap();
    assert_eq!(list.get_name().unwrap(), "(Not any of O)");
    assert!(list.has_match(&mol("CC")).unwrap());
    assert!(!list.has_match(&mol("CO")).unwrap());
    assert!(list.get_matches(&mol("CC")).unwrap().is_empty());
    assert!(is_runtime(list.add_pattern(&smarts("bad", "C("))));
    list.set_exclusion_patterns(&[&smarts("N", "[#7]"), &smarts("S", "[#16]")])
        .unwrap();
    assert_eq!(list.get_name().unwrap(), "(Not any of N S)");

    let m = FilterMatch::new(&c, &[(0, 1), (2, 3)]);
    assert_eq!(m.filter_match.get_name().unwrap(), "C");
    assert_eq!(m.atom_pairs, [(0, 1), (2, 3)]);
}

#[test]
fn hierarchies() {
    let mut root = FilterHierarchyMatcher::new(&smarts("C", "[#6]")).unwrap();
    assert_eq!(root.get_name().unwrap(), "C");
    assert!(root.is_valid());
    let co = FilterHierarchyMatcher::new(&smarts("CO", "[#6][#8]")).unwrap();
    let coc = FilterHierarchyMatcher::new(&smarts("COC", "[#6][#8][#6]")).unwrap();
    let cc = FilterHierarchyMatcher::new(&smarts("CC", "[#6][#6]")).unwrap();
    {
        let mut child = root.add_child(&co).unwrap();
        assert_eq!(child.get_name().unwrap(), "CO");
        child.add_child(&coc).unwrap();
    }
    root.add_child(&cc).unwrap();
    let expected = [
        ("COC", vec![("COC", vec![(0, 0), (1, 1), (2, 2)])]),
        (
            "CCO",
            vec![("CO", vec![(0, 1), (1, 2)]), ("CC", vec![(0, 0), (1, 1)])],
        ),
        ("C", vec![("C", vec![(0, 0)])]),
        ("O", vec![]),
    ];
    for (smiles, matches) in expected {
        let expected: Vec<_> = matches
            .into_iter()
            .map(|(name, pairs)| (name.to_string(), pairs))
            .collect();
        assert_eq!(names(&root.get_matches(&mol(smiles)).unwrap()), expected);
        assert_eq!(root.has_match(&mol(smiles)).unwrap(), smiles != "O");
    }
    // The tree holds copies: changing the original child changes nothing.
    let mut co = co;
    co.set_pattern(&smarts("N", "[#7]")).unwrap();
    assert_eq!(names(&root.get_matches(&mol("COC")).unwrap())[0].0, "COC");
    let entry = FilterCatalogEntry::new("tree", &root).unwrap();
    assert!(entry.has_filter_match(&mol("COC")).unwrap());

    // The functional groups are such trees.
    let groups = get_functional_group_hierarchy();
    assert_eq!(groups.get_num_entries(), 12);
    assert_eq!(
        names(&groups.get_filter_matches(&mol("c1ccccc1Cl")).unwrap()),
        [
            ("Halogen.Aromatic".to_string(), vec![(0, 6)]),
            ("Halogen.NotFluorine.Aromatic".to_string(), vec![(0, 6)])
        ]
    );
    assert_eq!(
        descriptions(&groups.get_matches(&mol("OCC(=O)O")).unwrap()),
        ["CarboxylicAcid", "Alcohol"]
    );
    let flat = get_flattened_functional_group_hierarchy();
    assert_eq!(flat.len(), 52);
    assert_eq!(
        mol_to_smarts(&flat["Alcohol.Aromatic"]).unwrap(),
        "[O&H1&$(O-&!@c)]"
    );
    assert_eq!(
        flat["BoronicAcid"].get_prop("Label").unwrap(),
        "Boronic Acid"
    );
    let normalized = get_flattened_functional_group_hierarchy_with_normalized(true);
    assert!(normalized.contains_key("halogen.bromine.aliphatic"));
    assert_eq!(normalized.len(), 52);
}

// A FilterHierarchyMatcher() has no pattern, which RDKit's IsValid,
// HasMatch, and GetMatches dereference; Python's crash.
#[test]
fn empty_hierarchy_nodes() {
    let mut root = FilterHierarchyMatcher::default();
    assert_eq!(root.get_name().unwrap(), "FilterMatcherHierarchy root");
    assert!(!root.is_valid());
    assert!(is_value(root.has_match(&mol("C"))));
    assert!(is_value(root.get_matches(&mol("C"))));
    // Nothing else takes one, since its methods would.
    let c = smarts("C", "[#6]");
    assert!(is_value(FilterCatalogEntry::new("root", &root)));
    assert!(is_value(And::new(&root, &c)));
    assert!(is_value(Or::new(&c, &root)));
    assert!(is_value(Not::new(&root)));
    assert!(is_value(FilterHierarchyMatcher::new(&root)));
    let mut list = ExclusionList::default();
    assert!(is_value(list.add_pattern(&root)));
    assert!(is_value(list.set_exclusion_patterns(&[&c, &root])));
    assert_eq!(list.get_name().unwrap(), "(Not any of)");
    let mut node = FilterHierarchyMatcher::new(&c).unwrap();
    assert!(is_value(node.set_pattern(&root)));
    // RDKit's own precondition, as in Python.
    assert!(is_runtime(node.add_child(&root)));

    // A pattern makes it a node like any other.
    root.set_pattern(&c).unwrap();
    assert!(root.is_valid());
    assert!(root.has_match(&mol("C")).unwrap());
    assert!(FilterCatalogEntry::new("root", &root).is_ok());
}

/// Matches molecules with at least `atoms` atoms, matching atom 0.
struct AtLeast {
    atoms: u32,
    calls: Arc<AtomicUsize>,
}

impl CustomFilterMatcher for AtLeast {
    fn has_match(&self, mol: &Mol) -> bool {
        self.calls.fetch_add(1, Ordering::SeqCst);
        mol.get_num_atoms() >= self.atoms
    }

    fn get_matches(&self, mol: &Mol, matches: &mut Vec<FilterMatch>) -> bool {
        let matched = self.has_match(mol);
        let me = FilterMatcher::new(
            "at least",
            AtLeast {
                atoms: self.atoms,
                calls: self.calls.clone(),
            },
        );
        matches.push(FilterMatch::new(&me, &[(0, 0)]));
        matched
    }
}

struct Panics;

impl CustomFilterMatcher for Panics {
    fn has_match(&self, _mol: &Mol) -> bool {
        panic!("no match for you")
    }

    fn get_matches(&self, _mol: &Mol, _matches: &mut Vec<FilterMatch>) -> bool {
        panic!("no matches for you")
    }
}

struct Invalid;

impl CustomFilterMatcher for Invalid {
    fn is_valid(&self) -> bool {
        false
    }

    fn has_match(&self, _mol: &Mol) -> bool {
        true
    }

    fn get_matches(&self, _mol: &Mol, _matches: &mut Vec<FilterMatch>) -> bool {
        true
    }
}

#[test]
fn custom_matchers() {
    let calls = Arc::new(AtomicUsize::new(0));
    let three = FilterMatcher::new(
        "three",
        AtLeast {
            atoms: 3,
            calls: calls.clone(),
        },
    );
    assert_eq!(three.get_name().unwrap(), "three");
    assert!(three.is_valid());
    assert!(three.has_match(&mol("CCO")).unwrap());
    assert!(!three.has_match(&mol("CO")).unwrap());
    assert_eq!(
        names(&three.get_matches(&mol("CCO")).unwrap()),
        [("at least".to_string(), vec![(0, 0)])]
    );
    // As in Python, the matches count only when it matches.
    assert!(three.get_matches(&mol("CO")).unwrap().is_empty());

    // In catalogs, with RDKit's matchers.
    let mut catalog = small_catalog();
    let both = And::new(&three, &smarts("O", "[#8]")).unwrap();
    catalog.add_entry(&FilterCatalogEntry::new("three and O", &both).unwrap());
    assert_eq!(
        descriptions(&catalog.get_matches(&mol("CCO")).unwrap()),
        ["carbon", "and", "three and O"]
    );
    assert_eq!(
        descriptions(&catalog.get_matches(&mol("CCN")).unwrap()),
        ["carbon"]
    );
    let before = calls.load(Ordering::SeqCst);
    let results =
        run_filter_catalog_with_num_threads(&catalog, &["CCO", "CO", "CCCO", "OCO"], 3).unwrap();
    let results: Vec<_> = results.iter().map(|row| descriptions(row)).collect();
    assert_eq!(
        results,
        [
            vec!["carbon", "and", "three and O"],
            vec!["and"],
            vec!["carbon", "and", "three and O"],
            vec!["and", "three and O"],
        ]
    );
    assert!(calls.load(Ordering::SeqCst) > before);
    // Boost knows only RDKit's matcher classes.
    assert!(is_value(catalog.serialize()));

    // A panic is caught at the bridge, on any thread.
    let panics = FilterMatcher::new("panics", Panics);
    match panics.has_match(&mol("C")) {
        Err(Error::Runtime { message }) => {
            assert!(message.contains("no match for you"), "{message}")
        }
        other => panic!("{other:?}"),
    }
    assert!(is_runtime(panics.get_matches(&mol("C"))));
    catalog.add_entry(&FilterCatalogEntry::new("panics", &panics).unwrap());
    assert!(is_runtime(catalog.has_match(&mol("N"))));
    assert!(is_runtime(run_filter_catalog_with_num_threads(
        &catalog,
        &["CCO", "CO", "C"],
        3
    )));

    // An invalid matcher's entry matches nothing.
    let invalid =
        FilterCatalogEntry::new("invalid", &FilterMatcher::new("invalid", Invalid)).unwrap();
    assert!(!invalid.is_valid());
    assert!(!invalid.has_filter_match(&mol("C")).unwrap());
}

#[test]
fn run_filter_catalogs() {
    let catalog = FilterCatalog::from_catalogs(FilterCatalogs::PAINS_A | FilterCatalogs::BRENK);
    let smiles = ["CC", "xx", HZONE];
    let expected = [
        vec![],
        vec!["no valid RDKit molecule"],
        vec![
            "hzone_phenol_A(479)",
            "imine_1",
            "Oxygen-nitrogen_single_bond",
        ],
    ];
    // RDKit counts threads in int, negating nonpositive counts (undefined
    // for i32::MIN); more threads than SMILES have nothing to do.
    for threads in [1, 2, 0, -1, 64, i32::MIN, i32::MAX] {
        let results = run_filter_catalog_with_num_threads(&catalog, &smiles, threads).unwrap();
        let results: Vec<_> = results.iter().map(|row| descriptions(row)).collect();
        assert_eq!(results, expected, "{threads}");
    }
    let results = run_filter_catalog(&catalog, &smiles).unwrap();
    let bad = &results[1][0];
    assert!(!bad.is_valid());
    assert_eq!(bad.get_prop_list().unwrap(), ["description"]);
    assert!(
        run_filter_catalog(&catalog, &[] as &[&str])
            .unwrap()
            .is_empty()
    );
    // A sanitization error is the call's, as in Python.
    assert!(matches!(
        run_filter_catalog_with_num_threads(&catalog, &["CC", "c1cc1"], 2),
        Err(Error::Kekulize { .. })
    ));
    // The unparsable entry is shared; changing it changes a copy.
    let mut first = run_filter_catalog(&catalog, &["xx"])
        .unwrap()
        .remove(0)
        .remove(0);
    first.set_description("mine");
    let again = run_filter_catalog(&catalog, &["xx"]).unwrap();
    assert_eq!(
        again[0][0].get_description().unwrap(),
        "no valid RDKit molecule"
    );
}

#[test]
fn serializing_checks_patterns() {
    // SmartsMatcher's save dereferences its pattern; Python's raises after
    // dereferencing null.
    let entry = FilterCatalogEntry::new("bad", &smarts("bad", "C(")).unwrap();
    assert!(is_value(entry.serialize()));
    let mut catalog = FilterCatalog::default();
    catalog.add_entry(&entry);
    assert!(is_value(catalog.serialize()));
    let nested = And::new(&smarts("C", "[#6]"), &SmartsMatcher::new("empty")).unwrap();
    let entry = FilterCatalogEntry::new("nested", &nested).unwrap();
    assert!(is_value(entry.serialize()));
}

// Archives RDKit wrote can still hold what Rust cannot make.
#[test]
fn loading_checks_matchers() {
    // SAFETY: bytes Python's Serialize wrote.
    let empty = unsafe { FilterCatalog::from_pickle(&hex(EMPTY_NODE_ARCHIVE)) };
    assert!(is_value(empty));
    // SAFETY: as above. RDKit's recursive methods would follow the cycle.
    let cycle = unsafe { FilterCatalog::from_pickle(&hex(CYCLE_ARCHIVE)) };
    assert!(is_value(cycle));

    // A pattern whose atom has monomer info of type PDBRESIDUE but no
    // residue fields, which RDKit 2024.09's pickler reads as PDB residue
    // info (U10), as Python's FilterCatalog does.
    let mut pattern = mol_from_smarts("[#6]").unwrap();
    pattern
        .get_atom_with_idx_mut(0)
        .unwrap()
        .set_monomer_info(&AtomMonomerInfo::new_with_name(AtomMonomerType::Other, "MARK").unwrap());
    let mut catalog = FilterCatalog::default();
    let matcher = SmartsMatcher::from_mol("marked", &pattern);
    catalog.add_entry(&FilterCatalogEntry::new("marked", &matcher).unwrap());
    let mut archive = catalog.serialize().unwrap();
    let at = archive
        .windows(9)
        .position(|w| w == b"\x04\0\0\0MARK\x02")
        .unwrap();
    // SAFETY: bytes serialize wrote.
    assert!(unsafe { FilterCatalog::from_pickle(&archive) }.is_ok());
    archive[at + 8] = 1;
    // SAFETY: bytes serialize wrote, but for the type byte of a molecule
    // pickle inside, which RDKit's pickler reads as it reads any.
    let catalog = unsafe { FilterCatalog::from_pickle(&archive) }.unwrap();
    assert_eq!(catalog.get_num_entries(), 1);
    assert!(catalog.has_match(&mol_from_smiles("C").unwrap()).unwrap());
}

/// Matches molecules with oxygen by calling functions that match: a
/// substructure search, a SmartsMatcher, and another Rust matcher.
struct Reentrant {
    oxygen: SmartsMatcher,
    other: FilterMatcher,
}

impl CustomFilterMatcher for Reentrant {
    fn has_match(&self, mol: &Mol) -> bool {
        let query = mol_from_smarts("[#8]").unwrap();
        let by_search = mol.has_substruct_match(&query).unwrap();
        let by_matcher = self.oxygen.has_match(mol).unwrap();
        assert_eq!(by_search, by_matcher);
        // A matcher that does not match, while the lock is already released.
        assert!(self.other.is_valid());
        self.other.has_match(mol).unwrap();
        by_matcher
    }

    fn get_matches(&self, mol: &Mol, matches: &mut Vec<FilterMatch>) -> bool {
        let found = self.oxygen.get_matches(mol).unwrap();
        let matched = !found.is_empty();
        matches.extend(found);
        matched
    }
}

// Matching holds control_c_mutex, for RDKit's SIGINT handler, and releases
// it while a Rust matcher runs, which may match itself. SubstructMatch writes
// a recursive query's matches into the query under its lock; copies of a
// pattern that others match are made under the same locks, after
// control_c_mutex. The functional groups' patterns are recursive and shared.
#[test]
fn threads_share_catalogs() {
    let reentrant = FilterMatcher::new(
        "reentrant",
        Reentrant {
            oxygen: smarts("O", "[#8]"),
            other: FilterMatcher::new(
                "at least",
                AtLeast {
                    atoms: 2,
                    calls: Arc::new(AtomicUsize::new(0)),
                },
            ),
        },
    );
    let mut catalog = FilterCatalog::from_catalogs(FilterCatalogs::PAINS_A);
    catalog.add_entry(&FilterCatalogEntry::new("reentrant", &reentrant).unwrap());
    let both = And::new(&reentrant, &smarts("C", "[#6]")).unwrap();
    catalog.add_entry(&FilterCatalogEntry::new("reentrant and C", &both).unwrap());
    let catalog = Arc::new(catalog);
    let groups = Arc::new(get_functional_group_hierarchy());
    let threads: Vec<_> = (0..4)
        .map(|i| {
            let catalog = catalog.clone();
            let groups = groups.clone();
            std::thread::spawn(move || {
                for _ in 0..10 {
                    assert_eq!(
                        descriptions(&catalog.get_matches(&mol("CCO")).unwrap()),
                        ["reentrant", "reentrant and C"]
                    );
                    assert!(!catalog.has_match(&mol("CC")).unwrap());
                    assert_eq!(
                        catalog
                            .get_first_match(&mol(HZONE))
                            .unwrap()
                            .unwrap()
                            .get_description()
                            .unwrap(),
                        "hzone_phenol_A(479)"
                    );
                    // RDKit's And replaces the matches gathered before it.
                    assert_eq!(
                        names(&catalog.get_filter_matches(&mol("CO")).unwrap()),
                        [
                            ("O".to_string(), vec![(0, 1)]),
                            ("C".to_string(), vec![(0, 0)])
                        ]
                    );
                    let rows =
                        run_filter_catalog_with_num_threads(&catalog, &["CCO", "CC", "O"], 4)
                            .unwrap();
                    let rows: Vec<_> = rows.iter().map(|row| descriptions(row)).collect();
                    assert_eq!(
                        rows,
                        [
                            vec!["reentrant", "reentrant and C"],
                            vec![],
                            vec!["reentrant"]
                        ]
                    );
                    if i % 2 == 0 {
                        assert_eq!(groups.get_matches(&mol("c1ccccc1Cl")).unwrap().len(), 1);
                    } else {
                        assert_eq!(get_flattened_functional_group_hierarchy().len(), 52);
                        let entry = groups.get_entry(0).unwrap();
                        assert_eq!(entry.get_description().unwrap(), "AcidChloride");
                    }
                }
            })
        })
        .collect();
    for thread in threads {
        thread.join().unwrap();
    }
}
