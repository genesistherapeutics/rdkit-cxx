// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! The recipes of the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)
//! in rdkit-cxx, one module per Cookbook chapter.
//!
//! Each recipe prints what the Cookbook prints and asserts the Cookbook's
//! output, so a recipe that stops matching RDKit fails. Run every recipe, or
//! name the ones to run by index:
//!
//! ```sh
//! cargo run --example cookbook
//! cargo run --example cookbook -- RDKitCB_3 RDKitCB_8
//! ```
//!
//! `cargo test --example cookbook` runs them all as a test.

mod bonds_and_bonding;
mod descriptor_calculations;
mod drawing_molecules;
mod error_messages;
mod manipulating_molecules;
mod miscellaneous_topics;
mod reactions;
mod rings_and_aromaticity;
mod stereochemistry;
mod substructure_matching;
mod writing_molecules;

use rdkit_cxx::Result;

/// A recipe: its Cookbook index, its title, and the function that runs it.
type Recipe = (&'static str, &'static str, fn() -> Result<()>);

/// Every chapter's recipes, in the Cookbook's order.
const CHAPTERS: &[&[Recipe]] = &[
    drawing_molecules::RECIPES,
    bonds_and_bonding::RECIPES,
    rings_and_aromaticity::RECIPES,
    stereochemistry::RECIPES,
    manipulating_molecules::RECIPES,
    substructure_matching::RECIPES,
    descriptor_calculations::RECIPES,
    writing_molecules::RECIPES,
    reactions::RECIPES,
    error_messages::RECIPES,
    miscellaneous_topics::RECIPES,
];

fn main() -> Result<()> {
    let selected: Vec<String> = std::env::args().skip(1).collect();
    for &(id, title, run) in CHAPTERS.iter().copied().flatten() {
        if selected.is_empty() || selected.iter().any(|name| name == id) {
            println!("== {id}: {title}");
            run()?;
            println!();
        }
    }
    Ok(())
}

#[test]
fn every_recipe_matches_the_cookbook() {
    for &(id, _, run) in CHAPTERS.iter().copied().flatten() {
        if let Err(err) = run() {
            panic!("{id}: {err}");
        }
    }
}
