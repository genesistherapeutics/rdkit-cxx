// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! On a Plot (RDKitCB_35): draw a molecule onto a plot.
//!
//! From the [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html#on-a-plot),
//! by Vincent Scalfani. The recipe puts the image on a matplotlib plot,
//! which has no counterpart here; this draws the image, which any plotting
//! library that places images can use.

mod common;

use common::{png_size, write_image};
use rdkit_cxx::chem::draw::mol_to_image;
use rdkit_cxx::chem::mol_from_smiles;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mol = mol_from_smiles("C1CNCCC1C(=O)C")?;
    let im = mol_to_image(&mol)?;
    write_image("on_a_plot.png", &im)?;
    assert_eq!(png_size(&im), (300, 300));
    Ok(())
}
