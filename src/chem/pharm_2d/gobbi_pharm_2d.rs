// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm2D.Gobbi_Pharm2D`, a Python module: the 2D
//! pharmacophores of Gobbi and Poppinger, Biotech. Bioeng. 61, 47-54
//! (1998).

use crate::chem::pharm_2d::sig_factory::{SigFactory, SigFactoryParams};
use crate::chem::rd_mol_chemical_features::build_feature_factory_from_string;

/// The feature definitions: lipophilic (`LH`), donor (`HD`), acceptor
/// (`HA`), aromatic (`AR`) and aliphatic (`RR`) attachment, unusual atom
/// (`X`), and basic (`BG`) and acidic (`AG`) group. Mirrors
/// `Gobbi_Pharm2D.fdef`.
pub fn fdef() -> &'static str {
    "
DefineFeature Hydrophobic [$([C;H2,H1](!=*)[C;H2,H1][C;H2,H1][$([C;H1,H2,H3]);!$(C=*)]),$(C([C;H2,H3])([C;H2,H3])[C;H2,H3])]
  Family LH
  Weights 1.0
EndFeature
DefineFeature Donor [$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]
  Family HD
  Weights 1.0
EndFeature
DefineFeature Acceptor [$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]
  Family HA
  Weights 1.0
EndFeature
DefineFeature AromaticAttachment [$([a;D3](@*)(@*)*)]
  Family AR
  Weights 1.0
EndFeature
DefineFeature AliphaticAttachment [$([A;D3](@*)(@*)*)]
  Family RR
  Weights 1.0
EndFeature
DefineFeature UnusualAtom [!#1;!#6;!#7;!#8;!#9;!#16;!#17;!#35;!#53]
  Family X
  Weights 1.0
EndFeature
DefineFeature BasicGroup [$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))]),$([N,n;X2;+0])]
  Family BG
  Weights 1.0
EndFeature
DefineFeature AcidicGroup [$([C,S](=[O,S,P])-[O;H1])]
  Family AG
  Weights 1.0
EndFeature
"
}

/// The distance bins: 2 to 8 bonds one by one, then 8 to 100. Mirrors
/// `Gobbi_Pharm2D.defaultBins`.
pub fn default_bins() -> Vec<(u32, u32)> {
    vec![(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 100)]
}

/// The Gobbi signature factory, of 2- and 3-point pharmacophores over
/// [`default_bins`], initialized. Python makes one, shared, on import;
/// this makes a new one each call. Mirrors `Gobbi_Pharm2D.factory`.
pub fn factory() -> SigFactory {
    let feat_factory =
        build_feature_factory_from_string(fdef()).expect("the Gobbi definitions parse");
    let params = SigFactoryParams {
        min_point_count: 2,
        max_point_count: 3,
        ..Default::default()
    };
    let mut factory = SigFactory::new_with_params(feat_factory, &params);
    factory
        .set_bins(&default_bins())
        .expect("the Gobbi factory initializes");
    factory
}
