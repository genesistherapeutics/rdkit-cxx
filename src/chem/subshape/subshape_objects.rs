// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Subshape.SubshapeObjects`, a Python module: the
//! shapes and skeleton points of subshapes.
//!
//! Python's `DisplaySubshape` and `DisplaySubshapeSkeleton` draw in PyMOL,
//! so they are not ported.

use crate::chem::feat_maps::feat_map_point::FeatMapPoint;
use crate::chem::feat_maps::feat_maps::FeatMap;
use crate::geometry::rd_geometry::{Point3D, UniformGrid3D};

/// A point of a shape's skeleton. Mirrors `SubshapeObjects.SkeletonPoint`.
#[derive(Clone, Debug, PartialEq)]
pub struct SkeletonPoint {
    /// The point. Mirrors `SkeletonPoint.location`.
    pub location: Point3D,
    /// The shape's moments around the point, largest first. Mirrors
    /// `SkeletonPoint.shapeMoments`.
    pub shape_moments: [f64; 3],
    /// The shape's principal directions around the point, in the moments'
    /// order; `None` until computed (Python's `[None, None, None]`).
    /// Mirrors `SkeletonPoint.shapeDirs`.
    pub shape_dirs: Option<[Point3D; 3]>,
    /// The families of the molecule's features near the point. Mirrors
    /// `SkeletonPoint.molFeatures`.
    pub mol_features: Vec<String>,
    /// The feature map's features near the point. Mirrors
    /// `SkeletonPoint.featmapFeatures`.
    pub featmap_features: Vec<String>,
    /// The fraction of the window's volume the shape fills. Mirrors
    /// `SkeletonPoint.fracVol`.
    pub frac_vol: f64,
}

impl Default for SkeletonPoint {
    fn default() -> SkeletonPoint {
        SkeletonPoint::from_location(Point3D::new(0.0, 0.0, 0.0))
    }
}

impl SkeletonPoint {
    /// A point at `location`. Python's has no location (`None`) without
    /// one. Mirrors `SkeletonPoint(location=location)`.
    pub fn from_location(location: Point3D) -> SkeletonPoint {
        SkeletonPoint {
            location,
            shape_moments: [0.0; 3],
            shape_dirs: None,
            mol_features: Vec::new(),
            featmap_features: Vec::new(),
            frac_vol: 0.0,
        }
    }
}

/// A shape's grid and skeleton. Mirrors `SubshapeObjects.ShapeWithSkeleton`.
#[derive(Clone, Debug, Default)]
pub struct ShapeWithSkeleton {
    /// The shape's encoding; `None` before one is made. Mirrors
    /// `ShapeWithSkeleton.grid`.
    pub grid: Option<UniformGrid3D>,
    /// The skeleton's points. Mirrors `ShapeWithSkeleton.skelPts`.
    pub skel_pts: Vec<SkeletonPoint>,
    /// The coarser grids the aligner adds to a target, as Python adds the
    /// attributes `medGrid` and `coarseGrid`.
    pub(crate) med_grid: Option<Box<ShapeWithSkeleton>>,
    pub(crate) coarse_grid: Option<Box<ShapeWithSkeleton>>,
}

/// A shape at several resolutions, with a feature map. Mirrors
/// `SubshapeObjects.SubshapeShape`.
#[derive(Clone, Debug, Default)]
pub struct SubshapeShape {
    /// The shape at increasing resolutions. Mirrors `SubshapeShape.shapes`.
    pub shapes: Vec<ShapeWithSkeleton>,
    /// Mirrors `SubshapeShape.featMap`.
    pub feat_map: Option<FeatMap>,
    /// Mirrors `SubshapeShape.keyFeat`.
    pub key_feat: Option<FeatMapPoint>,
}
