#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// A point in 3D. Mirrors `rdGeometry.Point3D`.
    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    struct Point3D {
        pub x: f64,
        pub y: f64,
        pub z: f64,
    }

    /// Atom or bond indices; nested lists cross as `Vec<IndexList>`.
    struct IndexList {
        indices: Vec<u32>,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/geometry.h");
    }

    // Other bridges return vectors of these.
    impl Vec<Point3D> {}
    impl Vec<IndexList> {}
}

impl ffi::IndexList {
    pub(crate) fn into_nested(lists: Vec<ffi::IndexList>) -> Vec<Vec<u32>> {
        lists.into_iter().map(|list| list.indices).collect()
    }
}
