#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        fn raise_for_tests(kind: &str) -> Result<()>;
    }
}
