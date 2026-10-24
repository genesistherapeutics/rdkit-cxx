#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/rd_base.h");
        fn rdkit_version() -> Result<UniquePtr<CxxString>>;
        fn boost_version() -> Result<UniquePtr<CxxString>>;
        fn disable_log(spec: &str) -> Result<()>;
        fn enable_log(spec: &str) -> Result<()>;
    }
}
