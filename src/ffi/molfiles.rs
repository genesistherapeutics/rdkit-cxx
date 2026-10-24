#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/molfiles.h");

        type RWMol = crate::ffi::mol::ffi::RWMol;

        fn mol_from_smiles(smiles: &str) -> Result<UniquePtr<RWMol>>;
        fn mol_to_smiles(mol: &RWMol) -> Result<UniquePtr<CxxString>>;
    }
}
