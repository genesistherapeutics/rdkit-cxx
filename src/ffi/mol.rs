#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/mol.h");

        type RWMol;

        fn mol_new() -> Result<UniquePtr<RWMol>>;
        fn mol_clone(mol: &RWMol) -> Result<UniquePtr<RWMol>>;
        fn mol_get_num_atoms(mol: &RWMol) -> Result<u32>;
        fn mol_get_num_atoms_with_only_explicit(mol: &RWMol, only_explicit: bool) -> Result<u32>;
        fn mol_get_num_heavy_atoms(mol: &RWMol) -> Result<u32>;
        fn mol_get_num_bonds(mol: &RWMol, only_heavy: bool) -> Result<u32>;
    }
}

// SAFETY: an RWMol owns all of its data and has no thread affinity, so it can
// move between threads. It is not Sync: RDKit writes computed properties
// through const molecules (CONVENTIONS.md S1).
unsafe impl Send for ffi::RWMol {}
