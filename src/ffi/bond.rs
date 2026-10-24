#[cxx::bridge(namespace = "rdkit_cxx")]
pub(crate) mod ffi {
    /// `RDKit::Bond::BondType`; cxx checks the values.
    #[repr(u8)]
    #[allow(clippy::upper_case_acronyms)]
    enum BondType {
        UNSPECIFIED = 0,
        SINGLE,
        DOUBLE,
        TRIPLE,
        QUADRUPLE,
        QUINTUPLE,
        HEXTUPLE,
        ONEANDAHALF,
        TWOANDAHALF,
        THREEANDAHALF,
        FOURANDAHALF,
        FIVEANDAHALF,
        AROMATIC,
        IONIC,
        HYDROGEN,
        THREECENTER,
        DATIVEONE,
        DATIVE,
        DATIVEL,
        DATIVER,
        OTHER,
        ZERO,
    }

    unsafe extern "C++" {
        include!("rdkit-cxx/errors.h");
        include!("rdkit-cxx/bond.h");

        type Bond;
        type BondType;

        fn bond_get_idx(bond: &Bond) -> Result<u32>;
        fn bond_get_bond_type(bond: &Bond) -> Result<BondType>;
        fn bond_get_bond_type_as_double(bond: &Bond) -> Result<f64>;
        fn bond_get_begin_atom_idx(bond: &Bond) -> Result<u32>;
        fn bond_get_end_atom_idx(bond: &Bond) -> Result<u32>;
        fn bond_get_other_atom_idx(bond: &Bond, this_idx: u32) -> Result<u32>;
        fn bond_get_is_aromatic(bond: &Bond) -> Result<bool>;
        fn bond_get_is_conjugated(bond: &Bond) -> Result<bool>;
    }
}
